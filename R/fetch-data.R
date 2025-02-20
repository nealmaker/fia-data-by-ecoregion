###############################
# STEP 0: LOAD REQUIRED LIBRARIES
###############################

# install.packages("sf")          # for spatial ops
# install.packages("dplyr")       # data wrangling
# install.packages("stringr")     # pattern matching
# install.packages("purrr")       # functional map
# install.packages("tidyr")       # data reshaping
# install.packages("readr")       # reading CSV
# install.packages("utils")       # for download.file (usually installed by default)

library(sf)
library(dplyr)
library(stringr)
library(purrr)
library(tidyr)
library(readr)

###############################
# HELPER: DOWNLOAD FIA CSVs
###############################

download_fia_csv <- function(states, tables = c("PLOT", "COND", "TREE"), out_dir = "FIA_Data") {
  # Create the output directory if it doesn't exist
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  base_url <- "https://apps.fs.usda.gov/fia/datamart/CSV"
  
  # For each state, for each table, try to download the ZIP and extract
  for (st in states) {
    for (tb in tables) {
      zip_file <- file.path(out_dir, paste0(st, "_", tb, ".zip"))
      csv_file <- file.path(out_dir, paste0(st, "_", tb, ".csv"))
      
      # If we haven't downloaded or if you want to refresh, proceed:
      # (You might want better checks e.g. file.exists(zip_file) or a user-defined "force = TRUE" param)
      file_url <- paste0(base_url, "/", st, "_", tb, ".zip")
      message("Attempting to download: ", file_url)
      
      utils::download.file(url = file_url, destfile = zip_file, mode = "wb", quiet = FALSE)
      
      # Unzip. The CSV inside typically has the same name as <STATE>_<TABLE>.csv
      utils::unzip(zipfile = zip_file, exdir = out_dir)
      
      # Optional: you could remove the ZIP after extraction
      # file.remove(zip_file)
    }
  }
}

###############################
# MAIN FUNCTION
###############################

create_FIA_dataset <- function(
    ecoregion_name,
    ecoregion_shapefile_path,
    states_to_download,
    # Directory where CSVs will be stored/loaded:
    csv_dir = "FIA_Data",
    # Do we want to attempt fresh downloads?
    do_download = TRUE
) {
  # 1) Optionally download the relevant FIA CSV files for each state
  if (do_download) {
    # You can add more tables if needed (e.g. POP_ESTN_UNIT, POP_EVAL, etc.)
    download_fia_csv(states = states_to_download, tables = c("PLOT", "COND", "TREE"), out_dir = csv_dir)
  }
  
  # 2) Read the CSVs for each state and combine
  # We'll read PLOT, TREE, and COND. 
  # Adjust column names if your version differs.
  
  # Helper function to read a single table for multiple states
  read_fia_table <- function(states, table_name, col_types = cols(.default = "c"), path_dir) {
    # Read all states, bind rows
    purrr::map_dfr(states, function(st) {
      csv_file <- file.path(path_dir, paste0(st, "_", table_name, ".csv"))
      # If needed, refine col_types to match known data types
      readr::read_csv(csv_file, col_types = col_types)
    })
  }
  
  # Read the combined tables
  plot_df <- read_fia_table(states_to_download, "PLOT", path_dir = csv_dir)
  cond_df <- read_fia_table(states_to_download, "COND", path_dir = csv_dir)
  tree_df <- read_fia_table(states_to_download, "TREE", path_dir = csv_dir)
  
  # 3) Load the shapefile for the ecoregions and filter to the requested L2 ecoregion
  ecoregions_sf <- st_read(ecoregion_shapefile_path, quiet = TRUE)
  # Suppose the shapefile column with Level II name is "NA_L2NAME" (often the case).
  # Adjust as needed if your shapefile is labeled differently.
  ecoregion_of_interest <- ecoregions_sf %>%
    filter(str_detect(NA_L2NAME, regex(ecoregion_name, ignore_case = TRUE))) %>%
    st_make_valid()
  
  # 4) Convert plot coordinates to sf and intersect with ecoregion
  #    The FIA columns might be "LAT" and "LON" or "LAT_DD" / "LON_DD".
  #    We'll assume "LAT" and "LON" exist. Adjust if needed.
  plot_df_sf <- plot_df %>%
    # Filter out rows with missing lat/long
    filter(!is.na(LAT) & !is.na(LON)) %>%
    # Convert to sf, WGS84 (EPSG:4326)
    st_as_sf(coords = c("LON", "LAT"), crs = 4326)
  
  # Intersect with the ecoregion
  plots_in_ecoregion_sf <- st_intersection(plot_df_sf, ecoregion_of_interest)
  
  # Grab these plot CNs
  plot_cns_in_ecoregion <- plots_in_ecoregion_sf$CN
  
  # 5) Filter the big data frames to only include those plot CNs
  plot_filtered <- plot_df %>%
    filter(CN %in% plot_cns_in_ecoregion)
  
  cond_filtered <- cond_df %>%
    filter(PLT_CN %in% plot_cns_in_ecoregion)
  
  tree_filtered <- tree_df %>%
    filter(PLT_CN %in% plot_cns_in_ecoregion)
  
  # 6) Build remeasurement intervals
  #    Each PLOT CN can appear in multiple inventory years (INVYR).
  #    We'll keep distinct combos of CN x INVYR (and MEASYEAR if present).
  
  # If your CSV uses "INVYR" or "YEAR", confirm the correct column name. 
  # We'll assume "INVYR" is correct.
  plot_measurements <- plot_filtered %>%
    select(CN, INVYR, MEASYEAR) %>% 
    distinct() %>%
    arrange(CN, INVYR)
  
  # Self-join to create all intervals for each plot
  intervals <- plot_measurements %>%
    rename(CN_1 = CN, INVYR_1 = INVYR, MEASYEAR_1 = MEASYEAR) %>%
    inner_join(
      plot_measurements %>% rename(CN_2 = CN, INVYR_2 = INVYR, MEASYEAR_2 = MEASYEAR),
      by = c("CN_1" = "CN_2")
    ) %>%
    filter(INVYR_2 > INVYR_1) %>%
    arrange(CN_1, INVYR_1, INVYR_2)
  
  # 7) Merge tree data from time 1 and time 2
  #    Usually you match on the same plot (PLT_CN), same SUBP, same TREE for the same tree.
  
  # We'll create an identifier for each measurement event
  tree_filtered_keyed <- tree_filtered %>%
    mutate(plot_meas_id = paste0(PLT_CN, "_", INVYR)) %>%
    select(plot_meas_id, PLT_CN, INVYR, SUBP, TREE, SPCD, DBH, STATUSCD)  # columns you need
  
  # Label time columns
  label_time <- function(df, suffix) {
    # rename DBH -> DBH_t1 or DBH_t2, etc.
    df %>%
      rename_with(
        .cols = c("DBH", "STATUSCD", "TREE", "SUBP", "SPCD"),
        .fn = ~ paste0(.x, "_", suffix)
      )
  }
  
  # Time1 trees
  time1 <- tree_filtered_keyed %>%
    rename(PLT_CN_1 = PLT_CN, INVYR_1 = INVYR) %>%
    mutate(plot_meas_id_1 = plot_meas_id) %>%
    select(-plot_meas_id) %>%
    label_time("t1")
  
  # Time2 trees
  time2 <- tree_filtered_keyed %>%
    rename(PLT_CN_2 = PLT_CN, INVYR_2 = INVYR) %>%
    mutate(plot_meas_id_2 = plot_meas_id) %>%
    select(-plot_meas_id) %>%
    label_time("t2")
  
  # Build final wide intervals
  # For each row in intervals, we merge the T1 and T2 trees by the same SUBP, TREE, etc.
  intervals_tree_level <- intervals %>%
    mutate(
      plot_meas_id_1 = paste0(CN_1, "_", INVYR_1),
      plot_meas_id_2 = paste0(CN_1, "_", INVYR_2)
    ) %>%
    # Join time1
    left_join(time1, by = c("plot_meas_id_1")) %>%
    # Then join time2
    left_join(
      time2,
      by = c("plot_meas_id_2" = "plot_meas_id_2",
             "CN_1" = "PLT_CN_2")  # often you'd also match SUBP_t1==SUBP_t2, TREE_t1==TREE_t2
    ) %>%
    # Possibly refine the join to ensure the same tree: 
    #   left_join(time2, 
    #             by = c("plot_meas_id_2", "SUBP_t1"="SUBP_t2", "TREE_t1"="TREE_t2"))
    select(CN_1, INVYR_1, MEASYEAR_1, CN_2, INVYR_2, MEASYEAR_2, everything())
  
  # The result is a data frame with all intervals for each plot in the ecoregion.
  # Each interval row can have many tree-level combos. (Plot-level intervals are repeated for each tree.)
  # You can compute growth, mortality, ingrowth, etc. from these wide columns or in subsequent steps.
  
  # Return the final assembled dataset
  return(intervals_tree_level)
}

###############################
# EXAMPLE USAGE
###############################
# 
# 1. Suppose you want ecoregion: "Mediterranean California" from the shapefile
#    that has a column "NA_L2NAME". 
# 2. You suspect it intersects CA, OR, NV.
# 3. You downloaded the "Level II Ecoregions of North America Shapefile" and unzipped
#    it to "data/Ecoregions/NA_CEC_Eco_Level2.shp".
# 4. Then you can do:

# states_of_interest <- c("CA", "OR", "NV")
# ecoregion_shp_path <- "data/Ecoregions/NA_CEC_Eco_Level2.shp"
# ecoregion_of_interest_name <- "Mediterranean California"

# final_dataset <- create_FIA_dataset(
#   ecoregion_name            = ecoregion_of_interest_name,
#   ecoregion_shapefile_path  = ecoregion_shp_path,
#   states_to_download        = states_of_interest,
#   csv_dir                   = "FIA_Data",
#   do_download               = TRUE   # set to FALSE if you already downloaded
# )
#
# head(final_dataset)

