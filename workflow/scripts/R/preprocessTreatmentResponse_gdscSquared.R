## ------------------- Parse Snakemake Object ------------------- ##
# Check if the "snakemake" object exists
# This snippet is run at the beginning of a snakemake run to setup the env
# Helps to load the workspace if the script is run independently or debugging
if (exists("snakemake")) {
    INPUT <- snakemake@input
    OUTPUT <- snakemake@output
    WILDCARDS <- snakemake@wildcards
    THREADS <- snakemake@threads

    # setup logger if log file is provided
    if (length(snakemake@log) > 0) {
        sink(
            file = snakemake@log[[1]],
            append = FALSE,
            type = c("output", "message"),
            split = TRUE
        )
    }

    # Assuming that this script is named after the rule
    # Saves the workspace to "resources/"preprocessTreatmentResponse_gdscSquared"
    file.path("resources", paste0(snakemake@rule, ".RData")) |>
        save.image()
} else {
    # If the snakemake object does not exist, load the workspace
    file.path("resources", "preprocessTreatmentResponse_gdscSquared.RData") |>
        load()
}
# IMPORTS
# pak::pkg_install('CancerRxGene/gdscmatrixanalyser')
required <- c("CancerRxGene/gdscmatrixanalyser", "CancerRxGene/gdscIC50", "data.table", "cli")
# if (!requireNamespace("gdscmatrixanalyser", quietly = TRUE)) {
#     remotes::install_github("CancerRxGene/gdscmatrixanalyser")
#     library(gdscmatrixanalyser)
# }
if (!requireNamespace("pak", quietly = TRUE)) {
    install.packages("pak")
}
pak::pkg_install(required)
library(gdscmatrixanalyser)
library(data.table)
library(cli)
library(dplyr)
library(tidyr)
library(purrr)
###############################################################################
# Load INPUT
###############################################################################
original_dir <- INPUT$original_dir

# assign the read in data to the filename (no extension)
# rawdata <- fs::dir_ls(original_dir, glob = "*raw_data*") |>
#     # assign clean file-stems as names
#     purrr::set_names(~ fs::path_ext_remove(fs::path_file(.x))) |>
#     # read each CSV; silence the spec messages if you like
#     purrr::map(readr::read_csv,
#         col_names = TRUE,
#         show_col_types = FALSE
#     )

rawdata <- fs::dir_ls(original_dir, glob = "**Original*raw_data*") |>
    data.table::fread(
        sep = ",",
        header = TRUE,
    )

validation <- fs::dir_ls(original_dir, glob = "**Validation*raw_data**") |>
    data.table::fread(
        sep = ",",
        header = TRUE,
    )
###############################################################################
# Main Script
###############################################################################
# steps <- c(
#     "Create drugset list",
#     "Create plate list",
# )


cli::cli_h1("Preprocess Treatment Response GDSC Squared")

cli::cli_inform("Creating drugset list")
drugset_list <- rawdata |> make_drugset_list()
# show done
cli::cli_alert_success("Drugset list created")


cli::cli_inform("Creating plate list")
plate_list <- rawdata |> make_plate_list()
cli::cli_alert_success("Plate list created")

# Process screen data to extract treatment information from drug set 269
# Extract, clean, and format position and tag data
process_one <- function(drugset_sub_dt) {
    ds_id <- unique(drugset_sub_dt$DRUGSET_ID)
    stopifnot(length(ds_id) == 1)
    cli::cli_alert_info("Processing drugset {ds_id}")
    processed <- drugset_sub_dt |>
        distinct(POSITION, TAG, DRUG_ID, CONC) |>
        mutate(CONC = round(CONC, digits = 5)) |>
        # Standardize TAG names
        mutate(
            TAG = ifelse(TAG == "NC-1", "NC1", TAG),
            TAG = ifelse(TAG == "NC-0", "NC0", TAG),
            TAG = ifelse(TAG == "UN-USED", "UNUSED", TAG)
        ) |>
        filter(TAG != "DMSO") |>
        arrange(POSITION, TAG) |>
        # Parse TAG components into separate columns
        separate(TAG, into = c("lib", "dose", "treatment"), sep = "-", extra = "merge", fill = "right") |>
        # Handle special case: when dose is "C" or "S", it's actually a treatment type (not a dose),
        # so move it to the treatment column and set dose to NA
        mutate(
            treatment = ifelse(dose %in% c("C", "S"), dose, treatment),
            dose = ifelse(dose %in% c("C", "S"), NA, dose)
        ) |>
        # Combine treatment information into a single identifier
        unite(col = "treatment", lib, dose, treatment, DRUG_ID, CONC, sep = "_") |>
        # Process multiple treatments per position
        # Some positions have multiple treatments (e.g., drug combinations)
        # We need to transform these multiple rows per position into columns
        # this will HOPEFULLY give TWO unique tments (anchor AND treatment)
        # A1_NA_C_1005_4      L1_D1_C_1058_4
        group_by(POSITION) |>
        mutate(tment = 1:n()) |> # Number each treatment per position (1,2,3...)
        ungroup() |>
        mutate(tment = paste("t", tment, sep = "")) |> # Create column names like "t1", "t2", etc.

        # Reshape from long to wide format - transform multiple treatments from rows to columns
        # This converts data from:
        #   POSITION | treatment
        #   A01      | treatment_info_1
        #   A01      | treatment_info_2
        # To:
        #   POSITION | t1              | t2
        #   A01      | treatment_info_1| treatment_info_2
        spread(key = tment, value = treatment) |>
        # Extract components from treatment columns
        separate(t1, into = c("lib1", "lib1_dose", "treatment1", "lib1_drug_id", "lib1_conc"), sep = "_") |>
        separate(t2, into = c("lib2", "lib2_dose", "treatment2", "lib2_drug_id", "lib2_conc"), sep = "_") |>
        # convert concentration columns to numeric
        mutate(
            lib1_conc = suppressWarnings(as.numeric(lib1_conc)),
            lib2_conc = suppressWarnings(as.numeric(lib2_conc))
        )
    # Convert to data.table format

    if (nrow(processed) == 0) {
        cli::cli_alert_warning("No data found for drugset {ds_id}")
        return(NULL)
    }
    if (processed %>% filter(treatment1 == "C" & treatment2 == "C") %>% nrow() > 0) {
        ds_layout <- processed %>%
            left_join(
                {
                    .
                } %>%
                    filter(treatment1 == "C" & treatment2 == "C") %>%
                    distinct(lib1, lib2) %>%
                    mutate(cmatrix = 1:n()),
                by = c("lib1", "lib2")
            )
    } else {
        ds_layout <- processed %>% mutate(cmatrix = NA)
    }
    return(drugset(
        drugset_id = ds_id,
        layout = ds_layout
    ))
}



# # Map over each drugset ID and process
drugset_list <- unique(rawdata$DRUGSET_ID) |>
    purrr::set_names(~ paste0("ds", .x)) |>
    purrr::map(\(id) {
        rawdata |>
            filter(DRUGSET_ID == id) |>
            process_one()
    })
# drugset_list
#
xx <- lapply(plate_list, set_plate_controls)


my_plate <- plate_list[[3]]
my_drugset <- gdscmatrixanalyser:::get_drugset(my_plate)
my_drugset
# set_plate_controls(my_plate)
set_plate_matrix_list(my_plate)
# get_monotherapies(my_plate, 'hi')

# # Create new drugset from processed data
# new_drugset <- processed_data

## ------------------------------------------------------------------------- ##
# Do something with the processed screen data

# For example, you might want to create a new drugset object

###############################################################################
# Save OUTPUT
###############################################################################
# ds_layout <- screen_data %>%
#     distinct(POSITION, TAG, DRUG_ID, CONC) %>%
#     mutate(CONC = round(CONC, digits = 8)) %>%
#     mutate(TAG = ifelse(TAG == "NC-1", "NC1", TAG)) %>%
#     mutate(TAG = ifelse(TAG == "NC-0", "NC0", TAG)) %>%
#     mutate(TAG = ifelse(TAG == "UN-USED", "UNUSED", TAG)) %>%
#     filter(TAG != "DMSO") %>%
#     arrange(POSITION, TAG) %>%
#     separate(TAG, into = c("lib", "dose", "treatment"), sep = "-", extra = "merge", fill = "right") %>%
#     mutate(treatment = ifelse(dose %in% c("C", "S"), dose, treatment)) %>%
#     mutate(dose = ifelse(dose %in% c("C", "S"), NA, dose)) %>%
#     unite(col = "treatment", lib, dose, treatment, DRUG_ID, CONC, sep = "_") %>%
#     group_by(POSITION) %>%
#     mutate(tment = 1:n()) %>%
#     ungroup() %>%
#     mutate(tment = paste("t", tment, sep = "")) %>%
#     spread(key = tment, value = treatment) %>%
#     separate(t1, into = c("lib1", "lib1_dose", "treatment1", "lib1_drug_id", "lib1_conc"), sep = "_") %>%
#     separate(t2, into = c("lib2", "lib2_dose", "treatment2", "lib2_drug_id", "lib2_conc"), sep = "_") %>%
#     mutate(
#         lib1_conc = suppressWarnings(as.numeric(lib1_conc)),
#         lib2_conc = suppressWarnings(as.numeric(lib2_conc))
#     ) %>%
#     # TEST THAT treatment1 == treatment2 unless NA...tbc - should be redundant
#     mutate(treatment = treatment1) %>%
#     select(-treatment1, -treatment2)
