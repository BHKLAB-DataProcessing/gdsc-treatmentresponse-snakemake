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
    # Saves the workspace to "resources/"buildTreatmentResponseExperiment"
    file.path("resources", paste0(snakemake@rule, ".RData")) |>
        save.image()
} else {
    # If the snakemake object does not exist, load the workspace
    file.path("resources", "buildTreatmentResponseExperiment.RData") |>
        load()
}


###############################################################################
# Load INPUT
###############################################################################
pre_raw <- data.table::fread(INPUT$preprocessed_raw)
pre_profiles <- data.table::fread(INPUT$preprocessed_profiles)


# TRE
raw_dt <- pre_raw
message(paste0("Loading TREDataMapper"))
# subset raw_dt to only the first 100 sampleids

treDataMapper <- CoreGx::TREDataMapper(rawdata = raw_dt)

# > names(raw_dt)
#  [1] "DRUG_ID_lib"      "MASTER_CELL_ID"   "SCAN_ID"          "BARCODE"
#  [5] "RESEARCH_PROJECT" "DATE_CREATED"     "DRUGSET_ID"       "CELL_LINE_NAME"
#  [9] "CELL_ID"          "COSMIC_ID"        "POSITION"         "Dose"
# [13] "INTENSITY"        "lib_drug"         "DoseNumber"       "treatment"
# [17] "NC"               "PC"               "Viability"        "norm_neg_pos"
# [21] "time_stamp"       "sw_version"       "GDSC.sampleid"    "GDSC.treatmentid"
# [25] "treatmentid"      "sampleid"
# >

# We need to define which columns identify the drug
# and which columns identify the sample
groups <- list(
    drugs_doses2 = c("treatmentid", "Dose", "lib_drug"),
    sample_ids = c("sampleid", "BARCODE"),
    assay2 = c("treatmentid", "Dose", "lib_drug", "sampleid", "BARCODE")
)

# subset out mapped?

subsets <- c(TRUE)

guess <- CoreGx::guessMapping(
    treDataMapper,
    groups = groups,
    subset = subsets
)

CoreGx::rowDataMap(treDataMapper) <- list(
    id_columns = guess$drugs_doses2$id_columns,
    mapped_columns = guess$drugs_doses2$mapped_columns
)

CoreGx::colDataMap(treDataMapper) <- list(
    id_columns = guess$sample_ids$id_columns,
    mapped_columns = guess$sample_ids$mapped_columns
)

CoreGx::assayMap(treDataMapper) <- list(
    sensitivity = list(
        id_columns = guess$assay2$id_columns,
        mapped_columns = guess$assay2$mapped_columns
    )
)

show(treDataMapper)

message("Running CoreGx::metaConstruct")
tre <- CoreGx::metaConstruct(treDataMapper)
CoreGx::metadata(tre) <- list(
    data_source = snakemake@config$treatmentResponse,
    annotation = "treatmentResponse",
    date = Sys.Date(),
    sessionInfo = capture.output(sessionInfo())
)
show(tre)
###############################################################################
# Save OUTPUT
###############################################################################
message("Saving the treatment response experiment object...")
saveRDS(
    tre_fit,
    file = OUTPUT$tre
)
