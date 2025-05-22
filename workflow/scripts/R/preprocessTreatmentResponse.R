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

    # # Assuming that this script is named after the rule
    # # Saves the workspace to "resources/"preprocessTreatmentResponse"
    # file.path("resources", paste0(snakemake@rule, ".RData")) |>
    #     save.image()
} else {
    # If the snakemake object does not exist, load the workspace
    file.path("resources", "preprocessTreatmentResponse.RData") |>
        load()
}

snakemake@source("utils.R")
library(data.table)
###############################################################################
# Load INPUT
###############################################################################
rawdata <- data.table::fread(INPUT$rawdata)
processed <- data.table::as.data.table(readxl::read_xlsx(INPUT$profiles))
treatmentMetadata <- data.table::fread(INPUT$treatmentMetadata)
sampleMetadata <- data.table::fread(INPUT$sampleMetadata)


# 2.0 Subset data
# ---------------

rawdata_s <- merge(
    rawdata,
    sampleMetadata,
    by.x = "COSMIC_ID",
    by.y = "GDSC.COSMIC_ID",
)

rawdata_st <- merge(
    rawdata_s,
    treatmentMetadata,
    by.x = "DRUG_ID",
    by.y = "GDSC.DRUG_ID",
    all.x = TRUE
)

subsetted_rawdata <- rawdata_st[(!is.na(DRUG_ID) & !is.na(GDSC.treatmentid)) | grepl("^(L|R|A|N|B)", TAG)]

# drop duplicate rows
subsetted_rawdata <- subsetted_rawdata[!duplicated(subsetted_rawdata)]
subsetted_rawdata

###############################################################################
# Main Script
###############################################################################

## ------------------- GDSC SPECIFIC STEP ------------------- ##
# In the essence of using the processed GDSC data as much as possible
# we will normalize the rawdata using the gdscIC50 package
# neg_control_TAGS <- subsetted_rawdata[grepl("^NC", TAG), unique(TAG)]
neg_control_tag <- ifelse(WILDCARDS$version == "GDSC1", "NC-0", "NC-1")

subsetted_rawdata <- gdscIC50::removeFailedDrugs(subsetted_rawdata)

subsetted_rawdata <- gdscIC50::removeMissingDrugs(subsetted_rawdata)

normData <- gdscIC50::normalizeData(subsetted_rawdata, trim = TRUE, neg_control = neg_control_tag, pos_control = "B")
## ------------------------------------------------------------------------- ##
# Do something
# the normalization removes the treatmentid and sampleid columns
normData <- merge(normData, unique(subsetted_rawdata[, .(MASTER_CELL_ID, GDSC.sampleid)]), by = "MASTER_CELL_ID", all.x = T)
normData <- merge(normData, unique(subsetted_rawdata[, .(DRUG_ID, GDSC.treatmentid)]), by.x = "DRUG_ID_lib", by.y = "DRUG_ID", all.x = T)
normData <- merge(normData, treatmentMetadata, by.x = c("GDSC.treatmentid","DRUG_ID_lib"), by.y = c("GDSC.treatmentid","GDSC.DRUG_ID"), all.x = T)

# remove anything with empty GDSC.sampleid, then rename CONC to Dose and normalized_intensity to Viability
x_assay <- normData[!is.na(GDSC.sampleid), ] |>
    data.table::setnames(
        old = c("CONC", "normalized_intensity", "dose"),
        new = c("Dose", "Viability", "DoseNumber")
    )

print(paste0("Merging treatmentMetadata with processed data"))
procData <- merge(processed, treatmentMetadata, by.x = "DRUG_ID", by.y = "GDSC.DRUG_ID", all.x = T)
procData <- procData[, GDSC.sampleid := cleanCharacterStrings(CELL_LINE_NAME)][!is.na(GDSC.sampleid) & GDSC.sampleid %in% sampleMetadata$GDSC.sampleid]


# 3.0 Subset data
# -------------------
print(paste0("Subsetting normData to not include missing values"))
subset_normData <- unique(x_assay[
    !is.na(Viability) &
        !is.na(GDSC.treatmentid) &
        !is.na(Dose) &
        !is.na(GDSC.sampleid)
])
subset_normData
# CoreGx Functions require that the main terminology is "treatmentid" and "sampleid"
# we will add these columns to the subset_normData

# rename the columns to match the CoreGx terminology
data.table::setnames(
    subset_normData,
    old = c("GDSC.treatmentid", "GDSC.sampleid"),
    new = c("treatmentid", "sampleid")
)

unique_treatments <- unique(subset_normData$treatmentid)
unique_samples <- unique(subset_normData$sampleid)
print(paste0("Number of unique treatments: ", length(unique_treatments)))
print(paste0("Number of unique samples: ", length(unique_samples)))

pre_raw <- OUTPUT$preprocessed_raw
pre_profiles <- OUTPUT$preprocessed_profiles

# ###############################################################################
# # # Save OUTPUT
# ###############################################################################
message("Saving the preprocessed raw data...")
data.table::fwrite(
    subset_normData,
    file = pre_raw,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
)
message("Saving the preprocessed profiles data...")
data.table::fwrite(
    procData,
    file = pre_profiles,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
)

# # TRE
# raw_dt <- subset_normData
# message(paste0("Loading TREDataMapper"))
# # subset raw_dt to only the first 100 sampleids

# treDataMapper <- CoreGx::TREDataMapper(rawdata = raw_dt)

# # > names(raw_dt)
# #  [1] "DRUG_ID_lib"      "MASTER_CELL_ID"   "SCAN_ID"          "BARCODE"
# #  [5] "RESEARCH_PROJECT" "DATE_CREATED"     "DRUGSET_ID"       "CELL_LINE_NAME"
# #  [9] "CELL_ID"          "COSMIC_ID"        "POSITION"         "Dose"
# # [13] "INTENSITY"        "lib_drug"         "DoseNumber"       "treatment"
# # [17] "NC"               "PC"               "Viability"        "norm_neg_pos"
# # [21] "time_stamp"       "sw_version"       "GDSC.sampleid"    "GDSC.treatmentid"
# # [25] "treatmentid"      "sampleid"
# # >

# # We need to define which columns identify the drug
# # and which columns identify the sample
# groups <- list(
#     drugs_doses2 = c("treatmentid", "Dose", "lib_drug"),
#     sample_ids = c("sampleid", "BARCODE"),
#     assay2 = c("treatmentid", "Dose", "lib_drug", "sampleid", "BARCODE")
# )

# # subset out mapped?

# subsets <- c(TRUE)

# guess <- CoreGx::guessMapping(
#     treDataMapper,
#     groups = groups,
#     subset = subsets
# )

# CoreGx::rowDataMap(treDataMapper) <- list(
#     id_columns = guess$drugs_doses2$id_columns,
#     mapped_columns = guess$drugs_doses2$mapped_columns
# )

# CoreGx::colDataMap(treDataMapper) <- list(
#     id_columns = guess$sample_ids$id_columns,
#     mapped_columns = guess$sample_ids$mapped_columns
# )

# CoreGx::assayMap(treDataMapper) <- list(
#     sensitivity = list(
#         id_columns = guess$assay2$id_columns,
#         mapped_columns = guess$assay2$mapped_columns
#     )
# )

# show(treDataMapper)

# message("Running CoreGx::metaConstruct")
# tre <- CoreGx::metaConstruct(treDataMapper)
# CoreGx::metadata(tre) <- list(
#     data_source = snakemake@config$treatmentResponse,
#     annotation = "treatmentResponse",
#     date = Sys.Date(),
#     sessionInfo = capture.output(sessionInfo())
# )
# show(tre)
# ###############################################################################
# # Save OUTPUT
# ###############################################################################
# message("Saving the treatment response experiment object...")
# saveRDS(
#     tre_fit,
#     file = OUTPUT$tre
# )
