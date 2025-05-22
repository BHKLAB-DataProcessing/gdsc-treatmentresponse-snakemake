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
    # Saves the workspace to "resources/"preprocessTreatmentResponse"
    file.path("resources", paste0(snakemake@rule, ".RData")) |>
        save.image()
} else {
    # If the snakemake object does not exist, load the workspace
    file.path("resources", "preprocessTreatmentResponse.RData") |>
        load()
}

snakemake@source("utils.R")
###############################################################################
# Load INPUT
###############################################################################
rawdata <- fread(INPUT$rawdata)
processed <- as.data.table(readxl::read_xlsx(INPUT$profiles))
treatmentMetadata <- fread(INPUT$treatmentMetadata)
sampleMetadata <- fread(INPUT$sampleMetadata)


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

# rawdata_s <- rawdata[, GDSC.sampleid := ..sampleMetadata[match(SANGER_MODEL_ID, CMP.model_id), GDSC.sampleid]]

# rawdata_st <- rawdata_s[, GDSC.treatmentid := ..treatmentMetadata[match(DRUG_ID, GDSC.DRUG_ID), GDSC.treatmentid]]

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
if (!require("gdscIC50", quietly = TRUE)) {
    pak::pkg_install("cancerrxgene/gdscIC50")
}


subsetted_rawdata <- gdscIC50::removeFailedDrugs(subsetted_rawdata)

subsetted_rawdata <- gdscIC50::removeMissingDrugs(subsetted_rawdata)

normData <- gdscIC50::normalizeData(subsetted_rawdata, trim = TRUE, neg_control = neg_control_tag, pos_control = "B")
## ------------------------------------------------------------------------- ##
# Do something
# the normalization removes the treatmentid and sampleid columns
normData <- merge(normData, unique(subsetted_rawdata[, .(MASTER_CELL_ID, GDSC.sampleid)]), by = "MASTER_CELL_ID", all.x = T)
normData <- merge(normData, unique(subsetted_rawdata[, .(DRUG_ID, GDSC.treatmentid)]), by.x = "DRUG_ID_lib", by.y = "DRUG_ID", all.x = T)

# remove anything with empty GDSC.sampleid, then rename CONC to Dose and normalized_intensity to Viability
x_assay <- normData[!is.na(GDSC.sampleid), ] |>
    data.table::setnames(
        old = c("CONC", "normalized_intensity", "dose"),
        new = c("Dose", "Viability", "DoseNumber")
    )

print(paste0("Merging treatmentMetadata with processed data"))
procData <- merge(processed, treatmentMetadata[, .(GDSC.DRUG_ID, GDSC.treatmentid)], by.x = "DRUG_ID", by.y = "GDSC.DRUG_ID", all.x = T)
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

# CoreGx Functions require that the main terminology is "treatmentid" and "sampleid"
# we will add these columns to the subset_normData
subset_normData[, "treatmentid" := GDSC.treatmentid]
subset_normData[, "sampleid" := GDSC.sampleid]

unique_treatments <- unique(subset_normData$GDSC.treatmentid)
unique_samples <- unique(subset_normData$GDSC.sampleid)
print(paste0("Number of unique treatments: ", length(unique_treatments)))
print(paste0("Number of unique samples: ", length(unique_samples)))


###############################################################################
# Save OUTPUT
###############################################################################
