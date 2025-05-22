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
    # Saves the workspace to "resources/"preprocessSampleMetadata"
    file.path("resources", paste0(snakemake@rule, ".RData")) |>
        save.image()
} else {
    # If the snakemake object does not exist, load the workspace
    file.path("resources", "preprocessMetadata") |>
        load()
}

# cleanCharacterStrings function from utils.R
snakemake@source("utils.R")

library(data.table)

# install readxl w pak if not installed
if (!requireNamespace("readxl", quietly = TRUE)) {
    install.packages("readxl")
}

# 1.0 Read Input Data
# -------------------
cat("Reading sample metadata file: ", INPUT$sampleMetadata, sep = "\n\t")
raw_sampleMetadata <- readxl::read_excel(INPUT$sampleMetadata, sheet = 1, col_names = TRUE, na = "NA")
raw_sampleMetadata <- as.data.table(raw_sampleMetadata)

cat("Reading treatment metadata file: ", INPUT$treatmentMetadata, sep = "\n\t")
treatmentMetadata <- fread(INPUT$treatmentMetadata)
names(treatmentMetadata) <- paste0("GDSC.", names(treatmentMetadata))
treatmentMetadata[, GDSC.treatmentid := AnnotationGx::cleanCharacterStrings(GDSC.DRUG_NAME)]


# 2.0 Preprocess Sample Metadata
# ------------------------------

cat("Cleaning sample metadata file", sep = "\n\t")
sampleMetadata <- copy(raw_sampleMetadata)
# clean column names using data.table
# remove \r and \n from column names
setnames(sampleMetadata, gsub("[\r\n]", "", names(sampleMetadata)))
setnames(sampleMetadata, gsub(" ", "_", names(sampleMetadata)))
# drop the row where  GDSC.Sample_Name == "TOTAL:""
sampleMetadata <- sampleMetadata[Sample_Name != "TOTAL:", ]

# Renames the columns of the sampleMetadata dataframe by adding a prefix "GDSC."
# to the column names that do not already start with "GDSC."
# The grepl() function is used to check if the column names start with "GDSC".
# The gsub() function is used to replace the "GDSC" prefix with "GDSC." if it already has.
names(sampleMetadata) <- ifelse(
    !grepl("^GDSC", names(sampleMetadata)),
    paste0("GDSC.", names(sampleMetadata)),
    gsub("^GDSC", "GDSC.", names(sampleMetadata))
)

# clean the sample names in the sampleMetadata dataframe
sampleMetadata[, GDSC.sampleid := cleanCharacterStrings(GDSC.Sample_Name)]

sampleMetadata[, GDSC.COSMIC_ID := as.numeric(GDSC.COSMIC_identifier)][, GDSC.COSMIC_identifier := NULL]
# CMP_sampleAnnotation[, GDSC.COSMIC_ID := as.numeric(GDSC.COSMIC_ID)]
# sampleMetadata <- merge(sampleMetadata, CMP_sampleAnnotation[, !c("GDSC.sampleid")], by = "GDSC.COSMIC_ID", all.x = TRUE)

# get all the column names, and reorder so that GDSC.sampleid is the first column and CCLE.sampleid is the second column
firstcols <- c("GDSC.sampleid", "GDSC.Sample_Name")
columnNames <- c(firstcols, names(sampleMetadata)[!names(sampleMetadata) %in% firstcols])
sampleMetadata <- sampleMetadata[, columnNames, with = FALSE]

cat("Writing sample metadata file: ", OUTPUT$sampleMetadata, sep = "\n\t")
fwrite(sampleMetadata, OUTPUT$sampleMetadata, sep = "\t", quote = FALSE, na = "NA", row.names = FALSE, col.names = TRUE)

# 3.0 Preprocess Treatment Metadata
# ------------------------------
# get all the column names, and reorder so that GDSC.treatment is the first column and CCLE.treatment is the second column
treatmentMetadata <- 
    treatmentMetadata[, 
    c("GDSC.treatmentid", setdiff(names(treatmentMetadata), c("GDSC.treatmentid"))),
    with = FALSE]

data.table::fwrite(
    treatmentMetadata,
    OUTPUT$treatmentMetadata,
    sep = "\t",
    quote = FALSE,
    na = "NA",
    row.names = FALSE,
    col.names = TRUE
)