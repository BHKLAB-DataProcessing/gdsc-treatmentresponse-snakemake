cleanCharacterStrings <- function(name){

    # # make sure name is a string
    name <- as.character(name)

    # if there is a colon like in "Cisplatin: 1 mg/mL (1.5 mM); 5 mM in DMSO"
    # remove everything after the colon
    name <- gsub(":.*", "", name)

    # remove ,  ;  -  +  *  $  %  #  ^  _  as well as any spaces
    name <- gsub("[\\,\\;\\+\\*\\$\\%\\#\\^\\_\\s]", "", name, perl = TRUE)

    # remove hyphen 
    name <- gsub(" ", "_", name)

    # # remove substring of round brackets and contents
    name <- gsub("\\s*\\(.*\\)", "", name)

    # # remove substring of square brackets and contents
    name <- gsub("\\s*\\[.*\\]", "", name)

    # # remove substring of curly brackets and contents
    name <- gsub("\\s*\\{.*\\}", "", name)

    # # convert entire string to uppercase
    name <- toupper(name)

    # dealing with unicode characters 
    name <- gsub("Unicode", "", iconv(name, "LATIN1", "ASCII", "Unicode"), perl=TRUE)

    name
}

# Function to calculate the mean of a specific tag in the data
# Parameters:
#   - myDat: The input data table
#   - tag_name: The name of the tag to calculate the mean for
#   - mean_col_name: The name of the column to store the calculated mean (default: "tag_mean")
# Returns:
#   - A data table with the calculated mean for each SCAN_ID
calcTagMean <- function(myDat, tag_name, mean_col_name = "tag_mean") {
    suppressMessages(
        check_for_tag <- left_join(
            myDat %>% select_(~SCAN_ID) %>% distinct(),
            myDat %>% group_by_(~SCAN_ID) %>% filter_(~TAG == tag_name) %>% count()
        )
    )

    e1 <- simpleError(paste(
        "calcTagMean:",
        tag_name,
        "is not present for some or all of the SCAN_IDs in your data.",
        sep = " "
    ))
    if (any(is.na(check_for_tag$n))) {
        stop(e1)
    }

    tag_means <- myDat %>%
        group_by_(~SCAN_ID) %>%
        filter_(~TAG == tag_name) %>%
        summarise_(tag_mean = ~mean(INTENSITY, na.rm = T))

    e2 <- simpleError(paste(
        "calcTagMean:",
        tag_name,
        "has a mean of NaN for some or all of the SCAN_IDs in your data.",
        sep = " "
    ))
    if (any(is.nan(tag_means$tag_mean))) {
        stop(e2)
    }

    tag_means <- tag_means %>% rename_(.dots = stats::setNames("tag_mean", mean_col_name))

    return(tag_means)
}


#' Normalize Data
#'
#' This function normalizes the input data.
#'
#' @param myDat The input data to be normalized.
#' @param trim Logical value indicating whether to trim the data.
#' @param neg_control The negative control sample used for normalization.
#' 
#' @return The normalized data.
#'
#' @examples
#' myData <- read.csv("data.csv")
#' normalizedData <- normalizeData(myData, trim = TRUE, neg_control = 'NC-1')
#' 
normalizeData <- function(myDat, trim = TRUE, neg_control = 'NC-1',
                                                    pos_control = 'B') {
    
    # Calculate mean for negative control
    nc1 <- calcTagMean(myDat, tag_name = neg_control, mean_col_name = "NC")
    
    # Calculate mean for positive control
    pc1 <- calcTagMean(myDat, tag_name = pos_control, mean_col_name = "PC")
    
    # Filter and select relevant columns
    #  

    # This code performs several operations on the 'myDat' data frame:
    # 1. Filters rows where the 'TAG' column matches the pattern "(A|L|R)\\d+(-D\\d+)?-(S|C)" using the 'grepl' function.
    # 2. Selects specific columns from the filtered data frame using the 'select_' function.
    #    The selected columns include 'SCAN_ID', 'BARCODE', 'RESEARCH_PROJECT', 'DATE_CREATED', 'DRUGSET_ID',
    #    'CELL_LINE_NAME', columns ending with "CELL_ID", 'COSMIC_ID', 'POSITION', 'TAG', 'DRUG_ID', 'CONC', and 'INTENSITY'.
    # 3. Modifies the data frame using the 'mutate_' function.
    #    - Creates a new column 'lib_drug' by extracting the pattern "((L|R)\\d+)(-D\\d+)?-(S|C)" from the 'TAG' column.
    #      If the extracted value starts with 'A', it is replaced with NA.
    #    - Creates a new column 'anchor' by extracting the pattern "(A\\d+)-(S|C)" from the 'TAG' column.
    #      If the extracted value starts with 'L' or 'R', it is replaced with NA.
    #    - Creates a new column 'dose' by extracting the pattern "(A|L|R)\\d+-?(D\\d+)?-(S|C)" from the 'TAG' column.
    #    - Creates a new column 'treatment' by extracting the pattern "((A|L|R)\\d+)(-D\\d+)?-(S|C)" from the 'TAG' column.
    #    The modified columns are added to the data frame.
    # 4. Removes the 'TAG' column from the data frame using the 'select_' function.
    # The resulting data frame is stored in the 'normalized_data' variable.
    normalized_data <- myDat %>% 
        filter_(~ grepl("(A|L|R)\\d+(-D\\d+)?-(S|C)", TAG)) %>%
        select_(~SCAN_ID, ~BARCODE,  ~RESEARCH_PROJECT, ~DATE_CREATED, ~DRUGSET_ID,
                        ~CELL_LINE_NAME, ~ends_with("CELL_ID"),  ~COSMIC_ID, ~POSITION,
                        ~TAG, ~DRUG_ID, ~CONC, ~INTENSITY) %>%
        mutate_(lib_drug = ~sub("((L|R)\\d+)(-D\\d+)?-(S|C)", "\\1", TAG),
                        lib_drug = ~ifelse(grepl("^A.+", lib_drug), yes = NA, no = lib_drug),
                        anchor = ~sub("(A\\d+)-(S|C)", "\\1", TAG),
                        anchor = ~ifelse(grepl("^(L|R).+", anchor), yes = NA, no = anchor),
                        dose = ~sub("(A|L|R)\\d+-?(D\\d+)?-(S|C)", "\\2", TAG),
                        treatment = ~sub("((A|L|R)\\d+)(-D\\d+)?-(S|C)", "\\4", TAG)
        ) %>%
        select_(~-TAG)
    
    # Create libraries dataframe
    libraries <- normalized_data %>% 
        filter_(~!is.na(lib_drug)) %>% 
        select_(~-anchor, DRUG_ID_lib = ~DRUG_ID, CONC = ~CONC)
    
    # Create anchors dataframe
    anchors <- normalized_data %>% 
        filter_(~!is.na(anchor)) %>% 
        select_(~-lib_drug, ~-dose, DRUG_ID_anch = ~DRUG_ID, CONC_anch = ~CONC)
    
    # Join libraries and anchors dataframes
    if (nrow(anchors) > 0) {
        suppressMessages(normalized_data <- full_join(libraries, anchors))
    } else {
        normalized_data <- libraries
    }

    # Left join with negative control mean
    normalized_data <- left_join(normalized_data, nc1, by = "SCAN_ID")
    
    # Left join with positive control mean
    normalized_data <- left_join(normalized_data, pc1, by = "SCAN_ID")
    
    # Calculate normalized intensity
    normalized_data <- normalized_data %>%
        mutate_(normalized_intensity = ~((INTENSITY - PC) / (NC - PC)))

    # Trim normalized intensity if trim is TRUE
    if (trim) {
        normalized_data <- normalized_data %>%
            mutate_(normalized_intensity = 
                                ~(ifelse(normalized_intensity > 1, 1, normalized_intensity))) %>%
            mutate_(normalized_intensity =
                                ~(ifelse(normalized_intensity < 0, 0, normalized_intensity)))
    }
    
    # Add norm_neg_pos column
    normalized_data <- normalized_data %>% 
        mutate_(norm_neg_pos = ~paste(neg_control, pos_control, sep = "+"))

    # Add time_stamp and sw_version columns
    normalized_data <- normalized_data %>% 
        mutate_(time_stamp = ~ Sys.time()) %>%
        mutate_(sw_version = ~ set_package_used())

    return(normalized_data)
}