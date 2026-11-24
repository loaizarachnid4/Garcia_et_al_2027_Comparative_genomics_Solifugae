library(biomaRt)
setwd("")
source("GoTermsFunctions.R")


####use this when only grepping for single files
##single with query searches matching fasta headers
head <- read.csv ("annotation.hogs2.absrel.csv")
#single file meta data after interpro
fatty <- read.csv ("Full_Metadata_cleaned_all_fatty.fasta.tsv.csv")

#Isolate branch names to search
targets <- head$Paralog

# use fixed = TRUE or escape them.
target_pattern <- paste(targets, collapse = "|")

# 4. Subset rows where V1 matches any of the targets
selected <- fatty[grepl(target_pattern, fatty$V1), ]

# 5. Write the individual 'selected' table to a file
write.table(selected, file = "ALL_fatty.txt", sep = "\t", 
            row.names = FALSE, col.names = TRUE, quote = FALSE)



#####this code subsets master interproscan list based on single column list with fasta headers###

##text file with branches
files <- list.files(full.names = TRUE, pattern = '*ps.txt')
file.list <- lapply(files, read.delim,header=F)
names(file.list) <- files


##tsv files from interproscan
tsvs <- list.files(full.names = TRUE, pattern = '*fasta.tsv')
tsv.list <- lapply(tsvs, read.delim, header =F)


#create empty list for subsetting
results_list <- list()

for (i in 1:length(file.list)) {
  x <- file.list[[i]]
  targets <- x$V2 #isolate branch names to search
  
  # Prepare search targets for grepl
  target_pattern <- paste(targets, collapse = "|")
  # Get the corresponding data from tsv.list
  current_tsv <- tsv.list[[i]]
  # Subset rows where V1 matches any of the targets
  selected <- current_tsv[grepl(target_pattern, current_tsv$V1), ]
  # Store in our results list using the original name
  results_list[[names(file.list)[i]]] <- selected
  # Define a unique filename for this specific file
  file_name <- paste0(names(file.list)[i], "_selected.txt")
  # Write the individual 'selected' table to a file
  write.table(selected, file = file_name, sep = "\t", 
              row.names = FALSE, col.names = TRUE, quote = FALSE)
}

# combine all into on dataframe
final_table <- do.call(rbind, results_list)
write.table(final_table, file = "all_positive_selected.txt", sep = "\t", 
            row.names = FALSE, col.names = TRUE, quote = FALSE)
  
#################################################################
##read files back in - this will only consider e-values with a numerical value
ip <- list.files(full.names = TRUE, pattern = '*.tsv')
ip.list <- lapply(ip, read.delim, header = F)

ip.list_cleaned <- lapply(seq_along(ip.list), function(i) {
  
  # Get the current data frame
  df <- ip.list[[i]]
  
  # this sorts by unique identifer and score taking smallest e value
  df$V9[df$V9 == "-"] <- "100"
  df$V9 <- as.numeric(df$V9)
  df <- df[order(df$V1, df$V9), ]
  df_unique <- df[!duplicated(df$V1), ]
  
  # Create a new filename (e.g., "original_name_cleaned.txt")
  out_name <- paste0("cleaned_", basename(ip[i]))
  
  write.table(df_unique, file = out_name, sep = "\t", 
              row.names = FALSE, quote = FALSE)
  
  # Still return the object so ip.list_cleaned is populated in R
  return(df_unique)
})

# no clean up
names(ip.list_cleaned) <- basename(ip)
#################GO annotations retrieval #################

#########################################################

# After doing some cleanup in excel, read files in again. sometimes multiple GO terms are retrieved
#so those must be seperated into different columns. Then the database has a (PANTHER) or (Interpro) match
#which must be deleted

clean.ip <- list.files(full.names = TRUE, pattern = "^cleaned.*")
ci.list <- lapply(clean.ip, read.delim, header = TRUE)

#iterate through all dataframes and change column names - must change starting column number based on how many go.terms columns
#this object will then be used in the next function
ci.list <- lapply(ci.list, function(df) {
  
  # find total number of go terms
  total_cols <- ncol(df)
  
  # If there are GO columns (starting at 14)
  if (total_cols >= 14) {
    # Generate names: GO.annotation.1, GO.annotation.2, etc.
    num_go_cols <- total_cols - 13
    go_names <- paste0("GO.annotation.", 1:num_go_cols)
    
    # Apply names to columns 14 through the end
    colnames(df)[14:total_cols] <- go_names
  }
  
  return(df)
})

## 

# iterate through each data frame in your ci.list
final_list <- lapply(seq_along(ci.list), function(idx) {
  
  gdor <- ci.list[[idx]]
  # Identify which columns are GO.annotation columns
  go_cols <- grep("GO.annotation", colnames(gdor), value = TRUE)
  
  # Extract and Fetch GO Data
  all_results <- list()
  
  for (col_name in go_cols) {
    # Extract the GO IDs for this specific column
    ids_vector <- gdor[[col_name]]
    
    # Apply your fetching function (get_go_data) from "GoTermsFunctions.R"
    go_results <- lapply(ids_vector, get_go_data)
    
    # Convert results to a clean data frame
    df <- do.call(rbind, lapply(go_results, function(x) {
      data.frame(
        id = ifelse(is.null(x$id), NA, x$id),
        name = ifelse(is.null(x$name), NA, x$name),
        definition = ifelse(is.null(x$definition), NA, x$definition$text),
        aspect = ifelse(is.null(x$aspect), NA, x$aspect),
        stringsAsFactors = FALSE
      )
    }))
    
    # Rename columns, e.g., GO.1_name, GO.1_aspect, etc.
    suffix <- gsub("GO.annotation.", "", col_name)
    colnames(df) <- paste0("GO.", suffix, "_", colnames(df))
    
    all_results[[col_name]] <- df
  }
  
  # Combine fetched data with the original file
  final_go_metadata <- do.call(cbind, all_results)
  full_table <- cbind(gdor, final_go_metadata)
  
  # Export the specific file
  # Uses the original filename from clean.ip to make a new CSV name
  output_filename <- paste0("Full_Metadata_", basename(clean.ip[idx]), ".csv")
  write.csv(full_table, file = output_filename, row.names = FALSE)
  
  return(full_table)
})


############################################################
#prepare interproscan output to only isolate low E values suggesting moderate to high confidence matches
gdor <-read.csv("Full_Metadata_cleaned_all_fatty.fasta.tsv.csv")
#convert e-score to numeric
gdor$V9 <- as.numeric(gdor$V9)

#keep those with low e-scores only
low.e <-gdor[gdor$E.value<=1,]

#remove those with NA
filt <- low.e[!is.na(low.e$E.value),]
go.terms <- filt[filt$GO.annotation !="-", ]

#use if low NA and blank filters were conducted externally
gdor <- gdor[gdor$GO.annotation.1 !="-", ]

##modify in excel to split multiple go terms into multiple columns
write.csv(go.terms, file="go.terms_hemero.csv", row.names = F)

#read in modified csv
gdor <- read.csv("go.terms_gdor.csv")


#######create vector of columns names depending on how many columns are present#####
#check amount of go term columns
colnames(go.terms)

#create column names for the number of columns (must do it manually in this code)
go.colms <- NULL
for (i in 1:9) {
  # Concatenate "GO.annotation" with the value of i depending on the amount of go terms hits
  go.colms[i] <- as.character(paste("GO.annotation", i, sep = "."))
}

#change columns names of the go terms
colnames(gdor)[15:23] <- go.colms


# Create an empty list to store the column names if needed
all.go <- list()

# Loop to iterate over the go.colms vector and extract columns from gdor
for (i in seq_along(go.colms)) {
  # Dynamically access each column of gdor using the go.colms vector element
  column_data <- gdor[[go.colms[i]]]  # Extract the column as a vector
  
  # Store the column data in all.go list
  all.go[[go.colms[i]]] <- column_data
}


