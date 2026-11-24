require(jsonlite)
library(ggtree)
library(ggplot2)
library(phytools)

setwd("")


###this first step was modified from https://github.com/stupornova33/parse_absREL_json###
###### we slightly modified the parse_significant_results.R script but maintained the same idea###
files <- list.files(full.names = TRUE, pattern = '*.json$')
#print(files)

#string to remove from file names - CHANGE ACCORDING TO YOUR FILE NAMES
string <- ""


#create an empty list for for an R table - this will be needed for downstream
pos.sel <- list()

#this isolates branches and nodes that were recovered as positively selected
#a table for each file that showed evidence for positive selection is exporded into working directory
for(file in files) {
  name_label <- gsub(string, "", basename(file)) #cut file name by string
  simplejson <- fromJSON(file)
  branchnames <- simplejson$`branch attributes`$`0`
  
  # Initialize an empty list or dataframe for this specific file's significant branches
  file_hits <- data.frame()
  
  for(idx in seq_along(branchnames)) {
    i <- branchnames[[idx]]
    branch <- names(branchnames)[idx]
    
    # Filter for significance
    if(i$`Corrected P-value` <= 0.05) {
      
      # Build local result row
      local_results <- data.frame(
        name = name_label, 
        branch = branch, 
        omega = as.numeric(i$`Rate Distributions`[1,1]),
        proportion = as.numeric(i$`Rate Distributions`[1,2]),
        omega2 = NA, # Initialize as NA
        proportion2 = NA,
        omega3 = NA,
        proportion3 = NA,
        corrected_pvalue=as.numeric(i$`Corrected P-value`),
        uncorrected_pvalue=as.numeric(i$`Uncorrected P-value`),
        LRT=as.numeric(i$`LRT`),
        stringsAsFactors = FALSE
      )
      
      # Fill based on classes
      if (i$`Rate classes` >= 2) {
        local_results$omega2   <- as.numeric(i$`Rate Distributions`[2,1])
        local_results$proportion2 <- as.numeric(i$`Rate Distributions`[2,2])
      } 
      if(i$`Rate classes` == 3) {
        local_results$omega3   <- as.numeric(i$`Rate Distributions`[3,1])
        local_results$percent3 <- as.numeric(i$`Rate Distributions`[3,2])
      }
      
      file_hits <- dplyr::bind_rows(file_hits, local_results)
    }
  }
  
  # Save individual files
  if (nrow(file_hits) > 0) {
    # Store in the list using the filename as the key
    pos.sel[[file]] <- file_hits
  }
}

####Combine the whole list into one giant table at the end and export to WD#####
final_table <- dplyr::bind_rows(pos.sel)
write.table(final_table, file = "all_positive_selected.txt", sep = "\t", 
            row.names = FALSE, col.names = TRUE, quote = FALSE)

###############################################################################
#### isolate trees from .json files that showed evidence for positive selection
###############################################################################
pos.trees <- list()

for (i in seq_along(files)) {
  # current_file is the full filename string
  current_file <- files[i]
  
  # Check if any name in your table is a partial match inside the filename
  # We use any() in case multiple patterns match, though usually it's 1:1
  match_idx <- which(sapply(final_table$name, function(x) grepl(x, current_file)))
  
  if (length(match_idx) > 0) {
    # If a match is found, extract the tree
    pos.trees[[i]] <- extract_newick(current_file)
    
    # Optional: Rename the list element for easier tracking
    names(pos.trees)[i] <- final_table$name[match_idx[1]]
  }
}
#final clean up step - necessary to remove empty and NA inputs
pos.trees <- pos.trees[!sapply(pos.trees, is.null)]

###############################################################################
#########isolate branches and nodes that are under positive selection 
#########from final_table result above to create a list of vectors######
###############################################################################
sel.branch <- list()

# Iterate through the names of your trees
for (tree_name in names(pos.trees)) {
  
  # Find the row index in final_table where the name matches
  match_row <- which(final_table$name == tree_name)
  
  # If a match is found, extract the branch vector
  if (length(match_row) > 0) {
    # Store the specific branch value(s) into your new list
    sel.branch[[tree_name]] <- final_table$branch[match_row]
  }
}

###############################################################################
#########plot trees: branches and nodes are separate###########
####this creates a combined pdf plot for all branches that showed signatures of positive selection##
###############################################################################
pdf("All_Branch_Trees.pdf", width = 12, height = 10)

for (tree_id in names(pos.trees)) {
  
  # Get the original vector
  current_targets <- sel.branch[[tree_id]]
  
  # Filter out items containing "Node" - this is in the next step
  # This keeps only the actual species names/labels
  clean_targets <- grep("Node", current_targets, value = TRUE, invert = TRUE, ignore.case = TRUE)
  
  # Check if the length is 0, it means the list only had "Node" items (or was empty)
  if (length(clean_targets) == 0) {
    message(paste("Skipping", tree_id, " - No tip labels found."))
    next  # Jump to the next tree_id in names(pos.trees)
  }
  
  # Highlight all branches in all trees what have evidence of positive selection
  p <- ggtree(pos.trees[[tree_id]]) + 
    geom_tree(aes(color = label %in% clean_targets), linewidth = 1.2) + #change line width
    scale_color_manual(values = c("black", "red")) +
    geom_tiplab() + hexpand(0.3, direction = 1) + labs(color = "Positive Selection") +
    labs(title = tree_id)
  
  print(p)
}
# Close and save the file
dev.off()


###############################################################################
####this creates a combined pdf plot for all nodes that showed signatures of positive selection##
###This also creates three separate text files with fasta headers####
###############################################################################

#initialize pdf generation
pdf("All_Node_Trees.pdf", width = 14, height = 10)

###MUST CHANGE THESE SEARCH STRINGS
#these will be the strings that will used for both naming the files and for searching for branch names
sols <- c("Hemer", "Gdor", "Ppal")

#initialize node searches with evidence of positive selection
for (tree_id in names(pos.trees)) {
  current_tree <- pos.trees[[tree_id]]
  current_targets <- sel.branch[[tree_id]]
  
  # Identify nodes that were showed evidence of positive selection
  clean_targets <- grep("Node", current_targets, value = TRUE, invert = FALSE, ignore.case = TRUE)
  node.idx <- match(clean_targets, current_tree$node.label)
  node.idx <- node.idx[!is.na(node.idx)]
  
  if (length(node.idx) == 0) next
  
  # Convert to absolute tree node IDs
  absolute_node_ids <- node.idx + length(current_tree$tip.label)
  
  # Tree generation to highlight nodes under positive selection - all trees get exported to pdf
  colored_tree <- paintSubTree(current_tree, node = length(current_tree$tip.label) + 1, 
                               state = "1", anc.state = "1")
  
  for (target_node in absolute_node_ids) {
    colored_tree <- paintSubTree(colored_tree, node = target_node, state = "2", stem = TRUE)
  }
  
  cols <- setNames(c("black", "red"), c("1", "2"))
  plotSimmap(colored_tree, colors = cols, lwd = 3)
  title(main = tree_id, line = -1)
  
  # now this extracts the tips for each of those nodes - the tips here correspond to fasta headers that will
  # be used for later extraction
  node.tips <- list()
  for (i in seq_along(absolute_node_ids)) {
    # Extract clade and grab tips directly
    tmp_clade <- extract.clade(current_tree, node = absolute_node_ids[i])
    node.tips[[i]] <- tmp_clade$tip.label
  }
  
  # unlist to be able to grep
  all_found_tips <- unlist(node.tips)
  
  # Now pull tip labels from branches and nodes using grep, then combine into a single dataframe
  for (s in sols) {
    # 1. Filter: Find ONLY matches for the species string 's'
    res1 <- grep(s, final_table$branch, value = TRUE, invert = FALSE, ignore.case = TRUE)
    res2 <- grep(s, all_found_tips, value = TRUE, invert = FALSE, ignore.case = TRUE)
    
    # Combine
    combined_headers <- c(res1, res2)
    
    # create file name based on the sols variable
    file_name <- paste0(s, "_headers.txt")
    
    # append into a single dataframe for each taxon
    cat(combined_headers, file = file_name, sep = "\n", append = TRUE)
  }
}

# close pdf call
dev.off()
