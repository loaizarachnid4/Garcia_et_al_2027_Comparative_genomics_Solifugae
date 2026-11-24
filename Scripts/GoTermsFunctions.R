# Load packages
library(httr)
library(jsonlite)


# Create a function to query AmiGO API for a single GO ID column
get_go_data_single <- function(go_id) {
  # QuickGO REST API URL
  url <- paste0("https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/", go_id)
  
  # Send GET request to the QuickGO API
  response <- GET(url)
  
  if (status_code(response) == 200) {
    # Parse the JSON response
    data <- fromJSON(content(response, "text", encoding = "UTF-8"))
    
    # Extract relevant information from the response
    id <- data$result$id
    name <- data$result$name
    definition <- data$result$definition
    aspect <- data$result$aspect
    
    # Return the data in a list format
    return(list(
      id= id,
      name = name,
      definition = definition,
      aspect = aspect
    ))
  } else {
    # Print a warning if the request failed and return NULL
    warning(paste("Failed to fetch data for GO term:", go_id, "with status code", status_code(response)))
    return(NULL)
  }
}


get_go_data <- function(go_id) {
  # Check if the go_id is blank or NULL
  if (is.null(go_id) || go_id == "") {
    warning("GO ID is blank or NULL. Skipping URL creation.")
    return(NULL)  # Skip the rest of the function and return NULL
  }
  
  # QuickGO REST API URL
  url <- paste0("https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/", go_id)
  
  # Send GET request to the QuickGO API
  response <- GET(url)
  
  if (status_code(response) == 200) {
    # Parse the JSON response
    data <- fromJSON(content(response, "text", encoding = "UTF-8"))
    
    # Extract relevant information from the response
    id <- data$result$id
    name <- data$result$name
    definition <- data$result$definition
    aspect <- data$result$aspect
    
    # Return the data in a list format
    return(list(
      id = id,
      name = name,
      definition = definition,
      aspect = aspect
    ))
  } else {
    # Print a warning if the request failed and return NULL
    warning(paste("Failed to fetch data for GO term:", go_id, "with status code", status_code(response)))
    return(NULL)
  }
}
