## Setup ----
library(dplyr)
library(tidyr)
library(gcplyr)

dup_rows <- function(mydf, mut_col, from_vals, to_vals) {
  #this function takes in a data.frame
  #finds all rows where df[, mut_col] == from_vals
  #deletes those rows
  #and adds rows that replicate the original but with to_vals in the
  # mut_col instead
  
  mydf <- as.data.frame(mydf)
  
  #Generate vector of rows to duplicate
  rowsrep <- ifelse(mydf[, mut_col] %in% from_vals,
                    length(to_vals), 1)
  idx <- rep(1:nrow(mydf), rowsrep)
  
  #Duplicate rows
  newdf <- mydf[idx, ]
  
  #replace from_vals with to_vals in duplicated rows
  newdf[which(idx %in% which(rowsrep > 1)), mut_col] <- to_vals
  
  return(newdf)
}

mynormalize <- function(x, a_input, b_input, a_output, b_output) {
  #normalize a vector of data, x, by scaling values such that 
  #a_input becomes a_output and b_input becomes b_output
  stopifnot(b_input > a_input, b_output > a_output)
  
  return((b_output-a_output)*(x-a_input)/(b_input - a_input) + a_output)
}

