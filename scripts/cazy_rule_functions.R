library(tidyverse)
library(lubridate)
#check for id in function_id column - T if id is in cell, F if id not in cell
checkID <- function(row, function_df, id){
  cell <- function_df %>%
    slice(row) %>%
    pull(function_ids)
  
  return(str_detect(cell, id))
}


#Make function/bin df
makeFunctBinDf <- function(columns, bin_row){
  df <- data.frame(row.names = bin_row)
  for(i in columns){
    c_add <- data.frame(i = "x")
    df <- cbind(df, c_add)
  }
  colnames(df) <- columns
  return(df)
}


#Make function boolean for bin
functionBool <- function(sub_names, column, function_df, dbcan){
  ids_bin <- dbcan[dbcan[,column] > 0, c("target_name", column)]
  df_out <- makeFunctBinDf(columns = sub_names, bin_row = column)
  
  for(i in sub_names){
    funct <- function_df %>%
      filter(str_detect(function_name, i))
    if(nrow(funct) == 1){
      detected <- 0
      for(j in ids_bin$target_name){
        detected <- detected + checkID(row = 1, function_df = funct, id=j)
      }
      df_out[column, i] <- detected > 0
    }
    if(nrow(funct) == 2){
      detected1 <- 0
      detected2 <- 0
      for(j in ids_bin$target_name){
        detected1 <- detected1 + checkID(row = 1, function_df = funct, id=j)
        detected2 <- detected2 + checkID(row = 2, function_df = funct, id=j)
      }
      df_out[column, i] <- detected1 > 0 & detected2 > 0
    }
  }
  
  return(df_out)
}


#create df with rownames = bins, colnames = function, TRUE/FALSE if present
buildProduct <- function(sub_names, function_df, dbcan){
  
  columns <- colnames(dbcan)[-1]
  df <- functionBool(sub_names,
                     column = columns[1],
                     function_df,
                     dbcan)
  #progressbar
  n_iter <- length(columns) - 1
  pb <- txtProgressBar(min = 0,
                       max = n_iter,
                       style = 3,
                       width = n_iter, 
                       char = "=") 
  
  init <- numeric(n_iter)
  end <- numeric(n_iter)
  
  for(i in  2:length(columns)){
    init[i] <- Sys.time()
    
    #code
    df <- rbind(df, functionBool(sub_names,
                                 column = columns[i],
                                 function_df,
                                 dbcan))
    
    end[i] <- Sys.time()
    
    #progress bar
    setTxtProgressBar(pb, i)
    time <- round(seconds_to_period(sum(end - init)), 0)
    est <- n_iter * (mean(end[end != 0] - init[init != 0])) - time
    remaining <- round(seconds_to_period(est), 0)
    cat(paste(" // Execution time:", time,
              " // Estimated time remaining:", remaining), "")
  }
  return(df)
}
