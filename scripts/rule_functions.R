library(tidyverse)

#functions to parse through rules

#build checkRule function from rule_df
makecheckRule <- function(rule_df) {
  switch_list <- list()
  switch_list <- append(switch_list, rule_df[["child"]])
  names(switch_list) <- rule_df[["parent"]]
  function(funct) {
    pluck(switch_list, funct)
  }
}

#evaluates string from rule structure and outputs rule list, used in evaluateGenome()
evaluateRule <- function(rule_df, rule_out){
  if(is.na(rule_out)){
    return(NA)
  }else{
    if(!hasIds(rule_out)){ 
      ids <- pullIds(rule_df, rule_out)
      rule_list <- map(ids, checkRule)
      names(rule_list) <- ids
    }else{
      rule_list <- as.list(rule_out)
      names(rule_list) <- rule_df[rule_df$child == rule_out, "parent"]
    }
    
    #handle ETC50
    for(rule in names(rule_list)){
      if(hasIds(rule_list[rule])){
        next
      }else{
        rule_list[rule] <- list(evaluateRule(rule_df, rule_list[rule]))
      }
    }
    
    #add original rule to list
    original_names <- names(rule_list)
    rule_list[length(rule_list)+1] <- rule_out
    names(rule_list) <- c(original_names, "org")
    
    return(rule_list)
  }
}

#replaces all instances of an id in a list with repl: used in buildLogicList()
replaceIdInList <- function(list, id, repl){
  list_out <- list
  if(listOfLists(list)){
    idx <- str_detect(map_chr(list, class), "list")
    for(elem in idx){
      if(elem){
        new_nested <- list_out[[which(idx)]]
        for(elem1 in seq_along(new_nested)){
          new_nested[[elem1]] <- gsub(pattern = id, 
                                      x = new_nested[[elem1]], 
                                      replacement = repl)
        }
        list_out[[which(idx)]] <- new_nested
      }else{
        for(elem1 in which(!idx)){
          list_out[[elem1]] <- gsub(pattern = id, 
                                    x = list_out[[elem1]], 
                                    replacement = repl)
        }
      }
    }
  }else{
    for(elem1 in seq_along(list_out)){
      list_out[[elem1]] <- gsub(pattern = id, 
                                x = list_out[[elem1]], 
                                replacement = repl)
    }
  }
  return(list_out)
}


#check if list has another list as an element: used in buildLogicList()
listOfLists <- function(list){
  return(sum(str_detect(map_chr(list, class), "list")) > 0)
}

#match EC number format EC-.-.-.- : used in hasIds(), pullIds(), buildLogicList()
matchEC <- function(string){
  x <- str_extract_all(string, 
                       "EC[\\d-]{1,2}\\.[\\d-]{1,2}\\.[\\d-]{1,2}\\.[\\d-]{1,2}")
  return(unlist(x))
}


#match KO number format K----- : used in hasIds(), pullIds(), buildLogicList()
matchKO <- function(string){
  x <- str_extract_all(string, 
                       "K[\\d]{5}")
  return(unlist(x))
}

#check if string has KO or EC ids: used in evaluateRule()
hasIds <- function(string){
  return(length(matchEC(string)) + length(matchKO(string)) > 0)
}

#pull individual Ids: used in buildLogicList(), evaluateRule()
pullIds <- function(rule_df, string){
  out_vect <- c()
  for(i in rule_df[["parent"]]){
    if(str_detect(string, i)){
      end_loc <- str_locate(string, i)[2]
      if(substr(string, end_loc+1, end_loc+1) %in% c("]", "&", "|", ",") |
         end_loc+1 > nchar(string)){
        out_vect <- c(out_vect, i)
      }else{
        next
      }
    }
  }
  if(is.null(out_vect)){
    out_vect <- c(out_vect, matchKO(string), matchEC(string))
    return(out_vect)
  }else{
    return(out_vect)
  }
}




#checks rule list against counts dataframe and provides T/F list for given genome: used in evaluateGenome()
buildLogicList <- function(dataframe, rule_df, rule_list, genome){
  
  #check for nestedness
  if(listOfLists(rule_list)){ 
    nested_index <- str_detect(map_chr(rule_list, class), "list")
    nested_list <- pluck(rule_list,which(nested_index))
    rule_list[[which(nested_index)]] <- buildLogicList(dataframe,
                                                       rule_df,
                                                       nested_list, 
                                                       genome)
  }
  
  #get unique IDs
  ids <- c()
  for(element in unlist(rule_list)){
    ids <- c(ids, pullIds(rule_df, element))
  }
  
  names(ids) <- NULL
  ids <- unique(c(matchKO(ids), matchEC(ids)))
  
  df <- dataframe %>%
    filter(bin.id == genome)
  
  for(id in ids){
    if(!(id %in% df[,"gene_id"])){
      rule_list <- replaceIdInList(rule_list, id, repl = "F")
    }else if(df[which(df[,"gene_id"] == id),"counts"] > 0){
      rule_list <- replaceIdInList(rule_list, id, repl = "T")
    }else if(df[which(df[,"gene_id"] == id),"counts"] == 0){
      rule_list <- replaceIdInList(rule_list, id, repl = "F")
    }
  }
  
  return(rule_list)
}



#build index of evaluation priority (run on individual elements of rule_list): used in evaluateFunction() and handleSublist()
buildIndex <- function(rule_string){
  index <- c()
  for(i in 1:nchar(rule_string)){
    if(str_sub(rule_string, i, i) == "["){
      index <- c(index, i)
    }else if(str_sub(rule_string, i, i) == "]"){
      index <- c(index, -i)
    }
  }
  return(index)
}

#replaces DBL instances in org string: used in buildPriorityList()
replaceDBLs <- function(priority_list, dbl_names){
  out_list <- priority_list
  for(i in dbl_names){
    dbl_repl <- str_extract(string = out_list[["org"]], pattern = "\\[.*\\]" )
    out_list[[i]] <- dbl_repl
    out_list[[i]] <- str_remove_all(out_list[[i]], "\\[")
    out_list[[i]] <- str_remove_all(out_list[[i]], "\\]")
    
    out_list[["org"]] <- str_replace(out_list[["org"]], 
                                     out_list[[i]], i)
    
    out_list[["org"]] <-str_remove_all(out_list[["org"]], "\\[")
    out_list[["org"]] <-str_remove_all(out_list[["org"]], "\\]")
  }
  return(out_list)
}

#build list for rule evaluation priority (run on individual elements of rule_list): used in evaluateFunction(), handleSublist()
buildPriorityList <- function(index, rule_string){
  priority_list <- list()
  dbl_idx_strt <- c()
  dbl_idx_end <- c()
  org <- rule_string
  
  if(is.null(index)){
    priority_list[length(priority_list)+1] <- rule_string
    names(priority_list) <- "org"
    return(priority_list)
  }else{
    let_cnt <- 1
    repl_idx <- abs(index)
    for(i in 2:length(index)){
      if(index[i] > 0 & index[i-1] > 0){
        dbl_idx_strt <- c(dbl_idx_strt, index[i-1])
      }else if(index[i] < 0 & index[i-1] > 0){
        x <- unlist(str_sub(rule_string, index[i-1], -index[i]))
        priority_list[length(priority_list)+1] <- str_remove_all(x, "\\[|\\]")
        str_sub(string = org, 
                start = repl_idx[i-1], 
                end = repl_idx[i]) <- letters[let_cnt]
        let_cnt <- let_cnt + 1
        repl_idx <- repl_idx - nchar(str_sub(string = org, 
                                             start = repl_idx[i-1], 
                                             end = repl_idx[i])) + 1 
      }else if(index[i] < 0 & index[i-1] < 0){
        dbl_idx_end <- c(dbl_idx_end, -index[i-1])
      }
    }
    names(priority_list) <- letters[1:length(priority_list)]
    first_names <- names(priority_list)
    
    dbl_names <- c()
    if(length(dbl_idx_strt) != 0 & length(dbl_idx_end) != 0){
      for(i in 1:length(dbl_idx_strt)){
        priority_list[length(priority_list)+1] <- str_sub(rule_string, dbl_idx_strt[i], dbl_idx_end[i])
        dbl_names <- c(dbl_names, paste0("DBL", as.character(i)))
      }
    }
    
    
    priority_list[length(priority_list)+1] <- org
    names(priority_list) <- c(first_names, dbl_names, "org")
    priority_list <- replaceDBLs(priority_list, dbl_names)
    
    return(priority_list)
  }
}


#collapses priority logic list: used in evaluateFunction(), handleSublist()
collapsePriority <- function(logic_list){
  collapsed_list <- logic_list
  elems_rm <- c()
  for(elem in names(collapsed_list)){ #itter over list names
    collapsed_list[[elem]] <- evalLogic(collapsed_list[[elem]]) #evaluate logic
    if(is.logical(collapsed_list[[elem]])){ #if list element is logical continue
      for(elem_idx in seq_along(collapsed_list)){ #itter over list elements after confirming logic
        if(str_detect(collapsed_list[[elem_idx]], elem)){ #T if list name is present in list object string
          repl_char <- ifelse(collapsed_list[[elem]], "T", "F") #create replacement character depending on logic value
          collapsed_list[[elem_idx]] <- str_replace(collapsed_list[[elem_idx]], 
                                                    elem, 
                                                    repl_char) #replace list name in list element string with T or F
          elems_rm <- c(elems_rm, elem) #add list element to be removed
        }
      }
    }else{
      next
    }
  }  
  for(i in elems_rm){
    collapsed_list[[i]] <- NULL
  }
  return(collapsed_list)
}



#handles logic vector with no spec.chrs: used in collapseVect()
collapseLogic <- function(string_vect, lgc_type){
  lgc <- string_vect
  while(length(lgc) > 1){
    if(lgc_type == "&"){
      new_val <- as.character(as.logical(lgc[1]) & as.logical(lgc[2]))
      lgc <- c(new_val, lgc[-c(1,2)])
    }else if(lgc_type == "|"){
      new_val <- as.character(as.logical(lgc[1]) | as.logical(lgc[2]))
      lgc <- c(new_val, lgc[-c(1,2)])
    }
  }
  return(lgc)
}

#collapse and logic vector to string: used in evalLogic()
collapseVect <- function(string_vect, lgc_type){
  alt_lgc <- ifelse(lgc_type == "&", "|", "&")
  num_elems <- length(string_vect)
  string_collapsed <- "" 
  if(length(which(nchar(string_vect) == 1)) == num_elems){
    string_collapsed <- collapseLogic(string_vect, lgc_type)
  }else{
    for(elem in 2:length(string_vect)){
      if(nchar(string_vect[elem-1]) == 1 & string_vect[elem-1] != alt_lgc 
         & nchar(string_vect[elem]) == 1 & string_vect[elem] != alt_lgc){
        if(lgc_type == "&"){
          add_str <- ifelse(as.logical(string_vect[elem-1]) & as.logical(string_vect[elem]),"T","F")
          string_collapsed <- paste0(string_collapsed, add_str)
        }else if(lgc_type == "|"){
          add_str <- ifelse(as.logical(string_vect[elem-1]) | as.logical(string_vect[elem]),"T","F")
          string_collapsed <- paste0(string_collapsed, add_str)
        }
      }else if(nchar(string_vect[elem-1]) == 1 & string_vect[elem-1] != alt_lgc 
               & nchar(string_vect[elem]) == 1 & string_vect[elem] == alt_lgc){
        string_collapsed <- paste0(string_collapsed, string_vect[elem])
      }else if(nchar(string_vect[elem-1]) == 1 & string_vect[elem-1] == alt_lgc 
               & nchar(string_vect[elem]) == 1 & string_vect[elem] != alt_lgc){
        next
      }else if(nchar(string_vect[elem-1]) > 1){
        string_collapsed <- paste0(string_collapsed, string_vect[elem-1])
      }else if(nchar(string_vect[elem]) == 1){
        next
      }else if(nchar(string_vect[elem]) > 1 & elem == length(string_vect)){
        string_collapsed <- paste0(string_collapsed, string_vect[elem])
      }
    }
  }
  return(string_collapsed)
}


#step through string and return characters adjacent to &: used in evalLogic()
stepThroughAnd <- function(string){
  if(str_detect(string, "\\&")){
    out_vect <- c()
    split_vect <- unlist(strsplit(string, "&", fixed = TRUE))
    num_elems <- length(split_vect)
    if(length(which(nchar(split_vect) == 1)) == num_elems){
      out_vect <- split_vect
    }else{
      for(i in 1:num_elems){
        if(i == 1 & nchar(split_vect[i]) > 1){
          out_vect <- c(substr(split_vect[i], 1, nchar(split_vect[i]) - 1),
                        substr(split_vect[i], nchar(split_vect[i]), nchar(split_vect[i])))
        }else if(i == 1 & nchar(split_vect[i]) == 1){
          out_vect <- c(split_vect[i])
        }else if(i > 1 & i < num_elems & nchar(split_vect[i]) > 1){
          out_vect <- c(out_vect,
                        substr(split_vect[i],1,1),
                        substr(split_vect[i],2,nchar(split_vect[i])-1),
                        substr(split_vect[i], nchar(split_vect[i]), nchar(split_vect[i])))
        }else if(i > 1 & i < num_elems & nchar(split_vect[i]) == 1){
          out_vect <- c(out_vect, split_vect[i])
        }else if(i == num_elems & nchar(split_vect[i]) > 1){
          out_vect <- c(out_vect,
                        substr(split_vect[i], 1, 1),
                        substr(split_vect[i], 2, nchar(split_vect[i])))
        }else if(i == num_elems & nchar(split_vect[i]) == 1){
          out_vect <- c(out_vect, split_vect[i])
        }
      }
    }
  }else{
    out_vect <- c(string)
  }
  return(out_vect)
}

#step through string and return characters adjacent to |: used in evalLogic()
stepThroughOr <- function(string){
  if(str_detect(string, "\\|")){
    out_vect <- c()
    split_vect <- unlist(strsplit(string, "|", fixed = TRUE))
    num_elems <- length(split_vect)
    for(i in 1:num_elems){
      if(i == 1 & nchar(split_vect[i]) > 1){
        out_vect <- c(substr(split_vect[i], 1, nchar(split_vect[i]) - 1),
                      substr(split_vect[i], nchar(split_vect[i]), nchar(split_vect[i])))
      }else if(i == 1 & nchar(split_vect[i]) == 1){
        out_vect <- c(split_vect[i])
      }else if(i > 1 & i < num_elems & nchar(split_vect[i]) > 1){
        out_vect <- c(out_vect,
                      substr(split_vect[i],1,1),
                      substr(split_vect[i],2,nchar(split_vect[i])-1),
                      substr(split_vect[i], nchar(split_vect[i]), nchar(split_vect[i])))
      }else if(i > 1 & i < num_elems & nchar(split_vect[i]) == 1){
        out_vect <- c(out_vect, split_vect[i])
      }else if(i == num_elems & nchar(split_vect[i]) > 1){
        out_vect <- c(out_vect,
                      substr(split_vect[i], 1, 1),
                      substr(split_vect[i], 2, nchar(split_vect[i])))
      }else if(i == num_elems & nchar(split_vect[i]) == 1){
        out_vect <- c(out_vect, split_vect[i])
      }
    }
  }else{
    out_vect <- c(string)
  }
  return(out_vect)
}



#evaluates | and & logic: used in collapsePriority(), handleSublist()
evalLogic <- function(logic_string){
  #check if pure logic string
  check_str <- logic_string
  for(i in c("&", "\\|", ",", "F", "T")){
    check_str <- str_remove_all(check_str, i)
  }
  
  if(check_str != ""){
    logic_string_out <- logic_string
  }else{
    #check for and logic (&) first
    if(str_detect(logic_string, "\\&")){ 
      x_and <- stepThroughAnd(logic_string) #splits and expressions, returns chr vect of length 1 if no & present
      and_string_collapsed <- collapseVect(x_and, "&") #collapses chr vect to string evaluating & logic
      #check and_string_collapsed for or logic (|)
      if(str_detect(and_string_collapsed, "\\|")){
        x_and_or <- stepThroughOr(and_string_collapsed) #splits or expressions, returns chr vect of length 1 if no | present
        logic_string_out <- collapseVect(x_and_or, "|") #collapses chr vect to string evaluating & logic
      }else{
        logic_string_out <- and_string_collapsed
      }
    }else if(str_detect(logic_string, "\\|")){
      x_or <- stepThroughOr(logic_string) #splits or expressions, returns chr vect of length 1 if no | present
      logic_string_out <- collapseVect(x_or, "|") #collapses chr vect to string evaluating & logic
    }else{
      logic_string_out <- logic_string
    }
  }
  
  #check for switch to logic type
  logic_string_out <- ifelse(is.na(as.logical(logic_string_out)), 
                             logic_string_out,
                             as.logical(logic_string_out))
  return(logic_string_out)
}



#returns percent true over percent false in comma separated logic string: used in evaluateFunction()
commaPercent <- function(string){
  vect <- unlist(strsplit(string, split = ","))
  log_vect <- as.logical(vect)
  trues <- length(which(log_vect))
  falses <- length(which(!log_vect))
  
  percT <- round(trues/length(log_vect)*100, 0)
  percF <- round(falses/length(log_vect)*100, 0)
  return(paste0(percT, '/', percF))
}

#returns number of T in comma separated logic string: used in evaluateFunction()
commaCount <- function(string){
  vect <- unlist(strsplit(string, split = ","))
  log_vect <- as.logical(vect)
  trues <- length(which(log_vect))
  
  return(trues)
}

#function will search org for "not" and switch logic of following letter: used in evaluateFunction()
switchNot <- function(string){
  loc_vect <- unlist(str_locate_all(string, "not"))
  ends <- loc_vect[((length(loc_vect)/2)+1):length(loc_vect)]
  
  #switch logic
  for(i in ends){
    if(substr(string, i+1, i+1) == "T"){
      substr(string, i+1, i+1) <- "F"
    }else if(substr(string, i+1, i+1) == "F"){
      substr(string, i+1, i+1) <- "T"
    }
  }
  
  #remove "not"
  string <- str_remove_all(string, "not")  
  return(string)
}

#function will change string with "TRUE" or "FALSE" to one with "T" or "F": used in handleSublist() and evaluateFunction()
truefalseToTF <- function(string){
  if(str_detect(string, "TRUE")){
    string <- str_replace_all(string, "TRUE", "T")
  }
  
  if(str_detect(string, "FALSE")){
    string <- str_replace_all(string, "FALSE", "F")
  }
  
  if(!str_detect(string, "TRUE|FALSE")){
    return(string)
  }
  return(string)
}

#match if 50% marker is present: used in handleSublist() and evaluateFunction()
check50 <- function(string){
  if(str_detect(string, "percent50")){
    loc_50 <- unlist(str_match_all(string, "percent50(.+?)[\\|,\\]]|percent50(.+?)$"))
    loc_50_1 <- loc_50[!str_detect(loc_50, "percent50")]
    return(loc_50_1[!is.na(loc_50_1)])
  }else{
    return(NA)
  }
}

#match if 1of  marker is present: used in evaluateFunction()
checkOne <- function(string){
  if(str_detect(string, "1of")){
    loc_1 <- unlist(str_match_all(string, "1of(.+?)[\\|,\\]]|1of(.*?)$"))
    loc_1_1 <- loc_1[!str_detect(loc_1, "1of")]
    return(loc_1_1[!is.na(loc_1_1)])
  }else{
    return(NA)
  }
}

#match if 2of  marker is present: used in evaluateFunction()
checkTwo <- function(string){
  if(str_detect(string, "2of")){
    loc_2 <- unlist(str_match_all(string, "2of(.+?)[\\|,\\]]|2of(.*?)$"))
    loc_2_1 <- loc_2[!str_detect(loc_2, "2of")]
    return(loc_2_1[!is.na(loc_2_1)])
  }else{
    return(NA)
  }
}

#match if true value is >= provided percent value: used in handleSublist() and evaluateFunction()
checkTruePerc <- function(string, perc_check){
  perc_t <- as.double(str_extract(string, "[^/]+"))
  return(perc_t >= perc_check)
} 

#function checks for a sublist and resolves it to T or F: used in evaluateFunction()
handleSublist <- function(logic_list){
  
  if(sum(map_lgl(logic_list, is.list)) > 0){
    sub_list <- logic_list[[which(map_lgl(logic_list, is.list))]] #save sub list as separate list
    sub_list_idx <- which(map_lgl(logic_list, is.list)) #save sub list position in main list
    first_elems <- logic_list[seq(1,sub_list_idx-1)] #save original list element index < sublist
    last_elems <- logic_list[seq(sub_list_idx+1, length(logic_list))] #save original list element index > sublist
    out_list_names <- names(logic_list) #save list names
    
    #evaluate sub list logic
    for(element in seq_along(sub_list)){
      if(names(sub_list[element]) == "org"){
        next
      }
      elem_index <- buildIndex(sub_list[[element]])
      prir_list <-buildPriorityList(elem_index, sub_list[[element]])
      sub_list[element] <- collapsePriority(prir_list)
    }
    
    #handle percent50
    if(!(NA %in% check50(sub_list[["org"]]))){
      out <- str_remove_all(sub_list[["org"]], "percent50")
      for(name in check50(sub_list[["org"]])){
        sub_list[name] <- checkTruePerc(commaPercent(sub_list[[name]]), 50)
        out <- str_replace_all(out, name, as.character(sub_list[[name]]))
        out <- truefalseToTF(out)
      }
      sub_list <- evalLogic(out)
    }
    
    #rebuild main list
    out_list <- map(c(first_elems, sub_list, last_elems), unlist)
    names(out_list) <- out_list_names
    
    out_list <- collapsePriority(out_list)
  }else{
    out_list <- logic_list
  }
  return(out_list)
}

#evaluates logic_list returning list of T or F
evaluateFunction <- function(logic_list){
  out_list <- handleSublist(logic_list) #resolve sublists
  
  if(length(out_list) != 1){
    #parse logic
    for(element in names(out_list)[which(names(out_list) != "org")]){
      elem_index <- buildIndex(out_list[[element]])
      prir_list <- buildPriorityList(elem_index, out_list[[element]])
      if(is.logical(unlist(collapsePriority(prir_list)))){
        out_list[element] <- collapsePriority(prir_list)
      }else{
        
        #handle percent50
        if(!(NA %in% check50(out_list[["org"]]))){
          out <- out_list[["org"]]
          for(name in check50(out_list[["org"]])){
            if(name == element){
              out_list[element] <- checkTruePerc(commaPercent(unlist(prir_list)), 50)
              out <- str_replace_all(out, str_glue("percent50{name}"), as.character(out_list[[element]]))
              out_list[["org"]] <- truefalseToTF(out)
            }
          }
        }
        
        #handle 2of
        if(!(NA %in% checkTwo(out_list[["org"]]))){
          twoOf <- TRUE
          out <- out_list[["org"]]
          for(name in checkTwo(out_list[["org"]])){
            if(name == element){
              working_logic <- ifelse(commaCount(unlist(prir_list)) >= 2, TRUE, FALSE)
              out <- str_replace_all(out, str_glue("2of{name}"), as.character(working_logic))
              out_list["org"] <- truefalseToTF(out)
            }
          }
        }else{
          twoOf <- FALSE
        }
        
        #handle 1of
        if(!(NA %in% checkOne(out_list[["org"]]))){
          oneOf <- TRUE
          out <- out_list[["org"]]
          for(name in checkOne(out_list[["org"]])){
            if(name == element){
              working_logic <- ifelse(commaCount(unlist(prir_list)) >= 1, TRUE, FALSE)
              out <- str_replace_all(out, str_glue("1of{name}"), as.character(working_logic))
              out_list["org"] <- truefalseToTF(out)
            }
          }
        }else{
          oneOf <- FALSE
        }
        
        if(oneOf | twoOf){
          out_list[element] <- working_logic
        }
      }
      out_list <- collapsePriority(out_list)
    }
    
    #clean up collapsed list elements
    if(length(out_list) > 1){
      remove <- c()
      for(i in which(names(out_list) != "org")){
        if(is.logical(out_list[[i]])){
          remove <- c(remove, i)
        }
      }
      out_list <- out_list[-remove]
    }
  }
  
  #check for not
  if(str_detect(out_list[["org"]], "not")){
    out_list[["org"]] <- switchNot(out_list[["org"]])
  }
  
  #collapse org
  if(!(is.logical(out_list[["org"]]))){
    elem_index <- buildIndex(out_list[["org"]])
    prir_list <- buildPriorityList(elem_index, out_list[["org"]])
    out_list["org"] <- collapsePriority(prir_list)
  }
  return(out_list)
}


#calls other functions to evaluate individual genome for presence of function: used in evaluateCountsDf()
evaluateGenome <- function(counts_df, rule_df, genome, omit = NA){
  out_list <- list()
  for(fun in unique(rule_df[['name']])){
    funct_name <- fun
    
    rule_list <- evaluateRule(rule_df, checkRule(fun))
    
    logic_list <- buildLogicList(counts_df, rule_df, rule_list, genome)
    
    fun_list <- evaluateFunction(logic_list)
    
    names(fun_list) <- funct_name
    
    out_list <- append(out_list, fun_list)
  }
  
  fun_names <- names(out_list)
  out_list <- append(out_list, genome)
  names(out_list) <- c(fun_names, "genome")
  
  if(!is.na(omit)){
    out_list[[omit]] <- NULL
  }
  return(out_list)
}


#evaluates function presence for all genomes in counts_df and all functions in rule_df
evaluateCountsDf <- function(counts_df, rule_df, omit = NA){
  out_df <- data.frame()
  
  #progressbar
  n_iter <- length(unique(counts_df[['bin.id']]))
  pb <- txtProgressBar(min = 0,
                       max = n_iter,
                       style = 3,
                       width = n_iter, 
                       char = "=") 
  
  init <- numeric(n_iter)
  end <- numeric(n_iter)
  
  i <- 1
  for(gen in unique(counts_df[['bin.id']])){
    init[i] <- Sys.time()
    
    #code
    add_row <- evaluateGenome(counts_df, rule_df, gen, omit = omit)
    out_df <- rbind(out_df, as.data.frame(add_row))
    
    end[i] <- Sys.time()
    
    #progress bar
    setTxtProgressBar(pb, i)
    time <- round(lubridate::seconds_to_period(sum(end - init)), 0)
    est <- n_iter * (mean(end[end != 0] - init[init != 0])) - time
    remaining <- round(lubridate::seconds_to_period(est), 0)
    cat(paste(" // Execution time:", time,
              " // Estimated time remaining:", remaining, " "))
    #i count
    i <- i+1
  }
  
  out_df <- out_df %>%
    select(genome, everything())
  
  cat("done\n")
  return(out_df)
}

