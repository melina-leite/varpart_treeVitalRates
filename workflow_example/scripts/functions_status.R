# Function for STATUS standardization


### TREE level status ----
# create tree level status information from stem level

# aggregating tree data from the stem data

# criteria to create the status of the tree
#- if at least one stem is alive, tree is A
#- if no stem is alive and there is at least one stem dead, tree is D
#- if no stem is alive or dead and there is at least one stem prior, tree is P
#- if tree is missing or NA in all the stem status -> M

treestat <-   function(x) {
  if (length(grep("A",x))>0) {"A"
  } else {
    if (length(grep("D",x))>0) {"D"
    } else {
      if(length(grep("P", x))>0) {"P"
      } else  {"M"}
    }
  }
}


### CORRECT STATUS ----

# correct sequence of status based on the discoveries of 1) A and 2) D
correction = function(x) {
  
  # trees that are recorded living at least once
  if (any(x == "A", na.rm = T)) {
    startA = min(which(x == "A"))
    endA = max(which(x == "A"))
    
    if (any(x[1:(startA-1)] == "D", na.rm = T)) {
      
      # trees that are falsely recorded as dead and then found alive
      startD = min(which(x == "D"))
      x[1:(startD-1)] = "P"                                # set everything before first D to P
      x[startD:endA] = "A"                                 # set everything between first D and last A to A
    } else {
      
      # trees that are not found dead before alive
      x[1:(startA-1)] = "P"                                # set everything before first A to P
      x[startA:endA] = "A"                                 # set everything between first and last A to A
    }
    
    if (endA < length(x)) {
      if (is.na(x[endA+1])) x[(endA+1):length(x)] = "M"  # set everything after last A to M if no D exists
      if (x[endA+1] == "D") x[(endA+1):length(x)] = "D"  # set everything after last A to D otherwise
    }
    
    # trees that are never recorded living
  } else {
    
    # trees are recorded dead at least once
    if (any(x == "D", na.rm = T)) {                      
      startD = min(which(x == "D"))
      x[startD:length(x)] = "D"                          # set everything after first D to D
      if (startD != 1) x[1:(startD-1)] = "P"             # set everything before first D to P
      
      # trees that are only recorded with M or NA
    } else {
      startM = min(which(x == "M"))
      x[startM:length(x)] = "M"                          # set everything after first M to M
      if (startM != 1) x[1:(startM-1)] = "P"             # set everything before first M to P
    }
  }
  x
}
#x = c("AP", "DP", "A")

# using the function above to generate again the correct status for the
# stem list object after accounting for all other issues.

correct_status <- function(x){
  
  # add arbitrary stemID if x is tree level data
  level = ifelse("stemID" %in% unique(unlist(lapply(x, function(x) names(x)))), "stem", "tree")
  if (level == "tree") {
   x = lapply(x, transform, stemID = 1) 
  }
  
  # combine all censuses
  dat <- bind_rows(x, .id="census") %>% 
    mutate(census = paste0("status_",substr(census, 5,5)))
  
  # make wide data format for status only
  status_wide <- dat %>% 
    pivot_wider(id_cols = c(treeID, stemID), 
                names_from = census,
                values_from = status) %>%  
    arrange(treeID, stemID)
  
  # get only status columns
  status <- status_wide %>% 
    select(starts_with("status_")) %>%
    mutate_all(as.character)
  
  # apply correction
  status_c <- t(apply(status, 1, correction)) %>% 
    as.data.frame() %>%
    mutate_all(as.character)
  
  
  # replace with corrected status columns
  for (i in 1:length(x)) {
    x[[i]]$status = status_c[match(paste(x[[i]]$treeID, x[[i]]$stemID),
                                   paste(status_wide$treeID, status_wide$stemID)), i]
  }
  
  # remove arbitrary stemID and hom if x is tree level data
  if (level == "tree") {
    x = lapply(x, select, -stemID)
  }
  
  return(x)
  
 
}
