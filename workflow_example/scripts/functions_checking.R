

# functions for checking ForestGEO data

library(crayon)


# function that allows running functions without showing cat() output (silent checks)
hush = function(code){
  sink("NUL") # use /dev/null in UNIX
  tmp = code
  sink()
  return(tmp)
}


# Check dimensions --------------------------------------------------------

# Do all censuses have the same dimensions, i.e. is the dataset appropriately expanded?
check_dim = function(data) {
  cat("\n")
  
  cat(underline("Check dimension"), sep = "\n")
  dims = as.data.frame(lapply(data, dim))
  res = apply(dims, 1, function(x) ifelse(length(unique(x))==1, 
                                          "dimensions are identical in censuses",
                                          "dimensions are NOT identical in censuses"))
  check = all(apply(dims, 1, function(x) length(unique(x))==1))
  cat(do.call(ifelse(check, "black", "red"), 
              list(paste(c("Rows:", "Columns:"), res))), sep = "\n")
  invisible(check)
}



# Check columns -----------------------------------------------------------

# Do all censuses have the correct columns?
check_col = function(data) {
  cat("\n")
  
  cat(underline("Check columns"), sep = "\n")
  # expected column names
  if (all(grepl("stem", names(data)))) {
    exp = c("treeID", "stemID", "sp", "gx", "gy", "dbh", "hom", 
            "homchange", "ba", "status", "date")
  } else {
    exp = c("treeID", "sp", "gx", "gy", "dbh",
                 "homchange", "ba", "nostems", "status", "date")
  }
  res = unlist(lapply(data, function(x) identical(names(x), exp)))
  check = all(res)
  message = ifelse(check, 
                   "column names are okay",
                   "column names are NOT okay")
  
  cat(do.call(ifelse(check, "black", "red"), list(message)),
      sep = "\n")
  invisible(check)
}



# Check column types ------------------------------------------------------

# Do all censuses have the correct columns?
check_types = function(data) {
  cat("\n")
  
  cat(underline("Check column types"), sep = "\n")
  # expected column types
  if (all(grepl("stem", names(data)))) {
    exp = c("character", "character", "character", "numeric", "numeric", "numeric", "numeric", 
            "logical", "numeric", "character", "Date")
  } else {
    exp = c("character", "character", "numeric", "numeric", "numeric",
            "logical", "numeric", "integer", "character", "Date")
  }

  res = unlist(lapply(data, function(x) identical(as.vector(sapply(x, class)), exp)))
  check = all(res)
  message = ifelse(check, 
                   "column types are okay",
                   "column types are NOT okay")
  
  cat(do.call(ifelse(check, "black", "red"), list(message)),
      sep = "\n")
  invisible(check)
}



# Check order -------------------------------------------------------------

# Is the order of stems and trees identical in each census?
check_order = function(data) {
  
  cat("\n")
  
  # run only if dimensions are okay
  if (hush(check_dim(data))) {
    cat(underline("Check order"), sep = "\n")
    IDs = as.data.frame(lapply(data, function(x) x$treeID))
    res = apply(IDs, 1, function(x) length(unique(x))==1) # same treeID in each row?
    check = all(res)
    message = ifelse(check, 
                     "treeIDs have the same order",
                     "treeIDs don't have the same order")
    cat(do.call(ifelse(check, "black", "red"), list(message)), 
        sep = "\n")
    
    if ("stemID" %in% unique(unlist(lapply(data, function(x) names(x))))) {
      IDs = as.data.frame(lapply(data, function(x) x$stemID))
      res = apply(IDs, 1, function(x) length(unique(x))==1) # same stemID in each row?
      check = all(res)
      message = ifelse(check, 
                       "stemIDs have the same order", 
                       "stemIDs don't have the same order")
      cat(do.call(ifelse(check, "black", "red"), list(message)), 
          sep = "\n")
    }
  } else {
    cat(underline("Check order of IDs NOT possible due to wrong dimensions."),
        sep = "\n")
  }
  invisible(check)
}




# Check status ------------------------------------------------------------

# Do all status entries follow the status convention?
# Stem level data: comparing corrected status with original one

check_status = function(data, orig.data) {
  options(dplyr.summarise.inform = FALSE)
  
  if(is.null(orig.data)){
    cat("\n")
    cat(underline("Check status"), sep = "\n")
    
    stat = unique(unlist(lapply(data, function(x) unique(x$status))))
    allowed = c("A", "D", "P", "M")
    not_allowed = sort(stat[!stat %in% allowed])
    
    check = length(not_allowed)==0
    message = ifelse(check, 
                     "status entries follow convention",
                     paste("the following status entries don't follow the convention:\n",
                           paste(not_allowed, collapse = ", ")))
    cat(do.call(ifelse(check, "black", "red"), list(message)), 
        sep = "\n")
  } else{
    cat("\n")
    cat(underline("Check status"), sep = "\n")
    
    stat = unique(unlist(lapply(data, function(x) unique(x$status))))
    allowed = c("A", "D", "P", "M")
    not_allowed = sort(stat[!stat %in% allowed])
    
    check = length(not_allowed)==0
    message = ifelse(check, 
                     "status entries follow convention",
                     paste("the following status entries don't follow the convention:\n",
                           paste(not_allowed, collapse = ", ")))
    cat(do.call(ifelse(check, "black", "red"), list(message)), 
        sep = "\n")
    
    # check the differences between original and corrected status
    
    orig <- bind_rows(orig.data, .id="census") %>% select(census,status)
    orig$status[which(is.na(orig$status))] <- "NA"
    correct <- bind_rows(data, .id="census") %>% select(census,status)
    correct$status[which(is.na(correct$status))] <- "NA"
    correct$orig.status <- orig$status
    correct$change <- correct$orig.status != correct$status
    correct$dif <- paste0(correct$orig.status,"_", correct$status)
    
    corr.t <- correct %>% group_by(census, dif) %>% summarise(n=sum(change)) %>%
      filter(n>0) %>%
      pivot_wider(names_from = dif, values_from=n) %>% as.data.frame()
    rownames(corr.t) <- corr.t$census
    
    corr.p <- round(corr.t[,-1, drop=FALSE]*100/dim(data[[1]])[1],6)
    
    n.change <- sum(orig$status != correct$status)
    
    message2 = ifelse(n.change==0, "No status was corrected",
                      paste0("There were ",n.change,
                             " status corrected. Corrections made:"))
    
    
    cat(do.call(ifelse(n.change==0, "black", "red"), list(message2)), 
        sep = "\n")
    if(n.change!=0) {
      print(corr.t[,-1, drop=FALSE])
      cat("Percentage of change from total census data", sep="\n")
      print(corr.p)
    }
    cat("", sep = "\n")
    
    
    # summary status values for each census
    sum_status = lapply(data, function(x) {
      x%>% tabyl(status) %>% select(-percent)
    })
    all_status <- bind_rows(sum_status, .id="census") %>%
      pivot_wider(names_from = status, values_from=n) %>% as.data.frame()
    rownames(all_status) <- all_status$census
    cat("summary of status N", sep = "\n")
    print(all_status[,-1])
    cat("", sep = "\n")
    
    cat("summary of status %", sep = "\n")
    print(round(all_status[,-1]*100/rowSums(all_status[,-1], na.rm = T)[1],2))
    cat("", sep = "\n")
  }
}




# Check units -------------------------------------------------------------

check_units = function(data) {
  
  cat("\n")
  
  cat(underline("Check units"), sep = "\n")
  
  par(mfrow = c(2, 1))
  
  if ("dbh" %in% unique(unlist(lapply(data, function(x) names(x))))) {
    plot(10000, 100000, xlim = c(1, 10000), ylim = c(0, 2), log = "x",
         xlab = "dbh in mm", ylab = "")
    polygon(c(10, 10, 1000, 1000), c(-10, 1000, 1000, -10), 
            col = "lightgrey", border = NA)
    sum_dbh = lapply(data, function(x) {
      x = x[x$status == "A", ]
      d = density((x$dbh[!is.na(x$dbh)]), bw = 0.05)
      d$x = (d$x)
      lines(d)
      c(min = min(x$dbh, na.rm = T), median = median(x$dbh, na.rm = T), 
        mean = mean(x$dbh, na.rm = T), max = max(x$dbh, na.rm = T))
    })
    box()
  }
  
  # basal area
  plot(10000, 100000, xlim = c(0.000001, 50), ylim = c(0, 0.4), log = "x",
       xlab = "ba in m2", ylab = "")

  polygon(c(pi*((10/1000)/2)^2, pi*((10/1000)/2)^2, 
            pi*((1000/1000)/2)^2, pi*((1000/1000)/2)^2), 
          c(-10, 1000, 1000, -10), 
          col = "lightgrey", border = NA)
  sum_ba = lapply(data, function(x) {
    x = x[x$status == "A", ]
    d = density(log((x$ba[!is.na(x$ba)])), bw = 0.5)
    d$x = exp(d$x)
    lines(d)
    c(min = min(x$ba, na.rm = T), median = median(x$ba, na.rm = T), 
      mean = mean(x$ba, na.rm = T), max = max(x$ba, na.rm = T))
  })
  box()
  
  if ("dbh" %in% unique(unlist(lapply(data, function(x) names(x))))) {
    cat("summary of dbh (should be in mm)", sep = "\n")
    print(t(as.data.frame(sum_dbh)))
    cat("", sep = "\n")
  }
  
  cat("summary of ba (should be in m2)", sep = "\n")
  print(round(t(as.data.frame(sum_ba)), 3))
}





# Check NAs ---------------------------------------------------------------

check_nas = function(data) {
  cat("\n")
  cat(underline("Check NAs"), sep = "\n")
  cat("percent of observations with NA per census", sep = "\n")
  
  # new column for querying only NA dates in living stems/trees
  data = lapply(data, function(x) {
    x$date_alive = x$date
    x$date_alive[x$status != "A" & is.na(x$date)] = Sys.Date() # sets NA to an arbitrary date for non-living stems/trees
    x
  })
  
  # the same for dbh
  if ("dbh" %in% unique(unlist(lapply(data, function(x) names(x))))) {
    data = lapply(data, function(x) {
      x$dbh_alive = x$dbh
      x$dbh_alive[x$status != "A" & is.na(x$dbh)] = 10 # sets NA to an arbitrary dbh for non-living stems/trees
      x
    })
  }
  
  # the same for hom
  data = lapply(data, function(x) {
    x$hom_alive = x$hom
    x$hom_alive[x$status != "A" & is.na(x$hom)] = 1.3 # sets NA to an arbitrary hom for non-living stems/trees
    x
  })
  
  if ("dbh" %in% unique(unlist(lapply(data, function(x) names(x))))) {
    columns = c("sp", "gx", "gy", "status", "date_alive", "dbh_alive", "hom_alive")
  } else {
    columns = c("sp", "gx", "gy", "status", "date_alive", "hom_alive")
  }
  na_proportion = t(as.data.frame(lapply(data, function(x) {
    apply(is.na(x[, columns]), 2, sum) / nrow(x)
  })))
  print(round(na_proportion *100, 3))
}



# Check coordinates boundaries ---------------------------------------------

# Do all stems/trees are inside plot boundaries?
check_coor = function(data) {
  cat("\n")
  
  cat(underline("Check stem/tree boundaries"), sep = "\n")
  dimX = as.data.frame(lapply(data, function(x) range(x$gx, na.rm=T)))
  lowX = dimX[1,] >= t(plot$plotdimension[1,1])
  upX  = dimX[2,] <= t(plot$plotdimension[1,2])
  dimY = as.data.frame(lapply(data, function(x) range(x$gy, na.rm=T)))
  resX <- ifelse(all(lowX, upX), 
                 "all trees/stems are inside X plot boundaries",
                 "there are trees/stems are outside X plot boundaries")
  dimY <= t(plot$plotdimension[2,])
  lowY = dimY[1,] >= t(plot$plotdimension[2,1])
  upY  = dimY[2,] <= t(plot$plotdimension[2,2])
  resY <- ifelse(all(lowY,upY), 
                 "all trees/stems are inside Y plot boundaries",
                 "there are trees/stems are outside Y plot boundaries")
  
  check = all(all(lowX, upX),all(lowY, upY))
  cat(do.call(ifelse(check, "black", "red"), 
              list(paste(c("X:", "Y:"), c(resX, resY)))), sep = "\n")
#  invisible(check)
}


# Check all ---------------------------------------------------------------

check_all = function(data, orig.data=NULL) {
  cat(bold(paste("Quality checks for", attributes(data)$site, 
                 substr(names(data)[1], 1, 4), "data")), sep = "\n")
  check_dim(data)
  check_col(data)
  check_types(data)
  check_order(data)
  check_status(data, orig.data)
  check_units(data)
  check_nas(data)
  check_coor(data)
}



