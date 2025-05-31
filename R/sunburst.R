# Code from Google Gemini

# https://alleninstitute.github.io/scrattch.example/Sunburst.html#
# Code by Cindy van Velthoven


as.sunburstDF <- function(cl.df, 
                          levels=c("class","subclass","cluster"),
                          valueCol = NULL, 
                          rootname="total") {
  
  ## select columns for plotting (levels)
  group_columns <- function(grouping) {
    list(id = paste0(grouping, "_id"),
         label = paste0(grouping, "_label"),
         color = paste0(grouping, "_color"))
  }
  
  group_cols <- list()
  
  for(i in 1:length(levels)){
    group_cols[[i]] <- group_columns(levels[i])
  }
  
  group_cols <- unlist(group_cols)
  
  sub.df <- cl.df[,group_cols]
  
  #add value column to level.df
  if(!is.null(valueCol)) {
    sub.df[,valueCol] <- cl.df[,valueCol]
  }
  
  ## set root (center)
  sub.df$root <- rootname
  sub.df <- data.frame(root = rootname, sub.df)
  
  
  ## add the first label set == all labels == root
  
  if(is.null(valueCol)) {
    hierarchyDF.base <- data.frame(cl.id ="0",
                                   labels=rootname,
                                   color= "white",
                                   parent = "",
                                   ids = rootname,
                                   values=nrow(sub.df))
  } else {
    hierarchyDF.base  <- data.frame(cl.id ="0",
                                    labels=rootname,
                                    color= "white",
                                    parent = "",
                                    ids = rootname,
                                    values=sum(sub.df[valueCol]))
  }
   
  hierarchyList <- list()
  
  for(i in seq_along(levels)){
    print(i)
    
    lev_labs <- paste0(levels, "_label") 
    currentCols <- c("root", lev_labs)
    
    parentCols <- currentCols[1:i]
    sub.df$parent <-do.call(paste, c(sub.df[parentCols], sep="-"))
    idCols <- currentCols[1:(i+1)]
    sub.df$ids <- do.call(paste, c(sub.df[idCols], sep="-"))
    
    lev_lab <- paste0(levels[i], "_label")
    lev_col <- paste0(levels[i], "_color")
    lev_id <- paste0(levels[i], "_id")
    par <- "parent"
    id <- "ids"
    
    if(is.null(valueCol)) {
      currentDF <- sub.df %>%
        group_by(.data[[lev_id]],
                 .data[[lev_lab]], 
                 .data[[lev_col]],
                 .data[[par]],
                 .data[[id]]) %>% 
        summarise(values=n()) %>%
        arrange(.data[[lev_id]]) 
    } else {
      currentDF <- sub.df %>%
        group_by(.data[[lev_id]],
                 .data[[lev_lab]], 
                 .data[[lev_col]],
                 .data[[par]],
                 .data[[id]]) %>% 
        summarise(values=sum(.data[[valueCol]])) %>%
        arrange(.data[[lev_id]])
    }
    
    colnames(currentDF) <- c("cl.id", "labels", "color","parent","ids", "values")
    hierarchyList[[i]] <- currentDF
  }
  
  hierarchyDF <- data.table::rbindlist(hierarchyList, use.names = TRUE, fill = TRUE)
  hierarchyDF <- rbind(hierarchyDF.base, hierarchyDF)
  return(hierarchyDF)
}



