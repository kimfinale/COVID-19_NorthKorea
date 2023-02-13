extract_dates <- function(string, month="월", day="일"){
  md <- stringr::str_extract_all(string, pattern=paste0("\\d+(?=", month, "|", day, ")"))
  df <- do.call('rbind', md)
  dates <- as.Date(paste0("2022-", df[,1], "-", df[,2]))
  
  return(dates)
}


