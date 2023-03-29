get_population_size <- function(region) {
  d <- readRDS("data/pop_by_province_20230213.rds")
  if (region == "overall") {
    pop <- sum(d$Total)
  } else{
    pop <- d[d$province_clean == region, ]$Total
    # pop <- d[province_clean == region, Total]
  }
  return (pop)
}
