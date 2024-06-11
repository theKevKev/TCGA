#' Basic Clean: Removes columns with 0 variance, and imputes all NA values
#'
#' @param data the data to be cleaned (columns are parameters, rows are samples)
#' @param removeZeros whether columns with too many zeros should be removed
#' @param zero_cutoff the cutoff for how many 0s imply a column to be removed
#' @param zerosAsNA whether 0s should be imputed as well as NA
#' @param knn whether to perform knn imputation (default TRUE)
#'
#' @return an unnormalized, but cleaner data set
#' @export
#'
#' @examples
basic_clean <- function(data,
                        zerosAsNA = TRUE,
                        removeZeros = TRUE,
                        zero_cutoff = 0.3,
                        knn = TRUE
) {
  ## Remove columns with var = 0, (removes 0-cols too)
  data <- data[, !apply(data, 2, function(x){var(x, na.rm = TRUE)}) == 0]
  
  if(removeZeros) {
    ## Remove columns that with more than x% of zeros
    data <- data[, colSums(data == 0, na.rm = TRUE)/nrow(data) < zero_cutoff]
  }
  
  if(zerosAsNA) {
    data[data == 0] <- NA
  }
  
  if (knn) {
    ## KNN imputation
    data <- t(data)
    data <- impute::impute.knn(data)
    data <- t(data$data)
  } else {
    stop("No Imputation Method Provided")
  }
  
  return(data)
}