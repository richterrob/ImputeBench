

#' Rounding the imputed data matrix
#'
#' @description Taking as an input the data matrix (with missing entries) and the imputed data
#' matrix this function tests which columns contain only positive integers and which contain only
#' 0's and 1's. Then it rounds the respective columns of the imputed data matrix accordingly.
#'
#' @param data A matrix with numeric entries or `NA` entries marking missing entries.
#' @param imputed_data A matrix with numeric entries of the same dimension as `data` with no `NA` entries.
#'
#' @return The appropriately rounded imputed data matrix.
#'
#' @export
#'



round_imputation = function(data, imputed_data){

  cb.clms = determine_count_binary(data)

  for(t in 1:ncol(imputed_data)){
    if(t %in% cb.clms[[1]]){
      imputed_data[,t] = pmax(round(imputed_data[,t]),0)
    }
    if(t %in% cb.clms[[2]]){
      imputed_data[,t] = pmin(imputed_data[,t],1)
    }
  }

  return(imputed_data)
}
