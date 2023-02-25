
#' Extracting an Example Imputation from the ImputeBench Protocol
#'
#' @description Returning a data matrix, a mask matrix and a set of imputed matrices for a specific run and scenario/parameter choice
#' of either `simulate_ImputeBench()` or `data_ImputeBench()`.
#'
#' @param Evaluation The output of `data_ImputeBench` or `simulation_ImputeBench`.
#' @param parameters A positive integer. Indicating the number of the scenario/parameter choice that should be extracted. Default is `1`.
#' @param run A positive integer. Indicating the run that should be extracted. Default is `1`.
#' @param data A numeric matrix with possible missing entries or `NULL`. In the case `Evaulation` stems from evaluation on real data, i.e. as the output of `data_ImputeBench()`,
#' the real data matrix needs to be provided via `data`. Default is `NULL`.
#' @param imputation_methods A vector of character strings chosen from `knn`, `MICE`, `missF`, `softImpute` and `median`. Indicating
#' the default methods that should be used to impute the data matrix. Default is `c("knn", "MICE", "missF", "softImpute", "baseline")`.
#' @param imputation_indices A vector of positive integers of the same length as `imputation_methods` of `NULL`. Needs to be provided
#' if `imputation_methods` is not `NULL`. The vector corresponds to the queue positions of the given default methods in the
#' benchmarking process. Default is `c(1,2,3,4,5)`.
#'
#'
#' @details This function extracts from the output of `simulate_ImputeBench()` or `data_ImputeBench()` a specific run of a specific
#' scenario or parameter choice. Integrated is also the imputation via default methods given the trained or fixed parameters in the
#' benchmarking protocol (in particular `k` of `knn` and `lambda` of `softImpute` are by default trained). If imputation is not wanted
#' one sets `imputation_methods` to `NULL`. Otherwise one needs to pass also the queue positions of the default methods. In detail,
#' if all default methods were set in the run of `simulate_ImputeBench()` or `data_ImputeBench()` and not additional methods were
#' used the queue positions are 1 for `knn`, 2 for `MICE`, 3 for `missF`, 4 for `softImpute`, 5 for the baseline method `baseline`.
#' In the case that additional methods were being benchmarked, the default methods are always at the end of the queue. Hence, their
#' positions are given by adding the number of methods considered. For example if all default methods are considered and 4 additional
#' non-default methods were being passed to say `simulate_ImputeBench()`, the queue position of `knn` is 5, that of `MICE` is 6, and
#' so on. Last, if a default method is not part of the benchmarking process, the "later" default methods are moving a place up. For
#' example if we consider all default methods but `missF` and no additional ones, the queue positions of `softImpute` and `baseline`
#' are now 4 and 5, respectively.
#'
#'
#' @return A ggplot2 plot.
#'
#' @import ggplot2
#' @import stringr
#'
#' @export
#'

example_ImputeBench = function(Evaluation, parameters = 1, run = 1, data = NULL,
                               imputation_methods = c("knn", "MICE", "missF", "softImpute", "baseline"),
                               imputation_indices = c(1,2,3,4,5)){

  if(is.null(data)){
    data = Evaluation$data.list[[parameters + 1]][[run]]$data
  }
  mask = Evaluation$data.list[[parameters + 1]][[run]]$mask

  miss.data = data
  miss.vec = c(miss.data)
  miss.mask = c(mask)
  miss.vec[which(miss.mask == 0)] =  NA
  miss.data = matrix(miss.vec, ncol = ncol(data))

  # results$data.list[[2]][[1]]$args[[4]]

  if(!is.null(imputation_methods)){
    args = Evaluation$data.list[[parameters + 1]][[run]]$args
    imputed_data = list()

    if("softImpute" %in% imputation_methods){
      imputed_data$softImpute = round_imputation(data = miss.data,
                                                 imputed_data = wsImpute(data = miss.data,
                                                                         args = args[[imputation_indices[which(imputation_methods == "softImpute")]]]))

    }
    if("baseline" %in% imputation_methods){
      imputed_data$median = round_imputation(data = miss.data,
                                             imputed_data =  baseline_imputation(data = miss.data,
                                                                                 args = args[[imputation_indices[which(imputation_methods == "baseline")]]]))
    }
    if("knn" %in% imputation_methods){
      imputed_data$knn =  round_imputation(data = miss.data,
                                           imputed_data = wKNN(data = miss.data,
                                                               args = args[[imputation_indices[which(imputation_methods == "knn")]]]))
    }
    if("missF" %in% imputation_methods){
      imputed_data$missF =  round_imputation(data = miss.data,
                                             imputed_data = wmissF(data = miss.data,
                                                                   args = args[[imputation_indices[which(imputation_methods == "missF")]]]))
    }
    if("MICE" %in% imputation_methods){
      imputed_data$MICE =  round_imputation(data = miss.data,
                                            imputed_data = wMICE(data = miss.data,
                                                                 args = args[[imputation_indices[which(imputation_methods == "MICE")]]]))
    }
  }

  if(is.null(imputation_methods)){
    example = list("data" = data,
                   "mask" = mask)
  } else{
    example = list("data" = data,
                   "mask" = mask,
                   "imputed_data" = imputed_data )
  }

  return(example)
}


