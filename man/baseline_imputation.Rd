% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/baseline_imputation.R
\name{baseline_imputation}
\alias{baseline_imputation}
\title{Baseline (Mean/Median) Imputation Method}
\usage{
baseline_imputation(data, args = list(method = "median"))
}
\arguments{
\item{data}{A numeric matrix with missing entries given by \code{NA}.}

\item{args}{A list of arguments consisting of a single entry \code{method} (possible entries are \code{"median"} (default) and \code{"mean"}).}
}
\value{
A numeric matrix with imputed missing entries.
}
\description{
Baseline (Mean/Median) Imputation Method
}
\details{
The function of the baseline imputation method of ImputeBench. Supported are mean and median imputation.
}
