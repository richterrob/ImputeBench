% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/log_odds.R
\name{log_odds}
\alias{log_odds}
\title{Computing the Log-Odds of an Imputed Matrix}
\usage{
log_odds(data, imputed_data, mask)
}
\arguments{
\item{data}{A matrix with binary entries or \code{NA} entries marking missing entries.}

\item{imputed_data}{A binary imputed matrix of \code{data}.}

\item{mask}{A matrix with binary entries (0 = missing, 1 = not missing), the (additional) missingness pattern on which the
imputed matrix is compared against the ground truth matrix.}
}
\value{
The log-odds error/loss of the imputed matrix.
}
\description{
Computing the log-odds error/loss for imputed binary matrices.
}
\details{
Hard coded is an \code{epsilon = 1e-15} to avoid division by 0.
}
