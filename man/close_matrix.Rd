% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/base_sim.R
\name{close_matrix}
\alias{close_matrix}
\title{Close-Ties matrix generator}
\usage{
close_matrix(close_matrix)
}
\arguments{
\item{close_matrix}{An object that can be coerced into a matrix. Values must be either 0 or 1 to denote connections.}
}
\value{
Returns to the global environment a close adjacency matrix labeled adj_close.
}
\description{
A renamed as.matrix function to take in an object that can be coerced into a matrix and returns it to the environment as a matrix ready for simulation.
}
