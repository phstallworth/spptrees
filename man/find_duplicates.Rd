% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/Paper_BW_Cross_Validation.R
\name{find_duplicates}
\alias{find_duplicates}
\title{Find the duplicate rows in a "final" matrix}
\usage{
find_duplicates(my_matrix)
}
\arguments{
\item{my_matrix}{A 6 row matrix with the first row validation set info,
the second and third row lattice points, the 4th and 5th row
nearest tree neighbor, and the 6th distance}
}
\value{
A vector with the row number of each non-duplicate entry
}
\description{
Find the duplicate rows in a "final" matrix
}

