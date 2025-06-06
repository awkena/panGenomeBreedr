% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{parent_het}
\alias{parent_het}
\title{Identify and subset loci with any heterozygous parent genotype.}
\usage{
parent_het(x, rp_row, dp_row, sep = ":", na_code = NA)
}
\arguments{
\item{x}{A data matrix or frame with markers as columns and samples as rows.}

\item{rp_row, dp_row}{An integer or character value indicating the row index
or name of Parent 1 and 2.}

\item{sep}{A character used as separator for genotype calls, default is a
colon.}

\item{na_code}{A value indicating missing data in \code{x}.}
}
\value{
A list object with the following components:
1) data frame of loci with at least one hetrozygous parent genotype, if present.
2) data frame of loci with all homozygous parent genotype.
}
\description{
Identify and subset loci with any heterozygous parent genotype.
}
\examples{
# example code
library(panGenomeBreedr)

# Marker data
dat <- data.frame(snp1 = c('C:A', 'A:A', 'C:A', 'C:A'),
                  snp2 = c('C:C', 'G:G', 'C:C', 'C:C'),
                  snp3 = c('C:T', 'C:C', 'C:T', 'C:T'),
                  snp4 = c('G:G', '-:-', 'G:-', 'G:G'),
                  snp5 = c('T:T', 'A:A', 'T:A', 'T:A'),
                  row.names = c('rp', 'dp', 'ind_1', 'ind_2'))

# Find loci with at least one heterozygous parent genotype
par_het <- parent_het(x = dat,
                      rp_row = 1,
                      dp_row = 2,
                      sep = ':')$par_het

}
