% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{find_indels}
\alias{find_indels}
\title{Identify and subset InDel markers from a marker panel.}
\usage{
find_indels(x, rp_row, dp_row, indel_sym = "-", sep = ":", na_code = NA)
}
\arguments{
\item{x}{A data matrix or frame with markers as columns and samples as rows.}

\item{rp_row, dp_row}{An integer or character value indicating the row index
or name of Parent 1 and 2.}

\item{indel_sym}{A character value that indicates the symbols for a deletion.}

\item{sep}{A character used as separator for genotype calls, default is a
colon.}

\item{na_code}{A value indicating missing data in \code{x}.}
}
\value{
A list object with the following components:
1) data frame of InDel loci in marker panel if present.
2) data frame of non-InDel markers in marker panel.
}
\description{
Identify and subset InDel markers from a marker panel.
}
\examples{
# example code
library(panGenomeBreedr)

# Marker data
dat <- data.frame(snp1 = c('C:C', 'A:A', 'C:A', 'C:A'),
                  snp2 = c('C:C', '-:-', 'C:C', 'C:C'),
                  snp3 = c('T:T', 'C:C', 'C:T', 'C:T'),
                  snp4 = c('G:G', '-:-', 'G:-', 'G:G'),
                  snp5 = c('T:T', 'A:A', 'T:A', 'T:A'),
                  row.names = c('rp', 'dp', 'art_het1', 'art_het2'))

# Find InDel loci
geno_indel <- find_indels(x = dat,
                          rp_row = 1,
                          dp_row = 2)$geno_indel


}
