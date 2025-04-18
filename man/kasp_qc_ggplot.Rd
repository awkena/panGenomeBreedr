% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kasp_validation.R
\name{kasp_qc_ggplot}
\alias{kasp_qc_ggplot}
\title{Make KASP marker genotyping QC plot.}
\usage{
kasp_qc_ggplot(
  x,
  FAM = "X",
  HEX = "Y",
  geno_call = "Call",
  color = "Color",
  snp_id = "SNPID",
  blank = "NTC",
  Group_id = NULL,
  scale = FALSE,
  pdf = TRUE,
  width = 6,
  height = 6,
  filename = "kasp_qc",
  expand_axis = 0.5,
  legend.pos.x = 0.6,
  legend.pos.y = 0.75,
  legend.box = "horizontal",
  legend.pos = "inside",
  alpha = 0.5,
  text_size = 12,
  ...
)
}
\arguments{
\item{x}{A list object of KASP genotype calls processed by the `kasp_color()`
function.}

\item{FAM}{A character indicating the column name for FAM fluorescence
coordinates in \code{x}.}

\item{HEX}{A character indicating the column name for HEX fluorescence
coordinates in \code{x}.}

\item{geno_call}{A character value indicating the column name for KASP genotype
calls in \code{x}.}

\item{color}{A character value indicating the column name for assigned colors
in \code{x}.}

\item{snp_id}{A character value indicating the column name for SNP IDs
in \code{x}.}

\item{blank}{A character value indicating `No Template Controls (NTC)`
genotype calls.}

\item{Group_id}{A character value for the column ID indicating the predictions
of the positive controls in \code{x}.}

\item{scale}{A logical value indicating whether to scale FAM and HEX axis
to values between 0 and 1.}

\item{pdf}{A logical value indicating whether to save plot as a pdf graphic
device when TRUE or output plot in R when FALSE.}

\item{width}{A numeric value for the width of pdf device.}

\item{height}{A numeric value for the height of pdf device.}

\item{filename}{A character value for path or file name for saving pdf.}

\item{expand_axis}{A numeric value to expand the axes for legend positioning.}

\item{legend.pos.x}{A numeric value representing the x coordinate for legend
placement.}

\item{legend.pos.y}{A numeric value representing the y coordinate for legend
placement.}

\item{legend.box}{A character value of either `horizontal` or `vertical`
legend placement.}

\item{legend.pos}{A character value for the position of legend; the
default value is `inside`.}

\item{alpha}{A numeric value between 0 and 1 for modifying the
opacity of colors.}

\item{text_size}{A numeric value for setting text size.}

\item{...}{Other valid arguments that can be passed to ggplot2.}
}
\value{
A graphic object or plot.
}
\description{
Make KASP marker genotyping QC plot.
}
\examples{
# example code
library(panGenomeBreedr)
\donttest{
# Assign KASP colors to plates
dat1 <- kasp_color(x = panGenomeBreedr::kasp_dat,
                   subset = 'MasterPlate',
                   sep = ':',
                   geno_call = 'Call',
                   uncallable = 'Uncallable',
                   unused = '?',
                   blank = 'NTC')

# KASP QC plot for Plate 12
kasp_qc_ggplot(x = dat1[12],
                    pdf = FALSE,
                    Group_id = 'Group',
                    scale = TRUE,
                    expand_axis = 0.6,
                    alpha = 0.5,
                    legend.pos.x = 0.6,
                    legend.pos.y = 0.8)
}

}
