% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/variant_annotation.R
\name{get_google_id}
\alias{get_google_id}
\title{Get the folder or file ID from a Google Drive shareable link.}
\usage{
get_google_id(drive_link, is.folder = TRUE)
}
\arguments{
\item{drive_link}{A character value indicating the shareable Google Drive link.}

\item{is.folder}{A logical value indicating if link is for a folder or file.
Set to `FALSE` if link is for a shareable file.

#' @examples
# example code
library(panGenomeBreedr)
folder_link <- "https://drive.google.com/drive/folders/1BotxaUb5emlrtgo473db3gDTUCLzKi70?usp=sharing"
folder_id <- get_google_id(drive_link = folder_link)}
}
\value{
A character object of Google Drive folder or file ID.
}
\description{
Get the folder or file ID from a Google Drive shareable link.
}
