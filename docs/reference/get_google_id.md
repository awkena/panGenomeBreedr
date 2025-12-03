# Get the folder or file ID from a Google Drive shareable link.

Get the folder or file ID from a Google Drive shareable link.

## Usage

``` r
get_google_id(drive_link, is.folder = TRUE)
```

## Arguments

- drive_link:

  A character value indicating the shareable Google Drive link.

- is.folder:

  A logical value indicating if link is for a folder or file. Set to
  \`FALSE\` if link is for a shareable file.

  \#' @examples \# example code library(panGenomeBreedr) folder_link \<-
  "https://drive.google.com/drive/folders/1BotxaUb5emlrtgo473db3gDTUCLzKi70?usp=sharing"
  folder_id \<- get_google_id(drive_link = folder_link)

## Value

A character object of Google Drive folder or file ID.
