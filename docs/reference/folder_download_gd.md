# Download files in a shared Google Drive folder without restrictions.

Download files in a shared Google Drive folder without restrictions.

## Usage

``` r
folder_download_gd(drive_link, output_path = tempdir(), is.folder = TRUE)
```

## Arguments

- drive_link:

  A character value indicating the shareable Google Drive link.

- output_path:

  A character value indicating the path to a directory for saving
  downloaded files.

- is.folder:

  A logical value indicating if link is for a folder or file. Set to
  \`FALSE\` if link is for a shareable file.

## Value

A list or vector containing the path to directory containing downloaded
files from Google Drive.

## Examples

``` r
# example code
# \donttest{
library(panGenomeBreedr)
f_link <- "https://drive.google.com/drive/folders/1BotxaUb5emlrtgo473db3gDTUCLzKi70?usp=sharing"
folder_path <- folder_download_gd(drive_link = f_link)
#> ℹ Not logged in as any specific Google user.
#> File downloaded:
#> • BR Field layout <id: 12JYxAz8OCmRLS1WIa0q9ohTB50wQAks-EclJja5l-4U>
#> Saved locally as:
#> • /var/folders/n_/swy48fpx1w76xyqp3qx2prz00000gn/T//Rtmpj16bNQ/Work Plan 2025/BR Field layout.xlsx
#> File downloaded:
#> • BR workplan for 2025 <id: 1ghJuIyhzRdVquSXepUh64I99xNcj2wvtH4tV5Qjsf_k>
#> Saved locally as:
#> • /var/folders/n_/swy48fpx1w76xyqp3qx2prz00000gn/T//Rtmpj16bNQ/Work Plan 2025/BR workplan for 2025.pptx
# }
```
