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
f_link <- "https://drive.google.com/drive/folders/1ZJ2cZcybCvUuwOdB_0kK5kg7StgayuPf?usp=drive_link"
folder_path <- folder_download_gd(drive_link = f_link)
#> File downloaded:
#> • VCFannotationformat_v1.0.pdf <id: 1m36OwCrwxYYGC016z807Dz_AdNbyCr0I>
#> Saved locally as:
#> • /var/folders/n_/swy48fpx1w76xyqp3qx2prz00000gn/T//RtmphzeVRl/VCF Resources/VCFannotationformat_v1.0.pdf
#> File downloaded:
#> • VCFv4.2.pdf <id: 1fEuDbVPB5xvI0WYHM27YCvNwvah8BUQC>
#> Saved locally as:
#> • /var/folders/n_/swy48fpx1w76xyqp3qx2prz00000gn/T//RtmphzeVRl/VCF Resources/VCFv4.2.pdf
# }
```
