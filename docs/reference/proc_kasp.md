# Process reshaped KASP genotype data for heatmap plotting

Process reshaped KASP genotype data for heatmap plotting

## Usage

``` r
proc_kasp(
  x,
  sample_id = "SubjectID",
  marker_start = 2,
  kasp_map,
  map_snp_id,
  chr = "chr",
  chr_pos = "pos"
)
```

## Arguments

- x:

  A data frame of genotype calls.

- sample_id:

  A string representing the column name of unique sample IDs.

- marker_start:

  An integer indicating the column index for the start of SNP calls.

- kasp_map:

  A data frame consisting of marker chromosome number and positions.

- map_snp_id:

  A character value indicating the column name for SNP IDs in
  `kasp_map`.

- chr:

  A character value for the column name for the chromosome number in
  `kasp_map`.

- chr_pos:

  A character value for the column name for the chromosome positions in
  `kasp_map`.

## Value

A data frame object of re-ordered SNPs based on chromosome numbers and
positions.

## Details

This function is experimental and should only be used with reshaped KASP
data. The data processing involves the following steps: 1. Transposing
the data 2. Matching marker names in KASP map file with names in KASP
genotype file 3. Column bind matched map for markers to the transposed
data 4. Sort data in ascending order of chromosome numbers 5. Split data
by chromosomes 6. Sort data in ascending order of physical positions per
chromosome 7. Row bind sorted data together 8. Re-transpose data to make
markers columns and samples as rows for plotting.

## Examples

``` r
# \donttest{
# example code
library(panGenomeBreedr)
dat1 <- panGenomeBreedr::beta_carotene
# Reshape KASP data for each master plate for the beta_carotene data
plate_wide <- kasp_reshape_wide(x = dat1,
                                subset = 'MasterPlate',
                                snp_id = 'SNPID',
                                geno_call = 'Call',
                                idvar = "SubjectID",
                                blank = 'NTC')

# Get Master Plate 1
plate1 <- plate_wide$`SE-24-1088_P01_d1`

# Generate a map for the beta_carotene KASP data
kasp_map <- data.frame(SNPID = unique(beta_carotene$SNPID),
                       SNPID_2 = paste0('snpSB', 1:4, '|', c(800, 803, 804, 805)),
                       chr = c(1:4),
                       pos = c(800, 803, 804, 805))

# Process Plate1 to re-order SNPs based on chrom. and position
proc_plate1 <- proc_kasp(x = plate1,
                         kasp_map = kasp_map,
                         map_snp_id = "SNPID",
                         sample_id = "SubjectID",
                         marker_start = 2,
                         chr = 'chr',
                         chr_pos = 'pos')
# }
```
