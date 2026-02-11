# Convert processed KASP data to numeric genotypes

Convert processed KASP data to numeric genotypes

## Usage

``` r
kasp_numeric(x, rp_row, dp_row, sep = ":", data_type = c("kasp", "agriplex"))
```

## Arguments

- x:

  A data frame of n rows of genotypes and p rows of SNPs; output from
  \`proc_kasp()\` function.

- rp_row:

  An integer or character value indicating the row index or name of the
  recurrent or Parent 1.

- dp_row:

  An integer or character value indicating the row index or name of
  donor or Parent 2.

- sep:

  A character used as separator for genotype calls, default is a colon.

- data_type:

  A character value indicating the data source; either \`kasp\` or
  \`agriplex\`.

## Value

A data frame of numeric codes for KASP genotype calls.

## Details

Re-coded as 1 if homozygous for Parent 1 allele; 0 if homozygous for
Parent 2 allele, and 0.5 if heterozygous. If any of the parents SNP call
is missing, it's coded as as -5. If the Parents 1 and 2 genotypes are
the same for any snp, it's coded as monomorphic. Loci with progeny
genotype calls different from RP and DP are coded as -2.

## Examples

``` r
# \donttest{
 # example code
library(panGenomeBreedr)

# Reshape KASP data for each master plate for the beta_carotene data
dat1 <- panGenomeBreedr::beta_carotene
plate_wide <- kasp_reshape_wide(x = dat1,
                                subset = 'MasterPlate',
                                snp_id = 'SNPID',
                                geno_call = 'Call',
                                idvar = "SubjectID",
                                blank = 'NTC')

# Get Master Plate 1
plate1 <- plate_wide$`SE-24-1088_P01_d1`

# Convert to numeric format for plotting
num_geno <- kasp_numeric(x = plate1[,-1],
                        rp_row = 1,
                        dp_row = 7,
                        data_type = 'kasp')

# }
```
