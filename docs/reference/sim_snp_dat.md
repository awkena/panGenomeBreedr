# Simulate raw SNP loci for any chromosome with or without LD.

Simulate raw SNP loci for any chromosome with or without LD.

## Usage

``` r
sim_snp_dat(
  nsnp = 10L,
  nobs = 100L,
  chr = 1L,
  start = 1000L,
  end = 20000L,
  sep = "/",
  add_LD = FALSE,
  LD_range = NULL
)
```

## Arguments

- nsnp:

  An integer value specifying the number of SNPs to simulate.

- nobs:

  An integer value specifying the number of individuals to simulate.

- chr:

  An integer value specifying the chromosome number.

- start:

  An integer value specifying the position of the first SNP.

- end:

  An integer value specifying the position of the last SNP.

- sep:

  A separator for deriving the genotypes.

- add_LD:

  A logical value indicating whether to simulate SNPs in LD or not.

- LD_range:

  A numeric vector of \`length = 2\` indicating the range of LD between
  SNPs, if `add_LD` = TRUE.

## Value

A data frame with the simulated SNPs as the columns and individuals as
rows.

## Examples

``` r
# example code
# Simulate data for 20 snp and 100 individuals on chr 1
geno_data <- sim_snp_dat(nsnp = 20,
                     nobs = 100,
                     chr = 1,
                     start = 1000,
                     end = 20000,
                     add_LD = TRUE,
                     LD_range = c(0.2, 1))
```
