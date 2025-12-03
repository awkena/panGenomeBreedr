# Design KASP markers based on causal variants.

Design KASP markers based on causal variants.

## Usage

``` r
kasp_marker_design(
  vcf_file = NULL,
  gt_df = NULL,
  variant_id_col = "variant_id",
  chrom_col = "chrom",
  pos_col = "pos",
  ref_al_col = "ref",
  alt_al_col = "alt",
  geno_start = 7,
  marker_ID,
  chr = NULL,
  genome_file,
  save_alignment = TRUE,
  plot_file = tempdir(),
  region_name = "loc_1",
  maf = 0.05
)
```

## Arguments

- vcf_file:

  Path to the vcf file containing identified variants.

- gt_df:

  A data frame or matrix containing the meta data of identified variants
  and sample VCF genotype calls. The variants are rows and samples as
  columns.

- variant_id_col, chrom_col, pos_col:

  A character value specifying the column names of variant IDs,
  chromosome, and positions in \`gt_df\` or \`vcf_file\`.

- ref_al_col, alt_al_col, :

  A character value specifying the column names of reference and
  alternate alleles, respectively in \`gt_df\` or \`vcf_file\`.

- geno_start:

  An integer value specifying the column index number of the start of
  the sample genotypes in \`gt_df\` or \`vcf_file\`.

- marker_ID:

  Designated name of variant for marker design. Name must be contained
  in \`gt_df\` or \`vcf_file\`.

- chr:

  A character value representation the chromosome description in the
  \`genome_file\`. Providing this helps to save memory in R.

- genome_file:

  Path to reference genome file in fasta format, either compressed (.gz)
  or decompressed.

- save_alignment:

  A logical value indicating whether to save the alignment plot to a
  directory specified by \`plot_file\`.

- plot_file:

  Path to save the sequence alignment if \`plot_draw = TRUE\`.

- region_name:

  A n optional character value for assigned region name.

- maf:

  A numeric value between 0 and 1 representing the minor allele
  frequency for variant subset.

## Value

A list object of two components comprising a data frame containing all
information required for KASP marker design, and a plot of DNA sequence
alignment to the reference genome.

## Details

This function provides the intertek sequence to be used for marker
development for the selected casual variants. It provides all the
information of allele and location to fill Intertek form.It needs a vcf
file of variants calls, and a genome sequence of the target crop in
fasta format.

## Examples

``` r
# \donttest{
# Example to design a KASP marker on a substitution variant
path <- tempdir()
path1 <- "https://raw.githubusercontent.com/awkena/panGB/main/Chr02.fa.gz"
path2 <-  system.file("extdata", "Sobic.002G302700_SNP_snpeff.vcf",
                     package = "panGenomeBreedr",
                     mustWork = TRUE)

ma1 <- kasp_marker_design(vcf_file = path2,
                          variant_id_col = 'ID',
                          chrom_col = 'CHROM',
                          pos_col = 'POS',
                          ref_al_col = 'REF',
                          alt_al_col = 'ALT',
                          genome_file = path1,
                          geno_start = 10,
                          marker_ID = "SNP_Chr02_69200443",
                          chr = "Chr02",
                          save_alignment = TRUE,
                          plot_file = path,
                          region_name = "ma1",
                          maf = 0.05)
#> using Gonnet

# View marker alignment output from temp folder
# path3 <- file.path(path, list.files(path = path, "alignment_"))
# system(paste0('open "', path3, '"'))

on.exit(unlink(path))
# }
```
