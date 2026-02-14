# Extract putative causal variants within a candidate gene from a tabix-indexed snpEff annotated VCF file.

Extract putative causal variants within a candidate gene from a
tabix-indexed snpEff annotated VCF file.

## Usage

``` r
extract_variant(
  cand_gene_id,
  gff_path,
  vcf_dir,
  vcf_file,
  output_path = tempdir(),
  outfile_suffix = "variants"
)
```

## Arguments

- cand_gene_id:

  A character value specifying the candidate gene ID.

- gff_path:

  A character value indicating the path to the GFF file, including the
  complete file name.

- vcf_dir:

  A character value indicating the path to directory containing snpEff
  annotated VCF files.

- vcf_file:

  A character value indicating the file name for snpEff annotated VCF
  file including the .vcf.gz extension.

- output_path:

  A character value indicating the path to directory for saving
  extracted variants.

- outfile_suffix:

  A character value indicating the file name to be used for saving
  extracted variants.

## Value

A VCF file containing extracted variants for `cand_gene_id`.

## Details

This wrapper function operates on tabix-indexed snpEff annotated VCF
files. However, if a tabix-indexed VCF file is not available, it can
create one from the inputted VCF file.

The file names of snpEff annotated VCF files are expected to consist of
three components: a common prefix, chromosome tag and a common suffix.

## Examples

``` r
# example code
# \donttest{
library(panGenomeBreedr)

# Work from the tempdir
vcf_dir <- tempdir()

# Google drive link to gff3 file
flink1 <- "https://drive.google.com/file/d/1XjYyJ2JLywbbniIU6oUIIxAmEBKfmHpz/view?usp=sharing"

# Download gff3 file to tempdir()
gff3 <- folder_download_gd(drive_link = flink1,
                           output_path = vcf_dir,
                           is.folder = FALSE)
#> File downloaded:
#> • Sbicolor_730_v5.1.gene.gff3 <id: 1XjYyJ2JLywbbniIU6oUIIxAmEBKfmHpz>
#> Saved locally as:
#> • /var/folders/n_/swy48fpx1w76xyqp3qx2prz00000gn/T//RtmphzeVRl/Sbicolor_730_v5.1.gene.gff3

# Google drive link to indel snpEff annotated vcf file on Chr05
flink2 <- "https://drive.google.com/file/d/1LiOeDsfIwbsCuHbw9rCJ1FLOZqICTrfs/view?usp=sharing"

# Download indel snpEff annotated vcf file to tempdir()
vcf_file_indel <- folder_download_gd(drive_link = flink2,
                                     output_path = vcf_dir,
                                     is.folder = FALSE)
#> File downloaded:
#> • Sorghum_d8.noduplicates.Chr05.indel._markernamesadded_imputed_snpeff.vcf.gz
#>   <id: 1LiOeDsfIwbsCuHbw9rCJ1FLOZqICTrfs>
#> Saved locally as:
#> • /var/folders/n_/swy48fpx1w76xyqp3qx2prz00000gn/T//RtmphzeVRl/Sorghum_d8.noduplicates.Chr05.indel._markernamesadded_imputed_snpeff.vcf.gz
# View downloaded files in tempdir
# list.files(vcf_dir)

# InDel variant extraction for lgs1 (Sobic.005G213600)
extract_variant(cand_gene_id = 'Sobic.005G213600',
                gff_path = gff3,
                vcf_dir = vcf_dir,
                vcf_file = basename(vcf_file_indel),
                output_path = vcf_dir,
                outfile_suffix = 'lgs_variants_indel')

# Clean tempdir after variant extraction
# unlink(vcf_dir, recursive = TRUE)
# }
```
