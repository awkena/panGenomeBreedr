% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/variant_annotation.R
\name{extract_variant}
\alias{extract_variant}
\title{Extract putative causal variants within a candidate gene from a tabix-indexed
snpEff annotated VCF file.}
\usage{
extract_variant(
  cand_gene_id,
  gff_path,
  vcf_dir,
  vcf_file,
  output_path = tempdir(),
  outfile_suffix = "variants"
)
}
\arguments{
\item{cand_gene_id}{A character value specifying the candidate gene ID.}

\item{gff_path}{A character value indicating the path to the GFF file,
including the complete file name.}

\item{vcf_dir}{A character value indicating the path to directory containing
snpEff annotated VCF files.}

\item{vcf_file}{A character value indicating the file name for snpEff
annotated VCF file including the .vcf.gz extension.}

\item{output_path}{A character value indicating the path to directory for
saving extracted variants.}

\item{outfile_suffix}{A character value indicating the file name to be used for
saving extracted variants.}
}
\value{
A VCF file containing extracted variants for \code{cand_gene_id}.
}
\description{
Extract putative causal variants within a candidate gene from a tabix-indexed
snpEff annotated VCF file.
}
\details{
This wrapper function operates on tabix-indexed snpEff annotated VCF files.
However, if a tabix-indexed VCF file is not available, it can create one from
the inputted VCF file.

The file names of snpEff annotated VCF files are expected to consist of three
components: a common prefix, chromosome tag and a common suffix.
}
\examples{
# example code
\donttest{
library(panGenomeBreedr)

# Work from the tempdir
vcf_dir <- tempdir()

# Google drive link to gff3 file
flink1 <- "https://drive.google.com/file/d/1XjYyJ2JLywbbniIU6oUIIxAmEBKfmHpz/view?usp=sharing"

# Download gff3 file to tempdir()
gff3 <- folder_download_gd(drive_link = flink1,
                           output_path = vcf_dir,
                           is.folder = FALSE)

# Google drive link to indel snpEff annotated vcf file on Chr05
flink2 <- "https://drive.google.com/file/d/1LiOeDsfIwbsCuHbw9rCJ1FLOZqICTrfs/view?usp=sharing"

# Download indel snpEff annotated vcf file to tempdir()
vcf_file_indel <- folder_download_gd(drive_link = flink2,
                                     output_path = vcf_dir,
                                     is.folder = FALSE)
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
}

}
