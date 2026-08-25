# Generate a Standardized Intertek Submission Table for KASP Markers

This function transforms marker data into a standardized data frame that
is compliant with the Intertek genotyping submission format. It
processes one or more markers, generates a unique Intertek-formatted SNP
name, and organizes all required information into a clean,
submission-ready table.

## Usage

``` r
make_intertek_table(
  marker_data,
  genome_version,
  region_name,
  trait,
  owner,
  publication_status = "Restricted"
)
```

## Arguments

- marker_data:

  A list containing a \`\$marker_data\` data frame. This data frame must
  have one or more rows, each with the following columns: \`SNP_Name\`,
  \`SNP\`, \`Marker_Name\`, \`Chromosome\`, \`Chromosome_Position\`,
  \`Sequence\`, \`ReferenceAllele\`, and \`AlternativeAllele\`.

- genome_version:

  A character string specifying the genome or assembly version (e.g.,
  "Sbv5.1").

- region_name:

  A character string for the QTL, gene, or region name associated with
  the markers.

- trait:

  A character string describing the trait associated with the markers.

- owner:

  A character string identifying the name of the person or lab
  responsible for the submission.

- publication_status:

  A character string for the publication status. Defaults to
  "Restricted".

## Value

A data frame with columns formatted for Intertek submission. If an error
occurs while processing a specific marker, the 'Comments' column for
that row will contain the error message.

## Details

The function constructs a unique \`SNP Name\*\` for each marker by
combining the genome version, chromosome, position, and a suffix
indicating the variant type. - \*\*Substitutions (SNPs):\*\* A suffix is
derived from the IUPAC ambiguity code for the allele pair (e.g., 'AG'
becomes 'R'). If the combination is not a standard IUPAC code, 'X' is
used. - \*\*Insertions/Deletions (INDELs):\*\* The suffix is always 'I'.

The final name format is \`genome_version_chr_pos_suffix\`.

## Examples

``` r
# \donttest{
marker_data <- list(
  marker_data = data.frame(
    SNP_Name = "SNP_Chr03_11361160",
    SNP = "Substitution",
    Marker_Name = "example",
    Chromosome = "Chr03",
    Chromosome_Position = 11361160,
    Sequence = "ACTG...[G/A]...CTGAAA",
    ReferenceAllele = "G",
    AlternativeAllele = "A",
    stringsAsFactors = FALSE
  )
)

make_intertek_table(marker_data = marker_data, genome_version = "Sbv5.1",
region_name = "qDT3.1", trait = "Stay-green", owner = "Cruet-Burgos")
#>   S/no            SNP Name* Alternative ID SNP*             Sequence*
#> 1      Sbv5.1_03_11361160_R                 G/A ACTG...[G/A]...CTGAAA
#>        Trait Gene/QTL Chromosome Number Chromosome Position Reference Allele
#> 1 Stay-green   qDT3.1                03            11361160                G
#>   Alternative Allele       Owner* Published/Restricted* Reference
#> 1                  A Cruet-Burgos            Restricted          
#>   SNP Geographical/Strain Relevance Comments
#> 1                                           
# }
```
