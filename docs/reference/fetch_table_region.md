# Fetch data from database tables by genomic coordinates

This function queries variant, annotation, genotype, or metadata tables
based on genomic coordinates (chromosome, start, end) or a gene name. It
can operate in either 'local' mode (connecting to a local DuckDB
database) or 'online' mode (fetching data from a remote API endpoint).

## Usage

``` r
fetch_table_region(
  con = NULL,
  table_name = c("variants", "annotations", "genotypes"),
  chrom,
  start = NULL,
  end = NULL,
  gene_name = NULL,
  connect_db_mode = c("local", "online")
)
```

## Arguments

- con:

  A DBI connection object to the local database. Required only when
  \`connect_db_mode\` is 'local'. Defaults to \`NULL\`.

- table_name:

  A character value specifying the target view to query. Must be one of
  \`"variants"\`, \`"annotations"\`, \`"genotypes"\`, or \`"metadata"\`.

- chrom:

  A character value specifying the target chromosome name (e.g.,
  \`"Chr05"\`).

- start:

  Integer. Optional start coordinate for the target window region.

- end:

  Integer. Optional end coordinate for the target window region.

- gene_name:

  A character value indicating the specific Sobic gene model ID.
  Utilized explicitly when subsetting the \`"annotations"\` table.

- connect_db_mode:

  A character string specifying the connection mode. Can be either
  \`'local'\` (default) or \`'online'\`.

## Value

A data frame containing the targeted genomic records. For the
\`"genotypes"\` table, it returns individual samples unpacked directly
into intuitive wide columns, along with calculated allele metrics.

## Examples

``` r
# \donttest{
# Load the package
library(panGenomeBreedr)

# Locate the package example database folder
my_db_folder <- system.file("extdata", "pangenome_scale_db",
                           package = "panGenomeBreedr",
                           mustWork = TRUE)

# Establish a virtual connection to the offline database engine
con_local <- connect_local_db(folder_path = my_db_folder)
#> Successfully connected to the local offline database! Pangenome-scale database  mounted safely.

# Extract functional annotations inside a candidate locus region from local DB
local_annota_region <- fetch_table_region(
  con = con_local,
  table_name = "annotations",
  chrom = "Chr05",
  start = 75104537,
  end = 75106403,
  gene_name = "Sobic.005G213600"
)
print(head(local_annota_region))
#>             variant_id allele          annotation   impact        gene_name
#> 1 INDEL_Chr05_75104541   TGAC 3_prime_UTR_variant MODIFIER Sobic.005G213600
#> 2   SNP_Chr05_75104557      T 3_prime_UTR_variant MODIFIER Sobic.005G213600
#> 3   SNP_Chr05_75104560      T 3_prime_UTR_variant MODIFIER Sobic.005G213600
#> 4 INDEL_Chr05_75104564     CA 3_prime_UTR_variant MODIFIER Sobic.005G213600
#> 5   SNP_Chr05_75104568      T 3_prime_UTR_variant MODIFIER Sobic.005G213600
#> 6   SNP_Chr05_75104569      A 3_prime_UTR_variant MODIFIER Sobic.005G213600
#>                 gene_id feature_type              feature_id transcript_biotype
#> 1 Sobic.005G213600.v5.1   transcript Sobic.005G213600.1.v5.1     protein_coding
#> 2 Sobic.005G213600.v5.1   transcript Sobic.005G213600.1.v5.1     protein_coding
#> 3 Sobic.005G213600.v5.1   transcript Sobic.005G213600.1.v5.1     protein_coding
#> 4 Sobic.005G213600.v5.1   transcript Sobic.005G213600.1.v5.1     protein_coding
#> 5 Sobic.005G213600.v5.1   transcript Sobic.005G213600.1.v5.1     protein_coding
#> 6 Sobic.005G213600.v5.1   transcript Sobic.005G213600.1.v5.1     protein_coding
#>   rank            hgvs_c hgvs_p chrom      pos
#> 1  2/2 c.*322_*324dupGTC        Chr05 75104541
#> 2  2/2         c.*309G>A        Chr05 75104557
#> 3  2/2         c.*306G>A        Chr05 75104560
#> 4  2/2        c.*301dupT        Chr05 75104564
#> 5  2/2         c.*298C>A        Chr05 75104568
#> 6  2/2         c.*297G>T        Chr05 75104569

# Extract matrix genotypes within the exact same coordinates window from local DB
local_gt_region <- fetch_table_region(
  con = con_local,
  table_name = "genotypes",
  chrom = "Chr05",
  start = 75104537,
  end = 75106403
)
print(local_gt_region[1:6, 1:12])
#>             variant_id chrom      pos ref  alt variant_type major_allele
#> 1 INDEL_Chr05_75104541 Chr05 75104541   T TGAC        INDEL            T
#> 2   SNP_Chr05_75104557 Chr05 75104557   C    T          SNP            C
#> 3   SNP_Chr05_75104560 Chr05 75104560   C    T          SNP            C
#> 4 INDEL_Chr05_75104564 Chr05 75104564   C   CA        INDEL            C
#> 5   SNP_Chr05_75104568 Chr05 75104568   G    T          SNP            G
#> 6   SNP_Chr05_75104569 Chr05 75104569   C    A          SNP            C
#>   minor_allele major_allele_freq minor_allele_freq IDMM ISGC
#> 1         TGAC           0.99881           0.00119  0|0  0|0
#> 2            T           0.89051           0.10949  0|0  0|0
#> 3            T           0.88962           0.11038  0|0  0|0
#> 4           CA           0.88544           0.11456  0|0  0|0
#> 5            T           0.99791           0.00209  0|0  0|0
#> 6            A           0.99791           0.00209  0|0  0|0

# Disconnect at the end of your session
disconnect_local_db(con_local)
#> Successfully disconnected from the local database. Memory cleared.
# }

if (FALSE) { # \dontrun{
# To query the public online resource instead:
online_variants_data <- fetch_table_region(
  table_name = "variants",
  chrom = "Chr05",
  start = 75104537,
  end = 75106403,
  connect_db_mode = 'online'
)
print(head(online_variants_data))

online_genotypes_data <- fetch_table_region(
  table_name = "genotypes",
  chrom = "Chr05",
  start = 75104537,
  end = 75106403,
  connect_db_mode = 'online'
)
print(online_genotypes_data[1:6, 1:12])
} # }
```
