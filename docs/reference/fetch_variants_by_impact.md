# Extract variants based on functional mutation impact

This function extracts variants based on their predicted functional
impact (e.g., HIGH, MODERATE, LOW, MODIFIER) within a specified genomic
region. It can operate in either 'local' mode (querying a local DuckDB
database) or 'online' mode (fetching data from a remote API endpoint).

## Usage

``` r
fetch_variants_by_impact(
  con = NULL,
  impact_level = c("HIGH", "MODERATE", "LOW", "MODIFIER"),
  chrom = NULL,
  start = NULL,
  end = NULL,
  connect_db_mode = c("local", "online")
)
```

## Arguments

- con:

  A DBI connection object to the local database. Required only when
  \`connect_db_mode\` is 'local'. Defaults to \`NULL\`.

- impact_level:

  A character vector specifying the variant impact types to retain.
  Allowed classifications are \`"HIGH"\`, \`"MODERATE"\`, \`"LOW"\`, and
  \`"MODIFIER"\`.

- chrom:

  A character value specifying the chromosome name (e.g., \`"Chr05"\`).

- start:

  Integer. Start coordinate for the genomic region.

- end:

  Integer. End coordinate for the genomic region.

- connect_db_mode:

  A character string specifying the connection mode. Can be either
  \`'local'\` (default) or \`'online'\`.

## Value

A data frame of variants showing their functional impact details and
positional coordinates.

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

# Extract high-impact variant lines inside a targeted coordinates locus window from local DB
local_high_impact_vars <- fetch_variants_by_impact(
  con = con_local,
  impact_level = "HIGH",
  chrom = "Chr05",
  start = 75104537,
  end = 75106403
)
print(head(local_high_impact_vars))
#>             variant_id chrom      pos   ref   alt qual filter variant_type
#> 1 INDEL_Chr05_75104881 Chr05 75104881     G GTCGA    .   PASS        INDEL
#> 2 INDEL_Chr05_75105587 Chr05 75105587    GC     G    .   PASS        INDEL
#> 3 INDEL_Chr05_75105598 Chr05 75105598    GT     G    .   PASS        INDEL
#> 4 INDEL_Chr05_75106156 Chr05 75106156 CGTAT     C    .   PASS        INDEL
#> 5 INDEL_Chr05_75106295 Chr05 75106295     A   ATC    .   PASS        INDEL
#> 6 INDEL_Chr05_75106325 Chr05 75106325     G   GTA    .   PASS        INDEL
#>   allele         annotation impact        gene_name               gene_id
#> 1  GTCGA frameshift_variant   HIGH Sobic.005G213600 Sobic.005G213600.v5.1
#> 2      G frameshift_variant   HIGH Sobic.005G213600 Sobic.005G213600.v5.1
#> 3      G frameshift_variant   HIGH Sobic.005G213600 Sobic.005G213600.v5.1
#> 4      C frameshift_variant   HIGH Sobic.005G213600 Sobic.005G213600.v5.1
#> 5    ATC frameshift_variant   HIGH Sobic.005G213600 Sobic.005G213600.v5.1
#> 6    GTA frameshift_variant   HIGH Sobic.005G213600 Sobic.005G213600.v5.1
#>   feature_type              feature_id transcript_biotype rank
#> 1   transcript Sobic.005G213600.1.v5.1     protein_coding  2/2
#> 2   transcript Sobic.005G213600.1.v5.1     protein_coding  2/2
#> 3   transcript Sobic.005G213600.1.v5.1     protein_coding  2/2
#> 4   transcript Sobic.005G213600.1.v5.1     protein_coding  2/2
#> 5   transcript Sobic.005G213600.1.v5.1     protein_coding  1/2
#> 6   transcript Sobic.005G213600.1.v5.1     protein_coding  1/2
#>               hgvs_c     hgvs_p
#> 1 c.1343_1344insTCGA p.Ser449fs
#> 2          c.637delG p.Ala213fs
#> 3          c.626delA p.Asp209fs
#> 4     c.65_68delATAC  p.His22fs
#> 5       c.38_39dupGA  p.Ser14fs
#> 6         c.8_9dupTA   p.Gln4fs

# Disconnect at the end of your session
disconnect_local_db(con_local)
#> Successfully disconnected from the local database. Memory cleared.
# }

if (FALSE) { # \dontrun{
# To query the public online resource instead:
online_high_impact_vars <- fetch_variants_by_impact(
  impact_level = "HIGH",
  chrom = "Chr05",
  start = 75104537,
  end = 75106403,
  connect_db_mode = 'online'
)
print(head(online_high_impact_vars))
} # }
```
