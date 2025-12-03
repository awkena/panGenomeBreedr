# kasp_dat

A sample KASP genotype data for marker QC visualization.

## Usage

``` r
data(kasp_dat)
```

## Format

A data frame with 1344 observations and 11 variables:

- `DaughterPlate`:

  `character` KASP Daughter Plate ID.

- `MasterPlate`:

  `character` KASP Master Plate ID.

- `MasterWell`:

  `character` KASP Master Well ID.

- `Call`:

  `character` KASP observed genotype calls.

- `X`:

  `double` FAM fluorescence values.

- `Y`:

  `double` HEX fluorescence values.

- `SNPID`:

  `character` KASP SNP ID.

- `SubjectID`:

  `character` KASP Subject ID.

- `CustomerID`:

  `character` Customer unique ID for samples.

- `DaughterWell`:

  `character` KASP Daughter Well ID.

- `Group`:

  `character` Predicted genotype for positive controls.
