# beta_carotene

A sample trait-predictive KASP genotype data for marker QC
visualization.

## Usage

``` r
data(beta_carotene)
```

## Format

A data frame with 768 observations and 11 variables:

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

- `DaughterWell`:

  `character` KASP Daughter Well ID.

- `Group`:

  `character` Predicted genotype for positive controls.

- `plates`:

  `character` Derived plate IDs for each 96-well KASP reaction plate.
