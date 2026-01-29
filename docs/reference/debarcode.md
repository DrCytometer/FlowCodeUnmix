# Debarcode

Debarcodes FlowCode data, identifying valid and invalid combinations.

## Usage

``` r
debarcode(unmixed, flowcode.fluors, flowcode.thresholds, combo.df)
```

## Arguments

- unmixed:

  Matrix of unmixed flow cytometry expression data. Cells in rows and
  fluorophores in columns. Fluorophores must include `flowcode.fluors`.

- flowcode.fluors:

  Named character vector of fluorophores used to identify the FlowCodes.
  Names should be FlowCode epitope tags, values should be the
  corresponding fluorophores.

- flowcode.thresholds:

  Named numeric vector of thresholds above which the cell (event) may be
  considered positive for a given FlowCode. Names should be
  fluorophores, values should be the threshold in raw, untransformed
  unmixed scale.

- combo.df:

  Data frame describing the valid combinations of FlowCodes. Structure:
  One row per combination. Columns are `Id`, `Procode.tag1`,
  `Procode.tag2` and `Procode.tag3`, describing the name (e.g., CRISPR
  target), and three epitopes for the combination, respectively.

## Value

A numeric vector, length `nrow(unmixed)`, linking each event to an `Id`
in `combo.df` by row number.
