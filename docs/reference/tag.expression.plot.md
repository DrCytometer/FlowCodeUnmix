# Plot Tag Expression

Creates ridgeline (density) plots of the FlowCode fluorophore expression
levels from unmixed data.

## Usage

``` r
tag.expression.plot(
  unmixed,
  flowcode.fluors,
  flowcode.ids,
  combo.df,
  asp,
  output.dir = "./flowcode_spectra",
  title = "Uncorrected",
  max.per.group = 2000,
  height = 20,
  width = 20
)
```

## Arguments

- unmixed:

  Matrix of unmixed flow cytometry expression data. Cells in rows and
  fluorophores in columns. Fluorophores must include `flowcode.fluors`.

- flowcode.fluors:

  Named character vector of fluorophores used to identify the FlowCodes.
  Names should be FlowCode epitope tags, values should be the
  corresponding fluorophores.

- flowcode.ids:

  Numeric vector containing the FlowCode combination assignments per
  cell (result of
  [`debarcode()`](https://drcytometer.github.io/FlowCodeUnmix/reference/debarcode.md)).

- combo.df:

  Data frame describing the valid combinations of FlowCodes. Structure:
  One row per combination. Columns are `Id`, `Procode.tag1`,
  `Procode.tag2` and `Procode.tag3`, describing the name (e.g., CRISPR
  target), and three epitopes for the combination, respectively.

- asp:

  The AutoSpectral parameter list.

- output.dir:

  File path determining where the plot will be saved. Default is
  `./flowcode_spectra`.

- title:

  Character, name to add to the saved plots. Default is `Uncorrected`.

- max.per.group:

  Numeric, maximum number of cells to plot per group (FlowCode
  combination). Larger numbers will take longer. Default is `2000`.

- height:

  Numeric, default `20`. Width (in inches) of the saved plot.

- width:

  Numeric, default `20`. Height (in inches) of the saved plot.

## Value

None. Saves the plot to the specified file path in `output.dir`.
