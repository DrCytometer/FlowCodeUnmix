# FlowCode Combo Plot

Creates side-by-side pseudocolor density dot plots (standard flow
cytometry) for the FlowCode markers. Plots the unmixed data on a
biexponential scale using
[`AutoSpectral::create.biplot()`](https://drcytometer.github.io/AutoSpectral/reference/create.biplot.html).

## Usage

``` r
flowcode.combo.plot(
  unmixed,
  corrected,
  flowcode.fluors,
  asp,
  output.dir = "./flowcode_spectra",
  filename = "Corrected_biplots.jpg",
  width = 5,
  height = 10
)
```

## Arguments

- unmixed:

  Matrix of unmixed flow cytometry expression data. Cells in rows and
  fluorophores in columns. Fluorophores must include `flowcode.fluors`.

- corrected:

  Corrected matrix of unmixed flow cytometry expression data after
  processing for FRET errors.

- flowcode.fluors:

  Named character vector of fluorophores used to identify the FlowCodes.
  Names should be FlowCode epitope tags, values should be the
  corresponding fluorophores.

- asp:

  The AutoSpectral parameter list. Prepare using
  `get.autospectral.param`.

- output.dir:

  File path determining where the plot will be saved. Default is
  `./flowcode_spectra`.

- filename:

  Character, name of the saved plot. Default is `Corrected_biplots.jpg`.

- width:

  Numeric, default `10`. Width (in inches) of the saved plot.

- height:

  Numeric, default `5`. Height (in inches) of the saved plot.

## Value

None. Saves the plot to the specified file path in `output.dir`.
