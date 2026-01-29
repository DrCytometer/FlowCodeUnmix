# Get FlowCode Spectra

Measure spectral unmixing errors per barcode in a sample stained only
with FlowCode tag antibodies (`flowcode.backbone.fcs`). Measures FRET
error and returns the calculated correction matrices and thresholds for
debarcoding.

## Usage

``` r
get.flowcode.spectra(
  backbone.rds,
  thresholds.file,
  asp,
  output.dir = "./flowcode_spectra",
  filename = "FlowCode_Spectra.rds",
  figures = TRUE,
  plot.corrections = TRUE,
  verbose = TRUE
)
```

## Arguments

- backbone.rds:

  File name and path to the full (large) rds file produced by running
  `unmix.backbone`. For this to work well, all combinations present in
  the experimental samples should be well represented in the backbone
  sample used.

- thresholds.file:

  File name and path to the thresholds CSV file produced using the
  threshold setting app.

- asp:

  The AutoSpectral parameter list. Prepare using
  `get.autospectral.param`.

- output.dir:

  Destination folder for output plots and saved spectra. Default is
  `./flowcode_spectra`.

- filename:

  Logical, whether to create and save plots of the data correction.
  Default is `TRUE`.

- figures:

  Logical, controls whether graphs of the process are produced. Default
  is `TRUE`.

- plot.corrections:

  Logical, default is `FALSE`. If `TRUE`, the FRET unmixing will be
  applied to `flowcode.backbone.fcs` and plots of the corrected data
  (FlowCode channels only) will be produced in `output.dir`.

- verbose:

  Logical, default is `TRUE`. Set to `FALSE` to suppress messaging.

## Value

A named list with two elements:

1.  Thresholds: numeric thresholds for debarcoding

2.  FRET: A list of matrices of FRET spectral variations per barcode
    combo

3.  Flowcode.fluors: Character vector of fluorophores linked to the
    FlowCodes

4.  Combos: The valid combinations, specified by the combination file.
