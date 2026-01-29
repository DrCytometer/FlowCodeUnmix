# Unmix FlowCode FCS Data

This function performs spectral unmixing on FCS data from FlowCode
samples.

## Usage

``` r
unmix.flowcode.fcs(
  fcs.file,
  spectra,
  asp,
  flow.control,
  af.spectra,
  spectra.variants,
  flowcode.spectra,
  thresholds.file = NULL,
  weights = NULL,
  cell.weighting = FALSE,
  cell.weight.regularize = TRUE,
  k = 10,
  output.dir = NULL,
  file.suffix = NULL,
  include.imaging = FALSE,
  parallel = TRUE,
  threads = if (parallel) 0 else 1,
  verbose = TRUE
)
```

## Arguments

- fcs.file:

  A character string specifying the path to the FCS file.

- spectra:

  A matrix containing the spectral data.

- asp:

  The AutoSpectral parameter list. Prepare using
  `get.autospectral.param`

- flow.control:

  A list containing flow cytometry control parameters.

- af.spectra:

  Spectral signatures of autofluorescences, normalized between 0 and 1,
  with fluorophores in rows and detectors in columns. Prepare using
  `get.af.spectra`. Required for `AutoSpectral` unmixing. Default is
  `NULL` and will thus provoke failure if no spectra are provided and
  `AutoSpectral` is selected.

- spectra.variants:

  Named list (names are fluorophores) carrying matrices of spectral
  signature variations for each fluorophore. Prepare using
  `get.spectral.variants`.

- flowcode.spectra:

  Structured output from
  [`get.flowcode.spectra()`](https://drcytometer.github.io/FlowCodeUnmix/reference/get.flowcode.spectra.md),
  which details the combination-level spectral unmixing errors due to
  FRET-like artefacts.

- thresholds.file:

  Optional, file name and path to the thresholds CSV file produced using
  the threshold setting app. Thresholds are provided by default as part
  of `flowcode.spectra`, but those will have been selected on the
  FlowCode backbone sample. If the thresholds for positivity of any of
  the fluorophores used to identify the FlowCodes differs in the fully
  stained samples being unmixed here, provide a new set of thresholds
  from the ThresholdApp. Run this on an unmixed sample (OLS from the
  instrument is fine), and refer to the new thresholds file here.
  Default is `NULL`, which will be ignored.

- weights:

  Optional numeric vector of weights (one per fluorescent detector).
  Default is `NULL`, in which case weighting will be done by channel
  means (Poisson variance).

- cell.weighting:

  Logical, whether to use cell-specific weighting for a more
  Poisson-like unmixing. Default is `FALSE`.

- cell.weight.regularize:

  Logical, whether to regularize cell-specific weights towards the bulk
  mean weighting set by `weights`. 50:50 averaging. Default is `TRUE`.
  Only active if `cell.weighting=TRUE`.

- k:

  Numeric, controls the number of variants tested for each fluorophore,
  autofluorescence and FRET spectrum. Default is `10`. Values up to `10`
  provide additional benefit in unmixing quality, `1` will be fastest.

- output.dir:

  A character string specifying the directory to save the unmixed FCS
  file. Default is `NULL`, which reverts to `asp$unmixed.fcs.dir`.

- file.suffix:

  A character string to append to the output file name. Default is
  `NULL`.

- include.imaging:

  A logical value indicating whether to include imaging parameters in
  the written FCS file. Default is `FALSE`.

- parallel:

  Logical, default is `TRUE`, which enables parallel processing for
  per-cell unmixing methods.

- threads:

  Numeric, default is `NULL`, in which case `asp$worker.process.n` will
  be used. `asp$worker.process.n` is set by default to be one less than
  the available cores on the machine. Multi-threading is only used if
  `parallel` is `TRUE`.

- verbose:

  Logical, whether to send messages to the console. Default is `TRUE`.

## Value

None. The function writes the unmixed FCS data to a file.
