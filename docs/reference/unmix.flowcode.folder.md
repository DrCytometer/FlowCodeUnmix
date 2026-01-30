# Unmix All FCS Files in a Directory

This function unmixes all FCS files in a specified directory using the
provided spectra and method, and saves the unmixed FCS files to an
output directory of the user's choice.

## Usage

``` r
unmix.flowcode.folder(
  fcs.dir,
  spectra,
  asp,
  flow.control,
  af.spectra,
  spectra.variants,
  flowcode.spectra,
  thresholds.file = NULL,
  weights = NULL,
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

- fcs.dir:

  Directory containing FCS files to be unmixed.

- spectra:

  Matrix containing spectra information.

- asp:

  The AutoSpectral parameter list. Prepare using
  `get.autospectral.param`

- flow.control:

  A list containing flow cytometry control parameters.

- af.spectra:

  Spectral signatures of autofluorescences, normalized between 0 and 1,
  with fluorophores in rows and detectors in columns. Prepare using
  `get.af.spectra`. Required for `FlowCodeUnmixe` unmixing.

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

  Optional numeric vector of weights: one per fluorescent detector.
  Default is `NULL`, in which case weighting will be done by channel
  means.

- k:

  Numeric, controls the number of variants tested for each fluorophore,
  autofluorescence and FRET spectrum. Default is `10`. Values up to `10`
  provide additional benefit in unmixing quality, `1` will be fastest.

- output.dir:

  Directory to save the unmixed FCS files (default is
  `asp$unmixed.fcs.dir`).

- file.suffix:

  A character string to append to the output file name. Default is
  `NULL`

- include.imaging:

  Logical indicating whether to include imaging data in the written FCS
  file: relevant for S8 and A8. Default is `FALSE`

- parallel:

  Logical, default is `TRUE`, which enables parallel processing for
  per-cell unmixing methods.

- threads:

  Numeric, defaults to all available cores if `parallel=TRUE`.

- verbose:

  Logical, controls messaging. Default is `TRUE`.

## Value

None. Saves the unmixed FCS files to the specified output directory.
