# Unmix FlowCode

Unmix FlowCode samples, correcting FRET errors and debarcoding the data.

## Usage

``` r
unmix.flowcode.test(
  raw.data,
  spectra,
  af.spectra,
  spectra.variants,
  flowcode.spectra,
  asp,
  weights = NULL,
  cell.weighting = FALSE,
  cell.weight.regularize = FALSE,
  k = 1,
  parallel = TRUE,
  threads = if (parallel) 0 else 1,
  verbose = TRUE
)
```

## Arguments

- raw.data:

  Expression data from raw FCS files. Cells in rows and detectors in
  columns. Columns should be fluorescent data only and must match the
  columns in spectra.

- spectra:

  Spectral signatures of fluorophores, normalized between 0 and 1, with
  fluorophores in rows and detectors in columns.

- af.spectra:

  Spectral signatures of autofluorescences, normalized between 0 and 1,
  with fluorophores in rows and detectors in columns. Prepare using
  `get.af.spectra`.

- spectra.variants:

  Named list (names are fluorophores) carrying matrices of spectral
  signature variations for each fluorophore. Prepare using
  `get.spectral.variants()`. Default is `NULL`.

- flowcode.spectra:

  Structured output from
  [`get.flowcode.spectra()`](https://drcytometer.github.io/FlowCodeUnmix/reference/get.flowcode.spectra.md),
  which details the combination-level spectral unmixing errors due to
  FRET-like artefacts.

- asp:

  description

- weights:

  Optional numeric vector of weights (one per fluorescent detector).
  Default is `NULL`, in which case weighting will be done by channel
  means (Poisson variance).

- cell.weighting:

  Logical, whether to use cell-specific weighting for a more
  Poisson-like unmixing. Default is `FALSE`.

- cell.weight.regularize:

  Logical, whether to regularize cell-specific weights towards the bulk
  mean weighting set by `weights`. 50:50 averaging. Default is `FALSE`.

- k:

  Numeric, controls the number of variants tested for each fluorophore,
  autofluorescence and FRET spectrum. Default is `1`, which will be
  fastest. Values up to `10` provide additional benefit in unmixing
  quality.

- parallel:

  Logical, whether to use parallel processing for the per-cell unmixing.
  Default is `FALSE`.

- threads:

  Numeric. Number of threads to use for parallel processing. Defaults to
  `1` for sequential processing, or `0` (all cores) if `parallel=TRUE`.

- verbose:

  Logical, whether to send messages to the console. Default is `TRUE`.

## Value

Unmixed data with cells in rows and fluorophores in columns.
