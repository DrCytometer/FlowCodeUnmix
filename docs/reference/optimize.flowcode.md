# Optimize FlowCode spectral unmixing

Parallel backend for FlowCode-based spectral optimization in R.

## Usage

``` r
optimize.flowcode(
  raw_data,
  unmixed,
  unmix_fun,
  combined_spectra,
  weights,
  pos_thresholds,
  af_idx,
  af_spectra,
  flowcode_ids,
  has_flowcode,
  combo_fret,
  fret_delta_list,
  fret_delta_norms,
  flowcode_combo_logical,
  flowcode_fluors,
  optimize_fluors,
  variants,
  delta_list,
  delta_norms,
  fluorophores,
  af_idx_in_spectra,
  asp,
  k = 10L,
  cell_weighting = FALSE,
  cell_weight_regularize = FALSE,
  nthreads = 1L,
  parallel = TRUE
)
```

## Arguments

- raw_data:

  Numeric matrix (n_cells x n_detectors)

- unmixed:

  Numeric matrix (n_cells x n_fluors)

- unmix_fun:

  Unmixing function (WLS)

- combined_spectra:

  Numeric matrix (n_fluors x n_detectors)

- weights:

  Numeric vector (n_detectors)

- pos_thresholds:

  Numeric vector (n_fluors)

- af_idx:

  Integer vector (n_cells), 1-based AF indices

- af_spectra:

  Numeric matrix (n_af_variants x n_detectors)

- flowcode_ids:

  Integer vector (n_cells), 0 = no flowcode

- has_flowcode:

  Integer/logical vector (n_cells)

- combo_fret:

  List of FRET spectra matrices

- fret_delta_list:

  List of delta matrices

- fret_delta_norms:

  List of norm vectors

- flowcode_combo_logical:

  Integer matrix (n_combos x n_flowcode_fluors)

- flowcode_fluors:

  Character vector

- optimize_fluors:

  Character vector

- variants:

  List of variant matrices per fluorophore

- delta_list:

  List of delta matrices per fluorophore

- delta_norms:

  List of delta norms per fluorophore

- fluorophores:

  Character vector of fluorophore names

- af_idx_in_spectra:

  Integer, 0-based index of AF in spectra

- asp:

  The AutoSpectral parameter list.

- k:

  Integer, number of variants to test

- cell_weighting:

  Logical, use cell-specific weights

- cell_weight_regularize:

  Logical, regularize cell weights

- nthreads:

  Integer, number of threads

- parallel:

  Logical, whether to use parallel processing

## Value

Unmixed data with cells in rows and fluorophores in columns.
