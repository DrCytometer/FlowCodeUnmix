# Launch the Thresholding Shiny App

A simple wrapper function to launch the manual threshold selection app.

## Usage

``` r
launch.threshold.app(
  flowcode.combo.file = NULL,
  flow.control = NULL,
  spectra = NULL
)
```

## Arguments

- flowcode.combo.file:

  Optional, default is `NULL`. Provide this when you wish to establish
  new thresholds using the app and will be working with an FCS file
  rather than loading the RDS object. When given: File name and path to
  the CSV file containing the information describing your FlowCode
  library. Describes the valid combinations of FlowCodes. Structure: One
  row per combination. Columns are `Id`, `Procode.tag1`, `Procode.tag2`
  and `Procode.tag3`, describing the name (e.g., CRISPR target), and
  three epitopes for the combination, respectively.

- flow.control:

  Optional, default is `NULL`. Provide this when you wish to establish
  new thresholds using the app and will be working with an FCS file
  rather than loading the RDS object. A list containing flow cytometry
  control parameters.

- spectra:

  Optional, default is `NULL`. Provide this when you wish to establish
  new thresholds using the app and will be working with an FCS file
  rather than loading the RDS object.Spectral signatures of
  fluorophores, normalized between 0 and 1, with fluorophores in rows
  and detectors in columns.
