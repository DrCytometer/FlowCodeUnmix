# Unmix Backbone

Unmix the backbone sample, stained only with anti-FlowCode epitope
antibodies, for downstream identification of FRET/unmixing errors on a
per-combo basis.

## Usage

``` r
unmix.backbone(
  flowcode.backbone.fcs,
  spectra,
  af.spectra,
  flow.control,
  asp,
  flowcode.combo.file,
  n.cells.gate = 100,
  output.dir = "./flowcode_spectra",
  filename = "FlowCode_Backbone.rds",
  verbose = TRUE
)
```

## Arguments

- flowcode.backbone.fcs:

  File name and path for the FlowCode backbone control FCS file. This
  should be a sample of cells, ideally the same cell source as your
  single-stained control samples. These cells should be stained with all
  the FlowCode epitope tag antibodies present in the fully stained
  sample and nothing else. All tag combinations should be present and
  well represented for best results. A minimum of 100 cells per
  combination is recommended; 2000+ is ideal.

- spectra:

  Spectral signatures of fluorophores, normalized between 0 and 1, with
  fluorophores in rows and detectors in columns.

- af.spectra:

  Spectral signatures of autofluorescences, normalized between 0 and 1,
  with fluorophores in rows and detectors in columns. Prepare using
  `get.af.spectra`.

- flow.control:

  A list containing flow cytometry control parameters.

- asp:

  The AutoSpectral parameter list. Prepare using
  `get.autospectral.param`.

- flowcode.combo.file:

  File name and path to the CSV file containing the information
  describing your FlowCode library. Describes the valid combinations of
  FlowCodes. Structure: One row per combination. Columns are `Id`,
  `Procode.tag1`, `Procode.tag2` and `Procode.tag3`, describing the name
  (e.g., CRISPR target), and three epitopes for the combination,
  respectively.

- n.cells.gate:

  Numeric, default `100`. The number of cells that will be used to
  define the gate. A gate region will be defined using
  [`define.gate()`](https://drcytometer.github.io/FlowCodeUnmix/reference/define.gate.md)
  around the brightest `n.cells.gate` for each fluorophore used to
  identify FlowCode epitope tags. Scatter coordinates of all brightest
  events will be combined to determine the optimal gating region.

- output.dir:

  Destination folder for output plots and saved spectra. Default is
  `./flowcode_spectra`.

- filename:

  Name of the output RDS file, default is `FlowCode_Backbone.rds`.

- verbose:

  Logical, whether to send messages to the console. Default is `TRUE`.

## Value

None. Saves RDS objects containing the processed data to disk.
