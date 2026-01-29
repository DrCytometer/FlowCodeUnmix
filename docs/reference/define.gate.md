# Define Gate

Defines a gate boundary on forward (FSC) and side (SSC) scatter
parameters. It does so by selecting the brightest events in the FlowCode
tag fluorescence channels, determining their location on scatter, and
defining a density/distance- based region around the centroid of those
events. The aim is to focus on the scatter region with the highest
signal, without any additional inputs or prior knowledge, as means of
determining where the cells are. This is important for identifying
positive events correctly, but also critical for matching the background
noise/autofluorescence/non-specific staining level on similar events
(cells rather than debris).

## Usage

``` r
define.gate(scatter.data, unmixed.data, flowcode.fluors, asp, n.cells = 200)
```

## Arguments

- scatter.data:

  Matrix of FSC and SSC data corresponding to the same events in
  `unmixed.data`. Used to define the scatter gating region.

- unmixed.data:

  Matrix of unmixed flow cytometry expression data. Cells in rows and
  fluorophores in columns. Fluorophores must include `flowcode.fluors`.

- flowcode.fluors:

  Named character vector of fluorophores used to identify the FlowCodes.
  Names should be FlowCode epitope tags, values should be the
  corresponding fluorophores.

- asp:

  The AutoSpectral parameter list.

- n.cells:

  Numeric, default is `200`. Number of brightest events per channel to
  use in defining the gate region.

## Value

Gate boundary
