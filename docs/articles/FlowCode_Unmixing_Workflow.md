# FlowCode Spectral Flow Unmixing Workflow

## Installation

You will both need `AutoSpectral` and `FlowCodeUnmix` to run this. If
this is your first time using these, run the installation code below,
which should take care of everything for you.

``` r
# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("flowWorkspace", "flowCore", "PeacoQC", "FlowSOM"))

# Install any missing CRAN packages that we'll need
cran.packages <- c(
  "shiny",
  "shinyFiles",
  "ggplot2",
  "scattermore",
  "dplyr"
)

missing.pkgs <- setdiff(
  cran.packages,
  rownames(installed.packages())
)

if (length(missing.pkgs) > 0) {
  install.packages(missing.pkgs, dependencies = TRUE)
}

# install `remotes` to allow you to install GitHub packages (devtools also works)
install.packages("remotes")
remotes::install_github("DrCytometer/AutoSpectral")
remotes::install_github("DrCytometer/FlowCodeUnmix")
```

## Basic AutoSpectral workflow

We start as we would with any spectral flow cytometry data set in
AutoSpectral. We will load the cytometer-specific parameters, load in
the single-stained control FCS data, clean-up the control data, extract
the spectral profiles of the fluorophores, identify autofluorescence
signatures and measure the spectral variation in the fluorophore output.

### Load the parameters

Here we are using data from the Cytek Aurora (5-laser system). For data
from a different spectral cytometer, call the appropriate cytometer as
“cytometer”. For more details on the AutoSpectral workflow, see
[AutoSpectral help](https://drcytometer.github.io/AutoSpectral/).

``` r
library(AutoSpectral)
autospectral.parameters <- AutoSpectral::get.autospectral.param(
  cytometer = "aurora",
  figure = TRUE
)
```

### Create and check the control file

The control file is a CSV file (spreadsheet) that describes your control
files for AutoSpectral. A draft of this will be created automatically by
AutoSpectral if you tell it where your single-stained control files are.
You will need to edit this and check it to be sure it is correct.
Details on how to do this are in this
[article](https://drcytometer.github.io/AutoSpectral/articles/02_Control_File_example.html).
There is also a [Shiny
helper](https://github.com/DrCytometer/AutoSpectralHelper) application
for this, which opens as an interactive webpage.

``` r
# where are your single-stained control files?
# replace the "my_control_folder" with the path to your files.
# Tip: if you start RStudio by double-clicking on this or another R document,
# that tells RStudio where you want to work. If you place a folder (say, called SSC)
# containing your control fcs files in that folder, you will be able to access 
# it with just "./SSC
control.folder <- "C://User/Documents/my_control_folder"

# or
control.folder <- "./SSC"

AutoSpectral::create.control.file(
  control.dir = control.folder,
  asp = autospectral.parameters
)
```

Once you have edited your control file, check it:

``` r
control.file <- "fcs_control_file.csv"
AutoSpectral::check.control.file(
  control.dir = control.folder,
  control.def.file = control.file,
  asp = autospectral.parameters
)
```

### Load the data

This step will likely take a few minutes. AutoSpectral will now take all
of the information in the control file, read in all the single-stained
FCS files, determine gates to use and extract the information in the
gates for downstream usage. This can go faster by setting
`parallel = TRUE`.

``` r
flow.control <- AutoSpectral::define.flow.control(
  control.dir = control.folder,
  control.def.file = control.file,
  asp = autospectral.parameters,
  parallel = FALSE
)
```

Check the plots in folder `figure_gate` to be sure the gating is working
appropriately for your samples. See the help pages if it does not look
good.

### Clean the controls

This will also take a few minutes, most likely. In this step, we are
trying to identify the most useful events in the single-stained controls
for determining the fluorophore spectra. This involves removing any
cells (in cell-based controls) that interfere with the measurement
(e.g., autofluorescent macrophages or dead cells), selecting the
brightest events and matching the background of those events in the
unstained sample.

``` r
cleaned.controls <- AutoSpectral::clean.controls(
  flow.control = flow.control,
  asp = autospectral.parameters
)
```

### Measure the fluorophore spectral profiles

Now we will extract the fluorophore spectral information from the
cleaned data. This should be relatively quick. Plots will appear in
`figure_spectra` and `figure_similarity_heatmap`.

``` r
fluorophore.spectra <- AutoSpectral::get.fluorophore.spectra(
  flow.control = cleaned.controls,
  asp = autospectral.parameters,
  use.clean.expr = TRUE
)
```

### Measure autofluorescence spectra

Autofluorescence is always present in our stained cells, including our
single-stained cell controls. This is a source of variability that
contributes to noise (spread) in our unmixed data. We need to identify
the variability in autofluorescence to counteract this. With FlowCode
data, we will often have different tissue sources, each of which will
have a distinct distribution of autofluorescence profiles. We will
eventually need to identify all of those. To start, though, we will
identify only the variability in the type of sample we are using for the
single-stained controls. For instance, with mouse studies, these will
likely be splenocytes, so we will use unstained spleen as the source of
our autofluorescence.

In this example, I have an unstained spleen sample as the “Unstained”
(so named by SpectroFlo) in my reference control group. This FCS file is
in the single-stained control folder, along with all of the other
controls.

``` r
spleen.autofluorescence <- AutoSpectral::get.af.spectra(
  unstained.sample = file.path( control.folder, "Unstained.fcs" ),
  asp = autospectral.parameters,
  spectra = fluorophore.spectra
)
```

### Determine spectral variation in the fluorophores

Now that we know the contribution of autofluorescence to the variation,
we can extract the variation from the fluorophore emissions themselves.

``` r
fluorophore.variation <- AutoSpectral::get.spectral.variants(
  control.dir = control.folder,
  control.def.file = control.file,
  asp = autospectral.parameters,
  spectra = fluorophore.spectra,
  af.spectra = spleen.autofluorescence,
  parallel = FALSE
)
```

That is the end of the normal AutoSpectral workflow for us. We now need
to do some FlowCode-specific steps before proceeding to the debarcoding
and unmixing.

## FlowCode Unmixing Workflow

For the FlowCode workflow, we need some additional pieces of
information: \* Which are the FlowCode epitope tag channels? \* What are
the valid FlowCode combinations (in the combination file)? \* What is
the name (guide/TCR) of each combination (in the combination file)? \*
Where are the cells (user-drawn gate)? \* What variation in fluorophore
output is created by the combinations due to FRET or a FRET-like
artifact? For this, we need your backbone control, stained only with the
FlowCode epitope antibodies. All combinations should be well represented
for best results.

### Where is everything?

Provide the location (file path) and name of your backbone FCS file.
Provide the name (and file path) of your combination CSV file.

``` r
flowcode.backbone <- "./Raw/FC control (Spleen)/F2 FC_FRET_02_Plate_002.fcs" 
combo.file <- "OTB01 Treg Tx factor combinations.csv"
```

Now we can unmix the backbone control, using a slight adaptation of the
AutoSpectral unmixing approach. This is a spleen-based sample, so we use
`spleen.autofluorescence`.

``` r
FlowCodeUnmix::unmix.backbone(
   flowcode.backbone.fcs = flowcode.backbone,
   spectra = fluorophore.spectra,
   af.spectra = spleen.autofluorescence,
   flow.control = cleaned.controls,
   asp = autospectral.parameters,
   flowcode.combo.file = combo.file,
   output.dir = "./flowcode_spectra",
   filename = "FlowCode_Backbone.rds"
)
```

This saves the output as an RDS file in the folder specified by
`output.dir`. If you change `output.dir`, you’ll need to change it later
in
[`get.flowcode.spectra()`](https://drcytometer.github.io/FlowCodeUnmix/reference/get.flowcode.spectra.md).
\*Note: change to returning object

### Setting thresholds

This part requires you to manually select the point on the graph where
the data start becoming “positive” for each of the FlowCode epitope
tags. That is, where is the threshold for having the FlowCode or not?
This is manual because it’s pretty hard to do it well in an automated
manner considering that at this stage we still have FRET-based unmixing
issues.

I have created a Shiny app to allow you to do this interactively, as in
Orian Bricard’s FlowCode debarcoding app.

``` r
FlowCodeUnmix::launch.threshold.app()
```

When you’re done, you should have a CSV file (spreadsheet), which is
pretty simple and just contains the numerical value corresponding to the
threshold for each FlowCode channel.

### Measuring FRET (unmixing) errors due to FlowCodes

In this step, we will take the FlowCode epitope-stained backbone data
(from the unmixed data in the RDS file), debarcode it to identify valid
combinations, and assess spectral unmixing errors per combination. This
gives us a set of potential corrections to apply to each cell for a
given barcode combination.

``` r
flowcode.spectra <- FlowCodeUnmix::get.flowcode.spectra(
  backbone.rds = "./flowcode_spectra/FlowCode_Backbone.rds",
  thresholds.file = "./flowcode_spectra/thresholds.csv",
  output.dir = "./flowcode_spectra",
  filename = "FlowCode_Spectra.rds",
  plot.corrections = TRUE
)
```

This call saves the output as an RDS file in the folder specified by
`output.dir`. If you select `plot.corrections = TRUE`, it will identify
the FRET errors on a cell-by-cell basis, correct them, and plot examples
of what the FlowCode epitope unmixing looks like on the backbone file
with and without corrections. That can be a bit slow, but is definitely
worth doing the first couple of times.

### FlowCode Unmixing

We now have all the data required for unmixing, or, at least we do for
the any files coming from spleen. We’ll get to other tissue sources in a
second.

To unmix a single file, call
[`unmix.flowcode.fcs()`](https://drcytometer.github.io/FlowCodeUnmix/reference/unmix.flowcode.fcs.md).

``` r
# for example:
fully.stained.spleen <- "./Raw/Spleen/Spleen_Mouse_001.fcs"

# then call:
FlowCodeUnmix::unmix.flowcode.fcs(
  fcs.file = fully.stained.spleen,
  spectra = fluorophore.spectra,
  asp = autospectral.parameters,
  flow.control = cleaned.controls,
  af.spectra = spleen.autofluorescence,
  spectra.variants = fluorophore.variation,
  flowcode.spectra = flowcode.spectra,
  weighted = TRUE, # for weighted least squares unmixing
  k = 1, # use k=1 for fastest unmixing, k = 10 for slower but better, k ~ 3 compromise
  parallel = TRUE
)
```

To unmix a folder containing many files, call `unmix.flowcode.folder()`.

``` r
# for example:
spleen.folder <- "./Raw/Spleen/"

# then call:
FlowCodeUnmix::unmix.flowcode.folder(
  sample.dir = spleen.folder,
  spectra = fluorophore.spectra,
  asp = autospectral.parameters,
  flow.control = cleaned.controls,
  af.spectra = spleen.autofluorescence,
  spectra.variants = fluorophore.variation,
  flowcode.spectra = "./flowcode_spectra/FlowCode_Spectra.rds",
  weighted = TRUE, # for weighted least squares unmixing (default)
  k = 1, # use k=1 (default) for fastest unmixing, k = 10 for slower but better, k ~ 3 compromise
  parallel = TRUE
)
```

If we want to look at other tissues sources, or any samples where the
autofluorescence may differ from the autofluorescence profiles we have
already extracted, we first need to call
[`AutoSpectral::get.af.spectra()`](https://drcytometer.github.io/AutoSpectral/reference/get.af.spectra.html)
on unstained samples representing those autofluorescence sources. We can
then unmix samples from those tissues, extracting the matching
autofluorescence correctly.

``` r
# Let's say we have brain samples

# this is our unstained brain sample (raw data)
unstained.brain <- "./Raw/Unstained controls/Unstained brain.fcs"

# we extract the AF profiles from it
brain.autofluorescence <- AutoSpectral::get.af.spectra(
  unstained.sample = unstained.brain,
  asp = autospectral.parameters,
  spectra = fluorophore.spectra
)

# now we want the raw data for the fully stained brain sample
fully.stained.brain <- "./Raw/Brain/Brain_Mouse_001.fcs"

# then call to unmix:
FlowCodeUnmix::unmix.flowcode.fcs(
  fcs.file = fully.stained.brain, # change the sample
  spectra = fluorophore.spectra,
  asp = autospectral.parameters,
  flow.control = cleaned.controls,
  af.spectra = brain.autofluorescence, # change the autofluorescence to match
  spectra.variants = fluorophore.variation,
  flowcode.spectra = "./flowcode_spectra/FlowCode_Spectra.rds",
  parallel = TRUE
)
```

In many cases, the thresholds for positivity may vary slightly by tissue
type or sample source. If that is the case, you will need to provide new
thresholds as part of the unmixing process. To start, unmix a single FCS
file from the source you want to use to set new thresholds. Do the
unmixing with AutoSpectral, extracting the AF per cell, since this will
get you closer to what the data will look like with the FlowCodeUnmix.
Then load this file into the ThresholdApp, select new thresholds for
each FlowCode marker, and save the results (you can and should rename
the CSV file). Now provide this new CSV file containing the updated
thresholds to the unmixing.

I suspect that in most cases there will be little need for
re-thresholding between different tissues. There may be a need to
re-threshold on a fully stained sample, though, since this will differ
from the backbone control we have used to establish the first set of
thresholds.

``` r

# this is our unstained brain sample (raw data)
unstained.brain <- "./Raw/Unstained controls/Unstained brain.fcs"

# we extract the AF profiles from it
brain.autofluorescence <- AutoSpectral::get.af.spectra(
  unstained.sample = unstained.brain,
  asp = autospectral.parameters,
  spectra = fluorophore.spectra
)

# now we want the raw data for the fully stained brain sample
fully.stained.brain <- "./Raw/Brain/Brain_Mouse_001.fcs"


launch.threshold.app(
  flowcode.combo.file = combo.file,
  flow.control = cleaned.controls,
  spectra = fluorophore.spectra
)

new.thresholds <- "./flowcode_spectra/New_thresholds.csv"

# then call to unmix:
FlowCodeUnmix::unmix.flowcode.fcs(
  fcs.file = fully.stained.brain, # change the sample
  spectra = fluorophore.spectra,
  asp = autospectral.parameters,
  flow.control = cleaned.controls,
  af.spectra = brain.autofluorescence, # change the autofluorescence to match
  spectra.variants = fluorophore.variation,
  flowcode.spectra = "./flowcode_spectra/FlowCode_Spectra.rds",
  thresholds.file = new.thresholds, # specify new thresholds
  parallel = TRUE
)
```

## FlowCode Debarcoding Workflow

With the new unmixing, the resulting FCS files are already debarcoded.
If you inspect the files in R or FlowJo, you should see a channel for
each of the tags (names, guides, TCRs) you have linked to each FlowCode
combination. The expression in these channels is the expression of that
combination of FlowCode epitopes in the cells. Additionally, you should
have a “FlowCode” channel, which is the expression level of the epiopes
in cells with valid FlowCode combinations. In other words, the
“FlowCode” channel measures transduction.

So, in FlowJo, you can gate on your cells and then plot a tag combo
versus whatever markers you have in your flow panel. This is useful for
checking CRISPR guide efficiency (for example, plotting GATA-3
expression in Gata3-targeted cells) or for phenotyping. If you want, you
can run dimensionality reduction and clustering approaches in FlowJo, R
or Python on the FCS files as they are or with some pre-gating.

For our analyses, we recommend pre-gating to specific biologically
relevant cell populations. I still need to write the part to do this in
R. At present, you should be able to treat the unmixed FCS files you get
from this workflow as you would previously, running them through Orian
Bricard’s
[FlowCodeDecoder](https://github.com/DrCytometer/FlowcodeDecoder). The
difference is that the unmixing will now be correct (hopefully), both
for the barcodes and for the phenotyping and gating markers.
