# Check RDS Contents

Checks the RDS file for expected elements.

## Usage

``` r
check.rds.contents(obj, required = list(), allow.extra = TRUE)
```

## Arguments

- obj:

  The loaded RDS object.

- required:

  List of required elements.

- allow.extra:

  Logical, whether to allow additional, unexpected elements. Default is
  `TRUE`.

## Value

Returns errors if any found; otherwise, returns `TRUE`.
