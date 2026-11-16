# Saved SIR simulations and synthetic field observations for calibration

A 6 by 6 grid with 26 observation times, including physical time zero.
Training uses 40 simulator inputs; ten further inputs are held out.
Responses and observations are on the log(I + 1) scale. This is a
different simulation design from dt_emulation: eta1 bounds are 2 to 4
here, whereas dt_emulation uses inputs near 1 to 2.

## Usage

``` r
dt_calibration
```

## Format

A list with training and validation lists (input matrices and
time-indexed response matrices); named lower and upper parameter bounds;
Nx, Ny, time, and grid metadata; a field list containing observations,
original-grid indices, discrepancy correlation, and prior settings; a
separate truth list used only for validation/plots; and provenance
containing input identifiers, response scale, seed, and source
checksums.

## Source

Saved SIR simulations following Equation (26) in Banerjee, Chen,
Frankenburg and Zhou (2025), JMLR 26(146). The field-generation/import
script is maintained in data-raw.

## See also

[`calibration_prepare()`](https://xiangchen-stat.github.io/spDBL/reference/calibration_prepare.md),
[`calibration_fit()`](https://xiangchen-stat.github.io/spDBL/reference/calibration_fit.md)

## Examples

``` r
data("dt_calibration")
d <- dt_calibration
rbind(lower = d$lower, upper = d$upper, truth = d$truth$eta)
#>       eta1 eta2 alpha1 alpha2 alpha3
#> lower    2  0.2      0    0.0      0
#> upper    4  0.4      0    0.2      0
#> truth    3  0.3      0    0.1      0
dim(d$field$z)
#> [1] 26 36
length(d$training$responses)
#> [1] 26
```
