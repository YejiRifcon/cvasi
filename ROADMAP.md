# Roadmap

Planned releases and changes.

## Unspecified version

FEATURES

 * Mixtox calculations
 
USABILITY

 * More robust handling of numerics with a self-correcting or *auto* option in
   case of numerical issues
 * Simpler workflows
   * Export of calibration sets as flat files


## Future release in v2.0.0

BREAKING CHANGES

 * `set_exposure()`: default value of argument `reset_times` set to `reset_times=FALSE`
 * `calibrate()` will only pass selected arguments to `optim()`, e.g. usage of
   `lower` and `upper` are prohibited.

## Released in v1.5.0 ✓

FEATURES

 * Log-likelihood profiling based on count data
 * Import of study data from flat files, automatic conversion to calibration sets

DEPRECATED
 * `CalibrationSet()` replaced by `caliset()`
 * `simulate_batch()` replaced by the more flexible `batch()`

