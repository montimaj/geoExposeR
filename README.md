# gwArsenicR
An R package for modeling and estimating arsenic exposure from groundwater, based on epidemiological studies and existing concentration models.

<!-- Package Info -->
[![Version](https://img.shields.io/badge/version-0.1.0-blue.svg)](https://gitlab.com/s.majumdar_mst.edu/gwArsenicR)
[![R](https://img.shields.io/badge/R-%E2%89%A5%204.4.0-blue.svg)](https://www.r-project.org/)
[![License](https://img.shields.io/badge/License-CC0%201.0-blue.svg)](https://creativecommons.org/publicdomain/zero/1.0/)

<!-- Build Status -->
[![Pipeline Status](https://code.usgs.gov/envirohealth/gwArsenicR/badges/main/pipeline.svg)](https://code.usgs.gov/envirohealth/gwArsenicR/-/pipelines)

<!-- Quality Metrics -->
[![Test Coverage](https://code.usgs.gov/envirohealth/gwArsenicR/badges/main/coverage.svg)](https://code.usgs.gov/envirohealth/gwArsenicR/-/commits/main)

<!-- Documentation and Usage -->
[![Documentation](https://img.shields.io/badge/docs-pkgdown-blue.svg)](https://gwarsenicr.code-pages.usgs.gov)
[![Vignettes](https://img.shields.io/badge/vignettes-available-blue.svg)](https://gwarsenicr.code-pages.usgs.gov/articles/)
[![GitLab](https://img.shields.io/badge/GitLab-repository-blue.svg)](https://code.usgs.gov/envirohealth/gwArsenicR/)

<!-- Development Status -->
[![Project Status: Active](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)

<!-- Code Quality -->
[![Code style: tidyverse](https://img.shields.io/badge/code%20style-tidyverse-blue.svg)](https://style.tidyverse.org)
[![Dependencies](https://img.shields.io/badge/dependencies-5-blue.svg)](https://code.usgs.gov/envirohealth/gwArsenicR/-/blob/main/DESCRIPTION)

<!-- Academic -->
[![DOI](https://img.shields.io/badge/DOI-pending-orange.svg)](https://gitlab.com/s.majumdar_mst.edu/gwArsenicR)
<!-- [![JOSS](https://joss.theoj.org/papers/10.21105/joss.XXXXX/status.svg)](https://joss.theoj.org/papers/10.21105/joss.XXXXX) -->

**Principal Investigator**: [Dr. Matthew O. Gribble](https://profiles.ucsf.edu/matthew.gribble) [matt.gribble@ucsf.edu]

**Co-Investigators**: [Dr. Sayantan Majumdar](https://www.dri.edu/directory/sayantan-majumdar/) [sayantan.majumdar@dri.edu], [Dr. Ryan G. Smith](https://www.engr.colostate.edu/ce/ryan-g-smith/) [ryan.g.smith@colostate.edu]

## Table of Contents
- [Statement of Need](#statement-of-need)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Data Requirements](#data-requirements)
- [Package Structure](#package-structure)
- [Core Functions](#core-functions)
- [Key Features](#key-features)
- [Documentation](#documentation)
- [Running Tests](#running-tests)
- [Modifying the R Package](#modifying-the-r-package)
- [Contributing](#contributing)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)
- [License](#license)

## Statement of Need

gwArsenicR is an R package that provides a streamlined workflow for epidemiologists and public health researchers to estimate arsenic exposure from private and public well water and to assess its association with various health outcomes. The package implements a sophisticated statistical approach that combines geospatial arsenic prediction models with a multiple imputation framework to account for exposure uncertainty. The core functionality integrates predicted arsenic concentration probabilities from U.S. Geological Survey (USGS) models for private wells with U.S. Environmental Protection Agency (EPA) data for public water systems. It then uses a weighting scheme based on the proportion of the population served by each water source to create a unified, county-level exposure probability distribution. By encapsulating these complex methodologies into a single, user-friendly function, gwArsenicR makes this type of analysis more accessible, reproducible, and standardized.

Existing methods from studies like _Bulka et al. (2022)_ and _Lombard et al. (2021)_ provide robust statistical frameworks but require complex implementation involving:
- Integration of multiple probabilistic models (USGS and EPA)
- Multiple imputation techniques for uncertain exposures
- Mixed-effects regression with proper pooling of results

**gwArsenicR solves this** by packaging these sophisticated methods into a single, user-friendly function. This enables researchers to:

- 🔬 **Focus on science**, not statistical programming
- 📊 **Ensure reproducibility** across studies  
- ⚡ **Accelerate research** on arsenic health effects
- 🎯 **Apply best practices** for uncertainty quantification

By democratizing access to these methods, gwArsenicR promotes more rigorous and comparable arsenic exposure research.

## Installation
### System Dependencies: Pandoc

To build the package vignettes, which provide detailed documentation and examples, you need to have Pandoc installed. RStudio bundles a version of Pandoc, so if you are an RStudio user, you may not need to install it separately.

If you do not have Pandoc, you can download and install it from the official website: [pandoc.org/installing.html](https://pandoc.org/installing.html).

### R Package Installation

**Prerequisites:** R (≥ 4.4.0) and [Pandoc](https://pandoc.org/installing.html) for building vignettes.

**Install from GitLab:**
```r
# Install devtools if not already installed
if (!require("devtools")) install.packages("devtools")

# Install gwArsenicR
devtools::install_git("https://code.usgs.gov/envirohealth/gwArsenicR.git", build_vignettes = TRUE)

# Load the package
library(gwArsenicR)
```

**For development and testing (includes coverage tools):**
```r
# Install additional development dependencies
install.packages(c("devtools", "testthat", "covr", "DT", "htmltools"))

# Install gwArsenicR
devtools::install_git("https://code.usgs.gov/envirohealth/gwArsenicR.git", build_vignettes = TRUE)
```

**Verify installation:**
```r
# Check package version
packageVersion("gwArsenicR")

# View help
?perform_sensitivity_analysis
```

**Package Documentation Website:**

The complete package documentation, including detailed function references, vignettes, and examples, is available at: **https://gwarsenicr.code-pages.usgs.gov**

This website is automatically generated and updated whenever changes are made to the main branch, ensuring you always have access to the latest documentation.

### Modifying the R Package
If you want to customize this package, run the following using R after you have modified the codes.
```
# Set the working directory to your package directory
setwd(<path to gwArsenicR>)

# Load devtools
library(devtools)

# Generate documentation
document()

# Build the package
build()

# Install the package
install()
```

## Quick Start

```r
library(gwArsenicR)

# Run sensitivity analysis with your data
results <- perform_sensitivity_analysis(
  ndraws = 10,
  as_usgs_prob_csv = "path/to/usgs_data.csv",
  as_epa_prob_csv = "path/to/epa_data.csv", 
  birth_data_txt = "path/to/birth_data.txt",
  regression_formula = "~ AsLevel + MAGE_R + rural + (1|FIPS)",
  output_dir = "results/",
  targets = c("OEGEST", "BWT")
)

# View results
print(results)
```

For detailed examples with dummy data, see the [vignette](vignettes/gwArsenicR-vignette.Rmd).

## Package Structure
The gwArsenicR package follows the standard structure for R packages with GitLab CI/CD integration:
```
gwArsenicR/
├── R/
│   ├── data-loading.R                  # Data loading and processing functions
│   ├── gwArsenic.R                     # Main exported function
│   ├── imputation.R                    # Multiple imputation functions
│   ├── regression.R                    # Mixed-effects regression analysis
│   └── utils.R                         # Utility and validation functions
├── man/
│   ├── perform_sensitivity_analysis.Rd
│   └── [other function documentation]
├── tests/
│   ├── testthat/
│   │   ├── helper-create_dummy_data.R  # Test helper functions
│   │   └── test-gwArsenic.R            # Comprehensive test suite
│   └── testthat.R                      # Test runner
├── scripts/
│   └── test_coverage_local.R           # Local coverage testing script
├── vignettes/
│   └── gwArsenicR-vignette.Rmd         # Package tutorial and examples
├── docs/                               # Generated by pkgdown (auto-created)
├── coverage/                           # Coverage reports (auto-created)
├── paper/                              # JOSS paper files
│   ├── paper.md                        # Main paper manuscript
│   ├── paper.bib                       # Bibliography
│   └── paper.html                      # Generated HTML version
├── .gitlab-ci.yml                      # GitLab CI/CD pipeline configuration
|── ._pkgdown.yml                       # Website configuration
├── CONTRIBUTING.md                     # Contribution guidelines
├── DESCRIPTION                         # Package metadata and dependencies
├── NAMESPACE                           # Package exports and imports
├── LICENSE.md                          # CC0 1.0 Universal License
├── NEWS.md                             # Change log and version history
├── README.md                           # README file
└── .gitignore                          # Git ignore rules
```

### R/ Directory Structure

The package uses a **modular architecture** where functionality is organized across multiple well-structured R files:

- **`gwArsenic.R`**: Main exported function `perform_sensitivity_analysis()`
- **`data-loading.R`**: Data loading and processing functions
  - `load_and_process_arsenic_data()` - Orchestrates USGS and EPA data loading
  - `load_usgs_data()` - Loads USGS probability data
  - `convert_epa_to_multinomial()` - Converts EPA lognormal to multinomial probabilities
  - `create_weighted_prob_matrix()` - Creates weighted combined probability matrix
  - `load_and_process_birth_data()` - Loads and processes birth/health outcome data
- **`imputation.R`**: Multiple imputation functions using MICE
  - `impute_arsenic_exposure()` - Creates multiple imputed datasets with arsenic exposure
  - `impute_additional_variables()` - Imputes missing covariates using MICE
  - `validate_imputed_datasets()` - Validates imputation results
- **`regression.R`**: Mixed-effects regression analysis functions
  - `regression_analysis()` - Performs mixed-effects regression on imputed data
  - `pool_estimates_by_term()` - Pools regression estimates using Rubin's Rules
  - `pool_single_estimate()` - Applies Rubin's Rules to individual parameters
- **`utils.R`**: Utility and validation functions
  - `validate_all_inputs()` - Input validation and error checking
  - `format_geographic_ids()` - Formats FIPS codes and geographic identifiers
  - Additional helper functions for data validation and formatting

This modular approach provides several advantages:
- **Clear separation of concerns**: Each file handles a specific aspect of the workflow
- **Easier maintenance**: Functions are logically grouped by purpose
- **Better code organization**: Related functions are kept together
- **Improved testability**: Each module can be tested independently
- **Enhanced readability**: Developers can focus on specific functionality

### Testing Infrastructure

The package includes comprehensive testing and coverage tools with GitLab CI/CD integration:

- **`tests/testthat/test-gwArsenic.R`**: Main test suite with >95% code coverage
- **`tests/testthat/helper-create_dummy_data.R`**: Creates synthetic data for testing
- **`tests/testthat.R`**: Standard R package test runner
- **`scripts/test_coverage_local.R`**: Local coverage analysis script
- **`.gitlab-ci.yml`**: GitLab CI/CD pipeline with automated testing
  - Multi-platform testing (Ubuntu with multiple R versions)
  - Automated coverage reporting and quality gates
  - Documentation validation
  - Security scanning
  - Automated pkgdown site deployment

### GitLab CI/CD Pipeline

The package uses GitLab CI/CD for continuous integration with the following stages:

1. **test**: Run R package tests on multiple versions/platforms
   - R release (Ubuntu latest)
   - R development version (Ubuntu latest)  
   - R oldrel-1 (Ubuntu with R 4.3.0)

2. **coverage**: Generate and analyze test coverage
   - Minimum threshold: 65% (fails pipeline if below)
   - Recommended threshold: 75% (warning if below)
   - Target threshold: 85% (optimal)

3. **docs**: Validate documentation is up to date
   - Checks that roxygen2 documentation is current
   - Validates package structure

4. **deploy**: Deploy to GitLab Pages and create releases
   - Automatic pkgdown site deployment
   - Release artifact creation on git tags

### Key Internal Functions

While users primarily interact with `perform_sensitivity_analysis()`, the package includes several internal functions organized by module:

**Data Loading (`data-loading.R`):**
- `load_and_process_arsenic_data()`: Orchestrates USGS and EPA data loading
- `load_usgs_data()`: Loads USGS probability data
- `convert_epa_to_multinomial()`: Converts EPA lognormal to multinomial probabilities
- `create_weighted_prob_matrix()`: Creates weighted combined probability matrix
- `load_and_process_birth_data()`: Loads and processes birth/health outcome data

**Multiple Imputation (`imputation.R`):**
- `impute_arsenic_exposure()`: Creates multiple imputed datasets with arsenic exposure
- `impute_additional_variables()`: Imputes missing covariates using MICE
- `validate_imputed_datasets()`: Validates imputation results

**Regression Analysis (`regression.R`):**
- `regression_analysis()`: Performs mixed-effects regression on imputed data
- `pool_estimates_by_term()`: Pools regression estimates using Rubin's Rules
- `pool_single_estimate()`: Applies Rubin's Rules to individual parameters

**Utilities:**
- `validate_all_inputs()`: Input validation and error checking
- `format_geographic_ids()`: Formats FIPS codes and geographic identifiers
- Additional helper functions for data validation and formatting

## Core Functions
The package's primary functionality is exposed through a single main function:

**`perform_sensitivity_analysis()`**: This is the main exported function of the package. It orchestrates the entire workflow from data loading through final analysis:

1. **Data Loading**: Loads and processes USGS probability data, EPA lognormal parameters, and birth/health outcome data
2. **Probability Integration**: Combines USGS and EPA models using weighted averages based on private well usage
3. **Multiple Imputation**: Creates multiple datasets with probabilistically assigned arsenic exposure levels
4. **Covariate Imputation**: Optionally imputes missing values in additional covariates using MICE
5. **Statistical Analysis**: Fits mixed-effects regression models to assess exposure-outcome relationships
6. **Results Pooling**: Applies Rubin's Rules to pool results across imputed datasets
7. **Output Generation**: Saves results and returns structured analysis output

The function is designed to handle the complex statistical methodology while providing a simple, user-friendly interface that requires minimal statistical programming expertise.

## Key Features

- ✅ **Automated workflow**: Single function handles entire analysis pipeline
- ✅ **Multiple imputation**: Accounts for uncertainty in arsenic exposure assignment
- ✅ **Flexible modeling**: Supports custom regression formulas and multiple outcomes
- ✅ **Robust statistics**: Implements Rubin's Rules for proper inference
- ✅ **Data integration**: Combines USGS and EPA models with population weighting
- ✅ **Missing data handling**: Optional MICE imputation for covariates
- ✅ **Reproducible results**: Seed control and comprehensive output saving
- ✅ **Extensive testing**: >95% code coverage with automated CI/CD

## Documentation

### Online Documentation
📖 **Complete documentation website**: https://gwarsenicr.code-pages.usgs.gov

The documentation website includes:
- **Function Reference**: Detailed documentation for all package functions
- **Vignettes**: Step-by-step tutorials with working examples
- **Getting Started Guide**: Quick introduction to package usage
- **Methodology**: Statistical background and implementation details
- **FAQ**: Common questions and troubleshooting

### Local Documentation
```r
# View package help
help(package = "gwArsenicR")

# View main function documentation
?perform_sensitivity_analysis

# Browse vignettes
browseVignettes("gwArsenicR")
```

## Data Requirements

The package requires three input files:

1. **USGS Probability Data** (`as_usgs_prob_csv`): CSV with arsenic concentration probabilities and geographic identifiers
2. **EPA Parameters** (`as_epa_prob_csv`): CSV with EPA lognormal distribution parameters  
3. **Health Outcome Data** (`birth_data_txt`): Text file with health outcomes and county identifiers

**Expected columns:**
- USGS data: `GEOID10`, `RFC3_C1v2`, `RFC3_C2v2`, `RFC3_C3v2`, `Wells_2010`
- EPA data: `EPA_AS_meanlog`, `PWELL_private_pct`
- Health data: `FIPS`, outcome variables (e.g., `BWT`, `OEGEST`), covariates

## Running Tests

The package includes comprehensive tests to ensure reliability and correctness. All tests use synthetic dummy data and run quickly for efficient development and CI/CD workflows.

**Run all tests:**
```bash
cd /path/to/gwArsenicR
Rscript -e "devtools::test()"
```

**Run specific test file:**
```bash
Rscript -e "testthat::test_file('tests/testthat/test-gwArsenic.R')"
```

**Run tests during development (faster - no package rebuild):**
```bash
Rscript -e "devtools::load_all(); devtools::test()"
```

**Run tests with clean output (suppress MICE warnings):**
```bash
Rscript -e "suppressWarnings(devtools::test())"
```

**Check package integrity:**
```bash
Rscript -e "devtools::check()"
```

**Interactive testing in R/RStudio:**
```r
# Install coverage dependencies if needed
if (!require("DT")) install.packages("DT")
if (!require("htmltools")) install.packages("htmltools")

# Load package for development
devtools::load_all()

# Run all tests
devtools::test()

# Run specific test with detailed output
testthat::test_file("tests/testthat/test-gwArsenic.R", reporter = "progress")

# Check test coverage (requires DT and htmltools)
Rscript scripts/test_coverage_local.R
```

### Testing Configuration

Tests are configured for speed and reliability:
- **Fast execution**: Uses `ndraws = 2` for quick multiple imputation
- **Synthetic data**: All tests use generated dummy data, no external dependencies
- **Comprehensive coverage**: Tests data loading, imputation, regression, and output generation
- **Warning suppression**: Expected MICE convergence warnings are suppressed for clean output
- **Automatic cleanup**: Temporary files are automatically removed after each test

### Expected Test Output

A successful test run will show:
```
✅ | F W  S  OK | Context
⠹ |          3 | Test perform_sensitivity_analysis with comprehensive coverage                                                                                                                                  
--- Pooled Analysis Results ---
$BWT
    term      q.mi    se.mi statistic  conf.low conf.high   p.value
1 As5-10  6.675316 32.49704 0.2054131 -60.42082  73.77145 0.8389941
2  As10+ 19.148237 32.77495 0.5842339 -45.79132  84.08780 0.5602382

$OEGEST
    term          q.mi     se.mi     statistic   conf.low conf.high   p.value
1 As5-10  3.618142e-02 0.1502064  0.2408779292 -0.2780558 0.3504186 0.8122120
2  As10+ -5.288651e-05 0.1693747 -0.0003122457 -0.3491999 0.3490941 0.9997534

-----------------------------
[1] "Checking results for target: BWT"
[1] "Checking results for target: OEGEST"
✅ |        189 | Test perform_sensitivity_analysis with comprehensive coverage [5.6s]                                                                                                                           

══ Results═══════════════════════════════════════════════════════════════════
Duration: 5.6 s

✅ [ FAIL 0 | WARN 0 | SKIP 0 | PASS 189 ]
```

### Continuous Integration

The package includes automated testing via GitLab CI/CD pipeline:
- **R CMD check**: Validates package structure and dependencies
- **Multiple R versions**: Tests on R release, devel, and oldrel-1
- **Coverage reporting**: Tracks test coverage with quality gates
- **Documentation checks**: Validates roxygen2 documentation completeness
- **Security scanning**: Optional vulnerability scanning with oysteR
- **Automated deployment**: pkgdown site deployment to GitLab Pages

The pipeline runs on:
- Pushes to main/develop branches
- Merge requests
- Git tags (for releases)
- Changes to R files, tests, DESCRIPTION, or CI configuration

[![Pipeline Status](https://code.usgs.gov/envirohealth/gwArsenicR/badges/main/pipeline.svg)](https://code.usgs.gov/envirohealth/gwArsenicR/-/pipelines)

### Troubleshooting Tests

**If tests fail with data.table errors:**
```bash
# Ensure data.table is properly loaded
Rscript -e "library(data.table); devtools::test()"
```

**If MICE warnings are excessive:**
```bash
# Run with warning suppression
Rscript -e "suppressWarnings(devtools::test())"
```

**If tests timeout:**
```bash
# Check system resources and reduce ndraws in test files if needed
# Default test configuration uses ndraws = 2 for speed

# For very slow systems, you can temporarily modify test parameters:
# Edit tests/testthat/test-gwArsenic.R and change ndraws = 2 to ndraws = 1

# Or run tests with extended timeout
Rscript -e "options(testthat.default_timeout = 600); devtools::test()"

# Skip slow integration tests and run only unit tests
Rscript -e "devtools::test(filter = 'unit')"
```

**Memory issues:**
```bash
# Clear R environment and restart
Rscript -e "rm(list=ls()); gc(); devtools::test()"

# Run tests with reduced memory usage
Rscript -e "options(testthat.progress.max_fails = 10); devtools::test()"
```

### Testing Best Practices

When developing or modifying the package:

1. **Run tests frequently** during development
2. **Add new tests** for new functionality
3. **Use small datasets** in tests for speed (ndraws = 2)
4. **Test edge cases** such as missing data or unusual parameter values
5. **Validate statistical correctness** of pooled results

## Contributing

We welcome contributions! Please see our [contribution guidelines](CONTRIBUTING.md) for details.

**Ways to contribute:**
- 🐛 Report bugs via [GitLab issues](https://code.usgs.gov/envirohealth/gwArsenicR/-/issues)
- 💡 Suggest features or improvements
- 📖 Improve documentation
- 🧪 Add tests or examples
- 🔧 Submit pull requests

**Development setup:**
```bash
git clone https://code.usgs.gov/envirohealth/gwArsenicR.git
cd gwArsenicR
Rscript -e "devtools::load_all(); devtools::test()"
```

## Troubleshooting

### RStudio Installation Issues

Some users may encounter missing dependency errors when installing or using gwArsenicR in RStudio, even after successful installation. This typically manifests as errors like:

```
Error in loadNamespace(i, c(lib.loc, .libPaths()), versionCheck = vI[[i]]) :
  there is no package called 'R6'
```

**Common missing packages:**
- R6
- digest  
- purrr
- vctrs
- glue

**Solution:**
Install these core dependencies manually before installing gwArsenicR:

```r
# Install missing dependencies
install.packages(c("R6", "digest", "purrr", "vctrs", "glue"))

# Then install gwArsenicR
devtools::install_git("https://code.usgs.gov/envirohealth/gwArsenicR.git")
```

**Why this happens:**
This issue typically occurs when:
- RStudio's package installation doesn't properly resolve all transitive dependencies
- System-level R and RStudio R installations differ
- Package library paths are misconfigured
- Previous incomplete installations left corrupted package states

**Alternative solutions:**

1. **Restart R session and try again:**
   ```r
   # In RStudio: Session -> Restart R
   .rs.restartR()
   ```

2. **Update all packages first:**
   ```r
   update.packages(ask = FALSE)
   ```

3. **Install from a clean state:**
   ```r
   # Remove any partial installations
   remove.packages("gwArsenicR")
   
   # Clear package cache
   .libPaths()  # Check library paths
   
   # Reinstall dependencies and package
   install.packages(c("devtools", "R6", "digest", "purrr", "vctrs", "glue"))
   devtools::install_git("https://code.usgs.gov/envirohealth/gwArsenicR.git")
   ```

4. **Check R version compatibility:**
   ```r
   # Ensure you're running R >= 4.4.0
   R.version.string
   ```

If issues persist, please [report them](https://code.usgs.gov/envirohealth/gwArsenicR/-/issues) with your R version and sessionInfo().

## Citation

If you use gwArsenicR in your research, please cite:

> Majumdar, S., Bartell, S. M., Lombard, M. A., Smith, R. G., & Gribble, M. O. (2025). 
> gwArsenicR: An R package for modeling and estimating arsenic exposure from groundwater. 
> *Journal of Open Source Software* (under review).

**BibTeX:**
```bibtex
@article{majumdar2025gwarsenic,
  title={gwArsenicR: An R package for modeling and estimating arsenic exposure from groundwater},
  author={Majumdar, Sayantan and Bartell, Scott M and Lombard, Melissa A and Smith, Ryan G and Gribble, Matthew O},
  journal={Journal of Open Source Software},
  year={2025},
  note={Under review},
  url={https://code.usgs.gov/envirohealth/gwArsenicR}
}
```

**DOI:** *Coming soon upon publication in JOSS*

## Acknowledgments

This work was supported by the National Heart, Lung, and Blood Institute (R21HL159574) and funding from the United States Geological Survey’s John Wesley Powell Center for Analysis and Synthesis.

## License

This project is licensed under the CC0 1.0 Universal License - see the [LICENSE](https://code.usgs.gov/envirohealth/gwArsenicR/-/blob/main/LICENSE.md) file for details.

**Version:** 0.1.0 | **Status:** Under active development

## References
Bulka, C. M., Bryan, M. S., Lombard, M. A., Bartell, S. M., Jones, D. K., Bradley, P. M., ... & Argos, M. (2022). Arsenic in private well water and birth outcomes in the United States. _Environment International, 163_, 107176. https://doi.org/10.1016/j.envint.2022.107176

Lombard, M. A., Bryan, M. S., Jones, D. K., Bulka, C., Bradley, P. M., Backer, L. C., ... & Ayotte, J. D. (2021). Machine learning models of arsenic in private wells throughout the conterminous United States as a tool for exposure assessment in human health studies. _Environmental science & technology, 55_(8), 5012-5023. https://doi.org/10.1021/acs.est.0c05239

Wickham, H., & Bryan, J. (2023). R Packages (2nd ed.). O'Reilly Media. https://r-pkgs.org/
