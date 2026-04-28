# geoExposeR 1.0.0

## First Stable Release

This is the first stable release of geoExposeR, a feature-complete R package for modeling health effects of environmental exposures. Originally designed for drinking water contaminants but broadly applicable to other environmental exposures.

### Highlights
* **100% test coverage** with 476 comprehensive tests
* **12 required dependencies** for robust functionality
* **7 exported functions** for data ingestion, analysis, and validation
* **Complete documentation** including vignettes and pkgdown site
* **CC0 1.0 Universal license** for unrestricted use

### New Features

#### Core Functionality
* **Main Function**: `perform_sensitivity_analysis()` - Complete workflow for exposure analysis (e.g., arsenic in drinking water)
* **Data Integration**: Combines modeled probability data with measured concentration data
* **Multiple Imputation**: Implements probabilistic exposure assignment across multiple datasets
* **Statistical Analysis**: Mixed-effects regression with proper pooling using Rubin's Rules
* **MICE Integration**: Optional multiple imputation for missing covariates

#### Data Loading and Processing
* **Probability Data Support**: Load and process modeled probability data from GeoTIFF rasters or tabular sources
* **Concentration Data Integration**: Convert measured concentration data to lognormal-derived multinomial probabilities
* **Weighted Combination**: Population-weighted integration of probability and concentration models
* **Flexible Identifier Support**: Works with any geographic or participant-level identifier (FIPS codes, census tracts, ZIP codes, or custom IDs)
* **Smart ID Formatting**: Automatic detection of identifier type—applies zero-padding for numeric FIPS codes (1-5 digits) while preserving alphanumeric identifiers unchanged. Controlled via `format_ids` parameter (enabled by default)
* **Birth Data Processing**: Specialized handling of health outcome datasets
* **Input Validation**: `validate_prepared_inputs()` for comprehensive pre-analysis data checks

#### Spatial Data Support
* **Raster Extraction**: Extract well water usage percentages from raster data using `add_well_percentages()`
* **Geographic Subsetting**: `subset_by_geography()` for filtering data by spatial boundaries
* **Coordinate Handling**: CRS transformation and centroid computation via sf and terra

#### Multiple Imputation Framework
* **Exposure Imputation**: Geographic-level probabilistic assignment based on combined models
* **Covariate Imputation**: MICE-based imputation for missing demographic and health variables
* **Validation Tools**: Comprehensive checks for imputation quality and convergence
* **Flexible Configuration**: Customizable imputation parameters and methods

#### Statistical Analysis
* **Mixed-Effects Models**: Support for complex nested data structures
* **Custom Formulas**: Flexible regression formula specification
* **Multiple Outcomes**: Simultaneous analysis of multiple health endpoints
* **Rubin's Rules**: Proper pooling of estimates across imputed datasets
* **Confidence Intervals**: Accurate uncertainty quantification

### Package Infrastructure

#### Testing and Quality Assurance
* **Comprehensive Test Suite**: 100% code coverage with automated testing
* **Synthetic Data Testing**: All tests use generated dummy data for reproducibility
* **GitLab CI/CD**: Automated testing across multiple R versions and operating systems
* **Performance Optimization**: Fast execution with `ndraws = 2` for development testing

#### Documentation and Usability
* **Extensive Documentation**: Complete roxygen2 documentation for all functions
* **Package Website**: Automated pkgdown site generation and deployment
* **Vignettes**: Step-by-step tutorials with working examples
* **README**: Comprehensive installation and usage instructions

#### Development Tools
* **Modular Code Organization**: Separate files for data ingestion, data loading, imputation, regression, and utilities
* **Internal Function Architecture**: Well-organized internal functions with clear responsibilities
* **Contribution Guidelines**: Detailed CONTRIBUTING.md with development setup instructions
* **Code Style**: Follows tidyverse style guidelines with consistent formatting

### Dependencies

#### Required Packages
* **Amelia**: Multiple imputation using bootstrapped EM algorithm
* **broom.mixed**: Tidy statistical output formatting
* **data.table**: High-performance data manipulation
* **dplyr**: Data transformation and summarization
* **exactextractr**: Fast raster value extraction for spatial analysis
* **Hmisc**: Weighted statistical functions
* **lme4**: Mixed-effects modeling
* **mice**: Multiple imputation by chained equations
* **rlang**: Tidy evaluation framework
* **rms**: Regression modeling strategies
* **sf**: Simple features for spatial vector data
* **terra**: Raster and spatial data processing

#### System Requirements
* **R** (>= 4.4.0): Minimum R version for compatibility
* **Pandoc**: Required for vignette building and documentation

### Performance and Scalability

#### Optimization Features
* **Vectorized Operations**: Efficient matrix operations for probability calculations
* **Memory Management**: Optimized data structures for large datasets
* **Parallel Processing**: Multi-core support for data loading operations
* **Caching**: Efficient dependency caching in CI/CD workflows

#### Practical Considerations
* **Geographic-Level Analysis**: Designed for exposure assessment at any geographic scale (counties, census tracts, ZIP codes, or custom regions)
* **Flexible Sample Sizes**: Supports datasets from small studies to large epidemiological cohorts
* **Configurable Parameters**: Adjustable imputation counts and convergence criteria
* **Output Management**: Comprehensive result saving and structured output formats

### Academic Integration

#### Methodological Foundation
* **Published Methods**: Implements approaches from Bulka et al. (2022) and Lombard et al. (2021)
* **Statistical Rigor**: Follows best practices for uncertainty quantification in exposure assessment
* **Reproducible Research**: Seed control and comprehensive output logging
* **Validation**: Extensive testing against known statistical properties

#### Research Applications
* **Birth Outcomes**: Specialized support for pregnancy and birth outcome studies
* **Environmental Epidemiology**: Designed for environmental health research workflows
* **Risk Assessment**: Tools for population-level exposure estimation
* **Policy Research**: Support for regulatory and public health decision-making

### Known Limitations

* **Geographic Scope**: Designed for US-based studies but supports any identifier format for international applications
* **Data Format Requirements**: Specific column naming conventions required for input data
* **Memory Usage**: Large datasets may require substantial memory for multiple imputation
* **Convergence**: MICE convergence warnings expected for complex imputation models

### Future Directions

* **International Support**: Expansion to non-US geographic coding systems
* **Additional Applications**: Expanding worked examples for other contaminants and exposure pathways
* **Advanced Modeling**: Integration of spatial autocorrelation and temporal trends
* **Visualization Tools**: Enhanced plotting and diagnostic visualization capabilities

---

## Development Team

Dr. Sayantan Majumdar (DRI), Dr. Scott M. Bartell (UC Irvine), Dr. Melissa A. Lombard (USGS), Dr. Ryan G. Smith (Colorado State University), Dr. Matthew O. Gribble (UCSF)

## Funding

This work was supported by the National Heart, Lung, and Blood Institute (R21HL159574) and funding from the United States Geological Survey’s John Wesley Powell Center for Analysis and Synthesis.

## Citation

Software release:

Majumdar, S., Bartell, S. M., Lombard, M. A., Smith, R. G., and Gribble, M. O., 2026, geoExposeR: An R package for modeling health effects of environmental exposures: U.S. Geological Survey software release, https://doi.org/10.5066/P1JGUKMD

Accompanying paper (under review at JOSS):

Majumdar, S., Bartell, S. M., Lombard, M. A., Smith, R. G., and Gribble, M. O., 2026, geoExposeR: An R package for modeling health effects of environmental exposures, Journal of Open Source Software* (under review).