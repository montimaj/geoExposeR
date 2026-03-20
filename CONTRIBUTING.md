# Contributing to geoExposeR

Thank you for your interest in contributing to geoExposeR! This document provides guidelines for contributing to the project.

## Table of Contents
- [Code of Conduct](#code-of-conduct)
- [Getting Started](#getting-started)
- [Development Setup](#development-setup)
- [How to Contribute](#how-to-contribute)
- [Reporting Issues](#reporting-issues)
- [Submitting Merge Requests](#submitting-merge-requests)
- [Development Guidelines](#development-guidelines)
- [Testing](#testing)
- [Documentation](#documentation)
- [Release Process](#release-process)

## Code of Conduct

We are committed to providing a welcoming and inclusive environment for all contributors. Please be respectful and professional in all interactions.

### Expected Behavior
- Use welcoming and inclusive language
- Be respectful of differing viewpoints and experiences
- Gracefully accept constructive criticism
- Focus on what is best for the community
- Show empathy towards other community members

### Unacceptable Behavior
- Harassment, discrimination, or offensive comments
- Personal attacks or trolling
- Publishing others' private information without permission
- Any conduct that would be inappropriate in a professional setting

## Getting Started

### Prerequisites
- R (≥ 4.4.0)
- VS Code (recommended)
- Git
- [Pandoc](https://pandoc.org/installing.html) for building vignettes

### Required R Packages
```r
# Development tools
install.packages(c("devtools", "roxygen2", "testthat", "covr"))

# Package dependencies (see DESCRIPTION file)
install.packages(c("data.table", "dplyr", "lme4", "mice", "broom.mixed"))
```

## Development Setup

1. **Fork the repository** on GitLab
2. **Clone your fork** locally:
   ```bash
   git clone https://code.usgs.gov/YOUR_USERNAME/geoExposeR.git
   cd geoExposeR
   ```

3. **Set up the development environment**:
   ```r
   # In R/RStudio
   devtools::load_all()
   devtools::test()
   devtools::check()
   ```

4. **Create a new branch** for your contribution:
   ```bash
   git checkout -b feature/your-feature-name
   # or
   git checkout -b fix/issue-number
   ```

## How to Contribute

We welcome several types of contributions:

### 🐛 Bug Reports
- Search existing issues first
- Use the bug report template
- Include reproducible examples
- Specify R version and package version

### 💡 Feature Requests
- Check if the feature already exists
- Clearly describe the use case
- Explain why it would be valuable
- Consider backward compatibility

### 📖 Documentation Improvements
- Fix typos or unclear explanations
- Add examples or clarify existing ones
- Improve function documentation
- Update vignettes

### 🧪 Code Contributions
- Bug fixes
- New features
- Performance improvements
- Test coverage improvements

### 📊 Examples and Use Cases
- Real-world application examples
- Educational materials
- Tutorial improvements

## Reporting Issues

When reporting issues, please include:

### For Bug Reports
1. **Describe the bug**: Clear description of what happened
2. **Reproduction steps**: Minimal code to reproduce the issue
3. **Expected behavior**: What you expected to happen
4. **System information**:
   ```r
   sessionInfo()
   packageVersion("geoExposeR")
   ```
5. **Data**: If possible, provide sample data that reproduces the issue

### For Feature Requests
1. **Feature description**: What functionality you'd like to see
2. **Use case**: Why this feature would be useful
3. **Proposed implementation**: If you have ideas about how to implement it

## Submitting Merge Requests

### Before Submitting
1. **Ensure your code passes all checks**:
   ```r
   devtools::load_all()
   devtools::test()
   devtools::check()
   ```

2. **Update documentation** if needed:
   ```r
   devtools::document()
   ```

3. **Add tests** for new functionality

4. **Update NEWS.md** with your changes

### Merge Request Process
1. **Commit your changes** with clear, descriptive messages:
   ```bash
   git add .
   git commit -m "Add: brief description of changes"
   ```

2. **Push to your fork**:
   ```bash
   git push origin feature/your-feature-name
   ```

3. **Create a merge request** on GitLab with:
   - Clear title describing the change
   - Detailed description of what was changed and why
   - Reference to any related issues
   - Confirmation that tests pass

4. **Respond to review feedback** promptly and professionally

### Merge Request Guidelines
- Keep changes focused and atomic
- Write clear commit messages
- Include tests for new functionality
- Update documentation as needed
- Ensure backward compatibility when possible

## Development Guidelines

### Code Style
Follow the [tidyverse style guide](https://style.tidyverse.org/):

```r
# Good
calculate_arsenic_exposure <- function(data, method = "weighted") {
  if (is.null(data)) {
    stop("Data cannot be NULL")
  }
  
  result <- data %>%
    filter(!is.na(arsenic_level)) %>%
    summarize(mean_exposure = mean(arsenic_level))
  
  return(result)
}

# Function names: snake_case
# Variable names: snake_case
# Constants: UPPER_SNAKE_CASE
# Use explicit returns
# Add input validation
```

### Package Structure
```
geoExposeR/
├── R/
│   ├── geoExposeR.R          # Main exported function
│   ├── data-loading.R        # Data loading functions
│   ├── imputation.R          # Imputation functions
│   ├── regression.R          # Analysis functions
│   └── utils.R               # Utility functions
├── tests/testthat/           # Test files
├── man/                      # Generated documentation
├── vignettes/               # Package vignettes
└── inst/                    # Additional package files
```

### Function Guidelines
1. **Single responsibility**: Each function should do one thing well
2. **Clear naming**: Function names should describe what they do
3. **Input validation**: Check arguments and provide helpful error messages
4. **Documentation**: Use roxygen2 for all exported functions
5. **Error handling**: Use informative error messages

### Internal Functions
- Use `@keywords internal` for functions not intended for end users
- Prefix with `.` if the function is truly internal (e.g., `.validate_input`)
- Keep internal functions focused and well-documented

## Testing

### Writing Tests
- Use `testthat` framework
- Test both success and failure cases
- Use descriptive test names
- Test edge cases and error conditions

```r
test_that("load_prob_data handles missing columns gracefully", {
  # Create test data missing required column
  test_data <- data.frame(wrong_col = 1:10)

  expect_error(
    load_prob_data(test_data, required_cols = "correct_col"),
    "missing columns"
  )
})

test_that("arsenic imputation produces valid probabilities", {
  results <- impute_arsenic_exposure(test_data, ndraws = 2)
  
  expect_true(all(results$probabilities >= 0))
  expect_true(all(results$probabilities <= 1))
  expect_equal(length(results$datasets), 2)
})
```

### Running Tests
```bash
# Run all tests
Rscript -e "devtools::test()"

# Run specific test file
Rscript -e "testthat::test_file('tests/testthat/test-imputation.R')"

# Check test coverage
Rscript -e "covr::package_coverage()"
```

### Test Data
- Use synthetic data for tests
- Keep test datasets small for speed
- Store test data in `tests/testthat/` directory
- Use `helper-*.R` files for test data generation functions

## Documentation

### Function Documentation
Use roxygen2 for all exported functions:

```r
#' Perform Arsenic Exposure Analysis
#'
#' This function performs a comprehensive analysis of arsenic exposure
#' using multiple imputation and mixed-effects modeling.
#'
#' @param data A data frame containing the input data
#' @param ndraws An integer specifying the number of imputed datasets (default: 10)
#' @param output_dir A character string specifying the output directory
#' @return A list containing analysis results and summary statistics
#' @examples
#' \dontrun{
#' results <- perform_sensitivity_analysis(
#'   data = my_data,
#'   ndraws = 5,
#'   output_dir = "results/"
#' )
#' }
#' @export
```

### Vignettes
- Update vignettes when adding new features
- Include realistic examples
- Explain the statistical methodology
- Show interpretation of results

### Package Documentation
- Update `README.md` for user-facing changes
- Update `NEWS.md` for all changes
- Keep `DESCRIPTION` file current

## Release Process

### Versioning
We use [Semantic Versioning](https://semver.org/):
- **MAJOR.MINOR.PATCH** (e.g., 1.2.3)
- **MAJOR**: Incompatible API changes
- **MINOR**: New functionality (backward compatible)
- **PATCH**: Bug fixes (backward compatible)

### Release Checklist
1. Update version in `DESCRIPTION`
2. Update `NEWS.md` with changes
3. Run comprehensive tests: `devtools::check()`
4. Update documentation: `devtools::document()`
5. Build and test package
6. Create GitLab release with release notes

## Getting Help

### Questions and Discussions
- **GitLab Discussions**: For questions about usage or development
- **GitLab Issues**: For bug reports and feature requests
- **Email**: Contact maintainers for sensitive issues

### Resources
- [R Packages book](https://r-pkgs.org/) by Hadley Wickham
- [Advanced R](https://adv-r.hadley.nz/) for advanced R programming
- [testthat documentation](https://testthat.r-lib.org/) for testing
- [roxygen2 documentation](https://roxygen2.r-lib.org/) for documentation

## Recognition

Contributors will be acknowledged in:
- `DESCRIPTION` file (for significant contributions)
- `README.md` contributors section
- Release notes for their contributions

Thank you for contributing to geoExposeR! Your efforts help advance arsenic exposure research and public health.

---

**Questions?** Feel free to open an issue or contact the maintainers:
- [Dr. Sayantan Majumdar](mailto:sayantan.majumdar@dri.edu)
- [Dr. Matthew O. Gribble](mailto:matt.gribble@ucsf.edu)