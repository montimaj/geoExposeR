#!/usr/bin/env Rscript

# Local Coverage Testing Script for gwArsenicR
# This script replicates the GitLab CI/CD coverage workflow locally
# Run this from the root directory of your gwArsenicR package

cat("=== Local Coverage Testing for gwArsenicR ===\n")
cat("Starting at:", as.character(Sys.time()), "\n\n")

# Set CRAN mirror to avoid repository errors
options(repos = c(CRAN = "https://cloud.r-project.org/"))

# Check if we're in the right directory
if (!file.exists("DESCRIPTION") ||
      !file.exists(file.path("R", "gwArsenic.R"))) {
  cat("❌ ERROR: Please run this script from the root of your ",
      "gwArsenicR package\n", sep = "")
  cat("Expected files: DESCRIPTION, R/gwArsenic.R\n")
  quit(status = 1)
}

# Load required packages
required_packages <- c("devtools", "covr", "testthat", "DT", "htmltools")
cat("Checking required packages...\n")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing missing package:", pkg, "\n")
    install.packages(pkg, repos = "https://cloud.r-project.org/")
  }
}

# Load all required libraries
suppressMessages({
  library(devtools)
  library(covr)
  library(testthat)
})

cat("✅ All required packages loaded\n\n")

# Load package first to make its functions available
cat("Loading package...\n")
tryCatch({
  devtools::load_all(".", quiet = FALSE)
  cat("✅ Package loaded successfully\n\n")
}, error = function(e) {
  cat("❌ ERROR loading package:\n")
  print(e)
  quit(status = 1)
})

# Generate coverage for the entire package
cat("Generating coverage for the entire package...\n")
start_time <- Sys.time()

tryCatch({
  full_cov <- covr::package_coverage(
    path = "R/gwArsenic.R",
    test_files = "tests/testhat/test-gwArsenic.R",
    quiet = FALSE,
    clean = FALSE,
    type = "tests"
  )


  end_time <- Sys.time()
  cat("✅ Coverage generation completed in",
      round(as.numeric(end_time - start_time, units = "secs"), 2),
      "seconds\n\n")
}, error = function(e) {
  cat("❌ ERROR generating coverage:\n")
  print(e)
  quit(status = 1)
})

# Calculate and print coverage percentage
cov <- full_cov
class(cov) <- class(full_cov)
attr(cov, "package") <- attr(full_cov, "package")
attr(cov, "relative") <- attr(full_cov, "relative")

# Generate HTML report for gwArsenic.R only
cat("=== Generating gwArsenic.R Coverage Report ===\n")
tryCatch({
  # Filter coverage to show only gwArsenic.R
  gwarsenic_only_cov <- Filter(function(x) {
    if (is.list(x) && !is.null(x$srcref)) {
      srcfile_name <- attr(x$srcref, "srcfile")$filename
      grepl("gwArsenic\\.R", basename(srcfile_name), ignore.case = TRUE)
    } else {
      FALSE
    }
  }, cov)

  # Recreate coverage object with proper class and attributes
  if (length(gwarsenic_only_cov) > 0) {
    class(gwarsenic_only_cov) <- class(cov)
    attr(gwarsenic_only_cov, "package") <- attr(cov, "package")
    attr(gwarsenic_only_cov, "relative") <- attr(cov, "relative")

    # Generate report with filtered coverage
    covr::report(gwarsenic_only_cov, file = "scripts/gwArsenic-coverage.html")
    cat("✅ HTML coverage report generated for gwArsenic.R only\n")
  } else {
    cat("⚠️ No gwArsenic.R coverage found, generating full report\n")
    covr::report(cov, file = "scripts/gwArsenic-coverage.html")
  }
}, error = function(e) {
  cat("⚠️ Warning: Could not generate HTML report:\n")
  print(e)
  cat("Trying fallback approach...\n")
  # Fallback: try to generate report with original coverage
  tryCatch({
    covr::report(cov, file = "scripts/gwArsenic-coverage.html")
    cat("✅ Fallback HTML report generated\n")
  }, error = function(e2) {
    cat("❌ Could not generate any HTML report\n")
  })
})

# Create a focused summary file
coverage_pct <- covr::percent_coverage(gwarsenic_only_cov)
cat(sprintf("📊 gwArsenic.R Coverage: %.1f%%\n\n", coverage_pct))
total_lines <- length(cov)
uncovered_lines <- ceiling(
  total_lines - coverage_pct / 100 * total_lines
)
coverage_summary <- data.frame(
  File = "gwArsenic.R",
  Coverage = sprintf("%.1f%%", coverage_pct),
  Status = if (coverage_pct >= 85) {
    "✅ EXCELLENT"
  } else if (coverage_pct >= 75) {
    "✅ GOOD"
  } else if (coverage_pct >= 65) {
    "⚠️ NEEDS IMPROVEMENT"
  } else {
    "❌ INSUFFICIENT"
  },
  Target = "≥85%",
  Uncovered_Lines = uncovered_lines,
  Total_Lines = total_lines,
  stringsAsFactors = FALSE
)

cat("\n=== gwArsenic.R Coverage Summary ===\n")
print(coverage_summary)

# Write summary to CSV
write.csv(
  coverage_summary,
  "scripts/gwArsenic-coverage-summary.csv",
  row.names = FALSE
)
cat("✅ Summary saved to: gwArsenic-coverage-summary.csv\n\n")

# Set quality gates and provide assessment
cat("=== Quality Gate Assessment ===\n")
if (coverage_pct < 65) {
  cat(sprintf(
    "❌ FAIL: gwArsenic.R coverage %.1f%% is below minimum threshold (65%%)\n",
    coverage_pct
  ))
  cat("🔧 Recommendation: Add more comprehensive tests\n")
  status_code <- 1
} else if (coverage_pct < 75) {
  cat(sprintf(
    paste0(
      "⚠️ WARNING: gwArsenic.R coverage %.1f%% is below recommended ",
      "threshold (75%%)\n"
    ),
    coverage_pct
  ))
  cat("🔧 Recommendation: Consider adding tests for edge cases\n")
  status_code <- 0
} else {
  cat(sprintf(
    "✅ SUCCESS: gwArsenic.R coverage %.1f%% meets quality standards\n",
    coverage_pct
  ))
  status_code <- 0
}

# Additional recommendations
if (coverage_pct < 85) {
  cat("\n📝 Suggestions to improve coverage:\n")
  cat("1. Add tests for error conditions and edge cases\n")
  cat("2. Test all branches of conditional statements\n")
  cat("3. Ensure all validation functions are tested\n")
  cat("4. Test both success and failure scenarios\n")
}

cat("\n=== Coverage Testing Complete ===\n")
cat("Files generated:\n")
if (file.exists("scripts/gwArsenic-coverage.html")) {
  cat("- gwArsenic-coverage.html (detailed HTML report)\n")
}
cat("- gwArsenic-coverage-summary.csv (summary data)\n")
cat("Completed at:", as.character(Sys.time()), "\n")

# Open HTML report if available
if (file.exists("scripts/gwArsenic-coverage.html")) {
  cat("\n🌐 Opening coverage report in browser...\n")
  tryCatch({
    if (.Platform$OS.type == "windows") {
      shell.exec("scripts/gwArsenic-coverage.html")
    } else if (Sys.info()["sysname"] == "Darwin") {
      system("open scripts/gwArsenic-coverage.html")
    } else {
      system("xdg-open scripts/gwArsenic-coverage.html")
    }
  }, error = function(e) {
    cat("Could not auto-open browser.\n")
    cat("Please manually open: scripts/gwArsenic-coverage.html\n")
  })
}

quit(status = status_code)
