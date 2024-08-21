# TCA-optimized
**TCA-optimized** is an enhanced and efficient version of the original TCA (Tensor Composition Analysis) R package, designed for the decomposition of bulk methylation array data into cell-specific sorted methylation data. This updated version of the tool brings significant improvements in performance and speed through extensive optimization, including code refactoring and integration of C++ implementations.

## Key Features
- **Cell Fraction Estimation**: Accurately estimate the proportions of different cell types in your bulk methylation data.
- **Cell-Specific Sorted Methylation Data**: Decompose bulk methylation data to obtain cell-specific methylation profiles.

## Installation
You can easily install the latest version of TCA-optimized directly from GitHub using the following devtools command in R:

r
```
devtools::install_github("yourusername/TCA-optimized")
# Make sure you have the devtools package installed:
install.packages("devtools")
```

## Usage
To use **TCA-optimized**, you will need to provide a reference panel data set. The tool utilizes this reference data to accurately decompose bulk methylation profiles into their cell-specific components.

```
# Load the package
library(TCAoptimized)

# Example usage
result <- TCAoptimized::decompose_methylation(bulk_data, reference_panel)

# Access the cell fraction estimation
cell_fractions <- result$cell_fractions

# Access the cell-specific methylation data
cell_specific_data <- result$cell_specific_data
```

## Optimizations
The **TCA-optimized** package has been meticulously improved for better performance:
**Code Refactoring**: The R code has been thoroughly refactored for better readability and maintainability.
**C++ Integration**: Critical parts of the computation have been re-implemented in C++ for faster execution.
**Efficient Data Handling**: Optimized algorithms and data structures to reduce computational overhead and improve memory management.

## Applications
This tool is ideal for researchers working with bulk methylation data who need to:
- Estimate the composition of different cell types within their samples.
- Obtain cell-specific methylation profiles from mixed samples.

## Contributions
Contributions to the TCA-optimized package are welcome! Please feel free to submit issues or pull requests on the GitHub repository.

