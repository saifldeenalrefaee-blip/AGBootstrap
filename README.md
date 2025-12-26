# AGBootstrap

## Overview
This package implements the **Adaptive Generative Bootstrap (AGB)** method for small sample inference. It is designed to correct coverage probabilities in multimodal or small datasets.

## Installation
You can install the development version from GitHub using:

```r
# install.packages("devtools")
devtools::install_github("saifldeenalrefaee-blip/AGBootstrap")
```

## Example
Here is a basic example of how to use the package:

```r
library(AGBootstrap)

# 1. Define sample data
data <- c(576, 635, 558, 578, 666, 580, 555, 661, 651, 605, 653, 575, 545, 572, 594)

# 2. Run the Adaptive Generative Bootstrap
result <- agb(data)

# 3. Print the results
print(result$ci)
```

## License
This package is licensed under the MIT License.
