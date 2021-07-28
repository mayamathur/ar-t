
# Overview

This repository contains all data, materials, and code required to reproduce all analyses in the paper:

Mathur MB & VanderWeele TJ (2022). Methods to address confounding and other biases in meta-analyses: Review and recommendations. *Annual Review of Public Health*. In press.

# Applied example

You can reproduce the analyses for the applied example by running the code `analyze.R` in the directory **Applied example/Code**. That file calls helper functions in `helper.R` and uses the dataset `kodama_prepped.csv`.

# Empirical benchmarks on agreement between RS and NRS

You can reproduce the re-analyses of the meta-meta-analyses by first running the code `prep.R` and then `analyze.R` in the directory **Empirical estimates of confounding bias**. Those files calls helper functions in `helper.R` and uses the datasets in **Raw data**.