# CGPSspatialmatch Version 1
This is an R package, CGPSspatialmatch

Manuscript: Adjustment for Unmeasured Spatial Confounding in Settings of Continuous Exposure Conditional on the Binary Exposure Status: Conditional Generalized Propensity Score-Based Spatial Matching

Authors: Honghyok Kim, Michelle Bell. Pre-print available at https://arxiv.org/abs/2202.00814; (the manuscript is under revision for a peer-reviewed publication)

To install this in R, devtools::install_github('HonghyokKim/CGPSspatialmatch')


You may want to use SampleCode.R in the "SampleCode" folder. Notes are appended for each of the lines.

An example dataset is provided. You can load this by "data(StrokePPR_CT)", which is in SampleCode.R.

This dataset was created by sampling the original dataset. As such, this sample code will not give you the same results of the paper.


This package includes a function, "estimateATT", that enables users to easily obtain effect estimates in settings of non-time-to-failure data (e.g., cross-sectional data) by automatically performing distance-matching, PS estimation, CGPS estimation, GPS estimation, GPS matching, and disease model-fitting.


For other designs such as longitudinal design and cohort design with survival analyses, users may have to

1) transform time-to-failure data to non-time-to-failure data.

2) use functions for distance-matching, PS estimation, CGPS estimation, GPS estimation, and GPS matching in this package.

"estimateATT" function  will be developed to include longitudinal and survival data with Cox models and mixed effect models in the future.

Error reporting and questions: 
honghyok.kim@yale.edu
https://hkimresearch.com


