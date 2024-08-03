# MultiarmSACE
This repository offers sample R code for the methods proposed in the article “Doubly Robust Estimation and Sensitivity Analysis with Outcomes Truncated by Death in Multi-arm Clinical Trials” by Tong et al.

ntp.rdata: The dataset used for analysis. 

SampleCodeFun.R contains R functions needed to run the code in SampleCode.R. Our sandwich variance estimates are calculated using the geex package. Please refer to their documentation for more details and examples.

SampleCode.R contains R code implementing the proposed methods.

SampleCodeFunPI.R contains R functions needed to run the code in SampleCode.R for sensitivity analysis when principal ignorability assumption is violated. 

SampleCodePI.R contains R code implementing the proposed methods for sensitivity analysis when principal ignorability assumption is violated.
