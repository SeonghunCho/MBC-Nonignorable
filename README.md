# Multiple Bias Calibration for Valid Statistical Inference under Nonignorable Nonresponse

# Author Contributions Checklist Form

## Data

### Abstract

The dataset used for the analyses in Section 7 of the paper is from the
National Health and Nutrition Examination Survey (NHANES) 2017–2018,
a program of studies designed to assess the health and nutritional status
of adults and children in the United States.

### Availability 

The dataset is available for download at <www.cdc.gov/nchs/nhanes.htm>.

### Description

The variable of interest $Y$ is the body fat percentage measured by *dxa*
(dual-energy X-ray absorptiometry, variable `DXDTOPF`). The demographics
data include *age*(variable `RIDAGEYR`, $8 \le \mathit{age} \le 59$) and
*gender*(variable `RIAGENDR`, 1 for male and 2 for female), and the
examination data include *bmi*(body mass index, variable `BMXBMI`). There
are $n = 5,055$ subjects and 1,411 (27.9\%) of which have missing $Y$ values.

## Code

### Abstract

All of the simulation study and the case study for this paper were done in R.
Separate code files are provided to conduct two simulation studies and analysis
on real data sets for the case study.

### Description

All of the R scripts used in this paper are available in this project.

# Instructions for Use

## Reproducibility

1. The folder named `code` contains all of the R scripts.
2. The R script `functions_sim.R` contains all of the R functions used in the simulation studies.
3. The R script `Simulation1.R` contains the R code to perform the simulation study in Section 6 of the paper. Its results are stored in `output/res_sim1/`.
4. The R script `Simulation1_summary.R` contains the R code to make Table 1 and to draw Figure 1.
5. The R script `Simulation2.R` contains the R code to perform the simulation study in Section A.1 of the supplementary material. Its results are stored in `output/res_sim2/`.
6. The R script `Simulation2_summary.R` contains the R code to make Table S1 and to draw Figure S1.
7. The R script `functions_case.R` contains all of the R functions used in the case study.
8. The R script `CaseStudy.R` contains the R code to perform the case study in Section 7 of the paper. Its results are stored in `output/res_case/`.
9. The R script `CaseStudy_summary.R` contains the R code to make Table 2 and 3 and to draw Figure S2 and S3.
10. The folder named `output` contains all of the results from two simulations studies and the case study.
11. The subfolder named `res_sim1` of `output` contains the results from the first simulation study.
12. The subfolder named `res_sim2` of `output` contains the results from the second simulation study.
13. The subfolder named `res_case` of `output` contains the results from the case study.
14. The subfolder named `summary` of `output` contains all of the tables and the figures presented in the paper.
