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

1. The subfolder named ``code''

1. The subfolder named “Figure1_rejectionBoundary” contains the R code to draw
the rejection boundaries as shown in Figure 1.
2. The subfolder named “Figure2_Rotation” contains the R code to draw Figure 2.
3. The subfolder named “Figure3_Table5_Table6” contains the R code to draw
Figure 3 and generate the numbers for Table 5 and Table 6.
4. The subfolder named “Merge_analysis” contains the R code to merge all the raw
GWAS summary statistics for the eight phenotypes used in the paper. We merged
those GWAS summary statistics by using their chromosome: position column, and then
performed PC based tests to the Z-scores of each phenotype across all the SNPs in the
merged data. The output file is a csv file which contains the p-values for each PC based
method for each SNP.
5. The subfolder named “FigureS2-FigureS8_QQplotCode” contains the R code to
draw QQ plots and compute genomic control factors using the output csv file.
6. The subfolder named “FigureS1” contains the R code to draw Figure S1.
7. The subfolder named “SimulationSize” contains the R code to perform simulation
studies for the type I error rates, and the output was summarized in Table 1.
8. The subfolder named “SimulationPower” contains the R code to perform
simulation studies for the power comparisons, and the results were summarized in
Table 2 (simulation setup) and Table 3 (simulation results).





