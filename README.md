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

The variable of interest $Y$ is the body fat percentage measured by \textit{dxa}(dual-energy X-ray absorptiometry, variable \texttt{DXDTOPF}). The demographics data include \textit{age} (variable \texttt{RIDAGEYR}, $8 \le \mathit{age} \le 59$) and \textit{gender} (variable \texttt{RIAGENDR}, 1 for male and 2 for female), and the examination data include \textit{bmi} (body mass index, variable \texttt{BMXBMI}). There are $n = 5,055$ subjects and 1,411 (27.9\%) of which have missing $Y$ values. Figure \ref{fig:hist_dxa} in the SM shows that the distribution of \textit{dxa} varies on \textit{gender}.

