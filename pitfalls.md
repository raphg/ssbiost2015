# Fifth Seattle Symposium in Biostatistics: Promises and pitfalls of omics experiments
Raphael Gottardo  
Friday, November 6  




## Motivation

- New high throughput technologies including gene and protein expression microarrays, flow and mass cytometry, next generation sequencing, etc.

- Experiments and protocols have become increasingly complex &rarr; can be very sensitive to specific settings

- Generate large datasets
    - Unique challenges (large p, small n)
    - Analyses of such data often require the use of many specialized tools

- Many challenges for data generation to analysis

## Retraction, retraction, retraction ...

<iframe src="http://retractionwatch.com/category/by-author/anil-potti-retractions/"></iframe>

## Confounding?

<iframe src="http://f1000research.com/articles/4-121/v1"></iframe>

## Outline

- Discuss some of the pitfalls encountered in “omics” data generation and analysis
    - Many of these are not new but have been amplified with omics data due to technical variation, high dimensionality, etc.
- Presents some strategies and guidelines that researchers should follow when undertaking such studies
- Use case: Gene expression analysis (microarray and RNA-seq)
- Example code will always be presented or made available (in the Rmd file)

## Short review on technologies

- Review Molecular [Biology 101 slides](./Biology_basics.html)
- Review [Microarray slides](./Microarrays.html) and [RNA-seq slides](https://www.biostat.wisc.edu/bmi776/lectures/rnaseq.pdf)

# Pitfall 1: Multiple testing!

## Differential gene expression

Let's assume that we are working with gene expression microarray data that have been normalized (and probes summarized).


Condition | 1 | -- | 1 | 2 | -- | 2 |
---|---|---|---|---|---|--- |
Replicate | 1 | -- | $R_1$ | 1 | -- | $R_2$ | 
Gene 1 | x | -- | x | y | -- | y |
Gene 2 | x | -- | x | y | -- | y |
Gene G | x | -- | x | y | -- | y |


 Our goal here is to find genes that are _differentially expressed_ between the two conditions. 

## Goal


**Note:** Here I will focus on oligo based arrays

- For each gene: Is gene $g$ differentially expressed between the two conditions?
- Is the mean expression level under condition 1 different from the mean expression level under condition 2?
- Test an hypothesis about the equality of the means of the two distributions

## Two-sample t-tests


- For each gene, are the mean log expression values equal?

Welch's t-test: $$t_g=(\bar{y}_{1g}-\bar{y}_{2g})/\sqrt{s^2_{1g}/R_1+s^2_{2g}/R_2}$$

If the means are equal, $t$ approximately follows a t-distribution with $R_1 + R_2 - 1$ degrees of freedom.

p-value $p = 2\cdot P(|t_{R_1+R_2-1}| > |t_g|)$

## Error rates

Called/Truth | Differentially expressed | Not differentially expressed|
---|---|--- |
Differentially expressed | 1-$\beta$ | $\alpha$ (Type I error)|
Not differentially expressed | $\beta$ (Type II error) | 1-$\alpha$|


## Multiple testing

- Fix the type I error rate (0.05)
- Minimize the type II
- This is what we do for each gene with a p-value cut off of 0.05
- Problem?
    - Look at many genes!

## Multiple testing

- 1000 t-tests, all null hypothesis are true ($\mu_1=\mu_2$) 
    - For one test, Pr of making an error is 0.05.
    - For 1000 tests, Pr of making at least one error is 1-(1-0.05)^1000 which is 1!

## Multiple testing


- The error probability is much greater when looking at many genes!
- We look at $G$ genes with $G$ very large! For each gene, $\alpha$ error probability
- Multiple testing is used to control an overall measure of error (FWER, FDR)

## Family Wise Error Rate


Controls the probability of making at least one type I error 


**Example:** Bonferroni multiple adjustment

$$\tilde p_g = G \cdot p_g$$

If $\tilde p_g \le \alpha$ then $FWER \le \alpha$


Many other (more powerful) FWER procedures exist (Holm's step-down, Hochberg's step-up).

## False Discovery Rate

Proportion of false positive among the genes called DE

First procedure introduced by Benjamini and Hochberg (1995)

- Order the p-values $p_{(1)} \le \dots \le p_{(g)} \le \dots \le p_{(G)}$

Let $k$ be the largest $g$ such that $p_{(g)} \le g/G\alpha$ then the 

FDR is controlled at $\alpha$

- Hypothesis need to be independent!
- Alternative approaches exist for dealing with the dependence at the cost of losing 
some power when the test are in fact independent. 

## False Discovery Rate


```r
library(ggplot2)
p <- rbeta(1000, 0.1, 0.1)
p_sorted <- sort(p)
qplot(1:1000, p_sorted) + geom_abline(intercept = 0, slope = 0.05/1000) + xlim(c(0, 500)) + ylim(c(0, 0.1))
```

```
## Warning: Removed 597 rows containing missing values (geom_point).
```

![](pitfalls_files/figure-html/unnamed-chunk-2-1.png) 

- Look at the `p.adjust` function in R!

## Recommended additional reading

- Benjamini, Yoav, and Daniel Yekutieli. 2001. “The Control of the False Discovery Rate in Multiple Testing under Dependency.” Annals of Statistics 29 (4). Institute of Mathematical Statistics: 1165–88.

- Gavrilov, Yulia, Yoav Benjamini, and Sanat K. Sarkar. 2009. “An Adaptive Step-down Procedure with Proven FDR Control under Independence.” Annals of Statistics 37 (2). Institute of Mathematical Statistics: 619–29.

- Storey, John D. 2002. “A Direct Approach to False Discovery Rates.” Journal of the Royal Statistical Society. Series B, Statistical Methodology 64 (3). Blackwell Publishers: 479–98.


# Pitfall 2: Small sample size!

## t-test - revisited

Microarray (and generally speaking) omics experiments are expensive, and as such the number of replicates is usually small. These can lead to the following issues:

- The test statistic is not normally distributed
- The variance estimates are noisy, with thousands of test, some of the estimated variances can be extremely small!


## Modified t-test


- Small variance problem: 

Regularized t-test: $$t_g=(\bar{y}_{1g}-\bar{y}_{2g})/\sqrt{s^2_{1g}/R_1+s^2_{2g}/R_2+c}$$

where $c$ is a positive constant used to regularize the variance estimate (e.g. 95% of all standard deviations $S_g$)

- Small sample size and distributional assumption:

Estimate the null distribution by permutation. Under the assumption of no differential expression we can permute the columns of the data matrix. 

This is the main idea behind SAM. 

Tusher, V. G., Tibshirani, R., & Chu, G. (2001). Significance analysis of microarrays applied to the ionizing radiation response. Proceedings of the National Academy of Sciences of the United States of America, 98(9), 5116–5121. doi:10.1073/pnas.091062498


## Linear Models for Microarray Data - LIMMA


LIMMA is popular Bioconductor package for the analysis of microarray data that provides a flexible linear modeling framework for assessing differential expression.

Smyth, G. K. (2004). Linear models and empirical bayes methods for assessing differential expression in microarray experiments. Statistical Applications in Genetics and Molecular Biology, 3(1), Article3. doi:10.2202/1544-6115.1027

## Linear Models for Microarray Data - LIMMA 


Let $\mathbf{y}^T_g=(y_{g1}, y_{g2},\dots, y_{gn})$ be the expression vector of gene $g$. 
The response would be log-ratios for 2-color microarrays or log-intensities for single-channel data. 

Smyth (2004) assumes that the mean expression value of $\mathbf{y}_g$ can be described through a linear model:

$$ \mathbb{E}(\mathbf{y}_g)=\mathbf{X}\boldsymbol{\alpha}_g$$ 

where $\mathbf{X}$ is a full column rank design matrix and $\boldsymbol{\alpha}_g$ is a coefficient vector. 

## Linear Models for Microarray Data - LIMMA 

It is futher assume that 

$$\mathrm{var}(\mathbf{y}_g)=\mathbf{W}_g\sigma_g^2$$ 

where $\mathbf{W}_g$ is a weight matrix. Certain constrasts $\mathbf{C}$ of the vector $\boldsymbol{\alpha}_g$ are of biological interest, and these are defined as:

$$\boldsymbol{\beta}_g=\mathbf{C}^T\boldsymbol{\alpha}_g.$$

The goal here will be to test if some of these contrats are equal to zero. 

## Linear Models for Microarray Data - LIMMA 

As an example of the difference between the design matrix $\mathbf{X}$ and the contrast matrix $\mathbf{C}$ consider a time course experiment for times $t_0$, $t_1$ and $t_2$ in which there are two replicates for each time point. A design matrix for this experiment would be:

$$\mathbf{X} = \left(\begin{array}{ccc}
        1 & 0 & 0\\
        1 & 0 & 0\\
        0 & 1 & 0\\
        0 & 1 & 0\\
        0 & 0 & 1\\
        0 & 0 & 1\end{array}\right)$$
        
If we are interested in the difference between $t_0$ and $t_1$, as well as the difference between $t_1$ and $t_2$, the transpose of the contrast matrix would be:

$$\mathbf{C^T} = \left(\begin{array}{ccc}
              -1 & 1 & 0 \\
               0 & -1 & 1 \end{array}\right)$$


## Linear Models for Microarray Data - LIMMA 


It is assumed that the linear model is fitted for each gene to obtain an estimator $\hat{\boldsymbol{\alpha}}_g$ of $\boldsymbol{\alpha}_g$, and estimator $s^2_g$ of $\sigma_g^2$ and estimated covariance matrix:

$$\mathrm{var}(\hat{\alpha}_g ) = \mathbf{V}_g s^2_g$$

where $\mathbf{V}_g$ is positive definite and does not depend on $s_g^2$. Then we have
$$ \mathrm{var}( \hat{\boldsymbol{\beta}}_g ) = \mathbf{C}^T \mathbf{V}_g \mathbf{C} s^2_g $$

**Note:** So far no distributional assumptions are made, and the fitting is not necessarily done by least-squares. However the contrast estimator will be assumed to be approximately normal with mean $\boldsymbol{\beta}_g$ and covariance $\mathbf{C}^T \mathbf{V}_g \mathbf{C} \sigma^2_g$


## Linear Models for Microarray Data - LIMMA 


Let $v_{gj}$ be the $j^{th}$ diagonal element of $\mathbf{C}^T \mathbf{V}_g \mathbf{C}$, then the distributional assumptions made are equivalent to: 
 
$$ \hat{\beta}_{gj} | \beta_{gj} , \sigma_g^2 \sim \mathrm{N}(\beta_{gj} , v_{gj} \sigma_g^2)$$
 
 and 
 
 $$ s^2_g|\sigma_g^2 \sim \frac{\sigma_g^2}{d_g}\chi^2_{d_g}$$
 
 where $d_g$ is the residual degrees of freedom for the linear model for gene $g$. Under these
assumptions the ordinary t-statistic

$$ t_{gj}=\frac{\hat{\beta}_{gj}}{s_g\sqrt{v_{gj}}}$$

follows an approximate t-distribution on $d_g$ degrees of freedom, which can be used to test the null hypothesis:
$H_0 : \beta_{gj} = 0.$

## Linear Models for Microarray Data - LIMMA 


This approach still suffers from the small sample size problem mentioned previously. One solution is to use a hierarchical model to borrow strength across genes. In particular, we will place a prior distribution on the inverse variances as follows:

$$\sigma_g^{-2}\sim \frac{1}{d_0 s_0^2}\chi^2_{d_0}$$ 

where ${d_0}$ and $s_0$ are fixed hyper-parameters. Similarly, a prior distribution can be placed on the unknown contrast parameters:

$$\beta_{gj} |\sigma_g^2,\beta_{gj}\ne 0 \sim \mathrm{N}(0,v_{0j}\sigma_g^2)$$

with $\mathrm{Pr}(\beta_{gj}\ne 0) = p_j$ where $v_{0j}$ and $p_j$ are given hyperparameters.

## LIMMA - Hierarchical model


Under the above hierarchical model, the posterior mean of $\sigma_g^2$ given $s_g^2$ is

$$\tilde{s}^2_g=\mathbb{E}(\sigma^2_g|s^2_g)= \frac{d_0s_0^2+d_gs_g^2}{d_0+d_g}$$

and we can define the following moderated t-statistics:

$$ \tilde{t}_{gj}=\frac{\hat{\beta}_{gj}}{\tilde{s}_g\sqrt{v_{gj}}}$$

The moderated t-statistics $\tilde{t}_{gj}$ and the residual sample variances $s^2$
are shown to be distributed independently. The moderated t is shown to follow a t-
distribution under the null hypothesis $H_0 : \beta_{gj}$ = 0 with degrees of freedom $d_g +d_0$.

## LIMMA - Estimation of hyperparameters


All parameters except $p_j$ are shared across genes and can easily be estimated using an empirical Bayes approach using all genes. The most difficult parameter to estimate is $p_j$, but this parameter is only used in the calculation of the posterior odds and is not required for inference via the moderated t-statistics. 

This is typical of frequentist inference where the alternative does not matter. 

## Other Bayesian approaches


Here are a few other Bayesian approaches that are available for the analysis of gene expression microarray data:
- Kendziorski, C. M., Newton, M. A., Lan, H., & Gould, M. N. (2003). On parametric empirical Bayes methods for comparing multiple groups using replicated gene expression profiles. Statistics in Medicine, 22(24), 3899–3914. doi:10.1002/sim.1548

- Gottardo, R., Raftery, A. E., Yeung, K. Y., & Bumgarner, R. E. (2006). Bayesian robust inference for differential gene expression in microarrays with multiple samples. Biometrics, 62(1), 10–18. doi:10.1111/j.1541-0420.2005.00397.x

- Lewin, A., Bochkina, N., & Richardson, S. (2007). Fully Bayesian mixture model for differential gene expression: simulations and model checks. Statistical Applications in Genetics and Molecular Biology, 6(1), Article36. doi:10.2202/1544-6115.1314

However, in my opinion, LIMMA provides the best user experience in terms of analysis in R and Bioconductor.

## The LIMMA package

Let's first install Limma:


```r
# You can skip this, if it's already installed on your machine
source("http://bioconductor.org/biocLite.R")
biocLite("limma")
```

Now we're ready to start using Limma


```r
library(limma)
library(Biobase)
library(data.table)
```

but we need some data!


## Getting some data with GEOquery


We're going to look at the dataset used in:

Nakaya, H. I., Wrammert, J., Lee, E. K., Racioppi, L., Marie-Kunze, S., Haining, W. N., et al. (2011). Systems biology of vaccination for seasonal influenza in humans. Nature Immunology, 12(8), 786–795. doi:10.1038/ni.2067


```r
library(GEOquery)
# Download the mapping information and processed data main serie #gds[[1]] = LAIV/TIV 0809, gds[[2]] = FACS, gds[[3]] =
# TIV 0708
dir.create("Data/GEO/", recursive = TRUE)
```

```
## Warning in dir.create("Data/GEO/", recursive = TRUE): 'Data/GEO' already
## exists
```

```r
gds <- getGEO("GSE29619", destdir = "Data/GEO/")
```

```
## ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE29nnn/GSE29619/matrix/
## Found 3 file(s)
## GSE29619-GPL13158_series_matrix.txt.gz
## Using locally cached version: Data/GEO//GSE29619-GPL13158_series_matrix.txt.gz
## Using locally cached version of GPL13158 found here:
## Data/GEO//GPL13158.soft 
## GSE29619-GPL3921_series_matrix.txt.gz
## Using locally cached version: Data/GEO//GSE29619-GPL3921_series_matrix.txt.gz
## Using locally cached version of GPL3921 found here:
## Data/GEO//GPL3921.soft 
## GSE29619-GPL570_series_matrix.txt.gz
## Using locally cached version: Data/GEO//GSE29619-GPL570_series_matrix.txt.gz
## Using locally cached version of GPL570 found here:
## Data/GEO//GPL570.soft
```

## Getting some data with GEOquery

but before we can use this, we need to clean up the pData a bit (see code in .Rmd file by clicking on the pencil icon above, which will bring you to this slide in the .Rmd file). 


```r
### Sanitize data and metadata
gds_new <- gds
sanitize_pdata <- function(pd) {
    keepCols <- c("characteristics_ch1.1", "characteristics_ch1.2", "description", "supplementary_file")
    pd <- pd[, keepCols]
    colnames(pd) <- c("ptid", "time", "description", "filename")
    pd$ptid <- gsub(".*: ", "", pd$ptid)
    pd$time <- gsub(".*: ", "", pd$time)
    pd$time <- gsub("Day", "D", pd$time)
    pd$description <- gsub("(-\\w*){2}$", "", pd$description)
    pd$filename <- basename(as.character(pd$filename))
    pd$filename <- gsub(".CEL.gz", "", pd$filename)
    pd
}
pData(gds_new[[1]]) <- sanitize_pdata(pData(gds_new[[1]]))
pData(gds_new[[2]]) <- sanitize_pdata(pData(gds_new[[2]]))
pData(gds_new[[3]]) <- sanitize_pdata(pData(gds_new[[3]]))
```

## Model set-up and estimation

Let's create seperate `ExpressionSet`s for the datasets of interests.


```r
TIV_08 <- gds_new[[1]][, grepl("2008-TIV", pData(gds_new[[1]])$description)]
LAIV_08 <- gds_new[[1]][, grepl("2008-LAIV", pData(gds_new[[1]])$description)]
TIV_07 <- gds_new[[3]][, grepl("2007-TIV", pData(gds_new[[3]])$description)]
```

TIV_08, LAIV_08 and TIV_07 are expression sets containing data from three time points (variable name is "time", with values D0, D3 and D7), for several probes (i.e., of form GSMXXXX) and patients (variable name "ptid"). 

We then use the limma R package to identify genes that are differentially expressed at D3 and D7 compared to baseline for each study. 



```r
mm_TIV_08 <- model.matrix(~ptid + time, TIV_08)  # design matrix
fit_TIV_08 <- lmFit(TIV_08, mm_TIV_08)  #Fit linear model for each gene given a series of arrays
ebay_TIV_08 <- eBayes(fit_TIV_08)  # compute moderated t-statistics, moderated F-statistic, and log-odds of differential expression
```

## Testing specific hypothesis

Let's first look at the estimated coefficients


```r
colnames(fit_TIV_08$coef)
```

```
##  [1] "(Intercept)" "ptid2"       "ptid29"      "ptid3"       "ptid32"     
##  [6] "ptid35"      "ptid38"      "ptid39"      "ptid4"       "ptid42"     
## [11] "ptid43"      "ptid44"      "ptid46"      "ptid47"      "ptid48"     
## [16] "ptid51"      "ptid53"      "ptid63"      "ptid65"      "ptid68"     
## [21] "ptid70"      "ptid72"      "ptid73"      "ptid74"      "ptid78"     
## [26] "ptid80"      "ptid83"      "ptid85"      "timeD3"      "timeD7"
```
In this case, the design matrix contains 1's and 0's, indicating which patient and time point matches up to a given measurement in the vector, $\mathbf{Y}$. There is no column for timeD0, since it is the reference point. When both timeD3 and timeD7 are zero, than we know that the measurement is from timeD0. 

Now we can test specific hypotheses.

## Testing specific hypothesis 

Here we look for genes differentially expressed at day 3 and day 7 wrt baseline:


```r
# Test t3=t0
topT3 <- topTable(ebay_TIV_08, coef = "timeD3", number = Inf, sort.by = "none")
# Test t7=t0
topT7 <- topTable(ebay_TIV_08, coef = "timeD7", number = Inf, sort.by = "none")
```


`topTable()` extracts a table of the top-ranked genes from a linear model fit and outputs a `data.frame` with the following columns:


```r
colnames(topT7)
```

```
##  [1] "ID"                               "GB_ACC"                          
##  [3] "SPOT_ID"                          "Species.Scientific.Name"         
##  [5] "Annotation.Date"                  "Sequence.Type"                   
##  [7] "Sequence.Source"                  "Target.Description"              
##  [9] "Representative.Public.ID"         "Gene.Title"                      
## [11] "Gene.Symbol"                      "ENTREZ_GENE_ID"                  
## [13] "RefSeq.Transcript.ID"             "Gene.Ontology.Biological.Process"
## [15] "Gene.Ontology.Cellular.Component" "Gene.Ontology.Molecular.Function"
## [17] "logFC"                            "AveExpr"                         
## [19] "t"                                "P.Value"                         
## [21] "adj.P.Val"                        "B"
```
as you can see it contains information about the probes contained in the `ExpressionSet` as well as values calculated by LIMMA. 

## MA plot d7 vs d0



```r
lm7 <- rowMeans(exprs(TIV_08)[, grepl("D7", pData(TIV_08)$time)])
lm0 <- rowMeans(exprs(TIV_08)[, grepl("D0", pData(TIV_08)$time)])
M <- lm7 - lm0
A <- (lm7 + lm0)/2
```


```r
dt <- data.table(A, M, abs_t = abs(topT7$t), p = topT7$adj.P.Val)
ggplot(dt, aes(x = A, y = M, color = abs_t, shape = p < 0.01)) + geom_point() + geom_point(data = dt[p < 0.01], aes(x = A, 
    y = M), color = "red")
```

![](pitfalls_files/figure-html/unnamed-chunk-10-1.png) 

## MA plot d7 vs d0 

Let's compare to ordinary t-statistics


```r
# Ordinary t-statistic
ordinary_t <- fit_TIV_08$coef/fit_TIV_08$stdev.unscaled/fit_TIV_08$sigma
ordinary_t <- ordinary_t[, "timeD7"]
# p-values based on normal approx with BH fdr adjustment
ordinary_p <- p.adjust(2 * pnorm(abs(ordinary_t), lower.tail = FALSE), method = "BH")
```


```r
dt <- data.table(A, M, abs_t = abs(ordinary_t), p = ordinary_p)
ggplot(dt[is.finite(abs_t)], aes(x = A, y = M, color = abs_t, shape = p < 0.01)) + geom_point() + geom_point(data = dt[p < 
    0.01], aes(x = A, y = M), color = "red")
```

![](pitfalls_files/figure-html/unnamed-chunk-12-1.png) 

## Setting up your own contrast


Suppose you want to look at the difference between timeD7 and timeD3. We need to create a contrast matrix that will get this information from the design matrix. This can easily be done using the `makeContrats` function as follows,


```r
cont_matrix <- makeContrasts(timeD7 - timeD3, levels = mm_TIV_08)
```

```
## Warning in makeContrasts(timeD7 - timeD3, levels = mm_TIV_08): Renaming
## (Intercept) to Intercept
```

```r
fit2 <- contrasts.fit(fit_TIV_08, cont_matrix)
```

```
## Warning in contrasts.fit(fit_TIV_08, cont_matrix): row names of contrasts
## don't match col names of coefficients
```

```r
fit2 <- eBayes(fit2)
topTable(fit2, adjust.method = "fdr")
```

```
##                            ID   GB_ACC SPOT_ID Species.Scientific.Name
## 211430_PM_s_at 211430_PM_s_at   M87789                    Homo sapiens
## 215946_PM_x_at 215946_PM_x_at AL022324                    Homo sapiens
## 214669_PM_x_at 214669_PM_x_at BG485135                    Homo sapiens
## 213502_PM_x_at 213502_PM_x_at AA398569                    Homo sapiens
## 215379_PM_x_at 215379_PM_x_at AV698647                    Homo sapiens
## 215121_PM_x_at 215121_PM_x_at AA680302                    Homo sapiens
## 214677_PM_x_at 214677_PM_x_at   X57812                    Homo sapiens
## 213182_PM_x_at 213182_PM_x_at   R78668                    Homo sapiens
## 216576_PM_x_at 216576_PM_x_at AF103529                    Homo sapiens
## 209138_PM_x_at 209138_PM_x_at   M87790                    Homo sapiens
##                Annotation.Date      Sequence.Type Sequence.Source
## 211430_PM_s_at    Aug 20, 2010  Exemplar sequence         GenBank
## 215946_PM_x_at    Aug 20, 2010 Consensus sequence         GenBank
## 214669_PM_x_at    Aug 20, 2010 Consensus sequence         GenBank
## 213502_PM_x_at    Aug 20, 2010 Consensus sequence         GenBank
## 215379_PM_x_at    Aug 20, 2010 Consensus sequence         GenBank
## 215121_PM_x_at    Aug 20, 2010 Consensus sequence         GenBank
## 214677_PM_x_at    Aug 20, 2010 Consensus sequence         GenBank
## 213182_PM_x_at    Aug 20, 2010 Consensus sequence         GenBank
## 216576_PM_x_at    Aug 20, 2010 Consensus sequence         GenBank
## 209138_PM_x_at    Aug 20, 2010  Exemplar sequence         GenBank
##                                                                                                                                                                                                                                                                                                                                                                                                                                            Target.Description
## 211430_PM_s_at                                                                                                                                                  gb:M87789.1 /DB_XREF=gi:185361 /FEA=FLmRNA /CNT=1 /TID=Hs.300697.0 /TIER=FL /STK=0 /UG=Hs.300697 /LL=3502 /UG_GENE=IGHG3 /DEF=Human (hybridoma H210) anti-hepatitis A IgG variable region, constant region, complementarity-determining regions mRNA, complete cds. /PROD=IgG /FL=gb:M87789.1
## 215946_PM_x_at gb:AL022324 /DB_XREF=gi:3702433 /FEA=DNA /CNT=2 /TID=Hs.296552.1 /TIER=ConsEnd /STK=0 /UG=Hs.296552 /LL=3545 /UG_GENE=IGLL3 /UG_TITLE=immunoglobulin lambda-like polypeptide 3 /DEF=Human DNA sequence from clone CTA-246H3 on chromosome 22 Contains the gene for IGLL1 (immunoglobulin lambda-like polypeptide 1, pre-B-cell specific), a pseudogene similar to LRP5 (Lipoprotein Receptor Related Protein.), ESTs, Genomic markers (D22S...
## 214669_PM_x_at                                                                                                                                                                                                                   gb:BG485135 /DB_XREF=gi:13417414 /DB_XREF=602503756F1 /CLONE=IMAGE:4617445 /FEA=mRNA /CNT=101 /TID=Hs.325722.1 /TIER=ConsEnd /STK=0 /UG=Hs.325722 /LL=28875 /UG_GENE=IGKV3D-15 /UG_TITLE=immunoglobulin kappa variable 3D-15
## 213502_PM_x_at                                                                                                                                                                                                                         gb:AA398569 /DB_XREF=gi:2051678 /DB_XREF=zt73g04.s1 /CLONE=IMAGE:728022 /FEA=DNA /CNT=47 /TID=Hs.296552.0 /TIER=Stack /STK=35 /UG=Hs.296552 /LL=3545 /UG_GENE=IGLL3 /UG_TITLE=immunoglobulin lambda-like polypeptide 3
## 215379_PM_x_at                                                                                                                                                                                                                                     gb:AV698647 /DB_XREF=gi:10300618 /DB_XREF=AV698647 /CLONE=GKCBJC12 /FEA=mRNA /CNT=4 /TID=Hs.289110.4 /TIER=ConsEnd /STK=0 /UG=Hs.289110 /LL=28831 /UG_GENE=IGLJ3 /UG_TITLE=immunoglobulin lambda joining 3
## 215121_PM_x_at                                                                                                                                                                                                                                     gb:AA680302 /DB_XREF=gi:2656270 /DB_XREF=ac83d05.s1 /CLONE=IMAGE:869193 /FEA=mRNA /CNT=18 /TID=Hs.181125.2 /TIER=ConsEnd /STK=1 /UG=Hs.181125 /LL=3535 /UG_GENE=IGL@ /UG_TITLE=immunoglobulin lambda locus
## 214677_PM_x_at                                                                                                                                                                                  gb:X57812.1 /DB_XREF=gi:33723 /GEN=immunoglobulin lambda light chain /FEA=mRNA /CNT=199 /TID=Hs.289110.2 /TIER=ConsEnd /STK=0 /UG=Hs.289110 /LL=28831 /UG_TITLE=immunoglobulin lambda joining 3 /DEF=Human rearranged immunoglobulin lambda light chain mRNA.
## 213182_PM_x_at                                                                                                                                                                                                                 gb:R78668 /DB_XREF=gi:854949 /DB_XREF=yi74c04.r1 /CLONE=IMAGE:144966 /FEA=EST /CNT=286 /TID=Hs.106070.2 /TIER=ConsEnd /STK=0 /UG=Hs.106070 /LL=1028 /UG_GENE=CDKN1C /UG_TITLE=cyclin-dependent kinase inhibitor 1C (p57, Kip2)
## 216576_PM_x_at                                                                                                                                                                  gb:AF103529.1 /DB_XREF=gi:4378387 /FEA=mRNA /CNT=1 /TID=Hs.247910.0 /TIER=ConsEnd /STK=0 /UG=Hs.247910 /DEF=Homo sapiens isolate donor N clone N88K immunoglobulin kappa light chain variable region mRNA, partial cds. /PROD=immunoglobulin kappa light chain variableregion
## 209138_PM_x_at                                                                                         gb:M87790.1 /DB_XREF=gi:185363 /FEA=FLmRNA /CNT=660 /TID=Hs.181125.0 /TIER=FL+Stack /STK=584 /UG=Hs.181125 /LL=3535 /UG_GENE=IGL@ /DEF=Human (hybridoma H210) anti-hepatitis A immunoglobulin lambda chain variable region, constant region, complementarity-determining regions mRNA, complete cds. /PROD=immunoglobulin lambda-chain /FL=gb:M87790.1
##                Representative.Public.ID
## 211430_PM_s_at                   M87789
## 215946_PM_x_at                 AL022324
## 214669_PM_x_at                 BG485135
## 213502_PM_x_at                 AA398569
## 215379_PM_x_at                 AV698647
## 215121_PM_x_at                 AA680302
## 214677_PM_x_at                   X57812
## 213182_PM_x_at                   R78668
## 216576_PM_x_at                 AF103529
## 209138_PM_x_at                   M87790
##                                                                                                                                                                                                         Gene.Title
## 211430_PM_s_at immunoglobulin heavy locus /// immunoglobulin heavy constant gamma 1 (G1m marker) /// immunoglobulin heavy constant mu /// immunoglobulin heavy variable 4-31 /// hypothetical protein LOC100290146
## 215946_PM_x_at                                   immunoglobulin lambda-like polypeptide 1 /// immunoglobulin lambda-like polypeptide 3 /// glucuronidase, beta/immunoglobulin lambda-like polypeptide 1 pseudogene
## 214669_PM_x_at                                                                       immunoglobulin kappa locus /// immunoglobulin kappa constant /// immunoglobulin kappa variable 3-20 /// similar to hCG1686089
## 213502_PM_x_at                                                                                                                             glucuronidase, beta/immunoglobulin lambda-like polypeptide 1 pseudogene
## 215379_PM_x_at                                                                                                                                                                 immunoglobulin lambda variable 1-44
## 215121_PM_x_at                                                                                                                                 cyclosporin A transporter 1 /// immunoglobulin lambda variable 1-44
## 214677_PM_x_at                                                                                                                                                                         Immunoglobulin lambda locus
## 213182_PM_x_at                                                                                                                                                    cyclin-dependent kinase inhibitor 1C (p57, Kip2)
## 216576_PM_x_at                      immunoglobulin kappa locus /// immunoglobulin kappa constant /// similar to Ig kappa chain V-I region HK102 precursor /// similar to Ig kappa chain V-I region HK102 precursor
## 209138_PM_x_at                                                                                                                                                                         Immunoglobulin lambda locus
##                                                          Gene.Symbol
## 211430_PM_s_at IGH@ /// IGHG1 /// IGHM /// IGHV4-31 /// LOC100290146
## 215946_PM_x_at                          IGLL1 /// IGLL3 /// LOC91316
## 214669_PM_x_at           IGK@ /// IGKC /// IGKV3-20 /// LOC100291682
## 213502_PM_x_at                                              LOC91316
## 215379_PM_x_at                                              IGLV1-44
## 215121_PM_x_at                                    CYAT1 /// IGLV1-44
## 214677_PM_x_at                                                  IGL@
## 213182_PM_x_at                                                CDKN1C
## 216576_PM_x_at             IGK@ /// IGKC /// LOC652493 /// LOC652694
## 209138_PM_x_at                                                  IGL@
##                                                ENTREZ_GENE_ID
## 211430_PM_s_at 100290146 /// 28396 /// 3492 /// 3500 /// 3507
## 215946_PM_x_at                       3543 /// 91316 /// 91353
## 214669_PM_x_at         100291682 /// 28912 /// 3514 /// 50802
## 213502_PM_x_at                                          91316
## 215379_PM_x_at                                          28823
## 215121_PM_x_at                            100290481 /// 28823
## 214677_PM_x_at                                           3535
## 213182_PM_x_at                                           1028
## 216576_PM_x_at           3514 /// 50802 /// 652493 /// 652694
## 209138_PM_x_at                                           3535
##                                                                RefSeq.Transcript.ID
## 211430_PM_s_at                                        XM_001718220 /// XM_002347483
## 215946_PM_x_at NM_001013618 /// NM_020070 /// NM_152855 /// NR_024448 /// NR_029395
## 214669_PM_x_at                                                         XM_002345544
## 213502_PM_x_at                                                            NR_024448
## 215379_PM_x_at                                                                     
## 215121_PM_x_at                                                         XM_002348112
## 214677_PM_x_at                                                                     
## 213182_PM_x_at                          NM_000076 /// NM_001122630 /// NM_001122631
## 216576_PM_x_at                                           XM_001724425 /// XM_942302
## 209138_PM_x_at                                                                     
##                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              Gene.Ontology.Biological.Process
## 211430_PM_s_at                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 0006955 // immune response // non-traceable author statement /// 0018298 // protein-chromophore linkage // inferred from electronic annotation
## 215946_PM_x_at                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              0005975 // carbohydrate metabolic process // inferred from electronic annotation /// 0006955 // immune response // non-traceable author statement
## 214669_PM_x_at                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   0006955 // immune response // non-traceable author statement
## 213502_PM_x_at                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               0005975 // carbohydrate metabolic process // inferred from electronic annotation
## 215379_PM_x_at                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   0006955 // immune response // non-traceable author statement
## 215121_PM_x_at                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   0006955 // immune response // non-traceable author statement
## 214677_PM_x_at                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   0006955 // immune response // non-traceable author statement
## 213182_PM_x_at 0000079 // regulation of cyclin-dependent protein kinase activity // traceable author statement /// 0000080 // G1 phase of mitotic cell cycle // traceable author statement /// 0000122 // negative regulation of transcription from RNA polymerase II promoter // inferred from electronic annotation /// 0007049 // cell cycle // inferred from electronic annotation /// 0007050 // cell cycle arrest // inferred from electronic annotation /// 0007050 // cell cycle arrest // traceable author statement /// 0008285 // negative regulation of cell proliferation // traceable author statement /// 0030511 // positive regulation of transforming growth factor beta receptor signaling pathway // inferred from mutant phenotype /// 0032582 // negative regulation of gene-specific transcription // inferred from direct assay /// 0033673 // negative regulation of kinase activity // inferred from direct assay /// 0042326 // negative regulation of phosphorylation // inferred from direct assay /// 0042551 // neuron maturation // inferred from electronic annotation /// 0050680 // negative regulation of epithelial cell proliferation // inferred from mutant phenotype
## 216576_PM_x_at                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   0006955 // immune response // non-traceable author statement
## 209138_PM_x_at                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
##                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            Gene.Ontology.Cellular.Component
## 211430_PM_s_at 0005576 // extracellular region // inferred from electronic annotation /// 0005576 // extracellular region // non-traceable author statement /// 0005624 // membrane fraction // non-traceable author statement /// 0005886 // plasma membrane // inferred from electronic annotation /// 0005887 // integral to plasma membrane // non-traceable author statement /// 0016020 // membrane // inferred from electronic annotation /// 0016021 // integral to membrane // inferred from electronic annotation
## 215946_PM_x_at                                                                                                                                                                                                                                                                                                                                                                             0005576 // extracellular region // inferred from electronic annotation /// 0016020 // membrane // non-traceable author statement
## 214669_PM_x_at                                                                                                                                                                                                                                                                                                                                                                                        0005576 // extracellular region // not recorded /// 0005576 // extracellular region // non-traceable author statement
## 213502_PM_x_at                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
## 215379_PM_x_at                                                                                                                                                                                                                                                                                                                                                                                                                                            0005576 // extracellular region // non-traceable author statement
## 215121_PM_x_at                                                                                                                                                                                                                                                                                                                                                                                                                                            0005576 // extracellular region // non-traceable author statement
## 214677_PM_x_at                                                                                                                                                                                                                                                                                                                                                                                                                                            0005576 // extracellular region // non-traceable author statement
## 213182_PM_x_at                                                                                                                                                                                                                                                                                 0005634 // nucleus // inferred from direct assay /// 0005634 // nucleus // inferred from electronic annotation /// 0005730 // nucleolus // inferred from direct assay /// 0005737 // cytoplasm // inferred from direct assay
## 216576_PM_x_at                                                                                                                                                                                                                                                                                                                                                                                        0005576 // extracellular region // not recorded /// 0005576 // extracellular region // non-traceable author statement
## 209138_PM_x_at                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
##                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    Gene.Ontology.Molecular.Function
## 211430_PM_s_at                                                                                                                                                                                          0003823 // antigen binding // traceable author statement /// 0003823 // antigen binding // inferred from electronic annotation /// 0003823 // antigen binding // non-traceable author statement /// 0004872 // receptor activity // inferred from electronic annotation /// 0005515 // protein binding // inferred from physical interaction /// 0008270 // zinc ion binding // inferred from electronic annotation /// 0046872 // metal ion binding // inferred from electronic annotation
## 215946_PM_x_at                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               0004553 // hydrolase activity, hydrolyzing O-glycosyl compounds // inferred from electronic annotation
## 214669_PM_x_at                                                                                                                                                                                                                                                                                                                                                                                                                                                                              0003823 // antigen binding // non-traceable author statement /// 0003823 // antigen binding // inferred from electronic annotation /// 0005515 // protein binding // inferred from physical interaction
## 213502_PM_x_at                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               0004553 // hydrolase activity, hydrolyzing O-glycosyl compounds // inferred from electronic annotation
## 215379_PM_x_at                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   0003823 // antigen binding // inferred from electronic annotation /// 0003823 // antigen binding // non-traceable author statement
## 215121_PM_x_at                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   0003823 // antigen binding // inferred from electronic annotation /// 0003823 // antigen binding // non-traceable author statement
## 214677_PM_x_at                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   0003823 // antigen binding // inferred from electronic annotation /// 0003823 // antigen binding // non-traceable author statement
## 213182_PM_x_at 0004860 // protein kinase inhibitor activity // inferred from electronic annotation /// 0004861 // cyclin-dependent protein kinase inhibitor activity // inferred from electronic annotation /// 0004861 // cyclin-dependent protein kinase inhibitor activity // traceable author statement /// 0005515 // protein binding // inferred from physical interaction /// 0005515 // protein binding // inferred from electronic annotation /// 0016301 // kinase activity // inferred from electronic annotation /// 0016563 // transcription activator activity // inferred from genetic interaction /// 0016564 // transcription repressor activity // inferred from mutant phenotype
## 216576_PM_x_at                                                                                                                                                                                                                                                                                                                                                                                                                                                                              0003823 // antigen binding // non-traceable author statement /// 0003823 // antigen binding // inferred from electronic annotation /// 0005515 // protein binding // inferred from physical interaction
## 209138_PM_x_at                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
##                     logFC   AveExpr         t      P.Value    adj.P.Val
## 211430_PM_s_at  2.7674612  8.427690  8.970877 2.624795e-12 1.436157e-07
## 215946_PM_x_at  1.0410176  6.039067  6.587058 1.881578e-08 3.468005e-04
## 214669_PM_x_at  0.8635100 11.327501  6.533285 2.300599e-08 3.468005e-04
## 213502_PM_x_at  1.0441970  9.142830  6.410658 3.637000e-08 3.468005e-04
## 215379_PM_x_at  1.1873632  9.664309  6.401504 3.763385e-08 3.468005e-04
## 215121_PM_x_at  1.0060274 10.416652  6.387251 3.968904e-08 3.468005e-04
## 214677_PM_x_at  1.0212244 11.183803  6.340074 4.732293e-08 3.468005e-04
## 213182_PM_x_at -0.7084388  5.203281 -6.321544 5.070646e-08 3.468005e-04
## 216576_PM_x_at  1.2227860  6.012631  6.282835 5.857118e-08 3.560803e-04
## 209138_PM_x_at  1.1115976 10.681856  6.241708 6.826052e-08 3.734874e-04
##                        B
## 211430_PM_s_at 16.626413
## 215946_PM_x_at  8.873759
## 214669_PM_x_at  8.695609
## 213502_PM_x_at  8.289495
## 215379_PM_x_at  8.259188
## 215121_PM_x_at  8.212009
## 214677_PM_x_at  8.055878
## 213182_PM_x_at  7.994573
## 216576_PM_x_at  7.866542
## 209138_PM_x_at  7.730570
```

## Your turn!

Ok, let's try to repeat what we've done with the TIV07 cohort. 






## RNA-seq

As opposed to microarrays, sequencing leads to count data. In this case, one could try to model counts over genes (or possibly genomic intervals) using discrete distributions such as Poisson and negative binomial:

1. Marioni, J. C., Mason, C. E., Mane, S. M., Stephens, M. & Gilad, Y. RNA-seq: an assessment of technical reproducibility and comparison with gene expression arrays. Genome Res. 18, 1509-1517 (2008).

2. Anders, S. et al. Count-based differential expression analysis of RNA sequencing data using R and Bioconductor. Nat Protoc 8, 1765-1786 (2013).

3. Lund, S. P., Nettleton, D., McCarthy, D. J. & Smyth, G. K. Detecting differential expression in RNA-sequence data using quasilikelihood with shrunken dispersion estimates. Stat Appl Genet Mol Biol 11, (2012).

4. Anders, S. & Huber, W. Differential expression analysis for sequence count data. Genome Biol. 11, R106 (2010).

Here I will focus on a slighlty different approach describe in this paper:

7. Law, C. W., Chen, Y., Shi, W. & Smyth, G. K. Voom: precision weights unlock linear model analysis tools for RNA-seq read counts. Genome Biol. 15, R29 (2014).



## Using a normal approximation?


- An alternative to the Poisson and NB models would be to find a data transformation that would make the count data approximately normal. 

- Law et al. (2014) propose to use the $\log_2(\mathrm{cpm}(\cdot))$ transformation. 

- The cpm transformation accounts for sequencing depth variability while the $\log_2$ transformation makes the data more normal. 

- However, the mean-variance relationship is quadratic, and a log transformation is not going to remove this dependence. As a consequence, it would be innapropriate to use a normal linear model with constant variance (even gene-wise). 

## Mean-variance trend estimation via voom

Law et al. (2014) propose to estimate the mean-variance trend, and then to use the inverse of the estimated standard deviation for each observation as weight in LIMMA. 
This is done by the `voom` function in LIMMA. 

Basically, `voom` fits the linear model without weights and uses the residuals of the model to estimate the weights, which are then pass onto a weighted LIMMA call for linear modeling. 

Law et al. (2014) show that this approach: 

1. Control the type I error rate
2. Is powerful among the methods that control the type I error rate
3. Has good FDR control
4. Is faster!

The `voom`+`limma` approach can also be used for gene set analysis, which is difficult to do with count based methods. 

## RNA-seq example



```r
# You should make sure the directory Data/GEO exist
gd <- getGEO("GSE45735", destdir = "Data/GEO/")
```

```
## ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE45nnn/GSE45735/matrix/
## Found 1 file(s)
## GSE45735_series_matrix.txt.gz
## Using locally cached version: Data/GEO//GSE45735_series_matrix.txt.gz
## Using locally cached version of GPL10999 found here:
## Data/GEO//GPL10999.soft
```

```r
pd <- pData(gd[[1]])
getGEOSuppFiles("GSE45735", makeDirectory = FALSE, baseDir = "Data/GEO/")
```

```
## ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE45nnn/GSE45735/suppl/
```

```r
# Note the regular expression to grep file names
files <- list.files(path = "Data/GEO/", pattern = "GSE45735_T.*.gz", full.names = TRUE)
# Read in gzip-compressed, tab-delimited files
file_list <- lapply(files, read.table, sep = "\t", header = TRUE)
# Subset to only those rows where Gene contains only non-space characters This addresses problems with T14 file
# containing 28 invalid rows at end of file
file_list <- lapply(file_list, function(file_list) subset(file_list, grepl("^[^[:space:]]+$", Gene)))
# Remove duplicated rows
file_list_unique <- lapply(file_list, function(x) {
    x <- x[!duplicated(x$Gene), ]
    x <- x[order(x$Gene), ]
    rownames(x) <- x$Gene
    x[, -1]
})
# Take the intersection of all genes
gene_list <- Reduce(intersect, lapply(file_list_unique, rownames))
file_list_unique <- lapply(file_list_unique, "[", gene_list, )
matrix <- as.matrix(do.call(cbind, file_list_unique))
# Clean up the pData
pd_small <- pd[!grepl("T13_Day8", pd$title), ]
pd_small$Day <- sapply(strsplit(gsub(" \\[PBMC\\]", "", pd_small$title), "_"), "[", 2)
pd_small$subject <- sapply(strsplit(gsub(" \\[PBMC\\]", "", pd_small$title), "_"), "[", 1)
colnames(matrix) <- rownames(pd_small)
```

## RNA-seq example 

Note that raw data files for sequencing experiments are available from the SRA database, which can be queried using the SRAdb package:

```r
source("http://bioconductor.org/biocLite.R")
biocLite("SRAdb")
```

The resulting files are usually very large!


## Using Voom


Let's first create an eSet we can use:

```r
# Note that I add one to the count
new_set <- ExpressionSet(assayData = matrix)
pData(new_set) <- pd_small
```

we now need to set-up our design matrix to estimate our weights:


```r
design <- model.matrix(~subject + Day, new_set)
new_set_voom <- voom(new_set, design = design)
```


```r
lm <- lmFit(new_set_voom, design)
eb <- eBayes(lm)
# Look at the other time-points
topTable(eb, coef = "DayDay1", number = 5)
```

```
##            logFC  AveExpr        t      P.Value   adj.P.Val        B
## TRIM25 0.3201758 7.440715 6.611904 6.898126e-08 0.001513311 7.861678
## IRF1   0.6714831 7.847838 6.089150 3.701578e-07 0.003992599 6.168622
## IRF9   0.4594858 6.717816 5.842538 8.188077e-07 0.003992599 5.468091
## ASPHD2 0.5445809 3.961196 5.744828 1.121315e-06 0.004099903 5.421792
## GBP2   0.5390285 7.976877 5.809736 9.099733e-07 0.003992599 5.271233
```


# Pitfall 3: Interpretability

## Going from genes to gene sets

- So far we have seen how to use microarrays or RNA-seq to derive a list of signiﬁcantly differentially expressed genes, while controlling for false discovery. 

- Sometimes it can be convenient to look at biological pathways, or more generally genesets to gain biological insights. 

## Goals of GSEA

Detecting changes in gene expression datasets can be hard due to
- the large number of genes/probes,
- the high variability between samples, and
- the limited number of samples.

The goal of GSEA is to enable the detection of modest but coordinated
changes in prespeciﬁed sets of related genes. Such a set might include all the genes in a speciﬁc pathway,
for instance, or genes that have been shown to be coregulated based on previously published studies. 

Most of what I'll be discussing here will be based on this paper:

- Wu, D. & Smyth, G. Camera: a competitive gene set test accounting for inter-gene correlation. 40, e133-e133 (2012).

## Competitive vs self-contained analyses

As explained in Wu & Smyth (2012)

Two approaches can be used to test the significance of a gene set:
1) ‘self-contained’ gene set tests examine a set of genes in their own right without reference to other genes in the genome (or array)
2) ‘competitive’ gene set tests compare genes in the test set relative to all other genes.

Competitive tests focus more on distinguishing the most important biological processes from those that are less important. Competitive tests are overwhelmingly more commonly used in the genomic literature. In particular, this is the approach used in the GSEA paper:

- Subramanian, A. et al. Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. Proc. Natl. Acad. Sci. U.S.A. 102, 15545-15550 (2005).

## Accounting for within set correlation


Most competitive gene set tests assume independence of genes, because they evaluate P-values by permutation of gene labels. However, these tests can be sensitive to inter-gene correlations. 

In this lecture we will talk about geneset analysis using the approach of Wu and Smyth that accounts for inter-gene correlations. 

## CAMERA

Camera, for competitive gene set test accounting for inter-gene correlation, makes heavy used of the limma framework. The same linear model is assumed for the mean gene expression level, namely, 

$$ \mathbb{E}(\mathbf{y}_g)=\mathbf{X}\boldsymbol{\alpha}_g$$ 

and we will also write $\mathrm{cov}(y_{gi},y_{g'i})=\rho_{gg'}$
Note that this correlation is the residual treatment effect, once any treatment effect has been removed. 

As with `limma`, we assumed that a specific contrast is of interest:

$$\beta_g=\sum_{j=1}^p c_j \alpha_{gj}$$

and we wish to test $H_0: \beta_g=0$, which can be done using the moderated $t$ statistics, $\tilde{t}$ that follows a t-statistics with $d+d_0$ degrees of freedom. 


## CAMERA 


Then W&S define a normalized version $z_g =F^{-1}F_t(\tilde{t}_g)$

**What is the distribution of $z_g$?**

Consider a set of m genewise statistics $z_1,\dots , z_m$. The variance of the mean of the statistics is

$$ \mathrm{var}(\overline{z})=1/m^2(\sum_{i} \tau_i^2 + 2\sum_{i\lt j}\tau_i\tau_j)$$

where $\tau_i$ is the standard deviation of $z_i$ and the $\rho_{ij}$ are the pairwise correlations. We can rewrite it as follows when all $\tau_i$'s are equal

$$\mathrm{var}(\overline{z})=\tau^2/m \mathrm{VIF}$$ 

where VIF is the variance inflation factor $1+(m-1)\overline{\rho}$.

## CAMERA - Testing


Now the idea of competitive gene-set analysis can be done by comparing two mean set statistics $\overline{z}_1$ and $\overline{z}_2$. Where $z_1$ is our set of interest and $z_2$ is the set composed of all other genes. The main idea behind CAMERA is to form a test-statistics that incorporates dependence between the genes in the first set. 

This can simply be done by forming the following T-statistics

$$
T=(\overline{z}_1-\overline{z}_2)/(s_p\sqrt{VIF/m_1+1/m_2})
$$

Now we just need a way to estimate the inter-gene correlation. W&S also propose a
modified version of the wilcoxon-rank-sum test, which we won't discuss here.

## CAMERA - Estimating covariances


Let's consider the QR decomposition of the design matrix $X=QR$, where $R$ is upper triangular ($n\times p$) and $Q$ ($n\times n$) is such that $Q^TQ=I$. An $m\times n$ matrix of independent residuals is obtained by $U = YQ_2$, where $Q_2$ represents the trailing $d$ columns of $Q$. The matrix $U$ is already available as a by-product of fitting genewise linear models to the expression values using standard numerical algorithms, so extracting it requires no extra computation in `limma`. 
The residual standard error $s_g$ for gene $g$ is equal to the root mean square of the corresponding row of $U$. We standardize each row of $U$ by dividing by $s_g$.
At this point, we could obtain the correlation matrix for the $m$ genes from $C = UU^T$; however, this is a numerically inefficient procedure if $m$ is large.
A numerically superior algorithm is to compute the column means $\overline{u}_{\cdot k}$ of $U$. Then we can form the following estimate of the $VIF$
$$
\widehat{\mathrm{VIF}}=\frac{m}{d}\sum_d \overline{u}^2_{\cdot k}
$$

If $m$ and $d$ are both reasonably large, and $\overline{\rho}$ is relatively small, then VIF is approximately distributed as VIF $\chi^2_d/d$. This can be used to find an asymptotic distribution for our test statistics. 

## Using CAMERA

Camera is readily available in the `limma` package. Let us go back to our RNA-seq example:


```r
gd <- getGEO("GSE45735", destdir = "Data/GEO/")
```

```
## ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE45nnn/GSE45735/matrix/
## Found 1 file(s)
## GSE45735_series_matrix.txt.gz
## Using locally cached version: Data/GEO//GSE45735_series_matrix.txt.gz
## Using locally cached version of GPL10999 found here:
## Data/GEO//GPL10999.soft
```

```r
pd <- pData(gd[[1]])
getGEOSuppFiles("GSE45735", makeDirectory = FALSE, baseDir = "Data/GEO/")
```

```
## ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE45nnn/GSE45735/suppl/
```

```r
# Note the regular expression to grep file names
files <- list.files(path = "Data/GEO/", pattern = "GSE45735_T.*.gz", full.names = TRUE)
# Read in gzip-compressed, tab-delimited files
file_list <- lapply(files, read.table, sep = "\t", header = TRUE)
```


```r
knitr::knit_exit()
```


## Using CAMERA


```r
# Subset to only those rows where Gene contains only non-space characters This addresses problems with T14 file
# containing 28 invalid rows at end of file
file_list <- lapply(file_list, function(file_list) subset(file_list, grepl("^[^[:space:]]+$", Gene)))
# Remove duplicated rows
file_list_unique <- lapply(file_list, function(x) {
    x <- x[!duplicated(x$Gene), ]
    x <- x[order(x$Gene), ]
    rownames(x) <- x$Gene
    x[, -1]
})
# Take the intersection of all genes
gene_list <- Reduce(intersect, lapply(file_list_unique, rownames))
file_list_unique <- lapply(file_list_unique, "[", gene_list, )
matrix <- as.matrix(do.call(cbind, file_list_unique))
# Clean up the pData
pd_small <- pd[!grepl("T13_Day8", pd$title), ]
pd_small$Day <- sapply(strsplit(gsub(" \\[PBMC\\]", "", pd_small$title), "_"), "[", 2)
pd_small$subject <- sapply(strsplit(gsub(" \\[PBMC\\]", "", pd_small$title), "_"), "[", 1)
colnames(matrix) <- rownames(pd_small)
```


## CAMERA - Hands on



```r
# Note that I add one to the count
new_set <- ExpressionSet(assayData = matrix)
pData(new_set) <- pd_small
```

we now need to set-up our design matrix to estimate our weights:


```r
design <- model.matrix(~subject + Day, new_set)
new_set_voom <- voom(new_set, design = design)
```

## CAMERA - Hands on


```r
lm <- lmFit(new_set_voom, design)
eb <- eBayes(lm)
# Look at the other time-points
topTable(eb, coef = "DayDay1", number = 5)
```

```
##            logFC  AveExpr        t      P.Value   adj.P.Val        B
## TRIM25 0.3201758 7.440715 6.611904 6.898126e-08 0.001513311 7.861678
## IRF1   0.6714831 7.847838 6.089150 3.701578e-07 0.003992599 6.168622
## IRF9   0.4594858 6.717816 5.842538 8.188077e-07 0.003992599 5.468091
## ASPHD2 0.5445809 3.961196 5.744828 1.121315e-06 0.004099903 5.421792
## GBP2   0.5390285 7.976877 5.809736 9.099733e-07 0.003992599 5.271233
```

## MSigDB

The [Molecular Signatures Database](http://www.broadinstitute.org/gsea/msigdb/index.jsp) (MSigDB) is a collection of annotated gene sets for use with GSEA analysis. 

These gene sets are available from download as gmt files, and can be read into R using `GSEAbase`. 

Let's first download and install the package that we need

```r
library(BiocInstaller)
biocLite("GSEABase")
```
note that `camera` is available as part of `limma`, so there is nothing to install.

## Getting started with GSEA analyses

We load the `GSEAbase` package for loading gene sets.

```r
library(GSEABase)
```

```
## Loading required package: annotate
## Loading required package: AnnotationDbi
## Loading required package: stats4
## Loading required package: GenomeInfoDb
```

```
## Warning: package 'GenomeInfoDb' was built under R version 3.2.2
```

```
## Loading required package: S4Vectors
```

```
## Warning: package 'S4Vectors' was built under R version 3.2.2
```

```
## Creating a generic function for 'nchar' from package 'base' in package 'S4Vectors'
## Loading required package: IRanges
```

```
## Warning: package 'IRanges' was built under R version 3.2.2
```

```
## Warning in .recacheSubclasses(def@className, def, doSubclasses, env):
## undefined subclass "externalRefMethod" of class "expressionORfunction";
## definition not updated
```

```
## Warning in .recacheSubclasses(def@className, def, doSubclasses, env):
## undefined subclass "externalRefMethod" of class "functionORNULL";
## definition not updated
```

```
## 
## Attaching package: 'IRanges'
## 
## The following object is masked from 'package:data.table':
## 
##     shift
## 
## Loading required package: XML
## Loading required package: graph
## 
## Attaching package: 'graph'
## 
## The following object is masked from 'package:XML':
## 
##     addNode
```
and convert the gene sets to gene indices

```r
c2_set <- getGmt("GSEA-sets/c2.all.v4.0.symbols.gmt")
gene_ids <- geneIds(c2_set)
# Camera requires gene-indices sets_indices <- ids2indices(gene_ids, rownames(new_set))
sets_indices <- ids2indices(gene_ids, rownames(new_set))
```

## Finding enriched gene sets

As with `limma`, we need to specify the contrast we wish to test at the set level:

```r
# Note that camera works on voom objects
cont_matrix <- makeContrasts("DayDay1", levels = design)
```

```
## Warning in makeContrasts("DayDay1", levels = design): Renaming (Intercept)
## to Intercept
```

```r
res <- camera(new_set_voom, sets_indices, design = design, cont_matrix)
res[1:10, ]
```

```
##                                                                   NGenes
## REACTOME_TRAF6_MEDIATED_NFKB_ACTIVATION                               21
## BIOCARTA_IL10_PATHWAY                                                 17
## REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM                         261
## ST_STAT3_PATHWAY                                                      11
## PETROVA_PROX1_TARGETS_DN                                              64
## REACTOME_INNATE_IMMUNE_SYSTEM                                        256
## REACTOME_IMMUNE_SYSTEM                                               880
## REACTOME_GROWTH_HORMONE_RECEPTOR_SIGNALING                            24
## RASHI_RESPONSE_TO_IONIZING_RADIATION_2                               127
## REACTOME_RIG_I_MDA5_MEDIATED_INDUCTION_OF_IFN_ALPHA_BETA_PATHWAYS     72
##                                                                     Correlation
## REACTOME_TRAF6_MEDIATED_NFKB_ACTIVATION                           -0.0100124047
## BIOCARTA_IL10_PATHWAY                                              0.0178927806
## REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM                       0.0218890159
## ST_STAT3_PATHWAY                                                  -0.0411621804
## PETROVA_PROX1_TARGETS_DN                                          -0.0003186575
## REACTOME_INNATE_IMMUNE_SYSTEM                                      0.0107912531
## REACTOME_IMMUNE_SYSTEM                                             0.0072485056
## REACTOME_GROWTH_HORMONE_RECEPTOR_SIGNALING                         0.0047726371
## RASHI_RESPONSE_TO_IONIZING_RADIATION_2                             0.0073250352
## REACTOME_RIG_I_MDA5_MEDIATED_INDUCTION_OF_IFN_ALPHA_BETA_PATHWAYS  0.0008893810
##                                                                   Direction
## REACTOME_TRAF6_MEDIATED_NFKB_ACTIVATION                                  Up
## BIOCARTA_IL10_PATHWAY                                                    Up
## REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM                             Up
## ST_STAT3_PATHWAY                                                         Up
## PETROVA_PROX1_TARGETS_DN                                                 Up
## REACTOME_INNATE_IMMUNE_SYSTEM                                            Up
## REACTOME_IMMUNE_SYSTEM                                                   Up
## REACTOME_GROWTH_HORMONE_RECEPTOR_SIGNALING                               Up
## RASHI_RESPONSE_TO_IONIZING_RADIATION_2                                   Up
## REACTOME_RIG_I_MDA5_MEDIATED_INDUCTION_OF_IFN_ALPHA_BETA_PATHWAYS        Up
##                                                                         PValue
## REACTOME_TRAF6_MEDIATED_NFKB_ACTIVATION                           1.746976e-05
## BIOCARTA_IL10_PATHWAY                                             3.632593e-05
## REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM                      5.064116e-05
## ST_STAT3_PATHWAY                                                  5.463043e-05
## PETROVA_PROX1_TARGETS_DN                                          5.700425e-05
## REACTOME_INNATE_IMMUNE_SYSTEM                                     5.733056e-05
## REACTOME_IMMUNE_SYSTEM                                            5.973135e-05
## REACTOME_GROWTH_HORMONE_RECEPTOR_SIGNALING                        9.082360e-05
## RASHI_RESPONSE_TO_IONIZING_RADIATION_2                            9.542540e-05
## REACTOME_RIG_I_MDA5_MEDIATED_INDUCTION_OF_IFN_ALPHA_BETA_PATHWAYS 1.013240e-04
##                                                                          FDR
## REACTOME_TRAF6_MEDIATED_NFKB_ACTIVATION                           0.04029306
## BIOCARTA_IL10_PATHWAY                                             0.04029306
## REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM                      0.04029306
## ST_STAT3_PATHWAY                                                  0.04029306
## PETROVA_PROX1_TARGETS_DN                                          0.04029306
## REACTOME_INNATE_IMMUNE_SYSTEM                                     0.04029306
## REACTOME_IMMUNE_SYSTEM                                            0.04029306
## REACTOME_GROWTH_HORMONE_RECEPTOR_SIGNALING                        0.04378187
## RASHI_RESPONSE_TO_IONIZING_RADIATION_2                            0.04378187
## REACTOME_RIG_I_MDA5_MEDIATED_INDUCTION_OF_IFN_ALPHA_BETA_PATHWAYS 0.04378187
```


## Finding enriched gene sets over time


```r
res <- vector("list", length = 10)
for (i in 1:10) {
    contrast <- paste0("DayDay", i)
    cont_matrix <- makeContrasts(contrast, levels = design)
    res[[i]] <- camera(new_set_voom, sets_indices, design = design, contrast = cont_matrix, sort = FALSE)
}
```

```
## Warning in makeContrasts(contrast, levels = design): Renaming (Intercept)
## to Intercept
```

```
## Warning in makeContrasts(contrast, levels = design): Renaming (Intercept)
## to Intercept
```

```
## Warning in makeContrasts(contrast, levels = design): Renaming (Intercept)
## to Intercept
```

```
## Warning in makeContrasts(contrast, levels = design): Renaming (Intercept)
## to Intercept
```

```
## Warning in makeContrasts(contrast, levels = design): Renaming (Intercept)
## to Intercept
```

```
## Warning in makeContrasts(contrast, levels = design): Renaming (Intercept)
## to Intercept
```

```
## Warning in makeContrasts(contrast, levels = design): Renaming (Intercept)
## to Intercept
```

```
## Warning in makeContrasts(contrast, levels = design): Renaming (Intercept)
## to Intercept
```

```
## Warning in makeContrasts(contrast, levels = design): Renaming (Intercept)
## to Intercept
```

```
## Warning in makeContrasts(contrast, levels = design): Renaming (Intercept)
## to Intercept
```

## Visualizing the results


```r
library(pheatmap)
PValue <- sapply(res, function(x) {
    ifelse(x$Direction == "Up", -10 * log10(x$PValue), 10 * log10(x$PValue))
})
rownames(PValue) <- rownames(res[[1]])
PValue_max <- rowMax(abs(PValue))
PValue_small <- PValue[PValue_max > 30, ]
anno <- data.frame(Time = paste0("Day", 1:10))
rownames(anno) <- colnames(PValue_small) <- paste0("Day", 1:10)
```

## Visualizing the results

![](pitfalls_files/figure-html/unnamed-chunk-31-1.png) 


## Using non MSigDB gene_sets


Any genesets can be used for a GSEA analysis. For example, we can use the sets published in:

Li, S. et al. Molecular signatures of antibody responses derived from a systems biology study of five human vaccines. Nat. Immunol. 15, 195–204 (2013).


```r
BTM_set <- getGmt("GSEA-sets/BTM_for_GSEA_20131008.gmt")
gene_ids <- geneIds(BTM_set)
# Camera requires gene-indices sets_indices <- ids2indices(gene_ids, rownames(new_set))
sets_indices <- ids2indices(gene_ids, rownames(new_set))
```


```r
res <- vector("list", length = 10)
for (i in 1:10) {
    contrast <- paste0("DayDay", i)
    cont_matrix <- makeContrasts(contrast, levels = design)
    res[[i]] <- camera(new_set_voom, sets_indices, design = design, contrast = cont_matrix, sort = FALSE)
}
```

```
## Warning in makeContrasts(contrast, levels = design): Renaming (Intercept)
## to Intercept
```

```
## Warning in makeContrasts(contrast, levels = design): Renaming (Intercept)
## to Intercept
```

```
## Warning in makeContrasts(contrast, levels = design): Renaming (Intercept)
## to Intercept
```

```
## Warning in makeContrasts(contrast, levels = design): Renaming (Intercept)
## to Intercept
```

```
## Warning in makeContrasts(contrast, levels = design): Renaming (Intercept)
## to Intercept
```

```
## Warning in makeContrasts(contrast, levels = design): Renaming (Intercept)
## to Intercept
```

```
## Warning in makeContrasts(contrast, levels = design): Renaming (Intercept)
## to Intercept
```

```
## Warning in makeContrasts(contrast, levels = design): Renaming (Intercept)
## to Intercept
```

```
## Warning in makeContrasts(contrast, levels = design): Renaming (Intercept)
## to Intercept
```

```
## Warning in makeContrasts(contrast, levels = design): Renaming (Intercept)
## to Intercept
```

## Visualizing the results



```r
PValue <- sapply(res, function(x) {
    ifelse(x$Direction == "Up", -10 * log10(x$PValue), 10 * log10(x$PValue))
})
rownames(PValue) <- rownames(res[[1]])
PValue_max <- rowMax(abs(PValue))
PValue_small <- PValue[PValue_max > 30, ]
anno <- data.frame(Time = paste0("Day", 1:10))
rownames(anno) <- colnames(PValue_small) <- paste0("Day", 1:10)
```

## Visualizing the results


```r
pheatmap(PValue_small, cluster_cols = FALSE)
```

![](pitfalls_files/figure-html/unnamed-chunk-35-1.png) 


## Conclusion

- Gene sets provide a convenient way to summarize gene activities over known pathways, which facilitate biological interpretation and replicability!

- Gene set analysis are complementary to gene based analyses, as gene sets might masked gene level signal

- You are not limited to using the predeﬁned gene sets. 
    - Specific applications might require specific gene sets (e.g. Immune gene sets). 
    - The GSEA gmt format provides a convenient way to do this.

# Pitfall 3: Confounding and batch effects!

## Reading

Before we start, you should read the following papers:

1. Leek, J. T. & Storey, J. D. Capturing Heterogeneity in Gene Expression Studies by Surrogate Variable Analysis. PLoS Genet 3, e161 (2007).
2. Leek, J. T. & Storey, J. D. A general framework for multiple testing dependence. Proc. Natl. Acad. Sci. U.S.A. 105, 18718–18723 (2008).
3. Leek, J. T. et al. Tackling the widespread and critical impact of batch effects in high-throughput data. Nature Reviews Genetics 11, 733-739 (2010).
4. Johnson, W. E., Li, C. & Rabinovic, A. Adjusting batch effects in microarray expression data using empirical Bayes methods. Biostatistics 8, 118-127 (2007).
5. Gagnon-Bartsch, J. A. & Speed, T. P. Using control genes to correct for unwanted variation in microarray data. Biostatistics 13, 539-552 (2012).

## Motivation

Batch effects are technical sources of variation that have been added to the samples during handling.
Example of batch variables: lot number, technician, instrument settings, etc. 

If not adjusted for, these batch variables can have a substantial effects on downstream analysis


## Normalization and batch effects

Unfortunately, normalization will not always correct for batch effects. Technical variation due to batch effects might only affect a subset of the genes.

<img src="http://www.nature.com/nrg/journal/v11/n10/images/nrg2825-f1.jpg" width=300>

"For a published bladder cancer microarray data set obtained using an Affymetrix platform, we obtained the raw data for only the normal samples. Here, green and orange represent two different processing dates."

## Adjusting for batch effects

**Two scenarios:**

1. You have information about the batch variable
    - Use your batch effect as a covariate in your analysis (e.g. limma)
2. You suspect a batch effect, but you don't know where it is coming from
    - The batch effect needs to be estimated first and then corrected for, by adding the estimated variables as co-variates


## Singular value decomposition

Let $X$ be a matrix of size $m\times n$ ($m \ge n$) and rank $r\le n$
then we can decompose $X$ as 

$$X=USV^T$$

- U is the matrix of left singular vectors (eigenassays)
- V is the matrix of right singular vectors (eigengenes)
- S is the matrix of singular values (eigenvalues)

$U^TU=VV^T=I$ (orthogonal vectors)


$S=diag(s_1, \dots, s_n)$ where $s_l\ge 0$ and $s_{r+1}=\dots=s_n=0$

$X_i=\sum_j u_{ij}s_j\mathbf{v}_j$, which can be interpreted as a change of coordinate

## Relationship to principal component analysis

$X=USV^T$, and we have $X^TX=VSU^TUSV^T=VS^2V^T$

What happens if the rows of X are scaled?


## Surrogate variable analysis

Let $X_{m\times n}=(x_1,..,x_m)^T$ be the matrix of normalized expression values, with $n$ arrays and $m$ genes. 
Let $y=(y_1,..,y_n)^T$ be a vector of length $n$ representing the primary variable of interest (e.g covariates, vector of unknown coefficients). Without loss of generality model $x_{ij}=\mu_i+f_i( y_j) + e_{ij}$, where $\mu_i$ is the baseline level of expression, $f_i(y_j)=\mathbb{E}(x_{ij} | y_j)-\mu_i$ gives the relationship between measured variable of interest and gene $i$, and $e_{ij}$ is random noise with mean zero.

Suppose there are $L$ biologically meaningful unmodeled factors, such as age, environmental exposure, genotype,
etc. Let $g_l = (g_{l1},...,g_{ln})$ be an _arbitrarily complicated function_ of
the lth factor across all $n$ arrays, for $l=1,2,...,L$. Our model becomes:

$$x_{ij}=\mu_i + f_i(y_j) +\sum_{l=1}^L \gamma_{l_i}g_{l_j} + e^*_{ij}$$

and if factor $l$ does not influence the expression of gene $i$, we have $\gamma_{l_i}=0$.

## Surrogate variable analysis

In practice it is impossible to estimate $\sum_{l=1}^L \gamma_{l_i}g_{l_j}$, so Leek and Storey propose to use singular value decomposition to approximate the matrix $(\sum_{l=1}^L \gamma_{l_i}g_{l_j})_{ij}$ by its singular value decomposition. Computationally, this is done in two steps:


1. Detect unmodeled factors
2. Construct surrogate variables

## Detect unmodel factor
  
The main idea is as follows:

- Compute the residual matrix $r_{ij} = x_{ij}- \hat{\mu}_i - \hat{f}_i(y_j)$
- Perform the SVD of $R=(r_{ij})$
- Permute the rows of the matrix $R$ to obtain $R^*$. Regress  $r^*_{ij} = x_{ij}- \hat{\mu}_i - \hat{f}_i(y_j)$ to get residual matrix $R_0$, and perform the SVD of $R_0$. Repeat this many times to generate a null distribution for the residuals, given that $y$ is accounted for. 
- Compare the observed eigenvalues to those generated from the null distribution to obtain significance p-values
- Record the $K$ significant variables

## Construct surrogate variables 

1. Compute the residual matrix $r_{ij} = x_{ij} - \hat{\mu}_i - \hat{f}_i(y_j)$
2. Perform the SVD of $R=(r_{ij})$

Let $e_k=(e_{k1},...,e_{kn})^T$ be the $k$-th column of $V$ (for $k=1,...,n$). These $e_k$ are the residual eigengenes and represent orthogonal residual signals independent of the signal due to the primary variable.

2. Regress $e_k$ on the $x_i$ to access the significance of the $k$-th factor on gene $i$
3. Use the selected genes to form a reduced expression matrix and repeat 1. The estimated factor will form the basis for the surrogate variables
4. In any subsequent analysis include these factors in your model

## Using the sva package

We're going to look at the dataset used in:

Nakaya, H. I., Wrammert, J., Lee, E. K., Racioppi, L., Marie-Kunze, S., Haining, W. N., et al. (2011). Systems biology of vaccination for seasonal influenza in humans. Nature Immunology, 12(8), 786795. doi:10.1038/ni.2067


```r
library(GEOquery)
# Download the mapping information and processed data
gds <- getGEO("GSE29619", destdir = "Data/GEO/")
```

```
## ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE29nnn/GSE29619/matrix/
## Found 3 file(s)
## GSE29619-GPL13158_series_matrix.txt.gz
## Using locally cached version: Data/GEO//GSE29619-GPL13158_series_matrix.txt.gz
## Using locally cached version of GPL13158 found here:
## Data/GEO//GPL13158.soft 
## GSE29619-GPL3921_series_matrix.txt.gz
## Using locally cached version: Data/GEO//GSE29619-GPL3921_series_matrix.txt.gz
## Using locally cached version of GPL3921 found here:
## Data/GEO//GPL3921.soft 
## GSE29619-GPL570_series_matrix.txt.gz
## Using locally cached version: Data/GEO//GSE29619-GPL570_series_matrix.txt.gz
## Using locally cached version of GPL570 found here:
## Data/GEO//GPL570.soft
```

```r
# main serie #gds[[1]] = LAIV/TIV 0809, gds[[2]] = FACS, gds[[3]] = TIV 0708
```

but before we can use this, we need to clean up the pData a bit (see code in .Rpres file).



## Using the sva package 

Let us estimate the (surrogate) factors as follows:

```r
library(sva)
```

```
## Loading required package: mgcv
## Loading required package: nlme
## 
## Attaching package: 'nlme'
## 
## The following object is masked from 'package:IRanges':
## 
##     collapse
## 
## This is mgcv 1.8-9. For overview type 'help("mgcv-package")'.
## Loading required package: genefilter
## 
## Attaching package: 'genefilter'
## 
## The following object is masked from 'package:base':
## 
##     anyNA
```

```r
TIV_08 <- gds_new[[1]][, grepl("2008-TIV", pData(gds_new[[1]])$description)]
mm_TIV_08 <- model.matrix(~ptid + time, TIV_08)
mm0_TIV_08 <- model.matrix(~ptid, TIV_08)
# Estimate the surrogate variables
sv_TIV_08 <- sva(exprs(TIV_08), mm_TIV_08, mod0 = mm0_TIV_08)
```

```
## Number of significant surrogate variables is:  12 
## Iteration (out of 5 ):1  2  3  4  5
```

## Using the limma with estimated SVs

Then we can use these variables in limma, as follows:

```r
library(limma)
# Add the surrogate variables to the design matrix
modSv <- cbind(mm_TIV_08, sv_TIV_08$sv)
# Use the new design matrix
fit_TIV_08 <- lmFit(TIV_08, modSv)
ebay_TIV_08 <- eBayes(fit_TIV_08)
topT7_sv <- topTable(ebay_TIV_08, coef = "timeD7", number = Inf)
# Compare to the old analysis
fit_TIV_08 <- lmFit(TIV_08, mm_TIV_08)
ebay_TIV_08 <- eBayes(fit_TIV_08)
topT7 <- topTable(ebay_TIV_08, coef = "timeD7", number = Inf)
```

## Using the limma with estimated SVs

The result from the adjusted analysis:


```r
topT7_sv[1:10, c("ID", "Gene.Symbol")]
```

```
##                            ID
## 211430_PM_s_at 211430_PM_s_at
## 232991_PM_at     232991_PM_at
## 219276_PM_x_at 219276_PM_x_at
## 239401_PM_at     239401_PM_at
## 239637_PM_at     239637_PM_at
## 227999_PM_at     227999_PM_at
## 226481_PM_at     226481_PM_at
## 216801_PM_at     216801_PM_at
## 241402_PM_at     241402_PM_at
## 235180_PM_at     235180_PM_at
##                                                          Gene.Symbol
## 211430_PM_s_at IGH@ /// IGHG1 /// IGHM /// IGHV4-31 /// LOC100290146
## 232991_PM_at                                                        
## 219276_PM_x_at                                               C9orf82
## 239401_PM_at                                                        
## 239637_PM_at                                                        
## 227999_PM_at                                                  PWWP2B
## 226481_PM_at                                                   VPRBP
## 216801_PM_at                                                        
## 241402_PM_at                                                  TSEN54
## 235180_PM_at                                                    STYX
```

## Using the limma with estimated SVs

The result from the un-adjusted analysis:


```r
topT7[1:10, c("ID", "Gene.Symbol")]
```

```
##                              ID
## 211430_PM_s_at   211430_PM_s_at
## 241824_PM_at       241824_PM_at
## 216576_PM_x_at   216576_PM_x_at
## 1559018_PM_at     1559018_PM_at
## 216207_PM_x_at   216207_PM_x_at
## 214669_PM_x_at   214669_PM_x_at
## 1568768_PM_s_at 1568768_PM_s_at
## 215946_PM_x_at   215946_PM_x_at
## 211645_PM_x_at   211645_PM_x_at
## 217157_PM_x_at   217157_PM_x_at
##                                                           Gene.Symbol
## 211430_PM_s_at  IGH@ /// IGHG1 /// IGHM /// IGHV4-31 /// LOC100290146
## 241824_PM_at                                                         
## 216576_PM_x_at              IGK@ /// IGKC /// LOC652493 /// LOC652694
## 1559018_PM_at                                                   PTPRE
## 216207_PM_x_at                                              IGKV1D-13
## 214669_PM_x_at            IGK@ /// IGKC /// IGKV3-20 /// LOC100291682
## 1568768_PM_s_at                                          LOC100302650
## 215946_PM_x_at                           IGLL1 /// IGLL3 /// LOC91316
## 211645_PM_x_at                                                       
## 217157_PM_x_at                            IGK@ /// IGKC /// LOC652493
```

What do you think?

# Pitfall 4: Replicability and reproducibility

## Life cycle of scientific discovery with omics technologies

<img src="./figures/experiment_overview.png" width="560">

## Pitfalls to avoid: summary

- Increased type I error rate and multiple testing

- Type II errors and statistical power
    - A bit more difficult for omics data due to multiplicity but many good software packages are available
    
- Confounding

- Experimental design

- Reproducibility

## Reproducibility (and replicability)

- Novel biomedical technologies generate large and high-dimensional datasets from individual experiments. Consequently, both experiments and analyses have become increasingly complex, rendering interpretation and replication of results difficult. 

- Highlighted the importance of defining precise study objectives, and implementing data management and analysis plans, as an integral part of experimental design

- Contribute significantly to the **reproducibility** and **replication** of an experiment or study

- Failing to follow the principles of good experimental design, researchers set themselves up for failure as such experiments often lead to “fishing expeditions” where data are tortured through a myriad of analysis pipelines and statistical tests until positive results (i.e. significant p-values) are found that are then reported in a manuscript. The tortuous data analysis procedures leading to those results are not described, resulting in an irreproducible and non-replicable study.

## Tools for reproducible research

In recent years, several open-source, community-based projects have emerged that enable researchers to construct and share complete and fully reproducible data analysis pipelines. Here are some examples:

- R/RStudio and Bioconductor for analysis in R
- GenomeSpace, GenePattern, Galaxy for building analysis pipelines using a GUI
- Git and GitHub for version control
- knitr and R markdown for authoring


