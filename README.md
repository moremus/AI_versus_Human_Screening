# AI versus Human Screening

<br>

## Human in the Loop Diagram
<img src="Images/Human_in_the_Loop.png" width="500">

<br>

## Overview

This repository contains the following files:

* CSV files from the screening process; and
* R or Latex programs.

File contents are summarized below.

<br>

## CSV Files

1. Catchii_Title_Abstract_Screening_Results.csv
2. Loon_Lens_Full_Text_Screening_Results.csv
3. Loon_Lens_Title_Abstract_Screening_Results.csv

These files contain the following information: article titles, authors, years of publication, journals, DOIs, whether the articles were included or excluded, reasons for inclusion or exclusion, and confidence levels for inclusion/exclusion decisions (Loon Lens only).

<br>

## R or Latex Programs

1. Class-imbalance-adjusted_Performance_Analysis: calculation of performance statistics obtained by chance alone and comparison with the point estimates in the manuscript; 
2. Confidence Concordance Regression: simple logistic regression of concordance on Loon Lens confidence scores;
3. Deduplication Script: use of R's ASySD library to remove duplicate citations obtained in the literature search;
4. Figure 3: generation of forest plots to show reliability statistics for title and abstract screening; 
5. Figure 4.1: confusion tables for performance statistics obtained by chance alone; and
6. Title and Abstract Screening Comparison: all analyses for the Loon Lens versus human title and abstract screening comparison, except for those analyses contained in the other R files.

*Please note:* for all screening analyses, we utilized the same code as shown in the Title_Abstract_Screening_Comparison.R file.

<br>

## Useful Links

[Catchii](https://catchii.org)
<br>
[Loon Inc.](https://loonbio.com)
