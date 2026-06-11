# AI versus Human Screening

<img src="Images/Appendix_D.png" width="500">

## Overview

This repository contains the following files:

* Appendix D PNG file;
* CSV files from the screening process; and
* R programs

File contents are summarized below.

<br>

## Appendix D

1. Human-in-the-loop workflow .png file that is Appendix D in the manuscript (shown above)

<br>

## CSV Files

1. Catchii_Title_Abstract_Screening_Results.csv
2. Loon_Lens_Full_Text_Screening_Results.csv
3. Loon_Lens_Title_Abstract_Screening_Results.csv

These files contain the following information: article titles, authors, years of publication, journals, DOIs, whether the articles were included or excluded, reasons for inclusion or exclusion, and confidence levels for inclusion/exclusion decisions (Loon Lens only)

<br>

## R Programs

1. Confidence Concordance Regression: simple logistic regression of concordance on Loon Lens confidence scores;
2. Deduplication Script: use of R's ASySD library to remove duplicate citations obtained in the literature search;
3. Figure 3: generation of forest plots to show reliability statistics for title and abstract screening; and
4. Title and Abstract Screening Comparison: all analyses for the Loon Lens versus human title and abstract screening comparison, except for those analyses contained in the other R files

*Please note:* for all screening analyses, we utilized the same code as shown in the Title_Abstract_Screening_Comparison.R file.

<br>

## Useful Links

[Catchii](https://catchii.org)
<br>
[Loon Inc.](https://loonbio.com)
