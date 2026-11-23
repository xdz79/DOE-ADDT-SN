# Code for "Optimal designs for skew-normal accelerated destructive degradation test models"

This repository contains the R code used to compute the optimal designs and numerical results presented in the paper:

> **Optimal designs for skew-normal accelerated destructive degradation test models**
> Yun-Juan Wang, Xiao-Dong Zhou, Min-Jue Zhang, Rong-Xian Yue

## 📁 Repository Structure

The repository includes the following R script files:

| File Name | Description |
| :--- | :--- |
| `merge_function.R` | Contains helper functions (e.g., `merge_design_points`, `merge_2d_design_points`) used to merge nearby design points and format the output. |
| `Polymer_Design_d_fixed.R` | Computes the results for **Table 1** and **Table 3** in the paper. (Polymer material example, with the inspection interval $d$ fixed). |
| `Polymer_Design.R` | Computes the results for **Table 2** and **Table 4** in the paper. (Polymer material example, with both temperature and inspection interval as design variables). |
| `BondB_Design_d_fixed.R` | Computes the results for **Table 5** and **Table 7** in the paper. (Adhesive Bond B example, with the inspection interval $d$ fixed). |
| `BondB_Design.R` | Computes the results for **Table 6** and **Table 8** in the paper. (Adhesive Bond B example, with both temperature and inspection interval as design variables). |
| `Polymer_Design_robust_m.R` | Computes the results for **Table 9** in the paper. (Robust design criteria for the polymer material example). |


## 🛠 Prerequisites

To run these scripts, you need to have R installed. The code relies on the following R packages:

* `sn` (for skew-normal distribution)
* `ICAOD` (for finding optimal approximate designs using Imperialist Competitive Algorithm)
* `xtable` (for generating LaTeX tables)
* `plot3D` (for 3D contour plots, if applicable)

You can install the required packages using:

```R
install.packages(c("sn", "ICAOD", "xtable", "plot3D"))
```
## 🚀 How to Run

1.Clone or download this repository.

2.Open R or RStudio and set the working directory to the folder containing these scripts.

3.Important: Make sure merge_function.R is in the same directory, as the main scripts source this file (source("merge_function.R")).

4.Open the R script corresponding to the table you want to reproduce.

5.Adjust parameters $l$ and $d$: Since a single script often generates multiple tables with different settings, you need to adjust the number of observation times $l$ and the inspection interval distance $d$ manually in the code before running.

 For example, in Polymer_Design_d_fixed.R to reproduce different rows of Table 1 and Table 3, change the lines:
 
```R
l <- 5  # the number of observation times: 5, 10, 15
d <- 6  # the distance between the inspections: 6, 12, 24, 36
```
6.Run the script. The script will output the optimal design points, weights, and criterion values directly in the console, and will generate LaTeX code for the tables using xtable.

## 📝 Notes

1.The scripts use nominal values (MLEs) reported in the original paper for the parameters.

2.Due to the stochastic nature of the optimization algorithm (ICAOD), the exact numerical results may vary slightly depending on the random seed and iteration settings. The default settings in the scripts are chosen to reproduce the results in the paper as closely as possible.

## 📧 Contact
For any questions regarding the code or the paper, please contact the corresponding author, Xiao-Dong Zhou.
