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

*(Note: The R code files have been renamed slightly for standard GitHub compatibility, removing the `(g)` label from the filenames.)*

## 🛠 Prerequisites

To run these scripts, you need to have R installed. The code relies on the following R packages:

* `sn` (for skew-normal distribution)
* `ICAOD` (for finding optimal approximate designs using Imperialist Competitive Algorithm)
* `xtable` (for generating LaTeX tables)
* `plot3D` (for 3D contour plots, if applicable)

You can install the required packages using:

```R
install.packages(c("sn", "ICAOD", "xtable", "plot3D"))
