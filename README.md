# Fast Step Graph

![image](https://github.com/juancolonna/FastStepGraph/assets/6243522/2c0f3f9c-e675-4bd9-bd3e-0f4ff197efdb)

Fast Step Graph is an optimized implementation in **R** of the Stepwise approach used for discovering high-dimensional Gaussian Graphical Models. It aims to accurately estimate the $\mathbf{\Omega}$ precision matrix when dealing with datasets where the number of features is significantly larger than the number of samples ($p >> n$), such as in genomics.

This implementation builds upon the original code available at [link](https://jdssv.org/index.php/jdssv/article/view/11), which accompanies the associated research paper. **Fast Step Graph** enhances the computational efficiency of the original code to handle much larger graph structures than previously reported, while reducing the training time. Several improvements have been made, including the elimination of redundant code, utilization of column-wise data structures (better for R), avoid list creation, manipulation and expansion within loops, and the integration of a faster subroutine for regression. Additionally, this implementation addresses a bug introduced in the original code.

Despite these enhancements, the primary bottleneck of Fast Step Graph lies in the requirement of substantial memory resources $\text{Memory} \propto \Big(\frac{p(p-1)}{2}\Big)$ for storing the entire graph, particularly when $p$ grows.

Clone this repository or simply download the .zip file and follow the instructions in [this link to see an example](https://github.com/juancolonna/FastStepGraph/blob/main/vignettes/How_to_use.pdf).

How to cite this repository?

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8233706.svg)](https://doi.org/10.5281/zenodo.8233706)
