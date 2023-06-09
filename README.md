# Fast Step Graph

Fast Step Graph is an optimized implementation in **R** of the Stepwise approach used for discovering high-dimensional Gaussian graphical models. It aims to accurately estimate the $\mathbf{\Omega}$ precision matrix when dealing with datasets where the number of features is significantly larger than the number of samples ($p >> n$), such as in genomics.

This implementation builds upon the original code available at [link](https://jdssv.org/index.php/jdssv/article/view/11), which accompanies the associated research paper. **Fast Step Graph** enhances the computational efficiency of the original code to handle much larger graph structures than previously reported, while reducing the training time. Several improvements have been made, including the elimination of redundant code, utilization of column-wise data structures (better for R), avoid list creation, manipulation and expansion within loops, and the integration of a faster subroutine for regression. Additionally, this implementation addresses a bug introduced in the original code.

Despite these enhancements, the primary bottleneck of Fast Step Graph lies in the requirement of substantial memory resources $memory \propto \Big(\frac{p(p-1)}{2}\Big)$ for storing the entire graph, particularly when $p$ grows.


Clone this repository or simply download the .zip and follow the instructions in [this link to run an example](instructions.md).
