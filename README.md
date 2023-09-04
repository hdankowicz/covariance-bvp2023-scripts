# covariance-bvp2023-scripts

This repository contains Matlab scripts for generating the figures in the manuscript

Zaid Ahsan, Harry Dankowicz, Christian Kuehn  (2023) *Adjoint-Based Projections for Uncertainty Quantification near Stochastically Perturbed Limit Cycles and Tori*

arXiv: 


## Content
This repository contains:

- `coco_2023August15`: August 15, 2023, release of COCO 
  [sourceforge.net/cocotools](https://sourceforge.net/projects/cocotools/)
  used for demos (note a change on line 35 of coco_save_full, adding '-v7.3'
  to the sequence of arguments)
- `Figure1.m`: script for generating Figure 1 in Section 2.2 of the paper
- `Figure2.m`: script for generating Figure 2 in Section 2.2 of the paper
- `Figure3.m`: script for generating Figure 3 in Section 2.2 of the paper
- `Figure3_validation.m`: script for solving the covariance boundary-value
  problem corresponding to Figure 3 in Section 2.2 of the paper
- `Figure4.m`: script for generating Figure 4 in Section 3.2 of the paper
- `Figure4_validation_nonautonomous.m`: script for solving the covariance
  boundary-value problem corresponding to Figure 4 in Section 3.2 of the
  paper using a non-autonomous encoding
- `Figure4_validation_autonomous.m`: script for solving the covariance
  boundary-value problem corresponding to Figure 4 in Section 3.2 of the
  paper using an autonomous encoding
- `Figure5.m`: script for generating Figure 5 in Section 3.2 of the paper
- `Figure6.m`: script for generating Figure 6 in Section 3.2 of the paper
- `Figures7and8.m`: script for generating Figures 7 and 8 in Section 3.2
  of the paper
- `projection.m`: zero problem for identifying coordinates of torus
  projection used to construct Figures 7 and 8
- `interpolant.m`: zero problem for straight-line interpolant used to
  construct Figures 7 and 8 
- `findPoints.m`: compute point on torus given its parameterization used to
  construct Figures 7 and 8 
- `findCovars.m`: compute covariance matrix for point on torus given its
  parameterization used to construct Figures 6, 7, and 8 


## Usage

1. execute script `startup.m` in folder `coco_2023August15/coco` to initialize `coco`
2. execute each of the scripts `Figure*.m` to generate the corresponding figure(s)
