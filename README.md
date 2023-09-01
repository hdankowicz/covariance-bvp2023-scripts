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
- `Section 2.2`: scripts for generating figures in Section 2.2 of the paper,
  as well as additional demos validating the encoding
- `Section 3.2`: scripts for generating figures in Section 3.2 of the paper,
  as well as additional utility functions used for data processing


## Usage

1. execute script `startup.m` in folder `coco_2023August15/coco` to initialize `coco`
2. execute each of the scripts `Figure*.m` to generate the corresponding figure(s)
