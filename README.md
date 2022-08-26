# ScRNA Analysis and Trajectory Inference
This workflow performs quality control, clustering, integration, and trajectory inference on the STARSolo outputs of 2 prostate cancer scRNA datasets. 

## Table of Contents 
* Installation Requirements 
* Set Up + Dataset 
* Author
* References 

## Installation Requirements
### Software 
* R 
* R studio

### Packages 
* Seurat 
* tidyverse
* stringr
* scales
* R.utils 
* SingleCellExperiment 
* slingshot 
* grDevices 
* GiNA

## Dataset + Set up
The dataset used in this workflow was STARSolo outputs that can be accessed through this [dropbox link](https://www.dropbox.com/scl/fo/0cvy5v1e30srupklrswni/h?dl=0&rlkey=7hzcn0846qdj5avzstp5zgo0j).

## Author 
Amara Alexander

## References 
Mary Piper, Meeta Mistry, Jihe Liu, William Gammerdinger, & Radhika Khetani. (2022, January 6). hbctraining/scRNA-seq_online: scRNA-seq Lessons from HCBC (first release). Zenodo. https://doi.org/10.5281/zenodo.5826256

Riemondy, K. (2019). Single-cell RNA-seq Workshop: Trajectory inference. Retrieved 26 August 2022, from https://rnabioco.github.io/cellar/previous/2019/docs/5_trajectories.html#:~:text=Slingshot%20is%20a%20Bioconductor%20package%20that%20draws%20curved,provides%20functionality%20for%20computing%20pseudotimes%20through%20multiple%20trajectories.
