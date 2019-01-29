
# Overview

This repo contains the python source code to run the method described in "Topological data analysis reveals principles ofchromosome structure throughout cellulardifferentiation" (Sauerwald, Shen and Kingsford), which apply persistent homology method to HiC data and anaylze structure changes across cell types. Specifically, the code will generate distance matrix from HiC contact matrix and use gudhi to output persistence pairs and all the simplices of the skeleton. Besides, it can be used for outputing persitence diagram, persistence barcode and other structural anaylsis. 


# Dependencies

To use this python source code, you must have installed GUDHI library first. GUDHI is a python package used for topological data analysis(TDA), and you can visit the webpage (http://gudhi.gforge.inria.fr/python/latest/) for more information. 

Install GUDHI with conda: 
conda install -c conda-forge gudhi


# Usage
```
python HiC_TDA.py -i input_file_name -o output_file_name -p output_path -r resolution
```
This will generate distance matrix of HiC data and output persistence pairs and all the simplices appeared during the persistent homology. 

Option Tag | Description
----------------------- | -----------------------------
-i \<inputfile>| the input file, a normalized HiC contact matrix
-p \<path> | the directory containing output files 
-o | the name of output files
-r | The resolution of HiC data

Here is an example:
```
python HiC_TDA.py -i example/input/RUES2_CM_combined_100000_iced_chr22.matrix -p example/output/ -o RUES2_CM_combined_100000_iced_chr22 -r 100000
```
Run this command line, and we will get three output files:

1: RUES2_CM_combined_100000_iced_chr22_distmat.txt
This file saves the distance matrix generated from original HiC contact matrix. 

2: RUES2_CM_combined_100000_iced_chr22_persisdiagram.txt
This file contains all the persistent diagram generated from persistent homology
- the first column: the dimension of a homology class
- the second column: birth time
- the third column: death time
- the fourth column: persitence pairs

3: RUES2_CM_combined_100000_iced_chr22_skeleton.txt
This file contains all the sinmplices generated from persistent homology. 
- the first column: a set of nodes representing one simplex
- the second column: the birth time of the simplex