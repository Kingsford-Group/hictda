
# Overview

This repo contains the python source code to run the method described in "Topological data analysis reveals principles of chromosome structure in cellular differentiation" (Sauerwald, Shen and Kingsford), which applies persistent homology method to HiC data and anaylze structure changes across cell types. Specifically, the code will generate distance matrix from HiC contact matrix and use gudhi to output persistence pairs and all the simplices of the skeleton. It can be used for outputing persitence diagram, persistence barcode and other structural anaylsis. The code also contains the implementation of loop trace back algorithm for generating representative loops of persistence homology modules. 


# Dependencies

To use this python source code, you must have installed GUDHI library first. GUDHI is a python package used for topological data analysis(TDA), and you can visit the webpage (http://gudhi.gforge.inria.fr/python/latest/) for more information. 

Install GUDHI with conda: 
conda install -c conda-forge gudhi

To run loop track back algorithm, you should install NetworkX. NetworkX is a python package used for manipulating graphs and networks. You can visit the webpage (https://networkx.github.io) for more information. 

Install NetworkX with conda:
conda install -c anaconda networkx


# Usage
```
python HiC_TDA.py -i input_file_name -o output_file_name -p output_path -r resolution
```
This will output distance matrix of HiC data, persistence pairs and all the simplices appeared during the persistent homology, and also generate loops corresponding to each persistence homology module. 

Option Tag | Description
----------------------- | -----------------------------
-i \<inputfile>| the input file, a normalized HiC contact matrix
-p \<path> | the directory containing output files 
-o | the name of output files
-r | The resolution of HiC data

Here is an example:
```
python HiC_TDA.py -i example/RUES2_CM_combined_100000_iced_chr22.matrix -p example/output/ -o RUES2_CM_combined_100000_iced_chr22 -r 100000
```
Run this command line, and we will get four output files:

1: RUES2_CM_combined_100000_iced_chr22_distmat.txt
This file saves the distance matrix generated from original HiC contact matrix. 

2: RUES2_CM_combined_100000_iced_chr22_persisdiagram.txt
This file contains all the persistent diagram generated from persistent homology.
- the first column: the dimension of a homology class
- the second column: birth time
- the third column: death time
- the fourth column: persitence pairs

3: RUES2_CM_combined_100000_iced_chr22_skeleton.txt
This file contains all the simplices generated from persistent homology. 
- the first column: a set of nodes representing one simplex
- the second column: the birth time of the simplex

4: RUES2_CM_combined_100000_iced_chr22_loop_information.txt
This file contains loops of all the persistence homoloy modules generated from persistent homology. Every two rows represent the information of one loop. 
- the first row: the first column of this row is the unique edge of a persistence homology module, of which the vertices are chromosome loci. The second column is the number of edges which are different from the unique edge but have the same length. 
- the second row: nodes of the representative loop, they are represented by chromosome loci. 