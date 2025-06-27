# Group-size-and-environmental-obstacles-drive-acoustic-call-properties-for-wild-gray-bats
Data and code for" Group size and environmental obstacles drive acoustic call properties for wild gray bats in flight: a data-driven analysis"

Data is found in the .csv files named "NF", "NM", "OF", and "OM" representing the environmental conditions: no obstacles few bats, no obstacles many bats, obstacles few bats, and obstacles many bats respectively. For each file, the columns are arranged by number of bats, calls per bat, acoustic power per bat, total calls, and total acoustic power. 

For this work, transfer entropy (TE) is calculated using the file "TE_wFigure.m", which uses the Java Information Dynamics Toolkit (JIDT) from Lizier found here: https://github.com/jlizier/jidt. Our code includes TE of permuted source arrays to compare the significance of the true TE via a z-score test. A figure of the TE of the original data vs a swarm plot of the TE of the permuted sources is also given by this code. 

Sparse identification of nonlinear dynamics (SINDy) is calculated using the file "SINDy_wFigure.m", which uses the MATLAB SINDy code base from Brunton found here: https://faculty.washington.edu/kutz/page26/. This file is divided into sections to make each step more organized. These sections are to be ran in order starting by loading data and adding a path to SINDy. This code includes spliting data and running test, train, and validate. The final section includes a figure to view your results as a color map. 
