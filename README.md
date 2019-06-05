# Bacterial multi-omic oscillatory networks
### Spatial periodicty and velocity of propagation in bacterial multi-omic network dynamics
This repository contains the manuscript mentioned in the title, and associated source code, the whole data set of signals and some functions to obtain fast results with several multi-omic combinations within and between bacteria. Should you need help running our code, please contact us.





#### Inter-organisms amplitude consensus (IOAC) procedure 
In order to compare the signals between different organisms, an amplitude discretization process was applied. The number of levels of discretization is decided through the following procedure, that it is called: Inter Organisms Amplitude Consensus (IOAC) [source code here](ioac_procedure/Data_norm.R). The relative dependencies to the RData datasets are explicitly indicated in the source code.

#### Median change point detector in order to search the periodicity windows ( Algorithm 1) 
The median change point detection algorithm for the estimation of the half periodicity search windows is calculated as described [here](SupplementaryAlgo1.R). 

#### Periodicity/Oscillation score and Velocity of propagation (Algorithm 2)
The optimized algorithm with which the scores of periodicity and propagation velocity of the multi-omic spatial signal were calculated is [here](SupplementaryAlgo2.R).


![image](figures/plot1_supp.png)

#### Whole dataset of the multi-omic signals
Note, in the paragraphs above we are providing a simple example, obviously these algorithms has been applied to all the signals generated for all the organisms and for all possible multi-omic combinations. To get the complete dataset with more than 2 million signals associated with all the functional classes, scores and various other statistics as it is described in the paper, it is possible to access the RData dataset [via this URL](https://thinfi.com/asxw). The password is: multi.org.2019.

In particular, this dataset presents these fields for each multi-omic signal:
1) "score" : the osc_s index, 
2) "m.s"   : the average value of osc_s
3) "med.s" : the median value of osc_s
4) "sd.s"  : the standard deviation of osc_s
5) "m.w"   : the average value of the period lengths
6) "med.w" : the median value of the period lengths
7) "sd.w"  : the standard deviation of the period lengths
8) "change.w"     : number of periodic oscillations
9) "v.change.w"   : velocity of propagation
10) "path.l"      : length of the signal
11) "n.path"      : KEGG pathway ID
12) "exp.cr"      : COLOMBOS condition contrast ID
13) "exp.ref"     : COLOMBOS treatment experiment ID
14) "exp.ctr"     : COLOMBOS control experiment ID
15) "kegg.id"     : KEGG organism ID
16) "code"        : ID of the multi-omic combination
17) "class"       : KEGG orthology level 1
18) "func"        : KEGG orthology level 2
19) "pathway_map" : KEGG orthology level 3

Note, for the field "code", these are the possible IDs:
1)  "n1" = CAI + Molecular Weight
2)  "n2" = CAI + mRNA CCs 
3)  "n3" = Molecular Weigth  + mRNA CCs 
4)  "n4" = CAI + Molecular Weigth  + mRNA CCs

While, for the operon compressed signals:
1)  "o1" = CAI + Molecular Weight
2)  "o2" = CAI + mRNA CCs 
3)  "o3" = Molecular Weigth  + mRNA CCs
4)  "o4" = CAI + Molecular Weigth  + mRNA CCs

#### Phase synchronizations and plots
