# Bacterial multi-omic oscillatory networks
### Spatial periodicty and velocity of propagation in bacterial multi-omic network dynamics
This repository contains the manuscript mentioned in the title, and associated source code, the whole data set of signals and some functions to obtain fast results with several multi-omic combinations within and between bacteria. Should you need help running our code, please contact us.





### Inter-organisms amplitude consensus (IOAC) procedure 
In order to compare the signals between different organisms, an amplitude discretization process was applied. The number of levels of discretization is decided through the following procedure, that it is called: Inter Organisms Amplitude Consensus (IOAC) [source code here](ioac_procedure/Data_norm.R). The relative dependencies to the RData datasets are explicitly indicated in the source code.

### Median change point detector in order to search the periodicity windows ( Algorithm 1) 
The median change point detection algorithm for the estimation of the half periodicity search windows is calculated as described [here](SupplementaryAlgo1.R). 

### Periodicity/Oscillation score and Velocity of propagation (Algorithm 2)
The optimized algorithm with which the scores of periodicity and propagation velocity of the multi-omic spatial signal were calculated is [here](SupplementaryAlgo2.R).


![image](figures/plot1_supp.png)

Note, in the paragraphs above we are providing a simple example, obviously these algorithms has been applied to all the signals generated for all the organisms and for all possible multi-omic combinations. To get the complete dataset with more than 2 million signals associated with all the functional classes, scores and various other statistics as it is described in the paper, it is possible to access the RData dataset [via this URL](https://thinfi.com/asxw). The password is: multi.org.2019 
