# Bacterial multi-omic oscillatory networks
### Spatial periodicty and velocity of propagation in bacterial multi-omic network dynamics
This repository contains the manuscript mentioned in the title, and associated source code, the whole data set of signals and some functions to obtain fast results with several multi-omic combinations within and between bacteria. Should you need help running our code, please contact us.





### Inter-organisms amplitude consensus (IOAC) procedure 
In order to compare the signals between different organisms, an amplitude discretization process was applied. The number of levels of discretization is decided through the following procedure, that it is called: Inter Organisms Amplitude Consensus (IOAC) [source code here](ioac_procedure/Data_norm.R). The relative dependencies to the RData datasets are explicitly indicated in the source code.

### Median change point detector in order to search the periodicity windows ( Algorithm 1) 

### Periodicity/Oscillation score and Velocity of propagation (Algorithm 2)
![image](figures/plot1_supp.png)

The optimized algorithm with which the scores of periodicity and propagation velocity of the multi-omic spatial signal were calculated is [here](SupplementaryAlgo2.R).
