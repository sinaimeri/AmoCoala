<img title="AmoCoale" alt="" src="/images/Logo_AmoCoala.png">

AmoCoala, is a tool for estimating the frequency of cophylogeny events based on an approximate Bayesian computation approach **in presence of multiple associations** between hosts and symbionts.  AmoCoala replaces the previous tool Coala and extends the classical 4 events
co-evolutionary model to include 2 additional **spreads events**: vertical and horizontal spreads that induce multiple associations. 


For more details about our methodology, please refer to our papers:

1. L. Urbini, B. Sinaimeri, M.-F. Sagot, and C. Matias. (2022). Cophylogeny Reconstruction Allowing for Multiple Associations Through Approximate Bayesian Computation. Submitted, [ArXiV preprint](https://arxiv.org/abs/2205.11084) 
 	- Note that a Supplementary Material for this paper is available in the [ArXiV preprint](https://arxiv.org/abs/2205.11084) version.
2. C. Baudet, B. Donati, B. Sinaimeri, P. Crescenzi, C. Gautier, C. Matias, and M.-F. Sagot. Cophylogeny reconstruction via an approximate Bayesian computation. *Systematic Biology*, 64(3): 416-431, 2015 [Journal link](https://doi.org/10.1093/sysbio/syu129).

<p align="center">
<img title="AmoCoale" alt="" src="docs/images/hr.png">
</p>

## Installation and requirements 

AmoCoala is released under the  [CeCILL](http://www.cecill.info/) license, compatible with [GNU GPL](http://www.gnu.org/licenses/gpl-3.0.html).

AmoCoala was developed in **Java** and tested in the Linux and Mac OS X environments. Additionally, in order to perform the clustering of the accepted vectors and produce graphical plots, AmoCoala requires the installation of R and of the **R packages dynamicTreeCut** and **ade4**.

Installation is through manual compilation and Jar file creation. You will need Java and the compiler javac.


1. Download the project folder from the provided source.

2. Open a command prompt or terminal on your system.

3. Navigate to the Project Directory (Use the `cd` command).
   
4. Execute the following command to compile the Java source code:

```  
javac -cp ".:./lib/commons-cli-1.2.jar:./lib/guava-14.0-rc3.jar:./lib/JSAT_r754.jar:./lib/junit-platform-console-standalone-1.10.0-M1.jar"  $(find ./* | grep .java) -Xlint
```

5. Once the code is successfully compiled, create the jar file with the following command:

```
jar cvfm AmoCoala.jar manifest.MF .
```

## Using AmoCoala

In a terminal, navigate to the directory where the JAR file is located and run the following command to execute the JAR file:
```
   java -jar AmoCoala.jar -input <nexusfile> [options]
```

**Examples**: The directory Datasets contains the 4 datasets (Nexus files) which we used in our study. You can process them with the command: 
```
    java -jar AmoCoala.jar -input Datasets/AP.nex -cluster -plot
```

<p align="center">
<img title="AmoCoale" alt="" src="images/hr.png">
</p>


## Input files 
AmoCoala receives a NEXUS file as input. The input file must contain a pair of trees (one host tree and one symbiont tree) and the (possibly multiple) association of their leaves. 
 
AmoCoala can read two types of Nexus files (.nex):

- Jane / TreeMap Nexus Format
- CoRe-Pa Nexus Format 

Please, verify that your file meets [the Nexus format description](nexus.md). 

<p align="center">
<img title="AmoCoale" alt="" src="images/hr.png">
</p>


## Options 
Here we list all options that are offered by the software:

	-a1 <value>	 
Defines the value of the constant $\alpha_1$ which is used to compute the metrics `LEAVES AND MAAC`, `EVENTS AND MAAC` and `EVENTS AND MAAC-MULTIPLE ASSOCIATIONS`.
Default value = 0.5 (for `LEAVES AND MAAC`) or 0.7 (for `EVENTS AND MAAC` and `EVENTS AND MAAC-MULTIPLE ASSOCIATIONS`).
If a value different from the default value is chosen, note that the 2 values from `-a1` and `-a2` options should sum to 1.

	-a2 <value>	
Defines the value of the constant $\alpha_2$ which is used to compute the metrics `LEAVES AND MAAC`, `EVENTS AND MAAC` and `EVENTS AND MAAC-MULTIPLE ASSOCIATIONS`.
Default value = 0.5 (for `LEAVES AND MAAC`) or 0.3 (for `EVENTS AND MAAC` and `EVENTS AND MAAC-MULTIPLE ASSOCIATIONS`).
If a value different from the default value is chosen, note that the 2 values from `-a1` and `-a2` options should sum to 1.


	-ac <value>	
Defines the value of $\alpha_{co-speciation}$ for the Dirichlet Process prior on the first round of the ABC-SMC process (sampling of the space of probability vectors). Default value = 1.0.

	-ad <value>	
Defines the value of $\alpha_{duplication}$ for the Dirichlet Process prior on the first round of the ABC-SMC process (sampling of the space of probability vectors). Default value = 1.0.

	-al <value>
Defines the value of $\alpha_{host-switch}$ for the Dirichlet Process prior on the first round of the ABC-SMC process (sampling of the space of probability vectors). Default value = 1.0.

	-as <value>	
Defines the value of $\alpha_{loss}$ for the Dirichlet Process prior on the first round of the ABC-SMC process (sampling of the space of probability vectors).
Default value = 1.0.

	-cluster	
If this option is present, the list of accepted vectors produced in the last round of the ABC-SMC process will be clustered using a hierarchical clustering procedure implemented by the R package `dynamicTreeCut`.

	-continue	
If this option is present, the software tries to continue a job that was started previously and interrupted in the middle. AmoCoala looks for the output files which are present in the output directory and starts the process from the last executed round. If no output file is found, AmoCoala starts the process from the beginning. Notice that, to use this option, the software must be configured with the same options of the job that was previously interrupted.

	-discard <value>
Defines the multiplier value that is used to compute the maximum number of vectors which can be discarded during the generation of the quantile populations (nb of rounds greater than 1). The maximum number of discarded vectors D is defined as: D = value × Q, where Q = N × threshold\_first\_round and N is the number of simulated vectors (option `-N`). Notice that D must be an integer greater than 1 (D > 1).
Default value = 10.

	-h
If this option is present, AmoCoala prints a help describing all available options and exits.

	-input <nexusfile>
Defines the path for the nexus file which contains the pair of host and parasite trees and their associations.

	-M <value>
Defines the number of trees which are going to be produced for each probability vector. Default value = 1000.

	-maxtree <value>
Defines the multiplier value that is used to compute the maximum number of trees which are going to be simulated in order to obtain the required number of trees M (trees that are more than 2 times bigger than the "real" symbiont tree are discarded). The maximum number of trees X is defined as: X = value × M. Notice that the value must be greater than 1 (value > 1).
Default value = 5.

	-metric <number>	
Defines the metric that will be used in the comparison between real and simulated trees. Currently, AmoCoala offers these 4 options: MAAC (1);     LEAVES AND MAAC (2);    EVENTS AND MAAC (3); EVENTS AND MAAC- MULTIPLE ASSOCIATIONS Metric (4).  Only the last one can handle trees with multiple associations.
Default metric = 4.

	-N <value>	
Defines the number of probability vectors which are going to be sampled in the first round of the ABC-SMC process.
Default value = 2000.

	-p <value>	
Defines a perturbation limit that is going to be applied to each element $p_i$  of a probability vector v in the refinement phases of the ABC-SMC process (rounds 2, 3, ...). During the perturbation routine, each probability $p_i$ receives an increment $delta_i$ that is uniformly sampled from the interval [-value,+value]. After that, the new vector v' is normalised such that the sum of all elements is equal to 1.
Default value = 0.01.
	
	-plot	
If this option is present, AmoCoala will produce plots with the results of each round.
	
	-R <value>
Defines the number of rounds of the ABC-SMC process.
Default value = 3. If a value different from the default value is chosen, the option `-t` must be specified.

	-root <value>	
Root mapping probability (real value in the interval [0.5,1.0]).
This probability value is used in the recursive process that chooses the starting position during the simulation of parasite trees. Starting from the root of the host tree, we generate a random number and compare it to the given probability value. If the random number is smaller than the probability value, the root of the host tree is chosen as starting point. Otherwise, we choose one of the two subtrees of the host root node to continue the recursion. Notice that to make a choice between the two subtrees, we attribute to each one a probability value which is proportional to the size of their leaf set.
Default value = 1.0

	-t <value>	
Defines a vector of tolerance values which are going to be used in each round of the ABC-SMC process. The list of tolerance values is composed by real numbers between 0 and 1, separated by commas (,). The size of this list must be equal to the number of rounds (option `-R`).
Default value = 0.10,0.25,0.25.

	-threads <value>	
Defines the number of threads that are going to be used to simulate symbiont  trees.
Default value = 1.


<p align="center">
<img title="AmoCoale" alt="" src="images/hr.png">
</p>


## Output files 
During its execution, AmoCoala produces some output files to register intermediate and final results. Given an input file `file.nex`, the program will produce the following files:

	file.nex.simul.round_X.csv

This file contains the list of probability vectors which were **simulated** during round X. It is a csv file which contains one line per probability vector and 7 (or 5) columns: probability of co-speciation, probability of duplication, probability of host-switch, probability of loss, number of vertical spreads (if metric=4), number of horizontal spreads (if metric=4) and observed distance.

	file.nex.accep.round_X.csv

This file contains the list of probability vectors which were **accepted** by the ABC rejection method at round X. It is a csv file which contains one line per probability vector and 7 (or 5) columns: probability of co-speciation, probability of duplication, probability of host-switch, probability of loss, number of vertical spreads (if metric=4), number of horizontal spreads (if metric=4) and distance observed.

	file.nex.plots.round_X.pdf

This PDF file contains a set of histograms which describe the results after round X. The plots show the distribution of the:

-  Probabilities of each event type, and number of vertical and horizontal spreads (if metric=4) among all **simulated** probability vectors (first row);

-  Probabilities of each event type, and number of vertical and horizontal spreads (if metric=4) among the **accepted** probability vectors (second row);

-  Distances observed among all simulated probability vectors and among the accepted probability vectors (third row).

This file is produced only if the option `-plot` is specified in the command line.
	
	file.nex.clust.round_X.Y.csv

This file contains a list of the probability vectors which were accepted during round X and were grouped together in the cluster Y. It is a csv file which has 6 columns: vector identifier, probability of co-speciation, probability of duplication, probability of host-switch, probability of loss, and observed distance. Additionally, the file contains statistics summaries (Min, Q1, Med, Mean, Q3, and Max) for each column and two proposals of representative probability vectors: one considering the average of each event probability (row `NMean`) and the other considering the median of each event probability (row `NMed`).

This file is produced only if the option `-cluster` is specified in the command line.

	file.nex.clusters.round_X.pdf

This file contains a plot which shows the projection of the clusters (of the list of accepted vectors during round X) on a plane.

This file is produced only if the options `-cluster` and `-plot` are specified in the command line.

**Observation**: All csv files use the TAB character (\t) as column separator.



<p align="center">
<img title="AmoCoale" alt="" src="images/hr.png">
</p>

## Contacts  
In case of doubts, comments, suggestions, or bug notifications, please write to bsinaimeri at luiss dot it. 

<p align="center">
<img title="AmoCoale" alt="" src="images/hr.png">
</p>

## License
All of the source code to this product is available under the [CeCILL](http://www.cecill.info/), compatible with [GNU GPL](http://www.gnu.org/licenses/gpl-3.0.html).
