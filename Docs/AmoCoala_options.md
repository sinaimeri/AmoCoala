<img title="AmoCoale" alt="" src="/images/Logo_AmoCoala.png">

# Options of the software AmoCoala

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
