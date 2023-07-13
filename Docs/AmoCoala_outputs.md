<img title="AmoCoale" alt="" src="/images/Logo_AmoCoala.png">

# Outputs files from AmoCoala

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
