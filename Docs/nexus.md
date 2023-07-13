<img title="AmoCoale" alt="" src="docs/images/Logo_AmoCoala.png">

This page describes the NEXUS formats (Jane / TreeMap) that you can use as input in AmoCoala. 


Jane [1] and TreeMap [2] use the same type of Nexus file (also known as tanglegram .tgl). This kind of file has the following structure:

	BEGIN HOST;
    	TREE HOST = ((A,B),C);
	ENDBLOCK;
	BEGIN PARASITE;
    	TREE PARASITE = ((A1,A2),(B1,C1));
	ENDBLOCK;
	BEGIN DISTRIBUTION;
    	RANGE
       C1 : C, B1 : B, A2 : A, A1 : A;
	ENDBLOCK;
        

As we can see in the example above, this Nexus file format is composed of three blocks of information.

## Block HOST

This block defines the host tree. This tree must be given in Newick format. Usually, in this kind of file the Newick format does not contain the labels of the internal nodes.

## Block PARASITE
This block defines the parasite (symbiont) tree. This tree must be given in Newick format. Usually, in this kind of file the Newick format does not contain the labels of the internal nodes.

## Block DISTRIBUTION

This block defines the association among host and parasite leaves. This list of assocations is composed by pairs \<Parasite\_Label : Host\_Label> separated by commas (,). In this file format, the separator is a colon character (:).


For technical reasons, AmoCoala creates label names for the internal nodes of the trees. These labels have the format #HX and #PX, respectively, for host and parasite trees, where X is an integer number. Thus, the user should avoid the use of labels like #H1, #H2, #P1, #P2, etc.


References:

   1. C. Conow, D. Fielder, Y. Ovadia, and R. Libeskind-Hadas. Jane: a new tool for the cophylogeny reconstruction problem. Algorithms for Molecular Biology, 5(16):10 pages, 2010.[DOI](https://doi.org/10.1186/1748-7188-5-16)

   2. R. D. M. Page. Parallel phylogenies: reconstructing the history of host-parasite assemblages. Cladistics, 10(2):155–173, 1994. [DOI](https://doi.org/10.1111/j.1096-0031.1994.tb00170.x)


    
