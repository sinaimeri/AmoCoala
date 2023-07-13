<img title="AmoCoale" alt="" src="docs/images/Logo_AmoCoala.png">

This page describes the different NEXUS formats (CoRe-Pa and Jane / TreeMap) that you can use as input in AmoCoala. 

# CoRe-Pa Nexus Format 
[CoRe-Pa](https://doi.org/10.1186/1471-2105-11-S1-S60) has a particular type of nexus file, which is composed by many blocks and sub-blocks of information to configure many parameters of the application. AmoCoala is able to read this kind of file, but it only extracts information from the blocks listed in the example below:


	BEGIN TREES;
    TRANSLATE
     	H0	H0,
		H1	H1,
		H2	A,
		H3	B,
		H4	C,
		P0	P0,
		P1	P2,
		P2	C1,
		P3	B1,
		P4	P1,
		P5	A2,
		P6	A1
		;
    TREE HOST = ((H2,H3)H1,H4)H0;
    TREE PARASITE = ((P2,P3)P1,(P5,P6)P4)P0;
	ENDBLOCK;

	BEGIN COPHYLOGENY;
    PHI
	C1	C,
	B1	B,
	A2	A,
	A1	A
	;
	ENDBLOCK;
        

## Block TREES

The block TREES has information about the pair of trees. The sub-block TRANSLATE contains pairs \<label\_newick label\_real\> which define the translation from the labels which appear in the Newick representation to the real labels of the nodes of the trees. For instance, in the example above the real label of the host tree node H2 is A. The labels of each pair are separated by a TAB character (\t) and each pair is separated by a comma character (,).

After the sub-block TRANSLATE, the file contains two lines with the Newick representation of host and parasite trees. Notice that, in this file format, the Newick representation contains the labels of the internal nodes of the trees.

## Block COPHYLOGENY

In the CoRe-Pa Nexus file, the block COPHYLOGENY has many sub-blocks which configure the parameters of the application. In this block, AmoCoala extracts information only from the sub-block PHI. This sub-block is composed by pairs of the type \<leaf\_parasite leaf\_host\> where the elements are separated by a TAB character (\t) and each pair is separated by a comma character (,).

Notice that, the leaf labels that are listed in the pairs refer to the real labels and not to the labels that appear in the newick representation of the trees.
Important observations

The input file must contain only one pair of host and parasite (symbiont) trees.


References:

1. D. Merkle, M. Middendorf, and N. Wieseke. A parameter-adaptive dynamic programming approach for inferring cophylogenies. BMC Bioinformatics, 11(Supplementary 1):10 pages, 2010. [DOI](https://doi.org/10.1186/1471-2105-11-S1-S60)
    
# Jane / TreeMap Nexus Format

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


    