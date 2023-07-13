<img title="AmoCoale" alt="" src="/images/Logo_AmoCoala.png">

AmoCoala, is a tool for estimating the frequency of cophylogeny events based on an approximate Bayesian computation approach **in presence of multiple associations** between hosts and symbionts.  AmoCoala replaces the previous tool Coala and extends the classical 4 events
co-evolutionary model to include 2 additional **spreads events**: vertical and horizontal spreads that induce multiple associations. 


For more details about our methodology, please refer to our papers:

1. L. Urbini, B. Sinaimeri, M.-F. Sagot, and C. Matias. (2022). Cophylogeny Reconstruction Allowing for Multiple Associations Through Approximate Bayesian Computation. Submitted, [ArXiV preprint](https://arxiv.org/abs/2205.11084) 
 	- Supplementary Material for this article is available [here](Docs/AmoCoala_SuppMat.pdf)
2. C. Baudet, B. Donati, B. Sinaimeri, P. Crescenzi, C. Gautier, C. Matias, and M.-F. Sagot. Cophylogeny reconstruction via an approximate Bayesian computation. *Systematic Biology*, 64(3): 416-431, 2015 [Journal link](https://doi.org/10.1093/sysbio/syu129).

<p align="center">
<img title="AmoCoale" alt="" src="images/hr.png">
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

6. Check the installation process by running this example:
```
     java -jar AmoCoala.jar -input Datasets/file.jane.nex  -metric 1 -N 10 -M 10 -R 1 -t 0.1
```


<p align="center">
<img title="AmoCoale" alt="" src="images/hr.png">
</p>

## Using AmoCoala

In a terminal, navigate to the directory where the JAR file is located and run the following command to execute the JAR file:
```
   java -jar AmoCoala.jar -input <nexusfile> [options]
```


<p align="center">
<img title="AmoCoale" alt="" src="images/hr.png">
</p>


## Input files 
AmoCoala receives a NEXUS file as input. The input file must contain a pair of trees (one host tree and one symbiont tree) and the (possibly multiple) association of their leaves. 
 
AmoCoala can read two types of Nexus files (.nex):

- Jane / TreeMap Nexus Format
- CoRe-Pa Nexus Format 

Please, verify that your file meets [the Nexus format description](Docs/nexus.md). 

<p align="center">
<img title="AmoCoale" alt="" src="images/hr.png">
</p>


## Options 
AmoCoala options are listed [here](Docs/AmoCoala_options.md). 


<p align="center">
<img title="AmoCoale" alt="" src="images/hr.png">
</p>


## Output files 
AmoCoala output files are described [here](Docs/AmoCoala_outputs.md). 



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
