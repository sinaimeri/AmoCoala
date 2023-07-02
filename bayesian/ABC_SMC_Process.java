/* 
############################################################################
#                                Coala                                     #
#         CO-evolution Assessment by a Likelihood-free Approach            #
############################################################################
#                                                                          #
# Copyright INRIA 2018                                                    #
#                                                                          #
# Contributors : Christian Baudet                                          #
#                Beatrice Donati                                           #
#                Blerina Sinaimeri                                         #
#                Pierluigi Crescenzi                                       #
#                Christian Gautier                                         #
#                Catherine Matias                                          #
#                Marie-France Sagot                                        #
#				 Laura Urbini											   #
#                                                                          #
# christian.baudet@inria.fr                                                #
# marie-france.sagot@inria.fr                                              #
# https://gforge.inria.fr/forum/forum.php?forum_id=11324                   #
#                                                                          #
# This software is a computer program whose purpose is to estimate         #
# the probability of events (co-speciation, duplication, host-switch and   #
# loss) for a pair of host and parasite trees through a Approximate        #
# Bayesian Computation - Sequential Monte Carlo (ABC-SMC) procedure.       #
#                                                                          #
# This software is governed by the CeCILL  license under French law and    #
# abiding by the rules of distribution of free software.  You can  use,    # 
# modify and/ or redistribute the software under the terms of the CeCILL   #
# license as circulated by CEA, CNRS and INRIA at the following URL        #
# "http://www.cecill.info".                                                # 
#                                                                          #
# As a counterpart to the access to the source code and  rights to copy,   #
# modify and redistribute granted by the license, users are provided only  #
# with a limited warranty  and the software's author,  the holder of the   #
# economic rights,  and the successive licensors  have only  limited       #
# liability.                                                               #
#                                                                          #
# In this respect, the user's attention is drawn to the risks associated   #
# with loading,  using,  modifying and/or developing or reproducing the    #
# software by the user in light of its specific status of free software,   #
# that may mean  that it is complicated to manipulate,  and  that  also    #
# therefore means  that it is reserved for developers  and  experienced    #
# professionals having in-depth computer knowledge. Users are therefore    #
# encouraged to load and test the software's suitability as regards their  #
# requirements in conditions enabling the security of their systems and/or #
# data to be ensured and,  more generally, to use and operate it in the    # 
# same conditions as regards security.                                     #
#                                                                          #
# The fact that you are presently reading this means that you have had     #
# knowledge of the CeCILL license and that you accept its terms.           #
############################################################################
 */

package  bayesian;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.NumberFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Random;

import  cycle.CyclicityTest;
import  threads.GenerateDistribution;
import  threads.PriorDistribution;
import  threads.QuantileDistribution;
import  threads.SynchronizedVectorManager;
import  trees.Tree;
import  trees.TreeNode;
import  util.NexusFileParser;
import  util.NexusFileParserException;
import  util.Statistics;

public class ABC_SMC_Process {

	/*
	 * -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:
	 * -:-:-:-:-:-:-:-
	 */
	/*
	 * Constants which define the default values for the input parameters of the
	 * process.
	 */
	public static final double DEFAULT_ALPHA_COSPECIATION = 1.0;

	public static final double DEFAULT_ALPHA_DUPLICATION = 1.0;

	public static final double DEFAULT_ALPHA_LOSS = 1.0;

	public static final double DEFAULT_ALPHA_SWITCH = 1.0;

	public static final double DEFAULT_LEAVES_AND_MAAC_ALPHA1 = GenerateDistribution.DEFAULT_LEAVES_AND_MAAC_ALPHA1;

	public static final double DEFAULT_LEAVES_AND_MAAC_ALPHA2 = GenerateDistribution.DEFAULT_LEAVES_AND_MAAC_ALPHA2;

	public static final double DEFAULT_EVENTS_AND_MAAC_ALPHA1 = GenerateDistribution.DEFAULT_EVENTS_AND_MAAC_ALPHA1;

	public static final double DEFAULT_EVENTS_AND_MAAC_ALPHA2 = GenerateDistribution.DEFAULT_EVENTS_AND_MAAC_ALPHA2;

	public static final int DEFAULT_MAXIMUM_NUMBER_OF_TREES_FACTOR = GenerateDistribution.DEFAULT_MAXIMUM_NUMBER_OF_TREES_FACTOR;

	public static final boolean DEFAULT_CLUSTERING = false;

	public static final boolean DEFAULT_CONTINUE = false;

	public static final int DEFAULT_METRIC = GenerateDistribution.METRIC_EVENTS_AND_MAAC_DISTANCE;

	public static final int DEFAULT_NUMBER_OF_ITERATIONS = 3;

	public static final int DEFAULT_MODEL = GenerateDistribution.DEFAULT_MODEL;

	public static final int DEFAULT_NUMBER_OF_TREES = 1000;

	public static final int DEFAULT_PRIOR_DISTRIBUTION_SIZE = 2000;

	public static final boolean DEFAULT_PLOT = false;

	public static final double[] DEFAULT_TOLERANCES = { 0.10, 0.25, 0.25 };

	public static final double DEFAULT_PERTURBATION = 0.01;

	public static final double DEFAULT_ROOT_MAPPING_PROBABILITY = 1.0;

	public static final int DEFAULT_NUMBER_OF_THREADS = 1;

	public static final int DEFAULT_CYCLICITY_TEST = CyclicityTest.DONATI;

	public static final int DEFAULT_MAXIMUM_NUMBER_OF_DISCARDED_VECTORS_FACTOR = 10;

	public static final int NCOLUMNS = 7;

	public static final int COL_COSPECIATION = 0;
	
	public static final int COL_DUPLICATION = 1;
	
	public static final int COL_HOSTSWITCH = 2;
	
	public static final int COL_LOSS = 3;
	
	public static final int COL_VERTICALSPREAD = 4;
	
	public static final int COL_HORIZONTALSPREAD = 5;
	
	public static final int COL_DISTANCE = 6;
	
	//public static final int COL_MAXJUMP = 5;
	
	//public static final int COL_MEANJUMP = 6;
	
	//public static final int COL_MEAN_MAXJUMP = 7;
	
	//public static final int COL_MEAN_MEANJUMP = 8;
	
	//public static final int COL_NB_MAXJUMP = 9;
	
	//public static final int COL_NB_MEANJUMP = 10;
	
	

	/*
	 * -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:
	 * -:-:-:-:-:-:-:-
	 */
	/* Input parameters of the process. */
	private double alphaCospeciation;

	private double alphaDuplication;

	private double alphaLoss;

	private double alphaSwitch;

	private double alpha1;

	private double alpha2;

	private int maximumNumberOfTreesFactor;

	private boolean clustering;

	private boolean continueProcess;

	private int metric;

	private int numberOfIterations;

	private String inputFileName;

	private int model;

	private int numberOfTrees;

	private int priorDistributionSize;

	private boolean plot;

	private double[] tolerances;

	private double perturbation;

	private int numberOfThreads;

	private int cyclicityTest;

	private int maximumNumberOfDiscardedVectorsFactor;

	private double rootMappingProbability;

	/*
	 * -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:
	 * -:-:-:-:-:-:-:-
	 */
	/* Work attributes */
	private Tree hostTree;

	private Tree parasiteTree;

	//private HashMap<String, String> mappingParasiteHost;
	private HashMap<String, List<String>> mappingParasiteHost;

	private String configurationDescription;

	public Random random;

	private static final SimpleDateFormat dateFormat = new SimpleDateFormat(
			"dd/MM/yyyy HH:mm:ss");

	private Calendar calendar;

	private NumberFormat format;

	/**
	 * Constructor: Create an ABC_SMC object filled with default parameters.
	 */
	public ABC_SMC_Process() {
		this.alphaCospeciation = DEFAULT_ALPHA_COSPECIATION;
		this.alphaDuplication = DEFAULT_ALPHA_DUPLICATION;
		this.alphaSwitch = DEFAULT_ALPHA_SWITCH;
		this.alphaLoss = DEFAULT_ALPHA_LOSS;
		this.alpha1 = DEFAULT_EVENTS_AND_MAAC_ALPHA1;
		this.alpha2 = DEFAULT_EVENTS_AND_MAAC_ALPHA2;
		this.maximumNumberOfTreesFactor = DEFAULT_MAXIMUM_NUMBER_OF_TREES_FACTOR;
		this.clustering = DEFAULT_CLUSTERING;
		this.continueProcess = DEFAULT_CONTINUE;
		this.metric = DEFAULT_METRIC;
		this.numberOfIterations = DEFAULT_NUMBER_OF_ITERATIONS;
		this.model = DEFAULT_MODEL;
		this.numberOfTrees = DEFAULT_NUMBER_OF_TREES;
		this.priorDistributionSize = DEFAULT_PRIOR_DISTRIBUTION_SIZE;
		this.plot = DEFAULT_PLOT;
		this.tolerances = DEFAULT_TOLERANCES;
		this.perturbation = DEFAULT_PERTURBATION;
		this.numberOfThreads = DEFAULT_NUMBER_OF_THREADS;
		this.cyclicityTest = DEFAULT_CYCLICITY_TEST;
		this.maximumNumberOfDiscardedVectorsFactor = DEFAULT_MAXIMUM_NUMBER_OF_DISCARDED_VECTORS_FACTOR;
		this.rootMappingProbability = DEFAULT_ROOT_MAPPING_PROBABILITY;
		this.random = new Random();
		this.calendar = Calendar.getInstance();
		this.format = NumberFormat.getInstance(Locale.US);
		this.format.setMinimumFractionDigits(4);
		this.format.setMaximumFractionDigits(4);
	}

	/*
	 * -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:
	 * -:-:-:-:-:-:-:-
	 */
	/* Set methods for process configuration. */

	public void setAlphaCospeciation(double alpha) {
		this.alphaCospeciation = alpha;
	}

	public void setAlphaDuplication(double alpha) {
		this.alphaDuplication = alpha;
	}

	public void setAlphaLoss(double alpha) {
		this.alphaLoss = alpha;
	}

	public void setAlphaSwitch(double alpha) {
		this.alphaSwitch = alpha;
	}

	public void setAlpha1(double alpha) {
		this.alpha1 = alpha;
	}

	public void setAlpha2(double alpha) {
		this.alpha2 = alpha;
	}

	public void setMaximumNumberOfTreesFactor(int value) {
		this.maximumNumberOfTreesFactor = value;
	}

	public void setClustering(boolean clustering) {
		this.clustering = clustering;
	}

	public void setContinueProcess(boolean continueProcess) {
		this.continueProcess = continueProcess;
	}

	public void setMetric(int metric) {
		this.metric = metric;
	}

	public void setNumberOfIterations(int iterations) {
		this.numberOfIterations = iterations;
	}

	public void setInputFile(String file) {
		this.inputFileName = file;
	}

	public void setModel(int model) {
		this.model = model;
	}

	public void setNumberOfTrees(int trees) {
		this.numberOfTrees = trees;
	}

	public void setMaximumNumberOfDiscardedVectorsFactor(int factor) {
		this.maximumNumberOfDiscardedVectorsFactor = factor;
	}

	public void setPriorDistributionSize(int size) {
		this.priorDistributionSize = size;
	}

	public void setPlot(boolean plot) {
		this.plot = plot;
	}

	public void setTolerances(double[] tolerances) {
		this.tolerances = tolerances;
	}

	public void setPerturbation(double perturbation) {
		this.perturbation = perturbation;
	}

	public void setRootMappingProbability(double probability) {
		this.rootMappingProbability = probability;
	}

	public void setNumberOfThreads(int threads) {
		this.numberOfThreads = threads;
	}

	public void setCyclicityTest(int cyclicityTest) {
		this.cyclicityTest = cyclicityTest;
	}

	/*
	 * -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:
	 * -:-:-:-:-:-:-:-
	 */
	/**
	 * Execute the ABC SMC process.
	 * 
	 * @throws IOException
	 *             When something goes wrong while reading the file.
	 * @throws NexusFileParserException
	 *             When something goes wrong while parsing the file.
	 */
	public void run() throws IOException, NexusFileParserException {
		readNexusFile(inputFileName);
		printConfiguration();
		runABC_SMC();
	}

	/**
	 * Read the nexus file and set up some parameters.
	 * 
	 * @param fileName
	 *            Nexus file name.
	 * @throws IOException
	 *             When something goes wrong while reading the file.
	 * @throws NexusFileParserException
	 *             When something goes wrong while parsing the file.
	 */
	private void readNexusFile(String fileName) throws IOException,
			NexusFileParserException {

		NexusFileParser parser = new NexusFileParser(fileName);

		hostTree = parser.getHostWithRHO();
		parasiteTree = parser.getParasite();

		/*
		 * Correct the leaf label of the parasite tree. In this process, we make
		 * assign to the parasite leaf the same label of the host leaf where it
		 * is mapped into.
		 */
		mappingParasiteHost = parser.getPhiNewickNames();
		/*for (TreeNode leaf : parasiteTree.getLeafNodes()) {
			leaf.setLabel(mappingParasiteHost.get(leaf.getLabel()));
		}*/
		
		
		for (TreeNode leaf : parasiteTree.getLeafNodes()) {
			leaf.setLabel(mappingParasiteHost.get(leaf.getLabel().get(0)));
		}

		
		/* Create the LCA tables for host and parasite trees. */
		parasiteTree.updateNumberOfMultipleAssociations();
		hostTree.createLCATable();
		parasiteTree.createLCATable();

		
		//System.out.println("hostTree " + hostTree.toString() );
		//System.out.println("parasiteTree " + parasiteTree.toString() );
		//System.out.println("mappingParasiteHost " + mappingParasiteHost.toString());
		
		
		/*Create the FREEZE table*/
		hostTree.createVERTICALSPREADAssociations(mappingParasiteHost, parasiteTree);
		//System.out.print("MapVerticalSpreadAssociations : ");
		//hostTree.printMapVerticalSpreadAssociations();
		//System.out.print("printNbParasitesSubTree : ");
		//hostTree.printNbParasitesSubTree();
		hostTree.createVERTICALSPREADTable(parasiteTree);
		//System.out.print("printVerticalSpreadMap : ");
		//hostTree.printVerticalSpreadMap();
		
		
		//System.out.println("---------------------");
		/*Create the HORIZONTAL SPREAD table*/
		hostTree.createHORIZONTALSPREADAssociations(mappingParasiteHost, parasiteTree);
		//System.out.println("printHorizontalSpreadMatrix : ");
		//hostTree.printHorizontalSpreadMatrix();
		
		/* Update the list of ancestors. */
		hostTree.updateAncestorList();
		parasiteTree.updateAncestorList();

		/* Update the list of descendants. */
		hostTree.updateDescendantList();
		parasiteTree.updateDescendantList();
		

		if (model != GenerateDistribution.DEFAULT_MODEL) {
			hostTree.getMapLabelToNode();
		}

	}

	/**
	 * Print the program configuration.
	 */
	private void printConfiguration() {

		int maximumNumberOfTrees = numberOfTrees * maximumNumberOfTreesFactor;

		String[] elements = inputFileName.split("/");
		String file = elements[elements.length - 1];

		String aux[] = new String[tolerances.length];
		for (int i = 0; i < tolerances.length; i++) {
			aux[i] = format.format(tolerances[i]);
		}
		String tol = Arrays.toString(aux).replaceAll("[\\[\\] ]", "");

		configurationDescription = "Nexus file: " + file + "\\n";
		configurationDescription += "Host/Parasite tree: "
				+ hostTree.getLeafNodes().length + "/"
				+ parasiteTree.getLeafNodes().length + " leaves\\n";
		configurationDescription += "Prior distribution: "
				+ priorDistributionSize + " vectors\\n";
		configurationDescription += "Number of trees: " + numberOfTrees
				+ " trees\\n";
		configurationDescription += "Maximum number of trees: "
				+ maximumNumberOfTrees + " trees\\n";
		configurationDescription += "Number of rounds: " + numberOfIterations
				+ " rounds\\n";
		configurationDescription += "Perturbation: "
				+ format.format(perturbation) + "\\n";
		configurationDescription += "Tolerances: " + tol + "\\n";

		System.out
				.println("---------------------------------------------------------------------");
		System.out.println("CONFIGURATION:");
		System.out.println("Nexus file              : " + file);
		System.out.println("Host tree               : "
				+ hostTree.getLeafNodes().length + " leaves");
		System.out.println("Parasite tree           : "
				+ parasiteTree.getLeafNodes().length + " leaves");
		System.out.println("Prior distribution      : " + priorDistributionSize
				+ " vectors");
		System.out.println("Number of trees         : " + numberOfTrees
				+ " trees");
		System.out.println("Maximum number of trees : " + maximumNumberOfTrees
				+ " trees");
		System.out.println("Number of rounds        : " + numberOfIterations
				+ " rounds");
		System.out.println("Perturbation            : "
				+ format.format(perturbation));
		System.out.println("Tolerances              : " + tol);
		System.out.println("Maximum discard factor  : "
				+ maximumNumberOfDiscardedVectorsFactor);

		if (model == GenerateDistribution.DEFAULT_MODEL) {
			System.out.println("Simulation Model        : " + model
					+ " - From the root to the leaves model");
			configurationDescription += "Simulation Model: " + model
					+ " - From the root to the leaves model\\n";
		} else if (model == GenerateDistribution.COALESCENCE_WITHOUT_DISTANCE) {
			System.out.println("Simulation Model        : " + model
					+ " - Coalescent model (Ignore host node distances).");
			configurationDescription += "Simulation Model: " + model
					+ " - Coalescent model (Ignore host node distances).\\n";
		} else {
			System.out.println("Simulation Model        : " + model
					+ " - Coalescent model (Consider host node distances).");
			configurationDescription += "Simulation Model: " + model
					+ " - Coalescent model (Consider host node distances).\\n";
		}

		if (cyclicityTest == CyclicityTest.STOLZER) {
			System.out.println("Cyclicity test          : "
					+ CyclicityTest.STOLZER + " - Stolzer et al., 2012.");
			configurationDescription += "Cyclicity test: "
					+ CyclicityTest.STOLZER + " - Stolzer et al., 2012.\\n";

		} else if (cyclicityTest == CyclicityTest.DONATI) {
			System.out.println("Cyclicity test          : "
					+ CyclicityTest.DONATI + " - Donati et al., 2014.");
			configurationDescription += "Cyclicity test: "
					+ CyclicityTest.DONATI + " - Donati et al., 2014.\\n";
		} else if (cyclicityTest == CyclicityTest.TOFIGH) {
			System.out.println("Cyclicity test          : "
					+ CyclicityTest.TOFIGH + " - Tofigh et al., 2011.");
			configurationDescription += "Cyclicity test: "
					+ CyclicityTest.TOFIGH + " - Tofigh et al., 2011.\\n";
		} else {
			System.out.println("Cyclicity test          : "
					+ CyclicityTest.TRANSFER_EDGES + " - Only transfer edges.");
			configurationDescription += "Cyclicity test: "
					+ CyclicityTest.TRANSFER_EDGES
					+ " - Only transfer edges\\n";
		}

		switch (metric) {
		case GenerateDistribution.METRIC_MAAC_DISTANCE:
			System.out.println("Metric                  : "
					+ GenerateDistribution.METRIC_MAAC_DISTANCE
					+ " - MAAC Metric");
			configurationDescription += "Metric: "
					+ GenerateDistribution.METRIC_MAAC_DISTANCE
					+ " - MAAC Distance\\n";
			break;
		case GenerateDistribution.METRIC_LEAVES_AND_MAAC_DISTANCE:
			System.out.println("Metric                  : "
					+ GenerateDistribution.METRIC_LEAVES_AND_MAAC_DISTANCE
					+ " - LEAVES AND MAAC Metric");
			System.out.println("(alpha1/alpha2): " + format.format(alpha1)
					+ "/" + format.format(alpha2));
			configurationDescription += "Metric: "
					+ GenerateDistribution.METRIC_LEAVES_AND_MAAC_DISTANCE
					+ " - LEAVES AND MAAC Metric\\n";
			configurationDescription += "(alpha1/alpha2): "
					+ format.format(alpha1) + "/" + format.format(alpha2)
					+ "\\n";
			break;
		case GenerateDistribution.METRIC_EVENTS_AND_MAAC_DISTANCE:
			System.out.println("Metric                  : "
					+ GenerateDistribution.METRIC_EVENTS_AND_MAAC_DISTANCE
					+ " - EVENTS AND MAAC Metric");
			System.out.println("(alpha1/alpha2)         : "
					+ format.format(alpha1) + "/" + format.format(alpha2));
			configurationDescription += "Metric: "
					+ GenerateDistribution.METRIC_EVENTS_AND_MAAC_DISTANCE
					+ " - EVENTS AND MAAC Metric\\n";
			configurationDescription += "(alpha1/alpha2): "
					+ format.format(alpha1) + "/" + format.format(alpha2)
					+ "\\n";
			break;
		case GenerateDistribution.METRIC_EVENTS_AND_MAAC_DISTANCE_MULTIPLEASSOCIATIONS:
			System.out.println("Metric                  : "
					+ GenerateDistribution.METRIC_EVENTS_AND_MAAC_DISTANCE_MULTIPLEASSOCIATIONS
					+ " - EVENTS AND MAAC-MULT. ASSOCIATIONS Metric");
			System.out.println("(alpha1/alpha2)         : "
					+ format.format(alpha1) + "/" + format.format(alpha2));
			configurationDescription += "Metric: "
					+ GenerateDistribution.METRIC_EVENTS_AND_MAAC_DISTANCE_MULTIPLEASSOCIATIONS
					+ " - EVENTS AND MAAC-MULT. ASSOCIATIONS Metric\\n";
			configurationDescription += "(alpha1/alpha2): "
					+ format.format(alpha1) + "/" + format.format(alpha2)
					+ "\\n";
			break;
		}

		if (model == GenerateDistribution.DEFAULT_MODEL) {
			System.out.println("Alpha cospeciation      : "
					+ format.format(alphaCospeciation));
			System.out.println("Alpha duplication       : "
					+ format.format(alphaDuplication));
			System.out.println("Alpha switch            : "
					+ format.format(alphaSwitch));
			System.out.println("Alpha loss              : "
					+ format.format(alphaLoss));
			configurationDescription += "Alpha (cosp.,dupl.,switch,loss): ("
					+ format.format(alphaCospeciation) + ","
					+ format.format(alphaDuplication) + ","
					+ format.format(alphaSwitch) + ","
					+ format.format(alphaLoss) + ")\\n";
			System.out.println("Root mapping probability: "
					+ format.format(rootMappingProbability));
			configurationDescription += "Root mapping probability: "
					+ format.format(rootMappingProbability);
		} else {
			System.out.println("Alpha cospeciation      : "
					+ format.format(alphaCospeciation));
			System.out.println("Alpha duplication       : "
					+ format.format(alphaDuplication));
			System.out.println("Alpha switch            : "
					+ format.format(alphaSwitch));
			configurationDescription += "Alpha (cosp.,dupl.,switch): ("
					+ format.format(alphaCospeciation) + ","
					+ format.format(alphaDuplication) + ","
					+ format.format(alphaSwitch) + ")";
		}

		

	}

	/**
	 * Read vectors (probabilities and distances) from the give file.
	 * 
	 * @param file
	 *            The file to be read.
	 * @return A list of vectors.
	 * @throws IOException
	 *             When something wrong happens while reading the file.
	 */
	private ArrayList<double[]> readVectorsFromFile(String file,
			ArrayList<double[]> vectors) throws IOException {
		File fileObject = new File(file);
		if (fileObject.exists()) {

			BufferedReader reader = new BufferedReader(new FileReader(
					fileObject));

			String line = "";
			while (line != null) {

				line = reader.readLine();

				if (line != null) {

					String[] elements = line.trim().split("\t");
					double[] vector = new double[NCOLUMNS];
					for (int i = 0; i < NCOLUMNS; i++) {
						vector[i] = Double.parseDouble(elements[i]);
					}
					vectors.add(vector);

				}

			}

			reader.close();

		}

		return vectors;

	}

	/**
	 * Run the ABC SMC process.
	 * 
	 * @throws IOException
	 *             When some IO operation fails.
	 */
	private void runABC_SMC() throws IOException {
		ArrayList<double[]> previousStep = new ArrayList<double[]>();
		ArrayList<double[]> lastStep = new ArrayList<double[]>();

		int lastIteration = 1;
		int populationQuantile = (int) Math.ceil(tolerances[0]
				* priorDistributionSize);
		
		if (continueProcess) {
			/* Discover where the process had stopped. */
			lastIteration = numberOfIterations;
			while (lastIteration > 1) {
				File file = new File(inputFileName + ".simul.round_"
						+ lastIteration + ".csv");
				if (file.exists()) {
					break;
				}
				lastIteration--;
			}

			/* Read the list of simulated vectors of the last round. */
			lastStep = readVectorsFromFile(inputFileName + ".simul.round_"
					+ lastIteration + ".csv", lastStep);

			if (lastIteration != 1) {
				/* Read the list of simulated vectors of the before last round. */
				readVectorsFromFile(inputFileName + ".simul.round_"
						+ (lastIteration - 1) + ".csv", previousStep);
			}

		}

		double epsilon = 1.0;
		ArrayList<double[]> quantile = null;
		ArrayList<double[]> acceptedVectors = null;

		if (lastIteration == 1) {

			int vectorsToGenerate = priorDistributionSize;
			if (lastStep != null) {
				vectorsToGenerate -= lastStep.size();
			}

			ArrayList<double[]> priorDistribution = lastStep;

			/* Check if the number of vectors are OK. */
			if (priorDistribution.size() > priorDistributionSize) {
				throw new RuntimeException(
						"Prior distribution computed previously has more vectors than the expected number ("
								+ priorDistribution.size()
								+ " > "
								+ priorDistributionSize + ")");
			}

			PrintWriter writer = openCSV(inputFileName + ".simul.round_1.csv");
			SynchronizedVectorManager manager = new SynchronizedVectorManager(
					priorDistributionSize, writer, format);
			for (double[] vector : priorDistribution) {
				manager.registerVector(vector);
			}

			printMessage("Running round 1 ...");
			printMessage("Generating prior distribution ...");
			if (continueProcess && priorDistribution.size() > 0) {
				printMessage("Already computed: " + priorDistribution.size()
						+ " vectors");
			}

			if (vectorsToGenerate > 0) {

				PriorDistribution[] threads = new PriorDistribution[numberOfThreads];
				for (int i = 0; i < numberOfThreads; i++) {
					/*threads[i] = new PriorDistribution(manager,
							priorDistributionSize, numberOfTrees, hostTree,
							parasiteTree, mappingParasiteHost, model, metric,
							maximumNumberOfTreesFactor, cyclicityTest,
							alphaCospeciation, alphaDuplication, alphaLoss,
							alphaSwitch, alpha1, alpha2, rootMappingProbability);
							*/
					threads[i] = new PriorDistribution(manager,
							priorDistributionSize, numberOfTrees, hostTree,
							parasiteTree, mappingParasiteHost, model, metric,
							maximumNumberOfTreesFactor, cyclicityTest,
							alphaCospeciation, alphaDuplication, alphaLoss,
							alphaSwitch, alpha1, alpha2, rootMappingProbability);
					threads[i].start();
				}

				for (int i = 0; i < numberOfThreads; i++) {
					try {
						threads[i].join();
					} catch (InterruptedException e) {
						throw new RuntimeException("Thread error: join " + i);
					}
				}

				priorDistribution = manager.getVectors();

				if (priorDistribution.size() < priorDistributionSize) {
					throw new RuntimeException(
							"Failure while generating prior distribution.");
				}

			}
			System.out.println();
			closeCSV(writer);
			

			printMessage("Executing ABC rejection method ...");
			acceptedVectors = new ArrayList<double[]>();
			epsilon = rejection(priorDistribution, acceptedVectors,
					populationQuantile);
			writer = openCSV(inputFileName + ".accep.round_1.csv");
			writeVectorToCSV(writer, acceptedVectors);
			closeCSV(writer);

			
			if (plot) {
				printMessage("Generating plots ...");
				plot(1, epsilon);
			}

			quantile = new ArrayList<double[]>();

			lastIteration = 2;

		} else {

			printMessage("Running round " + (lastIteration - 1) + " ...");
			printMessage("Executing ABC rejection method ...");

			/* Check if the number of vectors are OK. */
			if (lastIteration == 2
					&& previousStep.size() != priorDistributionSize) {
				throw new RuntimeException(
						"Prior distribution computed previously has more vectors than the expected number ("
								+ previousStep.size()
								+ " > "
								+ priorDistributionSize + ")");
			} else if (lastIteration > 2
					&& previousStep.size() != populationQuantile) {
				throw new RuntimeException(
						"Quantile distribution computed previously for round "
								+ (lastIteration - 1)
								+ " has more vectors than the expected number ("
								+ previousStep.size() + " > "
								+ populationQuantile + ")");
			}

			int populationToAccept = populationQuantile;
			if (lastIteration > 2) {
				populationToAccept = (int) Math
						.ceil(tolerances[lastIteration - 2]
								* populationQuantile);
			}

			acceptedVectors = new ArrayList<double[]>();
			epsilon = rejection(previousStep, acceptedVectors,
					populationToAccept);

			PrintWriter writer = openCSV(inputFileName + ".accep.round_"
					+ (lastIteration - 1) + ".csv");
			writeVectorToCSV(writer, acceptedVectors);
			closeCSV(writer);

			if (plot) {
				printMessage("Generating plots ...");
				plot(lastIteration - 1, epsilon);
			}

			quantile = lastStep;

			/* Just check if the number of vectors are ok. */
			if (quantile.size() > populationQuantile) {
				throw new RuntimeException(
						"Quantile distribution computed previously for round "
								+ lastIteration
								+ " has more vectors than the expected number ("
								+ quantile.size() + " > " + populationQuantile
								+ ")");
			}

		}

		int roundToCluster = numberOfIterations;
		for (int iteration = lastIteration; iteration <= numberOfIterations; iteration++) {

			printMessage("Running round " + (iteration) + " ...");
			printMessage("Generating population quantile ...");
			int maximumNumberOfDiscardedVectors = maximumNumberOfDiscardedVectorsFactor
					* populationQuantile;
			
			if (continueProcess && quantile.size() > 0) {
				printMessage("Already computed: " + quantile.size()
						+ " vectors");
				maximumNumberOfDiscardedVectors = maximumNumberOfDiscardedVectorsFactor
						* (populationQuantile - quantile.size());
			}

			PrintWriter writer = openCSV(inputFileName + ".simul.round_"
					+ iteration + ".csv");

			SynchronizedVectorManager manager = new SynchronizedVectorManager(
					populationQuantile, writer, format);

			for (double[] vector : quantile) {
				manager.registerVector(vector);
			}
			
			QuantileDistribution[] threads = new QuantileDistribution[numberOfThreads];
			for (int i = 0; i < numberOfThreads; i++) {
				/*threads[i] = new QuantileDistribution(manager, acceptedVectors,
						populationQuantile, numberOfTrees,
						maximumNumberOfDiscardedVectors, hostTree,
						parasiteTree, mappingParasiteHost, model, metric,
						maximumNumberOfTreesFactor, cyclicityTest, alpha1,
						alpha2, rootMappingProbability, perturbation, epsilon);
						*/
				threads[i] = new QuantileDistribution(manager, acceptedVectors,
						populationQuantile, numberOfTrees,
						maximumNumberOfDiscardedVectors, hostTree,
						parasiteTree, mappingParasiteHost, model, metric,
						maximumNumberOfTreesFactor, cyclicityTest, alpha1,
						alpha2, rootMappingProbability, perturbation, epsilon);
				threads[i].start();
			}

			for (int i = 0; i < numberOfThreads; i++) {
				try {
					threads[i].join();
				} catch (InterruptedException e) {
					throw new RuntimeException("Thread error: join " + i);
				}
			}

			quantile = manager.getVectors();
			int discarded = manager.getNumberOfRejectedVectors();
			boolean roundComplete = true;
			if (quantile.size() < populationQuantile) {
				if (discarded < maximumNumberOfDiscardedVectors) {
					throw new RuntimeException(
							"Unknown failure while generating quantile.");
				} else {
					roundComplete = false;
				}
			}

			System.out.println();
			closeCSV(writer);

			if (roundComplete) {
				/*
				 * In this case the round finished without any problem. we can
				 * run the rejection method without any problem.
				 */

				printMessage("Executing ABC rejection method ...");
				acceptedVectors = new ArrayList<double[]>();
				int populationToAccept = (int) Math
						.ceil(tolerances[iteration - 1] * populationQuantile);
				epsilon = rejection(quantile, acceptedVectors,
						populationToAccept);

				writer = openCSV(inputFileName + ".accep.round_" + iteration
						+ ".csv");
				writeVectorToCSV(writer, acceptedVectors);
				closeCSV(writer);

				if (plot) {
					printMessage("Generating plots ...");
					plot(iteration, epsilon);
				}

			} else {

				/*
				 * In this case, the round was aborted because it could not
				 * generate the necessary number of vectors in the population
				 * quantile.
				 */
				printMessage("[WARNING] Round " + (iteration) + " was aborted!");
				printMessage("[WARNING] Maximum number of discarded vectors was reached.");
				printMessage("[WARNING] Round failed to improve the sampled population (cannot improve epsilon).");
				if (iteration == 2) {
					printMessage("[WARNING] Failure occured at round 2.");
					printMessage("[WARNING] Please, use another set of parameters");
					printMessage("[WARNING] Try new population size (-N), number of trees (-M) or thresholds (-t)");
					roundToCluster = 0;
				} else {
					printMessage("[WARNING] Failure occured at round "
							+ (iteration) + ".");
					printMessage("[WARNING] Please, consider the results of round "
							+ (iteration - 1) + ".");
					roundToCluster = iteration - 1;
				}

				break;

			}

			quantile = new ArrayList<double[]>();

		}

		if (clustering) {
			if (roundToCluster != 0) {
				printMessage("Clustering accepted vectors of round "
						+ (roundToCluster) + "...");
				clusterAcceptedVectors(roundToCluster);
			} else {
				printMessage("[WARNING] Cannot run clustering procedure on the unfinished round 2.");
			}
		}

		printMessage("END!!!");

	}

	/**
	 * Run the ABC rejection method.
	 * 
	 * @param vectors
	 *            List of simulated vectors.
	 * @param accepted
	 *            Empty ArrayList that will receive the accepted vectors.
	 * @param numberOfVectorsToAccept
	 *            Number of vectors that must be accepted.
	 * @return Epsilon used in the rejection process.
	 */
	private double rejection(ArrayList<double[]> vectors,
			ArrayList<double[]> accepted, int numberOfVectorsToAccept) {

		int columnToEvaluate = COL_DISTANCE;

		double[] values = new double[vectors.size()];
		for (int i = 0; i < vectors.size(); i++) {
			double[] vector = vectors.get(i);
			values[i] = vector[columnToEvaluate];
		}

		double[] normalised = normalise(values);
		double[] normalisedAndSorted = Arrays.copyOf(normalised,
				normalised.length);
		Arrays.sort(normalisedAndSorted);

		double epsilon = 0.0;
		double normalisedEpsilon = normalisedAndSorted[numberOfVectorsToAccept - 1];
		for (int i = 0; i < normalised.length; i++) {
			if (normalised[i] == normalisedEpsilon) {
				epsilon = values[i];
				break;
			}
		}

		int index = 0;
		int[] indexes = new int[numberOfVectorsToAccept];
		for (int i = 0; i < normalised.length
				&& index < numberOfVectorsToAccept; i++) {
			if (normalised[i] < normalisedEpsilon) {
				indexes[index++] = i;
			}
		}

		for (int i = 0; i < normalised.length
				&& index < numberOfVectorsToAccept; i++) {
			if (normalised[i] == normalisedEpsilon) {
				indexes[index++] = i;
			}
		}

		for (int i : indexes) {
			accepted.add(vectors.get(i));
		}

		// now).

		return epsilon;
	}

	/**
	 * Normalise the vector values according with the Median Absolute Deviation.
	 * 
	 * @param values
	 *            Values to be normalised.
	 * @return The normalised values.
	 */
	private double[] normalise(double[] values) {
		double[] auxiliar = Arrays.copyOf(values, values.length);
		double median = Statistics.median(auxiliar);
		for (int i = 0; i < values.length; i++) {
			auxiliar[i] = Math.abs(values[i] - median);
		}
		double mad = Statistics.median(auxiliar);
		if (mad == 0) {
			return values;
		}
		for (int i = 0; i < values.length; i++) {
			auxiliar[i] = values[i] / mad;
		}
		return auxiliar;
	}

	/**
	 * Open a csv file for writing
	 * 
	 * @param fileName
	 *            File name
	 * @return A PrintWriter object.
	 * @throws IOException
	 *             When something wrong happens while creating the file.
	 */
	private PrintWriter openCSV(String fileName) throws IOException {
		return new PrintWriter(new File(fileName));
	}

	/**
	 * Close a csv file that was opened for writing.
	 * 
	 * @param writer
	 *            PrintWriter object.
	 */
	private void closeCSV(PrintWriter writer) {
		writer.close();
	}

	/**
	 * Write the list of vectors to the output file.
	 * 
	 * @param writer
	 *            PrintWriter object.
	 * @param vectors
	 *            List of vectors.
	 */
	private void writeVectorToCSV(PrintWriter writer,
			ArrayList<double[]> vectors) {
		for (double[] vector : vectors) {
			writer.println(formatVectorToPrint(vector));
		}
		writer.flush();
	}

	/**
	 * Format vector to print.
	 * 
	 * @param vector
	 *            Vector to be formatted.
	 * @return The String that represent the vector.
	 */
	private String formatVectorToPrint(double[] vector) {
		String[] str = new String[vector.length];
		for (int i = 0; i < vector.length; i++) {
			str[i] = format.format(vector[i]);
		}
		return Arrays.toString(str).replaceAll("[\\[\\],]", "")
				.replaceAll(" ", "\t");
	}

	/**
	 * Write the list of simulated and accepted vectors.
	 * 
	 * @param iteration
	 *            Iteration number
	 * @param simulated
	 *            List of simulated vectors.
	 * @param accepted
	 *            List of accepted vectors.
	 * @throws IOException
	 *             When something wrong happens while writing the files.
	 */
	private void plot(int iteration, double epsilon) throws IOException {

		String[] elements = inputFileName.split("/");
        String file = elements[elements.length - 1];

        String plotsFileName = inputFileName + ".plots.round_" + iteration + ".pdf";
        String simulationFileName = inputFileName + ".simul.round_" + iteration + ".csv";
        String acceptedFileName = inputFileName + ".accep.round_" + iteration + ".csv";
        String scriptFileName = inputFileName + ".script.round_" + iteration + ".R";
        String scriptLog = inputFileName + ".script.round_" + iteration + ".Rout";

        PrintWriter script = new PrintWriter(new File(scriptFileName));

        script.println("histPercent <- function(x, ...) {");
        script.println("  H <- hist(x, plot = FALSE, breaks=11)");
        script.println("  H$density <- with(H, 100 * density* diff(breaks)[1])");
        script.println("  labs <- paste(round(H$density), \"%\", sep=\"\")");
        script.println("  for (i in 1:length(labs)) {");
        script.println("    if (round(H$density[i]) == 0) {");
        script.println("      labs[i] = H$counts[i]");
        script.println("    }");
        script.println("  }");
        script.println("  plot(H, freq = FALSE, labels = labs, ylim=c(0, 1.08*max(H$density)),...)");
        script.println("}\n\n");

        script.println("simulation <- read.table(\"" + simulationFileName + "\", sep=\"\\t\", h=F)");
        script.println("accepted   <- read.table(\"" + acceptedFileName + "\", sep=\"\\t\", h=F)");

		if (model == GenerateDistribution.DEFAULT_MODEL) {
			script.println("names(simulation) <- c(\"p_cosp\", \"p_dup\", \"p_switch\", \"p_loss\", \"nb_vertical_Spread\", \"nb_horizontal_Spread\", \"dist\")");
			script.println("names(accepted) <- c(\"p_cosp\", \"p_dup\", \"p_switch\", \"p_loss\", \"nb_vertical_Spread\", \"nb_horizontal_Spread\", \"dist\")");
			//script.println("names(simulation) <- c(\"p_cosp\", \"p_dup\", \"p_switch\", \"p_loss\", \"dist\", \"maxMax_switch\", \"maxMean_switch\", \"meanMax_switch\", \"meanMean_switch\", \"nbMax_switch\", \"nbMean_switch\")");
			//script.println("names(accepted) <- c(\"p_cosp\", \"p_dup\", \"p_switch\", \"p_loss\", \"dist\", \"maxMax_switch\", \"maxMean_switch\", \"meanMax_switch\", \"meanMean_switch\", \"nbMax_switch\", \"nbMean_switch\")");
		} else {
			script.println("names(simulation) <- c(\"p_cosp\", \"p_dup\", \"p_switch\", \"f_loss\", \"dist\")");
			script.println("names(accepted) <- c(\"p_cosp\", \"p_dup\", \"p_switch\", \"f_loss\", \"dist\")");
		}

		script.println("pdf(file = \"" + plotsFileName
				+ "\", width = 16, height = 9, onefile = T)");

		script.println("par(mfrow=c(3,6))");

		String loopHeader = "for (i in 1:4) {";
		// if (model != GenerateDistribution.DEFAULT_MODEL) {
		// loopHeader = "for (i in c(1,2,3,5)) {";
		// }

		script.println(loopHeader);
		script.println("  histPercent(simulation[,i], main=\""
				+ (file + "\\nround " + iteration + " - simulation")
				+ "\", xlab=names(simulation)[i], cex.main=0.9, col=\"grey\")");
		script.println("}");

		// 
				
		script.println("if (length(unique(simulation$nb_vertical_Spread))==1) {"
				+ " hist(simulation$nb_vertical_Spread, breaks=c(unique(simulation$nb_vertical_Spread),unique(simulation$nb_vertical_Spread)+0.01), main=\""
				+ (file + "\\nround " + iteration + " - simulation")
				+ "\", xlab=\"vertical_Spread\", xlim=c(0,unique(simulation$nb_vertical_Spread)+1), cex.main=0.9, col=\"grey\")} else {");
		
		script.println("hist(simulation$nb_vertical_Spread,  main=\""
				+ (file + "\\nround " + iteration + " - simulation")
				+ "\", xlab=\"vertical_Spread\", cex.main=0.9, col=\"grey\")}");
		
		 // script.println("hist(simulation$nb_vertical_Spread,  main=\""
	    //+ (file + "\\nround " + iteration + " - simulation")
	    // + "\", xlab=\"vertical_Spread\", cex.main=0.9, col=\"grey\")");
		
		//script.println("hist(simulation$nb_horizontal_Spread,  main=\""
			//	+ (file + "\\nround " + iteration + " - simulation")
				//+ "\", xlab=\"horizontal_Spread\", cex.main=0.9, col=\"grey\")");
		
		script.println("if (length(unique(simulation$nb_horizontal_Spread))==1) {"
				+ " hist(simulation$nb_horizontal_Spread, breaks=c(unique(simulation$nb_horizontal_Spread),unique(simulation$nb_horizontal_Spread)+0.01), main=\""
				+ (file + "\\nround " + iteration + " - simulation")
				+ "\", xlab=\"horizontal_Spread\", xlim=c(0,unique(simulation$nb_horizontal_Spread)+1), cex.main=0.9, col=\"grey\")} else {");
		
		script.println("hist(simulation$nb_horizontal_Spread,  main=\""
				+ (file + "\\nround " + iteration + " - simulation")
				+ "\", xlab=\"horizontal_Spread\", cex.main=0.9, col=\"grey\")}");
		
		
		//script.println("hist(simulation$maxMax_switch,  main=\""
		//		+ (file + "\\nround " + iteration + " - simulation")
		//		+ "\", xlab=\"dist_jump\", breaks = length(seq(min(simulation$maxMax_switch),max(simulation$maxMax_switch))), cex.main=0.9, cex=0.5, col=\"grey\")");

		script.println(loopHeader);
		script.println("  histPercent(accepted[,i], main=\""
				+ (file + "\\nround " + iteration)
				+ (" - epsilon = " + format.format(epsilon))
				+ "\", xlab=names(accepted)[i], cex.main=0.9, col=\"grey\")");
		script.println("}");
		
		
		script.println("if (length(unique(accepted$nb_vertical_Spread))==1) {"
				+ " hist(accepted$nb_vertical_Spread, breaks=c(unique(accepted$nb_vertical_Spread),unique(accepted$nb_vertical_Spread)+0.01), main=\""
				+ (file + "\\nround " + iteration)
				+ "\", xlab=\"vertical_Spread\", xlim=c(0,unique(accepted$nb_vertical_Spread)+1), cex.main=0.9, col=\"grey\")} else {");
		
		script.println("hist(accepted$nb_vertical_Spread,  main=\""
				+ (file + "\\nround " + iteration)
				+ "\", xlab=\"vertical_Spread\", cex.main=0.9, col=\"grey\")}");
		
		
		
		//script.println("hist(accepted$nb_vertical_Spread,  main=\""
			//	+ (file + "\\nround " + iteration )
				//+ "\", xlab=\"vertical_Spread\", cex.main=0.9, col=\"grey\")");
		
		
		script.println("if (length(unique(accepted$nb_horizontal_Spread))==1) {"
				+ " hist(accepted$nb_horizontal_Spread, breaks=c(unique(accepted$nb_horizontal_Spread),unique(accepted$nb_horizontal_Spread)+0.01), main=\""
				+ (file + "\\nround " + iteration + " - simulation")
				+ "\", xlab=\"horizontal_Spread\", xlim=c(0,unique(accepted$nb_horizontal_Spread)+1), cex.main=0.9, col=\"grey\")} else {");
		
		script.println("hist(accepted$nb_horizontal_Spread,  main=\""
				+ (file + "\\nround " + iteration)
				+ "\", xlab=\"horizontal_Spread\", cex.main=0.9, col=\"grey\")}");
		
		
		//script.println("hist(accepted$nb_horizontal_Spread,  main=\""
			//	+ (file + "\\nround " + iteration )
				//+ "\", xlab=\"horizontal_Spread\", cex.main=0.9, col=\"grey\")");

		//script.println("hist(accepted$maxMax_switch,  main=\""
		//		+ (file + "\\nround " + iteration)
		//		+ (" - epsilon = " + format.format(epsilon))
		//		+ "\", xlab=\"dist_jump\", breaks = length(seq(min(accepted$maxMax_switch),max(accepted$maxMax_switch))), cex.main=0.9, cex=0.5, col=\"grey\")");
		

		script.println("histPercent(simulation$dist,  main=\""
				+ (file + "\\nround " + iteration + " - simulation")
				+ "\", xlab=\"Distance\", cex.main=0.9, cex=0.5, col=\"grey\")");
		script.println("histPercent(accepted$dist, main=\""
				+ (file + "\\nround " + iteration)
				+ (" - epsilon = " + format.format(epsilon))
				+ "\", xlab=\"Distance\", cex.main=0.9, cex=0.5, col=\"grey\")");
		

		script.println("plot(x=c(0,100), y=c(0,100), type=\"n\", xaxt=\"n\", yaxt=\"n\", ylab=\"\", xlab=\"\")");
		script.println("text(5,35, pos=4, cex=0.80, \""
				+ configurationDescription + "\")");

		script.println("dev.off()");

		script.flush();
		script.close();

		runR(scriptFileName, scriptLog);

	}

	/**
	 * Run the clustering procedure.
	 * 
	 * @param roundToCluster
	 *            Identifier of the round that must be clustered.
	 * @throws IOException
	 *             When something wrong happens while writing the files.
	 */
	private void clusterAcceptedVectors(int roundToCluster) throws IOException {

		String scriptFileName = inputFileName + ".script.cluster.R";
		String scriptLog = inputFileName + ".script.cluster.Rout";
		clusteringModel(scriptFileName, roundToCluster);
		runR(scriptFileName, scriptLog);

	}

	/**
	 * Call R to run the ABC script.
	 * 
	 * @param rFileName
	 *            Script file name.
	 * @param rFileOut
	 *            Script output file name.
	 * @throws IOException
	 *             When something goes wrong while writing the files.
	 */
	private void runR(String rFileName, String rFileOut) throws IOException {

		Process process = Runtime.getRuntime().exec(
				"R CMD BATCH --no-save --no-restore " + rFileName + " "
						+ rFileOut);
		try {
			process.waitFor();
		} catch (InterruptedException e) {
			e.printStackTrace();
		}

		File script = new File(rFileName);
		if (script.exists()) {
			script.delete();
		}

		File log = new File(rFileOut);
		if (log.exists()) {
			log.delete();
		}

	}

	/**
     * Create the script file for the clustering procedure.
     * 
     * @param scriptFileName
     *        Script file name.
     * @param roundToCluster
     *        Identifier of the round that must be clustered.
     * @throws IOException
     *         When something goes wrong while writing the files.
     */
    private void clusteringModel(String scriptFileName, int roundToCluster) throws IOException {

        String[] elements = inputFileName.split("/");
        String file = elements[elements.length - 1];
        String acceptedFileName = inputFileName + ".accep.round_" + roundToCluster + ".csv";
        String clusterPrefix = inputFileName + ".clust.round_" + roundToCluster + ".";
        String plotsFileName = inputFileName + ".clusters.round_" + roundToCluster + ".pdf";

        PrintWriter script = new PrintWriter(new File(scriptFileName));

        script.println("library(\"stats\")");
        script.println("library(\"dynamicTreeCut\")\n");

        script.println("chisq <- function(x, y) {");
        script.println("\ta <- (x + y) / 2");
        script.println("\tb <- (x - y) ** 2");
        script.println("\treturn (sum( b[a > 0] / a[a > 0] ))");
        script.println("}\n");

        script.println("chisqmatrix <- function(data) {");
        script.println("\tnr <- nrow(data)");
        script.println("\tnc <- ncol(data)");
        script.println("\tm <- NULL");
        script.println("\tfor (i in 1:nr) {");
        script.println("\t\tv <- apply(data[,1:nc], 1, function(x) { chisq(as.vector(t(data[i,1:nc])), as.vector(t(x[1:nc]))) })");
        script.println("\t\tif (is.null(m)) {");
        script.println("\t\t\tm <- matrix(v, ncol=nr)");
        script.println("\t\t} else {");
        script.println("\t\t\tm <- rbind(m, v)");
        script.println("\t\t}");
        script.println("\t}");
        script.println("\tm <- as.data.frame(m)");
        script.println("\tnames(m) <- 1:nr");
        script.println("\trownames(m) <- 1:nr");
        script.println("\treturn (as.dist(m))");
        script.println("}\n");

        script.println("accepted   <- read.table(\"" + acceptedFileName + "\", sep=\"\\t\", h=F)");

        if (model == GenerateDistribution.DEFAULT_MODEL) {
            clusteringDefaultModel(script, file, clusterPrefix, plotsFileName);
        } else {
            clusteringCoalescenttModel(script, file, clusterPrefix, plotsFileName);
        }

        script.flush();
        script.close();
    }


    private void clusteringDefaultModel(PrintWriter script,
                                        String file,
                                        String clusterPrefix,
                                        String plotsFileName) {
    	script.println("names(accepted) <- c(\"p_cosp\", \"p_dup\", \"p_switch\", \"p_loss\", \"nb_vertical_Spread\", \"nb_horizontal_Spread\", \"dist\")");
        //script.println("names(accepted) <- c(\"p_cosp\", \"p_dup\", \"p_switch\", \"p_loss\", \"dist\")");
    	//erase horizontal and vertical spread
    	script.println("accepted<-accepted[-5]");
    	script.println("accepted<-accepted[-5]");
        script.println("distanceMatrix <- chisqmatrix(accepted[,1:4])");
       
        script.println("cluster <- hclust(distanceMatrix)");
        
        script.println("groups <- as.vector(cutreeDynamic(cluster, distM=as.matrix(distanceMatrix), minClusterSize = 2))");

        script.println("header <- as.data.frame(matrix(c(\"-Cosp-\", \"-Dupl-\", \"-Swit-\", \"-Loss-\", \"-Dist-\"), nrow=1))");
        script.println("names(header) <- c(\"p_cosp\", \"p_dup\", \"p_switch\", \"p_loss\", \"dist\")");

        script.println("for (group in sort(unique(groups))) {");
        script.println("  aux <- accepted[groups == group, ]");
        script.println("  aux <- aux[order(aux$p_cosp), ]");
        script.println("  ids <- rownames(aux)");
        script.println("  co <- as.vector(summary(aux$p_cosp))");
        script.println("  du <- as.vector(summary(aux$p_dup))");
        script.println("  sw <- as.vector(summary(aux$p_switch))");
        script.println("  lo <- as.vector(summary(aux$p_loss))");
        script.println("  di <- as.vector(summary(aux$dist))");
        script.println("  aux <- format(aux, nsmall=4, digits=1, scientific=F)");
        script.println("  if (nrow(aux) > 1) {");
        script.println("    ids <- c(\"----\", ids, \"Summar\", \"Min\", \"Q1\", \"Med\", \"Mean\", \"Q3\", \"Max\", \"Normal\", \"Nmed\", \"Nmean\")");
        script.println("    aux <- rbind(aux, rep(\"------\", 5))");
        script.println("    for (i in 1:length(co)) {");
        script.println("      aux <- rbind(aux, format(c(co[i], du[i], sw[i], lo[i], di[i]), nsmall=4, digits=1, scientific=F))");
        script.println("    }");
        script.println("    aux <- rbind(aux, rep(\"------\", 5))");
        script.println("    x <- sum(c(co[3], du[3], sw[3], lo[3]))");
        script.println("    x <- format(c(co[3], du[3], sw[3], lo[3])/x, nsmall=4, digits=1, scientific=F)");
        script.println("    aux <- rbind(aux, c(x, \"------\"))");
        script.println("    x <- sum(c(co[4], du[4], sw[4], lo[4]))");
        script.println("    x <- format(c(co[4], du[4], sw[4], lo[4])/x, nsmall=4, digits=1, scientific=F)");
        script.println("    aux <- rbind(aux, c(x, \"------\"))");
        script.println("  } else {");
        script.println("    ids <- c( \"----\", ids)");
        script.println("  }");
        script.println("  aux <- rbind(header, aux)");
        script.println("  rownames(aux) <- ids");
        script.println("  fileName <- paste(\"" + clusterPrefix + "\", group, \".csv\", sep=\"\")");
        script.println("  write.table(aux, file=fileName, row.names=T, col.names=F, sep=\"\\t\", quote=F)");
        script.println("}");

        if (plot) {
            script.println("library(\"ade4\")");
            script.println("pdf(file = \"" + plotsFileName
                           + "\", width = 8, height = 11, onefile = T)");
            script.println("minD            <- min(accepted$dist)");
            script.println("maxD            <- max(accepted$dist)");
            script.println("description <- paste(\""
                           + file
                           + "\", \"\\nMin distance (red) = \", minD, \"\\nMax distance (yellow)= \", maxD, sep=\"\")");
            script.println("sorteddistances <- sort(unique(accepted$dist))");
            script.println("distancecolors  <- heat.colors(length(sorteddistances))");
            script.println("clustercolors   <- rainbow(length(unique(groups)))");
            script.println("dudipco         <- dudi.pco(distanceMatrix, scannf = FALSE)");
            script.println("s               <- s.class(dudipco$li, as.factor(groups), clabel=0, cstar=1, cellipse=1.5, cpoint=0, col=rep(\"black\", length(groups)))");
            script.println("if (length(sorteddistances) > 1) {");
            script.println("    description <- paste(\""
                           + file
                           + "\", \"\\nMin distance (red) = \", minD, \"\\nMax distance (yellow)= \", maxD, sep=\"\")");
            script.println("    for (i in 1:length(sorteddistances)) {");
            script.println("        d                <- sorteddistances[i]");
            script.println("        indexes          <- which(accepted$dist == d)");
            script.println("        pointcolors      <- c(\"#00000000\", distancecolors[i])");
            script.println("        factors          <- rep(0, length(groups))");
            script.println("        factors[indexes] <- 1");
            script.println("        s.class(dudipco$li, as.factor(factors), clabel=0, cstar=0, axesell=F, cellipse=0, cpoint=2.0, col=pointcolors, add.plot=T)");
            script.println("    }");
            script.println("} else {");
            script.println("    description <- paste(\""
                           + file
                           + "\", \"\\nMin distance (red) = \", minD, \"\\nMax distance (red)= \", maxD, sep=\"\")");
            script.println("    pointcolors      <- distancecolors[1]");
            script.println("    s.class(dudipco$li, as.factor(rep(1, length(groups))), clabel=0, cstar=0, axesell=F, cellipse=0, cpoint=2.0, col=pointcolors, add.plot=T)");
            script.println("}");
            script.println("s.class(dudipco$li, as.factor(groups), clabel=1.2, cstar=0, axesell=F, cellipse=0, cpoint=0, col=rep(\"black\", length(groups)), sub=description, add.plot=T, csub=1.25)");
            script.println("dev.off()");
        }

    }



	private void clusteringCoalescenttModel(PrintWriter script, String file,
			String clusterPrefix, String plotsFileName) {

		script.println("names(accepted) <- c(\"p_cosp\", \"p_dup\", \"p_switch\", \"f_loss\", \"dist\")");

		script.println("distanceMatrix <- chisqmatrix(accepted[,1:3])");
		script.println("cluster <- hclust(distanceMatrix)");

		script.println("groups <- as.vector(cutreeDynamic(cluster, distM=as.matrix(distanceMatrix), minClusterSize = 2))");

		script.println("header <- as.data.frame(matrix(c(\"-Cosp-\", \"-Dupl-\", \"-Swit-\", \"-Loss-\", \"-Dist-\"), nrow=1))");
		script.println("names(header) <- c(\"p_cosp\", \"p_dup\", \"p_switch\", \"f_loss\", \"dist\")");

		script.println("for (group in sort(unique(groups))) {");
		script.println("  aux <- accepted[groups == group, ]");
		script.println("  aux <- aux[order(aux$p_cosp), ]");
		script.println("  ids <- rownames(aux)");
		script.println("  co <- as.vector(summary(aux$p_cosp))");
		script.println("  du <- as.vector(summary(aux$p_dup))");
		script.println("  sw <- as.vector(summary(aux$p_switch))");
		script.println("  lo <- as.vector(summary(aux$f_loss))");
		script.println("  di <- as.vector(summary(aux$dist))");
		script.println("  aux <- format(aux, nsmall=4, digits=1, scientific=F)");
		script.println("  if (nrow(aux) > 1) {");
		script.println("    ids <- c(\"----\", ids, \"Summar\", \"Min\", \"Q1\", \"Med\", \"Mean\", \"Q3\", \"Max\", \"Normal\", \"Nmed\", \"Nmean\")");
		script.println("    aux <- rbind(aux, rep(\"------\", 5))");
		script.println("    for (i in 1:length(co)) {");
		script.println("      aux <- rbind(aux, format(c(co[i], du[i], sw[i], lo[i], di[i]), nsmall=4, digits=1, scientific=F))");
		script.println("    }");
		script.println("    aux <- rbind(aux, rep(\"------\", 5))");
		script.println("    x <- sum(c(co[3], du[3], sw[3]))");
		script.println("    x <- format(c(co[3], du[3], sw[3])/x, nsmall=4, digits=1, scientific=F)");
		script.println("    aux <- rbind(aux, c(x, \"------\", \"------\"))");
		script.println("    x <- sum(c(co[4], du[4], sw[4]))");
		script.println("    x <- format(c(co[4], du[4], sw[4])/x, nsmall=4, digits=1, scientific=F)");
		script.println("    aux <- rbind(aux, c(x, \"------\", \"------\"))");
		script.println("  } else {");
		script.println("    ids <- c( \"----\", ids)");
		script.println("  }");
		script.println("  aux <- rbind(header, aux)");
		script.println("  rownames(aux) <- ids");
		script.println("  fileName <- paste(\"" + clusterPrefix
				+ "\", group, \".csv\", sep=\"\")");
		script.println("  write.table(aux, file=fileName, row.names=T, col.names=F, sep=\"\\t\", quote=F)");
		script.println("}");

		if (plot) {
			script.println("library(\"ade4\")");
			script.println("pdf(file = \"" + plotsFileName
					+ "\", width = 8, height = 11, onefile = T)");
			script.println("minF            <- min(accepted$f_loss)");
			script.println("maxF            <- max(accepted$f_loss)");
			script.println("description <- paste(\""
					+ file
					+ "\", \"\\nMin loss frequency (red) = \", minF, \"\\nMax loss frequency (yellow)= \", maxF, sep=\"\")");
			script.println("sortedfrequencies <- sort(unique(accepted$f_loss))");
			script.println("frequencycolors   <- heat.colors(length(sortedfrequencies))");
			script.println("clustercolors     <- rainbow(length(unique(groups)))");
			script.println("dudipco           <- dudi.pco(distanceMatrix, scannf = FALSE)");
			script.println("s                 <- s.class(dudipco$li, as.factor(groups), clabel=0, cstar=1, cellipse=1.5, cpoint=0, col=rep(\"black\", length(groups)))");
			script.println("if (length(sortedfrequencies) > 1) {");
			script.println("    description <- paste(\""
					+ file
					+ "\", \"\\nMin loss frequency (red) = \", minF, \"\\nMax loss frequency (yellow)= \", maxF, sep=\"\")");
			script.println("    for (i in 1:length(sortedfrequencies)) {");
			script.println("        f                <- sortedfrequencies[i]");
			script.println("        indexes          <- which(accepted$f_loss == f)");
			script.println("        pointcolors      <- c(\"#00000000\", frequencycolors[i])");
			script.println("        factors          <- rep(0, length(groups))");
			script.println("        factors[indexes] <- 1");
			script.println("        s.class(dudipco$li, as.factor(factors), clabel=0, cstar=0, axesell=F, cellipse=0, cpoint=2.0, col=pointcolors, add.plot=T)");
			script.println("    }");
			script.println("} else {");
			script.println("    description <- paste(\""
					+ file
					+ "\", \"\\nMin loss frequency (red) = \", minF, \"\\nMax loss frequency (red)= \", maxF, sep=\"\")");
			script.println("    pointcolors      <- frequencycolors[1]");
			script.println("    s.class(dudipco$li, as.factor(rep(1, length(groups))), clabel=0, cstar=0, axesell=F, cellipse=0, cpoint=2.0, col=pointcolors, add.plot=T)");
			script.println("}");
			script.println("s.class(dudipco$li, as.factor(groups), clabel=1.2, cstar=0, axesell=F, cellipse=0, cpoint=0, col=rep(\"black\", length(groups)), sub=description, add.plot=T, csub=1.25)");
			script.println("dev.off()");
		}

	}

	/**
	 * Print the given message with a timestamp.
	 * 
	 * @param message
	 *            Message to be printed.
	 */
	private void printMessage(String message) {
		calendar.setTimeInMillis(System.currentTimeMillis());
		System.out.println(dateFormat.format(calendar.getTime()) + " - "
				+ message);
	}

}
