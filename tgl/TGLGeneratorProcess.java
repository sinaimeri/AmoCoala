/* 
############################################################################
#                                Coala                                     #
#         CO-evolution Assessment by a Likelihood-free Approach            #
############################################################################
#                                                                          #
# Copyright INRIA 2018                                                     #
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

package  tgl;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Random;

import cycle.CyclicityTest;
import generator.Association;
import generator.IParasiteGenerator;
import generator.IScenario;
import generator.defaultmodel.DefaultModel;
import trees.Tree;
import trees.TreeNode;
import util.NexusFileParserException;

public class TGLGeneratorProcess {

    /* -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:- */
    /* Constants which define the default values for the input parameters of the process. */
    public static final double DEFAULT_COSPECIATION = 0.7;

    public static final double DEFAULT_DUPLICATION = 0.1;

    public static final double DEFAULT_HOSTSWITCH = 0.1;

    public static final int DEFAULT_NUMBER_OF_TREES = 10;

    public static final double DEFAULT_EUCLIDEAN_DISTANCE_3_THRESHOLD = 1e-10;

    public static final double DEFAULT_EUCLIDEAN_DISTANCE_4_THRESHOLD = 0.05;

    public static final double DEFAULT_MAXIMUM_SIZE_FACTOR = 2.0;

    public static final int DEFAULT_CYCLICITY_TEST = CyclicityTest.DONATI;

    public static final double DEFAULT_ROOT_MAPPING_PROBABILITY = 1.0;

    private final char OPEN_BRACKET = '(';

    private final char CLOSE_BRACKET = ')';

    private final char CHILD_SEPARATOR = ',';

    /* -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:- */
    /* Input parameters of the process. */
    private double cospeciation;

    private double duplication;

    private double hostSwitch;

    private double cospeciationWithoutLoss;

    private double duplicationWithoutLoss;

    private double hostSwitchWithoutLoss;

    private double loss;

    private int numberOfTrees;

    private String inputFileName;

    private String outputDirectory;

    private String outputFileNamePrefix;

    private String formatString;

    private int padding;

    private int cyclicityTest;

    private boolean tabular;

    private double euclideanDistance3Threshold;

    private double euclideanDistance4Threshold;

    private double maximumSizeFactor;

    private double rootMappingProbability;

    /* -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:- */
    /* Work attributes */
    private Tree hostTree;

    private Tree hostTreeWithRHO;

    public Random random;

    private NumberFormat format;

    private int numberOfLeavesHostTree;

    private int minimumNumberOfLeaves;

    private int maximumNumberOfLeaves;


    /**
     * Constructor: Creates an TGLGeneratorProcess object filled with default parameters.
     */
    public TGLGeneratorProcess() {
        this.cospeciation = DEFAULT_COSPECIATION;
        this.duplication = DEFAULT_DUPLICATION;
        this.hostSwitch = DEFAULT_HOSTSWITCH;
        this.numberOfTrees = DEFAULT_NUMBER_OF_TREES;
        this.padding = 2;
        this.cyclicityTest = DEFAULT_CYCLICITY_TEST;
        this.inputFileName = null;
        this.outputFileNamePrefix = null;
        this.hostTree = null;
        this.hostTreeWithRHO = null;
        this.euclideanDistance3Threshold = DEFAULT_EUCLIDEAN_DISTANCE_3_THRESHOLD;
        this.euclideanDistance4Threshold = DEFAULT_EUCLIDEAN_DISTANCE_4_THRESHOLD;
        this.maximumSizeFactor = DEFAULT_MAXIMUM_SIZE_FACTOR;
        this.rootMappingProbability = DEFAULT_ROOT_MAPPING_PROBABILITY;
        this.random = new Random();
        this.format = NumberFormat.getInstance(Locale.US);
        this.format.setMinimumFractionDigits(4);
        this.format.setMaximumFractionDigits(4);
    }


    /* -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:- */
    /* Set methods for process configuration. */
    public void setCospeciation(double probability) {
        this.cospeciation = probability;
    }


    public void setDuplication(double probability) {
        this.duplication = probability;
    }


    public void setHostSwitch(double probability) {
        this.hostSwitch = probability;
    }


    public void setInputFile(String file) {
        this.inputFileName = file;
        File fileName = new File(inputFileName);
        File absolute = new File(fileName.getAbsolutePath());
        this.outputDirectory = absolute.getParent();
    }


    public void setOutputFileNamePrefix(String prefix) {
        this.outputFileNamePrefix = prefix;
    }


    public void setNumberOfTrees(int number) {
        this.numberOfTrees = number;
        this.padding = (int)Math.ceil(Math.log10(numberOfTrees + 0.1));
    }


    public void setCyclicityTest(int cyclicityTest) {
        this.cyclicityTest = cyclicityTest;
    }


    public void setTabular(boolean tabular) {
        this.tabular = tabular;
    }


    public void setEsuclideanDistance3Threshold(double threshold) {
        this.euclideanDistance3Threshold = threshold;
    }


    public void setEsuclideanDistance4Threshold(double threshold) {
        this.euclideanDistance4Threshold = threshold;
    }


    public void setMaximumSizeFactor(double factor) {
        this.maximumSizeFactor = factor;
    }


    public double getEventProbabilitiesSum() {
        return cospeciation + duplication + hostSwitch;
    }


    public void setRootMappingProbability(double probability) {
        this.rootMappingProbability = probability;
    }


    /* -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:- */
    /**
     * Execute the ABC SMC process.
     * 
     * @throws IOException
     *         When something goes wrong while reading the file.
     * @throws NexusFileParserException
     *         When something goes wrong while parsing the file.
     */
    public void run() throws IOException, NexusFileParserException {
        readNewickFile(inputFileName);
        generateTrees();
    }


    /**
     * Read the nexus file and set up some parameters.
     * 
     * @param fileName
     *        Nexus file name.
     * @throws IOException
     *         When something goes wrong while reading the file.
     * @throws NexusFileParserException
     *         When something goes wrong while parsing the file.
     */
    private void readNewickFile(String fileName) throws IOException {

        File file = new File(fileName);
        if (!file.exists()) {
            throw new IOException("File not found: " + fileName);
        }

        BufferedReader reader = new BufferedReader(new FileReader(file));
        try {
            String hostTreeString = reader.readLine();
            hostTree = parseNewickString(hostTreeString, "!H");
            hostTreeWithRHO = parseNewickString(hostTreeString, "!H");
        } finally {
            reader.close();
        }

        /* Add a fake root to the host tree. */
        TreeNode rho = hostTreeWithRHO.createNode();
        TreeNode root = hostTreeWithRHO.getRoot();
        hostTreeWithRHO.addChild(rho.getKey(), root.getKey(), 0);
        hostTreeWithRHO.setRoot(rho.getKey());

        /* Create the LCA tables for host and parasite trees. */
        hostTree.createLCATable();
        hostTreeWithRHO.createLCATable();

        /* Update the list of ancestors. */
        hostTree.updateAncestorList();
        hostTreeWithRHO.updateAncestorList();

        /* Update the list of descendants. */
        hostTree.updateDescendantList();
        hostTreeWithRHO.updateDescendantList();

        /* Save the number of leaves of the host tree. */
        numberOfLeavesHostTree = hostTree.getNumberOfLeafNodes();

    }


    /**
     * Generate the desired number of parasite trees.
     * 
     * @throws IOException
     *         When some IO operation fails.
     */
    private void generateTrees() throws IOException {

        double withoutLoss = cospeciation + duplication + hostSwitch;
        loss = 1.0 - withoutLoss;
        cospeciationWithoutLoss = cospeciation / withoutLoss;
        duplicationWithoutLoss = duplication / withoutLoss;
        hostSwitchWithoutLoss = hostSwitch / withoutLoss;

        formatString = outputDirectory + "/" + outputFileNamePrefix + "-%0" + padding + "d.tgl";

        maximumNumberOfLeaves = (int)Math.ceil(maximumSizeFactor * numberOfLeavesHostTree);

        IParasiteGenerator generator = new DefaultModel(hostTreeWithRHO,
        												null,
                                                        cospeciation,
                                                        duplication,
                                                        hostSwitch,
                                                        -1,
                                                        maximumNumberOfLeaves,
                                                        cyclicityTest,
                                                        true,
                                                        rootMappingProbability);

        int fileNumber = 1;
        int numberOfTooBigTrees = 0;
        int numberOfIterations = 0;
        while (fileNumber <= numberOfTrees) {
            IScenario scenario = generator.generateParasiteTree();
            if (!tabular) {
                if (writeTGLFile(scenario, fileNumber)) {
                    System.out.println(fileNumber);
                    fileNumber++;
                }
            } else {
                if (scenario != null) {
                    double eDistance = scenario.computeEuclideanDistance(cospeciation,
                                                                         duplication,
                                                                         hostSwitch,
                                                                         loss);
                    StringBuffer buffer = new StringBuffer();
                    buffer.append(format.format(eDistance));
                    buffer.append("\t");
                    buffer.append(scenario.getNumberOfLeaves());
                    buffer.append("\t");
                    buffer.append(format.format(scenario.getFrequencyOfCospeciations()));
                    buffer.append("\t");
                    buffer.append(format.format(scenario.getFrequencyOfDuplications()));
                    buffer.append("\t");
                    buffer.append(format.format(scenario.getFrequencyOfHostSwitches()));
                    buffer.append("\t");
                    buffer.append(format.format(scenario.getFrequencyOfLosses()));
                    buffer.append("\t");
                    buffer.append(scenario.getParasiteTree().toString(true, false));
                    System.out.println(buffer.toString());
                    fileNumber++;
                } else {
                    numberOfTooBigTrees++;
                }
            }
            numberOfIterations++;
            if (numberOfIterations > 10000 * numberOfTrees) {
                System.out.println("Limit reached.");
                break;
            }
        }

        if (tabular) {
            double percentage = (100.0 * numberOfTooBigTrees) / numberOfIterations;
            System.err.println(numberOfTooBigTrees + "\t" + numberOfIterations + "\t"
                               + format.format(percentage));
        }

    }


    /**
     * Write a TGL File for the given scenario.
     * 
     * @param scenario
     *        Scenario to be output into the TGL file.
     * @param fileNumber
     *        File number.
     * @return true if the scenario is considered valid to be writen into a file.
     * @throws IOException
     *         When something goes wrong while writing the file.
     */
    private boolean writeTGLFile(IScenario scenario, int fileNumber) throws IOException {

        if (scenario == null)
            return false;

        if (scenario.getNumberOfLeaves() < minimumNumberOfLeaves)
            return false;

        double eDistance3 = scenario.computeEuclideanDistance(cospeciationWithoutLoss,
                                                              duplicationWithoutLoss,
                                                              hostSwitchWithoutLoss);

        double eDistance4 = scenario.computeEuclideanDistance(cospeciation,
                                                              duplication,
                                                              hostSwitch,
                                                              loss);

        if (eDistance3 > euclideanDistance3Threshold) {
            return false;
        }

        if (eDistance4 > euclideanDistance4Threshold) {
            return false;
        }

        HashMap<String, Integer> labelCount = new HashMap<String, Integer>();
        for (TreeNode node: hostTree.getLeafNodes()) {
            //labelCount.put(node.getLabel(), 0);
        	labelCount.put(node.getLabel().toString(), 0);
        }

        HashMap<String, String> finalMapping = new HashMap<String, String>();
        for (Association a: scenario.getMapping()) {
            if (a.getParasite().isLeaf()) {
                TreeNode node = (TreeNode)a.getHostElement();
                int n = labelCount.get(node.getLabel()) + 1;
                //labelCount.put(node.getLabel(), n);
                labelCount.put(node.getLabel().toString(), n);
                String newLabel = node.getLabel() + "-" + n;
                List<String> newLabelList = new ArrayList<String>();
                newLabelList.add(newLabel);
                //a.getParasite().setLabel(newLabel);
                a.getParasite().setLabel(newLabelList);
                //finalMapping.put(newLabel, node.getLabel());
                finalMapping.put(newLabel, node.getLabel().toString());
            }
        }

        int nInternal = scenario.getParasiteTree().getLeafNodes().length - 1;

        int nC = scenario.getNumberOfCospeciations();
        int nD = scenario.getNumberOfDuplications();
        int nH = scenario.getNumberOfHostSwitches();
        int nL = scenario.getNumberOfLosses();
        int[] observed = {nC, nD, nH, nL};

        double nTotal = nInternal / (cospeciation + duplication + hostSwitch);

        String expC = String.format("%.3f", cospeciationWithoutLoss * nInternal);
        String expD = String.format("%.3f", duplicationWithoutLoss * nInternal);
        String expH = String.format("%.3f", hostSwitchWithoutLoss * nInternal);
        String expL = String.format("%.3f", nTotal * loss);
        String[] expected = {expC, expD, expH, expL};

        String pC = String.format("%.3f", cospeciation);
        String pD = String.format("%.3f", duplication);
        String pH = String.format("%.3f", hostSwitch);
        String pL = String.format("%.3f", loss);
        String[] probabilities = {pC, pD, pH, pL};

        String fCStr = String.format("%.3f", scenario.getFrequencyOfCospeciations());
        String fDStr = String.format("%.3f", scenario.getFrequencyOfDuplications());
        String fHStr = String.format("%.3f", scenario.getFrequencyOfHostSwitches());
        String fLStr = String.format("%.3f", scenario.getFrequencyOfLosses());
        String[] frequencies = {fCStr, fDStr, fHStr, fLStr};

        FileWriter writer = new FileWriter(new File(String.format(formatString, fileNumber)));
        writer.write("#NEXUS\n");
        writer.write("BEGIN HOST;\n");
        writer.write("\tTREE * Host1 = " + hostTree.toString(true, false) + ";\n");
        writer.write("ENDBLOCK;\n\n");
        writer.write("BEGIN PARASITE;\n");
        writer.write("\tTREE * Para1 = " + scenario.getParasiteTree().toString(true, false) + ";\n");
        writer.write("ENDBLOCK;\n\n");

        writer.write("BEGIN DISTRIBUTION;\n");
        writer.write("\tRANGE\n");

        boolean first = true;
        for (String key: finalMapping.keySet()) {
            if (!first)
                writer.write(",\n");
            first = false;
            writer.write("\t\t" + key + ": " + finalMapping.get(key));
        }
        writer.write("\n\t;\nEND;\n\n");

        writer.write("#Leaves        : " + (nInternal + 1) + "\n");
        writer.write("#Expected      : " + Arrays.toString(expected) + "\n");
        writer.write("#Observed      : " + Arrays.toString(observed) + "\n");
        writer.write("#Probabilities : " + Arrays.toString(probabilities) + "\n");
        writer.write("#Frequencies   : " + Arrays.toString(frequencies) + "\n");
        writer.write("#Euclidean3    : " + String.format("%.3f", eDistance3) + "\n");
        writer.write("#Euclidean4    : " + String.format("%.3f", eDistance4) + "\n");

        writer.close();

        return true;
    }


    /**
     * Parse the string to extract a tree.
     * 
     * @param treeString
     *        String to be parsed.
     * @param labelPrefix
     *        Label prefix for internal nodes.
     * @return The parsed tree.
     */
    private Tree parseNewickString(String treeString, String labelPrefix) {

        treeString = treeString.replaceAll(";", "");
        treeString = treeString.replaceAll(", +", ",");
        treeString = treeString.replaceAll("\\( +", "(");
        treeString = treeString.replaceAll("\\) +", ")");

        int cursor = 0;

        Tree tree = new Tree(true);
        TreeNode root = tree.getRoot();
        
        List<String> setLabelList = new ArrayList<String>();
        setLabelList.add(labelPrefix + Integer.toString(root.getKey()));
        //root.setLabel(labelPrefix + Integer.toString(root.getKey()));
        root.setLabel(setLabelList);

        int brackets = 0;

        TreeNode currentNode = root;
        while (cursor < treeString.length()) {
            if (treeString.charAt(cursor) == OPEN_BRACKET) {
                TreeNode child = tree.createNode();
                
                List<String> newLabelList = new ArrayList<String>();
                newLabelList.add(labelPrefix + Integer.toString(child.getKey()));
                //child.setLabel(labelPrefix + Integer.toString(child.getKey()));
                child.setLabel(newLabelList);
                
                tree.addChild(currentNode.getKey(), child.getKey(), 0);
                currentNode = child;
                cursor++;
                brackets++;
            } else if (treeString.charAt(cursor) == CHILD_SEPARATOR) {
                currentNode = currentNode.getParent();
                TreeNode child = tree.createNode();
                
                List<String> newLabelList = new ArrayList<String>();
                newLabelList.add(labelPrefix + Integer.toString(child.getKey()));
                //child.setLabel(labelPrefix + Integer.toString(child.getKey()));
                child.setLabel(newLabelList);
                
                tree.addChild(currentNode.getKey(), child.getKey(), 1);
                currentNode = child;
                cursor++;
            } else if (treeString.charAt(cursor) == CLOSE_BRACKET) {
                currentNode = currentNode.getParent();
                cursor++;
                brackets--;
                if (brackets < 0) {
                    throw new RuntimeException("Malformed tree string :\n" + treeString + "\n");
                }
            } else {
                String label = treeString.charAt(cursor) + "";
                cursor++;
                while (cursor < treeString.length() && treeString.charAt(cursor) != CHILD_SEPARATOR
                       && treeString.charAt(cursor) != CLOSE_BRACKET) {
                    label = label + treeString.charAt(cursor);
                    cursor++;
                }
                label = label.trim();
                
                List<String> labelList = new ArrayList<String>();
                labelList.add(label);
                //currentNode.setLabel(label);
                currentNode.setLabel(labelList);
            }

        }

        return tree;
    }

}
