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

package tgl;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Random;

import cycle.CyclicityTest;
import generator.IParasiteGenerator;
import generator.IScenario;
import generator.defaultmodel.DefaultModel;
import trees.Tree;
import trees.TreeNode;
import util.NexusFileParserException;

public class DistancesProcess {

    /* -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:- */
    /* Constants which define the default values for the input parameters of the process. */
    public static final double DEFAULT_COSPECIATION = 0.7;

    public static final double DEFAULT_DUPLICATION = 0.1;

    public static final double DEFAULT_HOSTSWITCH = 0.1;

    public static final int DEFAULT_NUMBER_OF_TREES = 10;

    public static final int DEFAULT_MAXIMUM_SIZE = 100;

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

    private double loss;

    private int numberOfTrees;

    private String hostFileName;

    private String parasiteFileName;

    private int cyclicityTest;

    private int maximumSize;

    private double rootMappingProbability;

    /** Expected number of co-speciations */
    private int expectedNumberOfCospeciationsHost;

    /** Expected number of duplications */
    private int expectedNumberOfDuplicationsHost;

    /** Expected number of host-switches */
    private int expectedNumberOfHostSwitchesHost;

    /** Expected number of losses */
    private int expectedNumberOfLossesHost;

    /** Expected number of co-speciations */
    private int expectedNumberOfCospeciationsParasite;

    /** Expected number of duplications */
    private int expectedNumberOfDuplicationsParasite;

    /** Expected number of host-switches */
    private int expectedNumberOfHostSwitchesParasite;

    /** Expected number of losses */
    private int expectedNumberOfLossesParasite;

    /* -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:- */
    /* Work attributes */
    private Tree hostTree;

    private Tree hostTreeWithRHO;

    private Tree parasiteTreeToCompare;

    public Random random;

    private NumberFormat format;

    private HashMap<String, Integer> hostTreeLeafSet;

    private int numberOfLeavesHostTree;

    private HashMap<String, Integer> parasiteTreeLeafSet;

    private int numberOfLeavesParasiteTree;


    /**
     * Constructor: Creates an TGLGeneratorProcess object filled with default parameters.
     */
    public DistancesProcess() {
        this.cospeciation = DEFAULT_COSPECIATION;
        this.duplication = DEFAULT_DUPLICATION;
        this.hostSwitch = DEFAULT_HOSTSWITCH;
        this.numberOfTrees = DEFAULT_NUMBER_OF_TREES;
        this.cyclicityTest = DEFAULT_CYCLICITY_TEST;
        this.rootMappingProbability = DEFAULT_ROOT_MAPPING_PROBABILITY;
        this.hostFileName = null;
        this.parasiteFileName = null;
        this.hostTree = null;
        this.hostTreeWithRHO = null;
        this.maximumSize = DEFAULT_MAXIMUM_SIZE;
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


    public void setHostFile(String file) {
        this.hostFileName = file;
    }


    public void setParasiteFile(String file) {
        this.parasiteFileName = file;
    }


    public void setNumberOfTrees(int number) {
        this.numberOfTrees = number;
    }


    public void setCyclicityTest(int cyclicityTest) {
        this.cyclicityTest = cyclicityTest;
    }


    public void setMaximumSize(int size) {
        this.maximumSize = size;
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
        readHostNewickFile(hostFileName);
        readParasiteNewickFile(parasiteFileName);

        double cds = cospeciation + duplication + hostSwitch;
        loss = 1.0 - cds;

        int nInternalHost = hostTree.getNumberOfLeafNodes() - 1;
        int nInternalParasite = parasiteTreeToCompare.getNumberOfLeafNodes() - 1;

        double pc = cospeciation / cds;
        double pd = duplication / cds;
        double ps = hostSwitch / cds;

        double nTotalHost = nInternalHost / cds;
        double nTotalParasite = nInternalParasite / cds;

        expectedNumberOfCospeciationsHost = (int)Math.round(nInternalHost * pc);
        expectedNumberOfDuplicationsHost = (int)Math.round(nInternalHost * pd);
        expectedNumberOfHostSwitchesHost = (int)Math.round(nInternalHost * ps);
        expectedNumberOfLossesHost = (int)Math.round(nTotalHost * loss);

        expectedNumberOfCospeciationsParasite = (int)Math.round(nInternalParasite * pc);
        expectedNumberOfDuplicationsParasite = (int)Math.round(nInternalParasite * pd);
        expectedNumberOfHostSwitchesParasite = (int)Math.round(nInternalParasite * ps);
        expectedNumberOfLossesParasite = (int)Math.round(nTotalParasite * loss);

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
    private void readHostNewickFile(String fileName) throws IOException {

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

        /* Create the LCA tables. */
        hostTree.createLCATable();
        hostTreeWithRHO.createLCATable();

        /* Update the list of ancestors. */
        hostTree.updateAncestorList();
        hostTreeWithRHO.updateAncestorList();

        /* Update the list of descendants. */
        hostTree.updateDescendantList();
        hostTreeWithRHO.updateDescendantList();

        /* Get the leaf set of the host tree. */
        hostTreeLeafSet = getLeafSet(hostTree);

        /* Save the number of leaves of the host tree. */
        numberOfLeavesHostTree = hostTree.getNumberOfLeafNodes();

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
    private void readParasiteNewickFile(String fileName) throws IOException {

        File file = new File(fileName);
        if (!file.exists()) {
            throw new IOException("File not found: " + fileName);
        }

        BufferedReader reader = new BufferedReader(new FileReader(file));
        try {
            String treeString = reader.readLine();
            parasiteTreeToCompare = parseNewickString(treeString, "!P");
        } finally {
            reader.close();
        }

        /* Create the LCA table. */
        parasiteTreeToCompare.createLCATable();

        /* Update the list of ancestors. */
        parasiteTreeToCompare.updateAncestorList();

        /* Update the list of descendants. */
        parasiteTreeToCompare.updateDescendantList();

        /* Get the leaf set of the host tree. */
        parasiteTreeLeafSet = getLeafSet(parasiteTreeToCompare);

        /* Save the number of leaves of the host tree. */
        numberOfLeavesParasiteTree = parasiteTreeToCompare.getNumberOfLeafNodes();

    }


    /**
     * Generate the desired number of parasite trees.
     * 
     * @throws IOException
     *         When some IO operation fails.
     */
    private void generateTrees() throws IOException {

        IParasiteGenerator generator = new DefaultModel(hostTreeWithRHO,
        												parasiteTreeToCompare,
                                                        cospeciation,
                                                        duplication,
                                                        hostSwitch,
                                                        -1,
                                                        maximumSize,
                                                        cyclicityTest,
                                                        true,
                                                        rootMappingProbability);

        int printed = 0;
        int numberOfIterations = 0;
        while (printed < numberOfTrees) {

            IScenario scenario = generator.generateParasiteTree();
            if (scenario != null) {
                double eDistance = scenario.computeEuclideanDistance(cospeciation,
                                                                     duplication,
                                                                     hostSwitch,
                                                                     loss);

                double c1H = computeComponent1Host(scenario);
                double[] c23H = computeComponents2And3Host(scenario);

                double c1P = computeComponent1Parasite(scenario);
                double[] c23P = computeComponents2And3Parasite(scenario);

                double[] d3H = new double[11];
                double[] d3P = new double[11];
                double[] d2H = new double[11];
                double[] d2P = new double[11];

                for (int i = 0; i < 11; i++) {
                    double w1 = 0.1 * i;
                    double w2 = 1.0 - w1;
                    double w23 = w2 / 2;
                    d3H[i] = compute3ComponentsDistance(c1H, c23H[0], c23H[1], w1, w23, w23);
                    d3P[i] = compute3ComponentsDistance(c1P, c23P[0], c23P[1], w1, w23, w23);
                    d2H[i] = compute2ComponentsDistance(c1H, c23H[1], w1, w2);
                    d2P[i] = compute2ComponentsDistance(c1P, c23P[1], w1, w2);
                }

                StringBuffer buffer = new StringBuffer();

                buffer.append(scenario.getNumberOfLeaves());
                buffer.append("\t");
                buffer.append(format.format(eDistance));
                buffer.append("\t");

                for (int i = 0; i < 11; i++) {
                    buffer.append(format.format(d3H[i]));
                    buffer.append("\t");
                }

                for (int i = 0; i < 11; i++) {
                    buffer.append(format.format(d3P[i]));
                    buffer.append("\t");
                }

                for (int i = 0; i < 11; i++) {
                    buffer.append(format.format(d2H[i]));
                    buffer.append("\t");
                }

                for (int i = 0; i < 11; i++) {
                    buffer.append(format.format(d2P[i]));
                    buffer.append("\t");
                }

                buffer.append(format.format(scenario.getFrequencyOfCospeciations()));
                buffer.append("\t");
                buffer.append(format.format(scenario.getFrequencyOfDuplications()));
                buffer.append("\t");
                buffer.append(format.format(scenario.getFrequencyOfHostSwitches()));
                buffer.append("\t");
                buffer.append(format.format(scenario.getFrequencyOfLosses()));

                // buffer.append("\t");
                // buffer.append(scenario.getParasiteTree().toString(true, false));

                System.out.println(buffer.toString());

                printed++;
            }

            numberOfIterations++;

            if (numberOfIterations > 10000 * numberOfTrees) {
                System.err.println("Aborting process: We reached the maximum number of interations!");
                break;
            }

        }

        System.err.println("Discarded = "
                           + format.format((numberOfIterations - numberOfTrees)
                                           / (double)numberOfIterations));
    }


    private double[] computeComponents2And3Host(IScenario scenario) {

        Tree tree = scenario.getParasiteTree();
        double maacValue = tree.computeMaac(hostTree).size();
        HashMap<String, Integer> treeLeafSet = getLeafSet(tree);
        double intersection = getIntersectionSize(treeLeafSet, hostTreeLeafSet);
        double numberOfLeaves = tree.getNumberOfLeafNodes();

        if (intersection > 0) {
            double[] toReturn = {
                                 (1.0 - (2.0 * intersection / (numberOfLeaves + numberOfLeavesHostTree))),
                                 (1.0 - (maacValue / intersection))};
            return toReturn;
        }

        double[] toReturn = {1.0, 1.0};
        return toReturn;
    }


    private double[] computeComponents2And3Parasite(IScenario scenario) {
        Tree tree = scenario.getParasiteTree();
        double maacValue = tree.computeMaac(parasiteTreeToCompare).size();
        HashMap<String, Integer> treeLeafSet = getLeafSet(tree);
        double intersection = getIntersectionSize(treeLeafSet, parasiteTreeLeafSet);
        double numberOfLeaves = tree.getNumberOfLeafNodes();
        if (intersection > 0) {
            double[] toReturn = {
                                 (1.0 - (2.0 * intersection / (numberOfLeaves + numberOfLeavesParasiteTree))),
                                 (1.0 - (maacValue / intersection))};
            return toReturn;
        }
        double[] toReturn = {1.0, 1.0};
        return toReturn;
    }


    private double computeComponent1Host(IScenario scenario) {

        double diffC = Math.abs(expectedNumberOfCospeciationsHost
                                - scenario.getNumberOfCospeciations());
        double diffD = Math.abs(expectedNumberOfDuplicationsHost
                                - scenario.getNumberOfDuplications());
        double diffS = Math.abs(expectedNumberOfHostSwitchesHost
                                - scenario.getNumberOfHostSwitches());
        double diffL = Math.abs(expectedNumberOfLossesHost - scenario.getNumberOfLosses());

        int maxC = expectedNumberOfCospeciationsHost;
        if (maxC < scenario.getNumberOfCospeciations())
            maxC = scenario.getNumberOfCospeciations();

        int maxD = expectedNumberOfDuplicationsHost;
        if (maxD < scenario.getNumberOfDuplications())
            maxD = scenario.getNumberOfDuplications();

        int maxS = expectedNumberOfHostSwitchesHost;
        if (maxS < scenario.getNumberOfHostSwitches())
            maxS = scenario.getNumberOfHostSwitches();

        int maxL = expectedNumberOfLossesHost;
        if (maxL < scenario.getNumberOfLosses())
            maxL = scenario.getNumberOfLosses();

        double partC = 0.0;
        double partD = 0.0;
        double partS = 0.0;
        double partL = 0.0;

        if (maxC > 0) {
            partC = diffC / maxC;
        }

        if (maxD > 0) {
            partD = diffD / maxD;
        }

        if (maxS > 0) {
            partS = diffS / maxS;
        }

        if (maxL > 0) {
            partL = diffL / maxL;
        }

        return (partC + partD + partS + partL) / 4.0;
    }


    private double compute3ComponentsDistance(double c1,
                                              double c2,
                                              double c3,
                                              double w1,
                                              double w2,
                                              double w3) {
        return (c1 * w1) + (c2 * w2) + (c3 * w3);
    }


    private double compute2ComponentsDistance(double c1, double c2, double w1, double w2) {
        return (c1 * w1) + (c2 * w2);
    }


    private double computeComponent1Parasite(IScenario scenario) {

        double diffC = Math.abs(expectedNumberOfCospeciationsParasite
                                - scenario.getNumberOfCospeciations());
        double diffD = Math.abs(expectedNumberOfDuplicationsParasite
                                - scenario.getNumberOfDuplications());
        double diffS = Math.abs(expectedNumberOfHostSwitchesParasite
                                - scenario.getNumberOfHostSwitches());
        double diffL = Math.abs(expectedNumberOfLossesParasite - scenario.getNumberOfLosses());

        int maxC = expectedNumberOfCospeciationsParasite;
        if (maxC < scenario.getNumberOfCospeciations())
            maxC = scenario.getNumberOfCospeciations();

        int maxD = expectedNumberOfDuplicationsParasite;
        if (maxD < scenario.getNumberOfDuplications())
            maxD = scenario.getNumberOfDuplications();

        int maxS = expectedNumberOfHostSwitchesParasite;
        if (maxS < scenario.getNumberOfHostSwitches())
            maxS = scenario.getNumberOfHostSwitches();

        int maxL = expectedNumberOfLossesParasite;
        if (maxL < scenario.getNumberOfLosses())
            maxL = scenario.getNumberOfLosses();

        double partC = 0.0;
        double partD = 0.0;
        double partS = 0.0;
        double partL = 0.0;

        if (maxC > 0) {
            partC = diffC / maxC;
        }

        if (maxD > 0) {
            partD = diffD / maxD;
        }

        if (maxS > 0) {
            partS = diffS / maxS;
        }

        if (maxL > 0) {
            partL = diffL / maxL;
        }

        return (partC + partD + partS + partL) / 4.0;

    }


    // /**
    // * Return the modified MAAC distance between the given tree and the base host tree.
    // *
    // * @param scenario
    // * Generated scenario.
    // * @return The modified MAAC distance between the given tree and the base host tree.
    // */
    // private double computeMMAACDistanceHost(IScenario scenario) {
    // Tree tree = scenario.getParasiteTree();
    // double maacValue = tree.computeMaac(hostTree).size();
    // HashMap<String, Integer> treeLeafSet = getLeafSet(tree);
    // double intersection = getIntersectionSize(treeLeafSet, hostTreeLeafSet);
    // double numberOfLeaves = tree.getNumberOfLeafNodes();
    // if (intersection > 0) {
    // return ALPHA1
    // * (1.0 - (2.0 * intersection / (numberOfLeaves + numberOfLeavesHostTree)))
    // + ALPHA2 * (1.0 - (maacValue / intersection));
    // }
    // return 1.0;
    // }

    // /**
    // * Return the modified MAAC distance between the given tree and the base host tree.
    // *
    // * @param scenario
    // * Generated scenario.
    // * @return The modified MAAC distance between the given tree and the base host tree.
    // */
    // private double computeMMAACDistanceParasite(IScenario scenario) {
    // Tree tree = scenario.getParasiteTree();
    // double maacValue = tree.computeMaac(parasiteTreeToCompare).size();
    // HashMap<String, Integer> treeLeafSet = getLeafSet(tree);
    // double intersection = getIntersectionSize(treeLeafSet, parasiteTreeLeafSet);
    // double numberOfLeaves = tree.getNumberOfLeafNodes();
    // if (intersection > 0) {
    // return ALPHA1
    // * (1.0 - (2.0 * intersection / (numberOfLeaves + numberOfLeavesParasiteTree)))
    // + ALPHA2 * (1.0 - (maacValue / intersection));
    // }
    // return 1.0;
    // }

    // private double computeDistanceHost(IScenario scenario) {
    //
    // double diffC = Math.abs(expectedNumberOfCospeciationsHost
    // - scenario.getNumberOfCospeciations());
    // double diffD = Math.abs(expectedNumberOfDuplicationsHost
    // - scenario.getNumberOfDuplications());
    // double diffS = Math.abs(expectedNumberOfHostSwitchesHost
    // - scenario.getNumberOfHostSwitches());
    // double diffL = Math.abs(expectedNumberOfLossesHost - scenario.getNumberOfLosses());
    //
    // int maxC = expectedNumberOfCospeciationsHost;
    // if (maxC < scenario.getNumberOfCospeciations())
    // maxC = scenario.getNumberOfCospeciations();
    //
    // int maxD = expectedNumberOfDuplicationsHost;
    // if (maxD < scenario.getNumberOfDuplications())
    // maxD = scenario.getNumberOfDuplications();
    //
    // int maxS = expectedNumberOfHostSwitchesHost;
    // if (maxS < scenario.getNumberOfHostSwitches())
    // maxS = scenario.getNumberOfHostSwitches();
    //
    // int maxL = expectedNumberOfLossesHost;
    // if (maxL < scenario.getNumberOfLosses())
    // maxL = scenario.getNumberOfLosses();
    //
    // double partC = 0.0;
    // double partD = 0.0;
    // double partS = 0.0;
    // double partL = 0.0;
    //
    // if (maxC > 0) {
    // partC = diffC / maxC;
    // }
    //
    // if (maxD > 0) {
    // partD = diffD / maxD;
    // }
    //
    // if (maxS > 0) {
    // partS = diffS / maxS;
    // }
    //
    // if (maxL > 0) {
    // partL = diffL / maxL;
    // }
    //
    // double component1 = (partC + partD + partS + partL) / 4.0;
    //
    // Tree tree = scenario.getParasiteTree();
    // double maacValue = tree.computeMaac(hostTree).size();
    // HashMap<String, Integer> treeLeafSet = getLeafSet(tree);
    // double intersection = getIntersectionSize(treeLeafSet, hostTreeLeafSet);
    // double numberOfLeaves = tree.getNumberOfLeafNodes();
    // if (intersection > 0) {
    // return (BETA1 * component1)
    // + (BETA2 * (1.0 - (2.0 * intersection / (numberOfLeaves + numberOfLeavesHostTree))))
    // + (BETA3 * (1.0 - (maacValue / intersection)));
    // }
    // return 1.0;
    // }
    //
    //
    // private double computeDistanceParasite(IScenario scenario) {
    //
    // double diffC = Math.abs(expectedNumberOfCospeciationsParasite
    // - scenario.getNumberOfCospeciations());
    // double diffD = Math.abs(expectedNumberOfDuplicationsParasite
    // - scenario.getNumberOfDuplications());
    // double diffS = Math.abs(expectedNumberOfHostSwitchesParasite
    // - scenario.getNumberOfHostSwitches());
    // double diffL = Math.abs(expectedNumberOfLossesParasite - scenario.getNumberOfLosses());
    //
    // int maxC = expectedNumberOfCospeciationsParasite;
    // if (maxC < scenario.getNumberOfCospeciations())
    // maxC = scenario.getNumberOfCospeciations();
    //
    // int maxD = expectedNumberOfDuplicationsParasite;
    // if (maxD < scenario.getNumberOfDuplications())
    // maxD = scenario.getNumberOfDuplications();
    //
    // int maxS = expectedNumberOfHostSwitchesParasite;
    // if (maxS < scenario.getNumberOfHostSwitches())
    // maxS = scenario.getNumberOfHostSwitches();
    //
    // int maxL = expectedNumberOfLossesParasite;
    // if (maxL < scenario.getNumberOfLosses())
    // maxL = scenario.getNumberOfLosses();
    //
    // double partC = 0.0;
    // double partD = 0.0;
    // double partS = 0.0;
    // double partL = 0.0;
    //
    // if (maxC > 0) {
    // partC = diffC / maxC;
    // }
    //
    // if (maxD > 0) {
    // partD = diffD / maxD;
    // }
    //
    // if (maxS > 0) {
    // partS = diffS / maxS;
    // }
    //
    // if (maxL > 0) {
    // partL = diffL / maxL;
    // }
    //
    // double component1 = (partC + partD + partS + partL) / 4.0;
    //
    // Tree tree = scenario.getParasiteTree();
    // double maacValue = tree.computeMaac(parasiteTreeToCompare).size();
    // HashMap<String, Integer> treeLeafSet = getLeafSet(tree);
    // double intersection = getIntersectionSize(treeLeafSet, parasiteTreeLeafSet);
    // double numberOfLeaves = tree.getNumberOfLeafNodes();
    // if (intersection > 0) {
    // return (BETA1 * component1)
    // + (BETA2 * (1.0 - (2.0 * intersection / (numberOfLeaves + numberOfLeavesParasiteTree))))
    // + (BETA3 * (1.0 - (maacValue / intersection)));
    // }
    // return 1.0;
    // }

    /**
     * Compute the size of the intersection of the leaf sets of both trees.
     * 
     * @param tree1
     *        Leaf set of tree 1.
     * @param tree2
     *        Leaf set of tree 2.
     * @return The size of the intersection of the leaf sets of both trees.
     */
    private int getIntersectionSize(HashMap<String, Integer> tree1, HashMap<String, Integer> tree2) {
        int sum = 0;
        for (String label: tree1.keySet()) {
            if (tree2.containsKey(label)) {
                int n1 = tree1.get(label);
                int n2 = tree2.get(label);
                if (n1 <= n2) {
                    sum += n1;
                } else {
                    sum += n2;
                }
            }
        }
        return sum;
    }


    /**
     * Create a Hashmap with the leaf set of the given tree.
     * 
     * @param tree
     *        Tree to compute the leaf set.
     * @return A HashMap where the key is the label and the value is the number of occurrences of
     *         the label in the tree.
     */
    private HashMap<String, Integer> getLeafSet(Tree tree) {
        HashMap<String, Integer> toReturn = new HashMap<String, Integer>();
        for (TreeNode node: tree.getLeafNodes()) {
            Integer n = toReturn.get(node.getLabel());
            if (n == null) {
                n = 1;
            } else {
                n++;
            }
            //toReturn.put(node.getLabel(), n);
            toReturn.put(node.getLabel().toString(), n);
        }
        return toReturn;
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
        List<String> setLabelListRoot = new ArrayList<String>();
        setLabelListRoot.add(labelPrefix + Integer.toString(root.getKey()));
        //root.setLabel(labelPrefix + Integer.toString(root.getKey()));
        root.setLabel(setLabelListRoot);

        int brackets = 0;

        TreeNode currentNode = root;
        while (cursor < treeString.length()) {
            if (treeString.charAt(cursor) == OPEN_BRACKET) {
                TreeNode child = tree.createNode();
                
                List<String> setLabelList = new ArrayList<String>();
                setLabelList.add(labelPrefix + Integer.toString(child.getKey()));
                //child.setLabel(labelPrefix + Integer.toString(child.getKey()));
                child.setLabel(setLabelList);
                
                tree.addChild(currentNode.getKey(), child.getKey(), 0);
                currentNode = child;
                cursor++;
                brackets++;
            } else if (treeString.charAt(cursor) == CHILD_SEPARATOR) {
                currentNode = currentNode.getParent();
                TreeNode child = tree.createNode();
                
                List<String> setLabelList = new ArrayList<String>();
                setLabelList.add(labelPrefix + Integer.toString(child.getKey()));
                //child.setLabel(labelPrefix + Integer.toString(child.getKey()));
                child.setLabel(setLabelList);
                
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
