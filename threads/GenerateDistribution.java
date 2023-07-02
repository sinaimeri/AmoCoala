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

package  threads;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import generator.IScenario;
import  trees.Tree;
import  trees.TreeNode;

public abstract class GenerateDistribution extends Thread {

    /**
     * Constant which indicates that the simulation model will use the a maximum size limit to abort
     * the generation of a tree.
     */
    public static final boolean MAXIMUM_SIZE_STOP_CRITERIUM = true;

    /**
     * Constant that identifies the default generation model.
     */
    public static final int DEFAULT_MODEL = 1;

    /**
     * Constant that identifies the coalescent model which does not use the information of distance
     * between host nodes.
     */
    public static final int COALESCENCE_WITHOUT_DISTANCE = 2;

    /**
     * Constant that identifies the coalescent model which does not use the information of distance
     * between host nodes.
     */
    public static final int COALESCENCE_WITH_DISTANCE = 3;

    /**
     * Constant that defines the default factor to compute the maximum number of trees that can be
     * simulated for a vector.
     */
    public static final int DEFAULT_MAXIMUM_NUMBER_OF_TREES_FACTOR = 5;

    /** Defines the maximum size factor for generating the simulated trees. */
    public static final double MAXIMUM_SIZE_FACTOR = 2.0;

    /** Alpha 1 defines the weight of the component 1 of the LEAVES AND MAAC metric. */
    public static final double DEFAULT_LEAVES_AND_MAAC_ALPHA1 = 0.5;

    /** Alpha 2 defines the weight of the component 2 of the LEAVES AND MAAC metric. */
    public static final double DEFAULT_LEAVES_AND_MAAC_ALPHA2 = 1.0 - DEFAULT_LEAVES_AND_MAAC_ALPHA1;

    /** Alpha 1 defines the weight of the component 1 of the EVENTS AND MAAC metric. */
    public static final double DEFAULT_EVENTS_AND_MAAC_ALPHA1 = 0.7;

    /** Alpha 2 defines the weight of the component 2 of the EVENTS AND MAAC metric. */
    public static final double DEFAULT_EVENTS_AND_MAAC_ALPHA2 = 1.0 - DEFAULT_EVENTS_AND_MAAC_ALPHA1;

    /**
     * Constant that identifies the MAAC metric.
     */
    public static final int METRIC_MAAC_DISTANCE = 1;

    /**
     * Constant that identifies the LEAVES AND MAAC metric. This distance has two components: the
     * first component deals with the leaf set of the trees and the second (MAAC) deals with the
     * topology.
     */
    public static final int METRIC_LEAVES_AND_MAAC_DISTANCE = 2;

    /**
     * Constant that identifies the EVENTS AND MAAC metric. This distance has two components: the
     * first component deals with the number of events and the second (MAAC) deals with the
     * topology.
     */
    public static final int METRIC_EVENTS_AND_MAAC_DISTANCE = 3;
    
    /**
     * Constant that identifies the EVENTS AND MAAC metric. This distance has two components: the
     * first component deals with the number of events and the second (MAAC) deals with the
     * topology.
     */
    public static final int METRIC_EVENTS_AND_MAAC_DISTANCE_MULTIPLEASSOCIATIONS = 4;

    /* ------------------------------------------------------------------------------------------ */

    /** Synchronized object which handles the list of generated vectors. */
    protected SynchronizedVectorManager vectorManager;

    /** Model which will be used during the simulation. */
    protected int model;

    /** Number of vectors to generate. */
    protected int numberOfVectors;

    /** Number of trees to generate for each vector. */
    protected int numberOfTrees;

    /** Maximum number of leaves for the generated trees. */
    protected int maximumNumberOfLeaves;

    /** Maximum number of trees that can be generated for a given vector. */
    protected int maximumNumberOfTrees;

    /** Metric that is going to be used to compare the trees. */
    protected int metric;

    /** Cyclicity test model */
    protected int cyclicityTestModel;

    /** Host tree. */
    protected Tree hostTree;

    /** Real parasite tree. */
    protected Tree parasiteTree;

    /** Mapping of the parasite leaf nodes into the host leaf nodes. */
    //protected HashMap<String, String> mappingParasiteHost;
    protected HashMap<String, List<String>> mappingParasiteHost;

    /** Number of leaves of the real parasite tree. */
    protected int numberOfLeavesParasiteTree;

    /** Real parasite leaf set. HashMap -> <label, number of occurences>. */
    protected HashMap<List<String>, Integer> parasiteLeafSet;

    /** Alpha 1 for distance metric : LEAVES AND MAAC or EVENTS AND MAAC */
    protected double alpha1;

    /** Alpha 2 for distance metric : LEAVES AND MAAC or EVENTS AND MAAC */
    protected double alpha2;

    protected double rootMappingProbability;

    /** Expected number of co-speciations */
    private int expectedNumberOfCospeciations;

    /** Expected number of duplications */
    private int expectedNumberOfDuplications;

    /** Expected number of host-switches */
    private int expectedNumberOfHostSwitches;

    /** Expected number of losses */
    private int expectedNumberOfLosses;
    
    /** Expected number of multiple associations */
    private int expectedNumberMultipleAssociations;
    
    protected static int nbFile = 1; 


    /**
     * Thread to generate vectors.
     * 
     * @param vectorManager
     *        Object that manage the access to the list of vectors in a synchronized way.
     * @param numberOfVectors
     *        Number of vectors to generate.
     * @param numberOfTrees
     *        Number of trees that must be generated for each vector.
     * @param hostTree
     *        Host tree.
     * @param parasiteTree
     *        Parasite tree.
     * @param mappingParasiteHost2
     *        Mapping of the parasite leaf nodes into the host leaf nodes.
     * @param model
     *        Model to be used during the simulation.
     * @param metric
     *        Metric to be used during the simulation.
     * @param maximumNumberOfTreesFactor
     *        Maximum number of trees factor.
     * @param cyclicityTestModel
     *        Cyclicity test model identifier.
     * @param alpha1
     *        Alpha 1 value for computing LEAVES AND MAAC or EVENTS AND MAAC metric.
     * @param alpha2
     *        Alpha 2 value for computing LEAVES AND MAAC or EVENTS AND MAAC metric.
     * @param rootMappingProbability
     *        Probability of mapping into the root node.
     */
    public GenerateDistribution(SynchronizedVectorManager vectorManager,
                                int numberOfVectors,
                                int numberOfTrees,
                                Tree hostTree,
                                Tree parasiteTree,
                                //HashMap<String, String> mappingParasiteHost,
                                HashMap<String, List<String>> mappingParasiteHost,
                                int model,
                                int metric,
                                int maximumNumberOfTreesFactor,
                                int cyclicityTestModel,
                                double alpha1,
                                double alpha2,
                                double rootMappingProbability) {
        super();
        
        if (model != DEFAULT_MODEL && model != COALESCENCE_WITHOUT_DISTANCE
            && model != COALESCENCE_WITH_DISTANCE) {
            throw new RuntimeException("Invalid model.");
        }

        if (metric != METRIC_MAAC_DISTANCE && metric != METRIC_LEAVES_AND_MAAC_DISTANCE
            && metric != METRIC_EVENTS_AND_MAAC_DISTANCE && metric != METRIC_EVENTS_AND_MAAC_DISTANCE_MULTIPLEASSOCIATIONS) {
            throw new RuntimeException("Invalid distance.");
        }

        this.vectorManager = vectorManager;
        this.numberOfVectors = numberOfVectors;
        this.numberOfTrees = numberOfTrees;

        this.hostTree = hostTree;

        this.parasiteTree = parasiteTree;
        this.mappingParasiteHost = mappingParasiteHost;

        this.model = model;
        this.metric = metric;

        this.cyclicityTestModel = cyclicityTestModel;
        this.maximumNumberOfTrees = numberOfTrees * maximumNumberOfTreesFactor;

        this.parasiteLeafSet = getLeafSet(parasiteTree);

        this.numberOfLeavesParasiteTree = parasiteTree.getLeafNodes().length;

        this.maximumNumberOfLeaves = (int)Math.round(numberOfLeavesParasiteTree
                                                     * MAXIMUM_SIZE_FACTOR);

        this.alpha1 = alpha1;
        this.alpha2 = alpha2;

        this.rootMappingProbability = rootMappingProbability;
       
    }


    /**
     * Return the MAAC Distance between the given tree and the real parasite tree.
     * 
     * @param scenario
     *        Generated scenario.
     * @return The MAAC Distance between the given tree and the real parasite tree.
     */
    protected double computeMAACDistanceFromRealParasite(IScenario scenario) {
        return scenario.getParasiteTree().computeMaacDistance(parasiteTree);
    }


    /**
     * Return the modified LEAVES AND MAAC metric between the given tree and the real parasite tree.
     * 
     * @param scenario
     *        Generated scenario.
     * @return The LEAVES AND MAAC metric between the given tree and the real parasite tree.
     */
    protected double computeLEAVES_AND_MAAC_Metric(IScenario scenario) {
        double maacValue = scenario.getParasiteTree().computeMaac(parasiteTree).size();
        HashMap<List<String>, Integer> treeLeafSet = getLeafSet(scenario.getParasiteTree());
        double intersection = getIntersectionSize(parasiteLeafSet, treeLeafSet);
        double numberOfLeaves = scenario.getNumberOfLeaves();
        if (intersection > 0) {
            return alpha1
                   * (1.0 - (2.0 * intersection / (numberOfLeaves + numberOfLeavesParasiteTree)))
                   + alpha2 * (1.0 - (maacValue / intersection));
        }
        return 1.0;
    }


    /**
     * Return the modified EVENTS AND MAAC metric between the given tree and the real parasite tree.
     * 
     * @param scenario
     *        Generated scenario.
     * @return The EVENTS AND MAAC metric between the given tree and the real parasite tree.
     */
    protected double computeEVENTS_AND_MAAC_Metric(IScenario scenario) {

        double diffC = Math.abs(expectedNumberOfCospeciations - scenario.getNumberOfCospeciations());
        double diffD = Math.abs(expectedNumberOfDuplications - scenario.getNumberOfDuplications());
        double diffS = Math.abs(expectedNumberOfHostSwitches - scenario.getNumberOfHostSwitches());
        double diffL = Math.abs(expectedNumberOfLosses - scenario.getNumberOfLosses());

        int maxC = expectedNumberOfCospeciations;
        if (maxC < scenario.getNumberOfCospeciations())
            maxC = scenario.getNumberOfCospeciations();

        int maxD = expectedNumberOfDuplications;
        if (maxD < scenario.getNumberOfDuplications())
            maxD = scenario.getNumberOfDuplications();

        int maxS = expectedNumberOfHostSwitches;
        if (maxS < scenario.getNumberOfHostSwitches())
            maxS = scenario.getNumberOfHostSwitches();

        int maxL = expectedNumberOfLosses;
        if (maxL < scenario.getNumberOfLosses())
            maxL = scenario.getNumberOfLosses();

        double partC = 0.0;
        if (maxC > 0)
            partC = diffC / maxC;

        double partD = 0.0;
        if (maxD > 0)
            partD = diffD / maxD;

        double partS = 0.0;
        if (maxS > 0)
            partS = diffS / maxS;

        double partL = 0.0;
        if (maxL > 0)
            partL = diffL / maxL;

        double component1 = (partC + partD + partS + partL) / 4.0;

        double maacValue = scenario.getParasiteTree().computeMaac(parasiteTree).size();
        HashMap<List<String>, Integer> treeLeafSet = getLeafSet(scenario.getParasiteTree());
        double intersection = getIntersectionSize(parasiteLeafSet, treeLeafSet);
     //   double union = getUnionLabels(parasiteLeafSet, treeLeafSet);
        if (intersection > 0) {
            return (alpha1 * component1) + (alpha2 * (1.0 - (maacValue / intersection)));
        }
        
        //if (union > 0) {
        	//return (alpha1 * component1) + (alpha2 * (1.0 - (maacValue / union)));
       // }		
        		

        return 1.0;
    }


    /**
     * Return the modified EVENTS AND MAAC metric between the given tree and the real parasite tree in the case of multiple associations.
     * 
     * @param scenario
     *        Generated scenario.
     * @return The EVENTS AND MAAC metric between the given tree and the real parasite tree in the case of multiple associations.
     */
    protected double computeEVENTS_AND_MAAC_Metric_MultipleAssociations(IScenario scenario, int expectedNumberMultipleAssociations) {
    	this.expectedNumberMultipleAssociations = expectedNumberMultipleAssociations;
    	
        double diffC = Math.abs(expectedNumberOfCospeciations - scenario.getNumberOfCospeciations());
        double diffD = Math.abs(expectedNumberOfDuplications - scenario.getNumberOfDuplications());
        double diffS = Math.abs(expectedNumberOfHostSwitches - scenario.getNumberOfHostSwitches());
        double diffL = Math.abs(expectedNumberOfLosses - scenario.getNumberOfLosses());
        double diffMA = Math.abs(expectedNumberMultipleAssociations - scenario.getParasiteTree().getNumberOfMultipleAssociations());

        int maxC = expectedNumberOfCospeciations;
        if (maxC < scenario.getNumberOfCospeciations())
            maxC = scenario.getNumberOfCospeciations();

        int maxD = expectedNumberOfDuplications;
        if (maxD < scenario.getNumberOfDuplications())
            maxD = scenario.getNumberOfDuplications();

        int maxS = expectedNumberOfHostSwitches;
        if (maxS < scenario.getNumberOfHostSwitches())
            maxS = scenario.getNumberOfHostSwitches();

        int maxL = expectedNumberOfLosses;
        if (maxL < scenario.getNumberOfLosses())
            maxL = scenario.getNumberOfLosses();
        
        int maxMA = expectedNumberMultipleAssociations;
        if (maxMA < scenario.getParasiteTree().getNumberOfMultipleAssociations())
            maxMA = scenario.getParasiteTree().getNumberOfMultipleAssociations();

        double partC = 0.0;
        if (maxC > 0)
            partC = diffC / maxC;

        double partD = 0.0;
        if (maxD > 0)
            partD = diffD / maxD;

        double partS = 0.0;
        if (maxS > 0)
            partS = diffS / maxS;

        double partL = 0.0;
        if (maxL > 0)
            partL = diffL / maxL;
        
        double partMA = 0.0;
        if (maxMA > 0)
        	partMA = diffMA / maxMA;

        double component1 = (partC + partD + partS + partL + partMA) / 5.0;
        //double component1 = (partC + partD + partS + partL) / 4.0;

        double maacValue = scenario.getParasiteTree().computeMaac(parasiteTree).size();
        TreeNode[] leaaafff = scenario.getParasiteTree().getLeafNodes();
        
        HashMap<List<String>, Integer> treeLeafSet = getLeafSet(scenario.getParasiteTree());
        //double union = getUnionLabels(parasiteLeafSet, treeLeafSet);
        double minLabels = getMinLabelBetweenTrees(parasiteLeafSet, treeLeafSet);
                
        //if (union > 0) {
        if (minLabels > 0) {
        	//return (alpha1 * component1) + (alpha2 * (1.0 - (maacValue / minLabels)));
        	
        }		
        		
        return (alpha1 * component1);
        //return 1.0;
    }

    
    /**
     * Create a Hashmap with the leaf set of the given tree.
     * 
     * @param tree
     *        Tree to compute the leaf set.
     * @return A HashMap where the key is the label and the value is the number of occurrences of
     *         the label in the tree.
     */
    private HashMap<List<String>, Integer> getLeafSet(Tree tree) {
    	//HashMap<String, Integer> toReturn = new HashMap<String, Integer>();
        HashMap<List<String>, Integer> toReturn = new HashMap<List<String>, Integer>();
        for (TreeNode node: tree.getLeafNodes()) {
            Integer n = toReturn.get(node.getLabel());
            if (n == null) {
                n = 1;
            } else {
                n++;
            }
            
            toReturn.put(node.getLabel(), n);
            
        }
        return toReturn;
    }


    /**
     * Compute the size of the intersection of the leaf sets of both trees.
     * 
     * @param tree1
     *        Leaf set of tree 1.
     * @param tree2
     *        Leaf set of tree 2.
     * @return The size of the intersection of the leaf sets of both trees.
     */
    private int getIntersectionSize(HashMap<List<String>, Integer> tree1, HashMap<List<String>, Integer> tree2) {
        int sum = 0;
        //for (String label: tree1.keySet()) {
        for (List<String> label: tree1.keySet()) {
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
     * Compute the size of the union of the labels leaf sets of both trees.
     * 
     * @param tree1
     *        Leaf set of tree 1.
     * @param tree2
     *        Leaf set of tree 2.
     * @return The size of the union of the labels leaf sets of both trees.
     */
    private int getUnionLabels(HashMap<List<String>, Integer> tree1, HashMap<List<String>, Integer> tree2) {
        int sum = 0;
        for (List<String> label: tree1.keySet()) {
            sum = sum + label.size();
        }
        
        for (List<String> label: tree2.keySet()) {
            sum = sum + label.size();
        }
        return sum;
    }
    
    
    /**
     * Get the size of the labels of tree with minimum number of labels sets of both trees.
     * 
     * @param tree1
     *        Leaf set of tree 1.
     * @param tree2
     *        Leaf set of tree 2.
     * @return The size of the labels of tree with minimum number of labels sets of both trees.
     */
    private int getMinLabelBetweenTrees(HashMap<List<String>, Integer> tree1, HashMap<List<String>, Integer> tree2){
    	
    	int sum1 = 0;
        for (List<String> label: tree1.keySet()) {
            sum1 = sum1 + label.size();
        }
        
        int sum2 = 0;
        for (List<String> label: tree2.keySet()) {
            sum2 = sum2 + label.size();
        }
        
        if(sum1 <= sum2)
        	return sum1; 
        
        return sum2;
    }


    /**
     * Update the expected number of events.
     * 
     * @param cospeciationProbability
     *        Co-speciation probability.
     * @param duplicationProbability
     *        Duplication probability.
     * @param hostSwitchProbability
     *        Host-switch probability.
     * @param lossProbability
     *        Loss probability.
     */
    protected void updateExpectedNumberOfEvents(double cospeciationProbability,
                                                double duplicationProbability,
                                                double hostSwitchProbability,
                                                double lossProbability) {
        int nInternal = parasiteTree.getNumberOfLeafNodes() - 1;
        double cds = cospeciationProbability + duplicationProbability + hostSwitchProbability;
        double pc = cospeciationProbability / cds;
        double pd = duplicationProbability / cds;
        double ps = hostSwitchProbability / cds;
        double nTotal = nInternal / cds;
        expectedNumberOfCospeciations = (int)Math.round(nInternal * pc);
        expectedNumberOfDuplications = (int)Math.round(nInternal * pd);
        expectedNumberOfHostSwitches = (int)Math.round(nInternal * ps);
        expectedNumberOfLosses = (int)Math.round(nTotal * lossProbability);
    }

}
