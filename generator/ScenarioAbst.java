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

package  generator;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import  trees.Tree;
import  trees.TreeNode;

public abstract class ScenarioAbst implements IScenario {

    /** Number of leaves. */
    private int numberOfLeaves;

    /** Fake parasite tree. */
    private Tree parasiteTree;

    /** Number of cospeciations. */
    private int cospeciations;

    /** Number of duplications. */
    private int duplications;

    /** Number of host switches. */
    private int hostSwitches;

    /** Number of losses. */
    private int losses;
    
    /** Number of vertical spread. */
    private int verticalSpread;
    
    /** Number of hertical spread. */
    private int horizontalSpread;
    
    /** Number of Multiple Associations. */
    private int multipleAssociations;

    /** List of associations. */
    private ArrayList<Association> associations;

    /** Hashmap which keeps the mapping parasite node to host node (key values). */
    protected HashMap<Integer, List<Integer>> parasiteToHost;

    /**
     * Maximum jump distance.
     */
    private double maxDistanceJump;
    
    /**
     * Add jump distance.
     */
    private double sumDistanceJump;

    /**
     * Store alle jumps of a tree.
     */
    private ArrayList<Double> jumpDistanceTot;
    
    /**
     * Add horizontal spread distance.
     */
    private double sumDistanceHorizontalSpread;
    
    
    /**
     * Creates an empty scenario: no associations and a parasite tree with only one node (root
     * node).
     */
    protected ScenarioAbst() {
        this.associations = new ArrayList<Association>();
        this.numberOfLeaves = 0;
        this.cospeciations = 0;
        this.duplications = 0;
        this.hostSwitches = 0;
        this.losses = 0;
        this.verticalSpread = 0;
        this.horizontalSpread = 0;
        this.parasiteToHost = new HashMap<Integer, List<Integer>>();
        this.maxDistanceJump = 0.00;
        this.sumDistanceJump = 0.00;
        this.sumDistanceHorizontalSpread = 0.00;
        jumpDistanceTot = new ArrayList<Double>();
    }


    /**
     * Performs the initialization of the parasite tree object.
     * 
     * @param createRoot
     *        If true, the tree is created with an initial root node. Otherwise, the tree will be
     *        empty.
     */
    protected void initParasiteTree(boolean createRoot) {
        if (parasiteTree == null) {
            if (createRoot) {
                numberOfLeaves = 1;
            }
            this.parasiteTree = new Tree(createRoot);
        } else {
            throw new RuntimeException("The tree object was already created before.");
        }
    }


    /**
     * Returns the parasite tree root node.
     * 
     * @return The parasite tree root node.
     */
    public TreeNode getRootNode() {
        return parasiteTree.getRoot();
    }


    /**
     * Returns the number of leaves.
     * 
     * @return The number of leaves.
     */
    public int getNumberOfLeaves() {
        return numberOfLeaves;
    }


    /**
     * Returns the number of cospeciations.
     * 
     * @return The number of cospeciations.
     */
    public int getNumberOfCospeciations() {
        return cospeciations;
    }


    /**
     * Returns the number of duplications.
     * 
     * @return The number of duplications.
     */
    public int getNumberOfDuplications() {
        return duplications;
    }


    /**
     * Returns the number of host switches.
     * 
     * @return The number of host switches.
     */
    public int getNumberOfHostSwitches() {
        return hostSwitches;
    }


    /**
     * Returns the number of losses.
     * 
     * @return The number of losses.
     */
    public int getNumberOfLosses() {
        return losses;
    }
    
    
    /**
     * Returns the number of multiple associations.
     * 
     * @return The number of multiple associations.
     */
    public int getNumberOfMultipleAssociations() {
        return multipleAssociations;
    }

    
    /**
     * Returns the number of vertical spread.
     * 
     * @return The number of vertical spread.
     */
    public int getNumberOfVerticalSpread() {
        return verticalSpread;
    }
    
    
    /**
     * Returns the number of horizontal spread.
     * 
     * @return The number of horizontal spread.
     */
    public int getNumberOfHorizontalSpread() {
        return horizontalSpread;
    }
    

    /**
     * Returns the generated parasite tree.
     * 
     * @return The generated parasite tree.
     */
    public Tree getParasiteTree() {
        return parasiteTree;
    }


    /**
     * Returns the mapping that describes how the parasite tree was generated while following the
     * host tree.
     * 
     * @return The mapping that describes how the parasite tree was generated while following the
     *         host tree.
     */
    public ArrayList<Association> getMapping() {
        return associations;
    }


    /**
     * Increases the loss counter in one unit.
     */
    protected void increaseCospeciationCounter() {
        cospeciations++;
    }


    /**
     * Increases the loss counter in one unit.
     */
    protected void increaseDuplicationCounter() {
        duplications++;
    }


    /**
     * Increases the loss counter in one unit.
     */
    protected void increaseHostSwitchCounter() {
        hostSwitches++;
    }


    /**
     * Increases the loss counter in one unit.
     */
    protected void increaseLossCounter() {
        losses++;
    }
    
    
    /**
     * Increases the multiple associations counter in one unit.
     */
    protected void increaseMultipleAssociationsCounter() {
        multipleAssociations++;
    }
    
    
    /**
     * Increases the verticalSpread counter in one unit.
     */
    protected void increaseVerticalSpreadCounter() {
    	verticalSpread++;
		
	}

    
    /**
     * Increases the horizontalSpread counter in one unit.
     */
    protected void increaseHorizontalSpreadCounter() {
    	horizontalSpread++;
		
	}


    /**
     * Increases the loss counter by the given number.
     */
    protected void increaseLossCounter(int numberOfLosses) {
        losses += numberOfLosses;
    }
    
    /**
     * Increases the multiple asociations counter by the given number.
     */
    protected void increaseMultipleAssociationsCounter(int numberOfMultipleAssociations) {
        losses += numberOfMultipleAssociations;
    }


    /**
     * Increases the number of leaves by the given number.
     */
    protected void increaseNumberOfLeaves(int leaves) {
        numberOfLeaves += leaves;
    }


    /**
     * Increases the number of leaves by one unit.
     */
    protected void increaseNumberOfLeaves() {
        numberOfLeaves++;
    }


    /**
     * Adds one assocation to the scenario.
     * 
     * @param association
     *        Association object that is going to be added.
     */
    protected void addAssociation(Association association) {
        associations.add(association);
        //System.out.println("associations: " + associations);
    }
    

    /**
     * Returns the frequency of cospeciations.
     * 
     * @return The frequency of cospeciations.
     */
    public double getFrequencyOfCospeciations() {
        double total = cospeciations + duplications + hostSwitches + losses;
        return cospeciations / total;
    }

	public double getFrequencyOfVerticalSpread() {
		return 0;
	}


	
	public double getSumDistanceVerticalSpread() {
		return 0;
	}


    /**
     * Returns the frequency of duplications.
     * 
     * @return The frequency of duplications.
     */
    public double getFrequencyOfDuplications() {
        double total = cospeciations + duplications + hostSwitches + losses;
        return duplications / total;
    }


    /**
     * Returns the frequency of host switches.
     * 
     * @return The frequency of host switches.
     */
    public double getFrequencyOfHostSwitches() {
        double total = cospeciations + duplications + hostSwitches + losses;
        return hostSwitches / total;
    }


    /**
     * Returns the frequency of losses.
     * 
     * @return The frequency of losses.
     */
    public double getFrequencyOfLosses() {
        double total = cospeciations + duplications + hostSwitches + losses;
        return losses / total;
    }



    /**
     * Computes the euclidean distance between the probability vector and the observed frequency
     * vector for this scenario.
     * 
     * @param pc
     *        Co-speciation probability used in the simulation model.
     * @param pd
     *        Duplication probability used in the simulation model.
     * @param ps
     *        Host-switch probability used in the simulation model.
     * @param pl
     *        Loss probability used in the simulation model.
     * @return The euclidean distance between the probability vector and the observed frequency
     *         vector for this scenario.
     */
    public double computeEuclideanDistance(double pc, double pd, double ps, double pl) {
        double total = cospeciations + duplications + hostSwitches + losses;
        double fc = cospeciations / total;
        double fd = duplications / total;
        double fs = hostSwitches / total;
        double fl = losses / total;
        double c = pc - fc;
        double d = pd - fd;
        double s = ps - fs;
        double l = pl - fl;
        return Math.sqrt((c * c) + (d * d) + (s * s) + (l * l));
    }


    /**
     * Computes the euclidean distance between the probability vector and the observed frequency
     * vector for this scenario.
     * 
     * @param pc
     *        Co-speciation probability used in the simulation model.
     * @param pd
     *        Duplication probability used in the simulation model.
     * @param ps
     *        Host-switch probability used in the simulation model.
     * @return The euclidean distance between the probability vector and the observed frequency
     *         vector for this scenario.
     */
    public double computeEuclideanDistance(double pc, double pd, double ps) {
        double total = cospeciations + duplications + hostSwitches;
        double fc = cospeciations / total;
        double fd = duplications / total;
        double fs = hostSwitches / total;
        double c = pc - fc;
        double d = pd - fd;
        double s = ps - fs;
        return Math.sqrt((c * c) + (d * d) + (s * s));
    }


    /**
     * Returns the hashmap which keeps the mapping parasite node to host node (key values).
     * 
     * @return The hashmap which keeps the mapping parasite node to host node (key values).
     */
    public HashMap<Integer, List<Integer>> getParasiteToHostMapping() {
        return parasiteToHost;
    }


    
    public void countMaxDistanceJump(double values) {
    	if (values > maxDistanceJump)
    		maxDistanceJump = values;
    }
    
    public void addDistanceJump(double values) {
    	sumDistanceJump = sumDistanceJump + values;
    }
    
    public void addDistanceVerticalSpread(double values) {
    	sumDistanceHorizontalSpread = sumDistanceHorizontalSpread + values;
    }
    
    public void jumpDistanceTot(double values){
    	jumpDistanceTot.add(values);
    }
    
    public double getMaxDistanceJump() {
    	return maxDistanceJump;
    }
    
    
    public double getSumDistanceJump() {
    	return sumDistanceJump;
    }
    
    public Double[] getJumpDistanceTot() {
    	return jumpDistanceTot.toArray(new Double[jumpDistanceTot.size()]);
    }
    
    public double getSumDistanceHorizontalSpread() {
    	return sumDistanceHorizontalSpread;
    }
}
