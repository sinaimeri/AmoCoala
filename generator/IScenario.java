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

import  trees.Tree;
import  trees.TreeNode;

public interface IScenario {

    /**
     * Return the parasite tree root node.
     * 
     * @return The parasite tree root node.
     */
    public TreeNode getRootNode();


    /**
     * Return the number of leaves.
     * 
     * @return The number of leaves.
     */
    public int getNumberOfLeaves();


    /**
     * Return the number of cospeciations.
     * 
     * @return The number of cospeciations.
     */
    public int getNumberOfCospeciations();


    /**
     * Return the number of duplications.
     * 
     * @return The number of duplications.
     */
    public int getNumberOfDuplications();


    /**
     * Return the number of host switches.
     * 
     * @return The number of host switches.
     */
    public int getNumberOfHostSwitches();


    /**
     * Return the number of losses.
     * 
     * @return The number of losses.
     */
    public int getNumberOfLosses();

    /**
     * Return the number of multiple Associations.
     * 
     * @return The number of multiple Associations.
     */
    public int getNumberOfMultipleAssociations();
    
    /**
     * Return the number of Vertical Spread.
     * 
     * @return The number of Vertical Spread.
     */
    public int getNumberOfVerticalSpread();
    
    
    /**
     * Return the number of Horizontal Spread.
     * 
     * @return The number of Horizontal Spread.
     */
    public int getNumberOfHorizontalSpread();
    

    /**
     * Return the generated parasite tree.
     * 
     * @return The generated parasite tree.
     */
    public Tree getParasiteTree();


    /**
     * Return the mapping that describes how the parasite tree was generated while following the
     * host tree.
     * 
     * @return The mapping that describes how the parasite tree was generated while following the
     *         host tree.
     */
    public ArrayList<Association> getMapping();


    /**
     * Returns the frequency of cospeciations.
     * 
     * @return The frequency of cospeciations.
     */
    public double getFrequencyOfCospeciations();


    /**
     * Returns the frequency of duplications.
     * 
     * @return The frequency of duplications.
     */
    public double getFrequencyOfDuplications();


    /**
     * Returns the frequency of host switches.
     * 
     * @return The frequency of host switches.
     */
    public double getFrequencyOfHostSwitches();


    /**
     * Returns the frequency of losses.
     * 
     * @return The frequency of losses.
     */
    public double getFrequencyOfLosses();
    
    /**
     * Returns the frequency of losses.
     * 
     * @return The frequency of losses.
     */
    public double getFrequencyOfVerticalSpread();

    
    /**
     * Returns the maximum distance of host switch between two host nodes.
     * 
     * @return The maximum distance of host switch between two host nodes.
     */
    public double getMaxDistanceJump();
    
    public double getSumDistanceJump();
    
    public Double[] getJumpDistanceTot();
    
    public double getSumDistanceVerticalSpread();
    
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
    public double computeEuclideanDistance(double pc, double pd, double ps, double pl);


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
    public double computeEuclideanDistance(double pc, double pd, double ps);

}
