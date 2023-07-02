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

import  graph.Edge;
import  trees.TreeNode;

public class Position {

    /** Parasite's tree node */
    private TreeNode parasiteNode;

    /** Host's tree edge */
    private Edge hostEdge;


    /**
     * Constructor: Creates a position that relates the given parasite node to the given host edge.
     * 
     * @param parasiteNode
     *        Parasite's tree node.
     * @param hostEdge
     *        Host's tree edge
     */
    public Position(TreeNode parasiteNode, Edge hostEdge) {
        this.parasiteNode = parasiteNode;
        this.hostEdge = hostEdge;
    }


    /**
     * Returns the parasite tree node.
     * 
     * @return The parasite tree node.
     */
    public TreeNode getParasiteNode() {
        return parasiteNode;
    }


    /**
     * Returns the host tree edge.
     * 
     * @return The host tree edge.
     */
    public Edge getHostEdge() {
        return hostEdge;
    }


    /**
     * Return a String which represents this object.
     * 
     * @return A String representation of this object.
     */
    public String toString() {
        return "{" + this.parasiteNode.getKey() + ";" + this.hostEdge.toString() + "}";
    }

}