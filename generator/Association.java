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

public class Association {

    /** Cospeciation event. */
    public static final int COSPECIATION = 1;

    /** Duplication event. */
    public static final int DUPLICATION = 2;

    /** Host switch event. */
    public static final int HOSTSWITCH = 3;

    /** Horizontal Spread event. */
    public static final int HORIZONTALSPREAD = 4;
    
    /** Vertical Spread event. */
    public static final int VERTICALSPREAD = 5;
    
    /** Stop event (parasite node mapped to a leaf node). */
    public static final int STOP = 6;

    /** Parasite Node. */
    private TreeNode parasite;

    /** Host node or Host edge. */
    private Object hostElement;

    /** Event that is associated to this association. */
    private int associatedEvent;


    /**
     * Constructor
     * 
     * @param parasite
     *        Parasite Node.
     * @param hostElement
     *        Host node or Host edge.
     * @param associatedEvent
     *        Event that is associated to this association.
     */
    public Association(TreeNode parasite, Object hostElement, int associatedEvent) {
        if (hostElement instanceof TreeNode || hostElement instanceof Edge) {
            this.parasite = parasite;
            this.hostElement = hostElement;
            switch (associatedEvent) {
                case COSPECIATION:
                case DUPLICATION:
                case HOSTSWITCH:
                case HORIZONTALSPREAD:
                case VERTICALSPREAD:
                case STOP:
                    this.associatedEvent = associatedEvent;
                    break;
                default:
                    throw new RuntimeException("Invalid event type!");
            }

        } else {
            throw new RuntimeException("Invalid object instance!");
        }
    }


    /**
     * Return the parasite node.
     * 
     * @return The parasite node.
     */
    public TreeNode getParasite() {
        return parasite;
    }


    /**
     * Return the host element: node or edge.
     * 
     * @return The host element: node or edge.
     */
    public Object getHostElement() {
        return hostElement;
    }


    /**
     * Return a string representation of this object.
     * 
     * @return A string representation of this object.
     */
    public String toString() {
        String event = "";
        switch (associatedEvent) {
            case COSPECIATION:
                event = "{C}";
                break;
            case DUPLICATION:
                event = "{D}";
                break;
            case HOSTSWITCH:
                event = "{H}";
                break;
            case HORIZONTALSPREAD:
            	event = "{HS}";
            	break;
            case VERTICALSPREAD:
                event = "{VS}";
                break;
            case STOP:
                event = "{S}";
                break;
        }
        return event + parasite + "@" + hostElement;
    }

}
