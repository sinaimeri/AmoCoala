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

package  generator.defaultmodel;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.HashMap;

import  generator.Association;
import  generator.ScenarioAbst;
import  trees.TreeNode;
import  trees.Tree;

public class DefaultModelScenario extends ScenarioAbst {

    /**
     * Creates an empty scenario: no associations and a parasite tree with only one node (root
     * node).
     */
    public DefaultModelScenario() {
        super();
        initParasiteTree(true);
    }


    /**
     * Increments in one unit the loss counter.
     */
    public void loss() {
        increaseLossCounter();
    }


    /**
     * Adds an association that represents a co-speciation event.
     * 
     * @param parasiteNode
     *        The parasite node.
     * @param hostNode
     *        The host node.
     */
    public void cospeciation(TreeNode parasiteNode, TreeNode hostNode) {
        addAssociation(new Association(parasiteNode, hostNode, Association.COSPECIATION));
        //System.out.println("COSP : P " + parasiteNode.getLabel() + " ---- H " + hostNode.getLabel() );
        List<Integer> nodeKey = new ArrayList<Integer>();
        nodeKey.add(hostNode.getKey());
        //parasiteToHost.put(parasiteNode.getKey(), hostNode.getKey());
        parasiteToHost.put(parasiteNode.getKey(), nodeKey);
        increaseCospeciationCounter();
    }


    /**
     * Adds an association that represents a duplication event.
     * 
     * @param parasiteNode
     *        The parasite node.
     * @param hostNode
     *        The host node.
     */
    public void duplication(TreeNode parasiteNode, TreeNode hostNode) {
        addAssociation(new Association(parasiteNode, hostNode, Association.DUPLICATION));
        //System.out.println("DUP : P " + parasiteNode.getLabel() + " ---- H " + hostNode.getLabel() );
        List<Integer> nodeKey = new ArrayList<Integer>();
        nodeKey.add(hostNode.getKey());
        //parasiteToHost.put(parasiteNode.getKey(), hostNode.getKey());
        parasiteToHost.put(parasiteNode.getKey(), nodeKey);
        increaseDuplicationCounter();
    }


    /**
     * Adds an association that represents a host switch event which involves the parasite edge
     * (u,v) and their image in the host tree.
     * 
     * @param parasiteNodeU
     *        The parasite node u.
     * @param parasiteNodeV
     *        The parasite node v.
     * @param hostNodeU
     *        The host node where the parasite node u is being mapped..
     * @param hostNodeV
     *        The host node where the parasite node v is being mapped..
     */
    public void hostSwitch(TreeNode parasiteNodeU,
                           TreeNode parasiteNodeV,
                           TreeNode hostNodeU,
                           TreeNode hostNodeV) {
        addAssociation(new Association(parasiteNodeU, hostNodeU, Association.HOSTSWITCH));
        List<Integer> nodeUKey = new ArrayList<Integer>();
        nodeUKey.add(hostNodeU.getKey());
        List<Integer> nodeVKey = new ArrayList<Integer>();
        nodeVKey.add(hostNodeV.getKey());
        //parasiteToHost.put(parasiteNodeU.getKey(), hostNodeU.getKey());
        //parasiteToHost.put(parasiteNodeV.getKey(), hostNodeV.getKey());
        parasiteToHost.put(parasiteNodeU.getKey(), nodeUKey);
        parasiteToHost.put(parasiteNodeV.getKey(), nodeVKey);
        
        //System.out.println("HS : P " + parasiteNodeU.getLabel() + " ---- H " + hostNodeU.getLabel() + "--CHILD : " + hostNodeU.getKey() + " " + hostNodeV.getKey() );
        increaseHostSwitchCounter();
        
       Tree hostTree;
    }

    /**
     * Adds an association that represents a host switch event which involves the parasite edge
     * (u,v) and their image in the host tree.
     * 
     * @param parasiteNodeU
     *        The parasite node u.
     * @param parasiteNodeV
     *        The parasite node v.
     * @param hostNodeU
     *        The host node where the parasite node u is being mapped..
     * @param hostNodeV
     *        The host node where the parasite node v is being mapped..
     */
    public void hostSwitch(TreeNode parasiteNodeU,
                           TreeNode parasiteNodeV,
                           TreeNode hostNodeU,
                           TreeNode hostNodeV, Tree hostTree) {
        addAssociation(new Association(parasiteNodeU, hostNodeU, Association.HOSTSWITCH));
        List<Integer> nodeUKey = new ArrayList<Integer>();
        nodeUKey.add(hostNodeU.getKey());
        List<Integer> nodeVKey = new ArrayList<Integer>();
        nodeVKey.add(hostNodeV.getKey());
        //parasiteToHost.put(parasiteNodeU.getKey(), hostNodeU.getKey());
        //parasiteToHost.put(parasiteNodeV.getKey(), hostNodeV.getKey());
        parasiteToHost.put(parasiteNodeU.getKey(), nodeUKey);
        parasiteToHost.put(parasiteNodeV.getKey(), nodeVKey);
        int lca_Dist = hostTree.getLCA(hostNodeU.getKey(), hostNodeV.getKey()).getHeight();
        int hostNodeU_Dist = hostNodeU.getHeight();
        int hostNodeV_Dist = hostNodeV.getHeight();
        double jumpDistance = (double)(hostNodeU_Dist-lca_Dist)+(hostNodeV_Dist-lca_Dist);  
        //System.out.println("HS : P " + parasiteNodeU.getLabel() + " ---- H " + hostNodeU.getLabel() + "--CHILD : " + hostNodeU.getKey() + " " + hostNodeV.getKey() );

        //countMaxDistanceJump(jumpDistance);
        addDistanceJump(jumpDistance);
        jumpDistanceTot(jumpDistance);
        increaseHostSwitchCounter();
    }

    
    public void  vertical_spread(TreeNode parasiteNode, TreeNode hostNode, boolean counter) {
        addAssociation(new Association(parasiteNode, hostNode, Association.VERTICALSPREAD));
      //parasiteToHost.put(parasiteNode.getKey(), hostNode.getKey());
        //System.out.println("FREZ : P " + parasiteNode.getLabel() + " ---- H " + hostNode.getLabel() );        
        List<Integer> hostKey = new ArrayList<Integer>();
       
        if(parasiteToHost.containsKey(parasiteNode.getKey())){
        	hostKey = parasiteToHost.get(parasiteNode.getKey());
        	hostKey.add(hostNode.getKey());
        }
        else{
        	hostKey.add(hostNode.getKey());
        }
        parasiteToHost.put(parasiteNode.getKey(), hostKey);
        if (counter)
        	increaseVerticalSpreadCounter();
    }
    
 
    

	/**
     * Adds an association that represents a stop event.
     * 
     * @param parasiteNode
     *        The parasite node.
     * @param hostNode
     *        The host node.
     */
    public void stop(TreeNode parasiteNode, TreeNode hostNode) {
        addAssociation(new Association(parasiteNode, hostNode, Association.STOP));
        parasiteNode.setLabel(hostNode.getLabel());
        List<Integer> nodeKey = new ArrayList<Integer>();
        nodeKey.add(hostNode.getKey());
        //parasiteToHost.put(parasiteNode.getKey(), hostNode.getKey());
        parasiteToHost.put(parasiteNode.getKey(), nodeKey);
    }
    
    /**
     * Adds a multiple-association that represents a stop event.
     * 
     * @param parasiteNode
     *        The parasite node.
     * @param hostNode
     *        The host node.
     */
    public void stop(TreeNode parasiteNode, TreeNode hostNode, Tree hostTree, boolean spread) {
    	
        addAssociation(new Association(parasiteNode, hostNode, Association.STOP));
        
    	//System.out.println("-------------");
    	//System.out.println("parasiteNode " + parasiteNode + " hostNode.getLabel():" + hostNode.getLabel());
    	//parasiteNode.setLabel(hostNode.getLabel());
        parasiteNode.setLabelSpread(hostNode.getLabel(), parasiteNode, hostTree);
        //System.out.println("parasiteNode " + parasiteNode);
        List<Integer> nodeKey = new ArrayList<Integer>();
        nodeKey.add(hostNode.getKey());
        //parasiteToHost.put(parasiteNode.getKey(), hostNode.getKey());
        parasiteToHost.put(parasiteNode.getKey(), nodeKey);
      
    }
    
    
   


    /**
     * Returns the hashmap which keeps the mapping parasite node to host node (key values).
     * 
     * @return The hashmap which keeps the mapping parasite node to host node (key values).
     */
    public HashMap<Integer, List<Integer>> getParasiteToHostMapping() {
        return parasiteToHost;
    }


    /**
     * Creates two children for the given node.
     * 
     * @param node
     *        The node that will receive two children.
     * @return Return an array containing the two new nodes.
     */
    public TreeNode[] createChildren(TreeNode node) {
        TreeNode child0 = getParasiteTree().createNode();
        TreeNode child1 = getParasiteTree().createNode();
        TreeNode[] children = {child0, child1};
        getParasiteTree().addChild(node.getKey(), child0.getKey(), 0);
        getParasiteTree().addChild(node.getKey(), child1.getKey(), 1);
        increaseNumberOfLeaves();
        /* Update the lists of descendants and ancestorss. */
        while (node != null) {
            child0.addAncestor(node);
            child1.addAncestor(node);
            node.addDescendant(child0);
            node.addDescendant(child1);
            node = node.getParent();
        }
        return children;
    }


    /**
     * Adds an association that represents a horizontal spread event which involves the parasite edge
     * (u) and its image in the host tree.
     * 
     * @param parasiteNodeU
     *        The parasite node u.
     * @param hostNode
     *        The host node where the parasite node u is being mapped..
     */
    public void horizontal_spread(TreeNode parasiteNode,
    		TreeNode hostNodeU,
    		TreeNode hostNodeV) {
    	addAssociation(new Association(parasiteNode, hostNodeU, Association.HORIZONTALSPREAD));
    	List<Integer> nodeUKey = new ArrayList<Integer>();
    	nodeUKey.add(hostNodeU.getKey());
    	List<Integer> nodeVKey = new ArrayList<Integer>();
    	nodeVKey.add(hostNodeV.getKey());
    	
    	List<Integer> hostKeyU = new ArrayList<Integer>();
    	if(parasiteToHost.containsKey(parasiteNode.getKey())){
        	hostKeyU = parasiteToHost.get(parasiteNode.getKey());
        	hostKeyU.add(hostNodeU.getKey());
        }
        else{
        	hostKeyU.add(hostNodeU.getKey());
        }
    	
    	List<Integer> hostKeyV = new ArrayList<Integer>();
    	if(parasiteToHost.containsKey(parasiteNode.getKey())){
        	hostKeyV = parasiteToHost.get(parasiteNode.getKey());
        	hostKeyV.add(hostNodeV.getKey());
        }
        else{
        	hostKeyV.add(hostNodeV.getKey());
        }
    	parasiteToHost.put(parasiteNode.getKey(), nodeUKey);
    	parasiteToHost.put(parasiteNode.getKey(), nodeVKey);

    	//System.out.println("HS : P " + parasiteNodeU.getLabel() + " ---- H " + hostNodeU.getLabel() + "--CHILD : " + hostNodeU.getKey() + " " + hostNodeV.getKey() );
    	increaseHorizontalSpreadCounter();

    	Tree hostTree;
    }


	
  
}
