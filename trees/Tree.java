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

package  trees;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map.Entry;
import java.util.Set;
import java.util.Map;

import com.google.common.collect.HashMultiset;
import com.google.common.collect.Multisets;


public class Tree {

	private final char OPEN_BRACKET = '(';

	private final char CLOSE_BRACKET = ')';

	private final char CHILD_SEPARATOR = ',';

	/**
	 * Variable which keeps a counter for generating key identifier for new
	 * nodes.
	 */
	private int key;

	/** Tree's height */
	private int height;

	/** Root of the tree. */
	private TreeNode root;

	/** List of nodes of the tree. */
	private ArrayList<TreeNode> nodes;

	/** List of leaf nodes of the tree. */
	private ArrayList<TreeNode> leaves;

	/**
	 * HashMap which keeps the mapping key -> TreeNode object for each created
	 * node.
	 */
	private HashMap<Integer, TreeNode> mapKeyToNode;

	/**
	 * HashMap which keeps the mapping label -> TreeNode object for each created
	 * node.
	 */
	private HashMap<String, TreeNode> mapLabelToNode;

	/** Table containing all LCA of every pair of leaf nodes of the tree. */
	private TreeNode[][] lcaTable;

	/** Array containing all pVerticalSpread of each node of the tree */
	private HashMap<TreeNode, Double> mapPVerticalSpread;

	/** HasMap containing the number of parasites in a subtree */
	private HashMap<TreeNode, HashMap<TreeNode, Integer>> nbParasitesSubTree;
	
	/** HasMap containing the associations between host and parasites */
	private HashMap<TreeNode, List<TreeNode>> mapVerticalSpreadAssociations;
	
	/** HasMap containing the probability of Vertical Spread for each host node */
	private HashMap<TreeNode, Double> verticalSpreadMap;
	
	/** matrix containing the probability of Horizontal Spread for each pair of nodes */
	private Double[][] horizontalSpreadMatrix;
	
	/** number of multiple associations for the tree*/
	private int nbMultipleAssociations;
	
	

	/**
	 * Creates a tree object.
	 * 
	 * @param createRoot
	 *            If true, creates a root node and sets it with height 1.
	 */
	public Tree(boolean createRoot) {
		this.key = 0;
		this.mapKeyToNode = new HashMap<Integer, TreeNode>();
		if (createRoot)
			this.root = createNode();
		else
			this.root = null;
		this.height = -1;
		this.nodes = null;
		this.leaves = null;
		this.lcaTable = null;
		this.mapLabelToNode = null;
		this.mapPVerticalSpread = null;

	}

	/**
	 * Returns the node which has the given key.
	 * 
	 * @param key
	 *            The node key.
	 * @return The node which has the given key.
	 */
	public TreeNode getNode(int key) {
		TreeNode node = mapKeyToNode.get(key);
		if (node == null)
			throw new RuntimeException("Invalid key identifier");
		return node;
	}
	
	/**
	 * Returns the node which has the given key.
	 * 
	 * @param key
	 *            The node key.
	 * @return The node which has the given key.
	 */
	public List<TreeNode> getNode(List<Integer> key) {
		List<TreeNode> nodefinal = new ArrayList<TreeNode>();
		for(Integer k: key){
			TreeNode node = mapKeyToNode.get(k);
			if (node == null)
				throw new RuntimeException("Invalid key identifier");
			nodefinal.add(node);
		}
		
		return nodefinal;
	}

	/**
	 * Returns the root of the tree.
	 * 
	 * @return The root of the tree.
	 */
	public TreeNode getRoot() {
		return root;
	}

	/**
	 * Set the root node of this tree.
	 * 
	 * @param key
	 *            The key of the new root node.
	 */
	public void setRoot(int key) {
		TreeNode node = mapKeyToNode.get(key);
		if (node == null)
			throw new RuntimeException("Invalid key identifier");
		if (node.getParent() != null)
			throw new RuntimeException(
					"The root of the tree cannot have a parent node.");
		this.root = node;
		this.nodes = null;
		this.height = -1;
		this.leaves = null;
		this.lcaTable = null;
		this.mapLabelToNode = null;
	}

	/**
	 * Creates a TreeNode object with unique identifier.
	 * 
	 * @return A TreeNode object with unique identifier.
	 */
	public TreeNode createNode() {
		TreeNode node = new TreeNode(key++);
		mapKeyToNode.put(node.getKey(), node);
		return node;
	}

	/**
	 * Creates a TreeNode object with unique identifier and the given label.
	 * 
	 * @return A TreeNode object with unique identifier and the given label.
	 */
	public TreeNode createNode(String label) {
		TreeNode node = new TreeNode(key++, label);
		mapKeyToNode.put(node.getKey(), node);
		return node;
	}

	/**
	 * Add a child into the given node at the give position.
	 * 
	 * @param node
	 *            Node that will receive the child.
	 * @param child
	 *            The new child of the node.
	 * @param position
	 *            0 - left child; 1 - right child.
	 */
	public void addChild(int node, int child, int position) {
		TreeNode parentNode = mapKeyToNode.get(node);
		TreeNode childNode = mapKeyToNode.get(child);
		if (parentNode != null && childNode != null)
			parentNode.addChild(childNode, position);
		this.nodes = null;
		this.height = -1;
		this.leaves = null;
		this.lcaTable = null;
		this.mapLabelToNode = null;
	}

	/**
	 * Returns the height of the tree.
	 * 
	 * @return The height of the tree.
	 */
	public int getHeight() {
		if (height == -1)
			updateAttributes();
		return height;
	}

	/**
	 * Returns the total number of nodes of the tree (excluding eventual fake
	 * root).
	 * 
	 * @return The total number of nodes of the tree (excluding eventual fake
	 *         root).
	 */
	public int getNumberOfNodes() {
		if (nodes == null)
			updateAttributes();
		return nodes.size();
	}

	/**
	 * Returns the total number of leaf nodes of the tree.
	 * 
	 * @return The total number of leaf nodes of the tree.
	 */
	public int getNumberOfLeafNodes() {
		if (leaves == null)
			updateAttributes();
		return leaves.size();
	}

	
	/**
	 * Returns the total number of multiple associations of the tree.
	 * 
	 * @return The total number of multiple associations of the tree.
	 */
	public int getNumberOfMultipleAssociations() {
		return nbMultipleAssociations;
	}
	
	
	/**
	 * Count the number of multiple associations of the tree.
	 */
	public void updateNumberOfMultipleAssociations() {
		nbMultipleAssociations = 0;	
		for (int u=0; u<getNumberOfLeafNodes(); u++){
			nbMultipleAssociations = nbMultipleAssociations + (getLeafNodes()[u].getLabel().size()-1);
		}
	}

	/**
	 * Traverses the tree in post-order and update the attributes of this tree
	 * (list of nodes, list of leaves and tree height), excluding the eventual
	 * fake root.
	 */
	private void updateAttributes() {
		this.nodes = new ArrayList<TreeNode>();
		this.leaves = new ArrayList<TreeNode>();
		this.height = -1;
		if (root == null)
			throw new RuntimeException("The tree does not have a root node.");
		if (root.getNumberOfChildren() == 1) {
			/* RHO. */
			updateAttributes(root, 0);
		} else {
			/* Real root. */
			updateAttributes(root, 1);
		}
	}
	
	/**
	 * Recursive function that gives the pre-order tree  
	 * 
	 * @param node
	 *            node of the subtree that is going to be traversed.
	 * @param s
	 *            arrayList that contain the pre-order tree.
	 * @return preOrderString
	 * 			  arrayList that contain the pre-order tree.          
	 */
	public ArrayList<TreeNode> getPreorderTree(TreeNode node, ArrayList<TreeNode> s){
		ArrayList <TreeNode> preOrderString = new ArrayList<TreeNode>(s);
		if(node != null){
			preOrderString.add(node);
			preOrderString = getPreorderTree(node.getChildren()[0], preOrderString);
			preOrderString = getPreorderTree(node.getChildren()[1], preOrderString);
		}
		return preOrderString;
	}
	
	
	/**
	 * Recursive function that helps in the task of updating the attributes of
	 * this tree.
	 * 
	 * @param subtreeRoot
	 *            Root of the subtree that is going to be traversed.
	 * @param height
	 *            The current height.
	 */
	private void updateAttributes(TreeNode subtreeRoot, int height) {
		subtreeRoot.setHeight(height);
		if (subtreeRoot.isLeaf()) {
			leaves.add(subtreeRoot);
			if (this.height < height)
				this.height = height;
		} else {
			if (subtreeRoot.getChild(0) != null)
				updateAttributes(subtreeRoot.getChild(0), height + 1);
			if (subtreeRoot.getChild(1) != null)
				updateAttributes(subtreeRoot.getChild(1), height + 1);
		}
		nodes.add(subtreeRoot);
	}

	/**
	 * Generate a string representing the traversing of the tree in post order.
	 * This method uses the key values of the nodes.
	 * 
	 * @param node
	 *            The root of the tree
	 * @param addInternalNodes
	 *            If <code>true</code>, the newick string includes the internal
	 *            node labels. Otherwise, only leaf labels are present in the
	 *            returned string.
	 * @return A string representing the traversing of the tree in post order
	 *         (using the key of the nodes).
	 */
	private String postOrderStringKey(TreeNode node, boolean addInternalNodes) {
		if (node.isLeaf())
			return Integer.toString(node.getKey());
		if (node.getNumberOfChildren() == 2) {
			String ls = postOrderStringKey(node.getChild(0), addInternalNodes);
			String rs = postOrderStringKey(node.getChild(1), addInternalNodes);
			if (addInternalNodes)
				return OPEN_BRACKET + ls + CHILD_SEPARATOR + rs + CLOSE_BRACKET
						+ Integer.toString(node.getKey());
			else
				return OPEN_BRACKET + ls + CHILD_SEPARATOR + rs + CLOSE_BRACKET;
		} else if (node.getChild(0) != null) {
			String ls = postOrderStringKey(node.getChild(0), addInternalNodes);
			if (addInternalNodes)
				return OPEN_BRACKET + ls + CLOSE_BRACKET
						+ Integer.toString(node.getKey());
			else
				return OPEN_BRACKET + ls + CLOSE_BRACKET;
		} else {
			String rs = postOrderStringKey(node.getChild(1), addInternalNodes);
			if (addInternalNodes)
				return OPEN_BRACKET + rs + CLOSE_BRACKET
						+ Integer.toString(node.getKey());
			else
				return OPEN_BRACKET + rs + CLOSE_BRACKET;
		}
	}

	/**
	 * Generate a string representing the traversing of the tree in post order.
	 * This method uses the label values of the nodes.
	 * 
	 * @param n
	 *            The root of the tree
	 * @param addInternalNodes
	 *            If <code>true</code>, the newick string includes the internal
	 *            node labels. Otherwise, only leaf labels are present in the
	 *            returned string.
	 * @return A string representing the traversing of the tree in post order
	 *         (using the label of the nodes).
	 */
	private String postOrderStringLabel(TreeNode node, boolean addInternalNodes) {
		if (node.isLeaf())
			//return node.getLabel();
			return node.getLabel().toString();
		if (node.getNumberOfChildren() == 2) {
			String ls = postOrderStringLabel(node.getChild(0), addInternalNodes);
			String rs = postOrderStringLabel(node.getChild(1), addInternalNodes);
			if (addInternalNodes)
				return OPEN_BRACKET + ls + CHILD_SEPARATOR + rs + CLOSE_BRACKET
						+ node.getLabel();
			else
				return OPEN_BRACKET + ls + CHILD_SEPARATOR + rs + CLOSE_BRACKET;
		} else if (node.getChild(0) != null) {
			String ls = postOrderStringLabel(node.getChild(0), addInternalNodes);
			if (addInternalNodes)
				return OPEN_BRACKET + ls + CLOSE_BRACKET + node.getLabel();
			else
				return OPEN_BRACKET + ls + CLOSE_BRACKET;
		} else {
			String rs = postOrderStringLabel(node.getChild(1), addInternalNodes);
			if (addInternalNodes)
				return OPEN_BRACKET + rs + CLOSE_BRACKET + node.getLabel();
			else
				return OPEN_BRACKET + rs + CLOSE_BRACKET;
		}
	}

	/**
	 * Returns a String representation of this tree.
	 * 
	 * @return A String representation of this tree.
	 */
	public String toString() {
		return postOrderStringKey(root, true);
	}

	/**
	 * Returns a String representation of this tree.
	 * 
	 * @param label
	 *            If <code>true</code>, returns a String that contains the
	 *            labels of the nodes. Otherwise, returns a String that contains
	 *            the key values of the nodes.
	 * @return A String representation of this tree.
	 */
	public String toString(boolean label) {
		if (label)
			return postOrderStringLabel(root, true);
		return postOrderStringKey(root, true);
	}

	/**
	 * Returns a String representation of this tree.
	 * 
	 * @param label
	 *            If <code>true</code>, returns a String that contains the
	 *            labels of the nodes. Otherwise, returns a String that contains
	 *            the key values of the nodes.
	 * @param addInternalNodes
	 *            If <code>true</code>, the newick string includes the internal
	 *            node labels. Otherwise, only leaf labels are present in the
	 *            returned string.
	 * @return A String representation of this tree.
	 */
	public String toString(boolean label, boolean addInternalNodes) {
		if (label)
			return postOrderStringLabel(root, addInternalNodes);
		return postOrderStringKey(root, addInternalNodes);
	}

	/**
	 * Update the list of ancestor nodes of all nodes of the tree.
	 */
	public void updateAncestorList() {
		updateAncestorList(this.root, new HashSet<TreeNode>());
	}

	/**
	 * Recursive method to update the list of ancestors of all nodes of the
	 * tree.
	 * 
	 * @param node
	 *            Node to be updated.
	 * @param ancestors
	 *            List of ancestors.
	 */
	private void updateAncestorList(TreeNode node, HashSet<TreeNode> ancestors) {
		node.emptyAncestors();
		node.addAncestor(node);
		node.addAncestors(ancestors);
		TreeNode[] children = node.getChildren();
		if (children[0] != null) {
			updateAncestorList(children[0], node.getDescendants());
		}
		if (children[1] != null) {
			updateAncestorList(children[1], node.getDescendants());
		}
	}

	/**
	 * Update the list of descendant nodes of all nodes of the tree.
	 */
	public void updateDescendantList() {
		updateDescendantList(this.root);
	}

	/**
	 * Recursive method to update the list of descendants of all nodes of the
	 * tree.
	 * 
	 * @param node
	 *            Node to be updated.
	 * @param ancestors
	 *            List of ancestors.
	 */
	private void updateDescendantList(TreeNode node) {
		node.emptyDescendants();
		node.addDescendant(node);
		TreeNode[] children = node.getChildren();
		if (children[0] != null) {
			updateDescendantList(children[0]);
			node.addDescendants(children[0].getDescendants());
		}
		if (children[1] != null) {
			updateDescendantList(children[1]);
			node.addDescendants(children[1].getDescendants());
		}
	}

	/**
	 * Returns the list of descendants of the given tree node.
	 * 
	 * @param key
	 *            Identifier of the tree node.
	 * @return The list of descendants of the given tree node.
	 */
	public HashSet<TreeNode> getDescendants(int key, boolean preprocessed) {

		TreeNode node = mapKeyToNode.get(key);
		if (node == null)
			throw new RuntimeException("Invalid key identifier");

		if (preprocessed)
			return node.getDescendants();

		return getDescendantsList(node);
	}

	/**
	 * Returns the list of descendants of the given tree node.
	 * 
	 * @param node
	 *            Root of the subtree.
	 * @return The list of descendants of the given tree node.
	 */
	private HashSet<TreeNode> getDescendantsList(TreeNode node) {
		HashSet<TreeNode> descendantsList = new HashSet<TreeNode>();
		if (node.isLeaf()) {
			descendantsList.add(node);
		} else {
			TreeNode child0 = node.getChild(0);
			TreeNode child1 = node.getChild(1);
			if (child0 != null)
				descendantsList.addAll(getDescendantsList(child0));
			if (child1 != null)
				descendantsList.addAll(getDescendantsList(child1));
			descendantsList.add(node);
		}
		return descendantsList;
	}

	/**
	 * Returns the list of ancestors of the given tree node.
	 * 
	 * @param key
	 *            Identifier of the tree node.
	 * @return The list of ancestors of the given tree node.
	 */
	public HashSet<TreeNode> getAncestors(int key, boolean preprocessed) {
		TreeNode node = mapKeyToNode.get(key);
		if (node == null)
			throw new RuntimeException("Invalid key identifier");
		if (preprocessed)
			return node.getAncestors();
		return getAncestorsHashSet(node);
	}

	/**
	 * Returns the list of descendants of the given tree node.
	 * 
	 * @param key
	 *            Identifier of the tree node.
	 * @return The list of ancestors of the given tree node.
	 */
	public ArrayList<TreeNode> getAncestorsArrayList(int key) {
		TreeNode node = mapKeyToNode.get(key);
		if (node == null)
			throw new RuntimeException("Invalid key identifier");
		ArrayList<TreeNode> ancestors = new ArrayList<TreeNode>();
		while (node != null) {
			ancestors.add(node);
			node = node.getParent();
		}
		return ancestors;
	}
	
	/**
	 * Returns the list of descendants of the given tree node.
	 * 
	 * @param key
	 *            Identifier of the tree node.
	 * @return The list of descendants of the given tree node.
	 */
	public HashSet<TreeNode> getDescendants(TreeNode node/*, boolean preprocessed*/) {

		
		if (node == null)
			throw new RuntimeException("Invalid key identifier");

		/*if (preprocessed)
			return node.getDescendants();
		 */
		return getDescendantsList(node);
	}
	
	/**
	 * Returns the common ancestor of two given tree nodes.
	 * 
	 * @param key1
	 *            Identifier of the tree node1.
	 * @param key2
	 *            Identifier of the tree node2.
	 * @return The common ancestor of two given tree nodes.
	 */
	public TreeNode getCommonAncestors(int key1, int key2) {
		TreeNode node1 = mapKeyToNode.get(key1);
		TreeNode node2 = mapKeyToNode.get(key2);
		if (node1 == null || node2 ==null )
			throw new RuntimeException("Invalid key identifier");
		
		
		ArrayList<TreeNode> ancestorsNode1 = getAncestorsArrayList(node1.getKey());
		ArrayList<TreeNode> ancestorsNode2 = getAncestorsArrayList(node2.getKey());
		Object[] ancestorsNode1Array = ancestorsNode1.toArray();
		TreeNode commonAncestor = null;

		//Collections.reverse(ancestorsNode1);
		
		for (int i = 0; i < ancestorsNode1.size(); i++){
			if(ancestorsNode2.contains(ancestorsNode1.get(i))){
				commonAncestor = ancestorsNode1.get(i);
				break;
			}
		}
		return commonAncestor;
	}

	/**
	 * Returns the list of descendants of the given tree node.
	 * 
	 * @param key
	 *            Identifier of the tree node.
	 * @return The list of descendants of the given tree node.
	 */
	private HashSet<TreeNode> getAncestorsHashSet(TreeNode node) {
		HashSet<TreeNode> ancestors = new HashSet<TreeNode>();
		while (node != null) {
			ancestors.add(node);
			node = node.getParent();
		}return ancestors;
	}

	
	/**
	 * Returns the table which contains the least common ancestor relationship
	 * for all pair of tree nodes.
	 * 
	 * @return The table which contains the least common ancestor relationship
	 *         for all pair of tree nodes.
	 */
	private TreeNode[][] getLCATable() {
		if (lcaTable == null)
			createLCATable();
		return lcaTable;
	}

	
	/**
	 * Returns the least common ancestor of the given pair of tree nodes.
	 * 
	 * @param u
	 *            Node u
	 * @param v
	 *            Node v
	 * @return The least common ancestor of u and v.
	 */
	public TreeNode getLCA(int u, int v) {
		TreeNode nodeU = mapKeyToNode.get(u);
		TreeNode nodeV = mapKeyToNode.get(v);
		if (nodeU == null)
			throw new RuntimeException("Invalid key identifier for node u");
		if (nodeV == null)
			throw new RuntimeException("Invalid key identifier for node v");
		TreeNode[][] table = getLCATable();
		int indexU = nodes.indexOf(nodeU);
		int indexV = nodes.indexOf(nodeV);
		if (indexU == -1 || indexV == -1) {
			return null;
		}
		return table[indexU][indexV];
	}

	/**
	 * Creates the LCA table of the tree. Notice that, only the nodes that are
	 * in the main tree are included in the table.
	 */
	public void createLCATable() {
		if (nodes == null)
			updateAttributes();
		lcaTable = new TreeNode[nodes.size()][];
		for (int i = 0; i < nodes.size(); i++) {
			lcaTable[i] = new TreeNode[nodes.size()];
			for (int j = 0; j < nodes.size(); j++) {
				lcaTable[i][j] = lca(nodes.get(i), nodes.get(j));
			}
		}
	}

	
	/**
	 * Returns the LCA of the two given nodes. Warning 1: This method does not
	 * use the lca table and performs the search for the LCA based on the height
	 * of the nodes. Warning 2: Only nodes which are in the main tree (whose
	 * root node is assigned to the attribute root of this class) are
	 * considered.
	 * 
	 * @param u
	 *            Node u
	 * @param v
	 *            Node v
	 * @return The LCA of the two given nodes.
	 */
	private TreeNode lca(TreeNode u, TreeNode v) {

		if (!nodes.contains(u) || !nodes.contains(v))
			return null;

		/* First put the nodes into the same height. */
		while (u.getHeight() != v.getHeight()) {
			if (u.getHeight() > v.getHeight()) {
				u = u.getParent();
			} else {
				v = v.getParent();
			}
		}

		if (u == v)
			return u;

		/* Keep looking for the least common ancestor. */
		while (u.getParent().getKey() != v.getParent().getKey()) {
			u = u.getParent();
			v = v.getParent();
		}

		return u.getParent();
	}
	
	/**
	 * Creates the associations between Host and Parasite nodes for the Horizontal spread.
	 */
	public void createHORIZONTALSPREADAssociations(
			HashMap<String, List<String>> mappingParasiteHost, Tree parasiteTree) {
		TreeNode[] nodesNi = this.getNodes();
		TreeNode[] nodesNj = this.getNodes();
		
		int m = this.getNumberOfNodes();
		horizontalSpreadMatrix = new Double[m][m];
		//initialize the horizontalSpreadMatrix
		for (int i=0; i<m; i++){
			for (int j=0; j<m; j++){
				horizontalSpreadMatrix[i][j]=-1.0;
			}
			
		}
		
		
		
		//insert the horizontal spread value in the matrix 
		for(TreeNode ni: nodesNi){
			double probab = 0;
			for(TreeNode nj: nodesNi){
				if (ni == nj /*|| ni.isLeaf()*/){
					horizontalSpreadMatrix[ni.getKey()][nj.getKey()]=0.0;
				}
				//calculate the probability of horizontal spread for two leaves
				else if (horizontalSpreadMatrix[ni.getKey()][nj.getKey()] == -1){
					//if ni and nj are comparable
					if (this.getAncestors(ni.getKey(), false).contains(nj) || this.getAncestors(nj.getKey(), false).contains(ni)){
						horizontalSpreadMatrix[ni.getKey()][nj.getKey()]=0.0;
						horizontalSpreadMatrix[nj.getKey()][ni.getKey()]=0.0;
					}
					else{
						double cm = commonParaNodes(ni, nj, mappingParasiteHost);
						double unionp = commonParaClades(ni, nj,mapVerticalSpreadAssociations);
						horizontalSpreadMatrix[ni.getKey()][nj.getKey()]=cm/unionp;
						probab = probab + (horizontalSpreadMatrix[ni.getKey()][nj.getKey()] * nj.getPverticalSpread() );
						
					}
				}//else if	
				
			}//for
			probab = probab * ni.getPverticalSpread();
			//the maximum horizontal-vertical probability has to be 1
			if(probab > 1)
					probab=1;
			ni.setShvProbability(probab);
            //else check for <0
		}//for
		
	}
	

	/**
	 * Give the number of common parasites between two nodes.
	 * 
	 * @param ni
	 *            node i.
	 * @param nj
	 *            node j.
	 * @param mappingParasiteHost
	 *            HashMap with the mapping between Parasites and Hosts tree.
	 * @return An integer of  of common parasites between two nodes.
	 */
	public int commonParaNodes(TreeNode ni, TreeNode nj, HashMap<String, List<String>> mappingParasiteHost){
		List<TreeNode> associationsNi = mapVerticalSpreadAssociations.get(ni);
		List<TreeNode> associationsNj = mapVerticalSpreadAssociations.get(nj);
		int commonPara = 0;
		if(associationsNi.size() <= associationsNj.size()){
			for(TreeNode p : associationsNi){
				if (associationsNj.contains(p))
					commonPara ++;
			}
		}
		else{
			for(TreeNode p : associationsNj){
				if (associationsNi.contains(p))
					commonPara ++;
			}
		}
		return commonPara;
		
	}
	
	
	/**
	 * Give the number of common parasites between two clades.
	 * 
	 * @param ci
	 *            clade i.
	 * @param cj
	 *            clade j.
	 * @param mapVerticalSpreadAssociations
	 *            HashMap with the mapping between Parasites and Hosts tree.
	 * @return An integer of  of common parasites between two clades.
	 */
	public int commonParaClades(TreeNode ci, TreeNode cj, HashMap<TreeNode, List<TreeNode>> mapVerticalSpreadAssociations){
		int result = mapVerticalSpreadAssociations.get(ci).size() + mapVerticalSpreadAssociations.get(cj).size();
		//if the clades has at least 1 parasite mapped
		if((mapVerticalSpreadAssociations.get(ci).size() != 0 ) && (mapVerticalSpreadAssociations.get(cj).size() != 0 )){
			if(mapVerticalSpreadAssociations.get(ci).size() <= mapVerticalSpreadAssociations.get(cj).size()){
				for(TreeNode n : mapVerticalSpreadAssociations.get(ci)){
					if(mapVerticalSpreadAssociations.get(cj).contains(n))
						result--;
				}
			}
			else{
				for(TreeNode n : mapVerticalSpreadAssociations.get(cj)){
					if(mapVerticalSpreadAssociations.get(ci).contains(n))
						result--;
				}
			}
		}
		return result;
	}
		
	
	
	/**
	 * Set the probability of horizontal-vertical spread for the incomparables nodes on ni.
	 * 
	 * @param ni
	 *            node i.
	 * @param incomparablesNodes
	 *            ArrayList with the incomparables nodes of ni.
	 *            
	 * @return An HashMap with the probabilities for the incomparables nodes to ni.
	 */
	public HashMap<TreeNode, Double> pHorizontalSpread_IncomparablesNodes(TreeNode ni, ArrayList<TreeNode> incomparablesNodes){
		HashMap<TreeNode, Double> pIncomparablesNodes = new HashMap<TreeNode, Double>();
		if (incomparablesNodes.isEmpty())
			return pIncomparablesNodes;
			
		double denominator = 0.0;
		for (TreeNode nj : incomparablesNodes){
			denominator = denominator + (horizontalSpreadMatrix[ni.getKey()][nj.getKey()] * verticalSpreadMap.get(nj));
		}
		denominator = denominator * verticalSpreadMap.get(ni) ;
		double control = 0.0;
		for (TreeNode nj : incomparablesNodes){
			double probab = 0;
			probab = horizontalSpreadMatrix[ni.getKey()][nj.getKey()] * verticalSpreadMap.get(ni) * verticalSpreadMap.get(nj);
			probab = probab/denominator;
			pIncomparablesNodes.put(nj, probab);
			control = control + probab;
		}
		
		//if(control!=1)
		//	System.out.println("ERROR p_invasion");
		
		return pIncomparablesNodes;
	}
	
	

	/**
	 * Creates the associations between Host and Parasite nodes for the Vertical Spread.
	 * 
	 * @param mappingParasiteHost
	 *            HashMap with the mapping between Parasites and Hosts tree.
	 * @param parasiteTree
	 * 			  Parasite Tree
	 */
	public void createVERTICALSPREADAssociations(
			HashMap<String, List<String>> mappingParasiteHost, Tree parasiteTree) {
		
		TreeNode[] leafNodesH = this.getLeafNodes();
		mapVerticalSpreadAssociations = new HashMap<TreeNode, List<TreeNode>>();
		nbParasitesSubTree = new HashMap<TreeNode, HashMap<TreeNode, Integer>>();
		TreeNode[] leafNodesP = parasiteTree.getLeafNodes();
		for (int i = 0; i < leafNodesH.length; i++){
			List<TreeNode> associatedParasite = new ArrayList<TreeNode>();
			HashMap<TreeNode, Integer> nbParasites = new HashMap<TreeNode, Integer>();
			
			
			for (int j = 0; j < leafNodesP.length; j++){
				//if host label is in the list of associated Labels
				if(leafNodesP[j].getLabel().size()>1){
					
					for (int k = 0; k < leafNodesP[j].getLabel().size(); k++) {
						String currentPLabel = leafNodesP[j].getLabel().get(k);
						if(leafNodesH[i].getLabel().get(0).equals(currentPLabel)){
							associatedParasite.add(leafNodesP[j]);
							nbParasites.put(leafNodesP[j], 1);
						}
					}
				}
				else {
					if(leafNodesH[i].getLabel().equals(leafNodesP[j].getLabel())){
						associatedParasite.add(leafNodesP[j]);
						nbParasites.put(leafNodesP[j], 1);
					}
				}
				
			}
			mapVerticalSpreadAssociations.put(leafNodesH[i], associatedParasite);
			nbParasitesSubTree.put(leafNodesH[i], nbParasites);
		}
		addVerticalSpread(this.getRoot());
	}

	
	/**
	 * Creates the associations between Host and Parasite nodes. Only for internal nodes
	 */
	public void addVerticalSpread(TreeNode node) {
		if (!node.isLeaf()) {

			TreeNode leftChild = node.getChild(0);
			TreeNode rightChild = node.getChild(1);

			if (leftChild != null && rightChild != null) {
				if (!mapVerticalSpreadAssociations.containsKey(leftChild))
					addVerticalSpread(leftChild);
				if (!mapVerticalSpreadAssociations.containsKey(rightChild))
					addVerticalSpread(rightChild);

				List<TreeNode> leftChildAssoc = mapVerticalSpreadAssociations
						.get(leftChild);
				List<TreeNode> rightChildAssoc = mapVerticalSpreadAssociations
						.get(rightChild);

				HashMap<TreeNode, Integer> leftNbParasites = nbParasitesSubTree
						.get(leftChild);
				HashMap<TreeNode, Integer> rightNbParasites = nbParasitesSubTree
						.get(rightChild);

				List<TreeNode> associatedParasite = new ArrayList<TreeNode>();
				HashMap<TreeNode, Integer> nbParasites = new HashMap<TreeNode, Integer>();
				

				for (TreeNode nodeAssoc : leftChildAssoc) {
					if (!associatedParasite.contains(nodeAssoc))
						associatedParasite.add(nodeAssoc);
						nbParasites.put(nodeAssoc,
								leftNbParasites.get(nodeAssoc));
				}
				
				for (TreeNode nodeAssoc : rightChildAssoc) {
					if (!associatedParasite.contains(nodeAssoc))
						associatedParasite.add(nodeAssoc);
						nbParasites.put(nodeAssoc,
								rightNbParasites.get(nodeAssoc));
				}

				/*control if children have the same parsite*/
				for (TreeNode nodeAssoc : rightChildAssoc) {
					if (leftChildAssoc.contains(nodeAssoc)) {
						//associatedParasite.add(nodeAssoc);
						nbParasites.put(nodeAssoc,
								leftNbParasites.get(nodeAssoc)
										+ rightNbParasites.get(nodeAssoc)/*+1*/ );
					}

				}
				
				mapVerticalSpreadAssociations.put(node, associatedParasite);
				nbParasitesSubTree.put(node, nbParasites);
				
				
				/*remove duplicate values from mapVerticalSpreadAssociations*/
				HashMap<List<TreeNode>, TreeNode> helpMap = new HashMap<List<TreeNode>, TreeNode>();
				Set<TreeNode> keys = mapVerticalSpreadAssociations.keySet(); 
				Iterator<TreeNode> keyIter = keys.iterator();

				while (keyIter.hasNext()) {
				    TreeNode key = keyIter.next();
				    List<TreeNode> value = mapVerticalSpreadAssociations.get(key);
				    helpMap.put(value, key);
				}
				
				//mapVerticalSpreadAssociations.clear();
				Set<List<TreeNode>> keys2 = helpMap.keySet(); 
				Iterator<List<TreeNode>> keyIter2 = keys2.iterator();
				while (keyIter2.hasNext()) {
					List<TreeNode> key = keyIter2.next();
				    TreeNode value = helpMap.get(key);
				    mapVerticalSpreadAssociations.put(value, key);	
				}
				/*end for remove duplicate values from mapVerticalSpreadAssociations */
				
				
			} else {
				addVerticalSpread(leftChild);
				
			}
		}
		
	}
	
	
	/**
	 * Creates the VERTICAL SPREAD table of the tree. Notice that, only the nodes that are
	 * in the main tree are included in the table.
	 */
	public void createVERTICALSPREADTable(Tree parasiteTree) {
		verticalSpreadMap = new HashMap<TreeNode, Double>();
		//for each node of the host Tree
		
		
		TreeNode [] hostNode = this.getNodes();
		for (int i=0; i<hostNode.length; i++){
			double probab = 0;
			
			
				if(hostNode[i].isLeaf() || hostNode[i].isRoot()){
					if(hostNode[i].isLeaf()){
						probab = 1.00;
						verticalSpreadMap.put(hostNode[i], probab);
						hostNode[i].setPVerticalSpread(probab);
					}
					else{
						probab = 0;
						verticalSpreadMap.put(hostNode[i], probab);
						hostNode[i].setPVerticalSpread(probab);
					}
				}
				else{
					double constant = 1/((double) mapVerticalSpreadAssociations.get(hostNode[i]).size());
					//double constatSubTree = hostNode[i].getNbNodesSubTree();
					double constatSubTree = hostNode[i].getLeafSet().size()-1;
					Map<TreeNode, Integer> infectedP = nbParasitesSubTree.get(hostNode[i]);
					for(Entry<TreeNode, Integer> entry : infectedP.entrySet()){
						int nbInfectedP = entry.getValue()-1;
						probab += nbInfectedP/constatSubTree;
						
					}
					probab = probab * constant;
					hostNode[i].setPVerticalSpread(probab);

					verticalSpreadMap.put(hostNode[i], probab);
					setInitialVERTICALSPREADactivated();
				}
			
			
		}
		
	}
	
	
	/**
	 * Set the default VERTICAL SPREAD Activation. TRUE if the probability is 1, FALSE otherwise
	 */
	public void setInitialVERTICALSPREADactivated(){
		for (Map.Entry<TreeNode, Double> entry : verticalSpreadMap.entrySet())
		{
		    if (entry.getValue()==1.0)
		    	entry.getKey().setVerticalSpreadActivated(true);
		    else
		    	entry.getKey().setVerticalSpreadActivated(false);
		}
	}
	
	
	
	/**
     * Active the Vertical Spread Activation for all the nodes of a subtree.
     * 
     * @param subTreeRoot
     * 		  Root of the subtree that is going to be traversed
     */
    public void setVerticalSpreadActivated_SubTree(TreeNode subTreeRoot)
    {
    	if (!subTreeRoot.isLeaf()){
    		subTreeRoot.setVerticalSpreadActivated(true);
    		setVerticalSpreadActivated_SubTree(subTreeRoot.getChildren()[0]);
    		setVerticalSpreadActivated_SubTree(subTreeRoot.getChildren()[1]);
    	}
    }

	/**
	 * Return an array containing all leaf nodes of the tree.
	 * 
	 * @return An array containing all leaf nodes of the tree.
	 */
	public TreeNode[] getLeafNodes() {
		if (leaves == null)
			updateAttributes();
		return leaves.toArray(new TreeNode[0]);
	}

	
	/**
	 * Return the ArrayList containing all leaf nodes of the tree.
	 * 
	 * @return The ArrayList containing all leaf nodes of the tree.
	 */
	public ArrayList<TreeNode> getLeafNodesArrayList(){
		return leaves;
	}
	
	
	/**
	 * Return an array containing all leaf nodes of the tree.
	 * 
	 * @return An array containing all leaf nodes of the tree.
	 */
	public TreeNode[] getNodes() {
		if (nodes == null)
			updateAttributes();
		return nodes.toArray(new TreeNode[0]);
	}

	/**
	 * Creates a map between each partition of the tree and their number of
	 * occurrences.
	 * 
	 * @param leavesOrder
	 *            An ordering of the leaves.
	 * @return A map between each partition of the tree and their number of
	 *         occurrences.
	 */
	private HashMap<ArrayList<Integer>, Integer> partitionMap(
			HashMap<String, Integer> leavesOrder) {
		HashMap<ArrayList<Integer>, Integer> partitionMap = new HashMap<ArrayList<Integer>, Integer>();
		if (nodes == null)
			updateAttributes();
		for (TreeNode n : nodes) {
			ArrayList<Integer> p = n.partition(leavesOrder);
			if (partitionMap.containsKey(p)) {
				int i = partitionMap.get(p);
				partitionMap.put(p, i + 1);
			} else {
				partitionMap.put(n.partition(leavesOrder), 1);
			}
		}
		return partitionMap;
	}

	/**
	 * Computes the (modifed)RF distance between this tree and the one given as
	 * argument using the order specified in the map.
	 * 
	 * @param other
	 *            Tree that will be compared to this tree.
	 * @param leavesOrder
	 *            An ordering of the leaves of this tree.
	 * @return The (modifed)RF distance between this tree and the one given as
	 *         argument.
	 */
	public int mRFDistance(Tree other, HashMap<String, Integer> leavesOrder) {
		int distance = 0;
		HashMap<ArrayList<Integer>, Integer> partitionMap = other
				.partitionMap(leavesOrder);
		for (TreeNode n : nodes) {
			ArrayList<Integer> partition = n.partition(leavesOrder);
			if (partitionMap.containsKey(partition)) {
				int i = partitionMap.get(partition);
				if (i == 0) {
					partitionMap.remove(partition);
					distance++;
				} else {
					partitionMap.put(partition, i - 1);
				}
			} else {
				distance++;
			}
		}
		for (int i : partitionMap.values()) {
			distance = distance + i;
		}
		return distance;
	}

	/**
	 * Computes the MAAC value between this tree and the one given as argument.
	 * 
	 * @param other
	 *            Tree that will be compared to this tree.
	 * @return The MAAC value between this tree and the one given as argument.
	 */
	public HashMultiset<String> computeMaac(Tree other) {
		if (nodes == null)
			updateAttributes();
		
		ArrayList<TreeNode> otherNodes = new ArrayList<TreeNode>(
				Arrays.asList(other.getNodes()));

		ArrayList<ArrayList<HashMultiset<String>>> maacTable = new ArrayList<ArrayList<HashMultiset<String>>>();
		
		
		for (int i = 0; i < nodes.size(); i++) {
			//System.out.println("maacTable1 " + maacTable);
			TreeNode v1 = nodes.get(i);
			maacTable.add(i,
					new ArrayList<HashMultiset<String>>(otherNodes.size()));
			//System.out.println("maacTable2 " + maacTable);
			for (int j = 0; j < otherNodes.size(); j++) {
				TreeNode v2 = otherNodes.get(j);

				if (v1.isLeaf() || v2.isLeaf()) {

					HashMultiset<String> seti = HashMultiset.create();
					HashMultiset<String> setj = HashMultiset.create();
					/*System.out.println("----------------------------");
					if(v1.getLabel().size() != 1){
						System.out.println("v1 : " + v1.getLabel().toString());
					}
					if(v2.getLabel().size() != 1){
						System.out.println("v2 : " + v2.getLabel().toString());
					}*/
					
					
					for (TreeNode leaf : v1.getLeafSet())
					{
						//seti.add(leaf.getLabel());
						seti.addAll(leaf.getLabel());
						//System.out.println("v1 : " + seti);
					}
						
					
					
					for (TreeNode leaf : v2.getLeafSet())
					{
						//setj.add(leaf.getLabel());
						setj.addAll(leaf.getLabel());
						//System.out.println("v2 : " + setj);
					}

					HashMultiset<String> setij = HashMultiset.create(Multisets
							.intersection(seti, setj));
					maacTable.get(i).add(j, setij);
					//System.out.println(maacTable.toString());

				} else {

					HashMultiset<String> match = match(v1, v2, maacTable,
							nodes, otherNodes);
					HashMultiset<String> diag = diag(v1, v2, maacTable, nodes,
							otherNodes);
					if (match.size() > diag.size()) {
						maacTable.get(i).add(j, match);
					} else {
						maacTable.get(i).add(j, diag);
					}
				}
			}
		}
		return maacTable.get(nodes.size() - 1).get(otherNodes.size() - 1);
	}

	/**
	 * Compute the MAAC distance between this tree and the one given as
	 * argument.
	 * 
	 * @param other
	 *            Tree that will be compared to this tree.
	 * @return The MAAC distance between this tree and the one given as
	 *         argument.
	 */
	public double computeMaacDistance(Tree other) {
		if (leaves == null)
			updateAttributes();
		double maacValue = computeMaac(other).size();
		double maxSize = Math.max(leaves.size(), other.getLeafNodes().length);
		System.out.println("maacValue: " + maacValue + " maxSize: " + maxSize + " distance: " + (maacValue / maxSize) );
		return 1.0 - (maacValue / maxSize);
		//return 0;
	}

	/**
	 * Compute the MAAC similarity between this tree and the one given as
	 * argument.
	 * 
	 * @param other
	 *            Tree that will be compared to this tree.
	 * @return The MAAC similarity between this tree and the one given as
	 *         argument.
	 */
	public double computeMaacSimilarity(Tree other) {
		return 1 - computeMaacDistance(other);
	}

	/**
	 * Auxiliar function for computing the MAAC between two trees.
	 */
	private HashMultiset<String> match(TreeNode v1, TreeNode v2,
			ArrayList<ArrayList<HashMultiset<String>>> maacTable,
			ArrayList<TreeNode> po1, ArrayList<TreeNode> po2) {

		TreeNode c11 = v1.getChild(0);
		TreeNode c12 = v1.getChild(1);
		TreeNode c21 = v2.getChild(1);
		TreeNode c22 = v2.getChild(0);

		int oc11 = po1.indexOf(c11);
		int oc12 = po1.indexOf(c12);

		int oc21 = po2.indexOf(c21);
		int oc22 = po2.indexOf(c22);

		HashMultiset<String> a = maacTable.get(oc11).get(oc21);
		HashMultiset<String> c = maacTable.get(oc12).get(oc22);

		HashMultiset<String> b = maacTable.get(oc11).get(oc22);
		HashMultiset<String> d = maacTable.get(oc12).get(oc21);

		HashMultiset<String> match;
		if ((a.size() + c.size()) < (b.size() + d.size())) {
			match = HashMultiset.create(b.size() + d.size());
			match.addAll(b);
			match.addAll(d);
		} else {
			match = HashMultiset.create(a.size() + c.size());
			match.addAll(a);
			match.addAll(c);
		}

		return match;
	}
	
	
	/**
	 * Computes the MAMST value between this tree and the one given as argument.
	 * 
	 * @param other
	 *            Tree that will be compared to this tree.
	 * @return The MAMST value between this tree and the one given as argument.
	 */
	public HashMultiset<String> computeMamst(Tree other) {
		
		//creation matrix for calculate the isomorphism
		int [][] isomMatrix = new int [this.getNumberOfNodes()][other.getNumberOfNodes()];
		//creation isomorphism between leaves
		
		return null;
		}
		
	
	/**
	 * Calculate the isomorphism between two node.
	 * 
	 * @param node1
	 *            node of original parasite tree
	 * @param
	 * 			  node of pther parasite tree, that will be compared to this tree.
	 * @return The MAMST value between this tree and the one given as argument.
	 */
	private int calculateIsomorphism (TreeNode node1, TreeNode node2){
		//assign at label1 the smaller node in term of length 
		List<String> labels1;
		List<String> labels2;
		if(node1.getLabel().size() <= node2.getLabel().size()){
			labels1 = node1.getLabel();
			labels2 = node2.getLabel();
		}
		else{
			labels2 = node1.getLabel();
			labels1 = node2.getLabel();
		}
		Iterator<String> labels1Iterator = labels1.iterator();
		int isomorpValue = 0;
		while (labels1Iterator.hasNext()) {
			if(labels2.contains(labels1Iterator))
				isomorpValue++;
		}//while
		return isomorpValue;
	}
	
	


	/**
	 * Auxiliar function for computing the MAAC between two trees.
	 */
	private HashMultiset<String> diag(TreeNode v1, TreeNode v2,
			ArrayList<ArrayList<HashMultiset<String>>> maacTable,
			ArrayList<TreeNode> po1, ArrayList<TreeNode> po2) {

		int oc11 = po1.indexOf(v1.getChild(0));
		int oc12 = po1.indexOf(v1.getChild(1));
		int ov1 = po1.indexOf(v1);
		int oc21 = po2.indexOf(v2.getChild(0));
		int oc22 = po2.indexOf(v2.getChild(1));
		int ov2 = po2.indexOf(v2);
		HashMultiset<String> opt;

		HashMultiset<String> opt1 = maacTable.get(ov1).get(oc21);
		HashMultiset<String> opt2 = maacTable.get(ov1).get(oc22);
		HashMultiset<String> opt3 = maacTable.get(oc11).get(ov2);
		HashMultiset<String> opt4 = maacTable.get(oc12).get(ov2);

		if (opt1.size() < opt2.size())
			opt = HashMultiset.create(opt2);
		else
			opt = HashMultiset.create(opt1);

		if (opt.size() < opt3.size())
			opt = HashMultiset.create(opt3);

		if (opt.size() < opt4.size())
			opt = HashMultiset.create(opt4);

		return opt;
	}
	
	
	

	/**
	 * Prints the LCA table of this tree.
	 */
	public String printLCATable() {

		TreeNode[][] table = getLCATable();

		/* Prepare the header of the table */
		StringBuffer buffer = new StringBuffer();
		for (int i = 0; i < nodes.size(); i++)
			buffer.append("\t" + Integer.toString(nodes.get(i).getKey()));

		buffer.append("\n--");

		for (int i = 0; i < nodes.size(); i++)
			buffer.append("\t--");

		buffer.append("\n");

		/* Add one line for each node. */
		for (int i = 0; i < nodes.size(); i++) {
			buffer.append(Integer.toString(nodes.get(i).getKey()) + ":");

			for (int j = 0; j < nodes.size(); j++)
				buffer.append("\t" + Integer.toString(table[i][j].getKey()));

			buffer.append("\n");
		}

		return buffer.toString();

	}
	
	/**
	 * Prints the VerticalSpreadAssociations Map of this tree.
	 */
	public void printMapVerticalSpreadAssociations(){
		System.out.println(mapVerticalSpreadAssociations.toString());
	}
	
	
	/**
	 * Prints the NbParasitesSubTree Map of this tree.
	 */
	public void printNbParasitesSubTree(){
		System.out.println(nbParasitesSubTree.toString());
	}

	
	/**
	 * Prints the Vertical Spread Map of this tree.
	 */
	public void printVerticalSpreadMap(){
		System.out.println(verticalSpreadMap.toString());
	}
	
	
	/**
	 * Prints the Horizontal Spread matrix of this tree.
	 */
	public void printHorizontalSpreadMatrix() {
		for (int i=0; i<horizontalSpreadMatrix.length; i++){
			for (int j=0; j<horizontalSpreadMatrix.length; j++){
				System.out.print(horizontalSpreadMatrix[i][j] + " ");
			}
			System.out.println();
		}
	}
	
	
	
	/**
	 * Returns the mapping between labels and nodes. WARNING!!! This method
	 * assumes that the labels are unique.
	 * 
	 * @return The mapping between labels and nodes.
	 */
	public HashMap<String, TreeNode> getMapLabelToNode() {
		if (mapLabelToNode == null) {

			if (nodes == null)
				updateAttributes();

			mapLabelToNode = new HashMap<String, TreeNode>();
			for (TreeNode node : nodes) {
				if (mapLabelToNode.containsKey(node.getLabel()))
					throw new RuntimeException(
							"This method should not be applied for tree which have non-unique labels.");
				//mapLabelToNode.put(node.getLabel(), node);
				mapLabelToNode.put(node.getLabel().toString(), node);
			}
		}
		return mapLabelToNode;
	}

	
	
	/**
	 * Returns all nodes that are in the path from the descendant to the
	 * ancestor.
	 * 
	 * @param descendant
	 *            Descendant node.
	 * @param ancestor
	 *            Ancestor node.
	 * @return All nodes that are in the path from the descendant to the
	 *         ancestor.
	 */
	public ArrayList<TreeNode> getPathFromDescendantToAncestor(int descendant,
			int ancestor) {

		TreeNode descendantNode = mapKeyToNode.get(descendant);
		TreeNode ancestorNode = mapKeyToNode.get(ancestor);

		ArrayList<TreeNode> toReturn = new ArrayList<TreeNode>();
		while (descendantNode != null && descendantNode != ancestorNode) {
			toReturn.add(descendantNode);
			descendantNode = descendantNode.getParent();
		}
		toReturn.add(ancestorNode);

		if (descendantNode == null)
			return null;

		return toReturn;
	}
	
	public  HashMap<TreeNode, List<TreeNode>> getMapVerticalSpreadAssociations(){
		return mapVerticalSpreadAssociations;
	}
	
	public HashMap<TreeNode, HashMap<TreeNode, Integer>> getNbParasitesSubTree(){
		return nbParasitesSubTree;
	}

	

}
