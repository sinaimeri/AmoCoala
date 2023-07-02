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

package  graph;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import  trees.Tree;
import  trees.TreeNode;

public class Graph {

    private static final int GRAY = 1;

    private static final int BLACK = 2;

    /** Set of edges of this graph object. */
    private HashMap<Integer, HashSet<Integer>> edges;

    /** Set of nodes of this graph object. */
    private HashSet<Integer> nodes;

    /** Set of temporary edges whose insertion on the graph can be commited or rolled back. */
    private HashSet<Edge> temporaryEdges;

    /** Set of temporary nodes whose insertion on the graph can be commited or rolled back. */
    private HashSet<Integer> temporaryNodes;

    /** HashMap for the cyclicity verification algorithm. */
    private HashMap<Integer, Integer> status;


    /**
     * Constructor which builds an empty graph.
     */
    public Graph() {
        nodes = new HashSet<Integer>();
        edges = new HashMap<Integer, HashSet<Integer>>();
        temporaryNodes = new HashSet<Integer>();
        temporaryEdges = new HashSet<Edge>();
    }


    /**
     * Constructor which builds a graph based on the given tree.
     * 
     * @param tree
     *        Tree that will be used in the construction of the graph.
     */
    public Graph(Tree tree) {

        nodes = new HashSet<Integer>();
        edges = new HashMap<Integer, HashSet<Integer>>();
        temporaryNodes = new HashSet<Integer>();
        temporaryEdges = new HashSet<Edge>();

        /* First add all nodes of the tree. */
        for (TreeNode n: tree.getNodes()) {
            this.addNode(n.getKey());
        }

        /* Then, add all edges. */
        for (TreeNode n: tree.getNodes()) {
            TreeNode c0 = n.getChild(0);
            TreeNode c1 = n.getChild(1);
            if (c0 != null) {
                Edge edge = new Edge(n.getKey(), c0.getKey());
                this.addEdge(edge);
            }
            if (c1 != null) {
                Edge edge = new Edge(n.getKey(), c1.getKey());
                this.addEdge(edge);
            }
        }

        /* Finally, commit the insertion. */
        commit();
    }


    /**
     * Adds the node to the graph, if it does not already exist. To commit (or rollback) the
     * addition of this node, the method commit() (rollback()) should be called.
     * 
     * @param node
     *        Node that is going to be inserted into the graph. Notice that, if the node already
     *        belongs to the graph, we do not add it (it does not receive the "temporary" status
     *        and, in a case of a call of the method rollback(), it will not be removed from the
     *        graph).
     */
    public void addNode(int node) {
        if (!nodes.contains(node) && !temporaryNodes.contains(node)) {
            /* It's a new node, we can add it in the list of temporary nodes. */
            temporaryNodes.add(node);
        }
    }


    /**
     * Adds the edge to the graph, if it does not already exist. To commit (or rollback) the
     * addition of this edge, the method commit() (rollback()) should be called.
     * 
     * @param edge
     *        Edge that is going to be inserted into the graph. Notice that, if the edge already
     *        belongs to the graph, we do not add it (it does not receive the "temporary" status
     *        and, in a case of a call of the method rollback(), it will not be removed from the
     *        graph).
     */
    public void addEdge(Edge edge) {

        int source = edge.getSource();
        int target = edge.getTarget();

        if ((nodes.contains(source) || temporaryNodes.contains(source))
            && (nodes.contains(target) || temporaryNodes.contains(target))) {
            /* We can add the edge, source and target nodes are registered into the graph. */

            if (edges.containsKey(source)) {
                HashSet<Integer> targets = edges.get(source);
                if (!targets.contains(target)) {
                    /*
                     * The edge does not exist, add it on the graph and on the list of temporary
                     * edges.
                     */
                    targets.add(target);
                    edges.put(source, targets);
                    temporaryEdges.add(edge);
                }
            } else {
                /* The edge does not exist, add it on the graph and on the list of temporary edges. */
                HashSet<Integer> targets = new HashSet<Integer>();
                targets.add(target);
                edges.put(source, targets);
                temporaryEdges.add(edge);
            }

        }

    }


    /**
     * Commit the insertion of nodes and edges (nodes and edges are going to lose the status
     * "temporary").
     */
    public void commit() {
        nodes.addAll(temporaryNodes);
        temporaryEdges.clear();
        temporaryNodes.clear();
    }


    /**
     * This method erases all insertions (nodes and edges) which were performed after the last call
     * of the confirm method.
     */
    public void roolback() {
        /* First cancel the edges. */
        for (Edge edge: temporaryEdges) {
            int source = edge.getSource();
            int target = edge.getTarget();
            HashSet<Integer> targets = edges.get(source);
            targets.remove(target);
            if (targets.size() == 0)
                edges.remove(source);
        }
        temporaryEdges.clear();
        /* Then, cancel the nodes. */
        temporaryNodes.clear();
    }


    /**
     * Check if the graph has cycles.
     * 
     * @return <code>TRUE</code> if the graph has cycles. Otherwise, returns <code>FALSE</code>.
     */
    public boolean isCyclic() {
        status = new HashMap<Integer, Integer>();
        for (Integer node: nodes) {
            if (!status.containsKey(node)) {
                if (isCyclic(node)) {
                    return true;
                }
            }
        }
        return false;
    }


    /**
     * Recursive function for checking graph cyclicity.
     * 
     * @param node
     *        Node to be verified.
     * @return <code>TRUE</code> if a cycle is found. Otherwise, returns <code>FALSE</code>.
     */
    private boolean isCyclic(int node) {
        if (!status.containsKey(node)) {
            /* The node is being visited for the first time. */
            /* Mark it as gray. */
            status.put(node, GRAY);
            /* Visit its targets */
            HashSet<Integer> targets = edges.get(node);
            if (targets != null)
                for (Integer target: targets) {
                    if (isCyclic(target)) {
                        return true;
                    }
                }
            /* We alredy visited all children. Mark it as black. */
            status.put(node, BLACK);
        } else if (status.get(node) == GRAY) {
            /* The node is grey. We detected a cycle. */
            return true;
        }
        return false;
    }


    /**
     * Returns the list of all nodes (confirmed and temporary).
     * 
     * @return The list of all nodes (confirmed and temporary).
     */
    public Integer[] getAllNodes() {
        ArrayList<Integer> all = new ArrayList<Integer>(nodes);
        all.addAll(temporaryNodes);
        return all.toArray(new Integer[0]);
    }


    /**
     * Returns true if the node belongs to the vertex set of the graph.
     * 
     * @param node
     *        Node to be verified for existence.
     * @return True if the node belongs to the vertex set of the graph and false, otherwise.
     */
    public boolean isInTheGraph(int node) {
        return nodes.contains(node) || temporaryNodes.contains(node);
    }


    /**
     * Returns a string represetation of this graph.
     * 
     * @return A string represetation of this graph.
     */
    public String toString() {
        StringBuffer buffer = new StringBuffer();
        buffer.append("Nodes: " + nodes.toString() + "\n");
        buffer.append("Temporary nodes: " + temporaryNodes.toString() + "\n");
        buffer.append("Edges:\n");
        for (Integer source: edges.keySet()) {
            HashSet<Integer> targets = edges.get(source);
            for (Integer target: targets) {
                buffer.append("(" + source + "," + target + ")\n");
            }
        }
        buffer.append("Temporary edges: " + temporaryEdges.toString() + "\n");
        return buffer.toString();
    }

}
