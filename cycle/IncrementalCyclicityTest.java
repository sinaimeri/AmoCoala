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
#                Marie-France Sagot    									   #
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
package  cycle;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Random;

import  graph.Edge;
import  graph.Graph;
import  trees.Tree;
import  trees.TreeNode;

public class IncrementalCyclicityTest implements CyclicityTest {

    /** Model identifier. */
    private int model;

    /** Graph. */
    private Graph graph;

    /** Host tree */
    private Tree hostTree;

    /** Parasite tree */
    private Tree parasiteTree;

    /** List of confirmed transfer edges (edges on the parasite tree). */
    private ArrayList<Edge> confirmedTransferEdges;

    /** Hashmap which keeps the mapping parasite node to host node (key values). */
    private HashMap<Integer, List<Integer>> mapping;

    /** HashMap which keeps for each host tree node, the list of incomparable nodes. */
    public HashMap<TreeNode, ArrayList<TreeNode>> incomparableNodesHostTree;

    /** HashMap which keeps for each host tree node, its list of descendants. */
    private HashMap<TreeNode, HashSet<TreeNode>> descendantsHostTree;

    /** HashMap which keeps for each host tree node, its list of ancestors. */
    private HashMap<TreeNode, HashSet<TreeNode>> ancestorsHostTree;


    /**
     * Creates a cyclicity test for the given model.
     * 
     * @param model
     *        Model identifier.
     * @param graph
     *        Graph that is going to be tested for cyclicity.
     */
    public IncrementalCyclicityTest(int model,
                                    Tree hostTree,
                                    Tree parasiteTree,
                                    HashMap<Integer, List<Integer>> hashMap) {

        switch (model) {
            case STOLZER:
                /* Original version : We start with an empty graph. */
                graph = new Graph();
                break;
            case TRANSFER_EDGES:
            case TOFIGH:
            case DONATI:
                /* We start with a graph with the host tree topology. */
                graph = new Graph(hostTree);
                break;
            default:
                throw new RuntimeException("Unexpected model identifier");
        }

        this.model = model;
        this.hostTree = hostTree;
        this.parasiteTree = parasiteTree;
        this.mapping = hashMap;
        this.confirmedTransferEdges = new ArrayList<Edge>();

        if (model == STOLZER) {

            /* Pre-processing for speed-up. */
            this.incomparableNodesHostTree = new HashMap<TreeNode, ArrayList<TreeNode>>();
            this.descendantsHostTree = new HashMap<TreeNode, HashSet<TreeNode>>();
            this.ancestorsHostTree = new HashMap<TreeNode, HashSet<TreeNode>>();

            for (TreeNode node: hostTree.getNodes()) {
                descendantsHostTree.put(node, hostTree.getDescendants(node.getKey(), true));
                ancestorsHostTree.put(node, hostTree.getAncestors(node.getKey(), true));
                ArrayList<TreeNode> incomparables = new ArrayList<TreeNode>();
                incomparableNodesHostTree.put(node, incomparables);
                for (TreeNode other: hostTree.getNodes()) {
                    TreeNode lca = hostTree.getLCA(node.getKey(), other.getKey());
                    if (lca != node && lca != other) {
                        incomparables.add(other);
                    }
                }
            }

        } else {

            /* Pre-processing for speed-up. */
            this.incomparableNodesHostTree = new HashMap<TreeNode, ArrayList<TreeNode>>();
            for (TreeNode node: hostTree.getNodes()) {
                ArrayList<TreeNode> incomparables = new ArrayList<TreeNode>();
                incomparableNodesHostTree.put(node, incomparables);
                for (TreeNode other: hostTree.getNodes()) {
                    TreeNode lca = hostTree.getLCA(node.getKey(), other.getKey());
                    if (lca != node && lca != other) {
                        incomparables.add(other);
                    }
                }
            }

        }
    }


    public Graph getGraph() {
        return graph;
    }


    /**
     * Selects a random incomparable node which allow the mapping of the parasite node into the host
     * node without creating cycles. Returns null if there is no available option.
     * 
     * @param parasiteNode
     *        The parasite node.
     * @param hostNode
     *        The host node.
     * @return A random incomparable node which allow the mapping of the parasite node into the host
     *         node without creating cycles. Returns null if there is no available option.
     */
    public TreeNode randomIncomparable(TreeNode parasiteNode, TreeNode hostNode) {

        ArrayList<TreeNode> candidates = null;

        if (confirmedTransferEdges.size() == 0) {
            /* Just to accelerate, if it is the first transfer edge, there is no need to test. */
            candidates = incomparableNodesHostTree.get(hostNode);
        } else {
            /*
             * For each possible incomparable node, we have to test if the effect of adding the
             * transfer edge generates a cycle.
             */
            candidates = new ArrayList<TreeNode>();
            for (TreeNode other: incomparableNodesHostTree.get(hostNode)) {
                if (testTransferEdge(parasiteNode, hostNode, other))
                    candidates.add(other);
            }
        }

        if (candidates.size() == 0)
            return null;

        return candidates.get((int)(Math.random() * candidates.size()));
    }


    /**
     * Returns true if the trasnfer edge passes in the cyclicity test (i.e., it does not create a
     * cycle) and false, otherwise.
     * 
     * @param u
     *        Parasite node u.
     * @param gammaU
     *        Host node where the parasite node u is mapped into.
     * @param gammaV
     *        Host node where the child of u will be mapped into.
     * @return True if the trasnfer edge passes in the cyclicity test (i.e., it does not create a
     *         cycle) and false, otherwise.
     */
    private boolean testTransferEdge(TreeNode u, TreeNode gammaU, TreeNode gammaV) {

        switch (model) {
            case STOLZER:
                stolzer(u, gammaU, gammaV);
                break;
            case TRANSFER_EDGES:
                /* Just add the transfer edge, the host nodes are already there. */
                graph.addEdge(new Edge(gammaU.getKey(), gammaV.getKey()));
                break;
            case TOFIGH:
                tofigh(u, gammaU, gammaV);
                break;
            case DONATI:
                donati(u, gammaU, gammaV);
                break;
        }

        boolean isCyclic = graph.isCyclic();

        graph.roolback();

        return !isCyclic;
    }


    /**
     * Adds a confirmed transfer edge.
     * 
     * @param u
     *        Parasite node u.
     * @param v
     *        Parasite node v.
     * @param gammaU
     *        Host node where the parasite node u is mapped into.
     * @param gammaV
     *        Host node where the parasite node v is mapped into.
     */
    public void addTransferEdge(TreeNode u, TreeNode v, TreeNode gammaU, TreeNode gammaV) {

        switch (model) {
            case IncrementalCyclicityTest.STOLZER:
                stolzer(u, gammaU, gammaV);
                break;
            case IncrementalCyclicityTest.TRANSFER_EDGES:
                /* Just add the transfer edge, the host nodes are already there. */
                graph.addEdge(new Edge(gammaU.getKey(), gammaV.getKey()));
                break;
            case IncrementalCyclicityTest.TOFIGH:
                tofigh(u, gammaU, gammaV);
                break;
            case IncrementalCyclicityTest.DONATI:
                donati(u, gammaU, gammaV);
                break;
        }

        confirmedTransferEdges.add(new Edge(u.getKey(), v.getKey()));
        graph.commit();

    }


    /**
     * Adds edges to the graph as reflect of the transfer edge (u',v') according with the Tofigh's
     * criterium.
     * 
     * @param uPrime
     *        Source of the trasnfer edge into the parasite tree.
     * @param gammaUprime
     *        Source of the transfer edge into the host tree (gamma(u')).
     * @param gammaVprime
     *        Target of the transfer edge into the host tree (gamma(v')).
     */
    private void tofigh(TreeNode uPrime, TreeNode gammaUprime, TreeNode gammaVprime) {
        HashSet<TreeNode> uPrimeDescendants = parasiteTree.getDescendants(uPrime.getKey(), true);
        Random random = new Random();
        for (Edge edge: confirmedTransferEdges) {
            /* edge = (u,v) */
            //TreeNode gammaU = hostTree.getNode(mapping.get(edge.getSource()));
        	List<TreeNode> gammaUList = hostTree.getNode(mapping.get(edge.getSource()));
            
            int randomIndex = random.nextInt(gammaUList.size());
            TreeNode gammaU = gammaUList.get(randomIndex);
            
            HashSet<TreeNode> descendantsGammaU = parasiteTree.getDescendants(edge.getSource(),
                                                                              true);
            if (descendantsGammaU.contains(uPrime)) {
                /* Get parent(gamma(u)) */
                TreeNode parentGammaU = gammaU.getParent();
                /* Add the edge (parent(gamma(u)),v'). */
                graph.addEdge(new Edge(parentGammaU.getKey(), gammaVprime.getKey()));
            } else if (uPrimeDescendants.contains(gammaU)) {
                /* Get parent(gamma(u')) */
                TreeNode parentGammaUprime = gammaUprime.getParent();
                /* Get v */
                TreeNode gammaV = hostTree.getNode(edge.getTarget());
                /* Add the edge (parent(gamma(u')),v). */
                graph.addEdge(new Edge(parentGammaUprime.getKey(), gammaV.getKey()));
            }
        }
    }


    /**
     * Adds edges to the graph as reflect of the transfer edge (u,v) according with the Stolzer's
     * criterium.
     * 
     * @param u
     *        Source of the trasnfer edge into the parasite tree.
     * @param gammaU
     *        Source of the transfer edge into the host tree (gamma(u)).
     * @param gammaV
     *        Target of the transfer edge into the host tree (gamma(v)).
     */
    private void stolzer(TreeNode u, TreeNode gammaU, TreeNode gammaV) {

        /* Add the nodes gamma(u), parent(gamma(u)) into the graph. */
        TreeNode parentGammaU = gammaU.getParent();
        graph.addNode(gammaU.getKey());
        graph.addNode(parentGammaU.getKey());

        /* Add the nodes gamma(v), parent(gamma(v)) into the graph. */
        TreeNode parentGammaV = gammaV.getParent();
        graph.addNode(gammaV.getKey());
        graph.addNode(parentGammaV.getKey());

        /* Create a list of addedNodes to speed-up the condition 3. */
        ArrayList<TreeNode> addedNodes = new ArrayList<TreeNode>();
        addedNodes.add(gammaU);
        addedNodes.add(parentGammaU);
        addedNodes.add(gammaV);
        addedNodes.add(parentGammaV);

        /* Add the edge (parent(gamma(u)),gamma(u)) into the graph. */
        graph.addEdge(new Edge(parentGammaU.getKey(), gammaU.getKey()));
        /* Add the edge (parent(gamma(v)),gamma(v)) into the graph. */
        graph.addEdge(new Edge(parentGammaV.getKey(), gammaV.getKey()));

        Random random = new Random();
        /* Condition 1 :-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:- */
        Integer[] graphNodes = graph.getAllNodes();
        ArrayList<TreeNode> gammaUdescendants = new ArrayList<TreeNode>(descendantsHostTree.get(gammaU));
        ArrayList<TreeNode> parentGammaUdescendants = new ArrayList<TreeNode>(descendantsHostTree.get(parentGammaU));
        gammaUdescendants.remove(gammaU);
        parentGammaUdescendants.remove(parentGammaU);

        ArrayList<TreeNode> gammaVdescendants = new ArrayList<TreeNode>(descendantsHostTree.get(gammaV));
        ArrayList<TreeNode> parentGammaVdescendants = new ArrayList<TreeNode>(descendantsHostTree.get(parentGammaV));
        gammaVdescendants.remove(gammaV);
        parentGammaVdescendants.remove(parentGammaV);

        for (int node: graphNodes) {

            TreeNode hostTreeNode = hostTree.getNode(node);
            ArrayList<TreeNode> nodeDescendants = new ArrayList<TreeNode>(descendantsHostTree.get(hostTreeNode));
            nodeDescendants.remove(hostTreeNode);

            if (nodeDescendants.contains(gammaU))
                graph.addEdge(new Edge(node, gammaU.getKey()));

            if (nodeDescendants.contains(gammaV))
                graph.addEdge(new Edge(node, gammaV.getKey()));

            if (nodeDescendants.contains(parentGammaU))
                graph.addEdge(new Edge(node, parentGammaU.getKey()));

            if (nodeDescendants.contains(parentGammaV))
                graph.addEdge(new Edge(node, parentGammaV.getKey()));

            if (gammaUdescendants.contains(hostTreeNode))
                graph.addEdge(new Edge(gammaU.getKey(), node));

            if (gammaVdescendants.contains(hostTreeNode))
                graph.addEdge(new Edge(gammaV.getKey(), node));

            if (parentGammaUdescendants.contains(hostTreeNode))
                graph.addEdge(new Edge(parentGammaU.getKey(), node));

            if (parentGammaVdescendants.contains(hostTreeNode))
                graph.addEdge(new Edge(parentGammaV.getKey(), node));

        }

        /* Condition 2 :-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:- */

        /* Confirmed Transfer Edge = (u,v) => (d,r) */
        /* New transfer edge = (u',v') => (d',r') */
        int d = 0;
        int r = 0;
        int dParent = 0;
        int rParent = 0;
        int dPrime = gammaU.getKey();
        int rPrime = gammaV.getKey();
        for (Edge confirmedTransferEdge: confirmedTransferEdges) {
            HashSet<TreeNode> descendants = parasiteTree.getDescendants(confirmedTransferEdge.getSource(),
                                                                        true);
            if (descendants.contains(u)) {
                //d = mapping.get(confirmedTransferEdge.getSource());
                //r = mapping.get(confirmedTransferEdge.getTarget());
                
                List<Integer> dList = mapping.get(confirmedTransferEdge.getSource());
                List<Integer> rList = mapping.get(confirmedTransferEdge.getTarget());
                
                int randomIndex = random.nextInt(dList.size());
                d = dList.get(randomIndex);
                randomIndex = random.nextInt(rList.size());
                r = rList.get(randomIndex);

                dParent = hostTree.getNode(d).getParent().getKey();
                rParent = hostTree.getNode(r).getParent().getKey();
                graph.addEdge(new Edge(dParent, dPrime));
                graph.addEdge(new Edge(dParent, rPrime));
                graph.addEdge(new Edge(rParent, dPrime));
                graph.addEdge(new Edge(rParent, rPrime));
            }
        }

        /*
         * This part is commented because, in the simulation, the only descendant of u is v. V is
         * being created now, so there is no other transfer edge whose source is descendant of u.
         */
        // /* New transfer edge = (u,v) => (d,r) */
        // /* Confirmed Transfer Edge = (u',v') => (d',r') */
        // d = gammaU.getKey();
        // r = gammaV.getKey();
        // dParent = gammaU.getParent().getKey();
        // rParent = gammaV.getParent().getKey();
        // dPrime = 0;
        // rPrime = 0;
        // ArrayList<TreeNode> uDescendants = parasiteTree.getDescendants(u.getKey());
        // for (Edge confirmedTransferEdge: confirmedTransferEdges) {
        // if (uDescendants.contains(parasiteTree.getNode(confirmedTransferEdge.getSource()))) {
        // dPrime = mapping.get(confirmedTransferEdge.getSource());
        // rPrime = mapping.get(confirmedTransferEdge.getTarget());
        // graph.addEdge(new Edge(dParent, dPrime));
        // graph.addEdge(new Edge(dParent, rPrime));
        // graph.addEdge(new Edge(rParent, dPrime));
        // graph.addEdge(new Edge(rParent, rPrime));
        // }
        // }

        /* Condition 3 :-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:- */
        for (Edge confirmedTransferEdge: confirmedTransferEdges) {

            //int gammaG = mapping.get(confirmedTransferEdge.getSource());
            //int gammaH = mapping.get(confirmedTransferEdge.getTarget());
            List<Integer> gammaGList = mapping.get(confirmedTransferEdge.getSource());
            List<Integer> gammaHList = mapping.get(confirmedTransferEdge.getTarget());
            
            int randomIndex = random.nextInt(gammaGList.size());
            int gammaG = gammaGList.get(randomIndex);
            randomIndex = random.nextInt(gammaHList.size());
            int gammaH = gammaHList.get(randomIndex);

            /* (g,h) => (gamma(g),gamma(h)) = (sk,sj) */
            TreeNode sk = hostTree.getNode(gammaG);
            ArrayList<TreeNode> ancestors = new ArrayList<TreeNode>(ancestorsHostTree.get(sk));
            ancestors.remove(sk);
            for (TreeNode si: addedNodes) {
                if (ancestors.contains(si)) {
                    graph.addEdge(new Edge(si.getKey(), gammaH));
                }
            }

            /* (g,h) => (gamma(g),gamma(h)) = (sj,sk) */
            sk = hostTree.getNode(gammaH);
            ancestors = new ArrayList<TreeNode>(ancestorsHostTree.get(sk));
            ancestors.remove(sk);
            for (TreeNode si: addedNodes) {
                if (ancestors.contains(si)) {
                    graph.addEdge(new Edge(si.getKey(), gammaG));
                }
            }

        }

        /* New transfer edge : (g,h) => (gamma(g),gamma(h)) = (gamma(u),gamma(v)) = (sk,sj) */
        ArrayList<TreeNode> ancestors = new ArrayList<TreeNode>(ancestorsHostTree.get(gammaU));
        ancestors.remove(gammaU);
        for (TreeNode si: ancestors) {
            if (graph.isInTheGraph(si.getKey()))
                graph.addEdge(new Edge(si.getKey(), gammaV.getKey()));
        }

        /* New transfer edge : (g,h) => (gamma(g),gamma(h)) = (gamma(u),gamma(v)) = (sj,sk) */
        ancestors = new ArrayList<TreeNode>(ancestorsHostTree.get(gammaV));
        ancestors.remove(gammaV);
        for (TreeNode si: ancestors) {
            if (graph.isInTheGraph(si.getKey()))
                graph.addEdge(new Edge(si.getKey(), gammaU.getKey()));
        }

    }


    /**
     * Adds edges to the graph as reflect of the transfer edge (u,v) according with the Stolzer's
     * criterium.
     * 
     * @param u
     *        Source of the trasnfer edge into the parasite tree.
     * @param gammaU
     *        Source of the transfer edge into the host tree (gamma(u)).
     * @param gammaV
     *        Target of the transfer edge into the host tree (gamma(v)).
     */
    private void donati(TreeNode u, TreeNode gammaU, TreeNode gammaV) {

        /* Get the parent of gamma(u) */
        TreeNode parentGammaU = gammaU.getParent();
        /* Get the parent of gamma(v) */
        TreeNode parentGammaV = gammaV.getParent();

        /* Add the edge (parent(gamma(u)),gamma(u')) = (parent(gamma(u)),gamma(u)) */
        /* (redundant) It is not necessary to add it */
        /* Add the edge (parent(gamma(u)),gamma(v')) = (parent(gamma(u)),gamma(v)) */
        graph.addEdge(new Edge(parentGammaU.getKey(), gammaV.getKey()));
        /* Add the edge (parent(gamma(v)),gamma(u')) = (parent(gamma(v)),gamma(u)) */
        graph.addEdge(new Edge(parentGammaV.getKey(), gammaU.getKey()));
        /* Add the edge (parent(gamma(v)),gamma(v')) = (parent(gamma(v)),gamma(v)) */
        /* (redundant) It is not necessary to add it */

        /* Now compare (u,v) with all the other confirmed transfer edges */
        Random random = new Random();
        for (Edge confirmedTransferEdge: confirmedTransferEdges) {

            /*
             * Here, we have to check only if u is the list of descendants of the sources of the
             * confirmed transfer edges. We do not need to verify if the source of the confirmed
             * transfer edges are in the list of descendants of u (in this list we have only v that
             * is a new node in the simulation, so there is no transfer edge that has source in v).
             */

            /*
             * In this case (u,v) = a confirmed tranfer edge and (u',v') is the tranfer edge
             * received by parameter.
             */
            HashSet<TreeNode> descendants = parasiteTree.getDescendants(confirmedTransferEdge.getSource(),
                                                                        true);
            if (descendants.contains(u)) {
            	//TreeNode testGammaU = hostTree.getNode(mapping.get(confirmedTransferEdge.getSource()));
            	//TreeNode testGammaV = hostTree.getNode(mapping.get(confirmedTransferEdge.getTarget()));
                List<TreeNode> testGammaUList = hostTree.getNode(mapping.get(confirmedTransferEdge.getSource()));
                List<TreeNode> testGammaVList = hostTree.getNode(mapping.get(confirmedTransferEdge.getTarget()));
                
                //System.out.println("testGammaUList: " + testGammaUList);
                //System.out.println("testGammaVList: " + testGammaVList);
                
            	int randomIndex = random.nextInt(testGammaUList.size());
                TreeNode testGammaU = testGammaUList.get(randomIndex);
                //System.out.println("testGammaU: " + testGammaU);
                randomIndex = random.nextInt(testGammaVList.size());
                TreeNode testGammaV = testGammaVList.get(randomIndex);
                //System.out.println("testGammaV: " + testGammaV);
                //System.out.println("\n \n ----------------------------");
                
            	

                TreeNode testParentGammaU = testGammaU.getParent();
                TreeNode testParentGammaV = testGammaV.getParent();

                /* Add the edge (parent(gamma(u)),gamma(u')) */
                graph.addEdge(new Edge(testParentGammaU.getKey(), gammaU.getKey()));
                /* Add the edge (parent(gamma(u)),gamma(v')) */
                graph.addEdge(new Edge(testParentGammaU.getKey(), gammaV.getKey()));
                /* Add the edge (parent(gamma(v)),gamma(u')) */
                graph.addEdge(new Edge(testParentGammaV.getKey(), gammaU.getKey()));
                /* Add the edge (parent(gamma(v)),gamma(v')) */
                graph.addEdge(new Edge(testParentGammaV.getKey(), gammaV.getKey()));

            }

        }

    }


    public ArrayList<Edge> confirmed() {
        return confirmedTransferEdges;
    }
}
