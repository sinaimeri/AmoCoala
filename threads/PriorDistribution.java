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

import java.text.NumberFormat;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Random;

import jsat.distributions.multivariate.Dirichlet;
import jsat.linear.DenseVector;
import jsat.linear.Vec;
import  generator.IParasiteGenerator;
import  generator.IScenario;
import  generator.defaultmodel.DefaultModel;
import  trees.Tree;
import  util.Statistics;

public class PriorDistribution extends GenerateDistribution {

	
	/**
     * Random generator.
     */
    private Random random;
    
    
	/** Alpha value for cospeciation (Dirichlet process). */
	protected double alphaCospeciation;

	/** Alpha value for cospeciation (Dirichlet process). */
	protected double alphaDuplication;

	/** Alpha value for cospeciation (Dirichlet process). */
	protected double alphaLoss;

	/** Alpha value for cospeciation (Dirichlet process). */
	protected double alphaSwitch;

	/**
	 * Thread to generate vectors for the prior distribution.
	 * 
	 * @param vectorManager
	 *            Object that manage the access to the list of vectors in a
	 *            synchronized way.
	 * @param numberOfVectors
	 *            Number of vectors to generate.
	 * @param numberOfTrees
	 *            Number of trees that must be generated for each vector.
	 * @param hostTree
	 *            Host tree.
	 * @param parasiteTree
	 *            Parasite tree.
	 * @param mappingParasiteHost
	 *            Mapping of the parasite leaf nodes into the host leaf nodes.
	 * @param model
	 *            Model to be used during the simulation.
	 * @param metric
	 *            Metric to be used during the simulation.
	 * @param maximumNumberOfTreesFactor
	 *            Maximum number of trees factor.
	 * @param cyclicityTestModel
	 *            Cyclicity test model identifier.
	 * @param alphaCospeciation
	 *            Alpha cospeciation (Dirichlet process).
	 * @param alphaDuplication
	 *            Alpha duplication (Dirichlet process).
	 * @param alphaLoss
	 *            Alpha loss (Dirichlet process).
	 * @param alphaSwitch
	 *            Alpha switch (Dirichlet process).
	 * @param alpha1
	 *            Alpha 1 value for computing LEAVES AND MAAC or EVENTS AND MAAC
	 *            metric.
	 * @param alpha2
	 *            Alpha 2 value for computing LEAVES AND MAAC or EVENTS AND MAAC
	 *            metric.
	 * @param rootMappingProbability
	 *            Probability of mapping into the root node.
	 */
	public PriorDistribution(SynchronizedVectorManager vectorManager,
			int numberOfVectors, int numberOfTrees, Tree hostTree,
			//Tree parasiteTree, HashMap<String, String> mappingParasiteHost,
			Tree parasiteTree, HashMap<String, List<String>> mappingParasiteHost,
			int model, int metric, int maximumNumberOfTreesFactor,
			int cyclicityTestModel, double alphaCospeciation,
			double alphaDuplication, double alphaLoss, double alphaSwitch,
			double alpha1, double alpha2, double rootMappingProbability) {
		super(vectorManager, numberOfVectors, numberOfTrees, hostTree,
				parasiteTree, mappingParasiteHost, model, metric,
				maximumNumberOfTreesFactor, cyclicityTestModel, alpha1, alpha2,
				rootMappingProbability);
		this.alphaCospeciation = alphaCospeciation;
		this.alphaDuplication = alphaDuplication;
		this.alphaLoss = alphaLoss;
		this.alphaSwitch = alphaSwitch;
	}

	public void run() {
		if (model == GenerateDistribution.DEFAULT_MODEL) {
			runDefaultModel();
		} else {
			//runCoalescentModel();
		}
	}

	private void runDefaultModel() {
		NumberFormat format = NumberFormat.getInstance();
		format = NumberFormat.getInstance(Locale.US);
		format.setMinimumFractionDigits(4);
		format.setMaximumFractionDigits(4);

		double[] vector = { alphaCospeciation, alphaDuplication, alphaSwitch,
				alphaLoss };

		this.random = new Random();
		DenseVector denseVector = new DenseVector(vector);
		Dirichlet dirichlet = new Dirichlet(denseVector);
		

		while (vectorManager.getNumberOfVectors() < numberOfVectors) {

			Vec sampledVector = dirichlet.sample(1, random).get(0);
			
			
			updateExpectedNumberOfEvents(sampledVector.get(0),
					sampledVector.get(1), sampledVector.get(2),
					sampledVector.get(3));

			IParasiteGenerator generator = new DefaultModel(hostTree, parasiteTree,
					sampledVector.get(0), sampledVector.get(1),
					sampledVector.get(2), parasiteTree.getNumberOfLeafNodes(),
					maximumNumberOfLeaves, cyclicityTestModel,
					MAXIMUM_SIZE_STOP_CRITERIUM, rootMappingProbability);

			int nScenarios = 0;

			double[] values = new double[numberOfTrees];
			double[] freqVerticalSpread = new double[numberOfTrees];
			double[] freqHorizontalSpread = new double[numberOfTrees];
			int[] nbMultipleAssociations = new int [numberOfTrees];
			int[] dimensionTrees = new int [numberOfTrees];
			
			Tree[] paraTree = new Tree[numberOfTrees];
            int[] sizeTree = new int[numberOfTrees];
            String[] mappi = new String[numberOfTrees];
			
			/*variable for print all data concerning jumps
			 * 
			 * double[] maxJump = new double[numberOfTrees];
			 * double[] meanJump = new double[numberOfTrees];
			 * double[] nbJump = new double[numberOfTrees]; 
			 * int[] nbHostSwitch = new int[numberOfTrees];
			 */
			//Double[] finale = null;
			//Double[] first = finale;
			for (int j = 0; nScenarios < numberOfTrees
					&& j < maximumNumberOfTrees; j++) {
				IScenario scenario = generator.generateParasiteTree();
				if (scenario != null) {
					double distance = 1.0;
					scenario.getParasiteTree().updateNumberOfMultipleAssociations();
					switch (metric) {
					case METRIC_MAAC_DISTANCE:
						distance = computeMAACDistanceFromRealParasite(scenario);
						break;
					case METRIC_LEAVES_AND_MAAC_DISTANCE:
						distance = computeLEAVES_AND_MAAC_Metric(scenario);
						break;
					case METRIC_EVENTS_AND_MAAC_DISTANCE:
						distance = computeEVENTS_AND_MAAC_Metric(scenario);
						break;
					case METRIC_EVENTS_AND_MAAC_DISTANCE_MULTIPLEASSOCIATIONS:
						distance = computeEVENTS_AND_MAAC_Metric_MultipleAssociations(scenario, parasiteTree.getNumberOfMultipleAssociations());
						break;
					}
					
					freqVerticalSpread[nScenarios] = scenario.getNumberOfVerticalSpread();
					freqHorizontalSpread[nScenarios] = scenario.getNumberOfHorizontalSpread();
					
					paraTree[nScenarios] = scenario.getParasiteTree();
                    sizeTree[nScenarios] = paraTree[nScenarios].getNumberOfNodes();
                    mappi[nScenarios] = scenario.getMapping().toString();
					
					
					int nb_multipleAssociations = 0;
					for (int u=0; u<scenario.getParasiteTree().getNumberOfLeafNodes(); u++){
						nb_multipleAssociations = nb_multipleAssociations + (scenario.getParasiteTree().getLeafNodes()[u].getLabel().size()-1);
					}
					dimensionTrees[nScenarios] = (2* scenario.getNumberOfLeaves())-1; 
					nbMultipleAssociations[nScenarios] = nb_multipleAssociations;
					values[nScenarios++] = distance;
					

				}	
				
			}

			
			double representativeDistance = 1.0;
			double representativeVerticalSpread = 0.0;
			double representativeHorizontalSpread = 0.0;
			//double mediana = 0;
			/*assign default values 
			 * 
			 *double maxMaxJumpTotal =-1;
			 *double maxMeanJumpTotal = -1;
			 *double meanMaxJumpTotal = -1;
			 *double meanMeanJumpTotal = -1;
			 *double nbMaxJumpTotal = -1;
			 *double nbMeanJumpTotal = -1;
			 */
			
			if (nScenarios == numberOfTrees) {
				representativeVerticalSpread = Statistics.mean(freqVerticalSpread);
				representativeHorizontalSpread = Statistics.mean(freqHorizontalSpread);
				representativeDistance = Statistics.mean(values);
				//System.out.println(Statistics.mean(dimensionTrees) + "\t" + Statistics.mean(nbMultipleAssociations));
				
				
			}
				
			double[] array = { sampledVector.get(0), sampledVector.get(1),
					sampledVector.get(2), sampledVector.get(3), representativeVerticalSpread, representativeHorizontalSpread, representativeDistance
					};
			vectorManager.registerVector(array);

		}

	}

}
