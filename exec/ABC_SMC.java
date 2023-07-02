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

package  exec;

import java.text.NumberFormat;
import java.util.Arrays;
import java.util.Locale;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.UnrecognizedOptionException;

import  bayesian.ABC_SMC_Process;
import  cycle.CyclicityTest;
import  threads.GenerateDistribution;
import  util.NexusFileParserException;

public class ABC_SMC {

    /* -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:- */
    public static void main(String[] args) {
        ABC_SMC_Process process = new ABC_SMC_Process();
        Options opts = createOptions();
        String programUsage = "java -jar Coala.jar -input <file> [opts]";
        try {
            if (processParameters(process, args, opts, programUsage)) {
                process.run();
            }
        } catch (UnrecognizedOptionException uoe) {
            System.out.println(uoe.getMessage());
            HelpFormatter f = new HelpFormatter();
            f.setWidth(120);
            f.printHelp(programUsage, opts);
        } catch (NexusFileParserException e) {
            System.out.println(e.getMessage());
            HelpFormatter f = new HelpFormatter();
            f.setWidth(120);
            f.printHelp(programUsage, opts);
        } catch (Exception e) {
            System.out.println("Unknown error:");
            e.printStackTrace();
            HelpFormatter f = new HelpFormatter();
            f.setWidth(120);
            f.printHelp(programUsage, opts);
        }
    }


    /* -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:- */
    @SuppressWarnings("static-access")
    private static Options createOptions() {

        NumberFormat format = NumberFormat.getInstance(Locale.US);
        format.setMinimumFractionDigits(4);
        format.setMaximumFractionDigits(4);

        String description = null;
        Options opts = new Options();

        /* -------------------------------------------------------------------------------------- */
        /* Alpha 1 for LEAVES AND MAAC or EVENTS AND MAAC metrics */
        description = "Constant alpha 1 for metrics LEAVES AND MAAC or EVENTS AND MAAC.\n"
                      + "Default value for LEAVES AND MAAC metric = "
                      + format.format(ABC_SMC_Process.DEFAULT_LEAVES_AND_MAAC_ALPHA1) + "\n"
                      + "Default value for EVENTS AND MAAC metric = "
                      + format.format(ABC_SMC_Process.DEFAULT_EVENTS_AND_MAAC_ALPHA1) + "\n-\n";
        Option a1 = OptionBuilder.withArgName("value").hasArgs(1).withValueSeparator().withDescription(description).create("a1");
        opts.addOption(a1);

        /* -------------------------------------------------------------------------------------- */
        /* Alpha 2 for LEAVES AND MAAC or EVENTS AND MAAC metrics */
        description = "Constant alpha 2 for metrics LEAVES AND MAAC or EVENTS AND MAAC.\n"
                      + "Default value for LEAVES AND MAAC metric = "
                      + format.format(ABC_SMC_Process.DEFAULT_LEAVES_AND_MAAC_ALPHA2) + "\n"
                      + "Default value for EVENTS AND MAAC metric = "
                      + format.format(ABC_SMC_Process.DEFAULT_EVENTS_AND_MAAC_ALPHA2) + "\n-\n";
        Option a2 = OptionBuilder.withArgName("value").hasArgs(1).withValueSeparator().withDescription(description).create("a2");
        opts.addOption(a2);

        /* -------------------------------------------------------------------------------------- */
        /* Alpha cospeciation */
        description = "Cospeciation alpha for Dirichlet Process.\n" + "Default value = "
                      + format.format(ABC_SMC_Process.DEFAULT_ALPHA_COSPECIATION) + "\n-\n";
        Option ac = OptionBuilder.withArgName("value").hasArgs(1).withValueSeparator().withDescription(description).create("ac");
        opts.addOption(ac);

        /* -------------------------------------------------------------------------------------- */
        /* Alpha duplication */
        description = "Duplication alpha for Dirichlet Process.\n" + "Default value = "
                      + format.format(ABC_SMC_Process.DEFAULT_ALPHA_DUPLICATION) + "\n-\n";
        Option ad = OptionBuilder.withArgName("value").hasArgs(1).withValueSeparator().withDescription(description).create("ad");
        opts.addOption(ad);

        /* -------------------------------------------------------------------------------------- */
        /* Alpha loss */
        // description = ("Loss alpha for Dirichlet Process.\n")
        // + ("Default value = " + ABC_SMC_Process.DEFAULT_ALPHA_LOSS + "\n"
        // + "Ignored by models " + GenerateDistribution.COALESCENCE_WITHOUT_DISTANCE
        // + " and " + GenerateDistribution.COALESCENCE_WITH_DISTANCE + "\n-\n");
        description = "Loss alpha for Dirichlet Process.\n" + "Default value = "
                      + format.format(ABC_SMC_Process.DEFAULT_ALPHA_LOSS) + "\n-\n";
        Option al = OptionBuilder.withArgName("value").hasArgs(1).withValueSeparator().withDescription(description).create("al");
        opts.addOption(al);

        /* -------------------------------------------------------------------------------------- */
        /* Alpha switch */
        description = "Switch alpha for Dirichlet Process.\n" + "Default value = "
                      + format.format(ABC_SMC_Process.DEFAULT_ALPHA_SWITCH) + "\n-\n";
        Option as = OptionBuilder.withArgName("value").hasArgs(1).withValueSeparator().withDescription(description).create("as");
        opts.addOption(as);

        /* -------------------------------------------------------------------------------------- */
        /* Clustering */
        description = "Cluster the vectors accepted in the last round.\n"
                      + "WARNING: To perform this task, the program needs R\n"
                      + "and the R package \"dynamicTreeCut\".\n-\n";
        opts.addOption("cluster", false, description);

        /* -------------------------------------------------------------------------------------- */
        /* Continue */
        description = "Continue the process from the point where it had stopped.\n-\n";
        opts.addOption("continue", false, description);

        /* -------------------------------------------------------------------------------------- */
        /* Cyclicity test */
        description = "Cyclicity test which is used to avoid time inconsistent scenarios:\n"
                      + CyclicityTest.STOLZER + " - Stolzer et al., 2012\n" + CyclicityTest.DONATI
                      + " - Donati et al., 2014\n" + CyclicityTest.TOFIGH
                      + " - Tofigh et al., 2011\n" + "Default cyclicity test = "
                      + ABC_SMC_Process.DEFAULT_CYCLICITY_TEST + "\n-\n";
        Option test = OptionBuilder.withArgName("number").hasArgs(1).withValueSeparator().withDescription(description).create("ct");
        opts.addOption(test);

        /* -------------------------------------------------------------------------------------- */
        /* Maximum number of discarded vectors factor. */
        description = "Integer which determines the multiplier factor that is used\n"
                      + "to compute the maximum number of vectors which can be discarded\n"
                      + "during the generation of the quantile populations of the rounds\n"
                      + "greater than 1. The maximum number of discarded vector D is defined\n"
                      + "as: D = factor * Q, where Q = N * threshold_first_round.\n"
                      + "Notice that factor must be greater than 1 (factor > 1).\n"
                      + ("Default value = "
                         + ABC_SMC_Process.DEFAULT_MAXIMUM_NUMBER_OF_DISCARDED_VECTORS_FACTOR + "\n-\n");
        Option d = OptionBuilder.withArgName("factor").hasArgs(1).withValueSeparator().withDescription(description).create("discard");
        opts.addOption(d);

        /* -------------------------------------------------------------------------------------- */
        /* Help */
        opts.addOption("h", false, "Print the help message.\n-\n");

        /* -------------------------------------------------------------------------------------- */
        /* Input */
        description = "Input nexus file.\n-\n";
        Option input = OptionBuilder.withArgName("file").hasArgs(1).withValueSeparator().withDescription(description).create("input");
        opts.addOption(input);

        /* -------------------------------------------------------------------------------------- */
        /* Population number of trees */
        description = "Number of trees for each vector of the population.\n"
                      + ("Default value = " + ABC_SMC_Process.DEFAULT_NUMBER_OF_TREES + "\n-\n");
        Option n = OptionBuilder.withArgName("number").hasArgs(1).withValueSeparator().withDescription(description).create("M");
        opts.addOption(n);

        /* -------------------------------------------------------------------------------------- */
        /* Maximum number of discarded vectors factor. */
        description = "Integer number which determines the multiplier factor that is used\n"
                      + "to compute the maximum number of trees which are going to be simulated\n"
                      + "in order to obtain the required number of trees M (trees that are more\n"
                      + "than 2 times bigger than the 'real' parasite are discarded).\n"
                      + "The maximum number of trees X is defined as: X = factor * M.\n"
                      + "Notice that factor must be greater than 1 (factor > 1).\n"
                      + ("Default value = "
                         + ABC_SMC_Process.DEFAULT_MAXIMUM_NUMBER_OF_TREES_FACTOR + "\n-\n");
        Option max = OptionBuilder.withArgName("factor").hasArgs(1).withValueSeparator().withDescription(description).create("maxtrees");
        opts.addOption(max);

        /* -------------------------------------------------------------------------------------- */
        /* Metric */
        //BS
        description = "Distance to be used:\n" + GenerateDistribution.METRIC_MAAC_DISTANCE
                      + " - MAAC Metric\n" + GenerateDistribution.METRIC_LEAVES_AND_MAAC_DISTANCE
                      + " - LEAVES AND MAAC Metric\n"
                      + GenerateDistribution.METRIC_EVENTS_AND_MAAC_DISTANCE
                      + " - EVENTS AND MAAC Metric\n" +  " - EVENTS AND MAAC Metric for MultipleAssociations\n" + GenerateDistribution.METRIC_EVENTS_AND_MAAC_DISTANCE_MULTIPLEASSOCIATIONS
                      +"Default metric = "
                      + ABC_SMC_Process.DEFAULT_METRIC + ".\n-\n";
        Option dist = OptionBuilder.withArgName("number").hasArgs(1).withValueSeparator().withDescription(description).create("metric");
        opts.addOption(dist);

        /* -------------------------------------------------------------------------------------- */
        /* Model */
        // TODO future release
        // description = "Simulation model to be used:\n"
        // + GenerateDistribution.DEFAULT_MODEL
        // + " - From Root to the Leaves model\n"
        // + GenerateDistribution.COALESCENCE_WITHOUT_DISTANCE
        // + " - Coalescent model which ignores the information about host node distances.\n"
        // + GenerateDistribution.COALESCENCE_WITH_DISTANCE
        // + " - Coalescent model which considers the information about host node distances.\n"
        // + ("Default value = " + ABC_SMC_Process.DEFAULT_MODEL + "\n-\n");
        // Option model =
        // OptionBuilder.withArgName("value").hasArgs(1).withValueSeparator().withDescription(description).create("model");
        // opts.addOption(model);

        /* -------------------------------------------------------------------------------------- */
        /* Prior distribuition size */
        description = "Prior distribution size (Number of vectors that\n"
                      + " are going to be simulated in the first ABC-SMC round).\n"
                      + ("Default value = " + ABC_SMC_Process.DEFAULT_PRIOR_DISTRIBUTION_SIZE + "\n-\n");
        Option p = OptionBuilder.withArgName("size").hasArgs(1).withValueSeparator().withDescription(description).create("N");
        opts.addOption(p);

        /* -------------------------------------------------------------------------------------- */
        /* Perturbation */
        description = "Pertubation value. Real number in the interval (0,0.2).\n"
                      + ("Default value = " + format.format(ABC_SMC_Process.DEFAULT_PERTURBATION) + "\n-\n");
        Option perturbation = OptionBuilder.withArgName("value").hasArgs(1).withValueSeparator().withDescription(description).create("p");
        opts.addOption(perturbation);

        /* -------------------------------------------------------------------------------------- */
        /* Plot */
        description = "Create plots with the results of each round.\n"
                      + "WARNING: To perform this task, the program needs R\n"
                      + "and the R package \"ade4\".\n-\n";
        opts.addOption("plot", false, description);

        /* -------------------------------------------------------------------------------------- */
        /* Iterations */
        description = "Number of rounds of the ABC-SMC process.\n"
                      + "If the number of rounds is different from the\n"
                      + "default value, the parameter -t must be specified.\n" + "Default value = "
                      + ABC_SMC_Process.DEFAULT_NUMBER_OF_ITERATIONS + "\n-\n";
        Option iterations = OptionBuilder.withArgName("number").hasArgs(1).withValueSeparator().withDescription(description).create("R");
        opts.addOption(iterations);

        /* -------------------------------------------------------------------------------------- */
        /* Root mapping probability */
        description = "Root mapping probability. Real number in the interval [0.5,1.0].\n"
                      + ("Default value = "
                         + format.format(ABC_SMC_Process.DEFAULT_ROOT_MAPPING_PROBABILITY) + "\n-\n");
        Option root = OptionBuilder.withArgName("value").hasArgs(1).withValueSeparator().withDescription(description).create("root");
        opts.addOption(root);

        /* -------------------------------------------------------------------------------------- */
        /* Tolerance for each round. */
        String aux[] = new String[ABC_SMC_Process.DEFAULT_TOLERANCES.length];
        for (int i = 0; i < ABC_SMC_Process.DEFAULT_TOLERANCES.length; i++) {
            aux[i] = format.format(ABC_SMC_Process.DEFAULT_TOLERANCES[i]);
        }
        String tol = Arrays.toString(aux).replaceAll("[\\[\\] ]", "");

        description = "Tolerance level for each ABC-SMC round.\n"
                      + "It must be a list of R real numbers in the interval\n"
                      + "(0,1), where R is the number of rounds defined by\n"
                      + "the parameter -R (number of rounds).\n"
                      + "The numbers must be separated by comma characters.\n" + "Default value = "
                      + tol + "\n-\n";
        Option tolerance = OptionBuilder.withArgName("vector").hasArgs(1).withValueSeparator().withDescription(description).create("t");
        opts.addOption(tolerance);

        /* -------------------------------------------------------------------------------------- */
        /* Threads */
        description = "Number of threads to execute the process.\n" + "Default value = "
                      + ABC_SMC_Process.DEFAULT_NUMBER_OF_THREADS + "\n-\n";
        Option threads = OptionBuilder.withArgName("value").hasArgs(1).withValueSeparator().withDescription(description).create("threads");
        opts.addOption(threads);

        return opts;
    }


    private static boolean processParameters(ABC_SMC_Process process,
                                             String[] args,
                                             Options opts,
                                             String usageString) throws UnrecognizedOptionException,
                                                                ParseException {

        BasicParser bp = new BasicParser();
        CommandLine cl = bp.parse(opts, args);
        if (cl.hasOption("h")) {
            HelpFormatter f = new HelpFormatter();
            f.setWidth(120);
            f.printHelp(usageString, opts);
            return false;
        }

        // TODO official release
        int model = ABC_SMC_Process.DEFAULT_MODEL;
        process.setModel(model);
        // TODO for future
        // if (cl.hasOption("model")) {
        // model = Integer.parseInt(cl.getOptionValue("model"));
        // if (model != GenerateDistribution.DEFAULT_MODEL
        // && model != GenerateDistribution.COALESCENCE_WITHOUT_DISTANCE
        // && model != GenerateDistribution.COALESCENCE_WITH_DISTANCE) {
        // throw new UnrecognizedOptionException("Invalid value for option -model");
        // }
        // process.setModel(model);
        // }

        /* Get metric */
        int metric = ABC_SMC_Process.DEFAULT_METRIC;
        if (cl.hasOption("metric")) {
            metric = Integer.parseInt(cl.getOptionValue("metric"));
            if (metric != GenerateDistribution.METRIC_MAAC_DISTANCE
                && metric != GenerateDistribution.METRIC_LEAVES_AND_MAAC_DISTANCE
                && metric != GenerateDistribution.METRIC_EVENTS_AND_MAAC_DISTANCE
                && metric != GenerateDistribution.METRIC_EVENTS_AND_MAAC_DISTANCE_MULTIPLEASSOCIATIONS) {
                throw new UnrecognizedOptionException("Invalid value for option -metric");
            }
            process.setMetric(metric);
        }

        /* Get alpha 1 */
        double alpha1 = ABC_SMC_Process.DEFAULT_EVENTS_AND_MAAC_ALPHA1;
        if (metric == GenerateDistribution.METRIC_LEAVES_AND_MAAC_DISTANCE) {
            alpha1 = ABC_SMC_Process.DEFAULT_LEAVES_AND_MAAC_ALPHA1;
        }

        if (cl.hasOption("a1")) {
            alpha1 = Double.parseDouble(cl.getOptionValue("a1"));
            if (alpha1 < 0 || alpha1 > 1) {
                throw new UnrecognizedOptionException("Invalid value for option -a1");
            }
            process.setAlpha1(alpha1);
        }

        /* Get alpha 2 */
        double alpha2 = ABC_SMC_Process.DEFAULT_EVENTS_AND_MAAC_ALPHA2;
        if (metric == GenerateDistribution.METRIC_LEAVES_AND_MAAC_DISTANCE) {
            alpha2 = ABC_SMC_Process.DEFAULT_LEAVES_AND_MAAC_ALPHA2;
        }

        if (cl.hasOption("a2")) {
            alpha2 = Double.parseDouble(cl.getOptionValue("a2"));
            if (alpha2 < 0 || alpha2 > 1) {
                throw new UnrecognizedOptionException("Invalid value for option -a2");
            }
            process.setAlpha2(alpha2);
        }

        if (metric != GenerateDistribution.METRIC_MAAC_DISTANCE && alpha1 + alpha2 != 1.0) {
            throw new UnrecognizedOptionException("Alpha 1 (-a1) and alpha 2 (-a2) must sum up to 1.0.");
        }

        /* Get cospeciation alpha */
        if (cl.hasOption("ac")) {
            double alphaCospeciation = Double.parseDouble(cl.getOptionValue("ac"));
            if (alphaCospeciation <= 0) {
                throw new UnrecognizedOptionException("Invalid value for option -ac");
            }
            process.setAlphaCospeciation(alphaCospeciation);
        }

        /* Get duplication alpha */
        if (cl.hasOption("ad")) {
            double alphaDuplication = Double.parseDouble(cl.getOptionValue("ad"));
            if (alphaDuplication <= 0) {
                throw new UnrecognizedOptionException("Invalid value for option -ad");
            }
            process.setAlphaDuplication(alphaDuplication);
        }

        /* Get switch alpha */
        if (cl.hasOption("as")) {
            double alphaSwitch = Double.parseDouble(cl.getOptionValue("as"));
            if (alphaSwitch <= 0) {
                throw new UnrecognizedOptionException("Invalid value for option -as");
            }
            process.setAlphaSwitch(alphaSwitch);
        }

        /* Get loss alpha */
        if (model != GenerateDistribution.DEFAULT_MODEL && cl.hasOption("al")) {
            double alphaLoss = Double.parseDouble(cl.getOptionValue("al"));
            if (alphaLoss <= 0) {
                throw new UnrecognizedOptionException("Invalid value for option -al");
            }
        }

        /* Get plot */
        process.setClustering(cl.hasOption("cluster"));

        /* Get continue */
        process.setContinueProcess(cl.hasOption("continue"));

        /* Get the cyclicity test. */
        int cyclicityTest = ABC_SMC_Process.DEFAULT_CYCLICITY_TEST;
        if (cl.hasOption("ct")) {
            cyclicityTest = Integer.parseInt(cl.getOptionValue("ct"));
            if (cyclicityTest != CyclicityTest.STOLZER && cyclicityTest != CyclicityTest.TOFIGH
                && cyclicityTest != CyclicityTest.DONATI) {
                throw new UnrecognizedOptionException("Invalid value for option -ct");
            }
            process.setCyclicityTest(cyclicityTest);
        }

        /* Get number of rounds */
        int numberOfIterations = ABC_SMC_Process.DEFAULT_NUMBER_OF_ITERATIONS;
        if (cl.hasOption("R")) {
            numberOfIterations = Integer.parseInt(cl.getOptionValue("R"));
            if (numberOfIterations <= 0) {
                throw new UnrecognizedOptionException("Invalid value for option -R");
            }
            process.setNumberOfIterations(numberOfIterations);
        }

        /* Get input file name */
        if (cl.hasOption("input")) {
            process.setInputFile(cl.getOptionValue("input"));
        } else {
            throw new UnrecognizedOptionException("Missing input file.");
        }

        /* Get number of trees */
        if (cl.hasOption("M")) {
            int numberOfTrees = Integer.parseInt(cl.getOptionValue("M"));
            if (numberOfTrees <= 0) {
                throw new UnrecognizedOptionException("Invalid value for option -M");
            }
            process.setNumberOfTrees(numberOfTrees);
        }

        /* Get prior distribution size */
        if (cl.hasOption("N")) {
            int size = Integer.parseInt(cl.getOptionValue("N"));
            if (size <= 0) {
                throw new UnrecognizedOptionException("Invalid value for option -N");
            }
            process.setPriorDistributionSize(size);
        }

        /* Maximum number of discarded vectors factor. */
        if (cl.hasOption("discard")) {
            int factor = Integer.parseInt(cl.getOptionValue("discard"));
            if (factor <= 1) {
                throw new UnrecognizedOptionException("Invalid value for option -discard");
            }
            process.setMaximumNumberOfDiscardedVectorsFactor(factor);
        }

        /* Maximum number of trees factor. */
        if (cl.hasOption("maxtrees")) {
            int factor = Integer.parseInt(cl.getOptionValue("maxtrees"));
            if (factor <= 1) {
                throw new UnrecognizedOptionException("Invalid value for option -maxtrees");
            }
            process.setMaximumNumberOfTreesFactor(factor);
        }

        /* Get plot */
        process.setPlot(cl.hasOption("plot"));

        /* Get tolerance levels */
        if (numberOfIterations != ABC_SMC_Process.DEFAULT_NUMBER_OF_ITERATIONS
            && !cl.hasOption("t")) {
            throw new UnrecognizedOptionException("Number of rounds different from default: missing parameter -t");
        }

        if (cl.hasOption("t")) {
            String[] elements = cl.getOptionValue("t").split(",");
            if (elements.length != numberOfIterations) {
                throw new UnrecognizedOptionException("Parameter -t: vector must have length equal to the number of rounds.");
            }
            double[] tolerances = new double[numberOfIterations];
            for (int i = 0; i < numberOfIterations; i++) {
                try {
                    tolerances[i] = Double.parseDouble(elements[i]);
                    if (tolerances[i] <= 0.0 || tolerances[i] >= 1.0) {
                        throw new UnrecognizedOptionException("Parameter -t: invalid tolerance value = "
                                                              + elements[i]);
                    }
                } catch (NumberFormatException e) {
                    throw new UnrecognizedOptionException("Parameter -t: invalid tolerance value = "
                                                          + elements[i]);
                }
            }
            process.setTolerances(tolerances);
        }

        /* Get threads */
        if (cl.hasOption("threads")) {
            int threads = Integer.parseInt(cl.getOptionValue("threads"));
            if (threads < 0) {
                throw new UnrecognizedOptionException("Invalid value for option -threads");
            }
            process.setNumberOfThreads(threads);
        }

        /* Get perturbation */
        if (cl.hasOption("p")) {
            double perturbation = Double.parseDouble(cl.getOptionValue("p"));
            if (perturbation <= 0.0 || perturbation >= 0.2) {
                throw new UnrecognizedOptionException("Invalid value for option -p");
            }
            process.setPerturbation(perturbation);
        }

        /* Get root mapping probability */
        if (cl.hasOption("root")) {
            double probability = Double.parseDouble(cl.getOptionValue("root"));
            if (probability < 0.5 || probability > 1.0) {
                throw new UnrecognizedOptionException("Invalid value for option -root");
            }
            process.setRootMappingProbability(probability);
        }

        return true;
    }
}
