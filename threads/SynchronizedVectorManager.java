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

import java.io.PrintWriter;
import java.text.Format;
import java.util.ArrayList;
import java.util.Arrays;

public class SynchronizedVectorManager {

    private int numberOfVectors;

    private int numberOfRejectedVectors;

    private ArrayList<double[]> vectors;

    private PrintWriter writer;

    private Format format;


    /**
     * Thread to generate vectors for the prior distribution.
     * 
     * @param numberOfVectors
     *        Maximum number of vectors that the array can receive.
     * @param writer
     *        PrintWriter object for the file that will receive the generated vectors.
     * @param format
     *        Object to format the real numbers.
     */
    public SynchronizedVectorManager(int numberOfVectors, PrintWriter writer, Format format) {
        this.numberOfVectors = numberOfVectors;
        this.numberOfRejectedVectors = 0;
        this.vectors = new ArrayList<double[]>();
        this.writer = writer;
        this.format = format;
    }


    /**
     * Writes the given vector to the output file.
     * 
     * @param vector
     *        Vector
     */
    public synchronized void registerVector(double[] vector) {
        if (vectors.size() < numberOfVectors) {
            vectors.add(vector);
            writer.println(formatVectorToPrint(vector));
            writer.flush();
        }
    }


    /**
     * Returns the number of vectors which were generated.
     * 
     * @return The number of vectors which were generated.
     */
    protected synchronized int getNumberOfVectors() {
        return vectors.size();
    }


    /**
     * Increases the number of rejected vectors
     */
    public synchronized void increaseNumberOfRejectedVectors() {
        this.numberOfRejectedVectors++;
    }


    /**
     * Returns the number of rejected vectors.
     * 
     * @return The number of rejected vectors.
     */
    public synchronized int getNumberOfRejectedVectors() {
        return numberOfRejectedVectors;
    }


    /**
     * Returns the list of vectors.
     * 
     * @return The list of vectors.
     */
    public ArrayList<double[]> getVectors() {
        return vectors;
    }


    /**
     * Formats vector to print.
     * 
     * @param vector
     *        Vector to be formatted.
     * @return The String that represent the vector.
     */
    private String formatVectorToPrint(double[] vector) {
        String[] str = new String[vector.length];
        for (int i = 0; i < vector.length; i++) {
            str[i] = format.format(vector[i]);
        }
        return Arrays.toString(str).replaceAll("[\\[\\],]", "").replaceAll(" ", "\t");
    }

}
