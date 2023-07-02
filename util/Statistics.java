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

package  util;

import java.util.Arrays;

public class Statistics {

    public static double mean(double[] l) {
        double sum = 0.0;
        int n = l.length;
        for (int e = 0; e < n; e++) {
            sum = sum + l[e];
        }
        return sum / n;
    }


    public static double mean(int[] l) {
        double sum = 0.0;
        int n = l.length;
        for (int e = 0; e < n; e++) {
            sum = sum + l[e];
        }
        return sum / n;
    }

    public static double max(double[] l) {
    	double maxVal;
    	if(l.length==0){
    		maxVal = 0.0;
    	}
    	else{
    		maxVal=l[0];
    		for (int i = 1; i < l.length; i++){
    			if (maxVal < l[i]){
    				maxVal = l[i];
    			}
    		}
    	}	
        return maxVal;
    }

    
    public static int max(int[] l) {
    	int maxVal;
    	if(l.length==0){
    		maxVal = 0;
    	}
    	else{
    		maxVal=l[0];
    		for (int i = 1; i < l.length; i++){
    			if (maxVal < l[i]){
    				maxVal = l[i];
    			}
    		}
    	}	
        return maxVal;
    }
    
    
    public static double median(double[] l) {
        int n = l.length;
        Arrays.sort(l);
        double median = l[n / 2];
        if (n % 2 == 0) {
            median = (median + l[n / 2 - 1]) / 2;
        }
        return median;
    }


    public static double median(int[] l) {
        int n = l.length;
        Arrays.sort(l);
        double median = l[n / 2];
        if (n % 2 == 0) {
            median = (median + l[n / 2 - 1]) / 2;
        }
        return median;
    }


    public static double pearson(double[] l) {
        int n = l.length;
        Arrays.sort(l);
        double median = l[n / 2];
        if (n % 2 == 0) {
            median = (median + l[n / 2 - 1]) / 2;
        }
        double sum = 0.0;
        double mean = mean(l);
        for (int e = 0; e < n; e++) {
            sum = sum + (mean - l[e]) * (mean - l[e]);
        }
        double std = Math.sqrt(sum / (n - 1));
        double pearson = 3 * (mean - median) / std;
        return pearson;
    }


    public static double pearson(int[] l) {
        int n = l.length;
        Arrays.sort(l);
        double median = l[n / 2];
        if (n % 2 == 0) {
            median = (median + l[n / 2 - 1]) / 2;
        }
        double sum = 0.0;
        double mean = mean(l);
        for (int e = 0; e < n; e++) {
            sum = sum + (mean - l[e]) * (mean - l[e]);
        }
        double std = Math.sqrt(sum / (n - 1));
        double pearson = 3 * (mean - median) / std;
        return pearson;
    }


    public static double sd(double[] l) {
        if (l.length > 1) {
            double mean = mean(l);
            double sum = 0.0;
            for (double value: l) {
                sum += Math.pow(value - mean, 2);
            }
            return Math.sqrt(sum / (l.length - 1));
        }
        return 0;
    }


    public static double sd(int[] l) {
        if (l.length > 1) {
            double mean = mean(l);
            double sum = 0.0;
            for (double value: l) {
                sum += Math.pow(value - mean, 2);
            }
            return Math.sqrt(sum / (l.length - 1));
        }
        return 0.0;
    }

}
