/*
 *  BEANDiscoPrior: functions for computing scores with priors
 *  
 *  Copyright 2014 Diane Oyen <doyen at cs.unm.edu>
 * 
 *  Modified from BEANDisco:
 *  Copyright 2011 Teppo Niinimäki <teppo.niinimaki(at)helsinki.fi>
 *  
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


#include <cmath>

#include "common.hpp"
#include "data.hpp"
#include "stacksubset.hpp"
#include "arc.hpp"

#ifndef SCORES_HPP
#define SCORES_HPP

double* logGammas = NULL;
void computeLogGammas(int n) {
	logGammas = new double[n + 1];
	for (int i = 0; i <= n; ++i)
		logGammas[i] = lgamma(i);
}

void freeLogGammas() {
	assert(logGammas != NULL);
	delete[] logGammas;
}

double BDELogScore(int nValues, int nParentValues, int* counts, int pseudocount) {
	double score = 0;
	for (int pv = 0; pv < nParentValues; ++pv) {
		int cumCount = 0;
		for (int v = 0; v < nValues; ++v) {
			int c = counts[pv * nValues + v];
			//score += lgamma(c + pseudocount) - lgamma(pseudocount);
			score += logGammas[c + pseudocount] - logGammas[pseudocount];
			cumCount += c;
		}
		//score += lgamma(nValues * pseudocount) - lgamma(cumCount + nValues * pseudocount);
		score += logGammas[nValues * pseudocount] - logGammas[cumCount + nValues * pseudocount];
	}
	return score;
}/**/

double computeScore(const Data& data, const StackSubset& parents, int node) {
	int nParentValues = 1;
	for (int i = 0; i < parents.size(); ++i)
		nParentValues *= data.arities[parents[i]];
	int nNodeValues = data.arities[node];
	int* counts = new int[nParentValues * nNodeValues];
	for (int i = 0; i < nParentValues * nNodeValues; ++i)
		counts[i] = 0;
	//memset(counts, 0, nParentValues * nNodeValues * sizeof(int));
	
	for (int j = 0; j < data.nSamples; ++j) {
		int index = 0;
		for (int i = 0; i < parents.size(); ++i)
			index = index * data.arities[parents[i]] + data(parents[i], j);
		index = index * data.arities[node] + data(node, j);
		++counts[index];
	}
	
	double score = BDELogScore(nNodeValues, nParentValues, counts, 1);
	delete[] counts;

	return score;
}/**/


/**
 * Energy function (just a piece of the prior) modularized per parentset or ideal
 **/ 
double computeEnergy(const StackSubset& parents, int node, ArcMap<Real>* priors) {
  double parentsetEnergy = 0;
  Arc arc;
  arc.head = node;
  for (int i = 0; i < priors->getNnodes(); ++i) {
    if (node == i)
      continue;
    arc.tail = i;
    if (parents.contains(i)) {
      parentsetEnergy += (1 - (*priors)[arc].getVal());
    } else {
      parentsetEnergy += (*priors)[arc].getVal();
    }
  }
  return parentsetEnergy;
}


#endif

