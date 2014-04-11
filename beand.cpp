/*
 *  BEANDiscoPrior: main program
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

#include <algorithm>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include "common.hpp"
#include "logger.hpp"
#include "lognum.hpp"
#include "timer.hpp"
#include "data.hpp"
#include "stacksubset.hpp"
#include "scores.hpp"
#include "parentsetmap.hpp"
#include "bucketorder.hpp"
#include "parbucketorder.hpp"
#include "arc.hpp"

//#define NDEBUG
#include <cassert>


#define BEAND_VERSION_STRING  "1.0.0"


// create a logger
Logger logger;


/**
 * Computes K2 scores for all node-parentset pairs
 */
void computeScores(const Data& data, ParentsetMap<Real>& scores, ParentsetMap<Real>* parentPriors = NULL) {
	computeLogGammas(2 * data.nSamples);
	StackSubset parents(scores.maxParents);
	for (int node = 0; node < scores.nNodes; ++node) {
		parents.clear();
		do {
			if (parents.contains(node))
				continue;
			double logScore = computeScore(data, parents, node);
			if (parentPriors != NULL) {
			  logScore += logScore + (*parentPriors)(node, parents).getLog();
			}
			//printf("%g ", logScore);
			Lognum<double> tmp;
			tmp.setLog(logScore);
			scores(node, parents) = to<Real>(tmp);
		} while (parents.next(0, scores.nNodes, scores.maxParents));
	}
	freeLogGammas();
}

/** 
 * Computes priors per parentset from arc priors
 */
void computePriors(ArcMap<Real>*arcPriors, double beta, ParentsetMap<Real>& parentPriors) {
  StackSubset parents(parentPriors.maxParents);
  for (int node = 0; node < parentPriors.nNodes; ++node) {
    parents.clear();
    double norm = 0;
    do {
      if (parents.contains(node))
	continue;
      double logPrior = -1 * beta * computeEnergy(parents, node, arcPriors);
      Lognum<double> tmp;
      tmp.setLog(logPrior);
      parentPriors(node, parents) = to<Real>(tmp);
      norm += pow(2, parentPriors.nNodes - 1 - parents.size()) * parentPriors(node, parents).getVal();
    } while (parents.next(0, parentPriors.nNodes, parentPriors.maxParents));
    // Need to divide out the normalization constant
    if (norm != 0) { // avoid divide by zero
      do {
	if (parents.contains(node))
	  continue;
	parentPriors(node, parents).setLog(parentPriors(node, parents).getLog() - log(norm));
      } while (parents.next(0, parentPriors.nNodes, parentPriors.maxParents));
    }
  }
}


/**
 * Adjusts scores to incorporate prior (assuming it was not factored in previously).
 */
void adjustScores(ParentsetMap<Real>& scores, ParentsetMap<Real>* parentPriors) {
  StackSubset parents(scores.maxParents);
  for (int node = 0; node < scores.nNodes; ++node) {
    parents.clear();
    do {
      if (parents.contains(node))
	continue;
      // score (without prior) already computed, just factor in the prior
      scores(node, parents) += (*parentPriors)(node, parents).getLog();
    } while (parents.next(0, scores.nNodes, scores.maxParents));
  }
}


template <class T>
T binom(int n, int k) {
	return round(exp(lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1)));
}

template <>
Lognum<double> binom<Lognum<double> >(int n, int k) {
	Lognum<double> res;
	res.setLog(lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1));
	return res;
}


Real* preCalcInvBinoms(int n) {
	Real* invBinoms = new Real[n];
	for (int k = 0; k < n; ++k)
		invBinoms[k] = Real(1.0) / binom<Real>(n, k);
	return invBinoms;
}

/**
 * Divides scores by the number of parent sets with same size.
 */
void weightScores(ParentsetMap<Real>& scores) {
	Real* invBinoms = preCalcInvBinoms(scores.nNodes - 1);
	StackSubset pa(scores.maxParents);
	for (int v = 0; v < scores.nNodes; ++v) {
		//pa.length = 0;
		do {
			scores(v, pa) *= invBinoms[pa.size()];
		} while (pa.next(0, scores.nNodes, scores.maxParents));
	}
	delete[] invBinoms;
}



/**
 * Computes tail sums of scores in given partial order.
 */
template <class POF>
void calcTailSums(
		const POF& pof,
		const typename POF::Order& po,
		ParentsetMap<Real>& scores,
		int node,
		Arc feature,
		typename POF::template IdealMap<Real>& tailSums
		) {
	tailSums.setAll(0.0);
	StackSubset x(scores.maxParents);
	StackSubset xt(scores.maxParents); // translated x
	typename POF::Ideal xId(pof);
	do {
		translateSubset(po.getOrder(), x, xt);
		//size_t xti = scores.getParentsetIndex(xt);
		xId.setSuperOf(x);
		if (feature.holds(node, xt))
			tailSums[xId] += scores(node, xt);
	} while (x.next(0, scores.nNodes, scores.maxParents));
}/**/


/**
 *
 */
template <class POF>
Real getRho(
	    const POF& pof,
	    const typename POF::Order& po,
	    int v,
	    typename POF::Ideal ideal,
	    ArcMap<Real>* priors,
	    double beta
	    ) {
        StackSubset x(pof.n);
	StackSubset xt(pof.n);
	ideal.getElements(x);
	translateSubset(po.getOrder(), x, xt);
        //std::cout << xt << std::endl;

	double logRho = -1 * beta * computeEnergy(xt, po[v], priors);
	Lognum<double> tmp;
	tmp.setLog(logRho);
	Real rho = to<Real>(tmp);
	return rho;
}



/**
 * Computes alphas from scores.
 */
template <class POF>
void calcAlphas(
		const POF& pof,
		const typename POF::Order& po,
		ParentsetMap<Real>& scores,
		Arc arc,
		ArcMap<Real>* priors,
		double beta,
		std::vector<typename POF::template IdealMap<Real> >& alphas
		) {
	for (int v = 0; v < pof.n; ++v) {
		calcTailSums(pof, po, scores, po[v], arc, alphas[v]);
		alphas[v].fastSparseZetaTransform();

		// order priors
		if (priors) {
		  typename POF::Ideal ideal(pof);
		  ideal.setFirstExpandableWith(v);
		  do {
		    //std::cout << "calcAlphas: v=" << v << std::endl;
		    alphas[v][ideal] *= getRho(pof, po, v, ideal, priors, beta);
		  } while (ideal.nextExpandableWith(v));
		}
	}
}



/**
 * Computes gammas from forward and backward sums.
 */
template <class POF>
void calcGammas(
		const POF& pof,
		const typename POF::Order& po,
		ArcMap<Real>* priors,
		double beta,
		typename POF::template IdealMap<Real>& fp,
		typename POF::template IdealMap<Real>& bp,
		std::vector<typename POF::template IdealMap<Real> >& gamma
		) {
	// for each variable
	for (int v = 0; v < pof.n; ++v) {
		// iterate over all ideals y
		typename POF::Ideal y(pof);
		do {
			if (y.isShrinkableWith(v)) {
				Real tmp = bp[y];
				y.shrinkWith(v);
				if (priors)
				  gamma[v][y] = fp[y] * tmp * getRho(pof, po, v, y, priors, beta);
				else 
				  gamma[v][y] = fp[y] * tmp;
				y.expandWith(v);
			}
		} while(y.next());
		gamma[v].fastSparseUpZetaTransform();
	}
}/**/


/**
 * Computes the final (unnormalized) probabilities for each arc from gammas and local scores.
 */
template <class POF>
void addParentsetSums(
		const POF& pof,
		ParentsetMap<Real>& scores,
		std::vector<typename POF::template IdealMap<Real> >& gammas,
		const typename POF::Order& po,
		ArcMap<Real>& sums
		) {
	StackSubset pa(scores.maxParents);
	StackSubset pat(scores.maxParents);
	typename POF::Ideal paId(pof);
	do {
		translateSubset(po.getOrder(), pa, pat);
		//size_t pai = scores.getParentsetIndex(pa);
		size_t pati = scores.getParentsetIndex(pat);
		paId.setSuperOf(pa);
		
		Arc arc;
		for (int i = 0; i < pat.size(); ++i) {
			arc.tail = pat[i];
			for (int headt = 0; headt < pof.n; ++headt) {
				arc.head = po.getOrder()[headt];
				if (arc.head == arc.tail)
					continue;
				sums[arc] += scores(arc.head, pati) * gammas[headt][paId];
			}
		}
	} while (pa.next(0, scores.nNodes, scores.maxParents));
	
}/**/


/**
 * Computes the (unnormalized) probability of given arc in given partial order.
 */
template <class POF>
Real calcUnnormProb(
		const POF& pof,
		ParentsetMap<Real>& scores,
		const typename POF::Order& po,
		ArcMap<Real>* priors,
		double beta,
		Arc arc
	) {
	
	// compute alphas
	std::vector<typename POF::template IdealMap<Real> >
			alphas(pof.n, typename POF::template IdealMap<Real>(pof));
	calcAlphas(pof, po, scores, arc, priors, beta, alphas);
	
	// compute the probability
	typename POF::template IdealMap<Real> fp(pof);
	fp.sparseForwardSum(alphas);

	Real p = fp.getFull();
	
	return p;
}


/**
 * Computes the (unnormalized) probabilities of all arc simultaneously in given partial order.
 */
template <class POF>
void calcUnnormArcProbs(
		const POF& pof,
		ParentsetMap<Real>& scores,
		const typename POF::Order& po,
		ArcMap<Real>* priors,
		double beta,
		ArcMap<Real>& probs
	) {
	// compute alphas for null feature
	std::vector<typename POF::template IdealMap<Real> >
			nullAlphas(pof.n, typename POF::template IdealMap<Real>(pof));
	calcAlphas(pof, po, scores, NullArc, priors, beta, nullAlphas);
	
	// compute forward and backward functions
	typename POF::template IdealMap<Real> fp(pof);
	fp.sparseForwardSum(nullAlphas);
	typename POF::template IdealMap<Real> bp(pof);
	bp.sparseBackwardSum(nullAlphas);
	
	// compute gammas
	std::vector<typename POF::template IdealMap<Real> >
			gammas(pof.n, typename POF::template IdealMap<Real>(pof));
	calcGammas(pof, po, priors, beta, fp, bp, gammas);
	
	// compute all arc probs at once
	probs.setAll(0.0);
	addParentsetSums(pof, scores, gammas, po, probs);
}



template <class POF>
class ExactArcProbComputer {
private:
	ParentsetMap<Real>& scores_;
	const POF& pof_;
        ArcMap<Real>* orderPriors_;
        double beta_;
	
public:
	ExactArcProbComputer(ParentsetMap<Real>& scores, const POF& pof) :
		scores_(scores), pof_(pof)
	{
	  orderPriors_ = NULL;
	}
	
	~ExactArcProbComputer() {
	}
	
        void setOrderPriors(ArcMap<Real>* priors) {
	        orderPriors_ = priors;
	}

        void setPriorBeta(double beta) {
	        beta_ = beta;
	}

	double calcProb(Arc arc) {
		//ParentsetMap<Real> transScores;
		
		Real cumMarginalLikelihood = 0;
		Real cumArcLikelihood = 0;
		typename POF::OrderEnumerator poe(pof_);
		do {
		        Real lhPO = calcUnnormProb(pof_, scores_, poe.getOrder(), orderPriors_, beta_, NullArc);
			Real lhFPO = calcUnnormProb(pof_, scores_, poe.getOrder(), orderPriors_, beta_, arc);
			cumMarginalLikelihood += lhPO;
			cumArcLikelihood += lhFPO;
		} while(poe.next());
		
		return to<double>(cumArcLikelihood / cumMarginalLikelihood);
	}
	
	void printAllProbs(std::ostream& resStream) {
		Arc arc; arc.setFirst();
		do {
			double p = calcProb(arc);
			resStream << arc << "   " << p << std::endl;
		} while (arc.next(pof_.n));
	}

	void printArcProbs(std::ostream& resStream) {
		Real cumMarginalProb = 0;
		ArcMap<Real> cumArcProbs(pof_.n);
		ArcMap<Real> probs(pof_.n);
		cumArcProbs.setAll(0.0);
		typename POF::OrderEnumerator poe(pof_);
		do {
		        Real marginalProb = calcUnnormProb(pof_, scores_, poe.getOrder(), orderPriors_, beta_, NullArc);
			calcUnnormArcProbs(pof_, scores_, poe.getOrder(), orderPriors_, beta_, probs);
		        cumMarginalProb += marginalProb;
			Arc arc; arc.setFirst();
			do {
				cumArcProbs[arc] += probs[arc];
			} while (arc.next(pof_.n));
		} while(poe.next());
		
		Arc arc; arc.setFirst();
		do {
			double p = to<double>(cumArcProbs[arc] / cumMarginalProb);
			resStream << arc << "   " << p << std::endl;
		} while (arc.next(pof_.n));
	}

};



template <class POF>
class MCMCArcProbComputer {
private:
	ParentsetMap<Real>& scores_;
	const POF& pof_;
        ArcMap<Real>* orderPriors_;
        double beta_;
	
	Real marginUnnormProb_;
	typename POF::Order po_;
	
	int nAccepts_;
	int nSteps_;
	
	std::ostream* marginStream_;
	
public:
  MCMCArcProbComputer(ParentsetMap<Real>& scores, const POF& pof, ArcMap<Real>* orderPriors, double beta) :
		scores_(scores), pof_(pof), po_(pof), marginStream_(NULL)
	{
		logger.println(1, "  Initialize starting state (random permutation)...");
		po_.rand();
		
		logger.println(1, "  Compute initial probability...");
		orderPriors_ = orderPriors;
		beta_ = beta;
		marginUnnormProb_ = calcUnnormProb(pof_, scores_, po_, orderPriors_, beta_, NullArc);
		
		resetStats();
	}
	
	~MCMCArcProbComputer() {
	}

	void resetStats() {
		nAccepts_ = 0;
		nSteps_ = 0;
	}
	
	double getAcceptRatio() {
		return nAccepts_ / (double) nSteps_;
	}
	
	void setMarginStream(std::ostream& targetStream) {
		marginStream_ = &targetStream;
	}
	
        void setOrderPriors(ArcMap<Real>* priors) {
	        orderPriors_ = priors;
	}

        void setPriorBeta(double beta) {
	        beta_ = beta;
	}

	void mcmcStep(int nSwaps = 1) {
		typename POF::Order poNew(pof_);
		poNew = po_;
		for (int i = 0; i < nSwaps; ++i)
			poNew.randSwap();
		Real pnew = calcUnnormProb(pof_, scores_, poNew, orderPriors_, beta_, NullArc);
		++nSteps_;
		if (randu() < to<double>(pnew / marginUnnormProb_)) {
			po_ = poNew;
			marginUnnormProb_ = pnew;
			++nAccepts_;
		}
		if (marginStream_)
			(*marginStream_) << to<double>(log(marginUnnormProb_)) << std::endl;
	}
	
	
	void temperedIdleRun(int nSteps, int nSwaps = 1) {
		for (int i = 1; i <= nSteps; ++i) {
			int ns = 1 + nSwaps * (nSteps - i) / nSteps;
			mcmcStep(ns);
		}
	}
	
	void idleRun(int nSteps) {
		for (int i = 0; i < nSteps; ++i) {
			mcmcStep();
		}
	}
	
	double calcProb(int nSamples, int nStepsPerSample, Arc arc, double* probs = NULL) {
		
		double psum = 0.0;
		for (int i = 0; i < nSamples; ++i) {
			for (int j = 0; j < nStepsPerSample; ++j) {
				mcmcStep();
			}
			Real arcUnnormProb = calcUnnormProb(pof_, scores_, po_, orderPriors_, beta_, arc);
			double pi = to<double>(arcUnnormProb / marginUnnormProb_);
			if (probs)
				probs[i] = pi;
			psum += pi;
		}
		double p = psum / nSamples;
		
		return p;
	}
	
	void printAllProbs(std::ostream& resStream, int nSamples, int nStepsPerSample,
			bool printSamples = false) {
		double* samples = NULL;
		if (printSamples)
			samples = new double[nSamples];
		Arc arc; arc.setFirst();
		do {
			double p = calcProb(nSamples, nStepsPerSample, arc, samples);
			resStream << arc << "   " << p;
			if (printSamples)
				for (int i = 0; i < nSamples; ++i)
					resStream << "  " << samples[i];
			resStream << std::endl;
		} while (arc.next(pof_.n));
		if (printSamples)
			delete[] samples;
	}
	
	void printArcProbs(std::ostream& resStream, int nSamples, int nStepsPerSample,
			bool printSamples = false) {
		
		ArcMap<double> cumArcProbs(pof_.n);
		ArcMap<Real> arcUnnormProbs(pof_.n);
		
		cumArcProbs.setAll(0.0);
		for (int i = 0; i < nSamples; ++i) {
			for (int j = 0; j < nStepsPerSample; ++j) {
				mcmcStep();
			}
			
			calcUnnormArcProbs(pof_, scores_, po_, orderPriors_, beta_, arcUnnormProbs);
			
			Arc arc; arc.setFirst();
			do {
				double p = to<double>(arcUnnormProbs[arc] / marginUnnormProb_);
				if (printSamples)
					resStream << "  " << p;
				cumArcProbs[arc] += p;
			} while (arc.next(pof_.n));
			
			if (printSamples) {
				resStream << std::endl;
			}
		}
		
		if (!printSamples) {
			Arc arc; arc.setFirst();
			do {
				double p = cumArcProbs[arc] / nSamples;
				resStream << arc << "   " << p << std::endl;
			} while (arc.next(pof_.n));
		}
	}
};



/**
 * Writes local scores to file.
 */
void writeScores(std::ostream& file, const ParentsetMap<Real>& scores) {
	file << scores.nNodes << std::endl;
	file << scores.maxParents << std::endl;
	file.precision(16);
	StackSubset parents(scores.maxParents);
	for (int node = 0; node < scores.nNodes; ++node) {
		parents.clear();
		do {
			if (parents.contains(node))
				continue;
			file << log(scores(node, parents)) << " ";
		} while (parents.next(0, scores.nNodes, scores.maxParents));
		file << std::endl;
	}
}


/**
 * Reads local scores from file.
 */
ParentsetMap<Real>* readScores(std::istream& file) {
	int nNodes, maxParents;
	file >> nNodes;
	file >> maxParents;
	ParentsetMap<Real>* scores = new ParentsetMap<Real>(nNodes, maxParents);
	StackSubset parents(maxParents);
	for (int node = 0; node < nNodes; ++node) {
		parents.clear();
		do {
			if (parents.contains(node))
				continue;
			double tmp;
			file >> tmp;
			if (file.fail())
				throw Exception("File corrupted; could not read all scores.");
			Lognum<double> tmp2;
			tmp2.setLog(tmp);
			(*scores)(node, parents) = to<Real>(tmp2);
		} while (parents.next(0, nNodes, maxParents));
	}
	file >> std::ws;
	if (!file.eof())
		throw Exception("File corrupted; contains more data than expected.");
	return scores;
}


/**
 * Read prior probabilities from file.
 */
void readPriorFile(std::istream& file, ArcMap<Real> *priors) {
  int nNodes = priors->getNnodes();
  Arc arc;

  double tmp;
  // Rows are parents (tails), columns are children (heads)
  for (int nodeTail = 0; nodeTail < nNodes; ++nodeTail) {
    for (int nodeHead = 0; nodeHead < nNodes; ++nodeHead) {
      file >> tmp;
      if (file.fail())
	throw Exception("File corrupted; could not read all priors.");
      if (nodeHead == nodeTail) {
	continue;
      } else {
	Lognum<double> tmp2;
	tmp2 = tmp;
	arc.head = nodeHead;
	arc.tail = nodeTail;
	(*priors)[arc] = to<Real>(tmp2);
      }
    }
  }

  file >> std::ws;
  if (!file.eof())
    throw Exception("File <prior> corrupted; contains more data than expected.");

}


using namespace std;


/**
 * Compute and print exact probabilities.
 */
template <class POF>
void printExact(std::ostream& resStream, std::ostream& logStream,
		const POF& pof, ParentsetMap<Real>& scores,
		ArcMap<Real>* priors, double beta) {
	logger.println(1, "Initialize...");
	ExactArcProbComputer<POF> eapc(scores, pof);
	if (priors) {
	  eapc.setOrderPriors(priors);
	  eapc.setPriorBeta(beta);
	}
	logger.printfln(1, "Actual computation...");
	eapc.printArcProbs(resStream);
}


/**
 * Compute and print MCMC approximated probabilities.
 */
template <class POF>
void printMCMC(std::ostream& resStream, std::ostream& marginStream, std::ostream& logStream,
		const POF& pof, ParentsetMap<Real>& scores, int nBurnInSteps, int nSamples,
	       int nStepsPerSample, int nBurnOutSteps, bool printSamples,
	       ArcMap<Real>* priors, double beta) {
	logger.println(1, "Initialize...");
	MCMCArcProbComputer<POF> mapc(scores, pof, priors, beta);
	mapc.setMarginStream(marginStream);
	
	// burn-in
	if (nBurnInSteps > 0) {
		Timer burninTimer; burninTimer.start();
		logger.printfln(1, "Burn-in (%d steps)...", nBurnInSteps);
		mapc.temperedIdleRun(nBurnInSteps, 1);
		double burninTime = burninTimer.elapsed();
		logStream << "burnin_time = " << burninTime << endl;
		logger.printfln(1, "  Elapsed %.2f s.", burninTime);
	}
	
	// actual sampling
	if (nSamples > 0) {
		logger.printfln(1, "Actual computation (%d samples x %d steps)...",
				nSamples, nStepsPerSample);
		mapc.resetStats();
		Timer samplingTimer; samplingTimer.start();
		mapc.printArcProbs(resStream, nSamples, nStepsPerSample, printSamples);
		double samplingTime = samplingTimer.elapsed();
		logStream << "sampling_time = " << samplingTime << endl;
		logger.printfln(1, "  Elapsed %.2f s.", samplingTime);
		double acceptRatio = mapc.getAcceptRatio();
		logger.printfln(1, "  Acceptance ratio was %.3g.", acceptRatio);
		logStream << "sampling_acceptance_ratio = " << acceptRatio << endl;
	}
	
	// burn-out
	if (nBurnOutSteps > 0) {
		Timer burnoutTimer; burnoutTimer.start();
		logger.printfln(1, "Burn-out (%d steps)...", nBurnOutSteps);
		mapc.idleRun(nBurnOutSteps);
		double burnoutTime = burnoutTimer.elapsed();
		logStream << "burnout_time = " << burnoutTime << endl;
		logger.printfln(1, "  Elapsed %.2f s.", burnoutTime);
	}
}



/*
 * Main program.
 */

#include <string>
#include <iomanip>
#include <boost/program_options.hpp>

namespace opts = boost::program_options;

int main(int argc, char** argv) {

	string inFilename;
	string outFilename;
	string scoreFilename;
	string marginOutFilename;
	string logFilename;
	string priorFilename;
	
	int nVariables;
	int maxIndegree;
	int nDataSamples;
	
	string poType;
	int maxBucketSize;
	int nChains;
	
	int nSamples;
	int nStepsPerSample;
	int nBurnInSteps;
	int nBurnOutSteps;

	double beta;

	unsigned int rngSeed;
	
	int verbosity;
	
	opts::options_description desc("Options");
	desc.add_options()
		("help,h",          "produce help message")
		("exact,e",         "use exact computation instead of MCMC")
		("test-conv",       "just test MCMC convergence")
		("verbose,v",       opts::value<int>(&verbosity)->default_value(0)->implicit_value(1),
		                    "set verbosity level")
		("quiet,q",         "use quiet mode, does not print anything unnecessary")
		("num-rows,r",      opts::value<int>(&nDataSamples)->default_value(0),
		                    "set number of data rows (samples)")
		("num-variables,n", opts::value<int>(&nVariables)->default_value(0),
		                    "set number of variables")
		("max-indegree,m",  opts::value<int>(&maxIndegree)->default_value(0),
		                    "set maximum indegree")
		("order-type",      opts::value<string>(&poType)->default_value("bo"),
		                    "set partial order type, possible values: bo, pbo")
		("bucket-size,b",   opts::value<int>(&maxBucketSize)->default_value(1),
		                    "set (maximum) bucket size")
		("num-chains,c",    opts::value<int>(&nChains)->default_value(1),
		                    "set number of bucket chains")
		("num-samples,s",   opts::value<int>(&nSamples)->default_value(0),
		                    "set number of samples to draw")
		("sample-steps,S",  opts::value<int>(&nStepsPerSample)->default_value(1),
		                    "set number of steps per sample")
		("burnin-steps,B",  opts::value<int>(&nBurnInSteps)->default_value(0),
		                    "set number of burn-in steps")
		("burnout-steps",   opts::value<int>(&nBurnOutSteps)->default_value(0),
		                    "set number of burn-out steps")
		("seed",            opts::value<unsigned int>(&rngSeed)->default_value(time(0)),
		                    "set seed for random number generator")
		("input-file,i",    opts::value<string>(&inFilename)->default_value(""),
		                    "set input file for data")
		("score-file",      opts::value<string>(&scoreFilename)->default_value(""),
		                    "set score file (for output if data file given and for input otherwise)")
		("output-file,o",   opts::value<string>(&outFilename)->default_value("-"),
		                    "set output file for feature probabilities")
		("margin-file",     opts::value<string>(&marginOutFilename)->default_value(""),
		                    "set output file for marginal probabilities")
		("log-file",        opts::value<string>(&logFilename)->default_value(""),
		                    "set log file to write statistics about computation")
	        ("prior-file",      opts::value<string>(&priorFilename)->default_value(""),
	                            "set input file for prior probabilities of arcs")
	        ("print-samples",   "output instead all p-value samples")
	        ("beta",            opts::value<double>(&beta)->default_value(1),
		                    "set prior temperature, beta (default=1)")
	        ("order-prior",     "Use priors as ancestor relationships (order priors) instead of parent priors")
		;

	opts::positional_options_description pdesc;
	pdesc.add("input-file", 1);
	pdesc.add("output-file", 1);
	
	opts::variables_map vm;
	
	try {
		opts::store(opts::command_line_parser(argc, argv).options(desc).positional(pdesc).run(), vm);
		opts::notify(vm);
	} catch (opts::error& err) {
		logger.println(-1, "Error: ", err.what());
		logger.println(-1, "Aborting.");
		return 1;
	}
	
	if (vm.count("help")) {
		logger.println(-1, "BEANDiscoPrior - Bayesian Exact and Approximate Network Discovery with Priors");
		logger.println(-1, "Version " BEAND_VERSION_STRING);
		logger.println(-1);
		logger.println(-1, "Usage:");
		logger.printfln(-1, "  %s [options] [infile [outfile]]", argv[0]);
		logger.println(-1);
		logger.println(-1, desc);
		return 1;
	}
	
	bool exact = vm.count("exact");
	bool testConv = vm.count("test-conv");
	
	bool printSamples = vm.count("print-samples");
	bool orderPrior = vm.count("order-prior");
	
	if (vm.count("quiet"))
		logger.setVerbosity(-1);
	else
		logger.setVerbosity(verbosity);
	
	
	logger.println(2, "Parameters:");
	logger.printfln(2, "  computation type = %s", exact ? "Exact" : "MCMC");
	logger.printfln(2, "  data file = %s", inFilename.c_str());
	logger.printfln(2, "  prior file = %s", priorFilename.c_str());
	logger.printfln(2, "  outFilename = %s", outFilename.c_str());
	logger.printfln(2, "  nVariables = %d", nVariables);
	logger.printfln(2, "  nDataSamples = %d", nDataSamples);
	logger.printfln(2, "  maxIndegree = %d", maxIndegree);
	logger.printfln(2, "  nMcmcSamples = %d", nSamples);
	logger.printfln(2, "  maxBucketSize = %d", maxBucketSize);
	
	
	// open log stream for statistics
	ofstream logFile;
	ostream logStream(0);
	if (!logFilename.empty()) {
		if (logFilename == "-") {
			logStream.rdbuf(cout.rdbuf());
		} else {
			logFile.open(logFilename.c_str());
			if (!logFile) {
				logger.printfln(-1, "Error: couldn't open file '%s' for writing.", logFilename.c_str());
				return 1;
			}
			logStream.rdbuf(logFile.rdbuf());
		}
	}
	logStream.setf(ios::fixed);
	logStream.precision(2);
	
	logStream << "data_file = " << inFilename << endl;
	logStream << "score_file = " << scoreFilename << endl;
	logStream << "prior_file = " << priorFilename << endl;
	logStream << "output_file = " << outFilename << endl;
	logStream << "margin_output_file = " << marginOutFilename << endl;
	logStream << "variables = " << nVariables << endl;
	logStream << "rows = " << nDataSamples << endl;
	logStream << "maximum_indegree = " << maxIndegree << endl;
	logStream << "bucket_size = " << maxBucketSize << endl;
	logStream << "chains = " << nChains << endl;
	logStream << "samples = " << nSamples << endl;
	logStream << "steps = " << nStepsPerSample << endl;
	logStream << "burnin_steps = " << nBurnInSteps << endl;
	logStream << "burnout_steps = " << nBurnOutSteps << endl;
	logStream << "seed = " << rngSeed << endl;

	// initialize rng
	rng.seed(rngSeed);
	
	// start global timer
	Timer timer;
	timer.start();
	
	// map for local scores
	ParentsetMap<Real>* scores;
	// map for arc priors
	ArcMap<Real> *priors = NULL;
	// map for local priors (per parentset)
	ParentsetMap<Real>* parentPriors;
	
	// if data file given, read the data and compute the scores
	if (!inFilename.empty()) {
		if (maxIndegree <= 0) {
			logger.println(-1, "Error: The maximum in-degree not given.");
			logger.println(-1, "Aborting.");
			return 1;
		}
		
		logger.println(1, "Reading data...");
		Data data;
		istream inStream(0);
		ifstream inFile;
		if (inFilename == "-") {
			inStream.rdbuf(cin.rdbuf());
		} else {
			inFile.open(inFilename.c_str());
			if (!inFile) {
				logger.printfln(-1, "Error: Could not open file '%s' for reading.", inFilename.c_str());
				return 1;
			}
			inStream.rdbuf(inFile.rdbuf());
		}
		try {
			if (nDataSamples > 0 || nVariables > 0)
				data.read(inFile, nVariables, nDataSamples);
			else
				data.read(inStream);
		} catch (Exception& e) {
			logger.printfln(-1, "Error: While reading data file '%s': %s", inFilename.c_str(), e.what());
			return 1;
		}
		if (inFile.is_open())
			inFile.close();
		
		nVariables = data.nVariables;
		
		// if prior file given, read it in
		if (!priorFilename.empty()) {
		  priors = new ArcMap<Real>(nVariables);
		  ifstream inPriorFile;
		  inPriorFile.open(priorFilename.c_str());
		  if (!inPriorFile) {
		    logger.printfln(-1, "Error: Could not open file '%s' for reading.", priorFilename.c_str());
		    return 1;
		  }
		  try {
		    readPriorFile(inPriorFile, priors);
		  } catch (Exception& e) {
		    logger.printfln(-1, "Error: While reading prior file '%s': %s", priorFilename.c_str(), e.what());
		    return 1;
		  }
		  if (inPriorFile.is_open())
		    inPriorFile.close();
		}

		// If priors are parent priors, calculate priors
		if (!priorFilename.empty() && !orderPrior) {
		  logger.println(1, "Computing priors...");
		  Timer priorTimer; priorTimer.start();
		  parentPriors = new ParentsetMap<Real>(nVariables, maxIndegree);
		  computePriors(priors, beta, *parentPriors);
		  double priorTime = priorTimer.elapsed();
		  logStream << "prior_computation_time = " << priorTime << endl;
		  logger.printfln(1, "  Elapsed %.2f s.", priorTime);
		} else {
		  parentPriors = NULL;
		}
		
		logger.println(1, "Computing scores...");
		Timer scoreTimer; scoreTimer.start();
		logger.println(2, "  Allocating a score map...");
		scores = new ParentsetMap<Real>(nVariables, maxIndegree);
		logger.println(2, "  Computing the actual scores...");
		computeScores(data, *scores, parentPriors);
		double scoreTime = scoreTimer.elapsed();
		logStream << "score_computation_time = " << scoreTime << endl;
		logger.printfln(1, "  Elapsed %.2f s.", scoreTime);
		
		// optionally write scores to file
		if (!scoreFilename.empty()) {	
			logger.println(1, "Writing scores...");
			if (scoreFilename == "-") {
				writeScores(cout, *scores);
			} else {
				ofstream scoreFile(scoreFilename.c_str());
				if (!scoreFile) {
					logger.printfln(-1, "Error: Could not open file '%s' for writing.",
							scoreFilename.c_str());
					return 1;
				}
				writeScores(scoreFile, *scores);
				scoreFile.close();
			}
		}
	
	// if data file not given, read the scores
	} else if (!scoreFilename.empty()) {
		logger.println(1, "Reading scores...");
		if (scoreFilename == "-") {
			scores = readScores(cin);
		} else {
			ifstream scoreFile(scoreFilename.c_str());
			if (!scoreFile) {
				logger.printfln(-1, "Error: Could not open file '%s' for reading.",
						scoreFilename.c_str());
				return 1;
			}
			try {
				scores = readScores(scoreFile);
			} catch (Exception& e) {
				logger.printfln(-1, "Error: While reading score file '%s': %s",
						scoreFilename.c_str(), e.what());
				return 1;
			}
			scoreFile.close();
		}
		nVariables = scores->nNodes;
		maxIndegree = scores->maxParents;

		// if prior file given, read it in
		if (!priorFilename.empty()) {
		  priors = new ArcMap<Real>(nVariables);
		  ifstream inPriorFile;
		  inPriorFile.open(priorFilename.c_str());
		  if (!inPriorFile) {
		    logger.printfln(-1, "Error: Could not open file '%s' for reading.", priorFilename.c_str());
		    return 1;
		  }
		  try {
		    readPriorFile(inPriorFile, priors);
		  } catch (Exception& e) {
		    logger.printfln(-1, "Error: While reading prior file '%s': %s", priorFilename.c_str(), e.what());
		    return 1;
		  }
		  if (inPriorFile.is_open())
		    inPriorFile.close();

		  // if priors are parent priors, re-compute scores
		  if (!orderPrior) {
		    logger.println(1, "Computing priors...");
		    Timer priorTimer; priorTimer.start();
		    parentPriors = new ParentsetMap<Real>(nVariables, maxIndegree);
		    computePriors(priors, beta, *parentPriors);
		    double priorTime = priorTimer.elapsed();
		    logStream << "prior_computation_time = " << priorTime << endl;
		    logger.printfln(1, "  Elapsed %.2f s.", priorTime);
		    
		    logger.println(1, "Computing scores with priors...");
		    Timer scoreTimer; scoreTimer.start();
		    logger.println(2, "  Adjusting scores for priors...");
		    adjustScores(*scores, parentPriors);
		    double scoreTime = scoreTimer.elapsed();
		    logStream << "score_computation_time = " << scoreTime << endl;
		    logger.printfln(1, "  Elapsed %.2f s.", scoreTime);
		  } else {
		    parentPriors = NULL;
		  }
		}
		// either way, if priorfilename was given, priors will be assigned (can be used as order priors)
	
	// complain if neither data nor scores was given
	} else {
		logger.println(-1, "Error: either data or score file should be provided.");
		logger.printfln(-1, "Type '%s --help' to see help message.", argv[0]);
		return 1;
	}
	
	// weight scores
	logger.println(1, "Weighting scores...");
	weightScores(*scores);
	
	// open result stream for writing
	ofstream resFile;
	ostream resStream(0);
	if (outFilename == "-") {
		resStream.rdbuf(cout.rdbuf());
	} else {
		resFile.open(outFilename.c_str());
		if (!resFile) {
			logger.printfln(-1, "Error: couldn't open file '%s' for writing.", outFilename.c_str());
			return 1;
		}
		resStream.rdbuf(resFile.rdbuf());
	}
	
	// open margin stream if needed
	ofstream marginFile;
	ostream marginStream(0);
	if (marginOutFilename == "-") {
		marginStream.rdbuf(cout.rdbuf());
	} else if (marginOutFilename != "") {
		marginFile.open(marginOutFilename.c_str());
		if (!marginFile) {
			logger.printfln(-1, "Error: couldn't open file '%s' for writing.", marginOutFilename.c_str());
			return 1;
		}
		marginStream.rdbuf(marginFile.rdbuf());
	}
	marginStream.setf(std::ios::scientific);
	marginStream.precision(10);

	
	// actual computation
	
	// version for bucket order
	if (poType == "bo") {
		if (nChains != 1)
			logger.println(0, "Warning: partial order type 'bo' does not support multiple chains.");
		BucketOrderFamily pof(nVariables, maxBucketSize);
		if (exact) {
			// Exact arc probabilities
		        if (orderPrior) {
			        printExact(resStream, logStream, pof, *scores, priors, beta);
			} else {
			        printExact(resStream, logStream, pof, *scores, NULL, 1);
			}
		} else {
			// MCMC sampling => approximated probabilities
		        if (orderPrior) {
			        printMCMC(resStream, marginStream, logStream, pof, *scores, nBurnInSteps,
					  nSamples, nStepsPerSample, nBurnOutSteps, printSamples,
					  priors, beta);
			} else {
			        printMCMC(resStream, marginStream, logStream, pof, *scores, nBurnInSteps,
					  nSamples, nStepsPerSample, nBurnOutSteps, printSamples,
					  NULL, 1);
			}
		}
	}
	// version for parallel bucket order
	//	else if (poType == "pbo") {
	//	ParBucketOrderFamily pof(nVariables, maxBucketSize, nChains);
	//	if (exact) {
	//		// Exact arc probabilities, TODO: implement
	//		logger.println(-1, "Exact computation is currently not supported for Parallel Bucket Orders.");
	//		return 1;
	//	} else {
	//		// MCMC sampling => approximated probabilities
	//		printMCMC(resStream, marginStream, logStream, pof, *scores, nBurnInSteps,
	//			  nSamples, nStepsPerSample, nBurnOutSteps, printSamples,
	//			  NULL, 1); // priors not supported
	//	}
	//}
	else {
		logger.println(-1, "Error: unknown partial order type.");
		return 1;
	}
	
	// actual computation ends
	
	
	// print timing
	double elapsedTime = timer.elapsed();	
	logger.printfln(0, "Elapsed %g s total.", elapsedTime);
	logStream << "time = " << elapsedTime << endl;
	
	// release resources
	delete scores;
	
	if (resFile.is_open())
		resFile.close();

	if (marginFile.is_open())
		marginFile.close();
	
	if (logFile.is_open())
		logFile.close();
	
	return 0;
}


