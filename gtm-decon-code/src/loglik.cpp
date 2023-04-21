#include "JCVB0.h"

using namespace std;
using namespace arma;

// marginal likelihood estimated by the
// variational evidence lower bound (elbo)
// NB: bc in JCVB0, variance is not estimated
// this entails not including teh qlog(q) term
// in the original ELBO
double JCVB0::trainLogLik() {
	int i;
	int D_train = trainRSSamples->size();

	//cout << "loglik 15" << endl;

	// theta elbo
	double elbo = D_train * (lgamma(alphaSum) - accu(lgamma(alpha)));

	//cout << "loglik 20" << endl;

	for(std::vector<SampleRS>::iterator RSSamplej = trainRSSamples->begin(); RSSamplej!=trainRSSamples->end(); RSSamplej++) {

		rowvec alpha_i = zeros<rowvec>(RSSamplej->Ki);
		for (unordered_map<int, double>::iterator iter = RSSamplej->topicMap.begin(); iter != RSSamplej->topicMap.end(); iter++){
			i = std::distance(RSSamplej->topicMap.begin(), iter);
			alpha_i[i] = alpha[iter->first] + iter->second;
		}

		elbo += accu(lgamma(alpha_i + RSSamplej->metagene)) - lgamma(alphaSum + accu(RSSamplej->metagene));
	}

	//cout << "loglik 27" << endl;

	//	cout << "theta elbo: " << elbo << endl;

	// phi elbo
	double betaPrior = 0;

	for(unordered_map<int,double>::iterator iter=betaSum.begin(); iter!=betaSum.end(); iter++) {

		int t = iter->first;

		//		cout << "betaSum[t]: " << betaSum[t] << endl;
		//		cout << "sumLogGammaBeta[t]: " << sumLogGammaBeta[t] << endl;

		betaPrior += lgamma(betaSum[t]) - sumLogGammaBeta[t];
	}

	//cout << "loglik 44" << endl;

	double betaLike = 0;

	for(unordered_map<pair<int,int>, GeneParams*>::iterator genePar = geneParams.begin(); genePar != geneParams.end(); genePar++) {

		betaLike += accu(lgamma(genePar->second->beta + genePar->second->phi));
	}

	//cout << "loglik 53" << endl;

	double betaNorm = 0;

	for(unordered_map<int,double>::iterator iter=betaSum.begin(); iter!=betaSum.end(); iter++) {

		int t = iter->first;

		for(int k = 0; k < K; k++) {

			betaNorm += lgamma(betaSum[t] + topicCountsPerGeneType[t](k));
		}
	}

	//cout << "loglik 67" << endl;

	elbo += K * betaPrior + betaLike - betaNorm;

	//	cout << "phi elbo: " << K * betaPrior + betaLike - betaNorm << endl;
	//	cout << "betaPrior: " << betaPrior << endl;
	//	cout << "betaLike: " << betaLike << endl;
	//	cout << "betaNorm: " << betaNorm << endl;
	//	throw::runtime_error("");

	vector<SampleRS>::iterator RSSample0 = trainRSSamples->begin();

	int j = 0;

	vec logEq = zeros<vec>(trainRSSamples->size());

	// Eq[logq(z,h)]
#pragma omp parallel for shared(j)
	for(j=0; j < (int) trainRSSamples->size(); j++) {

		vector<SampleRS>::iterator RSSamplej = RSSample0 + j;

		for(unordered_map<pair<int, int>, int>::iterator iter = RSSamplej->geneDict.begin(); iter != RSSamplej->geneDict.end(); iter++) {

			pair<int,int> geneId = iter->first;

			logEq(j) += accu(RSSamplej->gamma[geneId] % log(RSSamplej->gamma[geneId]));
		}

		// free up the memory allocated for the RSSample gamma hash
		RSSamplej->gamma.clear();

	}

	//cout << "loglik 171" << endl;

	elbo -= accu(logEq);

	return elbo/C_train;
}


// first predict theta using x% of test terms
// then predict the 1-x% target terms using theta
double JCVB0::predictLogLik() {

	int j = 0;

	vec llk_vec = zeros<vec>(testRSSamples->size());

	vector<SampleRS>::iterator RSSample0 = testRSSamples->begin();

#pragma omp parallel for shared(j)
	for(j=0; j < (int) testRSSamples->size(); j++) {

		vector<SampleRS>::iterator RSSample = RSSample0 + j;

//		cout << "testRSSample " << RSSample->RSSampleId << endl;

		inferRSSampleParams(RSSample, inferTestRSSampleMetagene_maxiter_intermediateRun);

//		cout << "done" << endl;

		// infer RSSample target genes
		for(unordered_map<pair<int,int>, int>::iterator iter = RSSample->geneDict.begin(); iter != RSSample->geneDict.end(); iter++) {

			pair<int,int> geneId = iter->first;

			if(RSSample->isTestGene[geneId]) {

				double llk_ij = accu(RSSample->metagene_normalized % geneParams[geneId]->phi_normalized);

				llk_vec(j) += iter->second * log(llk_ij);
			}
		}

//		cout << "gene loglik done" << endl;
	}

	return accu(llk_vec)/C_tar;
}



// DEBUG oNLY
rowvec JCVB0::trainLogLik_breakdowns() {

	rowvec elbo_bkdw = zeros<rowvec>(4+(!mar));

	int D_train = trainRSSamples->size();
	int i = 0;

	// theta elbo
	double thetaELBO = D_train * (lgamma(alphaSum) - accu(lgamma(alpha)));

	for(std::vector<SampleRS>::iterator RSSamplej = trainRSSamples->begin(); RSSamplej!=trainRSSamples->end(); RSSamplej++) {

		thetaELBO += accu(lgamma(alpha + RSSamplej->metagene)) - lgamma(alphaSum + accu(RSSamplej->metagene));
	}

	elbo_bkdw(i) = thetaELBO;
	i++;


	// phi elbo
	double betaPrior = 0;

	for(unordered_map<int,double>::iterator iter=betaSum.begin(); iter!=betaSum.end(); iter++) {

		int t = iter->first;

		betaPrior += lgamma(betaSum[t]) - sumLogGammaBeta[t];
	}

	double betaLike = 0;

	for(unordered_map<pair<int,int>, GeneParams*>::iterator genePar = geneParams.begin(); genePar != geneParams.end(); genePar++) {

		betaLike += accu(lgamma(genePar->second->beta + genePar->second->phi));
	}

	double betaNorm = 0;

	for(unordered_map<int,double>::iterator iter=betaSum.begin(); iter!=betaSum.end(); iter++) {

		int t = iter->first;

		for(int k = 0; k < K; k++) {

			betaNorm += lgamma(betaSum[t] + topicCountsPerGeneType[t](k));
		}
	}

	elbo_bkdw(i) = K * betaPrior + betaLike - betaNorm;
	i++;


	double elbo_logq = 0;

	// Eq[logq(z,h)]
	for(std::vector<SampleRS>::iterator RSSamplej = trainRSSamples->begin(); RSSamplej!=trainRSSamples->end(); RSSamplej++) {

		for(unordered_map<pair<int, int>, int>::iterator iter = RSSamplej->geneDict.begin(); iter != RSSamplej->geneDict.end(); iter++) {

			pair<int,int> geneId = iter->first;

			elbo_logq += accu(RSSamplej->gamma[geneId] % log(RSSamplej->gamma[geneId]));
		}

	}

	elbo_bkdw(i) = elbo_logq;

	return elbo_bkdw/C_train;
}



























