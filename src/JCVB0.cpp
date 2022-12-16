#include "JCVB0.h"
#include "GeneParams.h"

using namespace std;
using namespace arma;

void JCVB0::initialize(
		unordered_map<int,int> geneCnt,
		int numOfTopics, int numTopicsPerCellType, int maxiter,
		unordered_map<pair<int,int>, GeneParams*> geneParamsMap
		)
{
	K = numOfTopics;
	numTopicsCellType = numTopicsPerCellType;

	numOfGenes = unordered_map<int,int>();

	C_train=0;
	C_tar=0;

	iterations = maxiter;

	svi = false;
	mar = false;
	imp = false;

	trainRSSamples = new vector<SampleRS>();
	testRSSamples = new vector<SampleRS>();
	imputeTargetRSSamples = new vector<SampleRS>();

	// initialize alpha
	alpha = zeros<rowvec>(K);
	alpha.fill(0.1);
	alphaSum = accu(alpha);

	geneParams = geneParamsMap;

	betaSum = unordered_map<int,double>();

	sumLogGammaBeta = unordered_map<int,double>();

	topicCountsPerGeneType = unordered_map<int, rowvec>();

	inferTestRSSampleMetagene_maxiter_intermediateRun = 2;

	// update by gene types
	for(unordered_map<int,int>::iterator iter=geneCnt.begin(); iter!=geneCnt.end(); iter++) {

		int t = iter->first;

		numOfGenes[t] = geneCnt[t];

		topicCountsPerGeneType[t] = zeros<rowvec>(K);

		betaSum[t] = 0;

		sumLogGammaBeta[t] = 0;
	}

	for(unordered_map<pair<int,int>, GeneParams*>::iterator iter = geneParams.begin(); iter != geneParams.end(); iter++) {

		int t = iter->first.first;

		betaSum[t] += iter->second->beta;

		sumLogGammaBeta[t] += lgamma(iter->second->beta);
	}

	updateParamSums();

	normalizeParams();

	sortID();
}


JCVB0::~JCVB0() {

	geneParams.clear();

	betaSum.clear();
	sumLogGammaBeta.clear();

	for(vector<SampleRS>::iterator RSSample = trainRSSamples->begin(); RSSample != trainRSSamples->end(); RSSample++) {
		RSSample->~SampleRS();
	}
	trainRSSamples->clear();

	for(vector<SampleRS>::iterator RSSample = testRSSamples->begin(); RSSample != testRSSamples->end(); RSSample++) {
		RSSample->~SampleRS();
	}
	testRSSamples->clear();
}


// Update variational parameters by JCVB0
void JCVB0::train(bool updateHyper) {

	cout << "E-step" << endl;

	// E-step
	inferAllRSSampleParams();

	cout << "M-step" << endl;

	// M-step
	updateMainParams();

	cout << "EB-step" << endl;

	if(updateHyper) updateHyperParams();
}


// E-step
void JCVB0::inferAllRSSampleParams(int trainIter) {
	int j = 0;
	vector<SampleRS>::iterator RSSample0 = trainRSSamples->begin();

#pragma omp parallel for shared(j)
	for(j=0; j < (int) trainRSSamples->size(); j++) {
		inferRSSampleParams(RSSample0 + j, trainIter);
	}
}



void JCVB0::inferTestRSSampleParams_finalRun() {
	normalizeParams();
	int j = 0;
	vector<SampleRS>::iterator RSSample0 = testRSSamples->begin();

#pragma omp parallel for shared(j)
	for(j=0; j < (int) testRSSamples->size(); j++) {
		inferRSSampleParams(RSSample0 + j, inferTestRSSampleMetagene_maxiter_finalRun);
	}
}

// E-step
void JCVB0::inferRSSampleParams(vector<SampleRS>::iterator RSSample, int maxiter) {
	int i;
	vector<int> mapping(RSSample->Ki,0);
	rowvec alpha_i = zeros<rowvec>(RSSample->Ki);
	rowvec phi_i = zeros<rowvec>(RSSample->Ki);
	rowvec ztot_i = zeros<rowvec>(RSSample->Ki);

	//rowvec prior_prob_i = zeros<rowvec>(RSSample->Ki);

	for (unordered_map<int, double>::iterator iter = RSSample->topicMap.begin(); iter != RSSample->topicMap.end(); iter++){
		i = std::distance(RSSample->topicMap.begin(), iter);
		mapping[i] = iter->first;
		alpha_i[i] = alpha[iter->first] + iter->second;
	}

	for(unordered_map<pair<int, int>, int>::iterator iter = RSSample->geneDict.begin(); iter != RSSample->geneDict.end(); iter++) {

		if(!RSSample->isTestGene[iter->first]) {

			RSSample->gamma[iter->first] = zeros<rowvec>(RSSample->Ki); // 1 x K
		}
	}

	double diff = 1;
	int iter = 1;
	while(diff > 1e-2 && iter <= maxiter) {

		rowvec RSSampleMetagene_j_prev = RSSample->metagene;

		// iterate latents for ONLY the observed RSSample-gene
		for(unordered_map<pair<int,int>, int>::iterator RSSampleiter = RSSample->geneDict.begin(); RSSampleiter != RSSample->geneDict.end(); RSSampleiter++) {

			pair<int,int> geneId = RSSampleiter->first;

			for(i=0;i<RSSample->Ki;i++){
				phi_i[i] = geneParams[geneId]->phi[mapping[i]];
				ztot_i[i] = topicCountsPerGeneType[geneId.first][mapping[i]];				
			}

			if(!RSSample->isTestGene[geneId]) {

				// removing the current RSSample-gene token
				rowvec geneWeight_ij = RSSampleiter->second * RSSample->gamma[geneId];

				RSSample->gamma[geneId] = (RSSample->metagene - geneWeight_ij + alpha_i) %
						(phi_i - geneWeight_ij + geneParams[geneId]->beta) /
						(ztot_i - geneWeight_ij + betaSum[geneId.first]);

				if(any(RSSample->gamma[geneId] < 0)) {

					if(RSSample->isTestRSSample || svi || imp) {

						RSSample->gamma[geneId] = (RSSample->metagene + alpha_i) %
								(phi_i + geneParams[geneId]->beta) /
								(ztot_i + betaSum[geneId.first]);

					} else {

						// DEBUG S impossible for training unless there is a bug

						throw::runtime_error("inferRSSampleParams(): RSSample->gamma[geneId] has negative values for training RSSamples");
						// DEBUG E
					}
				}

				// not needed because normalizing over gene tokens produce much better results (below)
				RSSample->gamma[geneId] = RSSample->gamma[geneId]/accu(RSSample->gamma[geneId]);
			}
		}

		RSSample->metagene.zeros();

		// infer RSSample metagene
		// update gene params
		// iterate latents for ONLY the observed RSSample-gene
		for(unordered_map<pair<int, int>, int>::iterator iter = RSSample->geneDict.begin(); iter != RSSample->geneDict.end(); iter++) {

			pair<int,int> geneId = iter->first;

			if(!RSSample->isTestGene[geneId]) {

				RSSample->metagene += iter->second * RSSample->gamma[geneId];
			}
		}

		diff = accu(abs(RSSample->metagene - RSSampleMetagene_j_prev))/RSSample->Ki;

		iter++;

	}

	RSSample->metagene_normalized = alpha_i + RSSample->metagene;
	RSSample->metagene_normalized = RSSample->metagene_normalized/accu(RSSample->metagene_normalized);
}



// E-step
void JCVB0::inferRSSampleParamsUnsupervised(vector<SampleRS>::iterator RSSample, int maxiter) {


	for(unordered_map<pair<int, int>, int>::iterator iter = RSSample->geneDict.begin(); iter != RSSample->geneDict.end(); iter++) {

		if(!RSSample->isTestGene[iter->first]) {

			RSSample->gamma[iter->first] = zeros<rowvec>(K); // 1 x K
		}
	}

	double diff = 1;
	int iter = 1;

	while(diff > 1e-2 && iter <= maxiter) {

		rowvec RSSampleMetagene_j_prev = RSSample->metagene;

		// iterate latents for ONLY the observed RSSample-gene
		for(unordered_map<pair<int,int>, int>::iterator RSSampleiter = RSSample->geneDict.begin(); RSSampleiter != RSSample->geneDict.end(); RSSampleiter++) {

			pair<int,int> geneId = RSSampleiter->first;

			if(!RSSample->isTestGene[geneId]) {

				// removing the current RSSample-gene token
				rowvec geneWeight_ij = RSSampleiter->second * RSSample->gamma[geneId];


				/*cout << "Error after: 1" << endl;
				cout << "RSSampleiter->second" << endl;
				cout << RSSampleiter->second << endl;
				cout << "RSSample->gamma[geneId]" << endl;
				cout << RSSample->gamma[geneId] << endl;
				cout << "geneParams[geneId]->phi" << endl;
				cout << geneParams[geneId]->phi << endl;
				cout << "geneWeight_ij" << endl;
				cout << geneWeight_ij << endl;
				cout << "RSSample->metagene" << endl;
				cout << RSSample->metagene << endl;
				cout << "topicCountsPerGeneType" << endl;
				cout << topicCountsPerGeneType[geneId.first] << endl;*/
				RSSample->gamma[geneId] = (RSSample->metagene - geneWeight_ij + alpha) %
						(geneParams[geneId]->phi - geneWeight_ij + geneParams[geneId]->beta) /
						(topicCountsPerGeneType[geneId.first] - geneWeight_ij + betaSum[geneId.first]);
				//cout << "Error after: 2" << endl;

				if(any(RSSample->gamma[geneId] < 0)) {

					if(RSSample->isTestRSSample || svi || imp) {

						RSSample->gamma[geneId] = (RSSample->metagene + alpha) %
								(geneParams[geneId]->phi + geneParams[geneId]->beta) /
								(topicCountsPerGeneType[geneId.first] + betaSum[geneId.first]);

					} else {

						// DEBUG S impossible for training unless there is a bug
						cout << "RSSample->RSSampleId: " << RSSample->RSSampleId << "; typeId: " <<
								geneId.first << "; geneId: " << geneId.second << endl;

						cout << "iter: " << iter << endl;

						cout << "RSSample->gamma[geneId]: " << RSSample->gamma[geneId] << endl;

						cout << "geneWeight_ij: " << geneWeight_ij << endl;
						cout << "RSSample->metagene - geneWeight_ij + alpha: " << RSSample->metagene - geneWeight_ij + alpha << endl;


						cout << "geneParams[geneId]->phi: " << geneParams[geneId]->phi << endl;
						cout << "geneWeight_ij: " << geneWeight_ij << endl;

						cout << "geneParams[geneId]->beta: " << geneParams[geneId]->beta << endl << endl;


						cout << "(geneParams[geneId]->phi - geneWeight_ij + geneParams[geneId]->beta)" <<
								(geneParams[geneId]->phi - geneWeight_ij + geneParams[geneId]->beta) << endl;

						cout << "(topicCountsPerGeneType[geneId.first] - geneWeight_ij + betaSum[geneId.first]): "<<
								(topicCountsPerGeneType[geneId.first] - geneWeight_ij + betaSum[geneId.first]) << endl;

						throw::runtime_error("inferRSSampleParams(): RSSample->gamma[geneId] has negative values for training RSSamples");
						// DEBUG E
					}
				}

				// not needed because normalizing over gene tokens produce much better results (below)
				RSSample->gamma[geneId] = RSSample->gamma[geneId]/accu(RSSample->gamma[geneId]);
			}
		}

		RSSample->metagene.zeros();

		// infer RSSample metagene
		// update gene params
		// iterate latents for ONLY the observed RSSample-gene
		for(unordered_map<pair<int, int>, int>::iterator iter = RSSample->geneDict.begin(); iter != RSSample->geneDict.end(); iter++) {

			pair<int,int> geneId = iter->first;

			if(!RSSample->isTestGene[geneId]) {

				RSSample->metagene += iter->second * RSSample->gamma[geneId];
			}
		}

		diff = accu(abs(RSSample->metagene - RSSampleMetagene_j_prev))/K;

		iter++;
	}

	RSSample->metagene_normalized = alpha + RSSample->metagene;
	RSSample->metagene_normalized = RSSample->metagene_normalized/accu(RSSample->metagene_normalized);

}
