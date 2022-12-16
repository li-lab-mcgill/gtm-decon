#include "JCVB0.h"
#include "GeneParams.h"

using namespace std;
using namespace arma;

namespace bm = boost::math;

void JCVB0::updateAlpha() {

	// update alpha by Minka's fixed point iteration
	rowvec alpha_numer = zeros<rowvec>(K);

	double alpha_denom = 0;

	for(vector<SampleRS>::iterator RSSample = trainRSSamples->begin(); RSSample != trainRSSamples->end(); RSSample++) {

		alphaSum = 0;

		for (unordered_map<int, double>::iterator iter = RSSample->topicMap.begin(); iter != RSSample->topicMap.end(); iter++){
			int i = std::distance(RSSample->topicMap.begin(), iter);
			alpha_numer(iter->first) += bm::digamma(RSSample->metagene(i) + alpha(iter->first));
			alphaSum += alpha(iter->first);
		}

		alpha_denom += bm::digamma(accu(RSSample->metagene) + alphaSum);
	}

	for(int k=0; k<K; k++) {
		alpha_numer(k) -= D_train * bm::digamma(alpha(k));
	}

	alpha_denom -= D_train * bm::digamma(alphaSum);

	// having 10 and 100 imply a gamma(2,20) with mean at 0.1
	alpha = alpha % alpha_numer/(10+alpha_denom);

//	cout << "alpha: " << alpha << endl;

	alphaSum = accu(alpha);
}


void JCVB0::updateBeta() {

	unordered_map<int, double> beta_denom = unordered_map<int, double>();

	for(unordered_map<int,double>::iterator iter=betaSum.begin(); iter!=betaSum.end(); iter++) {

		int t = iter->first;

		beta_denom[t] = 0;

		betaSum[t] = 0;
		sumLogGammaBeta[t] = 0;
	}

	for(int k=0; k<K; k++) {

		unordered_map<int, double> beta_denom_k = unordered_map<int, double>();

		for(unordered_map<int,double>::iterator iter=betaSum.begin(); iter!=betaSum.end(); iter++) {

			beta_denom_k[iter->first] = 0;
		}

		for(unordered_map<pair<int,int>, GeneParams*>::iterator iter = geneParams.begin(); iter != geneParams.end(); iter++) {

			int t = iter->first.first;

			beta_denom_k[t] += iter->second->beta + iter->second->phi(k);
		}

		for(unordered_map<int,double>::iterator iter=betaSum.begin(); iter!=betaSum.end(); iter++) {

			int t = iter->first;

			beta_denom[t] += bm::digamma(beta_denom_k[t]);
		}

		beta_denom_k.clear();
	}



	unordered_map<int, double> beta_denom_t = unordered_map<int, double>();

	for(unordered_map<int,double>::iterator iter=betaSum.begin(); iter!=betaSum.end(); iter++) {

		beta_denom_t[iter->first] = 0;
	}


	for(unordered_map<pair<int,int>, GeneParams*>::iterator iter = geneParams.begin(); iter != geneParams.end(); iter++) {

		int t = iter->first.first;

		beta_denom_t[t] += iter->second->beta;
	}


	for(unordered_map<int,double>::iterator iter=betaSum.begin(); iter!=betaSum.end(); iter++) {

		int t = iter->first;

		beta_denom[t] -= K * bm::digamma(beta_denom_t[t]);
	}

	beta_denom_t.clear();


	for(unordered_map<pair<int,int>, GeneParams*>::iterator iter = geneParams.begin(); iter != geneParams.end(); iter++) {

		int t = iter->first.first;

		double beta_numer = 0;

		for(int k=0; k<K; k++) {

			beta_numer += bm::digamma(iter->second->beta + iter->second->phi(k));
		}

		beta_numer -= K*bm::digamma(iter->second->beta);

		iter->second->beta = (1 + iter->second->beta * beta_numer) / (100 + beta_denom[t]); // fixed point update

		if(iter->second->beta==0) {
			throw::runtime_error("iter->second->beta become zeros");
		}

		betaSum[t] += iter->second->beta;

		sumLogGammaBeta[t] += lgamma(iter->second->beta);
	}

	beta_denom.clear();
}


void JCVB0::updateHyperParams() {

	updateAlpha();

//	cout << "alpha updated" << endl;

	updateBeta();

//	cout << "beta updated" << endl;

}


void JCVB0::updateParamSums() {

	// update topicCountsPerGeneType
	unordered_map<int, rowvec> topicCountsPerGeneType_new = unordered_map<int, rowvec>();

	for(unordered_map<int, rowvec>::iterator iter=topicCountsPerGeneType.begin(); iter!=topicCountsPerGeneType.end(); iter++)
		topicCountsPerGeneType_new[iter->first] = zeros<rowvec>(K);

	for(unordered_map<pair<int,int>, GeneParams*>::iterator iter = geneParams.begin(); iter != geneParams.end(); iter++)
		topicCountsPerGeneType_new[iter->first.first] += iter->second->phi;

	for(unordered_map<int, rowvec>::iterator iter=topicCountsPerGeneType.begin(); iter!=topicCountsPerGeneType.end(); iter++)
		topicCountsPerGeneType[iter->first] = topicCountsPerGeneType_new[iter->first];

	topicCountsPerGeneType_new.clear();
}


void JCVB0::normalizeParams() {

	//cout << "num topics per cell type:  " << JCVB0::numTopicsCellType << endl;
	for(unordered_map<pair<int,int>, GeneParams*>::iterator iter = geneParams.begin(); iter != geneParams.end(); iter++) {

		iter->second->phi_normalized = (iter->second->beta + iter->second->phi)/

				(betaSum[iter->first.first] + topicCountsPerGeneType[iter->first.first]);
	}
}


void JCVB0::updateMainParams() {
	int i;

	for(unordered_map<pair<int,int>, GeneParams*>::iterator iter = geneParams.begin(); iter != geneParams.end(); iter++) {
		iter->second->phi.zeros(); // 1 x K
	}

	// main update
	for(vector<SampleRS>::iterator RSSample = trainRSSamples->begin(); RSSample != trainRSSamples->end(); RSSample++) {

		vector<int> mapping(RSSample->Ki,0);
		for (unordered_map<int, double>::iterator iter = RSSample->topicMap.begin(); iter != RSSample->topicMap.end(); iter++){
			i = std::distance(RSSample->topicMap.begin(), iter);
			mapping[i] = iter->first;
		}

		// update gene params
		// iterate latents for ONLY the observed RSSample-gene
		for(unordered_map<pair<int, int>, int>::iterator iter = RSSample->geneDict.begin(); iter != RSSample->geneDict.end(); iter++) {

			pair<int,int> geneId = iter->first;

			for (unordered_map<int, double>::iterator iter = RSSample->topicMap.begin(); iter != RSSample->topicMap.end(); iter++){
				int i = std::distance(RSSample->topicMap.begin(), iter);
				geneParams[geneId]->phi[iter->first] += iter->second * RSSample->gamma[geneId][i];
			}
		}
	}

	updateParamSums();

	normalizeParams();

}


// sorted id for outputs
void JCVB0::sortID() {

	// sorted gene id
	vector<int> typeIds;

	for(unordered_map<pair<int,int>, GeneParams*>::iterator genePar = geneParams.begin(); genePar != geneParams.end(); genePar++) {

		int t = genePar->first.first;
		int w = genePar->first.second;

		if(geneIds.find(t)==geneIds.end()) {

			typeIds.push_back(t);

			geneIds[t].push_back(w);

		} else {

			geneIds[t].push_back(w);
		}
	}

	sort(typeIds.begin(), typeIds.end());

	for(vector<int>::iterator iter=typeIds.begin(); iter!=typeIds.end(); iter++) {

		sort(geneIds[*iter].begin(), geneIds[*iter].end());
	}
}


void JCVB0::showParams() {

	cout << endl << "<<Gene parameters>>" << endl;

	for(map<int, vector<int>>::iterator t = geneIds.begin(); t != geneIds.end(); t++) {

		int typeId = t->first;

		for(vector<int>::iterator w = geneIds[t->first].begin(); w != geneIds[t->first].end(); w++) {

			int geneId = *w;

			cout << "--------------------" << endl;
			printf("typeId %d; geneId %d:\n", typeId, geneId);
			cout << "--------------------" << endl;

			cout << geneParams[make_pair(typeId, geneId)]->phi;

			cout << "--------------------" << endl;
			cout << "beta: " << geneParams[make_pair(typeId, geneId)]->beta << endl;
		}
	}

	cout << endl << "--------------------" << endl;
	cout << "topicCountsPerGeneType: " << endl;
	for(map<int, vector<int>>::iterator t = geneIds.begin(); t != geneIds.end(); t++) {
		cout << topicCountsPerGeneType[t->first];
	}

	cout << "--------------------" << endl;
	cout << "alpha: " << endl << alpha << endl;
	cout << "--------------------" << endl;


	cout << endl << "--------------------" << endl << endl;
}



