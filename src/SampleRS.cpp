#include <armadillo>

#include "SampleRS.h"

using namespace arma;

SampleRS::SampleRS(int id,
		unordered_map<pair<int, int>, int> geneMap,
		int K)
{
	RSSampleId = id;
	Ki = 0;
	geneDict = geneMap;
	tarGeneFrac = 0.5;
	topicMap = *(new unordered_map<int, double>());
//	rowvec prior;
	metagene = randu<rowvec>(K);
	metagene_normalized = metagene/accu(metagene);

	Cj = 0;

	isTestRSSample = false;

	for(unordered_map<pair<int, int>, int>::iterator iter = geneMap.begin(); iter != geneMap.end(); iter++) {

		isTestGene[iter->first] = false;

		Cj += iter->second;
	}

	Cj_train = Cj;
	Cj_tar = 0;
}

void SampleRS::assignTargetView(int targetViewId) {

	if(geneDict.size() > 0) {

		for(unordered_map<pair<int, int>, int>::iterator iter = geneDict.begin(); iter != geneDict.end(); iter++) {

			if(iter->first.first == targetViewId) {

				isTestGene[iter->first] = true;

				Cj_tar += iter->second;

			} else {

				isTestGene[iter->first] = false;
			}
		}
	}
}


void SampleRS::assignTargetGenenotypes() {

	if(geneDict.size() > 0) {

		vector<pair<int,int>> obsGene;

		for(unordered_map<pair<int, int>, int>::iterator iter = geneDict.begin(); iter != geneDict.end(); iter++) {

			obsGene.push_back(iter->first);

			isTestGene[iter->first] = false;
		}

		std::random_shuffle(obsGene.begin(), obsGene.end());

		vector<pair<int,int> >::const_iterator first_tarGene = obsGene.begin();

		vector<pair<int,int> >::const_iterator last_tarGene = obsGene.begin() + floor(tarGeneFrac * obsGene.size());

		vector<pair<int,int> > tmp(first_tarGene, last_tarGene);

		for(vector<pair<int, int>>::iterator iter = tmp.begin(); iter != tmp.end(); iter++) {

			Cj_tar += geneDict[*iter];

			isTestGene[*iter] = true;
		}

		Cj_train -= Cj_tar;
	}
}

rowvec SampleRS::Decompress(rowvec vec, int K){
	rowvec target = zeros<rowvec>(K);
	for (unordered_map<int, double>::iterator iter = topicMap.begin(); iter != topicMap.end(); iter++){
		int i = std::distance(topicMap.begin(), iter);
		target[iter->first] = vec[i];
	}
	return target;
}


SampleRS::~SampleRS() {

	geneDict.clear();
	gamma.clear();
	isTestGene.clear();
}


