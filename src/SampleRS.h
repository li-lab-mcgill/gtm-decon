#ifndef PATIENT_H_
#define PATIENT_H_

#include <map>
#include <armadillo>

#include "pairkeyhash.h"
#include "GeneParams.h"

using namespace std;
using namespace arma;

class SampleRS {
public:

	int RSSampleId;
	int Ki;

	bool isTestRSSample;

	rowvec metagene; // 1 x K
	rowvec metagene_normalized; // 1 x K
	rowvec prior;

	// key: <typeId,geneId>
	// value: freq
	unordered_map<pair<int, int>, int> geneDict;

	unordered_map<int, double> topicMap;

	// key: <typeId, geneId>
	// value: 1 x K
	unordered_map<pair<int, int>, rowvec> gamma;


	SampleRS(int id,
			unordered_map<pair<int, int>, int> geneMap,
			int K);

	~SampleRS();

	void assignTargetGenenotypes();
	void assignTargetView(int targetViewId);
	rowvec Decompress(rowvec vec, int K);

	double tarGeneFrac;
	int Cj;
	int Cj_tar;
	int Cj_train;

	unordered_map<pair<int,int>, bool> isTestGene;
};

#endif /* PATIENT_H_ */

