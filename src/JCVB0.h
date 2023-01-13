#ifndef JCVB0_H_
#define JCVB0_H_

#include <map>
#include <vector>
#include <iostream>
#include <armadillo>

#include "SampleRS.h"
#include "GeneParams.h"


using namespace std;
using namespace arma;

class JCVB0 {
public:

	int K; //Num of Topics
	int numTopicsCellType;

	bool mar;
	bool svi;
	bool imp;

	unordered_map<int,int> numOfGenes;

	int iterations;

	int C_train; // total no of EHR tokens in training RSSample set
	int C_tar; // total no of *target* tokens in test RSSample set

	// key: typeId,geneId
	// value: GenenoParams struct
	unordered_map<pair<int,int>, GeneParams*> geneParams;


	// key: typeId
	// value: 1 x K for K topics
	unordered_map<int, rowvec> topicCountsPerGeneType; // T x K total topic counts per gene type

	// sum of hyperparameters saved for latter computations
	rowvec alpha;
	double alphaSum;


	unordered_map<int, double> betaSum;
	unordered_map<int, double> sumLogGammaBeta;

	int D_train;
	int D_test;

	vector<SampleRS> *trainRSSamples;
	vector<SampleRS> *testRSSamples;

	void initialize(
			unordered_map<int,int> geneCnt,
			int numOfTopics, int numTopicsPerCellType, int maxiter,
			unordered_map<pair<int,int>, GeneParams*> geneParamsMap);

	int inferTestRSSampleMetagene_maxiter_intermediateRun; // during training for monitoring purpose (fast the better)
	int inferTestRSSampleMetagene_maxiter_finalRun; // after training for evaluation purpose (accurate the better)

	void inferRSSampleParams(vector<SampleRS>::iterator RSSample, int maxiter);
	void inferRSSampleParamsUnsupervised(vector<SampleRS>::iterator RSSample, int maxiter);

	void inferAllRSSampleParams(int trainIter=1);
	virtual void updateParamSums();
	virtual void updateMainParams();

	void inferTestRSSampleParams_finalRun();

	// update hyperparams
	virtual void updateAlpha();
	virtual void updateBeta();
	virtual void updateHyperParams();


	map<int, vector<int>> geneIds;
	void sortID(); // sorted id for outputs

	void showParams();
	void normalizeParams();

	virtual ~JCVB0();
	virtual void train(bool updateHyper);
	virtual double trainLogLik();

	virtual rowvec trainLogLik_breakdowns(); // model diagnostics

	double predictLogLik();
};

#endif /* JCVB0_H_ */







