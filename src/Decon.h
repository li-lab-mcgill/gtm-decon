#ifndef Decon_H_
#define Decon_H_

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <stdlib.h>
#include <vector>
#include <armadillo>
#include <sstream>
#include <random>
#include <iterator>
#include <algorithm>
#include <omp.h>

using namespace std;
using namespace arma;

#include "GeneParams.h"
#include "SampleRS.h"
#include "JCVB0.h"


class Decon {
public:

	int numOfRSSamples;
	int numOfTopics;
	int numOfIters;
	int numOfGeneTypes;
	int k_nearest_neighbors;
	int numTopicsPerCellType;

	unordered_map<int,int> numOfGenes;

	string metaFile;
	string trainDataFile;
	string testDataFile;
	string trainPriorFile;
	string testPriorFile;
	string imputeTargetsFile; // for impute targets in new RSSamples
	string imputeRSSampleDataFile; // new RSSamples data for imputation

	unordered_map<pair<int,int>, GeneParams*>* geneParamsMap;

	vector<pair<int,int>> geneImputeTargets;

	string newDatafile;

	string trainedModelPrefix;

	string filePrefixDE;

	bool inferNewSampleRSMetagene;

	bool outputIntermediates;

	bool imputeNewSampleRSData; // impute new RSSample data

	int inferRSSampleParams_maxiter;

	bool inferTrainRSSampleMetagene_only; // infer train RSSample mix for imputation
	string trainRSSampleMetageneFile; // file to save train RSSample mix
	string trainRSSampleIdFile; // train RSSample ID to match with the trainRSSampleMetagene matrix rows

	// inference method
	string inference;

	// evaluation
	double testRSSamplesFrac;

	bool mar;

	int targetTypeId;
	bool evalTargetTypeOnly;

	vec logTrainLik;
	mat logTrainLik_breakdowns; // (model diagnostics) elbo breakdowns to check convergence issues

	vec logPredLik;

	vec trainTime;

	string outPrefix_trainData;
	string outPrefix_testData;

	string output_dir;

	int maxcores;

	double loglikdiff_thres;

	Decon(string datafile_train,
			string datafile_test,
			string prior_train,
			string prior_test,
			string metafile,
			string prefixDE,
			int topics,
			int numcelltypes,
			int maxiter, double diff_thres,
			string inferMethod,
			double testSetFrac,
			string newDataRSSampleh, string modelPrefix,
			bool inferNewSampleRSParams,

			bool inferTrainSampleRSMix,
			string trainRSSampleMixFile,
			string trainRSSampleIdFile,

			string imputeTargetsFileRSSampleh,
			string imputeRSSampleDataFileRSSampleh,
			int knn,
			bool saveIntermediates,
			bool missingAtRandom,
			int targetViewId,
			bool evalTargetViewOnly,
			int inferTestRSSampleThetaMaxiter,
			string output_RSSampleh,
			int maxthreads);

	void parseMetaInfo();
	JCVB0* initialize_infer();
	void parseTrainData(JCVB0* jcvb0); // change JCVB0 within function
	void parseTestData(JCVB0* jcvb0); // change JCVB0 within function
	void parseTrainPrior(JCVB0* jcvb0); // Train prior
	void parseTestPrior(JCVB0* jcvb0); // Test prior

	void parseImputeTargetsFile();
	void parseImputeRSSampleDataFile(JCVB0* jcvb0);

	// parse model files
	void parsePhi();

	// parse hyperparameter files
	rowvec parseAlpha();
	void parseBeta();

	//	void trainOnline();
	//
	void parsePhiFixed();

	void exportResults(JCVB0* jcvb0, int iter, bool verbose);

	void exportLogLik(int iter);

	void exportTrainTime(int iter);

	void exportLogLik_breakdowns(); // output all elbo breakdown elements per iteration

	void exportTestRSSampleData(JCVB0* jcvb0);

	JCVB0* parseTrainedModelFiles();

	void parseImputeTargetList(); // get a list of impute targets (i.e., a subset of variables in metainfo)

	JCVB0* parseNewData();

	void inferNewRSSampleMetagene(JCVB0* jcvb0, bool output_to_file=false);
	void inferNewRSSampleMetageneUnsupervised(JCVB0* jcvb0);

	void inferTrainRSSampleMetagene(); // infer and save the train RSSample mix for imputing test RSSamples
	void imputeNewGeneData(JCVB0* jcvb0, int nearestNeighborK);
	void imputeNewRSSampleData(int nearestNeighborK);
};

#endif
















