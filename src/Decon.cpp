#include "Decon.h"

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include <sys/types.h>
#include <sys/stat.h>

using namespace std;
using namespace arma;

namespace bo = boost::iostreams;

Decon::Decon(string datafile_train,
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
		string trainRSSampleIDFile,

		string imputeTargetsFileRSSampleh,
		string imputeRSSampleDataFileRSSampleh,
		int impute_knn,
		bool saveIntermediates,
		bool missingAtRandom,
		int targetViewId,
		bool evalTargetViewOnly,
		int inferTestRSSampleThetaMaxiter,
		string output_RSSampleh,
		int maxthreads)
{
	numOfRSSamples = 0;

	numOfGeneTypes = 0;

	trainDataFile = datafile_train;

	testDataFile = datafile_test;

	trainPriorFile = prior_train;

	testPriorFile = prior_test;

	metaFile = metafile;

	filePrefixDE = prefixDE;

	geneParamsMap = new unordered_map<pair<int,int>, GeneParams*>();

	numOfTopics = topics;

	numTopicsPerCellType = numcelltypes;

	numOfIters = maxiter;

	loglikdiff_thres = diff_thres;

	inference = inferMethod;

	testRSSamplesFrac = testSetFrac;

	logPredLik = zeros<vec>(maxiter);

	logTrainLik = zeros<vec>(maxiter);

	trainTime = zeros<vec>(maxiter);

	logTrainLik_breakdowns = zeros<mat>(maxiter, 5);

	trainedModelPrefix = modelPrefix;

	newDatafile = newDataRSSampleh;

	inferNewSampleRSMetagene=inferNewSampleRSParams;


	// infer RSSample mix for imputation
	inferTrainRSSampleMetagene_only = inferTrainSampleRSMix;
	trainRSSampleMetageneFile = trainRSSampleMixFile;
	trainRSSampleIdFile = trainRSSampleIDFile;


	imputeNewSampleRSData = false;

	if(imputeTargetsFileRSSampleh.length() > 0) {

		imputeTargetsFile = imputeTargetsFileRSSampleh;

		if(imputeRSSampleDataFileRSSampleh.length() > 0) {
			imputeRSSampleDataFile = imputeRSSampleDataFileRSSampleh;
		} else {
			throw::runtime_error("imputeRSSampleDataFileRSSampleh is missing");
		}

		imputeNewSampleRSData = true;
	}

	inferRSSampleParams_maxiter = inferTestRSSampleThetaMaxiter;

	k_nearest_neighbors = impute_knn;

	outputIntermediates = saveIntermediates;

	mar = missingAtRandom;

	evalTargetTypeOnly = evalTargetViewOnly;

	targetTypeId = targetViewId;

	size_t lastindex = trainDataFile.find_last_of(".");

	outPrefix_trainData = trainDataFile.substr(0, lastindex);

	lastindex = testDataFile.find_last_of(".");

	outPrefix_testData = testDataFile.substr(0, lastindex);


	if(mar) {
		outPrefix_trainData = outPrefix_trainData + "_" + inference + "_mar" + "_K" + to_string(numOfTopics);
		outPrefix_testData = outPrefix_testData + "_" + inference + "_mar" + "_K" + to_string(numOfTopics);
	} else {
		outPrefix_trainData = outPrefix_trainData + "_" + inference + "_nmar" + "_K" + to_string(numOfTopics);
		outPrefix_testData = outPrefix_testData + "_" + inference + "_nmar" + "_K" + to_string(numOfTopics);
	}


	if(newDatafile.length() > 0) {

		ifstream ifile(newDatafile.c_str());

		string trainedModelFile = trainedModelPrefix + "_phi.csv";

		ifstream ifile2(trainedModelFile.c_str());

		if(ifile && ifile2) {

			ifile.close();

			ifile2.close();

		} else {
			cout << "newDatafile: " << newDatafile << endl;
			cout << "trainedModelFile: " << trainedModelFile << endl;
			throw::runtime_error("File does not exists!");
		}
	}

	output_dir = output_RSSampleh;

	numOfGenes = unordered_map<int,int>();

	maxcores = maxthreads;
}


JCVB0* Decon::initialize_infer() {

	JCVB0 *infer = new JCVB0();

	infer->initialize(
			numOfGenes, 
			numOfTopics, numTopicsPerCellType, numOfIters,
			*geneParamsMap 
			);

	infer->mar = mar;

	infer->inferTestRSSampleMetagene_maxiter_finalRun = inferRSSampleParams_maxiter;

	return infer;
}


void Decon::exportResults(JCVB0* jcvb0, int iter, bool verbose) {

	// model parameters
	// W x K x V
	string outfile_phi = outPrefix_trainData + "_iter" + to_string(iter) + "_phi.csv";
	string outfile_phi_normalized = outPrefix_trainData + "_iter" + to_string(iter) + "_phi_normalized.csv";

	// hyper-parameters
	string outfile_alpha = outPrefix_trainData + "_iter" + to_string(iter) + "_alpha.csv";
	string outfile_beta = outPrefix_trainData + "_iter" + to_string(iter) + "_beta.csv";

	if(verbose) {

		cout << "Save unnormalized model parameters:" << endl;
		cout << "Save phi to " << outfile_phi << endl;
		cout << "Save phi_normalized to " << outfile_phi_normalized << endl;

		cout << "Save hyper-parameters: " << endl;
		cout << "Save alpha to " << outfile_alpha << endl;
		cout << "Save beta to " << outfile_beta << endl;
	}


	// main param
	ofstream outfile_stream_phi;
	ofstream outfile_stream_phi_normalized;

	// hyper
	ofstream outfile_stream_alpha;
	ofstream outfile_stream_beta;

	// main param
	outfile_stream_phi.open(outfile_phi);
	outfile_stream_phi_normalized.open(outfile_phi_normalized);
	
	// hyper param
	outfile_stream_alpha.open(outfile_alpha);
	outfile_stream_beta.open(outfile_beta);

	// output matrix ordered by typeId, geneId, stateId from the unordered_map using an ordered map
	for(map<int, vector<int>>::iterator t = jcvb0->geneIds.begin(); t != jcvb0->geneIds.end(); t++) {

		int typeId = t->first;

		for(vector<int>::iterator w = jcvb0->geneIds[t->first].begin(); w != jcvb0->geneIds[t->first].end(); w++) {

			int geneId = *w;

			// caution: hashing missing keys will create an entry
			// export beta
			outfile_stream_beta << typeId << "," << geneId << "," <<
					jcvb0->geneParams[make_pair(typeId, geneId)]->beta << endl;

			outfile_stream_phi << typeId << "," << geneId;
			outfile_stream_phi_normalized << typeId << "," << geneId;


			for(int k=0; k<numOfTopics; k++) {

				outfile_stream_phi << "," << jcvb0->geneParams[make_pair(typeId, geneId)]->phi(k);
				outfile_stream_phi_normalized << "," << jcvb0->geneParams[make_pair(typeId, geneId)]->phi_normalized(k);
				
			}

			outfile_stream_phi << endl;
			outfile_stream_phi_normalized << endl;
		}
	}

	// export hyperparameter alpha
	outfile_stream_alpha << jcvb0->alpha(0);
	for(int k=1; k<numOfTopics; k++) {
		outfile_stream_alpha << "," << jcvb0->alpha(k);
	}
	outfile_stream_alpha << endl;

	// main params
	outfile_stream_phi.close();
	outfile_stream_phi_normalized.close();

	// hyper
	outfile_stream_alpha.close();

	//export metagene
	string outfile_metagene = outPrefix_trainData + "_iter" + to_string(iter) + "_metagene.csv";
	ofstream outfile_stream_metagene;
	outfile_stream_metagene.open(outfile_metagene);

	for(vector<SampleRS>::iterator RSSample = jcvb0->trainRSSamples->begin(); RSSample != jcvb0->trainRSSamples->end(); RSSample++) {
	        rowvec metagene = RSSample->Decompress(RSSample->metagene,numOfTopics);
		// output theta
		outfile_stream_metagene << metagene(0);

                for(int k=1; k<numOfTopics; k++) {
	                outfile_stream_metagene << "," << metagene(k);
                }
                outfile_stream_metagene << endl;
                metagene.clear();
        }
        outfile_stream_metagene.close();
}


void Decon::exportLogLik(int iter) {

	string outfile_logTrainLik = outPrefix_trainData + "_iter" + to_string(iter) + "_logTrainLik.txt";

	ofstream outfile_stream_logTrainLik;

	outfile_stream_logTrainLik.open(outfile_logTrainLik);

	vec loglik_train = logTrainLik.head(iter);

	for(int i=0; i<(int)loglik_train.n_elem; i++) {

		outfile_stream_logTrainLik << loglik_train(i) << endl;
	}

	outfile_stream_logTrainLik.close();

	// test predictive log lik
	if(testDataFile.compare("")!=0) {

		string outfile_logPredLik;

		outfile_logPredLik = outPrefix_testData + "_iter" + to_string(iter) + "_logPredLik.txt";

		ofstream outfile_stream_logPredLik;

		outfile_stream_logPredLik.open(outfile_logPredLik);

		vec loglik_pred = logPredLik.head(iter);

		for(int i=0; i<(int)loglik_pred.n_elem; i++) {

			outfile_stream_logPredLik << loglik_pred(i) << endl;
		}

		outfile_stream_logPredLik.close();
	}
}

void Decon::exportTrainTime(int iter) {

	string outfile_trainTime = outPrefix_trainData + "_iter" + to_string(iter) + "_trainTime.txt";

	ofstream outfile_stream_trainTime;

	outfile_stream_trainTime.open(outfile_trainTime);

	vec loglik_train = logTrainLik.head(iter);

	for(int i=0; i<(int)trainTime.n_elem; i++) {

		outfile_stream_trainTime << trainTime(i) << endl;
	}

	outfile_stream_trainTime.close();
}


void Decon::exportLogLik_breakdowns() {

	string outfile_logTrainLik = outPrefix_trainData + "_logTrainLikBKDW.csv";

	ofstream outfile_stream_logTrainLik;

	outfile_stream_logTrainLik.open(outfile_logTrainLik);

	for(int i=0; i<(int)logTrainLik_breakdowns.n_rows; i++) {

		outfile_stream_logTrainLik << i;

		for(int j=0; j<(int)logTrainLik_breakdowns.n_cols; j++) {

			outfile_stream_logTrainLik << "," << logTrainLik_breakdowns(i,j);
		}

		outfile_stream_logTrainLik << endl;
	}

	outfile_stream_logTrainLik.close();
}



// output test RSSample data for evaluation
void Decon::exportTestRSSampleData(JCVB0* jcvb0) {

	string outfile_testRSSample_obsGene = outPrefix_testData + "_testRSSample_obsGene.csv";
	string outfile_testRSSample_tarGene = outPrefix_testData + "_testRSSample_tarGene.csv";


	ofstream outfile_stream_testRSSample_obsGene;
	ofstream outfile_stream_testRSSample_tarGene;

	outfile_stream_testRSSample_obsGene.open(outfile_testRSSample_obsGene);
	outfile_stream_testRSSample_tarGene.open(outfile_testRSSample_tarGene);

	for(vector<SampleRS>::iterator RSSample = jcvb0->testRSSamples->begin(); RSSample != jcvb0->testRSSamples->end(); RSSample++) {

		for(unordered_map<pair<int,int>, int>::iterator iter = RSSample->geneDict.begin(); iter != RSSample->geneDict.end(); iter++) {

			// obs gene used for inferring theta
			if(RSSample->isTestGene[iter->first]) {

				// output format (same as input): RSSampleId,typeId,geneId,stateId,freq
				outfile_stream_testRSSample_tarGene <<
						RSSample->RSSampleId << "," <<
						iter->first.first << "," <<
						iter->first.second << "," <<
						1 << "," <<
						iter->second << endl;

			} else {


				outfile_stream_testRSSample_obsGene <<
						RSSample->RSSampleId << "," <<
						iter->first.first << "," <<
						iter->first.second << "," <<
						1 << "," <<
						iter->second << endl;
			}
		}
	}
}


void Decon::inferTrainRSSampleMetagene() {

	// parse model files
	JCVB0* jcvb0 = parseTrainedModelFiles();

	parseTrainData(jcvb0); // parse train RSSample data
	parseTrainPrior(jcvb0);

	// create directory to save the output files
	struct stat buffer;
	if (stat (output_dir.c_str(), &buffer) != 0) {
		const int status = mkdir(output_dir.c_str(), S_IRWXU);
		if (status == -1) {
			cout << "Error creating directory: " << output_dir;
			exit(1);
		}
	}

	// infer train RSSample mix
	cout << "Infer training RSSample mix" << endl;

	jcvb0->inferAllRSSampleParams(inferRSSampleParams_maxiter);

	cout << "Training mix inference completed" << endl;

	size_t lastindex = trainRSSampleMetageneFile.find_last_of(".");
	string prefix = trainRSSampleMetageneFile.substr(0, lastindex);
	string outfile_train_RSSampleid = output_dir + "/" + prefix + "_RSSampleId.csv";

	ofstream outfile_stream_train_RSSampleid;
	outfile_stream_train_RSSampleid.open(outfile_train_RSSampleid);

	mat train_RSSample_mix = zeros<mat>(jcvb0->trainRSSamples->size(), numOfTopics);

	int j = 0;
	for(vector<SampleRS>::iterator RSSample = jcvb0->trainRSSamples->begin(); RSSample != jcvb0->trainRSSamples->end();RSSample++) {
		rowvec metagene = RSSample->Decompress(RSSample->metagene,numOfTopics);
		train_RSSample_mix.row(j) = metagene;
		metagene.clear();
		outfile_stream_train_RSSampleid << RSSample->RSSampleId << endl;
		j++;
	}

	string outfile = output_dir + "/" + trainRSSampleMetageneFile;

	// save train RSSample mix
	train_RSSample_mix.save(outfile, csv_ascii);

	outfile_stream_train_RSSampleid.close();
}




// infer expectation of RSSamples' theta variables only
void Decon::inferNewRSSampleMetagene(JCVB0* jcvb0, bool output_to_file) {

	string inputFile;

	if (inferNewSampleRSMetagene)
		inputFile = newDatafile;
	else
		inputFile = trainDataFile;

	size_t lastindex0 = trainedModelPrefix.find_last_of("/");
	string trainedModel =  trainedModelPrefix.substr(lastindex0+1, trainedModelPrefix.length());
	size_t lastindex = inputFile.find_last_of(".");
	string prefix = inputFile.substr(0, lastindex);

	string outfile_testRSSample_theta = prefix + "_" + trainedModel + "_metagene.csv";

	cout << "Saving inferred RSSample metagene: " << outfile_testRSSample_theta << endl;

	ofstream outfile_stream_testRSSample_theta;

	outfile_stream_testRSSample_theta.open(outfile_testRSSample_theta);

	int j = 0;

	vector<SampleRS>::iterator RSSample0 = jcvb0->testRSSamples->begin();

	int D = jcvb0->testRSSamples->size();

#pragma omp parallel for shared(j)
	for(j=0; j<D; j++) {

		vector<SampleRS>::iterator RSSamplej = RSSample0 + j;

		//cout << j << endl;
		//jcvb0->inferRSSampleParams(RSSamplej, inferRSSampleParams_maxiter);
		jcvb0->inferRSSampleParamsUnsupervised(RSSamplej, inferRSSampleParams_maxiter);

		// free up the memory allocated for the RSSample gamma hash
		RSSamplej->gamma.clear();

	}

	/*
	for(vector<SampleRS>::iterator RSSample = jcvb0->testRSSamples->begin(); RSSample != jcvb0->testRSSamples->end(); RSSample++) {
		//rowvec metagene = RSSample->Decompress(RSSample->metagene,numOfTopics);

		// output theta
		cout << RSSample->metagene_normalized(0);

		for(int k=1; k<numOfTopics; k++) {
			cout << "," << RSSample->metagene_normalized(k);
		}
		cout << endl;
	}*/

	
	for(vector<SampleRS>::iterator RSSample = jcvb0->testRSSamples->begin(); RSSample != jcvb0->testRSSamples->end(); RSSample++) {
		//rowvec metagene = RSSample->Decompress(RSSample->metagene,numOfTopics);

		// output theta
		outfile_stream_testRSSample_theta << RSSample->metagene_normalized(0);

		for(int k=1; k<numOfTopics; k++) {

			outfile_stream_testRSSample_theta << "," << RSSample->metagene_normalized(k);
		}
		outfile_stream_testRSSample_theta << endl;
	}

	outfile_stream_testRSSample_theta.close();
}

void Decon::imputeNewGeneData(JCVB0* jcvb0, int nearestNeighborK) {

	int t = 0;
	umat target_gene_true = zeros<umat>(jcvb0->imputeTargetRSSamples->size(), geneImputeTargets.size());
	mat target_gene_pred = zeros<mat>(jcvb0->imputeTargetRSSamples->size(), geneImputeTargets.size());

	for(vector<pair<int,int>>::iterator tar = geneImputeTargets.begin(); tar != geneImputeTargets.end(); tar++) {

		int tar_typeId = tar->first;
		int tar_geneId = tar->second;

		cout << "Predict target code: " << tar_typeId << ", " << tar_geneId << endl;

		int j = 0;
		vector<SampleRS>::iterator RSSample0 = jcvb0->imputeTargetRSSamples->begin();

		mat test_RSSample_mix = zeros<mat>(jcvb0->imputeTargetRSSamples->size(), numOfTopics);

#pragma omp parallel for shared(j)
		for(j=0; j < (int) jcvb0->imputeTargetRSSamples->size(); j++) {

			pair<int,int> targeneid = make_pair(tar_typeId, tar_geneId);

			vector<SampleRS>::iterator tar_RSSample = RSSample0 + j;

			rowvec metagene = tar_RSSample->Decompress(tar_RSSample->metagene,numOfTopics);
			rowvec metagene_normalized = tar_RSSample->Decompress(tar_RSSample->metagene_normalized,numOfTopics);

			test_RSSample_mix.row(j) = metagene;

			unordered_map<pair<int,int>, int>::const_iterator hasit = tar_RSSample->geneDict.find(targeneid);

			double targene_freq = 0;

			if(hasit != tar_RSSample->geneDict.end()) { // actual positive
				targene_freq = tar_RSSample->geneDict[targeneid];
				tar_RSSample->geneDict.erase(hasit);
				target_gene_true(j,t) = 1;
			} // else actual positive

			// infer target RSSamples' mix without the target label
			jcvb0->inferRSSampleParams(tar_RSSample, inferRSSampleParams_maxiter);

			// save test RSSample mix for export
			test_RSSample_mix.row(j) = metagene_normalized;

			vec dist = zeros<vec>(jcvb0->trainRSSamples->size());

			// find most similar training RSSamples based on their mixes
			int i = 0;
			for(vector<SampleRS>::iterator RSSample = jcvb0->trainRSSamples->begin(); RSSample != jcvb0->trainRSSamples->end(); RSSample++) {
				dist(i) = accu(square(metagene_normalized - metagene_normalized)); // 1 x K
				i++;
			}

			// pick the top K most similar training RSSamples
			uvec top_trainRSSample_indices = sort_index(dist);
			top_trainRSSample_indices = top_trainRSSample_indices.head(nearestNeighborK);

			uvec bot_trainRSSample_indices = sort_index(dist, "descend");
			bot_trainRSSample_indices = bot_trainRSSample_indices.head(nearestNeighborK);


			// DEBUG BEGINS
//			cout << "test RSSample: " << j << endl;
//			cout << "dist(top_trainRSSample_indices): " << dist(top_trainRSSample_indices) << endl;
//			cout << "dist(bot_trainRSSample_indices): " << dist(bot_trainRSSample_indices) << endl;
			// DEBUG ENDS


			// predict target code based on the average of the top K most similar training RSSamples' target label
			double RSSampleSum = 0;
			double totalSum = 0;
			for(int i=0; i<nearestNeighborK; i++) {
				int top_RSSample_index = top_trainRSSample_indices(i);
				vector<SampleRS>::iterator top_trainRSSample = jcvb0->trainRSSamples->begin() + top_RSSample_index;
				if(top_trainRSSample->geneDict.find(targeneid) != top_trainRSSample->geneDict.end()) {
					double weighted_count = top_trainRSSample->geneDict[targeneid];
					RSSampleSum += weighted_count;
					totalSum += weighted_count;
				} else {
					totalSum++;
				}

			}
			target_gene_pred(j,t) = RSSampleSum/totalSum;

			// restore the erased observation
			if(target_gene_true(j,t)==1) { // actual positive
				tar_RSSample->geneDict[targeneid] = targene_freq;
			}

			metagene.clear();
			metagene_normalized.clear();
		}

		// save test RSSample mix matrix in csv format
		//		string outfile_tar_RSSample_mix = output_dir + "/target_RSSample_mix_" + to_string(tar_typeId) + "_" + to_string(tar_geneId) + ".csv";
		//		test_RSSample_mix.save(outfile_tar_RSSample_mix, csv_ascii);

		t++;
	}

	// save impute target id
	string outfile_target_geneid = output_dir + "/target_geneid.csv";
	ofstream outfile_stream_target_geneid;
	outfile_stream_target_geneid.open(outfile_target_geneid);
	for(vector<pair<int,int>>::iterator tar = geneImputeTargets.begin(); tar != geneImputeTargets.end(); tar++) {
		int tar_typeId = tar->first;
		int tar_geneId = tar->second;
		outfile_stream_target_geneid << tar_typeId << "," << tar_geneId << endl;
	}
	outfile_stream_target_geneid.close();

	// save prediction and true labels in matrix csv format
	string outfile_target_pred = output_dir + "/target_gene_pred.csv";
	string outfile_target_true = output_dir + "/target_gene_true.csv";

	target_gene_pred.save(outfile_target_pred, csv_ascii);
	target_gene_true.save(outfile_target_true, csv_ascii);
}


void Decon::imputeNewRSSampleData(int nearestNeighborK) {

	// parse model files
	JCVB0* jcvb0 = parseTrainedModelFiles();

	// infer train RSSamples' mix
	parseImputeTargetsFile();

	parseTrainData(jcvb0); // parse train RSSample data

	parseImputeRSSampleDataFile(jcvb0);  // parse new RSSample data for imputation

	// create directory to save the output files
	struct stat buffer;
	if (stat (output_dir.c_str(), &buffer) != 0) {
		const int status = mkdir(output_dir.c_str(), S_IRWXU);
		if (status == -1) {
			cout << "Error creating directory: " << output_dir;
			exit(1);
		}
	}

	// infer train RSSample mix
//	cout << "Infer training RSSample mix" << endl;
//	jcvb0->inferAllRSSampleParams(inferRSSampleParams_maxiter);
//	cout << "Training mix inference completed" << endl;

	// insert previously inferred topic mixture to each train RSSample
	mat train_RSSample_mix;
	train_RSSample_mix.load(trainRSSampleMetageneFile, csv_ascii);
	vec train_RSSample_id;
	train_RSSample_id.load(trainRSSampleIdFile, csv_ascii);

	for(int j=0; j<(int)train_RSSample_id.n_rows; j++) {

		vector<SampleRS>::iterator RSSample = jcvb0->trainRSSamples->begin() + j;

		if(RSSample->RSSampleId == train_RSSample_id(j)) {

			for (unordered_map<int, double>::iterator iter = RSSample->topicMap.begin(); iter != RSSample->topicMap.end(); iter++){
				int i = std::distance(RSSample->topicMap.begin(), iter);
				RSSample->metagene[i] = train_RSSample_mix(j,iter->first);
			}
			RSSample->metagene_normalized = RSSample->metagene/accu(RSSample->metagene);

		} else {

			cout << "RSSample->RSSampleId" << RSSample->RSSampleId;
			cout << "train_RSSample_id(j)" << train_RSSample_id(j);
			throw::runtime_error("Train RSSamples are in same order train RSSample mix");
		}
	}

	imputeNewGeneData(jcvb0, nearestNeighborK);

	// save target RSSample id
	string outfile_target_RSSampleid = output_dir + "/target_RSSampleid.csv";
	ofstream outfile_stream_target_RSSampleid;
	outfile_stream_target_RSSampleid.open(outfile_target_RSSampleid);
	for(vector<SampleRS>::iterator RSSampleiter = jcvb0->imputeTargetRSSamples->begin();
			RSSampleiter != jcvb0->imputeTargetRSSamples->end(); RSSampleiter++) {
		outfile_stream_target_RSSampleid << RSSampleiter->RSSampleId << endl;
	}
	outfile_stream_target_RSSampleid.close();
}


// infer expectation of RSSamples' theta variables only
void Decon::inferNewRSSampleMetageneUnsupervised(JCVB0* jcvb0) {

	string inputFile;

	if (inferNewSampleRSMetagene)
		inputFile = newDatafile;
	else
		inputFile = trainDataFile;

	int j = 0;

	vector<SampleRS>::iterator RSSample0 = jcvb0->testRSSamples->begin();

	int D = jcvb0->testRSSamples->size();


#pragma omp parallel for shared(j)
	for(j=0; j<D; j++) {

		vector<SampleRS>::iterator RSSamplej = RSSample0 + j;

		//jcvb0->inferRSSampleParams(RSSamplej, inferRSSampleParams_maxiter);
		jcvb0->inferRSSampleParamsUnsupervised(RSSamplej, inferRSSampleParams_maxiter);


		//free up the memory allocated for the RSSample gamma hash
		RSSamplej->gamma.clear();

	}


	size_t lastindex0 = trainedModelPrefix.find_last_of("/");
	string trainedModel =  trainedModelPrefix.substr(lastindex0+1, trainedModelPrefix.length());
	size_t lastindex = inputFile.find_last_of(".");
	string prefix = inputFile.substr(0, lastindex);

	string outfile_testRSSample_theta = prefix + "_" + trainedModel + "_metagene.csv";

	cout << "Saving inferred RSSample metagene: " << outfile_testRSSample_theta << endl;

	ofstream outfile_stream_testRSSample_theta;

	outfile_stream_testRSSample_theta.open(outfile_testRSSample_theta);
	
	for(vector<SampleRS>::iterator RSSample = jcvb0->testRSSamples->begin(); RSSample != jcvb0->testRSSamples->end(); RSSample++) {

		// output theta
		outfile_stream_testRSSample_theta << RSSample->metagene_normalized(0);

		for(int k=1; k<numOfTopics; k++) {

			outfile_stream_testRSSample_theta << "," << RSSample->metagene_normalized(k);
		}
		outfile_stream_testRSSample_theta << endl;
	}

	outfile_stream_testRSSample_theta.close();
}
