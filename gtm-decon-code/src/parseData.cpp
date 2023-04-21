#include "Decon.h"
#include "SampleRS.h"
#include "JCVB0.h"
#include "GeneParams.h"

using namespace std;
using namespace arma;


void Decon::parseMetaInfo() {

	// parse meta data file
	ifstream metafile(metaFile.c_str());

	int typeId,geneId,stateNum;

	numOfGeneTypes = 0;

	int oldTypeId = -1;

	while(metafile >> typeId >> geneId >> stateNum) {

		if(oldTypeId!=typeId) { // new data type found

			numOfGeneTypes++;

			oldTypeId=typeId;
		}

		if(oldTypeId!=typeId) { // new data type found
			numOfGeneTypes++;
		}

		GeneParams *genePar = new GeneParams();

		genePar->phi = randu<rowvec>(numOfTopics);

		(*geneParamsMap)[make_pair(typeId, geneId)] = genePar;

		if(numOfGenes.find(typeId)==numOfGenes.end()) {

			numOfGenes[typeId] = 1;

		} else {

			numOfGenes[typeId]++;
		}

		oldTypeId=typeId;
	}

	metafile.close();
}


void Decon::parseTrainData(JCVB0* jcvb0) {

	// parse RSSample data file
	ifstream datafile(trainDataFile.c_str());

	int typeId,geneId;

	int RSSampleId,stateId,freq;

	bool eof = false;

	if(!(datafile >> RSSampleId >> typeId >> geneId >> stateId >> freq)) {
		eof = true;
	}

	while(!eof) {

		unordered_map<pair<int, int>, int>* geneMap = new unordered_map<pair<int, int>, int>();

		int oldRSSampleId = RSSampleId;

		while(RSSampleId == oldRSSampleId) {

			(*geneMap)[make_pair(typeId,geneId)] = freq;


			if(!(datafile >> RSSampleId >> typeId >> geneId >> stateId >> freq)) {

				eof = true;
				break;
			}

			//printf("RSSampleId: %d; typeId: %d; geneId: %d; stateId: %d; freq: %d\n", RSSampleId,typeId,geneId,stateId,freq);
		}

		SampleRS* newRSSample = new SampleRS(oldRSSampleId,
				*geneMap,
				numOfTopics);

		jcvb0->trainRSSamples->push_back(*newRSSample);

		jcvb0->C_train += newRSSample->Cj_train;

		// printf("%.3f\n", (double)(jcvb0->trainRSSamples->size() + jcvb0->testRSSamples->size())/numOfRSSamples);
	}

	datafile.close();

	jcvb0->D_train = jcvb0->trainRSSamples->size();

	cout << "numOfDataTypes: " << numOfGeneTypes << endl;
	cout << "numOfGenes: " << jcvb0->geneParams.size() << endl;
	cout << "numOfRSSamples: " << jcvb0->D_train << endl;
	cout << "C_train: " << jcvb0->C_train << endl;
	printf("Training data file parsing completed.\n");
	cout << "--------------------" << endl;

//	throw::runtime_error("parseData");
}


void Decon::parseTestData(JCVB0* jcvb0) {

	// parse RSSample data file
	ifstream datafile(testDataFile.c_str());

	int RSSampleId,typeId,geneId,stateId,freq,obs;

	bool eof = false;


	if(!(datafile >> RSSampleId >> typeId >> geneId >> stateId >> freq >> obs)) {
		eof = true;
	}

//	printf("RSSampleId: %d; typeId: %d; geneId: %d; stateId: %d; freq: %d; obs: %d\n", RSSampleId,typeId,geneId,stateId,freq,obs);

	while(!eof) {

		unordered_map<pair<int, int>, int>* geneMap = new unordered_map<pair<int, int>, int>();

		int oldRSSampleId = RSSampleId;

		while(RSSampleId == oldRSSampleId) {

			(*geneMap)[make_pair(typeId,geneId)] = freq;

			if(!(datafile >> RSSampleId >> typeId >> geneId >> stateId >> freq >> obs)) {

				eof = true;
				break;
			}

//			printf("RSSampleId: %d; typeId: %d; geneId: %d; stateId: %d; freq: %d; obs: %d\n", RSSampleId,typeId,geneId,stateId,freq,obs);
		}

		SampleRS* newRSSample = new SampleRS(oldRSSampleId,
				*geneMap,
				numOfTopics);


		newRSSample->isTestRSSample = true;

		if(evalTargetTypeOnly) {
			newRSSample->assignTargetView(targetTypeId);
		} 
		jcvb0->testRSSamples->push_back(*newRSSample);
		jcvb0->C_tar += newRSSample->Cj_tar;
	}

	datafile.close();

	jcvb0->D_test = jcvb0->testRSSamples->size();

	printf("testRSSamples: %d\n", (int)jcvb0->D_test);
	cout << "C_tar: " << jcvb0->C_tar << endl;
	printf("Testing data file parsing completed.\n");
	cout << "--------------------" << endl;
}


void Decon::parseTrainPrior(JCVB0* jcvb0) {

	// parse RSSample data file
	ifstream datafile(trainPriorFile.c_str());

	int RSSampleId,topId;
	double prob;

	if(!(datafile >> RSSampleId >> topId >> prob)) {
		throw runtime_error("Prior file empty");
	}

	//printf("RSSampleId: %d; typeId: %d; geneId: %d; stateId: %d; freq: %d\n", RSSampleId,typeId,geneId,stateId,freq);
	//printf("RSSampleId: %d; topId: %d; prob: %.3f\n",RSSampleId,topId,prob);

	for(vector<SampleRS>::iterator RSSample = jcvb0->trainRSSamples->begin(); RSSample != jcvb0->trainRSSamples->end(); RSSample++) {

		if (RSSample->RSSampleId != RSSampleId){
			throw runtime_error("ID mismatch: train ID = " + to_string(RSSample->RSSampleId) + ", prior ID = " + to_string(RSSampleId));
		}

		//cout << "RSSample:" << RSSample->RSSampleId << "RSSampleId:" << RSSampleId << "prob:" << prob << endl;
		//RSSample->Ki = 0;
		while(RSSample->RSSampleId == RSSampleId) {

			(RSSample->topicMap)[topId] = prob;
			//printf("RSSampleId: %d; topId: %d; prob: %.3f; topicmap: %.3f; Ki: %d\n",RSSampleId,topId,prob,RSSample->topicMap[topId],RSSample->Ki);
			RSSample->Ki += 1;

			if(!(datafile >> RSSampleId >> topId >> prob)) {
				break;
			}

//			printf("RSSampleId: %d; typeId: %d; geneId: %d; stateId: %d; freq: %d\n", RSSampleId,typeId,geneId,stateId,freq);
		}

		RSSample->metagene = randu<rowvec>(RSSample->Ki);
		RSSample->metagene_normalized = RSSample->metagene/accu(RSSample->metagene);

		// printf("%.3f\n", (double)(jcvb0->trainRSSamples->size() + jcvb0->testRSSamples->size())/numOfRSSamples);
	}

	datafile.close();

	printf("Prior data file parsing completed.\n");
	cout << "--------------------" << endl;

//	throw::runtime_error("parseData");
}


void Decon::parsePhi() {

	// parse trained state eta parameters
	string trainedModelState = trainedModelPrefix + "_phi.csv";

	cout << "--------------------" << endl;
	cout << "Parsing phi: " << trainedModelState << " ... ";

	// parse trained phi file
	ifstream modelfile(trainedModelState.c_str());

	if(!modelfile.good()) {
		cout << trainedModelState << endl;
		throw runtime_error("file does not exist");
	}

	string tmp;

	getline(modelfile, tmp, ',');
	int typeId = atoi(tmp.c_str());
	getline(modelfile, tmp, ',');
	int geneId = atoi(tmp.c_str());

	while(modelfile.good()) {

		for(int k=0; k<numOfTopics-1; k++) {

			getline(modelfile, tmp, ',');

			(*geneParamsMap)[make_pair(typeId,geneId)]->phi(k) = atof(tmp.c_str());

//			if(typeId == 1 && geneId == 30774) {
//				cout << "typeId: " << typeId << "; " << "geneId: " << geneId << "; " << "k: " <<
//						k << "; " << (*geneParamsMap)[make_pair(typeId,geneId)]->phi(k) << endl;
//			}
		}

//		if(typeId == 1 && geneId == 30774)
//			throw runtime_error("stop here");

		getline(modelfile, tmp, '\n');

		(*geneParamsMap)[make_pair(typeId,geneId)]->phi(numOfTopics-1) = atof(tmp.c_str());

//		cout << "typeId: " << typeId << "; " << "geneId: " << geneId << "; " << "k: " <<
//				numOfTopics-1 << "; " << (*geneParamsMap)[make_pair(typeId,geneId)]->phi(numOfTopics-1) << endl;

		getline(modelfile, tmp, ',');
		typeId = atoi(tmp.c_str());
		getline(modelfile, tmp, ',');
		geneId = atoi(tmp.c_str());
	}

	modelfile.close();

	cout << "done." << endl;
}

void Decon::parsePhiNormalized() {

	// parse trained state eta parameters
	string trainedModelState = trainedModelPrefix + "_phi_normalized.csv";

	cout << "--------------------" << endl;
	cout << "Parsing phi normalized: " << trainedModelState << " ... ";

	// parse trained phi file
	ifstream modelfile(trainedModelState.c_str());

	if(!modelfile.good()) {
		cout << trainedModelState << endl;
		throw runtime_error("file does not exist");
	}

	string tmp;

	getline(modelfile, tmp, ',');
	int typeId = atoi(tmp.c_str());
	getline(modelfile, tmp, ',');
	int geneId = atoi(tmp.c_str());

	while(modelfile.good()) {

		for(int k=0; k<numOfTopics-1; k++) {

			getline(modelfile, tmp, ',');

			(*geneParamsMap)[make_pair(typeId,geneId)]->phi(k) = atof(tmp.c_str());

//			if(typeId == 1 && geneId == 30774) {
//				cout << "typeId: " << typeId << "; " << "geneId: " << geneId << "; " << "k: " <<
//						k << "; " << (*geneParamsMap)[make_pair(typeId,geneId)]->phi(k) << endl;
//			}
		}

//		if(typeId == 1 && geneId == 30774)
//			throw runtime_error("stop here");

		getline(modelfile, tmp, '\n');

		(*geneParamsMap)[make_pair(typeId,geneId)]->phi(numOfTopics-1) = atof(tmp.c_str());

//		cout << "typeId: " << typeId << "; " << "geneId: " << geneId << "; " << "k: " <<
//				numOfTopics-1 << "; " << (*geneParamsMap)[make_pair(typeId,geneId)]->phi(numOfTopics-1) << endl;

		getline(modelfile, tmp, ',');
		typeId = atoi(tmp.c_str());
		getline(modelfile, tmp, ',');
		geneId = atoi(tmp.c_str());
	}

	modelfile.close();

	cout << "done." << endl;
}

rowvec Decon::parseAlpha() {

	string tmp;

	// parse trained alpha parameter file
	string trainedModelAlpha = trainedModelPrefix + "_alpha.csv";

	cout << "--------------------" << endl;
	cout << "Parsing alpha: " << trainedModelAlpha << " ... ";

	ifstream modelfile(trainedModelAlpha.c_str());

	if(!modelfile.good()) {
		cout << trainedModelAlpha << endl;
		throw runtime_error("file does not exist");
	}

	rowvec alpha = zeros<rowvec>(numOfTopics);

	for(int k=0; k<numOfTopics-1; k++) {

		getline(modelfile, tmp, ',');

		alpha(k) = atof(tmp.c_str());
	}

	getline(modelfile, tmp, '\n');

	alpha(numOfTopics-1) = atof(tmp.c_str());

	modelfile.close();

	cout << "done." << endl;

	return alpha;
}


void Decon::parseBeta() {

	// parse trained alpha parameter file
	string trainedModelBeta= trainedModelPrefix + "_beta.csv";

	cout << "--------------------" << endl;
	cout << "Parsing beta: " << trainedModelBeta << " ... ";

	ifstream modelfile(trainedModelBeta.c_str());

	if(!modelfile.good()) {
		cout << trainedModelBeta << endl;
		throw runtime_error("file does not exist");
	}

	string tmp;

	getline(modelfile, tmp, ',');
	int typeId = atoi(tmp.c_str());
	getline(modelfile, tmp, ',');
	int geneId = atoi(tmp.c_str());

	while(modelfile.good()) {

		getline(modelfile, tmp, '\n');

		(*geneParamsMap)[make_pair(typeId,geneId)]->beta = atof(tmp.c_str());

		getline(modelfile, tmp, ',');
		typeId = atoi(tmp.c_str());
		getline(modelfile, tmp, ',');
		geneId = atoi(tmp.c_str());
	}

	modelfile.close();

	cout << "done." << endl;
}


JCVB0* Decon::parseTrainedModelFiles() {

	if(trainedModelPrefix.length() == 0) {
		throw runtime_error("trainedModelPrefix is undefined");
	}

	//parsePhi();
	parsePhiNormalized();

	parseBeta();

	JCVB0* jcvb0 = new JCVB0();

	jcvb0->initialize(
			numOfGenes, 
			numOfTopics, numTopicsPerCellType, numOfIters,
			*geneParamsMap);

	jcvb0->alpha = parseAlpha();

	jcvb0->updateParamSums();

	jcvb0->normalizeParams();

	jcvb0->mar = mar;

	return jcvb0;
}


JCVB0* Decon::parseNewData() {

	JCVB0* jcvb0 = parseTrainedModelFiles();

	cout << "--------------------" << endl;
	printf("Trained model data files parsing completed.\n");
	cout << "--------------------" << endl;

//	jcvb0->showParams();

	// parse RSSample data file
	ifstream datafile(newDatafile.c_str());

	if(!datafile.is_open()) {
		cout << newDatafile << endl;
		throw::runtime_error("Cannot open file");
	}

	int RSSampleId,typeId,geneId,stateId,freq;

	bool eof = false;

	if(!(datafile >> RSSampleId >> typeId >> geneId >> stateId >> freq)) {
		eof = true;
	}

//	printf("RSSampleId: %d; typeId: %d; geneId: %d; stateId: %d; freq: %d\n", RSSampleId,typeId,geneId,stateId,freq);

	while(!eof) {

		unordered_map<pair<int, int>, int>* geneMap = new unordered_map<pair<int, int>, int>();

		int oldRSSampleId = RSSampleId;

		while(RSSampleId == oldRSSampleId) {

			(*geneMap)[make_pair(typeId,geneId)] = freq;

			if(!(datafile >> RSSampleId >> typeId >> geneId >> stateId >> freq)) {

				eof = true;
				break;
			}
		}

		SampleRS* newRSSample = new SampleRS(oldRSSampleId,
				*geneMap,
				numOfTopics);

		newRSSample->isTestRSSample = true;

		jcvb0->testRSSamples->push_back(*newRSSample);

		jcvb0->C_train += newRSSample->Cj_train;
	}

	datafile.close();

	cout << "--------------------" << endl;
	printf("New RSSample data files parsing completed.\n");
	cout << "--------------------" << endl;

	cout << "numOfGeneParamTypes: " << numOfGeneTypes << endl;
	cout << "numOfGenes: " << jcvb0->geneParams.size() << endl;

	printf("testRSSamples: %d\n", (int)jcvb0->testRSSamples->size());
	cout << "C_train: " << jcvb0->C_train<< endl;

	cout << "--------------------" << endl;

	return jcvb0;
}


void Decon::parsePhiFixed() {

	// parse trained state phi parameters
        string trainedModelState = filePrefixDE + "_phi.csv";
        string trainedModelStateNorm = filePrefixDE + "_phi_normalized.csv";

	cout << "--------------------" << endl;
	cout << "Parsing phi for celltypes from single cell topics: " << trainedModelState << " ... ";

	// parse trained phi file
	ifstream modelfile(trainedModelState.c_str());

	if(!modelfile.good()) {
		cout << trainedModelState << endl;
		throw runtime_error("file does not exist");
	}

	string tmp;

	getline(modelfile, tmp, ',');
	int typeId = atoi(tmp.c_str());
	getline(modelfile, tmp, ',');
	int pheId = atoi(tmp.c_str());

	while(modelfile.good()) {

		for(int k=0; k<numOfTopics-1; k++) {

			getline(modelfile, tmp, ',');

			(*geneParamsMap)[make_pair(typeId,pheId)]->phi(k) = atof(tmp.c_str());

//			if(typeId == 1 && pheId == 30774) {
//				cout << "typeId: " << typeId << "; " << "pheId: " << pheId << "; " << "k: " <<
//						k << "; " << (*geneParamsMap)[make_pair(typeId,pheId)]->phi(k) << endl;
//			}
		}

//		if(typeId == 1 && pheId == 30774)
//			throw runtime_error("stop here");

		getline(modelfile, tmp, '\n');

		(*geneParamsMap)[make_pair(typeId,pheId)]->phi(numOfTopics-1) = atof(tmp.c_str());

//		cout << "typeId: " << typeId << "; " << "pheId: " << pheId << "; " << "k: " <<
//				numOfTopics-1 << "; " << (*geneParamsMap)[make_pair(typeId,pheId)]->phi(numOfTopics-1) << endl;

		getline(modelfile, tmp, ',');
		typeId = atoi(tmp.c_str());
		getline(modelfile, tmp, ',');
		pheId = atoi(tmp.c_str());
	}

	modelfile.close();

	cout << "done." << endl;

	// parse trained phi normalized file
	ifstream modelnormfile(trainedModelStateNorm.c_str());
	if(!modelnormfile.good()) {
		cout << trainedModelStateNorm << endl;
		throw runtime_error("file does not exist");
	}
	cout << "Parsing phiNorm for celltypes from single cell topics: " << trainedModelStateNorm << " ... ";

	getline(modelnormfile, tmp, ',');
	typeId = atoi(tmp.c_str());
	getline(modelnormfile, tmp, ',');
	pheId = atoi(tmp.c_str());

	while(modelnormfile.good()) {

		for(int k=0; k<numOfTopics-1; k++) {

			getline(modelnormfile, tmp, ',');

			(*geneParamsMap)[make_pair(typeId,pheId)]->phi_normalized(k) = atof(tmp.c_str());

//			if(typeId == 1 && pheId == 30774) {
//				cout << "typeId: " << typeId << "; " << "pheId: " << pheId << "; " << "k: " <<
//						k << "; " << (*geneParamsMap)[make_pair(typeId,pheId)]->phi(k) << endl;
//			}
		}

//		if(typeId == 1 && pheId == 30774)
//			throw runtime_error("stop here");

		getline(modelnormfile, tmp, '\n');

		(*geneParamsMap)[make_pair(typeId,pheId)]->phi_normalized(numOfTopics-1) = atof(tmp.c_str());

//		cout << "typeId: " << typeId << "; " << "pheId: " << pheId << "; " << "k: " <<
//				numOfTopics-1 << "; " << (*geneParamsMap)[make_pair(typeId,pheId)]->phi(numOfTopics-1) << endl;

		getline(modelnormfile, tmp, ',');
		typeId = atoi(tmp.c_str());
		getline(modelnormfile, tmp, ',');
		pheId = atoi(tmp.c_str());
	}
 	cout << "done." << endl;
}
