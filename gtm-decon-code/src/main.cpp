#include "Decon.h"
#include <math.h>

using namespace std;
using namespace arma;

Decon* parseCmdLine(int argc, char *argv[]) {

	string datafile_train = "";
	string datafile_test = "";
	string prior_train = "";
	string prior_test = "";
	string metafile = "";
	string prefixDE = "";

	string inferMethod="JCVB0";

	string newDatafile="";
	string trainedModelPrefix="";

	string trainRSSampleMetageneFile="";
	string trainRSSampleIdFile = "";

	string output_dir=".";

	// default parameters
	int iterations=5;
	double loglikdiff_thres = 1e-7;
	int numOfTopics=10;
	int numTopicsPerCellType=10;
	double testSetFrac=0.1;

	int targetViewId = 1;
	bool evalTargetViewOnly = false;

	bool inferNewSampleRSMetagene = false;

	bool inferTrainSampleRSMetagene = false;

	bool outputIntermediates = false;

	bool mar = false;

	int inferRSSampleParams_maxiter = 10;

	int maxcores = omp_get_max_threads(); // @suppress("Function cannot be resolved")

	string helpMsg = "./decon -f examples/toydata.txt -m1 examples/toymeta_gene.txt -i 10 -k 10";

	if (argc < 2) {

		cout << helpMsg << endl;

		exit(1);

	} else {

		for(int i = 1; i < argc; i++){

			string argComp = argv[i];

			if(argComp.compare("-h")== 0 || argComp.compare("--help")== 0) {

				cout << helpMsg << endl;

				exit(0);

			} else if(argComp.compare("-f")== 0 || argComp.compare("--trainDataFile")== 0) {

				datafile_train = string(argv[i+1]);

			} else if(argComp.compare("-t")== 0 || argComp.compare("--testDataFile")== 0) {

				datafile_test = string(argv[i+1]);

			} else if(argComp.compare("-trp")== 0 || argComp.compare("--trainPriorFile")== 0) {

				prior_train = string(argv[i+1]);

			} else if(argComp.compare("-tep")== 0 || argComp.compare("--testPriorFile")== 0) {

				prior_test = string(argv[i+1]);

			} else if(argComp.compare("-m")== 0 || argComp.compare("--metaFile")== 0) {

				metafile = string(argv[i+1]);

			} else if(argComp.compare("-de")== 0 || argComp.compare("--prefixDE")== 0) {

				prefixDE = string(argv[i+1]);

			} else if(argComp.compare("-i") == 0 || argComp.compare("--iter")== 0){

				iterations = atoi(argv[i+1]);

			} else if(argComp.compare("--convergenceThreshold") == 0){

				loglikdiff_thres = atof(argv[i+1]);

			} else if(argComp.compare("-k") == 0 || argComp.compare("--topics")== 0) {

				numOfTopics = atoi(argv[i+1]);

				if(numOfTopics < 1) throw::runtime_error("numOfTopics must be at least 1!");


			} else if(argComp.compare("-c") == 0 || argComp.compare("--numTopicsPerCellType")== 0) {

				numTopicsPerCellType = atoi(argv[i+1]);

			} else if(argComp.compare("-n") == 0 || argComp.compare("--inferenceMethod")== 0) {

				inferMethod = string(argv[i+1]);

			} else if(argComp.compare("--mar") == 0) {

				mar = true;

			} else if(argComp.compare("--targetTypeId")== 0) { // evaluating other imputation

				targetViewId = atoi(argv[i+1]);
				evalTargetViewOnly = true;

			} else if(argComp.compare("--newRSSamplesData")== 0) {

				newDatafile = string(argv[i+1]);

			} else if(argComp.compare("--trainedModelPrefix")== 0) {

				trainedModelPrefix =string(argv[i+1]);

			} else if(argComp.compare("--inferNewSampleRSMetagene")== 0) {

				inferNewSampleRSMetagene = true;

			} else if (argComp.compare("--inferTrainSampleRSMetagene")== 0) {

				inferTrainSampleRSMetagene = true;

			} else if(argComp.compare("--trainRSSampleMetageneFile")== 0) {

				trainRSSampleMetageneFile = string(argv[i+1]);

			} else if(argComp.compare("--trainRSSampleIdFile")== 0) {

				trainRSSampleIdFile = string(argv[i+1]);

			} else if(argComp.compare("--outputIntermediates")== 0) {

				outputIntermediates = true;

			} else if(argComp.compare("--inferRSSampleParams_maxiter")== 0) {

				inferRSSampleParams_maxiter = atoi(argv[i+1]);

			} else if(argComp.compare("-o")== 0 || argComp.compare("--output_dir")== 0) {

				output_dir = string(argv[i+1]);

			} else if(argComp.compare("--maxcores")== 0) {

				maxcores = atoi(argv[i+1]);

				if(maxcores > omp_get_max_threads()) // @suppress("Function cannot be resolved")
					maxcores = omp_get_max_threads(); // @suppress("Function cannot be resolved")
			}
		}
	}

	cout << "--------------------" << endl;
	cout << "Input arguments: " << endl;
	cout << "trainDataFile: " << datafile_train << endl;
	cout << "testDataFile: " << datafile_test << endl;
	cout << "trainPriorFile: " << prior_train << endl;
	cout << "testPriorFile: " << prior_test << endl;
	cout << "metaFile: " << metafile << endl;
	cout << "numTopics#: " << numOfTopics << endl;
	cout << "iter#: " << iterations << endl;
	cout << "convergenceThreshold: " << loglikdiff_thres << endl;
	cout << "inference method: " << inferMethod << endl;
	cout << "NMAR inference enabled: " << !mar << endl;
	cout << "maxcores: " << maxcores << endl;
	cout << "--------------------" << endl;


	return new Decon(datafile_train,
			datafile_test,
			prior_train,
			prior_test,
			metafile,
			prefixDE,
			numOfTopics,
			numTopicsPerCellType,
			iterations,
			loglikdiff_thres,
			inferMethod,
			testSetFrac,
			newDatafile,
			trainedModelPrefix,
			inferNewSampleRSMetagene,

			inferTrainSampleRSMetagene,
			trainRSSampleMetageneFile,
			trainRSSampleIdFile,

			outputIntermediates,
			mar, 
			targetViewId,
			evalTargetViewOnly,
			inferRSSampleParams_maxiter,
			output_dir,
			maxcores);
}

int main(int argc, char *argv[]) {

	arma_rng::set_seed(123);

	Decon* decon = parseCmdLine(argc, argv);

//	omp_set_num_threads(omp_get_max_threads()); // @suppress("Type cannot be resolved")

	omp_set_num_threads(decon->maxcores); // @suppress("Function cannot be resolved")

	// parse meta information of the clinical variables
	decon->parseMetaInfo();

	JCVB0 *infer = decon->initialize_infer();

	double loglikdiff = 0;
	double logprddiff = 0;

	double tStart = omp_get_wtime(); // @suppress("Function cannot be resolved")

	if(decon->inferNewSampleRSMetagene) {

		cout << "Use the trained model " << decon->trainedModelPrefix << endl;
		cout << "to infer new RSSample meta-genotypes from " << decon->newDatafile << endl;

		double tStart_parse = omp_get_wtime(); // @suppress("Function cannot be resolved")

		infer = decon->parseNewData();

		double tEnd_parse = omp_get_wtime(); // @suppress("Function cannot be resolved")

		printf("Data import time taken: %.2fs\n", (double) tEnd_parse - tStart_parse);

		double tStart_infer = omp_get_wtime(); // @suppress("Function cannot be resolved")

		decon->inferNewRSSampleMetagene(infer);
		//decon->inferNewRSSampleMetageneUnsupervised(infer);

		double tEnd_infer = omp_get_wtime(); // @suppress("Function cannot be resolved")

		printf("Meta-genotypes inference time taken: %.2fs\n", (double) tEnd_infer - tStart_infer);

	} else { // train

		cout << "decon->outputIntermediates " << decon->outputIntermediates << endl;

		int myints[] = {decon->numOfIters};

		std::vector<int> iter2print(myints, myints + sizeof(myints) / sizeof(int) );

		if(decon->outputIntermediates) {

			iter2print.clear();
			int myints2[] = {1,2,3,4,5,10,15,20,25,30,35,40,45,50,60,70,80,90,100,150,200,300,500,1000};
			int moreiters = decon->numOfIters/1e3;
			std::vector<int> iter2print2(myints2, myints2 + sizeof(myints2) / sizeof(int));

			for(std::vector<int>::iterator it = iter2print2.begin(); it != iter2print2.end(); it++) {
				if(*it < decon->numOfIters) {
					iter2print.push_back(*it);
				}
			}

			// output every 500 iterations from now on
			if(decon->numOfIters > 1000) {
				for(int i = 2; i <= moreiters; i++) {
					iter2print.push_back(i * 500);
				}
			}

			if(decon->numOfIters > (moreiters * 1e3)) iter2print.push_back(decon->numOfIters);
			cout << "Intermediate results will be saved from the following iterations: " << endl;
			for(std::vector<int>::iterator it = iter2print.begin(); it != iter2print.end(); it++)
				cout << *it << endl;
		}

		double tStart_parse = omp_get_wtime(); // @suppress("Function cannot be resolved")

		if(decon->inference.compare("JCVB0")==0) {

			decon->initialize_infer();

			cout << "Initialized" << endl;

			decon->parseTrainData(infer);

			cout << "Train data parsed" << endl;

			decon->parseTrainPrior(infer);

			cout << "Train prior parsed" << endl;

			if(decon->filePrefixDE != ""){
			      cout << "Using preset topics from " << decon->filePrefixDE << endl;
			      decon->parsePhiFixed();
			}

		} else {
			cout << decon->inference << endl;
			throw::runtime_error("invalid inference method");
		}
		if(decon->testDataFile.compare("")!= 0) {
			decon->parseTestData(infer);
		}

		double tEnd_parse = omp_get_wtime(); // @suppress("Function cannot be resolved")

		printf("Data import time taken: %.2fs\n", (double) tEnd_parse - tStart_parse);

		double trainStart = omp_get_wtime(); // @suppress("Function cannot be resolved")

		//cout << "main line 397" << endl;

		//		cout << endl << "before training" << endl << endl;
		//		infer->showParams();
		//		cout << endl;

		int iter = 0;

		//cout << "main line 405" << endl;

//		cout << "exportResults" << endl;
		if(decon->outputIntermediates) decon->exportResults(infer, iter, false);

//		cout << "trainLogLik" << endl;
		decon->logTrainLik(iter) = 0;

//		cout << "logPredLik" << endl;
		if(infer->testRSSamples->size() > 0) decon->logPredLik(iter) = infer->predictLogLik();

		printf("%d: logTrainLik: %.8f; logTrainLik diff: %.8f; logPredLik: %.8f; logPredLik diff: %.8f\n",
				iter+1,
				decon->logTrainLik(iter),loglikdiff,
				decon->logPredLik(iter), logprddiff);

		for (iter = 1; iter < decon->numOfIters; iter++) {

			infer->train(true);
			// infer->showParams();
			//cout << "main line 421" << endl;

			decon->logTrainLik(iter) = infer->trainLogLik();

			//cout << "main line 425" << endl;

			double trainSoFar = omp_get_wtime(); // @suppress("Function cannot be resolved")

			decon->trainTime(iter) = (double) (trainSoFar - trainStart);

			//cout << "main 431" << endl;

			//			decon->logTrainLik_breakdowns.row(iter) = infer->trainLogLik_breakdowns();
			//			cout << "logPredLik" << endl;

			if(infer->testRSSamples->size() > 0) decon->logPredLik(iter) = infer->predictLogLik();

			loglikdiff = decon->logTrainLik(iter) - decon->logTrainLik(iter-1);
			logprddiff = decon->logPredLik(iter) - decon->logPredLik(iter-1);

			printf("%d: logTrainLik: %.8f; logTrainLik diff: %.8f; logPredLik: %.8f; logPredLik diff: %.8f\n",
					iter+1,
					decon->logTrainLik(iter),loglikdiff,
					decon->logPredLik(iter), logprddiff);

			if(binary_search(iter2print.begin(), iter2print.end(), iter+1) && decon->outputIntermediates) {

				decon->exportResults(infer, iter, false);
				decon->exportLogLik(iter);
				decon->exportTrainTime(iter);
			}

			if(abs(loglikdiff) < decon->loglikdiff_thres && iter > 100 && !infer->svi) {
				break;
			}
		}

		cout << "Training completed after " << iter << " iterations" << endl;

		double trainEnd = omp_get_wtime(); // @suppress("Function cannot be resolved")

		printf("Training time taken: %.2fs\n", (double) (trainEnd - trainStart));

		if(iter < decon->numOfIters) { // converged

			decon->logTrainLik = decon->logTrainLik.head(iter);

			decon->logPredLik = decon->logPredLik.head(iter);

			decon->trainTime = decon->trainTime.head(iter);

			//			decon->logTrainLik_breakdowns = decon->logTrainLik_breakdowns.head_rows(iter+1);
		}

		if(infer->testRSSamples->size() > 0) {

			cout << "Inferring meta-genotypes of test RSSamples for the final round" << endl;

			infer->inferTestRSSampleParams_finalRun();
		}

		decon->exportResults(infer, iter, false);

		decon->exportLogLik(iter);

		decon->exportTrainTime(iter);

		//		decon->exportLogLik_breakdowns();
		//		decon->exportTestRSSampleData(infer);
	}

	double tEnd = omp_get_wtime(); // @suppress("Function cannot be resolved")

	printf("Total time taken: %.2fs\n", (double) (tEnd - tStart));

	return 0;
}




