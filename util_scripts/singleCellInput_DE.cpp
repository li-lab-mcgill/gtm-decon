//This script generates input for single cell or bulkRS datasets for DE genes - where it requires both number of celltypes, and number of topics to be trained per category
//Usage: 
//g++ SingleCellInput.cpp -o SingleCellInput -larmadillo
//./SingleCellInput ../../../scDECON_CANCER/data/MP/MP_counts_matrix.tab ../../../scDECON_CANCER/data/MP/ETM/cell_labels_oh.csv ../scMP/scMP_trainData.txt ../scMP/scMP_priorData.txt ../scMP/scMP_meta.txt ../scMP/scMP_genes.txt
//./SingleCellInput ../../../scDECON_CANCER/data/MP/MP_counts_matrix_train.tab ../../../scDECON_CANCER/data/MP/ETM/cell_labels_oh_train.csv ../scMP/scMP_trainData.txt ../scMP/scMP_priorData.txt ../scMP/scMP_meta.txt ../scMP/scMP_genes.txt
//./SingleCellInput ../../../scDECON_CANCER/data/MP/MP_counts_matrix_pp_train.tab ../../../scDECON_CANCER/data/MP/ETM/cell_labels_oh_train.csv ../scMP/scMP_pp_trainData.txt ../scMP/scMP_pp_priorData.txt ../scMP/scMP_pp_meta.txt ../scMP/scMP_pp_genes.txt
//./SingleCellInput ../../../RA/other_datasets/HBC_NSR/counts_matrix_pp_train.tab ../../../RA/other_datasets/HBC_NSR/cell_labels_oh_pp_train.csv ../scHB/scHB_pp_trainData.txt ../scHB/scHB_pp_priorData.txt ../scHB/scHB_pp_meta.txt ../scHB/scHB_pp_genes.txt
//./SingleCellInput ../../../RA/other_datasets/HBC_NSR/labels_7celltypes/counts_matrix_pp_train.tab ../../../RA/other_datasets/HBC_NSR/labels_7celltypes/cell_labels_oh_pp_train.csv ../scHB_labels_7celltypes_with_trainmetaphe/scHB_pp_trainData.txt ../scHB_labels_7celltypes_with_trainmetaphe/scHB_pp_priorData.txt ../scHB_labels_7celltypes_with_trainmetaphe/scHB_pp_meta.txt ../scHB_labels_7celltypes_with_trainmetaphe/scHB_pp_genes.txt

#include <bits/stdc++.h> 
#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <armadillo>

using namespace std;
using namespace arma;

//Convert the C x G matrix (cells x genes) from single cell data into train_Data file format for MixEHR
//Convert the cells_ohe (cells - one hot encoding) file for the single cell data into prior_Data file format for MixEHR

int main(int argc, char* argv[]) {
	string line;
	ifstream single_cell_file (argv[1]);
	ifstream cell_ohe_file (argv[2]);

        string output_path = argv[3];

	int num_topics = atoi(argv[4]);
	int num_celltypes = atoi(argv[5]);

        string train_file = output_path + "trainData.txt";
        string prior_file = output_path + "priorData.txt";
        string meta_file = output_path + "meta.txt";
        string genes_file = output_path + "genes.txt";

        ofstream ofs_train_file;
	ofstream ofs_prior_file;
	ofstream ofs_meta_file;
	ofstream ofs_genes_file;

        ofs_train_file.open(train_file);
        ofs_prior_file.open(prior_file);
        ofs_meta_file.open(meta_file);
        ofs_genes_file.open(genes_file);

	// Vector of string to save tokens 
	vector <string> genes;
	int cell_id = 1;

	//Generate the training data file from single_cell_file
	if (single_cell_file.is_open() && ofs_train_file.is_open() && ofs_meta_file.is_open() && ofs_genes_file.is_open())
	{
		//Use the first line to store list of genes
		getline(single_cell_file,line);
	
		//stringstream class check1 
		stringstream check1(line); 
		string intermediate; 
		//Tokenizing w.r.t. tab  
		while(getline(check1, intermediate,'\t')) 
		{ 
			if(intermediate != "")
				genes.push_back(intermediate); 
		} 

		//Printing the token vector 
		for(int i = 0; i < genes.size(); i++)
			ofs_genes_file << genes[i] << '\n';

		while (getline(single_cell_file,line))
		{
			//stringstream class check1 
			stringstream check1(line); 
			string intermediate; 
			vector <string> tokens;
			//Tokenizing w.r.t. tab  
			while(getline(check1, intermediate,'\t')) 
			{ 
				tokens.push_back(intermediate); 
			} 

			//Generating train data for this cell as "Cell ID, Gene ID, Counts"
			//From 2nd cell on it is cell counts for each gene for that cell
			//Starting cell_IDs from 1
			//Starting gene_IDs from 1
			for(int i = 1; i < tokens.size(); i++)
				if(stoi(tokens[i]) > 0)
					ofs_train_file << cell_id << ' ' << '1' << ' ' << i << ' ' << '0' << ' ' << tokens[i] << '\n';

			if(cell_id == 1){
				for(int i = 1; i < tokens.size(); i++)
					ofs_meta_file << '1' << ' ' << i << ' ' << '0' << '\n';
			}

			cell_id++;
		}
		single_cell_file.close();
		ofs_train_file.close();
	}
	else cout << "Unable to open file"; 

	cell_id = 1;
	if(cell_ohe_file.is_open() && ofs_prior_file.is_open())
	{
		while (getline(cell_ohe_file,line))
		{
			//stringstream class check1 
			stringstream check1(line); 
			string intermediate; 
			vector <string> tokens;
			//Tokenizing w.r.t. ,  
			while(getline(check1, intermediate,',')) 
			{ 
				tokens.push_back(intermediate); 
			} 

			//Generating prior data for this cell as "Cell ID, Celltype ID, prior prob"
			//Starting cell_IDs from 1
			//Cell type IDs are the indices with a value of 1 in one hot encoding
			//Setting prior prob as 1 for those cell types with prob == 1
			//
			/*Option 1:
			double random_prob = 1.0/(tokens.size()-1);
			double highest_prob = 0.9;
			for(int i = 0; i < tokens.size(); i++)
				if(stoi(tokens[i]) > 0)
					prior_file << cell_id << ' ' << i << ' ' << highest_prob << '\n';
				else
					prior_file << cell_id << ' ' << i << ' ' << random_prob << '\n';
			*/
			

			/*Option 2:
			rowvec random_prob = randu<rowvec>(tokens.size());
			double highest_prob = 0.9;
			for(int i = 0; i < tokens.size(); i++)
				if(stoi(tokens[i]) > 0)
					prior_file << cell_id << ' ' << i << ' ' << highest_prob << '\n';
				else
					prior_file << cell_id << ' ' << i << ' ' << random_prob[i] << '\n';
			*/
			

			//Option 3: Random numbers between 0.001 and 0.01 for each topic
			//
			//
			double highest_prob = 0.0;
			double randNum = 0.0;
			if (num_topics == 1){

				//double highest_prob = 5;
				//double highest_prob = 10;
				//double highest_prob = 20;
				//double highest_prob = 50;
				//double highest_prob = 100;
				double highest_prob = 0.9;
				//double highest_prob = 0.8;
				//double highest_prob = 0.7;
				//double highest_prob = 0.95;
				//double highest_prob = 0.18;
				int j = 0;
				if(num_celltypes == 1){
					for(int i = 0; i < tokens.size(); i++)
					{
						if(stoi(tokens[i]) > 0)
							ofs_prior_file << cell_id << ' ' << j << ' ' << highest_prob << '\n';
						else{
							randNum = ((double)rand()/RAND_MAX)* (0.01-0.001) + 0.001;
							ofs_prior_file << cell_id << ' ' << j << ' ' << randNum << '\n';
						}
					}
					j++;
				}
				else if(num_celltypes == 2){
					for(int i = 0; i < tokens.size(); i++)
					{
						if(stoi(tokens[i]) > 0){
							ofs_prior_file << cell_id << ' ' << j << ' ' << highest_prob << '\n';
							ofs_prior_file << cell_id << ' ' << j+1 << ' ' << highest_prob << '\n';
						}
						else{
							randNum = ((double)rand()/RAND_MAX)* (0.01-0.001) + 0.001;
							ofs_prior_file << cell_id << ' ' << j << ' ' << randNum << '\n';
							ofs_prior_file << cell_id << ' ' << j+1 << ' ' << randNum << '\n';
						}
						j=j+2;
					}
				}
				else if(num_celltypes == 3){
					for(int i = 0; i < tokens.size(); i++)
					{
						if(stoi(tokens[i]) > 0){
							ofs_prior_file << cell_id << ' ' << j << ' ' << highest_prob << '\n';
							ofs_prior_file << cell_id << ' ' << j+1 << ' ' << highest_prob << '\n';
							ofs_prior_file << cell_id << ' ' << j+2 << ' ' << highest_prob << '\n';
						}
						else{
							randNum = ((double)rand()/RAND_MAX)* (0.01-0.001) + 0.001;
							ofs_prior_file << cell_id << ' ' << j << ' ' << randNum << '\n';
							ofs_prior_file << cell_id << ' ' << j+1 << ' ' << randNum << '\n';
							ofs_prior_file << cell_id << ' ' << j+2 << ' ' << randNum << '\n';
						}
						j=j+3;
					}
				}
				else if(num_celltypes == 4){
					for(int i = 0; i < tokens.size(); i++)
					{
						if(stoi(tokens[i]) > 0){
							ofs_prior_file << cell_id << ' ' << j << ' ' << highest_prob << '\n';
							ofs_prior_file << cell_id << ' ' << j+1 << ' ' << highest_prob << '\n';
							ofs_prior_file << cell_id << ' ' << j+2 << ' ' << highest_prob << '\n';
							ofs_prior_file << cell_id << ' ' << j+3 << ' ' << highest_prob << '\n';
						}
						else{
							randNum = ((double)rand()/RAND_MAX)* (0.01-0.001) + 0.001;
							ofs_prior_file << cell_id << ' ' << j << ' ' << randNum << '\n';
							ofs_prior_file << cell_id << ' ' << j+1 << ' ' << randNum << '\n';
							ofs_prior_file << cell_id << ' ' << j+2 << ' ' << randNum << '\n';
							ofs_prior_file << cell_id << ' ' << j+3 << ' ' << randNum << '\n';
						}
						j=j+4;
					}
				}
				else if(num_celltypes == 5){
					for(int i = 0; i < tokens.size(); i++)
					{
						if(stoi(tokens[i]) > 0){
							ofs_prior_file << cell_id << ' ' << j << ' ' << highest_prob << '\n';
							ofs_prior_file << cell_id << ' ' << j+1 << ' ' << highest_prob << '\n';
							ofs_prior_file << cell_id << ' ' << j+2 << ' ' << highest_prob << '\n';
							ofs_prior_file << cell_id << ' ' << j+3 << ' ' << highest_prob << '\n';
							ofs_prior_file << cell_id << ' ' << j+4 << ' ' << highest_prob << '\n';
						}
						else{
							randNum = ((double)rand()/RAND_MAX)* (0.01-0.001) + 0.001;
							ofs_prior_file << cell_id << ' ' << j << ' ' << randNum << '\n';
							ofs_prior_file << cell_id << ' ' << j+1 << ' ' << randNum << '\n';
							ofs_prior_file << cell_id << ' ' << j+2 << ' ' << randNum << '\n';
							ofs_prior_file << cell_id << ' ' << j+3 << ' ' << randNum << '\n';
							ofs_prior_file << cell_id << ' ' << j+4 << ' ' << randNum << '\n';
						}
						j=j+5;
					}
				}
			}
			else if (num_topics == 2){
	                        highest_prob = 0.45;
	                        int j=0;
	                        for(int i = 0; i < tokens.size(); i=i+1)
	                        {
	                               if(stoi(tokens[i]) > 0)
				       {
					       //highest_prob = ((double)rand()/RAND_MAX)* (0.46-0.44) + 0.44;
				               ofs_prior_file << cell_id << ' ' << j << ' ' << highest_prob << '\n';
					       //highest_prob = ((double)rand()/RAND_MAX)* (0.46-0.44) + 0.44;
				               ofs_prior_file << cell_id << ' ' << j+1 << ' ' << highest_prob << '\n';
				       }
				       else
				       {
					       //randNum = (0.1)/(tokens.size()*2 - 2);
					       randNum = ((double)rand()/RAND_MAX)* (0.01-0.001) + 0.001;
				               ofs_prior_file << cell_id << ' ' << j << ' ' << randNum << '\n';
					       //randNum = ((double)rand()/RAND_MAX)* (0.01-0.001) + 0.001;
					       ofs_prior_file << cell_id << ' ' << j+1 << ' ' << randNum << '\n';
				       }
				       j=j+2;
				}
			}
			else if (num_topics == 3){
	                        highest_prob = 0.3;
	                        int j=0;
	                        for(int i = 0; i < tokens.size(); i=i+1)
	                        {
	                               if(stoi(tokens[i]) > 0)
				       {
					       //highest_prob = ((double)rand()/RAND_MAX)* (0.34-0.32) + 0.32;
				               ofs_prior_file << cell_id << ' ' << j << ' ' << highest_prob << '\n';
					       //highest_prob = ((double)rand()/RAND_MAX)* (0.34-0.32) + 0.32;
				               ofs_prior_file << cell_id << ' ' << j+1 << ' ' << highest_prob << '\n';
					       //highest_prob = ((double)rand()/RAND_MAX)* (0.34-0.32) + 0.32;
				               ofs_prior_file << cell_id << ' ' << j+2 << ' ' << highest_prob << '\n';
				       }
				       else
				       {
					       //randNum = (1.0 - 0.9)/(tokens.size()*3 - 3);
					       randNum = ((double)rand()/RAND_MAX)* (0.01-0.001) + 0.001;
				               ofs_prior_file << cell_id << ' ' << j << ' ' << randNum << '\n';
					       //randNum = ((double)rand()/RAND_MAX)* (0.01-0.001) + 0.001;
					       ofs_prior_file << cell_id << ' ' << j+1 << ' ' << randNum << '\n';
					       //randNum = ((double)rand()/RAND_MAX)* (0.01-0.001) + 0.001;
				               ofs_prior_file << cell_id << ' ' << j+2 << ' ' << randNum << '\n';
				       }
				       j=j+3;
				}
			}
			else if (num_topics == 4){
	                        highest_prob = 0.225;
	                        int j=0;
	                        for(int i = 0; i < tokens.size(); i=i+1)
	                        {
	                               if(stoi(tokens[i]) > 0)
				       {
					       //highest_prob = ((double)rand()/RAND_MAX)* (0.235-0.215) + 0.215;
				               ofs_prior_file << cell_id << ' ' << j << ' ' << highest_prob << '\n';
					       //highest_prob = ((double)rand()/RAND_MAX)* (0.235-0.215) + 0.215;
				               ofs_prior_file << cell_id << ' ' << j+1 << ' ' << highest_prob << '\n';
					       //highest_prob = ((double)rand()/RAND_MAX)* (0.235-0.215) + 0.215;
				               ofs_prior_file << cell_id << ' ' << j+2 << ' ' << highest_prob << '\n';
					       //highest_prob = ((double)rand()/RAND_MAX)* (0.235-0.215) + 0.215;
				               ofs_prior_file << cell_id << ' ' << j+3 << ' ' << highest_prob << '\n';
				       }
				       else
				       {
					       //randNum = (1.0 - 0.9)/(tokens.size()*4 - 4);
					       randNum = ((double)rand()/RAND_MAX)* (0.01-0.001) + 0.001;
				               ofs_prior_file << cell_id << ' ' << j << ' ' << randNum << '\n';
					       //randNum = ((double)rand()/RAND_MAX)* (0.01-0.001) + 0.001;
					       ofs_prior_file << cell_id << ' ' << j+1 << ' ' << randNum << '\n';
					       //randNum = ((double)rand()/RAND_MAX)* (0.01-0.001) + 0.001;
				               ofs_prior_file << cell_id << ' ' << j+2 << ' ' << randNum << '\n';
					       //randNum = ((double)rand()/RAND_MAX)* (0.01-0.001) + 0.001;
				               ofs_prior_file << cell_id << ' ' << j+3 << ' ' << randNum << '\n';
				       }
				       j=j+4;
				}
			}
			else if (num_topics == 5){
	                        highest_prob = 0.18;
	                        int j=0;
	                        for(int i = 0; i < tokens.size(); i=i+1)
	                        {
	                               if(stoi(tokens[i]) > 0)
				       {
					       //highest_prob = ((double)rand()/RAND_MAX)* (0.19-0.17) + 0.17;
				               ofs_prior_file << cell_id << ' ' << j << ' ' << highest_prob << '\n';
					       //highest_prob = ((double)rand()/RAND_MAX)* (0.19-0.17) + 0.17;
				               ofs_prior_file << cell_id << ' ' << j+1 << ' ' << highest_prob << '\n';
					       //highest_prob = ((double)rand()/RAND_MAX)* (0.19-0.17) + 0.17;
				               ofs_prior_file << cell_id << ' ' << j+2 << ' ' << highest_prob << '\n';
					       //highest_prob = ((double)rand()/RAND_MAX)* (0.19-0.17) + 0.17;
				               ofs_prior_file << cell_id << ' ' << j+3 << ' ' << highest_prob << '\n';
					       //highest_prob = ((double)rand()/RAND_MAX)* (0.19-0.17) + 0.17;
				               ofs_prior_file << cell_id << ' ' << j+4 << ' ' << highest_prob << '\n';
				       }
				       else
				       {
					       //randNum = (1.0 - 0.9)/(tokens.size()*5 - 5);
					       randNum = ((double)rand()/RAND_MAX)* (0.01-0.001) + 0.001;
					       //randNum = ((double)rand()/RAND_MAX)* (0.01-0.001) + 0.001;
				               ofs_prior_file << cell_id << ' ' << j << ' ' << randNum << '\n';
					       //randNum = ((double)rand()/RAND_MAX)* (0.01-0.001) + 0.001;
					       ofs_prior_file << cell_id << ' ' << j+1 << ' ' << randNum << '\n';
					       //randNum = ((double)rand()/RAND_MAX)* (0.01-0.001) + 0.001;
				               ofs_prior_file << cell_id << ' ' << j+2 << ' ' << randNum << '\n';
					       //randNum = ((double)rand()/RAND_MAX)* (0.01-0.001) + 0.001;
				               ofs_prior_file << cell_id << ' ' << j+3 << ' ' << randNum << '\n';
					       //randNum = ((double)rand()/RAND_MAX)* (0.01-0.001) + 0.001;
				               ofs_prior_file << cell_id << ' ' << j+4 << ' ' << randNum << '\n';
				       }
				       j=j+5;
				}
			}
			else if (num_topics == 10){
	                        highest_prob = 0.09;
	                        int j=0;
	                        for(int i = 0; i < tokens.size(); i=i+1)
	                        {
	                               if(stoi(tokens[i]) > 0)
				       {
				               ofs_prior_file << cell_id << ' ' << j << ' ' << highest_prob << '\n';
				               ofs_prior_file << cell_id << ' ' << j+1 << ' ' << highest_prob << '\n';
				               ofs_prior_file << cell_id << ' ' << j+2 << ' ' << highest_prob << '\n';
				               ofs_prior_file << cell_id << ' ' << j+3 << ' ' << highest_prob << '\n';
				               ofs_prior_file << cell_id << ' ' << j+4 << ' ' << highest_prob << '\n';
				               ofs_prior_file << cell_id << ' ' << j+5 << ' ' << highest_prob << '\n';
				               ofs_prior_file << cell_id << ' ' << j+6 << ' ' << highest_prob << '\n';
				               ofs_prior_file << cell_id << ' ' << j+7 << ' ' << highest_prob << '\n';
				               ofs_prior_file << cell_id << ' ' << j+8 << ' ' << highest_prob << '\n';
				               ofs_prior_file << cell_id << ' ' << j+9 << ' ' << highest_prob << '\n';
					}
				       else
					{
					       randNum = ((double)rand()/RAND_MAX)* (0.001-0.0001) + 0.0001;
				               ofs_prior_file << cell_id << ' ' << j << ' ' << randNum << '\n';
				               ofs_prior_file << cell_id << ' ' << j+1 << ' ' << randNum << '\n';
				               ofs_prior_file << cell_id << ' ' << j+2 << ' ' << randNum << '\n';
				               ofs_prior_file << cell_id << ' ' << j+3 << ' ' << randNum << '\n';
				               ofs_prior_file << cell_id << ' ' << j+4 << ' ' << randNum << '\n';
				               ofs_prior_file << cell_id << ' ' << j+5 << ' ' << randNum << '\n';
				               ofs_prior_file << cell_id << ' ' << j+6 << ' ' << randNum << '\n';
				               ofs_prior_file << cell_id << ' ' << j+7 << ' ' << randNum << '\n';
				               ofs_prior_file << cell_id << ' ' << j+8 << ' ' << randNum << '\n';
				               ofs_prior_file << cell_id << ' ' << j+9 << ' ' << randNum << '\n';
					}
			       		j=j+10;	       
				}
			}
			else{
				cout << "Error: Please enter topic number between 1 - 5" << endl;
			}
			cell_id++;

		}
		cell_ohe_file.close();
		ofs_prior_file.close();
	}


    return 0;
}
