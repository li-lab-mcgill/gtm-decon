//Usage: 
//g++ SingleCellInput_TestData.cpp -o SingleCellInput_TestData
//./SingleCellInput_TestData ../../../scDECON_CANCER/data/MP/MP_counts_matrix_test.tab ../scMP/scMP_testData.txt

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
	ofstream test_data_file (argv[2]);
	ofstream ofs_meta_file (argv[3]);
	ofstream ofs_genes_file (argv[4]);
	// Vector of string to save tokens 
	vector <string> genes;
	int cell_id = 1;

	//Generate the training data file from single_cell_file
	if (single_cell_file.is_open() && test_data_file.is_open() && ofs_meta_file.is_open())
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
					test_data_file << cell_id << ' ' << '1' << ' ' << i << ' ' << '0' << ' ' << tokens[i] << '\n';
					//test_data_file << cell_id << ' ' << '1' << ' ' << i << ' ' << '0' << ' ' << tokens[i] << '\n';
				//else
				//	test_data_file << cell_id << ' ' << '1' << ' ' << i << ' ' << '0' << ' ' << tokens[i] << ' ' << '0' << '\n';
				//
                        if(cell_id == 1){
                                for(int i = 1; i < tokens.size(); i++)
                                        ofs_meta_file << '1' << ' ' << i << ' ' << '0' << '\n';
                        }

			cell_id++;
		}
		single_cell_file.close();
		test_data_file.close();
		ofs_meta_file.close();
	}
	else cout << "Unable to open file"; 


    return 0;
}
