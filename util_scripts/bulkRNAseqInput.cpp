//Usage: 
//g++ BulkRNAseqInput.cpp -o BulkRNAseqInput
//./BulkRNAseqInput ~/projects/RA/main_datasets/MH.SARDs.RNASeq.DC.csv  ../scHB_labels_7celltypes/MH.SARDs.RNASeq_BulkData.txt

#include <bits/stdc++.h> 
#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <armadillo>

using namespace std;
using namespace arma;

//Convert the C x G matrix (cells x genes) from bulk RNAseq data into New Patients Data file format for MixEHR

int main(int argc, char* argv[]) {
	string line;
	ifstream bulk_RNAseq_file (argv[1]);
	ofstream test_data_file (argv[2]);
	// Vector of string to save tokens 
	vector <string> genes;
	int cell_id = 1;

	//Generate the training data file from bulk_RNAseq_file
	if (bulk_RNAseq_file.is_open() && test_data_file.is_open())
	{
		//Use the first line to store list of genes
		getline(bulk_RNAseq_file,line);
	
		//stringstream class check1 
		stringstream check1(line); 
		string intermediate; 
		//Tokenizing w.r.t. tab  
		while(getline(check1, intermediate,'\t')) 
		{ 
			if(intermediate != "")
				genes.push_back(intermediate); 
		} 

		while (getline(bulk_RNAseq_file,line))
		{
			//stringstream class check1 
			stringstream check1(line); 
			string intermediate; 
			vector <string> tokens;
			//Tokenizing w.r.t. tab  
			while(getline(check1, intermediate,'\t')) 
			{ 
				if(intermediate != "")
					tokens.push_back(intermediate); 
				else
					tokens.push_back("0");
			}

			//Generating train data for this cell as "Cell ID, Gene ID, Counts"
			//From 2nd cell on it is cell counts for each gene for that cell
			//Starting cell_IDs from 1
			//Starting gene_IDs from 1
			for(int i = 1; i < tokens.size(); i++)
			{
				if(stoi(tokens[i]) > 0)
					test_data_file << cell_id << ' ' << '1' << ' ' << i << ' ' << '0' << ' ' << tokens[i] << '\n';
			}
			cell_id++;
		}
		bulk_RNAseq_file.close();
		test_data_file.close();
	}
	else cout << "Unable to open file"; 


    return 0;
}
