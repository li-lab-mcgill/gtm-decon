//Usage: 
//g++ Accuracy.cpp -o Accuracy 
//./Accuracy ../scHB_labels_ordered/scHB_pp_testData_scHB_pp_trainData_JCVB0_nmar_K32_iter500_metaphe.csv ../../../RA/other_datasets/HBC_NSR/labels_ordered/cell_labels_oh_pp_test.csv
//./Accuracy ../scHB_labels_7celltypes/scHB_pp_testData_scHB_pp_trainData_JCVB0_nmar_K7_iter500_metaphe.csv ../../../RA/other_datasets/HBC_NSR/labels_7celltypes/cell_labels_oh_pp_test.csv

#include <bits/stdc++.h> 
#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <armadillo>

using namespace std;
using namespace arma;

//Calculate the accuracy between "_Metaphe.csv" file and "_cell_ohe_test.csv" files by comparing whether the indices of the max elements are the same - ie. does the index of the topic with the highest value correspond to the index of the cell type (in a guided topic model version)

int main(int argc, char* argv[]) {
	string line;
	ifstream metaphe_file (argv[1]);
	ifstream cell_ohe_file (argv[2]);

	vector <int> indices_metaphe;
	vector <int> indices_cell_ohe;

	vector <string> metaphe;

	//Generate list of max indices for metaphe_file 
	if (metaphe_file.is_open())
	{
		while (getline(metaphe_file,line))
		{
			metaphe.push_back(line);

			//stringstream class check1 
			stringstream check1(line); 
			string intermediate; 
			vector <string> tokens;
			//Tokenizing w.r.t. comma  
			while(getline(check1, intermediate,',')) 
			{ 
				tokens.push_back(intermediate); 
			}

			//identify the index of the element with the highest value
			int max_index = -1;
			double max_val = -1.0;
			for(int i = 0; i < tokens.size(); i++)
			{
				if(stod(tokens[i]) > max_val)
				{
					max_index = i;
					max_val = stod(tokens[i]);
				}
			}

			//Insert the max index of the cell into the vector
			indices_metaphe.push_back(max_index);
		}
		metaphe_file.close();
	}
	else cout << "Unable to open file"; 

	//Generate list of max indices for metaphe_file 
	if (cell_ohe_file.is_open())
	{
		while (getline(cell_ohe_file,line))
		{
			//stringstream class check1 
			stringstream check1(line); 
			string intermediate; 
			vector <string> tokens;
			//Tokenizing w.r.t. comma  
			while(getline(check1, intermediate,',')) 
			{ 
				tokens.push_back(intermediate); 
			}

			//Identify the index of the element with the highest value
			int max_index = -1;
			double max_val = -1.0;
			for(int i = 0; i < tokens.size(); i++)
			{
				if(stod(tokens[i]) > max_val)
				{
					max_index = i;
					max_val = stod(tokens[i]);
				}
			}

			//Insert the max index of the cell into the vector
			indices_cell_ohe.push_back(max_index);
		}
		cell_ohe_file.close();
	}
	else cout << "Unable to open file"; 

	//Calculate how similar the indices_metaphe and indices_cell_ohe are.. to calculate accuracy of prediction (considering indices_cell_ohe as ground truth
	
	if ( indices_cell_ohe.size() != indices_metaphe.size() )
	{
		cout << "Something wrong.. the number of elements in metaphe file NOT same as cell_ohe file. Please check" << endl;
		exit(0);
	}

	int matched_indices = 0;
	for (int i = 0; i < indices_cell_ohe.size(); i++)
	{
		if ( indices_cell_ohe[i] == indices_metaphe[i] )
		{
			//cout << indices_metaphe[i] << endl;
			matched_indices++;
			//cout << metaphe[i] << endl;
		}
	}
	cout << "Matched indices: " << matched_indices << endl;
	cout << "Total_number_of_cells: " << indices_cell_ohe.size() << endl;
	cout << "Accuracy: " << (matched_indices*100)/indices_cell_ohe.size() << "%" << endl;

    return 0;
}
