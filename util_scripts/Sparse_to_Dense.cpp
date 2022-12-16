//Convert a sparse matrix to a dense matrix
//Input in MatrixMarket matrix format <genes> <cells> <counts>
//Output in <cells x genes format> as a normal matrix
// ../../gtm-decon/Sparse_to_Dense genes.read.txt cells.read.new.txt counts.read.txt meta_celltypes.txt > chk.txt
//
//Usage: 

#include <bits/stdc++.h> 
#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <armadillo>

using namespace std;
using namespace arma;

int main(int argc, char* argv[]) {
	string line;
	ifstream gene_file(argv[1]);
	ifstream cell_file(argv[2]);
	ifstream MM_file(argv[3]);
	ifstream meta_file(argv[4]);

	vector<string> cell_names;
	vector<string> cell_types;
	vector<string> cell_names_types;
        vector<vector<int>> dense_matrix(50000, vector<int> (50000));

	//Parse gene file and output on a single line
	cout << "Unnamed: 0\t" ;
	if (gene_file.is_open()){
		while(getline(gene_file,line)){
			stringstream check1(line); 
			string intermediate; 
			string gene_symbol;
			while(getline(check1,intermediate,'_')){
				gene_symbol = intermediate;
			}
			cout << gene_symbol << "\t";
		}
	}
	cout << "\n";

	//Parse cell names and output values
	if(cell_file.is_open()){
		while(getline(cell_file,line)){
			cell_names.push_back(line);
		}
	}

	//Parse cell types and output values
	if(meta_file.is_open()){
		while(getline(meta_file,line)){
			stringstream check1(line); 
			string intermediate; 
			vector <string> tokens;
			while(getline(check1,intermediate,',')){
				tokens.push_back(intermediate);
			}
			cell_names_types.push_back(tokens[0]);
			cell_types.push_back(tokens[1]);
		}
	}

	/*for(int i=0; i<cell_names.size(); i++){
		cout << cell_names[i] << endl;
	}*/

	//Parse MM file and output all values in matrix format
	//
	
	int num_genes = 0;
	int num_cells = 0;
	if (MM_file.is_open())
	{
		//Skip first line as it contains comment info 
		getline(MM_file,line);
		//2nd line contains info about number of genes and number of cells
		getline(MM_file,line);
		stringstream check1(line); 
		string intermediate;
		vector <int> tokens;
	       	
		while(getline(check1, intermediate,' ')){
			tokens.push_back(stoi(intermediate));
		}

		num_genes = tokens[0]; 
		num_cells = tokens[1]; 
		//cout << num_genes << "," << num_cells << endl;

		int gene_num = 1;
		int cell_num = 0;
		int meta_num = 0;
		int prev_cell_num = 0;
		while (getline(MM_file,line))
		{
			//stringstream class check1 
			stringstream check1(line); 
			string intermediate; 
			vector <int> tokens;
			//Tokenizing w.r.t. space  
			while(getline(check1, intermediate,' ')) 
			{ 
				tokens.push_back(stoi(intermediate)); 
			}
			
			if (tokens[1] > prev_cell_num){
				//cout << "Here" << endl;
				if(prev_cell_num > 0){
					for(int j = gene_num; j <= num_genes; j++)
						cout << "0\t";
					if(cell_names[prev_cell_num-1] == cell_names_types[meta_num]){
						cout << cell_types[meta_num] << endl;
						//cout << cell_names[prev_cell_num-1] << "-" << cell_names_types[meta_num] << endl;
						meta_num++;
					}
					else{
						cout << "Not_specified" << endl;
						//cout << cell_names[prev_cell_num-1] << "," << cell_names_types[meta_num] << endl;
					}
				}
				cout << cell_names[cell_num] << "\t";
				cell_num++;
				prev_cell_num = cell_num;
				gene_num = 1;
			}

			for(int j = gene_num; j < tokens[0]; j++)
				cout << "0\t";
			cout << tokens[2] << "\t";
			gene_num = tokens[0] + 1;
			
		}
		//For the last cell
		for(int j = gene_num; j <= num_genes; j++)
			cout << "0\t";
		if(cell_names[prev_cell_num-1] == cell_names_types[meta_num]){
			cout << cell_types[meta_num] << endl;
			//cout << cell_names[prev_cell_num-1] << "-" << cell_names_types[meta_num] << endl;
		}
		else{
			cout << "Not_specified" << endl;
			//cout << cell_names[prev_cell_num-1] << "," << cell_names_types[meta_num] << endl;
		}
	}
    return 0;
}
