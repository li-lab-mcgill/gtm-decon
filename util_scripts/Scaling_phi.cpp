//Usage: 
//g++ Scaling_phi.cpp -o Scaling_phi

//This script identifies what the scaling factor should be for "immune topics" vs "tissue topics" and scales the immune  topics' phi values accordingly
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
	ifstream phi_file (argv[1]);
	int tissue_topics = atoi(argv[2]);
	ofstream new_phi_file (argv[3]);

        vector<vector<double>> phi_matrix(100000, vector<double> (500));

	//Parse phi_file and input all the results
	int i=0;
	int j=2;
	if (phi_file.is_open())
	{
		while (getline(phi_file,line))
		{
			//stringstream class check1 
			stringstream check1(line); 
			string intermediate; 
			vector <double> tokens;
			//Tokenizing w.r.t. comma  
			while(getline(check1, intermediate,',')) 
			{ 
				tokens.push_back(stod(intermediate)); 
			}

                        for (j = 2; j < tokens.size(); j++) //skip first 2 cols as they are not phi values
                        {
                                phi_matrix[i][j-2]=tokens[j];
			}

			i++;
		}
	}
	else cout << "Unable to open phi file"; 
	
	int rows_phi = i;
	int total_topics = j-2;

	double overall_immune_m = 0.0;
	for(int x = tissue_topics; x < total_topics; x++)
	{
		double m_col = 0.0;	//Calculate the average ie. m of each column
		for(int y = 0; y < rows_phi; y++)
		{	
			m_col =  m_col + phi_matrix[y][x];
		}
		
		cout << x << "," << m_col/rows_phi << endl;

		overall_immune_m = overall_immune_m + (m_col/rows_phi);
	}

	overall_immune_m = overall_immune_m / (total_topics - tissue_topics);

	cout << total_topics - tissue_topics << endl;
	cout << "Average immune topics:" << overall_immune_m  << endl;

	double overall_tissue_m = 0.0;
	for(int x = 0; x < tissue_topics; x++)
	{
		double m_col = 0.0;	//Calculate the average ie. m of each column
		for(int y = 0; y < rows_phi; y++)
		{	
			m_col =  m_col + phi_matrix[y][x];
		}
		
		cout << x << "," << m_col/rows_phi << endl;

		overall_tissue_m = overall_tissue_m + (m_col/rows_phi);
	}

	overall_tissue_m = overall_tissue_m / (tissue_topics);

	cout << tissue_topics << endl;
	cout << "Average tissue topics:" << overall_tissue_m  << endl;

	double scaling_factor = overall_tissue_m / overall_immune_m;

	cout << "Scaling factor:" << scaling_factor << endl;

	//Regenerate phi file based on the scaling factor - Multiply immune topics by that scaling factor

	i = 1;	
	for(int y = 0; y < rows_phi; y++)
	{
		new_phi_file << "1," << i << ",";
		/*for(int x = 0; x < tissue_topics; x++)
		{
			new_phi_file << phi_matrix[y][x] * scaling_factor <<",";
		}
		for(int x = tissue_topics; x < total_topics-1; x++)
		{
			new_phi_file << phi_matrix[y][x] << ",";
		}
		new_phi_file << phi_matrix[y][total_topics-1] << endl;*/

		for(int x = 0; x < tissue_topics; x++)
		{
			new_phi_file << phi_matrix[y][x] <<",";
		}
		for(int x = tissue_topics; x < total_topics-1; x++)
		{
			new_phi_file << phi_matrix[y][x] * scaling_factor << ",";
		}
		new_phi_file << phi_matrix[y][total_topics-1] * scaling_factor << endl;

		i++;
	}
	//cout << total_topics << endl;

    return 0;
}
