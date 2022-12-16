//Usage: 
//../../../../../MixEHR-sureLDA-GTM/Phi_3topics artificial_bulkData_trainData_JCVB0_nmar_K42_iter9_metaphe.csv > one

#include <bits/stdc++.h> 
#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <armadillo>

using namespace std;
using namespace arma;

//Generate a new Phi file from the Phi file for 3 topics by summing the values of 3 consecutive topics (which corresponds to a cell type)


int main(int argc, char* argv[]) {
	string line;
	ifstream metaphe_file (argv[1]);
	int num_topics = stoi(argv[2]);

	//Generate list of max indices for metaphe_file 
	if (metaphe_file.is_open())
	{
		while (getline(metaphe_file,line))
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

			//Calculate summe of 2 consecutive topics - corresponding to each cell type, and then identify max index 

			double tot_val = 0;
			vector<double> avg_vals;
			if(num_topics == 2){
				for(int i = 0; i < tokens.size(); i=i+2)
				{
					tot_val = (stod(tokens[i]) + stod(tokens[i+1]));
					//avg_vals.push_back(tot_val/2.0);
					avg_vals.push_back(tot_val);
				}
			}
			else if(num_topics == 3){
				for(int i = 0; i < tokens.size(); i=i+3)
				{
					tot_val = (stod(tokens[i]) + stod(tokens[i+1]) + stod(tokens[i+2]));
					//avg_vals.push_back(tot_val/3.0);
					avg_vals.push_back(tot_val);
				}
			}
			else if(num_topics == 4){
				for(int i = 0; i < tokens.size(); i=i+4)
				{
					tot_val = (stod(tokens[i]) + stod(tokens[i+1]) + stod(tokens[i+2]) + stod(tokens[i+3]));
					//avg_vals.push_back(tot_val/4.0);
					avg_vals.push_back(tot_val);
				}
			}
			else if(num_topics == 5){
				for(int i = 0; i < tokens.size(); i=i+5)
				{
					tot_val = (stod(tokens[i]) + stod(tokens[i+1]) + stod(tokens[i+2]) + stod(tokens[i+3]) + stod(tokens[i+4]));
					//avg_vals.push_back(tot_val/5.0);
					avg_vals.push_back(tot_val);
				}
			}
			else if(num_topics == 10){
				for(int i = 0; i < tokens.size(); i=i+10)
				{
					tot_val = (stod(tokens[i]) + stod(tokens[i+1]) + stod(tokens[i+2]) + stod(tokens[i+3]) + stod(tokens[i+4]) + stod(tokens[i+5]) + stod(tokens[i+6]) + stod(tokens[i+7]) + stod(tokens[i+8]) + stod(tokens[i+9]));
					//avg_vals.push_back(tot_val/10.0);
					avg_vals.push_back(tot_val);
				}
			}
			else{
                                cout << "Check your number of topics: Should be 2,3,4,5,10" << endl;
                        }

			for(int i = 0; i < avg_vals.size(); i=i+1)
			{
				//cout << avg_vals[i]/tot_avg_vals << "," ;
				cout << avg_vals[i] << "," ;
			}
			cout << endl;
		}
		metaphe_file.close();
	}
	else cout << "Unable to open file"; 

    return 0;
}
