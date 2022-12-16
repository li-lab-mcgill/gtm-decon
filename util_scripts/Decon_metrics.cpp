//Usage: 
//g++ Deconvolution_metrics.cpp -o Deconvolution_metrics
//./Deconvolution_metrics <simulated_cell_proportions> ../scHB_labels_7celltypes/scHB_pp_testData_scHB_pp_trainData_JCVB0_nmar_K7_iter500_metaphe.csv
//src/Deconvolution_metrics results/simulated_bulkRS/scHB/scHB_test_origcellprop_all_proportions_sets.csv results/scHB/scHB_tvt_pp_7celltypes/scHB_pp_origcellprop_simbulkRS_scHB_pp_trainData_JCVB0_nmar_K7_iter1_metaphe.csv

#include <bits/stdc++.h> 
#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <armadillo>

using namespace std;
using namespace arma;

//Calculate two deconvolution metrics to compare cell-type proportions - ground truth and _metaphe.csv (normalized) file ie. theta
//Pearson Correlation Coefficient (PCC), Root Mean Square Estimation (RMSE)

//function that returns correlation coefficient. 
//Taken from https://www.geeksforgeeks.org/program-find-correlation-coefficient/
double correlationCoefficient(vector<double> X, vector<double> Y, int n) 
{ 
	double sum_X = 0, sum_Y = 0, sum_XY = 0; 
      	double squareSum_X = 0, squareSum_Y = 0; 
        
	for (int i = 0; i < n; i++) 
        { 
   		//cout << X[i] << "," << Y[i] << endl;
       		// sum of elements of array
		sum_X = sum_X + X[i]; 
                         
		// sum of elements of array Y
		sum_Y = sum_Y + Y[i]; 
                                           
		// sum of X[i] * Y[i]. 
		sum_XY = sum_XY + X[i] * Y[i]; 
                                                                           
                // sum of square of array elements. 
                squareSum_X = squareSum_X + X[i] * X[i]; 
                squareSum_Y = squareSum_Y + Y[i] * Y[i]; 
        }
       	//cout << endl;	
                                                                                                         
        // use formula for calculating correlation coefficient. 
        double corr = (float)(n * sum_XY - sum_X * sum_Y) / sqrt((n * squareSum_X - sum_X * sum_X) * (n * squareSum_Y - sum_Y * sum_Y));

        return corr; 
}

//function that returns mean absolute deviation

double meanAbsoluteDeviation(vector<double> X, vector<double> Y, int n)
{
	double diff = 0.0;
	for(int i = 0; i < n; i++)
	{
		diff = diff + abs(X[i] - Y[i]);
	}

	return diff/n;
}

//function that returns root mean square deviation  (RMSD)
double rootMeanSquareDeviation(vector<double> X, vector<double> Y, int n)
{
	double diff = 0.0;
	for(int i = 0; i < n; i++)
	{
		diff = diff + (X[i] - Y[i]) * (X[i] - Y[i]);
	}
	return sqrt(diff/n);
}

// Function for calculating median 
double findMedian(vector<double> a) 
{ 
     // First we sort the array 
     sort(a.begin(),a.end()); 
     int n = a.size();

     // check for even case 
     if (n % 2 != 0) 
 	    return (double)a[n / 2]; 
                             
     return (double)(a[(n - 1) / 2] + a[n / 2]) / 2.0; 
} 

// Function for calculating min 
double findMin(vector<double> a) 
{ 
     // First we sort the array 
     sort(a.begin(),a.end()); 
     int n = a.size();

     return (double)a[0]; 
}
                             
int main(int argc, char* argv[]) {
	string line;
	ifstream cell_proportions_file (argv[1]);
	ifstream metaphe_file (argv[2]);

        vector<vector<double>> proportions_metaphe(1000, vector<double> (50));
        vector<vector<double>> proportions_ground_truth(1000, vector<double> (50));

	int topics=0;

	//Parse ground_truth_file and input all the results
	int i=0;
	if (cell_proportions_file.is_open())
	{
		//Skip first line as it contains celltype names
		getline(cell_proportions_file,line);
		while (getline(cell_proportions_file,line))
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

                        for (int j = 0; j < tokens.size(); j++)
                        {
                                proportions_ground_truth[i][j]=tokens[j];
			}

			i++;
		}
	}
	int rows_proportions = i;

	//Parse metaphe_file and input all the results
	i=0;
	if (metaphe_file.is_open())
	{
		while (getline(metaphe_file,line))
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

                        for (int j = 0; j < tokens.size(); j++)
                        {
                                proportions_metaphe[i][j]=tokens[j];
                        }

			i++;
                        topics = tokens.size();
	         }
		 metaphe_file.close();
	}
	else cout << "Unable to open metaphe file"; 

	int rows_metaphe = i;

	//Calculate how similar the indices_metaphe and indices_cell_ohe are.. to calculate accuracy of prediction (considering indices_cell_ohe as ground truth

	if ( rows_proportions != rows_metaphe )
	{
		cout << "Something wrong.. the number of elements in metaphe file NOT same as ground_truth file. Please check" << endl;
		exit(0);
	}

	vector<double> all_corr;
	vector<double> all_mad;
	vector<double> all_rmsd;

	double avg_corr = 0.0;
	double avg_mad = 0.0;
	double avg_rmsd = 0.0;

	for(int x = 0; x < rows_metaphe; x++)
	{
		std::vector<std::vector<double>>::iterator row=proportions_metaphe.begin() + x;
		std::vector<double>::iterator start=row->begin();
		std::vector<double>::iterator end=row->end(); 
		std::vector<double> proportions_metaphe_1D(start,end);

		std::vector<std::vector<double>>::iterator row1=proportions_ground_truth.begin() + x;
		std::vector<double>::iterator start1=row1->begin();
		std::vector<double>::iterator end1=row1->end(); 
		std::vector<double> proportions_ground_truth_1D(start1,end1);

	        //Function call to correlationCoefficient. 
                float corr = correlationCoefficient(proportions_metaphe_1D, proportions_ground_truth_1D, topics); 
                float mad = meanAbsoluteDeviation(proportions_metaphe_1D, proportions_ground_truth_1D, topics); 
                float rmsd = rootMeanSquareDeviation(proportions_metaphe_1D, proportions_ground_truth_1D, topics); 

		cout << corr << "," << mad << "," << rmsd << endl;

		all_corr.push_back(corr);
		all_mad.push_back(mad);
		all_rmsd.push_back(rmsd);

		avg_corr = avg_corr + corr;
		avg_mad = avg_mad + mad;
		avg_rmsd = avg_rmsd + rmsd;
	}

	avg_corr = avg_corr / rows_metaphe;
	avg_mad = avg_mad / rows_metaphe;
	avg_rmsd = avg_rmsd / rows_metaphe;

	double median_corr = findMedian(all_corr);
	double min_corr = findMin(all_corr);

	double median_mad = findMedian(all_mad);
	double min_mad = findMin(all_mad);

	double median_rmsd = findMedian(all_rmsd);
	double min_rmsd = findMin(all_rmsd);

	cout << "Average:" << avg_corr << "," << avg_mad << "," << avg_rmsd << endl;
	cout << "Median:" << median_corr << "," << median_mad << "," << median_rmsd << endl;
	cout << "Min:" << min_corr << "," << min_mad << "," << min_rmsd << endl;

    return 0;
}
