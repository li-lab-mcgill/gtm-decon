#Step 1: Generate input files suitable for nested guided topic modeling

mkdir $output_dir/DEgenes/
./gtm-decon/singleCellInput_DE $input_dir/counts_matrix_train.tab $input_dir/cell_labels_oh_train.csv $output_dir/DEgenes/ 1 5 
./gtm-decon/singleCellInput_TestData $input_dir/counts_matrix_test.tab $output_dir/DEgenes/testData.txt

#Step 2: Generating combined phi matrix for initialization
#There are 7 celltypes for breast tissue - lets get the 5 characterised ones
#mportant: get the genes matching to the one in bulkRS data - and in the same order
awk -F "," '{print $3","$4","$5","$6","$7}' $output_dir/BN/pp/trainData_JCVB0_nmar_K7_iter4_phi.csv > one
paste -d , $output_dir/BN/pp/genes.txt one > two
python3 ./python/prepare_input_phi.py --path_input two --preprocessed_genes $output_dir/DEgenes/genes.txt --path_save three
#Rearrange the ids, drop gene names etc. to mimic phi format
awk -F "," '{print $2","$1}' three.csv > one 
awk -F "," '{print $4","$5","$6","$7","$8}' three.csv > two
paste -d , one two two > $output_dir/DEgenes/input_phi.csv

#For phi_normalized
awk -F "," '{print $3","$4","$5","$6","$7}' $output_dir/BN/pp/trainData_JCVB0_nmar_K7_iter4_phi_normalized.csv > one
paste -d , $output_dir/BN/pp/genes.txt one > two
python3 ./python/prepare_input_phi.py --path_input two --preprocessed_genes $output_dir/DEgenes/genes.txt --path_save three
awk -F "," '{print $2","$1}' three.csv > one 
awk -F "," '{print $4","$5","$6","$7","$8}' three.csv > two
paste -d , one two two > $output_dir/DEgenes/input_phi_normalized.csv

sed -i 's/\.0,/,/g' $output_dir/DEgenes/input_phi.csv
sed -i 's/\.0,/,/g' $output_dir/DEgenes/input_phi_normalized.csv

sed -i 's/\.0$//' $output_dir/DEgenes/input_phi.csv
sed -i 's/\.0$//' $output_dir/DEgenes/input_phi_normalized.csv

path=../../results_gtm-decon/TCGA/BRCA/DEgenes/
scdata=$path/trainData.txt
scmeta=$path/meta.txt
scprior=$path/priorData.txt
preset_path=$path/input
testdata=$path/testData.txt

K=10

niter=100

# nmar model
./gtm-decon --outputIntermediates -f $scdata -m $scmeta -trp $scprior -k $K -i $niter --inferenceMethod JCVB0 --maxcores 10 --prefixDE $preset_path

#Inference
for file in `ls $path/trainData_JCVB0_nmar_K10_iter* | awk -F "_iter" '{print $2}' |  awk -F "_" '{print $1}' | sort -u`
do
	  echo $file
       ./gtm-decon -m $scmeta -n JCVB0 --newRSSamplesData $scdata \
                             --trainedModelPrefix $path/trainData_JCVB0_nmar_K10_iter$file -k $K --inferNewSampleRSMetagene \
	                                           --inferRSSampleParams_maxiter 100
       ./gtm-decon -m $scmeta -n JCVB0 --newRSSamplesData $testdata \
                             --trainedModelPrefix $path/trainData_JCVB0_nmar_K10_iter$file -k $K --inferNewSampleRSMetagene \
                                                   --inferRSSampleParams_maxiter 100
done
