#!/bin/bash
:<<END
#1: Convert the training files into input format suitable for gtm-decon
input_dir=<path_to_data>/pancreas/train_valid_test/
output_dir=<path_to_results/pancreas_dataset/

mkdir $output_dir/
mkdir $output_dir/1topic/
mkdir $output_dir/2topics/
mkdir $output_dir/3topics/
mkdir $output_dir/4topics/
mkdir $output_dir/5topics/

./singleCellInput $input_dir/counts_matrix_pp_train.tab $input_dir/cell_labels_oh_train.csv $output_dir/1topic/ 1
./singleCellInput $input_dir/counts_matrix_pp_train.tab $input_dir/cell_labels_oh_train.csv $output_dir/2topics/ 2
./singleCellInput $input_dir/counts_matrix_pp_train.tab $input_dir/cell_labels_oh_train.csv $output_dir/3topics/ 3
./singleCellInput $input_dir/counts_matrix_pp_train.tab $input_dir/cell_labels_oh_train.csv $output_dir/4topics/ 4
./singleCellInput $input_dir/counts_matrix_pp_train.tab $input_dir/cell_labels_oh_train.csv $output_dir/5topics/ 5

./singleCellInput_TestData $input_dir/counts_matrix_pp_test.tab $output_dir/1topic/testData.txt
cp $output_dir/1topic/testData.txt $output_dir/5topics/
END

##2: GTM-decon: Training 
path=pancreas_dataset/
scdata=$path/trainData.txt
scmeta=$path/meta.txt
scprior=$path/priorData.txt

K=14 ## K=number of topics; for modeling each celltype as 2topics, set K to be 28 for this set; for modeling each celltype as 3topics, set K to be 42 for this set; etc.

niter=5
:<<END
# nmar model
./gtm-decon --outputIntermediates -f $scdata -m $scmeta -trp $scprior -k $K -i $niter --inferenceMethod JCVB0 --maxcores 10

##3: Preparing bulk RNAseq data in the same gene order as training dataset
python3 python/prepare_TCGA_RNAseq.py --path_input ../../../data/TCGA/PAAD/gdac.broadinstitute.org_PAAD.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0/PAAD.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt --path_save ../../../data/TCGA/PAAD/gdac.broadinstitute.org_PAAD.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0/PAAD.csv
python3 python/prepare_bulkRNAseq_input.py --path_input ../../../data/TCGA/PAAD/gdac.broadinstitute.org_PAAD.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0/PAAD.csv --path_save ../../../data/TCGA/PAAD/gdac.broadinstitute.org_PAAD.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0/PAAD.pancreas.pp --preprocessed_genes pancreas_dataset/genes.txt

sed -i 's/\.0//g' ../../../data/TCGA/PAAD/gdac.broadinstitute.org_PAAD.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0/PAAD.pancreas.pp.tab
sed -i 's/\.0//g' ../../../data/TCGA/PAAD/gdac.broadinstitute.org_PAAD.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0/PAAD.pancreas.pp_normr.tab
sed -i 's/\.0//g' ../../../data/TCGA/PAAD/gdac.broadinstitute.org_PAAD.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0/PAAD.pancreas.pp_log1p_normr.tab

./bulkRNAseqInput ../../../data/TCGA/PAAD/gdac.broadinstitute.org_PAAD.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0/PAAD.pancreas.pp.tab $output_dir/1topic/TCGA_PAAD_bulkData.txt
cp $output_dir/1topic/TCGA_PAAD_bulkData.txt $output_dir/5topics/TCGA_PAAD_bulkData.txt
END
# Infer metagene mixing proportions for all RS samples

artificialbulkdata=$path/artificial_bulkData.txt
tcgabulkdata=$path/TCGA_PAAD_bulkData.txt
for file in `ls $path/trainData_JCVB0_nmar_K14_iter* | awk -F "_iter" '{print $2}' |  awk -F "_" '{print $1}' | sort -u`
do
	echo $file
	./gtm-decon -m $scmeta -n JCVB0 --newRSSamplesData $artificialbulkdata \
		                --trainedModelPrefix $path/trainData_JCVB0_nmar_K14_iter$file -k $K --inferNewSampleRSMetagene \
	        	                        --inferRSSampleParams_maxiter 100
	./gtm-decon -m $scmeta -n JCVB0 --newRSSamplesData $tcgabulkdata \
		                --trainedModelPrefix $path/trainData_JCVB0_nmar_K14_iter$file -k $K --inferNewSampleRSMetagene \
	        	                        --inferRSSampleParams_maxiter 100
done
exit

#3:  Normalize "metagene mixing proportions for number of topics > 1
for file in `ls $path/2topics/*metagene.csv`
do
       ofile=`echo $file | sed 's/metagene.csv/metagene_normalized.csv/g'`
       ./gtm-decon/Metagene_topics_avg $file 2 > $ofile
done

for file in `ls $path/3topics/*metagene.csv`
do
       ofile=`echo $file | sed 's/metagene.csv/metagene_normalized.csv/g'`
       ./gtm-decon/Metagene_topics_avg $file 3 > $ofile
done

for file in `ls $path/4topics/*metagene.csv`
do
       ofile=`echo $file | sed 's/metagene.csv/metagene_normalized.csv/g'`
       ./gtm-decon/Metagene_topics_avg $file 4 > $ofile
done

for file in `ls $path/5topics/*metagene.csv`
do
        ofile=`echo $file | sed 's/metagene.csv/metagene_normalized.csv/g'`
        ./gtm-decon/Metagene_topics_avg $file 5 > $ofile
done


