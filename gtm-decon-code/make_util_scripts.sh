g++ ../util_scripts/singleCellInput_MCLR.cpp -o singleCellInput_MCLR -larmadillo
g++ ../util_scripts/singleCellInput.cpp -o singleCellInput -larmadillo
g++ ../util_scripts/singleCellInput_TestData.cpp -o singleCellInput_TestData -larmadillo
g++ ../util_scripts/Accuracy.cpp -o Accuracy -larmadillo
g++ ../util_scripts/Decon_metrics.cpp -o Decon_metrics -larmadillo
g++ ../util_scripts/bulkRNAseqInput.cpp -o bulkRNAseqInput -larmadillo
g++ ../util_scripts/Metaphe_topics_avg.cpp -o Metaphe_topics_avg -larmadillo
g++ ../util_scripts/Scaling_phi.cpp -o Scaling_phi -larmadillo
g++ ../util_scripts/singleCellInput_BulkData.cpp -o singleCellInput_BulkData -larmadillo
g++ ../util_scripts/singleCellInput_DE.cpp -o singleCellInput_DE -larmadillo

chmod u+x singleCellInput
chmod u+x singleCellInput_TestData
chmod u+x Accuracy
chmod u+x Decon_metrics

