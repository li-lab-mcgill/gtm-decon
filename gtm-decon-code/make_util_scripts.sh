g++ util_scripts/singleCellInput_MCLR.cpp -o singleCellInput_MCLR -larmadillo
g++ util_scripts/singleCellInput.cpp -o singleCellInput -larmadillo
g++ util_scripts/singleCellInput_TestData.cpp -o singleCellInput_TestData -larmadillo
g++ util_scripts/Accuracy.cpp -o Accuracy -larmadillo
g++ util_scripts/Accuracy_BulkData.cpp -o Accuracy_BulkData -larmadillo
g++ util_scripts/Deconvolution_metrics.cpp -o Deconvolution_metrics -larmadillo
g++ util_scripts/bulkRNAseqInput.cpp -o bulkRNAseqInput -larmadillo
g++ util_scripts/singleCellInput.cpp -o singleCellInput -larmadillo
g++ util_scripts/Accuracy_3topics.cpp -o Accuracy_3topics -larmadillo
g++ util_scripts/Metaphe_2topics.cpp -o Metaphe_2topics -larmadillo
g++ util_scripts/Metaphe_3topics.cpp -o Metaphe_3topics -larmadillo
g++ util_scripts/Metaphe_4topics.cpp -o Metaphe_4topics -larmadillo
g++ util_scripts/Metaphe_5topics.cpp -o Metaphe_5topics -larmadillo
g++ util_scripts/Metaphe_10topics.cpp -o Metaphe_10topics -larmadillo
g++ util_scripts/Metaphe_topics_avg.cpp -o Metaphe_topics_avg -larmadillo
g++ util_scripts/Metaphe_topics_avg_denovo.cpp -o Metaphe_topics_avg_denovo -larmadillo
g++ util_scripts/Phi_2topics.cpp -o Phi_2topics -larmadillo
g++ util_scripts/Phi_3topics.cpp -o Phi_3topics -larmadillo
g++ util_scripts/Phi_4topics.cpp -o Phi_4topics -larmadillo
g++ util_scripts/Phi_5topics.cpp -o Phi_5topics -larmadillo
g++ util_scripts/Phi_10topics.cpp -o Phi_10topics -larmadillo
g++ util_scripts/Perplexity.cpp -o Perplexity -larmadillo
g++ util_scripts/Perplexity_log.cpp -o Perplexity_log -larmadillo
g++ util_scripts/Perplexity_log_phi.cpp -o Perplexity_log_phi -larmadillo
g++ util_scripts/Accuracy_10topics.cpp -o Accuracy_10topics -larmadillo
exit


g++-8 -std=c++17 -I/usr/local/boost-1.68.0/include -I/usr/local/include -L/usr/local/lib -g -c -Wall util_scripts/singleCellInput.cpp -O2 -larmadillo -o singleCellInput

g++-8 -std=c++17 -I/usr/local/boost-1.68.0/include -I/usr/local/include -L/usr/local/lib -g -c -Wall util_scripts/singleCellInput_TestData.cpp -O2 -larmadillo -o singleCellInput_TestData

g++-8 -std=c++17 -I/usr/local/boost-1.68.0/include -I/usr/local/include -L/usr/local/lib  -g -c -Wall util_scripts/Accuracy.cpp -O2 -o larmadilo Accuracy 

g++-8 -std=c++17 -I/usr/local/boost-1.68.0/include -I/usr/local/include -L/usr/local/lib  -g -c -Wall util_scripts/Accuracy_BulkData.cpp -O2 -larmadillo -o Accuracy_BulkData 

g++-8 -std=c++17 -I/usr/local/boost-1.68.0/include -I/usr/local/include -L/usr/local/lib -g -c -Wall util_scripts/Deconvolution_metrics.cpp -O2 -larmadillo -o Deconvolution_metrics

chmod u+x singleCellInput
chmod u+x singleCellInput_TestData
chmod u+x Accuracy
chmod u+x Accuracy_BulkData
chmod u+x Deconvolution_metrics

