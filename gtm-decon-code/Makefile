gpp := g++-8 -std=c++17 -I/usr/local/boost-1.68.0/include -I/usr/local/include -L/usr/local/lib
all: gtm-decon
gtm-decon: main.o Decon.o JCVB0.o SampleRS.o parseData.o updateParams.o loglik.o
	$(gpp) -fopenmp -o gtm-decon -g main.o Decon.o JCVB0.o SampleRS.o parseData.o updateParams.o loglik.o -O2 -larmadillo

main.o: src/main.cpp src/Decon.h src/JCVB0.h
	$(gpp) -fopenmp -g -c -Wall src/main.cpp -O2 -larmadillo

Decon.o: src/Decon.cpp src/Decon.h
	$(gpp) -fopenmp -g -c -Wall src/Decon.cpp -O2 -larmadillo

JCVB0.o: src/JCVB0.cpp src/JCVB0.h src/SampleRS.h src/GeneParams.h 
	$(gpp) -fopenmp -g -c -Wall src/JCVB0.cpp -O2 -larmadillo
	
loglik.o: src/loglik.cpp src/JCVB0.h
	$(gpp) -fopenmp -g -c -Wall src/loglik.cpp -O2 -larmadillo
	
updateParams.o: src/updateParams.cpp src/JCVB0.h
	$(gpp) -fopenmp -g -c -Wall src/updateParams.cpp -O2 -larmadillo

SampleRS.o: src/SampleRS.cpp src/SampleRS.h
	$(gpp) -fopenmp -g -c -Wall src/SampleRS.cpp -O2 -larmadillo
	
parseData.o: src/parseData.cpp src/Decon.h src/GeneParams.h src/JCVB0.h src/SampleRS.h
	$(gpp) -fopenmp -g -c -Wall src/parseData.cpp -O2 -larmadillo

clean:
	rm -f *.o

clean2:
	rm -f examples/*CVB0*

