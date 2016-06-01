#include<fstream>
#include <cassert>
#include<iostream>
#include <stdlib.h>

#include "MG.h"
#include "Timer.h"
#include <sys/time.h>


using namespace std;

//#define DEBUG
#define CLUSTER

#ifdef DEBUG
	#include "TEST.h"
#endif

int main(int argc, char *argv[]){
	(void) argc; //to suppress Warnings about unused argc
	assert(argc>0);
	int l = atoi(argv[1]);
	int n = 10;//atoi(argv[2]);
	//cerr<<l<<'\n';  
	//cerr<<l<<'\n'; 

	#ifdef DEBUG
		TEST d(l,n);
		//do some test
	
	#else

		MG solver(l,n); 
		#ifdef CLUSTER
			solver.initiate();
			solver.writeGnuFile("init.dat");
			std::cout<<"Your Alias: "<<"teamAlias"<<std::endl;
			struct timeval t0, t;
			gettimeofday(&t0, NULL);	
			
		#else
			siwir::Timer timer;
	    #endif

		solver.solve();
	    
		
	 	#ifdef CLUSTER
	 		gettimeofday(&t, NULL);
	 		std::cout << "Wall clock time of MG execution: " <<
			((int64_t)(t.tv_sec - t0.tv_sec) * (int64_t)1000000 +
			(int64_t)t.tv_usec - (int64_t)t0.tv_usec) * 1e-3
			<< " ms" << std::endl;
			solver.writeGnuFile("solution.dat");
		#else
			double time = timer.elapsed();
			cout<<"Total time:"<<time<<'\n';
	    #endif
	    


   	#endif


}
