#include<fstream>
#include <cassert>
#include<iostream>
#include <stdlib.h>

#include "MG.h"
#include "Timer.h"


using namespace std;

//#define DEBUG

#ifdef DEBUG
	#include "TEST.h"
#endif

int main(int argc, char *argv[]){
	(void) argc; //to suppress Warnings about unused argc
	assert(argc>1);
	int l = atoi(argv[1]);
	int n = atoi(argv[2]);
	//cerr<<l<<'\n';  
	//cerr<<l<<'\n'; 

	#ifdef DEBUG
		TEST d(l,n);
		//do some test
	
	#else

		MG solver(l,n);

		siwir::Timer timer;
	    
		solver.solve();
	    
		double time = timer.elapsed();
	


		//solver.print_gnuplot();
		cout<<"Total time:"<<time<<'\n';

   	#endif


}
