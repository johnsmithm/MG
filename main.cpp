#include <iostream>
#include <math.h>
#include <cstdlib>
#include "matrixread.h"

const double pi = M_PI; //3.1415926535897932;

void initBC(Matrix& src , const double dx, const int yMax, int xCells , int yCells){

    std::cout << "Initializing the boundary conditions" << std::endl;
    const int yIndex = yCells;

    const double sinh_y = sinh(pi * yMax );
//    std::cout << "sin h at y max " << sinh_y << std::endl;
    for ( int xIndex = 0; xIndex < xCells ; ++xIndex){
        // not iterating over last xIndex as we know that sin(2*pi*xMax) = 0.0 as xMax= 2(this is taken care in initialization to zero )
       src(yIndex,xIndex) = sin(pi * dx * xIndex) * sinh_y;
    }
}

int main(int argc, char *argv[]){
	(void) argc; //to suppress Warnings about unused argc
	assert(argc>1);
	const int l = atoi(argv[1]);
	//const int n = atoi(argv[2]);    //For use in multi grid implementation but not here at bc initialization
   

   //input parameters
   const int xMin = 0;
   const int xMax = 1;
   //const int yMin = 0;
   const int yMax = 1;

//Assigning the matrix sizes
//src(yIndex,xIndex) = 0;
   const int xCells = pow(2,l); 
   const int yCells = pow(2,l);
   const int n = yCells+1;
   const int m= xCells +1;
   //const int iter = atoi(argv[3]);

   // Declarations of user defined variables
   // src stores value of u at grid poitns and rhs stores value of f at gridpoints

   //siwir::Timer timer;
   //double time;
   Matrix src(n , m);
   Matrix rhs(n , m);
   //double residualNorm = 0;


//Initialising the Rhs matrix

   for (int yIndex=0; yIndex < yCells+1; yIndex++)
   {
       for (int xIndex=0; xIndex < xCells+1; xIndex++)
       {
           rhs(yIndex,xIndex) = 0;
	   src(yIndex,xIndex) = 0;
       }
   }


   //calucating constants used in calculations
   const double dx = double((xMax - xMin)) / xCells;
   //const double dy = double((yMax - yMin)) / yCells;
   //const double invX 	 = 1.0/(dx*dx);
   //const double invY 	 = 1.0/(dy*dy);
   //const double weightC  = 2.0*invX + 2.0*invY + (4.0 * pi * pi);
   //const double invC 	 = 1.0/weightC;


   initBC(src , dx , yMax, xCells ,yCells);
   src.writeFile(argv[2],n,m);
   //rhs.writeFile(argv[3],n,m);
   //initRHS(rhs , dx , dy , xCells ,yCells);

return 0;

}
