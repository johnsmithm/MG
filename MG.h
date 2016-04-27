#include<math.h>
#include<iostream>
#include<fstream>

using namespace std;

class MG{
	public:
		MG(int l_, int n_):l(l_),n(n_),N((1<<l)+1){
		//do something
		}

		~MG(){}

	private:
		void initiate(){
			grids = new double[N*2*N];
			f = new double[N * (N/2) + 1];
			//initiate grid NXN
			//initiate f first N entries
		}

		void smooth(double * a, int nr, double * f2){
			//some code using GS method
			*a = 0 + nr + *f2;
		}

		void downsampling(double * a, int nr){
			*a = nr;
		}

		void interpolation(double * a, int nr){
			*a = nr;
		}

		//a - pointer to top-left corner of the grid
		//nr - size of a matrix, note that the width is 2N when propragate through
		//f - ponter to the left hand side the left hand side of size nr
		void recoursionMG(double * a, int nr, double * f2){

			for(int i=0;i<2;i++)
				smooth(a,nr,f2);

			downsampling(a,nr);

			//calculate f

			if(nr == 1) {
				*a = -(*a);
			}else{
				//set to zero nr/2 times nr/2 entries on (a+nr) matrix
				recoursionMG( (a+nr), (nr/2), (f2+nr));
			}

			interpolation(a,nr);
	
			for(int i=0;i<2;i++)
				smooth(a,nr,f2);
		}
	public:	
		void solve(){
			initiate();
			//perform MG n times
			for(int i=0;i<n;++i)
				recoursionMG(grids,N, f);	
		}

		void print_gnuplot(){
		//print NXN grid 
		}
	private:
		int l,n,N;//levels, # v-cycles, # grid points including boundary

		double *grids;// pointer to the N times 2N matrix that will be storing the grids
		// Each grid left-top corner will be: *grids, *(grids+N), *(grids + N/2) and so on
		// Note that the lenght of the matrix id 2N 

		double * f; // right hand side -vector of vectors N
		//each vector for a size(the begining): *v , *(v+N), *(v+N+N/2) ans so on 
};
