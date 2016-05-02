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
		void initiateExtra(){
			grids = new double[N*N*2];// we will have a lot of non usefull entries
			f = new double[N*N*2];
			//other boundary condition
		}

		void initiate(){
			grids = new double[N*2*N];// we will have a lot of non usefull entries
			f = new double[N * (N/2) + 1];
			//initiate grid NXN
			//initiate f first N entries
		}

		void smooth(double * a, int nr, double * f2){
			//(red-black) Gauss-Seidel for relaxation
			//make two for loops for red and black nodes, it is enough for the task
			*a = 0 + nr + *f2;
		}

		void downsampling(double * a, int nr){//from nr to nr/2, a is nr*nr
			//full weighting for restriction
			//for
			*a = nr;
		}

		void interpolation(double * a, int nr){//from nr/2 to nr, (a) is (nr/2)*(nr/2)
			//bi-linear interpolation
			for(int i=1;i<nr/2;++i){
				for(int j=1;j<nr/2;++j){
					//Todo choose corect the indices!!!
					//spread to the neighbours
					if(i+1 < nr/2-1)//check if we are not at the boundary
						a[(i*2)*N+j*2] = (1/2)*a[(i)*N+j+1+nr];//nr is here offset
					//amd so on more 9 times, or can use the stencil and one more for loop
				}
			}
			
		}

		//a - pointer to top-left corner of the grid
		//nr - size of a matrix, note that the width is 2N when propragate through
		//f - ponter to the left hand side the left hand side of size nr
		// note 2N is the lenght of the total grid, so the second line of the grid (will be 2*N + a)
		void recoursionMG(double * a, int nr, double * f2){

			for(int i=0;i<2;i++)
				smooth(a,nr,f2);

			downsampling(a,nr);

			//calculate f

			if(nr == 1) {//use a solver
				*a = -(*a);
			}else{
				//set to zero nr/2+1 times nr/2+1 entries on (a+nr) matrix
				// nr/2+1 is because nr is odd, and inner grid begin with 1
				recoursionMG( (a+nr), (nr/2), (f2+nr*nr));//Todo choose corect the indices!!!
			}

			interpolation(a+nr/2+1,nr);
	
			for(int i=0;i<2;i++)
				smooth(a,nr,f2);


		}

		double Lnorm(){
			double norm = 0.;

			return norm;
		}

		double residualNorm(){
			double norm = 0.;
			//ToDo calculate the rezidual of the grid from the squere with the top-left corner at 'a'
			return norm;
		}

	public:	
		void solve(){
			initiate();
			//perform MG n times
			for(int i=0;i<n;++i){
				recoursionMG(grids,N, f);
				cout<<"Step:"<<i<<"\n";
				cout<<"Lnorm:"<<Lnorm()<<"\n";
				cout<<"residualNorm:"<<residualNorm()<<"\n";
			}	
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
