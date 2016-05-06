#include<math.h>
#include<iostream>
#include<fstream>

using namespace std;

class MG{
	public:
		MG(int l_, int n_):l(l_),n(n_),N((1<<l)+1){
		//do something
			//cerr<<l<<'\n';  
		}

		~MG(){}

	private:
		
		void initiateExtra(){
			
			//other boundary condition
		}

		void initiate(){
			grids = new double*[l];
			f = new double*[l];
			for(int i=l;i;--i){// note we have l matrixes put one after another, no lost of memory!!!
				grids[l-i] = new double[((1<<i)+1)*((1<<i)+1)];
				f[l-i] = new double[((1<<i)+1)*((1<<i)+1)];				
			}
			double h = 1./(1<<l);
			for(int i=0;i<N;++i){
				grids[0][i*N+N-1] = sin(pi)*sinh(i*pi*h);//(i,1) - column
				grids[0][N*(N-1)+i] = sin(pi*i*h)*sinh(pi);// (1,1) - row
			}
		}

		void smooth(double * a, int nr, double * f2){
			//(red-black) Gauss-Seidel for relaxation
			//make two for loops for red and black nodes, it is enough for the task
			
			double TB = 1./((nr-1)*(nr-1)), MIJ = 2*TB, LR = TB; // here we have the stencil

			//black nodes
			for(int i=1;i<nr-1;i=i+2)
				for(int j=(i%2?1:2);j<nr-1;j+=2)
					a[i*nr+j] = TB*(a[(i-1)*nr+j]+a[(i+1)*nr]+j) + LR*(a[i+nr+j+1]+a[i*nr+j-1]) + MIJ*f2[i*nr+j];


			//red nodes
			for(int i=1;i<nr-1;i=i+2)
				for(int j=(i%2?2:1);j<nr-1;j+=2)
					a[i*nr+j] = TB*(a[(i-1)*nr+j]+a[(i+1)*nr]+j) + LR*(a[i+nr+j+1]+a[i*nr+j-1]) + MIJ*f2[i*nr+j];

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
		void recoursionMG(int lev){

			for(int i=0;i<1;i++)
				smooth(grids[lev],(1<<(l-lev))+1,f[lev]);

			return;
			downsampling(grids[lev],(1<<(l-lev))+1);

			//calculate f

			if(lev == l) {//use a solver
				//*a = -(*a);
			}else{
				//set to zero nr/2+1 times nr/2+1 entries on (a+nr) matrix
				// nr/2+1 is because nr is odd, and inner grid begin with 1
				recoursionMG(l+1);//Todo choose corect the indices!!!
			}

			interpolation(grids[lev],(1<<(l-lev))+1);
	
			for(int i=0;i<2;i++)
				smooth(grids[lev],(1<<(l-lev))+1,f[lev]);

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

			cerr<<((1<<l)+1)<<"x"<<((1<<l)+1)<<"\n";
			test_print(grids[0],((1<<l)+1));

			//perform MG n times
			for(int i=0;i<1;++i){
				recoursionMG(0);
				cout<<"Step:"<<i<<"\n";
				cout<<"Lnorm:"<<Lnorm()<<"\n";
				cout<<"residualNorm:"<<residualNorm()<<"\n";
			}	

			cerr<<((1<<l)+1)<<"x"<<((1<<l)+1)<<"\n";
			test_print(grids[0],((1<<l)+1));
		}

		void print_gnuplot(){
		//print NXN grid 
		}

		void test_print(double a[], int nr){
			for(int i=0;i<nr;++i){
				for(int j=0;j<nr;++j)cerr<<a[i*nr+j]<<" ";
				cerr<<'\n';
			}
		}


	private:
		int l,n,N;//levels, # v-cycles, # grid points including boundary


		double **grids;// pointer to the N times 2N matrix that will be storing the grids
		// Each grid left-top corner will be: *grids, *(grids+N), *(grids + N/2) and so on
		// Note that the lenght of the matrix id 2N 

		double ** f; // right hand side -vector of vectors N
		//each vector for a size(the begining): *v , *(v+N), *(v+N+N/2) ans so on 

		double pi = 3.14;
};
