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
				grids[0][i*N+N-1] = sin(pi)*sinh(i*pi*h);//(i,1) -last column
				grids[0][N*(N-1)+i] = sin(pi*i*h)*sinh(pi);// (1,i) -first row
			}
		}

		void smooth(double * a, int nr, double * f2){
			//(red-black) Gauss-Seidel for relaxation
			//make two for loops for red and black nodes, it is enough for the task
			
			double TB = 1./((nr-1)*(nr-1)), MIJ = 2*TB, LR = TB; // here we have the stencil
			//todo calculate the stencil corect
			
			//black nodes
			for(int i=1;i<nr-1;++i)
				for(int j=(i%2?1:2);j<nr-1;j+=2)
					a[i*nr+j] = TB*(a[(i-1)*nr+j]+a[(i+1)*nr]+j) + LR*(a[i+nr+j+1]+a[i*nr+j-1]) + MIJ*f2[i*nr+j];


			//red nodes
			for(int i=1;i<nr-1;++i)
				for(int j=(i%2?2:1);j<nr-1;j+=2)
					a[i*nr+j] = TB*(a[(i-1)*nr+j]+a[(i+1)*nr]+j) + LR*(a[i+nr+j+1]+a[i*nr+j-1]) + MIJ*f2[i*nr+j];

		}

		//restriction
		void downsampling(int lev){//from nr to nr/2, a is nr*nr
			//full weighting for restriction
			//for
			int nr = (1<<(l-lev))+1;
			double TB = 1./((nr-1)*(nr-1)), MIJ = 2*TB, LR = TB; // here we have the stencil

			//black nodes
			for(int i=1;i<nr-1;++i)
				for(int j=1;j<nr-1;j+=1){
					//calculate rezidual
					double rezidialCell =  - (TB*(grids[lev][(i-1)*nr+j]+grids[lev][(i+1)*nr]+j)
					+ LR*(grids[lev][i+nr+j+1]+grids[lev][i*nr+j-1]) + MIJ*grids[lev][i*nr+j]) + f[lev][i*nr+j];
					
					//from small matrix from rezidual grids[lev+1]
					//grid[lev+1][indices2] += residialCell * somesclaler(1/2,1,1/4); todo
				}
		}

		void interpolation(int lev){//from nr/2 to nr, (a) is (nr/2)*(nr/2)
			//bi-linear interpolation
			int nr = (1<<(l-lev))+1;
			for(int i=1;i<nr-1;++i){
				for(int j=1;j<nr-1;++j){
					//construct each cell
					//amd so on more 9 times, or can use the stencil and one more for loop
					//grid[lev][indices2] += grid[lev+1][indices1] * somesclaler(1/2,1,1/4); todo
				}
			}
			
		}

		//level
		void recoursionMG(int lev){

			for(int i=0;i<1;i++)
				smooth(grids[lev],(1<<(l-lev))+1,f[lev]);

			return;
			downsampling(lev);


			if(lev == l) {//use a solver
				//todo calculate the solution
			}else{
				recoursionMG(l+1);
			}

			interpolation(lev);
	
			for(int i=0;i<2;i++)
				smooth(grids[lev],(1<<(l-lev))+1,f[lev]);

		}

		double Lnorm(){
			double norm = 0.;
			//ToDo calculate the rezidual of the grid 
			return norm;
		}

		double residualNorm(){
			double norm = 0.;
			//ToDo calculate the rezidual of the grid 
			return norm;
		}

	public:	
		void solve(){
			initiate();
			
			//debug
			cerr<<((1<<l)+1)<<"x"<<((1<<l)+1)<<"\n";
			test_print(grids[0],((1<<l)+1));

			//perform MG n times
			for(int i=0;i<1;++i){
				recoursionMG(0);
				cout<<"Step:"<<i<<"\n";
				cout<<"Lnorm:"<<Lnorm()<<"\n";
				cout<<"residualNorm:"<<residualNorm()<<"\n";
			}	

			//debug
			cerr<<((1<<l)+1)<<"x"<<((1<<l)+1)<<"\n";
			test_print(grids[0],((1<<l)+1));
		}

		void print_gnuplot(){
		//print NXN grid todo
		}

		void test_print(double a[], int nr){
			for(int i=0;i<nr;++i){
				for(int j=0;j<nr;++j)cerr<<a[i*nr+j]<<" ";
				cerr<<'\n';
			}
		}


	private:
		int l,n,N;//levels, # v-cycles, # grid points including boundary


		double **grids;// Vector od vector(matrix)
		//Each grid  will be: (1^l)+1 + (1^l)+1 size and so on
		

		double ** f; // right hand side - Vector od vector(matrix)
		//Each grid  will be: (1^l)+1 + (1^l)+1 size and so on

		double pi = 3.14;
};
