#include<math.h>
#include<iostream>
#include<fstream>
#include<vector>

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
			grids = new double*[l];
			f = new double*[l];
			for(int i=l;i;--i){// note we have l matrixes put one after another, no lost of memory!!!
				grids[l-i] = new double[((1<<i)+1)*((1<<i)+1)];
				f[l-i] = new double[((1<<i)+1)*((1<<i)+1)];				
			}

			for(int i=0;i<N;++i)
				for(int j=0;j<N;++j)f[0][i*N+j]=2;

			double h = 1./(1<<l);
			for(int i=0;i<N;++i){
				grids[0][i*N] = i*h*(i*h-1);//(i,1) -first column
				grids[0][i] = -1;// (1,i) -first row
				grids[0][i*N+N-1] = i*h*(i*h-1);//(i,1) -last column
				grids[0][N*(N-1)+i] = -1;// (1,i) -last row
				//to do : corect boundaries condition for rows
			}

		}

		void initiate(){
			grids = new double*[l];
			f = new double*[l];
			resvector= new double[n];
			for(int i=l;i;--i){// note we have l matrixes put one after another, no lost of memory!!!
				grids[l-i] = new double[((1<<i)+1)*((1<<i)+1)];
				f[l-i] = new double[((1<<i)+1)*((1<<i)+1)];				
			}
			double h = 2./(1<<l);
			double yTB = -1., xLR = -1;
			for(int i=0;i<N;++i){
				//r = ~/ (x^2+y^2) teta = tan^-1 (y/x)
				// g(x,y) = r^1/2 * sin(teta/2)
				grids[0][i*N] = sqrt(sqrt(1+yTB*yTB))*sin((atan(yTB/(-1.))+(yTB > 0 ? 2*pi : 2*pi))/2.); //sqrt((sqrt(1+yTB*yTB)-1)/2.);//(i,N) -first column
				grids[0][i*N+N-1] = sqrt(sqrt(1+yTB*yTB))*sin((atan(yTB/(1.))+(yTB > 0 ? 0 : 4*pi))/2.); //sqrt((sqrt(1+yTB*yTB)-1)/2.);//(i,1) -last column
				grids[0][N*(N-1)+i] = sqrt(sqrt(1+xLR*xLR))*sin((atan((-1.)/xLR)+(xLR > 0 ? 4*pi : 2*pi))/2.); // sqrt((sqrt(1+xLR*xLR)-(xLR > 0 ? xLR: -xLR))/2.);// (1,i) -first row
				grids[0][i] = sqrt(sqrt(1+xLR*xLR))*sin((atan(1/xLR)+(xLR > 0 ? 0 : 2*pi))/2.); //sqrt((sqrt(1+xLR*xLR)-(xLR > 0 ? xLR: -xLR))/2.);// (N,i) -last row
				yTB+=h;
				xLR+=h;

				//I	Use the calculator value
				//II	Add 180° to the calculator value
				//III	Add 180° to the calculator value
				//IV	Add 360° to the calculator value
			}
			//cerr << "x="<<xLR<<" y="<<yTB<<"\n";
		}

		void smooth(double * a, int nr, double * f2){
			//(red-black) Gauss-Seidel for relaxation
			//make two for loops for red and black nodes, it is enough for the task
			double h = 2./(nr-1);
			double st = 1./(h*h), co = 4./(h*h); // here we have the stencil
			//todo calculate the stencil corect
			
			//exclude from smoothing inner boundary
			//int yIndexInnerBoundary = nr / 2;
			//int xIndexInnerBoundary = nr / 2 ;

			//black nodes
			for(int i=1;i<nr-1;++i)
				for(int j=(i%2?1:2);j<nr-1;j+=2)
					if(!(i == (nr / 2) && j >= (nr / 2)))
						a[i*nr+j] = (st*(a[(i-1)*nr+j]+a[(i+1)*nr+j]) + st*(a[i*nr+j+1]+a[i*nr+j-1]) + f2[i*nr+j])/co;


			//red nodes
			for(int i=1;i<nr-1;++i)
				for(int j=(i%2?2:1);j<nr-1;j+=2)
					if(!(i == (nr / 2) && j >= (nr / 2)))
						a[i*nr+j] = (st*(a[(i-1)*nr+j]+a[(i+1)*nr+j] + a[i*nr+j+1]+a[i*nr+j-1]) + f2[i*nr+j])/co;

		}

		bool notBoundary(int x, int nr){
			if(((x >= (((nr / 2)-1)*nr + (nr/2)) )&& (x >= ((nr / 2)*nr - 1) ) )) {return false;}
			if(x<0 || x>=nr*nr){return false;}
			return true;
			if(x<nr)return false;
			if(x >= nr*(nr-1))return false;
			if(x%nr == 1) return false;
			if(x%nr == nr-1)return false;
			return true;
			//todo check if correct
		}

		//restriction
		void downsampling(int lev){//from nr to nr/2, a is nr*nr
			//full weighting for restriction
			//for
			int nr = (1<<(l-lev))+1;
			int nr1 = (1<<(l-lev-1))+1;
			double h = 1./(nr-1);          
            double st = -1./(h*h), co = 4./(h*h);
			
			//  set f[lev+1] to 0
			for(int i=0;i<nr1;++i)
				for(int j=0;j<nr1;++j)
					f[lev+1][i*nr1+j]=0;

			for(int i=1;i<nr-1;++i){
				for(int j=1;j<nr-1;j+=1){
					if((i == (nr / 2) && j >= (nr / 2)))
						continue;
					//calculate rezidual
					double rezidialCell = f[lev][i*nr+j]  - (st*(grids[lev][(i-1)*nr+j]+grids[lev][(i+1)*nr+j])
					+ st*(grids[lev][i*nr+j+1]+grids[lev][i*nr+j-1]) + co*grids[lev][i*nr+j]) ;
					//rezidialCell = sqrt(rezidialCell*rezidialCell); 
					//cout<<i<<" j="<<j<<" r="<<rezidialCell<<" \n";
					//cout<<rezidialCell<<" ";
					//from small matrix from rezidual grids[lev+1]
					//f[lev+1][indices2] += residialCell * somesclaler(1/2,1,1/4);
					if(i%2==0 && j%2==0){//the point to be restricted
						f[lev+1][(i/2)*nr1+j/2] += rezidialCell/4.;
						//cout<<"the cell\n";
					}else{
						if(i%2==1 && j%2==1){//we add to the diagonal
							//cout<<"diagonal"<<"\n";
							if(notBoundary(((i+1)/2)*(nr1)+(j+1)/2,nr1))
								f[lev+1][((i+1)/2)*(nr1)+(j+1)/2] += rezidialCell/16.;//right top							
							if(notBoundary(((i-1)/2)*(nr1)+(j+1)/2,nr1))
								f[lev+1][((i-1)/2)*(nr1)+(j+1)/2] += rezidialCell/16.;//right bottom
							if(notBoundary(((i+1)/2)*(nr1)+(j-1)/2,nr1))
								f[lev+1][((i+1)/2)*(nr1)+(j-1)/2] += rezidialCell/16.;//left top
							if(notBoundary(((i-1)/2)*(nr1)+(j-1)/2,nr1))
								f[lev+1][((i-1)/2)*(nr1)+(j-1)/2] += rezidialCell/16.;//right bottom

						}else if(i%2==1){//we add to top and bottom
							//cout<<"tb\n";
							if(notBoundary(((i+1)/2)*(nr1)+(j)/2,nr1))
								f[lev+1][((i+1)/2)*(nr1)+(j)/2] += rezidialCell/8.;// top							
							if(notBoundary(((i-1)/2)*(nr1)+(j)/2,nr1))
								f[lev+1][((i-1)/2)*(nr1)+(j)/2] += rezidialCell/8.;// bottom

						}else{//we add to left and right		
							//cout<<"rl\n";
							if(notBoundary(((i)/2)*(nr1)+(j+1)/2,nr1))
								f[lev+1][((i)/2)*(nr1)+(j+1)/2] += rezidialCell/8.;//right 							
							if(notBoundary(((i)/2)*(nr1)+(j-1)/2,nr1))
								f[lev+1][((i)/2)*(nr1)+(j-1)/2] += rezidialCell/8.;//right 
						}
					}
				}
				//cout<<'\n';
			}
			//set boundaries to zeor
				for(int j=0;j<nr1;++j){
					f[lev+1][(nr1-1)*nr1+j]=0;
					f[lev+1][(0)*nr1+j]=0;
					f[lev+1][(j)*nr1+0]=0;
					f[lev+1][(j)*nr1+nr1-1]=0;					
				}
		}

		void interpolation(int lev){//from nr/2 to nr, (a) is (nr/2)*(nr/2)
			//bi-linear interpolation
			int nr = (1<<(l-lev))+1;
			int nr1 = (1<<(l-lev-1))+1;

			for(int i=1;i<nr-1;++i){
				for(int j=1;j<nr-1;++j){
					if((i == (nr / 2) && j >= (nr / 2)))
						continue;
					//construct each cell
					//amd so on more 9 times, or can use the stencil and one more for loop
					//grids[lev][indices2] += grids[lev+1][indices1] * somesclaler(1/2,1,1/4); 
					if(i%2==0 && j%2==0){//the point to be restricted
						grids[lev][(i)*nr+j] += grids[lev+1][(i/2)*nr1+j/2];
						//cout<<"the cell\n";
					}else{
						if(i%2==1 && j%2==1){//we add to the diagonal
							//cout<<"diagonal"<<"\n";
							if(notBoundary(((i+1)/2)*(nr1)+(j+1)/2,nr1))
								grids[lev][(i)*nr+j] += grids[lev+1][((i+1)/2)*(nr1)+(j+1)/2]/4.;//right top							
							if(notBoundary(((i-1)/2)*(nr1)+(j+1)/2,nr1))
								grids[lev][(i)*nr+j] += grids[lev+1][((i-1)/2)*(nr1)+(j+1)/2]/4.;//right bottom
							if(notBoundary(((i+1)/2)*(nr1)+(j-1)/2,nr1))
								grids[lev][(i)*nr+j] += grids[lev+1][((i+1)/2)*(nr1)+(j-1)/2]/4.;//left top
							if(notBoundary(((i-1)/2)*(nr1)+(j-1)/2,nr1))
								grids[lev][(i)*nr+j] += grids[lev+1][((i-1)/2)*(nr1)+(j-1)/2]/4.;//right bottom

						}else if(i%2==1){//we add to top and bottom
							//cout<<"tb\n";
							if(notBoundary(((i+1)/2)*(nr1)+(j)/2,nr1))
								grids[lev][(i)*nr+j] += grids[lev+1][((i+1)/2)*(nr1)+(j)/2]/2.;// top							
							if(notBoundary(((i-1)/2)*(nr1)+(j)/2,nr1))
								grids[lev][(i)*nr+j] += grids[lev+1][((i-1)/2)*(nr1)+(j)/2]/2.;// bottom

						}else{//we add to left and right		
							//cout<<"rl\n";
							if(notBoundary(((i)/2)*(nr1)+(j+1)/2,nr1))
								grids[lev][(i)*nr+j] += grids[lev+1][((i)/2)*(nr1)+(j+1)/2]/2.;//right 							
							if(notBoundary(((i)/2)*(nr1)+(j-1)/2,nr1))
								grids[lev][(i)*nr+j] += grids[lev+1][((i)/2)*(nr1)+(j-1)/2]/2.;//right 
						}
					}
				}
			}
			
		}

		//level
		void recoursionMG(int lev){
			//for testing with GS method just ancoment the return and set 2 to 20 in the for loop
			for(int i=0;i<2;i++)
				smooth(grids[lev],(1<<(l-lev))+1,f[lev]);
			//return;
			//debug
			//cerr<<((1<<(l-lev))+1)<<"x"<<((1<<(l-lev))+1)<<" solution-before downs\n";
			//test_print(grids[lev],((1<<(l-lev))+1));
			//cerr<<"residual\n";
			downsampling(lev);
			//debug
			//cerr<<((1<<(l-lev-1))+1)<<"x"<<((1<<(l-lev-1))+1)<<" f2- after downs\n";
			//test_print(f[lev+1],((1<<(l-lev-1))+1));
			

			if(lev+2 == l) {//use a solver
				//todo calculate the solution
				double h = 1./2;
				double co = 4./(h*h), st = -1./(h*h);
				int nr = 3;

				for(int i =0;i<3;++i)
					for(int j=0;j<3;++j)grids[lev+1][i*3+1]=0;

				grids[lev+1][1*3+1] = (st*(grids[lev+1][(1-1)*nr+1]+grids[lev+1][(1+1)*nr+1])
				 + st*(grids[lev+1][1*nr+1+1]+grids[lev+1][1*nr+1-1]) + f[lev+1][1*nr+1])/co;//for testing
				//cout<<"sol1:"<<grids[lev+1][1*3+1]<<"\n";
				//grids[lev+1][1*3+1] = 0.1;
			}else{
				int nr = (1<<(l-lev-1))+1;
				for(int i=0;i<nr;++i)
					for(int j=0;j<nr;++j)grids[lev+1][i*nr+j]=0;
				recoursionMG(lev+1);
			}


			//debug
			//cerr<<((1<<(l-lev))+1)<<"x"<<((1<<(l-lev))+1)<<" before\n";
			//test_print(grids[lev],((1<<(l-lev))+1));

				//cout<<"before residualNorm:"<<residualNorm(lev)<<"\n";

			interpolation(lev);

				//cout<<"after residualNorm:"<<residualNorm(lev)<<"\n";

			//debug
			//cerr<<((1<<(l-lev))+1)<<"x"<<((1<<(l-lev))+1)<<" after\n";
			//test_print(grids[lev],((1<<(l-lev))+1));
	
			for(int i=0;i<2;i++)
				smooth(grids[lev],(1<<(l-lev))+1,f[lev]);

		}

		double Lnorm(){
		  
			double error = 0, temp = 0;
			int nr = (1<<l) +1 ;
			
			double h = 2./(nr-1);
			for ( int i=1; i<nr-1;i+=1){
				for (int j=1;j<nr-1;j+=1)
					{
						temp =  grids[0][(i*nr)+j] - sin(pi*j*h)*sinh(pi*i*h);
						error += temp*temp;
					}
						    }
			double norm = sqrt(error);
			return norm;
		}

		double residualNorm(int lev){
			double temp = 0, residuum=0;
			int nr = (1<<(l-lev)) +1 ;
			int dom = (1<<(l-lev))-1;    // Size of the domain where we calculate the norm>> internal grid points
			//int lev=0; 
			double h = 2./(nr-1);            // To check if this gives the updated fine matrix after interpolation or the original matrix
            double st = -1./(h*h), co = 4./(h*h);
			for ( int i=1;i< nr-1 ; i+=1)
				{
				for (int j=1;j<nr-1;j+=1)
				     {
					temp = f[lev][i*nr+j] - ( st*(grids[lev][(i-1)*nr+j]+grids[lev][(i+1)*nr+j])
					+ st*(grids[lev][i*nr+j+1]+grids[lev][i*nr+j-1]) + co*grids[lev][i*nr+j]);

					residuum+=temp*temp;

				      }
				}
 					double norm= sqrt(residuum/(dom*dom));
					return norm;

		}

	public:	
		void solve(){
			initiate();
			//initiateExtra();
			
			//debug
			//cerr<<((1<<l)+1)<<"x"<<((1<<l)+1)<<"\n";
			//test_print(grids[0],((1<<l)+1));

		

			//perform MG n times
			for(int i=0;i<n;++i){
				recoursionMG(0);
				cout<<"Step:"<<i<<"\n";
			  	//cout<<"Lnorm:"<<Lnorm()<<"\n";
					cout<<"residualNorm:"<<residualNorm(0)<<"\n";
					resvector[i]= residualNorm(0);	
					//cout<<"Error Norm:"<<Lnorm()<<"\n";	
					if(i>0){
						double convergencerate = resvector[i]/resvector[i-1];
						cout<< "convergence rate at step :  "<<" "<<i<<"="<< convergencerate<<"\n" ;
					}
				}
				
					

			//debug
			cerr<<((1<<l)+1)<<"x"<<((1<<l)+1)<<"\n";
			test_print(grids[0],((1<<l)+1));
			//writeGnuFile("solution.txt");
			//writeGnuFile1("realsol.txt");
			//cerr<<((1<<(l-1))+1)<<"x"<<((1<<(l-1))+1)<<"\n";
			//test_print(f[1],((1<<(l-1))+1));
		}

bool writeGnuFile(const std::string& name){

    std::cout << "Solution file being written " << std::endl;
    std::ofstream file(name,std::ios::out);
    //double hy = 1./l;
    int nr = (1<<l)+1;
    double h =  1./(nr-1);
    if (file.is_open()) {
        file << "#" << "x" << "\t" <<  "y" << "\t" << "u" << "\n";
        for (int i=1; i<nr-1; i+=1)
        {
            for (int j=1; j <nr-1; j+=1)
            {
               file<< j*h << "\t" << i*h << "\t" <<grids[0][i*nr+j]<<"\n";
            }
        }
     }

    file.close();
    return false;
}

bool writeGnuFile1(const std::string& name){

    std::cout << "Solution file being written " << std::endl;
    std::ofstream file(name,std::ios::out);
    //double hy = 1./l;
    int nr = (1<<l)+1;
    double h =  1./(nr-1);
    if (file.is_open()) {
        file << "#" << "x" << "\t" <<  "y" << "\t" << "u" << "\n";
		for(int i=1; i<nr-1;i+=1){
			for(int j=1;j<nr-1;j+=1){
				file<< j*h << "\t" << i*h << "\t" <<sin(pi*j*h)*sinh(pi*i*h)<<"\n";
}}}

file.close();
return false;
}

//writeGnuFile1("realsol.txt");


		/*void print_gnuplot(){
		//print NXN grid todo
		}*/

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

		double *resvector;  // Vector to store the residual norms


		double pi = 3.14159265359;
};
