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
		
		

		double polar(double x, double y)
		{ 
			double r = sqrt((pow(x,2))+(pow(y,2)));

		    double theta = atan2(y,x);
		    if(y<0)
		    	theta+=pi*2;

			return sqrt(r)*sin(theta/2.);
		}

		double getNr(int i, int j,int lev,int nr){
			// black (0,0)
			if((i%2==0 && j%2==0)||(i%2==1&&j%2==1)){
				if(i%2==0)
					return gridsB[lev][(i/2)*nr + j/2] ;
				else
					return gridsB[lev][(i*nr)/2 + j/2 - nr/2] ;
			}else{
				if(i%2==0)
					return gridsR[lev][(i/2)*nr + j/2] ;
				else
					return  gridsR[lev][(i*nr)/2 + j/2 - nr/2] ;
			}
			return 0;
		}

		void setNr(int i,int j,int lev,int nr,int val){
			// black (0,0)
			if((i%2==0 && j%2==0)||(i%2==1&&j%2==1)){
				if(i%2==0)
					gridsB[lev][(i/2)*nr + j/2] = val;
				else
					gridsB[lev][(i*nr)/2 + j/2 - nr/2] = val;
			}else{
				if(i%2==0)
					gridsR[lev][(i/2)*nr + j/2] = val;
				else
					gridsR[lev][(i*nr)/2 + j/2 - nr/2] = val;
			}
		}

		void initiate(){
			gridsR = new double*[l];
			gridsB = new double*[l];
			f = new double*[l];
			re = new double*[l];
			for(int i=l;i;--i){// note we have l matrixes put one after another, no lost of memory!!!
				gridsB[l-i] = new double[((1<<i)+1)*((1<<i)+1)/2 + 1];
				gridsR[l-i] = new double[((1<<i)+1)*((1<<i)+1)/2 + 1];
				f[l-i] = new double[((1<<i)+1)*((1<<i)+1)];
				re[l-i] = new double[((1<<i)+1)*((1<<i)+1)];
			}
			double h = 2./(1<<l);
			double yTB = 1., xLR = -1;
			int bC,nr = (1<<l)+1;//backCells, halfRow, smallIndex;

			//first touch

			//black nodes
			for(int i=0;i<nr;++i){
				for(int j=(i%2?1:2);j<nr;j+=2){//not sure condition j<nr
					bC = (i*nr)/2 + j/2;
					if(i == 0 || j == 0 || i == nr-1 || j==nr-1){
						xLR = -1. + (double)j*h;
						yTB = 1. - (double)i*h;
						gridsB[0][bC] = polar(xLR,yTB) ;
					}else gridsB[0][bC] = 0;
				}
			}


			//red nodes
			for(int i=0;i<nr;++i){
				for(int j=(i%2?2:1);j<nr;j+=2){
					bC = (i*nr)/2 + j/2;
					if(i == 0 || j == 0 || i == nr-1 || j==nr-1){	
												
						xLR = -1. + (double)j*h;
						yTB = 1. - (double)i*h;
						gridsR[0][bC] = polar(xLR,yTB);
					}else gridsR[0][bC] = 0;
				}
			}
		}

		void smooth(double * aB, double * aR, int nr, double * f2){
			//(red-black) Gauss-Seidel for relaxation
			double h = 2./(nr-1);// precalculate ToDo
			double co = 4./(h*h); 
			int bC,hR=nr/2;//backCells, halfRow, smallIndex;

			//black nodes
			for(int i=1;i<nr-1;++i)
				for(int j=(i%2?1:2);j<nr-2;j+=2)
					if(!(i == (nr / 2) && j >= (nr / 2))){		
						bC = (i*nr)/2+j/2;
						aB[bC] = 0.25*(aR[bC-hR]+aR[bC+hR] + aR[bC]+aR[bC+((i&1)?1:-1)])
						+ f2[i*nr+j]/co;
					}


			//red nodes
			for(int i=1;i<nr-1;++i)
				for(int j=(i%2?2:1);j<nr-2;j+=2)
					if(!(i == (nr / 2) && j >= (nr / 2))){
						bC = (i*nr)/2 + j/2;
						aR[bC] = 0.25*(aB[bC-hR]+aB[bC+hR] + aB[bC]+aB[bC+((i&1)?-1:1)])
						+ f2[i*nr+j]/co;
					}
					//else cerr<<'-';
			//cerr<<"\n";

			//for(int i=nr/2;i<nr-1;++i)
				//a[(nr/2)*nr + i] = 0;

		}

		bool notBoundary(int x, int nr){
			//return true;
			//if(((x >= (((nr / 2))*nr + (nr/2)) )&& (x >= ((nr / 2 + 1)*nr - 1) ) )) {return false;}
			if(x<0 || x>=nr*nr){return false;}
			return true;
			if(x<nr)return false;
			if(x >= nr*(nr-1))return false;
			if(x%nr == 1) return false;
			if(x%nr == nr-1)return false;
			return true;
			//todo check if correct
		}

		void clasifyD(int lev, double rezidialCell, int i, int j){
			int nr1 = (1<<(l-lev-1)) +1;
			if(i%2==0 && j%2==0){//the point to be restricted
						f[lev+1][(i/2)*nr1+j/2] += rezidialCell/4.;
			}else{
				if(i%2==1 && j%2==1){//we add to the diagonal
					f[lev+1][((i+1)/2)*(nr1)+(j+1)/2] += rezidialCell/16.;//right top	
					f[lev+1][((i-1)/2)*(nr1)+(j+1)/2] += rezidialCell/16.;//right bottom
					f[lev+1][((i+1)/2)*(nr1)+(j-1)/2] += rezidialCell/16.;//left top
					f[lev+1][((i-1)/2)*(nr1)+(j-1)/2] += rezidialCell/16.;//right bottom

				}else if(i%2==1){//we add to top and bottom
					f[lev+1][((i+1)/2)*(nr1)+(j)/2] += rezidialCell/8.;// top							
					f[lev+1][((i-1)/2)*(nr1)+(j)/2] += rezidialCell/8.;// bottom

				}else{//we add to left and right		
					f[lev+1][((i)/2)*(nr1)+(j+1)/2] += rezidialCell/8.;//right 							
					f[lev+1][((i)/2)*(nr1)+(j-1)/2] += rezidialCell/8.;//right 
				}
			}
		}


		//restriction
		void downsampling(int lev){//from nr to nr/2, a is nr*nr

			double temp = 0;
			int nr = (1<<(l-lev)) +1 ; 
			double h = 2./(nr-1);            // To check if this gives the updated fine matrix after interpolation or the original matrix
            double st = -1./(h*h), co = 4./(h*h);
			int bC,hR=nr/2;//backCells, halfRow, smallIndex;
			double *f2 = f[lev], *aB = gridsB[lev], *aR = gridsR[lev];

			//black nodes
			for(int i=1;i<nr-1;++i)
				for(int j=(i%2?1:2);j<nr-2;j+=2)
					if(!(i == (nr / 2) && j >= (nr / 2))){		
						bC = (i*nr)/2+j/2;
						temp =f2[i*nr+j] - (aB[bC]*co + st*(aR[bC-hR]+aR[bC+hR] + aR[bC]+aR[bC+((i&1)?1:-1)]));
						clasifyD(lev,temp,i,j);
					}


			//red nodes
			for(int i=1;i<nr-1;++i)
				for(int j=(i%2?2:1);j<nr-2;j+=2)
					if(!(i == (nr / 2) && j >= (nr / 2))){	
						bC = (i*nr)/2 + j/2;
						temp = f2[i*nr+j] -(aR[bC]*co + st*(aB[bC-hR]+aB[bC+hR] + aB[bC]+aB[bC+((i&1)?-1:1)]));
						clasifyD(lev,temp,i,j);
					}

			int nr1 = (1<<(l-lev-1)) +1;
			//set boundaries to zero
				for(int j=0;j<nr1;++j){
					f[lev+1][(nr1-1)*nr1+j]=0;
					f[lev+1][(0)*nr1+j]=0;
					f[lev+1][(j)*nr1+0]=0;
					f[lev+1][(j)*nr1+nr1-1]=0;					
				}

				for(int j=nr1/2;j<nr1;++j)
					f[lev+1][(nr1/2)*nr1+j]=0;
		}

		void clasifyI(double * A, double * B, int lev, int i, int j){//??
			int nr1 = (1<<(l-lev-1)) +1, nr = (1<<(l-lev)) +1;
			int id = (i)*(nr/2)+j/2;
			if(i%2==0 && j%2==0){//the point to be restricted
				A[id] += B[(i/2)*(nr1/2)+j/4];
			}else{
				if(i%2==1 && j%2==1){//we add to the diagonal
					A[id] += B[(((i+1)/2)*(nr1/2))+(j+1)/4]/4.;//right top							
					A[id] += B[(((i-1)/2)*(nr1/2))+(j+1)/4]/4.;//right bottom
					A[id] += B[(((i+1)/2)*(nr1/2))+(j-1)/4]/4.;//left top
					A[id] += B[(((i-1)/2)*(nr1/2))+(j-1)/4]/4.;//right bottom

				}else if(i%2==1){//we add to top and bottom
					A[id] += B[(((i+1)/2)*(nr1/2))+(j)/4]/2.;// top							
					A[id] += B[(((i-1)/2)*(nr1/2))+(j)/4]/2.;// bottom

				}else{//we add to left and right		
					A[id] += B[(((i)/2)*(nr1/2))+(j+1)/4]/2.;//right 							
					A[id] += B[(((i)/2)*(nr1/2))+(j-1)/4]/2.;//right 
				}
			}
		}



		void interpolation(int lev){//from nr/2 to nr, (a) is (nr/2)*(nr/2)

			int nr = (1<<(l-lev)) +1 ; 
			double  *aB = gridsB[lev], *aR = gridsR[lev];

			//black nodes
			for(int i=1;i<nr-1;++i)
				for(int j=(i%2?1:2);j<nr-2;j+=2)
					if(!(i == (nr / 2) && j >= (nr / 2))){		
						//bC = (i*nr)/2+j/2;
						//temp =f2[i*nr+j] - (aB[bC]*co + st*(aR[bC-hR]+aR[bC+hR] + aR[bC]+aR[bC+((i&1)?1:-1)]));
						clasifyI(aB,gridsB[lev+1],lev,i,j);
					}


			//red nodes
			for(int i=1;i<nr-1;++i)
				for(int j=(i%2?2:1);j<nr-2;j+=2)
					if(!(i == (nr / 2) && j >= (nr / 2))){
						clasifyI(aR,gridsR[lev+1],lev,i,j);
					}


			//corect inner boundary
			for(int j=nr/2;j<nr;j+=2){
				gridsR[lev][((nr/2)*nr)/2+j/2]=0;
				gridsB[lev][((nr/2)*nr)/2+j/2]=0;
			}
		}

		void interpolation1(int lev){//from nr/2 to nr, (a) is (nr/2)*(nr/2)
			//bi-linear interpolation
			int nr = (1<<(l-lev))+1;
			int nr1 = (1<<(l-lev-1))+1;

			for(int i=1;i<nr-1;++i){
				for(int j=1;j<nr-1;++j){
					//construct each cell
					//amd so on more 9 times, or can use the stencil and one more for loop
					//grids[lev][indices2] += grids[lev+1][indices1] * somesclaler(1/2,1,1/4); 
					if(i%2==0 && j%2==0){//the point to be restricted
						//grids[lev][(i)*nr+j] += grids[lev+1][(i/2)*nr1+j/2];
						setNr(i,j,lev,nr,getNr(i/2,j/2,lev+1,nr1));
						//cout<<"the cell\n";
					}else{
						if(i%2==1 && j%2==1){//we add to the diagonal
							//cout<<"diagonal"<<"\n";
							//if(notBoundary(((i+1)/2)*(nr1)+(j+1)/2,nr1))
								//grids[lev][(i)*nr+j] += grids[lev+1][((i+1)/2)*(nr1)+(j+1)/2]/4.;//right top
								setNr(i,j,lev,nr,getNr((i+1)/2,(j+1)/2,lev+1,nr1)/4.);							
							//if(notBoundary(((i-1)/2)*(nr1)+(j+1)/2,nr1))
								//grids[lev][(i)*nr+j] += grids[lev+1][((i-1)/2)*(nr1)+(j+1)/2]/4.;//right bottom
								setNr(i,j,lev,nr,getNr((i-1)/2,(j+1)/2,lev+1,nr1)/4.);
							//if(notBoundary(((i+1)/2)*(nr1)+(j-1)/2,nr1))
								//grids[lev][(i)*nr+j] += grids[lev+1][((i+1)/2)*(nr1)+(j-1)/2]/4.;//left top
								setNr(i,j,lev,nr,getNr((i+1)/2,(j-1)/2,lev+1,nr1)/4.);
							//if(notBoundary(((i-1)/2)*(nr1)+(j-1)/2,nr1))
								//grids[lev][(i)*nr+j] += grids[lev+1][((i-1)/2)*(nr1)+(j-1)/2]/4.;//right bottom
								setNr(i,j,lev,nr,getNr((i-1)/2,(j-1)/2,lev+1,nr1)/4.);

						}else if(i%2==1){//we add to top and bottom
							//cout<<"tb\n";
							//if(notBoundary(((i+1)/2)*(nr1)+(j)/2,nr1))
								//grids[lev][(i)*nr+j] += grids[lev+1][((i+1)/2)*(nr1)+(j)/2]/2.;// top
								setNr(i,j,lev,nr,getNr((i+1)/2,(j)/2,lev+1,nr1)/2.);							
							//if(notBoundary(((i-1)/2)*(nr1)+(j)/2,nr1))
								//grids[lev][(i)*nr+j] += grids[lev+1][((i-1)/2)*(nr1)+(j)/2]/2.;// bottom
								setNr(i,j,lev,nr,getNr((i-1)/2,(j)/2,lev+1,nr1)/2.);

						}else{//we add to left and right		
							//cout<<"rl\n";
						//	if(notBoundary(((i)/2)*(nr1)+(j+1)/2,nr1))
								//grids[lev][(i)*nr+j] += grids[lev+1][((i)/2)*(nr1)+(j+1)/2]/2.;//right 	
								setNr(i,j,lev,nr,getNr((i)/2,(j+1)/2,lev+1,nr1)/2.);						
						//	if(notBoundary(((i)/2)*(nr1)+(j-1)/2,nr1))
								//grids[lev][(i)*nr+j] += grids[lev+1][((i)/2)*(nr1)+(j-1)/2]/2.;//right 
								setNr(i,j,lev,nr,getNr((i)/2,(j-1)/2,lev+1,nr1)/2.);
						}
					}
				}
			}
			
		}

		//level
		void recoursionMG(int lev){
			//for testing with GS method just ancoment the return and set 2 to 20 in the for loop
			for(int i=0;i<2;i++)
				smooth(gridsB[lev],gridsR[lev],(1<<(l-lev))+1,f[lev]);
			//return;
			if(lev + 1 == l)
				return;
			//debug
			//cerr<<((1<<(l-lev))+1)<<"x"<<((1<<(l-lev))+1)<<" solution-before downs\n";
			//test_print(grids[lev],((1<<(l-lev))+1));
			//cerr<<"residual\n";
		  downsampling(lev);
		  //return;
			//debug
			//cerr<<((1<<(l-lev-1))+1)<<"x"<<((1<<(l-lev-1))+1)<<" f2- after downs\n";
			//test_print(f[lev+1],((1<<(l-lev-1))+1));
			

			int nr = (1<<(l-lev-1))+1;
			for(int i=0;i<nr;++i)
				for(int j=0;j<nr;j+=2)
					{
						gridsR[lev+1][i*(nr/2)+j/2]=0;
						gridsB[lev+1][i*(nr/2)+j/2]=0;
					}
			recoursionMG(lev+1);
			


			//debug
			//cerr<<((1<<(l-lev))+1)<<"x"<<((1<<(l-lev))+1)<<" before\n";
			//test_print(grids[lev],((1<<(l-lev))+1));

				//cout<<"before residualNorm:"<<residualNorm(lev)<<"\n";

			interpolation1(lev);

				//cout<<"after residualNorm:"<<residualNorm(lev)<<"\n";

			//debug
			//cerr<<((1<<(l-lev))+1)<<"x"<<((1<<(l-lev))+1)<<" after\n";
			//test_print(grids[lev],((1<<(l-lev))+1));
	
			for(int i=0;i<2;i++)
				smooth(gridsB[lev],gridsR[lev],(1<<(l-lev))+1,f[lev]);

		}

		double Lnorm(){
		  
			double error = 0, temp = 0;
			int nr = (1<<l) +1 ;
			
			double h = 2./(nr-1);
			double yTB = 1., xLR = -1;
			int bC;//backCells, halfRow, smallIndex;

			//black nodes
			for(int i=1;i<nr-1;++i){
				for(int j=(i%2?1:2);j<nr-2;j+=2){
					if(!(i == (nr / 2) && j >= (nr / 2))){		
						bC = (i*nr)/2+j/2;
						xLR = -1. + (double)j*h;
						yTB = 1. - (double)i*h;
						temp = gridsB[0][bC] - polar(xLR,yTB) ;

						error+=temp*temp;
					}
				}
			}


			//red nodes
			for(int i=1;i<nr-1;++i){
				for(int j=(i%2?2:1);j<nr-2;j+=2){
					if(!(i == (nr / 2) && j >= (nr / 2))){	
						bC = (i*nr)/2 + j/2;						
						xLR = -1. + (double)j*h;
						yTB = 1. - (double)i*h;
						temp = gridsR[0][bC] - polar(xLR,yTB);
						error+=temp*temp;
					}
				}
			}


			//cerr<<nr<<" "<<x1<<" "<<y1<<" "<<tM<<"-\n";
			double norm = sqrt(error/(nr*nr));
			return norm;
		}

		double residualNorm(int lev){
			double temp = 0, residuum=0;
			int nr = (1<<(l-lev)) +1 ;
			int dom = (1<<(l-lev))-1;    // Size of the domain where we calculate the norm>> internal grid points
			//int lev=0; 
			double h = 2./(nr-1);            // To check if this gives the updated fine matrix after interpolation or the original matrix
            double st = -1./(h*h), co = 4./(h*h);
			int bC,hR=nr/2;//backCells, halfRow, smallIndex;
			double *f2 = f[lev], *aB = gridsB[lev], *aR = gridsR[lev];

			//black nodes
			for(int i=1;i<nr-1;++i)
				for(int j=(i%2?1:2);j<nr-2;j+=2)
					if(!(i == (nr / 2) && j >= (nr / 2))){		
						bC = (i*nr)/2+j/2;
						temp =f2[i*nr+j] - (aB[bC]*co + st*(aR[bC-hR]+aR[bC+hR] + aR[bC]+aR[bC+((i&1)?1:-1)]));
						re[lev][i*nr+j] = temp;
						residuum+=temp*temp;
					}


			//red nodes
			for(int i=1;i<nr-1;++i)
				for(int j=(i%2?2:1);j<nr-2;j+=2)
					if(!(i == (nr / 2) && j >= (nr / 2))){	
						bC = (i*nr)/2 + j/2;
						temp = f2[i*nr+j] -(aR[bC]*co + st*(aB[bC-hR]+aB[bC+hR] + aB[bC]+aB[bC+((i&1)?-1:1)]));
						residuum+=temp*temp;
						re[lev][i*nr+j] = temp;
					}

 			double norm= sqrt(residuum/(dom*dom));
			return norm;

		}

	public:	
		void solve(){
			initiate();

			//debug
			//cerr<<((1<<l)+1)<<"x"<<((1<<l)+1)<<"\n";
			//test_print(grids[0],((1<<l)+1));

		

			//perform MG n times
			for(int i=0;i<n;++i){
				recoursionMG(0);
				cout<<"Step:"<<i<<"\n";
			  	//cout<<"Lnorm:"<<Lnorm()<<"\n";
					cout<<"residualNorm:"<<residualNorm(0)<<"\n";
					double resvectorC= residualNorm(0);	
					cout<<"Error Norm:"<<Lnorm()<<"\n";	
					if(i>0){
						double convergencerate = resvectorC/resvector;						
						cout<< "convergence rate at step :  "<<" "<<i<<"="<< convergencerate<<"\n" ;
					}
					resvector = resvectorC;
				}
				
					

			//debug
			//cerr<<((1<<l)+1)<<"x"<<((1<<l)+1)<<"\n";
			//test_print(grids[0],((1<<l)+1));
			writeGnuFile("solution.txt");
			writeGnuFile1("realsol.txt");
			//cerr<<((1<<(l-1))+1)<<"x"<<((1<<(l-1))+1)<<"\n";
			//test_print(f[1],((1<<(l-1))+1));
		}

		bool writeGnuFile(const std::string& name){

		    std::cout << "Solution file being written " << std::endl;
		    std::ofstream file(name,std::ios::out);
		    //double hy = 1./l;
		    int nr = (1<<l)+1;
		    double h =  2./(nr-1);
		    double yTB = 1., xLR = -1;
		    int bC;//backCells, halfRow, smallIndex;

			//black nodes
			for(int i=0;i<nr;++i){
				for(int j=(i%2?1:2);j<nr;j+=2){	
						bC = (i*nr)/2+j/2;
						xLR = -1. + (double)j*h;
						yTB = 1. - (double)i*h;
						file<< xLR << "\t" << yTB << "\t" <<gridsB[0][bC]<<"\n";
					
				}
			}


			//red nodes
			for(int i=0;i<nr;++i){
				for(int j=(i%2?2:1);j<nr;j+=2){	
						bC = (i*nr)/2 + j/2;						
						xLR = -1. + (double)j*h;
						yTB = 1. - (double)i*h;
						file<< xLR << "\t" << yTB << "\t" <<gridsB[0][bC]<<"\n";
					
				}
			}

		    file.close();
		    return false;
		}

		bool writeGnuFile1(const std::string& name){

		    std::cout << "Real Solution file being written " << std::endl;
		    std::ofstream file(name,std::ios::out);
		    //double hy = 1./l;
		    int nr = (1<<l)+1;
		    double h =  2./(nr-1);
		    double yTB = 1., xLR = -1;
		    if (file.is_open()) {
		        file << "#" << "x" << "\t" <<  "y" << "\t" << "u" << "\n";
				for(int i=0; i<nr;i+=1){
					xLR = -1.;
					for(int j=0;j<nr;j+=1){
						cerr<<"";
						file<< xLR << "\t" << yTB << "\t" <<polar(xLR,yTB)<<"\n";
								
								xLR+=h;
					}
					yTB -= h;
				}
			}

		file.close();
		return false;
		}



		void test_print(double a[], int nr){
			for(int i=0;i<nr;++i){
				for(int j=0;j<nr;++j)cerr<<a[i*nr+j]<<" ";
				cerr<<'\n';
			}
		}


	private:
		int l,n,N;//levels, # v-cycles, # grid points including boundary


		double **gridsR,**gridsB;// Red-Black
		

		double ** f, **re; // right hand side - Vector od vector(matrix)

		double resvector;  // Vector to store the residual norms


		double pi = 3.14159265359;
};
