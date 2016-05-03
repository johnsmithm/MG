#include "matrixread.h"
#include<algorithm>
// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ////////////////////////////////////////////////////   Class Member Functions      /////////////////////////////////////////////////////////////
// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Matrix::Matrix(int rows,int columns){
    rows_ = rows;
    columns_ = columns;
    arrayPtr_ = new double __attribute__((aligned(64))) [ rows_*columns_];
    //for (int i =0; i<rows_*columns_;i++)
    //    arrayPtr_[i] = 0.0;
}

bool Matrix::readFile(const std::string& name){
    int dummyR,dummyC;
    std::ifstream file(name,std::ios::in);
    if (file.is_open()) {
        file>>dummyR;
        file>>dummyC;

        for (int i = 0 ; i < dummyR ; ++i ){
                for ( int j = 0 ; j < dummyC ; ++j ){
                    file>>arrayPtr_[i*columns_+j];
                }
            }
        }
    file.close();
    return false;
}

bool Matrix::writeFile(const std::string& name,int a1, int b2){
	std::ofstream file(name,std::ios::out);
	if (file.is_open()) {
        file<<a1<<" ";
        file<<b2;
	file<<"\n";
        for (int i = 0 ; i < a1 ; ++i ){
                for ( int j = 0 ; j <b2 ; ++j ){
                    file<<arrayPtr_[i*columns_+j]<<"\t";
				
                }
                file << "\n";
            }
        }
    file.close();
    return false;
}
	



void Matrix::print(){

       for(int i = columns_ -1; i >= 0  ; --i ) {
           for(int j = 0; j < rows_ ; ++j){
               std::cout<<std::left<<" "<<arrayPtr_[j*columns_+i]<<"  ";
           }
           std::cout<<"\n";
        }
}


Matrix::~Matrix(){
    //std::cout<<"destructor"<<std::endl;
    delete [] arrayPtr_;
    arrayPtr_ = NULL;
}

// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ////////////////////////////////////////////////////  Non Class Member Functions   /////////////////////////////////////////////////////////////
// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool writeGnuplotFile(const std::string& name, Matrix &src, int xCells,int yCells){

    std::cout << "Writing solution file..." << std::endl;
    std::ofstream file(name,std::ios::out);
    double hx =  2.0/xCells;
    double hy = 1.0/yCells;
    if (file.is_open()) {
        file << "#" << "x" << "\t" <<  "y" << "\t" << "u(x,y)" << "\n";
        for (int yIndex=0; yIndex <yCells+1; yIndex++)
        {
            for (int xIndex=0; xIndex <xCells+1; xIndex++)
            {
                file << xIndex*hx << "\t" << yIndex*hy << "\t" << src(yIndex, xIndex) << "\n";
            }
        }
     }
    file.close();
    return false;
}

int readDimensions(const std::string& name,int choice){
    //choice 1 = number of rows
    //choice 2 = number of columns
    int rows,columns;
    std::ifstream file(name,std::ios::in);
    if (file.is_open()) {
        file>>rows;
        file>>columns;
    }
    file.close();
    if( choice == 1)
        return rows;
    else if(choice == 2)
        return columns;
    else
        std::cout<<"Wrong choice. Try again with valid choice for reading dimension in readdimension()"<<std::endl;
        return 0;

}




