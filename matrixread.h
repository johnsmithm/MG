#ifndef MATRIXREAD_H
#define MATRIXREAD_H

#include<iostream>
#include<vector>
#include<assert.h>
#include <iomanip>
#include<fstream>
#include<sstream>
#include<vector>

class Matrix{
public:
    Matrix(int rows,int columns);
    inline double& operator()(int,int);
    void print();
    bool readFile(const std::string& name);
    bool writeFile(const std::string& name, int a1,int b2);
    ~Matrix();

private:
int rows_ ,columns_ ;
double* arrayPtr_;
};



inline double& Matrix::operator ()(int i, int j ){
    assert(i<rows_ && j<columns_);
    return  arrayPtr_[i*columns_+j];
}

int readDimensions(const std::string& name,int choice);
bool writeGnuplotFile(const std::string& name, Matrix &src, int xCells, int yCells);
#endif // MATMULT_H
