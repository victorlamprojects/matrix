#ifndef MATRIX_H_
#define MATRIX_H_
#include<vector>
#include<string>
using namespace std;
template<typename T>
class Matrix{
	public:
		Matrix();	
		Matrix(int m,int n,int defaultVal=0);
		Matrix(vector<vector<T>>&);
		Matrix(const Matrix&);
		Matrix(Matrix&&);	
		Matrix<T>& operator=(const Matrix<T>&);
		Matrix<T>& operator=(Matrix&&);	
		Matrix<T>& operator+=(const Matrix<T>&);
		Matrix<T>& operator-=(const Matrix<T>&);
		Matrix<T>& operator*=(const Matrix<T>&);	
		Matrix<T> operator+(const Matrix<T>&);
		Matrix<T> operator-(const Matrix<T>&);
		Matrix<T> operator*(const Matrix<T>&);
		T det();	
		int LUPDecompose(T**,int, int*);
		Matrix<T> multiply(Matrix<T>&&,Matrix<T>&&);//only for 2 square matrices
		Matrix<T> pseudo_multiply(const Matrix<T>&,const Matrix<T>&,int,int,int,int); 	
		Matrix<T> top_left(const Matrix<T>&);	//only for 2n*2m matrix
		Matrix<T> top_right(const Matrix<T>&);	//...
		Matrix<T> bottom_left(const Matrix<T>&);//...
		Matrix<T> bottom_right(const Matrix<T>&);//...	
		void printErrorMessage(string msg) const;
		void print() const;
		~Matrix();
	private:
		int m, n;
		T** matrix;
};

#include "matrix.cpp"
#endif
