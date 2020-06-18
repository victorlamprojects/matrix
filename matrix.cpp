#ifndef MATRIX_TPP
#define MATRIX_TPP
#include<vector>
#include<iostream>
#include<iomanip>
#include<string>
using namespace std;

template <typename T>
Matrix<T>::Matrix():m(0), n(0), matrix(nullptr){}

template <typename T>
Matrix<T>::Matrix(int m, int n, int defaultVal){
	this->m = m;
	this->n = n;
	this->matrix = new T* [m];
	for(int i=0;i<m;i++){
		this->matrix[i] = new T[n];	
		for(int j=0;j<n;j++){
			this->matrix[i][j] = defaultVal;	
		}
	}
}

template<typename T>
Matrix<T>::Matrix(vector<vector<T>>& v){
	if(v.empty())
		return;
	this->m = v.size();
	this->matrix = new T* [this->m];
	for(int i=0;i<this->m;i++){
		this->n = v[i].size();
		this->matrix[i] = new T [this->n];	
		for(int j=0;j<this->n;j++){
			this->matrix[i][j] = v[i][j];
		}
	}
}

template<typename T>
Matrix<T>::Matrix(const Matrix& mat){
	this->n = 0;
	this->m = 0;
	*this = mat;
}
template<typename T>
Matrix<T>::Matrix(Matrix&& mat){
	this->n = mat.n;
	this->m = mat.m;
	this->matrix = mat.matrix;
	mat.n = 0;
	mat.m = 0;
	mat.matrix = nullptr;
}

template<typename T>
Matrix<T>::~Matrix(){
	if(this->matrix != nullptr){
		for(int i=0;i<this->m;i++){
			delete[] this->matrix[i];
			this->matrix[i] = nullptr;
		}
		delete[] this->matrix;
		this->matrix = nullptr;
	}
}

template<typename T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T>& mat){
	if(this == &mat)
		return *this;		
	T** temp = new T* [mat.m];
	for(int i=0;i<mat.m;i++){
		temp[i] = new T [mat.n];	
		for(int j=0;j<mat.n;j++)
			temp[i][j] = mat.matrix[i][j];	
	}
	if(this->matrix != nullptr && this->n != 0){
		for(int i=0;i<this->m;i++){
			delete[] this->matrix[i];
			this->matrix[i] = nullptr;
		}
		delete[] this->matrix;
		this->matrix = nullptr;
	}
	this->matrix = temp;
	this->m = mat.m;
	this->n = mat.n;	
	return *this;
}
template<typename T>
Matrix<T>& Matrix<T>::operator=(Matrix&& mat){
	if(this == &mat)
		return *this;		
	if(this->matrix != nullptr && this->n != 0){
		for(int i=0;i<this->m;i++){
			delete[] this->matrix[i];
			this->matrix[i] = nullptr;
		}
		delete[] this->matrix;
		this->matrix = nullptr;
	}
	this->matrix = mat.matrix;
	this->m = mat.m;
	this->n = mat.n;	
	mat.matrix = nullptr;	
	mat.n = 0;
	mat.m = 0;	
	return *this;
}

template<typename T>
Matrix<T>& Matrix<T>::operator+=(const Matrix<T>& mat){
	if(this->m == mat.m && this->n == mat.n)
		for(int i=0;i<this->m;i++)
			for(int j=0;j<this->n;j++)
				this->matrix[i][j] += mat.matrix[i][j];
	else
		this->printErrorMessage("Matrix size does not match.");
	return *this;
}

template<typename T>
Matrix<T>& Matrix<T>::operator-=(const Matrix<T>& mat){
	if(this->m == mat.m && this->n == mat.n)
		for(int i=0;i<this->m;i++)
			for(int j=0;j<this->n;j++)
				this->matrix[i][j] -= mat.matrix[i][j];
	else
		this->printErrorMessage("Matrix size does not match.");
	return *this;
}

template<typename T>
Matrix<T>& Matrix<T>::operator*=(const Matrix<T>& mat){
	if(this->m == 0 || this->n == 0 || mat.m == 0 || mat.n == 0)
		this->printErrorMessage("Matrix is empty!");
	else if(this->n != mat.m)
		this->printErrorMessage("Matrix size does not match");
	else
		*this = forward<Matrix<T>>((*this) * mat);	
	return *this;
}

template<typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T>& mat){
	if(this->m != mat.m){
		this->printErrorMessage("Matrices' row number does not match!");
		return Matrix<T>();	
	}
	if(this->n != mat.n){
		this->printErrorMessage("Matrices' column number does not match!");
		return Matrix<T>();
	}
	Matrix<T> n_mat(this->m, this->n);
	
	for(int i=0;i<this->m;i++)
		for(int j=0;j<this->n;j++)
			n_mat.matrix[i][j] = this->matrix[i][j] + mat.matrix[i][j];	
	return n_mat;	
}

template<typename T>
Matrix<T> Matrix<T>::operator-(const Matrix<T>& mat){
        if(this->m != mat.m){           
                this->printErrorMessage("Matrices' row number does not match!");
		return Matrix<T>();	
	}
        if(this->n != mat.n){
                this->printErrorMessage("Matrices' column number does not match!");
		return Matrix<T>();
        }
	Matrix<T> n_mat(this->m, this->n);    
        for(int i=0;i<this->m;i++)
                for(int j=0;j<this->n;j++)
                        n_mat.matrix[i][j] = this->matrix[i][j] - mat.matrix[i][j];
        return n_mat;
}
template<typename T>
Matrix<T> Matrix<T>::operator*(const Matrix<T>& mat){
	if(this->n != mat.m){
		this->printErrorMessage("Matrices' size does not match!");
		return Matrix<T>();
	}
	int size = max(max(mat.n, mat.m), max(this->n, this->m));
	
	int n = 2;
	while(size > n)
		n*=2;
	Matrix<T> m1(n,n,0);
	Matrix<T> m2(n,n,0);
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			if(i<this->m && j<this->n)
				m1.matrix[i][j] = this->matrix[i][j];
			else
				m1.matrix[i][j] = 0;
			if(i<mat.m && j<mat.n)
				m2.matrix[i][j] = mat.matrix[i][j];
			else
				m2.matrix[i][j] = 0;
		}
	}	
	
	Matrix<T> result = multiply(forward<Matrix<T>>(m1),forward<Matrix<T>>(m2));
	Matrix<T> product(this->m, mat.n, 0);
	for(int i=0;i<this->m;i++)
		for(int j=0;j<mat.n;j++)
			product.matrix[i][j] = result.matrix[i][j];
	return product;
}

template<typename T>
T Matrix<T>::det(){
	if(this->n != this->m){
		this->printErrorMessage("This is not a square matrix");
		return -9999;
	}
	if(this->n == 1)
		return this->matrix[0][0];	
	if(this->n == 2)
		return this->matrix[0][0]*this->matrix[1][1] - this->matrix[0][1]*this->matrix[1][0];
	T det;
	Matrix<T> X(*this);
	int *p=new int [this->n+1];
	LUPDecompose(X.matrix, X.n, p);
	det = X.matrix[0][0];
	for(int i=1;i<this->n;i++)
		det *= X.matrix[i][i];
	if((p[this->n]-this->n)%2 == 0)
		return det;
	return -1*det;
}	

template<typename T>
int Matrix<T>::LUPDecompose(T** A, int N, int *P){
	int imax;
	T maxA, absA;
	for(int i=0;i<=N;i++)
		P[i] = i;
	for(int i=0;i<N;i++){
		maxA = 0;
		imax = i;
		for(int k=i;k<N;k++){
			absA = A[k][i] < 0? -1*A[k][i]: A[k][i];
			if(absA > maxA){
				maxA = absA;
				imax = k;	
			}	
		}
		//if(maxA < Tol) return 0;//matrix is degenerate
		if(imax != i){
			//pivoting P
			int t = P[i];
			P[i] = P[imax];
			P[imax] = t;	
			//pivoting rows of A
			T* ptr = A[i];
			A[i] = A[imax];
			A[imax] = ptr;
			//counting pivots starting from N(for determinant)
			P[N]++;
		}
		for(int j = i+1; j<N;j++){
			A[j][i] /= A[i][i];
			
			for(int k=i+1;k<N;k++){
				A[j][k] -= A[j][i] * A[i][k];
			}	
		}
	}
	return 1;
}

template<typename T>
Matrix<T> Matrix<T>::multiply(Matrix<T>&& m1, Matrix<T>&& m2){
	int size = m1.n;
	Matrix<T> ans(size,size,0);
	if(size == 2){
		Matrix<T> m(2,2,0);
		m.matrix[0][0] = m1.matrix[0][0]*m2.matrix[0][0]+m1.matrix[0][1]*m2.matrix[1][0];
		m.matrix[0][1] = m1.matrix[0][0]*m2.matrix[0][1]+m1.matrix[0][1]*m2.matrix[1][1];
		m.matrix[1][0] = m1.matrix[1][0]*m2.matrix[0][0]+m1.matrix[1][1]*m2.matrix[1][0];
		m.matrix[1][1] = m1.matrix[1][0]*m2.matrix[0][1]+m1.matrix[1][1]*m2.matrix[1][1];
		return m;
	}
	Matrix<T> A11 = top_left(m1);
	Matrix<T> A12 = top_right(m1);
	Matrix<T> A21 = bottom_left(m1);
	Matrix<T> A22 = bottom_right(m1);
	Matrix<T> B11 = top_left(m2);
	Matrix<T> B12 = top_right(m2);
	Matrix<T> B21 = bottom_left(m2);
	Matrix<T> B22 = bottom_right(m2);	
	Matrix<T> M1 = multiply(A11 + A22 ,B11 + B22);
	Matrix<T> M2 = multiply(A21 + A22,forward<Matrix<T>>(B11));
	Matrix<T> M3 = multiply(forward<Matrix<T>>(A11), B12-B22);
	Matrix<T> M4 = multiply(forward<Matrix<T>>(A22), B21-B11);
	Matrix<T> M5 = multiply(A11 + A12,forward<Matrix<T>>(B22));
	Matrix<T> M6 = multiply(A21 - A11, B11 + B12);
	Matrix<T> M7 = multiply(A12 - A22, B21 + B22);
	Matrix<T> C11 = M1 + M4 - M5 + M7;
	Matrix<T> C12 = M3 + M5;
	Matrix<T> C21 = M2 + M4;
	Matrix<T> C22 = M1 - M2 + M3 + M6;
	for(int i=0;i<size/2;i++){
		for(int j=0;j<size/2;j++){
			ans.matrix[i][j] = C11.matrix[i][j];
			ans.matrix[i][j+size/2] = C12.matrix[i][j];
			ans.matrix[i+size/2][j] = C21.matrix[i][j];
			ans.matrix[i+size/2][j+size/2] = C22.matrix[i][j];
		}
	}
	return ans;
}
template<typename T>
Matrix<T> Matrix<T>::pseudo_multiply(const Matrix<T>& m1, const Matrix<T>& m2, int r1, int r2, int c1, int c2){
	int size = r2-r1;	
	if(r2-r1 == 1){
		Matrix<T> m(2,2,0);
		m.matrix[r1][c1] = m1.matrix[r1][c1]*m2.matrix[r1][c1]+m1.matrix[r1][c2]*m2.matrix[r2][c1];
		m.matrix[r1][c2] = m1.matrix[r1][c1]*m2.matrix[r1][c2]+m1.matrix[r1][c2]*m2.matrix[r2][c2];
		m.matrix[r2][c1] = m1.matrix[r2][c1]*m2.matrix[r1][c1]+m1.matrix[r2][c2]*m2.matrix[r2][c1];
		m.matrix[r2][c2] = m1.matrix[r2][c1]*m2.matrix[r1][c2]+m1.matrix[r2][c2]*m2.matrix[r2][c2];
		return m;
	}
}

template<typename T>
Matrix<T> Matrix<T>::top_left(const Matrix<T>& m){
	int size = m.n/2;
	Matrix<T> v(size,size,0);
	for(int i=0;i<size;i++)
		for(int j=0;j<size;j++)
			v.matrix[i][j] = m.matrix[i][j];
	return v;	
}
template<typename T>
Matrix<T> Matrix<T>::top_right(const Matrix<T>& m){
	int size = m.n/2;
	Matrix<T> v(size,size,0);
	for(int i=0;i<size;i++)
		for(int j=0;j<size;j++)
			v.matrix[i][j] = m.matrix[i][j+size];	
	return v;
}
template<typename T>
Matrix<T> Matrix<T>::bottom_left(const Matrix<T>& m){
	int size = m.n/2;
	Matrix<T> v(size,size,0);
	for(int i=0;i<size;i++)
		for(int j=0;j<size;j++)
			v.matrix[i][j] = m.matrix[i+size][j];	
	return v;
}
template<typename T>
Matrix<T> Matrix<T>::bottom_right(const Matrix<T>& m){
	int size = m.n/2;
	Matrix<T> v(size,size,0);
	for(int i=0;i<size;i++)
		for(int j=0;j<size;j++)
			v.matrix[i][j] = m.matrix[i+size][j+size];	
	return v;
}

template<typename T>
void Matrix<T>::printErrorMessage(string msg) const{
	cout<<"ERROR: "<<msg<<endl;
}
template<typename T>
void Matrix<T>::print() const{
	cout<<endl;	
	for(int i=0;i<this->m;i++){
		for(int j=0;j<this->n;j++){
			cout<<setw(5)<<this->matrix[i][j];	
		}
		cout<<endl;
	}
	cout<<endl;
}

#endif
