#include<iostream>
#include<vector>
#include "matrix.h"
using namespace std;
int main(){
	vector<vector<double>> v = {{2,3,4},{5,6,7},{8,9,1}};
	vector<vector<double>> v2 = {{1,2},{3,4},{5,6}};	
	Matrix<double> mat(v);
	cout<<mat.det()<<endl;
	Matrix<double> mat2(v2);
	mat *= mat2;	
	mat.print();
	return 0;
}
