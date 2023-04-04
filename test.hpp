#ifndef  INDEX_HPP
#define  INDEX_HPP
#include <vector>

using namespace std;

double normal(double mu , double sigma, double x);

double uniform(double a, double b);

double sinfun1(double q);

double sinfunn(double * Q);

int C(int n, int m);

vector<double> weight(double a, double b, vector<double> x);

double lagrangeBasis(vector<double> X, double x, int i);

int index(const int d, int N, int X, int Y);

int m(int i);

int sum(int **x, int p, int q);

int index_num(const int d, int N);

vector<double> ccpoints(int n, double a, double b);



#endif
