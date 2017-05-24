// basicfuncs.h

//#include <math.h>
//#include <vector>
//#include <stdio.h>

#ifndef _BASICFUNCS_H_
#define _BASICFUNCS_H_

using namespace std;

double norm(double* v);
double norm(vector<double> v);

void standardizeVector(double* V_, int n);

double dot(double* v1, double* v2);
double dot(vector<double> v1, vector<double> v2);
double dot(double* v1, vector<double> v2);

void getRow(double* V_, int nn_, double* v, int idx);
void getRow(double* V_, int nn_, std::vector<double>& v, int idx);

void cross(double* v1, double* v2, double* v);
void cross(vector<double> v1, vector<double> v2, double* v);
void cross(double* v1, vector<double> v2, double* v);
void cross(vector<double> v1, double* v2, double* v);

void minus_(double* v1, double* v2, double* v);
void minus_(vector<double> v1, vector<double> v2, double* v);

void mul_(double* v, double s);
bool cross_circles(double *c1, double r1, double *u, double *v,
	double *c2, double r2, double *n2, double *p);
void get_uv(double* n, double r, double* u, double* v);
bool isPossibleCrossCircles(double *c1, double r1, double* n1,
	double* c2, double r2, double* n2);
double det(double* A);
bool intersectionPoint(double* p,
	double *c1, double r1, double* n1,
	double *c2, double r2, double* n2,
	double *c3, double r3, double* n3);


#endif // _BASICFUNCS_H_