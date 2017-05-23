// basicfuncs.h

//#include <math.h>
//#include <vector>
//#include <stdio.h>

#ifndef _BASICFUNCS_H_
#define _BASICFUNCS_H_

template<class T> double norm(T v);
void standardizeVector(double* V_, int n);
template<class T1, class T2> double dot(T1 v1, T2 v2);
void getRow(double* V_, int nn_, double* v, int idx);
void getRow(double* V_, int nn_, std::vector<double>& v, int idx);
template<class T1, class T2> void cross(T1 v1, T2 v2, double* v);
template<class T1, class T2> void minus_(T1 v1, T2 v2, double* v);
template<class T1, class T2> void minus_(T1 v1, T2 v2, std::vector<double>& v);
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