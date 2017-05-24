// basicfuncs.cpp

#include <math.h>
#include <vector>
#include <stdio.h>

#include "basicfuncs.h"

double norm(double* v)
{
	return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}
double norm(vector<double> v)
{
	return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

void standardizeVector(double* V_, int n)
{
	int j;
	double ns;

	for (j = 0; j < n; j++)
	{
		ns = sqrt(V_[j] * V_[j] + V_[j + n] * V_[j + n] + V_[j + n * 2] * V_[j + n * 2]);

		V_[j] = V_[j] / ns;
		V_[j + n] = V_[j + n] / ns;
		V_[j + n * 2] = V_[j + n * 2] / ns;
	}
}

double dot(double* v1, double* v2)
{
	double val = 0;

	for (int i = 0; i < 3; i++)
	{
		val += v1[i] * v2[i];
	}

	return val;
}
double dot(vector<double> v1, vector<double> v2)
{
	double val = 0;

	for (int i = 0; i < 3; i++)
	{
		val += v1[i] * v2[i];
	}

	return val;
}
double dot(double* v1, vector<double> v2)
{
	double val = 0;

	for (int i = 0; i < 3; i++)
	{
		val += v1[i] * v2[i];
	}

	return val;
}

void getRow(double* V_, int nn_, double* v, int idx)
{
	v[0] = V_[idx];
	v[1] = V_[idx + nn_];
	v[2] = V_[idx + nn_ * 2];
}

void getRow(double* V_, int nn_, std::vector<double>& v, int idx)
{
	v[0] = V_[idx];
	v[1] = V_[idx + nn_];
	v[2] = V_[idx + nn_ * 2];
}

void cross(double* v1, double* v2, double* v)
{
	v[0] = v1[1] * v2[2] - v1[2] * v2[1];
	v[1] = v1[2] * v2[0] - v1[0] * v2[2];
	v[2] = v1[0] * v2[1] - v1[1] * v2[0];
}
void cross(vector<double> v1, vector<double> v2, double* v)
{
	v[0] = v1[1] * v2[2] - v1[2] * v2[1];
	v[1] = v1[2] * v2[0] - v1[0] * v2[2];
	v[2] = v1[0] * v2[1] - v1[1] * v2[0];
}
void cross(double* v1, vector<double> v2, double* v)
{
	v[0] = v1[1] * v2[2] - v1[2] * v2[1];
	v[1] = v1[2] * v2[0] - v1[0] * v2[2];
	v[2] = v1[0] * v2[1] - v1[1] * v2[0];
}
void cross(vector<double> v1, double* v2, double* v)
{
	v[0] = v1[1] * v2[2] - v1[2] * v2[1];
	v[1] = v1[2] * v2[0] - v1[0] * v2[2];
	v[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

void minus_(double* v1, double* v2, double* v)
{
	v[0] = v1[0] - v2[0];
	v[1] = v1[1] - v2[1];
	v[2] = v1[2] - v2[2];
}
void minus_(vector<double> v1, vector<double> v2, double* v)
{
	v[0] = v1[0] - v2[0];
	v[1] = v1[1] - v2[1];
	v[2] = v1[2] - v2[2];
}

template<class T1, class T2>
void minus_(T1 v1, T2 v2, std::vector<double>& v)
{
	v[0] = v1[0] - v2[0];
	v[1] = v1[1] - v2[1];
	v[2] = v1[2] - v2[2];
}

void mul_(double* v, double s)
{
	v[0] *= s;
	v[1] *= s;
	v[2] *= s;
}

bool cross_circles(double *c1, double r1, double *u, double *v,
	double *c2, double r2, double *n2, double *p)
{
	double c[3];
	double nc, nu, nv, a1, a2, a3, d, d_;
	double t1, t2;

	minus_(c1, c2, c);
	nc = dot(n2, c);
	nu = dot(n2, u);
	nv = dot(n2, v);

	a1 = nu*nu + nv*nv;
	a2 = nc*nu;
	a3 = nc*nc - nv*nv;

	d = a2*a2 - a1*a3;
	if (d < 0) return false;

#define CHECK_CROSS()								\
													\
	p[0] = c1[0] + u[0] * t1 + v[0] * t2;			\
	p[1] = c1[1] + u[1] * t1 + v[1] * t2;			\
	p[2] = c1[2] + u[2] * t1 + v[2] * t2;			\
	minus_(p, c2, c);								\
	if (abs(dot(n2, c)) < 1e-12) {					\
		d_ = dot(c, c);								\
		if (d_ < r2*r2) return true;				\
	}												\


	t1 = (-a2 + sqrt(d)) / a1;
	if (abs(t1) <= 1) {
		t2 = sqrt(1 - t1*t1);
		CHECK_CROSS()
			t2 = -t2;
		CHECK_CROSS()
	}
	t1 = (-a2 - sqrt(d)) / a1;
	if (abs(t1) <= 1) {
		t2 = sqrt(1 - t1*t1);
		CHECK_CROSS()
			t2 = -t2;
		CHECK_CROSS()
	}

	return false;
}

void get_uv(double* n, double r, double* u, double* v)
{
	if (n[2] != 0) {
		u[0] = 0; u[1] = 1;
		u[2] = -n[1] / n[2];
	}
	else if (n[1] != 0) {
		u[0] = 0; u[2] = 1;
		u[1] = -n[2] / n[1];
	}
	else {
		u[1] = 0; u[2] = 1;
		u[0] = -n[2] / n[0];
	}
	standardizeVector(u, 1);
	cross(n, u, v);
	mul_(u, r);
	mul_(v, r);
}

bool isPossibleCrossCircles(double *c1, double r1, double* n1,
	double* c2, double r2, double* n2)
{
	double w1[3], w2[3], c[3];
	double d21, d12, rp21, rp12;

	// distance d21 from the j2-th circle to the j1-th circle plane
	minus_(c2, c1, c);
	d21 = abs(dot(n1, c));

	cross(n1, n2, w1);

	cross(w1, n2, w2);
	standardizeVector(w2, 1);
	mul_(w2, r2);
	rp21 = abs(dot(w2, n1));

	d21 -= rp21;
	if (d21 >= 0) return false;

	// distance d12 from the j1-th circle to the j2-th circle plane
	minus_(c1, c2, c);
	d12 = abs(dot(n2, c));

	cross(w1, n1, w2);
	standardizeVector(w2, 1);
	mul_(w2, r1);

	rp12 = abs(dot(w2, n2));
	d12 -= rp12;
	if (d12 >= 0) return false;

	return true;
}

double det(double* A)
{
#define A_(j,i)	A[j*3+i]
	double d;
	d = A_(0, 0) * A_(1, 1) * A_(2, 2) + A_(0, 1) * A_(1, 2) * A_(2, 0) + A_(0, 2) * A_(1, 0) * A_(2, 1)
		- A_(0, 2) * A_(1, 1) * A_(2, 0) - A_(0, 0) * A_(1, 2) * A_(2, 1) - A_(0, 1) * A_(1, 0) * A_(2, 2);
	return d;
}

bool intersectionPoint(double* p,
	double *c1, double r1, double* n1,
	double *c2, double r2, double* n2,
	double *c3, double r3, double* n3)
{
	double d, D_[3][3], T_[3][3];
	double b1, b2, b3;
	double p1[3];

	memcpy(&D_[0][0], n1, sizeof(double) * 3);
	memcpy(&D_[1][0], n2, sizeof(double) * 3);
	memcpy(&D_[2][0], n3, sizeof(double) * 3);

	d = det((double*)D_);
	if (abs(d) < 1e-15) return false;

	b1 = dot(n1, c1);
	b2 = dot(n2, c2);
	b3 = dot(n3, c3);

	memcpy(T_, D_, sizeof(T_));
	T_[0][0] = b1; T_[1][0] = b2; T_[2][0] = b3;
	p[0] = det((double*)T_) / d;

	memcpy(T_, D_, sizeof(T_));
	T_[0][1] = b1; T_[1][1] = b2; T_[2][1] = b3;
	p[1] = det((double*)T_) / d;

	memcpy(T_, D_, sizeof(T_));
	T_[0][2] = b1; T_[1][2] = b2; T_[2][2] = b3;
	p[2] = det((double*)T_) / d;

	minus_(p, c1, p1);
	if (dot(p1, p1) > r1*r1) return false;
	minus_(p, c2, p1);
	if (dot(p1, p1) > r2*r2) return false;
	minus_(p, c3, p1);
	if (dot(p1, p1) > r3*r3) return false;

	return true;
}

