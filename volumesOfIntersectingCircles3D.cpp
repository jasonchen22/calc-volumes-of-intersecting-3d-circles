// volumesOfIntersectingCircles3D.cpp


#include <mex.h>

#include <math.h>
#include <matrix.h>
#include <vector>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <time.h>

#include "basicfuncs.h"

using namespace std;

#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))
static double mpi = 3.14159265358979323846;

#define TRIPLE_NUM_UNIT	(0x800000 / 2)

struct stCrossInfo {
	float	crossPts[3];
	__int32 rev;
	__int64 nIdx;
};

stCrossInfo crossInfo[TRIPLE_NUM_UNIT];
stCrossInfo crossInfo1[TRIPLE_NUM_UNIT];

double *C_, *R, *N_;

int		nn_;
int		Triple[TRIPLE_NUM_UNIT][3];
int		Triple1[TRIPLE_NUM_UNIT][3];
int		Pair_[TRIPLE_NUM_UNIT][2];

__int64 numTriple_;
__int64 numPair_;
__int64 countVolume_;

unsigned char* bpair_;
vector<float> thrpt(3);

char	fnametri[256], fnamepts[256];

vector<vector<int> > vertextris;
vector<vector<double> > vertexpts;
vector<vector<int> > patches;
vector<int> jcirclePatch;
vector<vector<double> > patchdirs;
vector<int> cirpatchflag;
vector<int> vertlapcnt;

vector<vector<int> > vertextristmp;
vector<vector<double> > vertexptstmp;
vector<vector<int> > patchestmp;
vector<int> jcirclePatchtmp;
vector<vector<double> > patchdirstmp;
vector<int> cirpatchflagtmp;
vector<int> vertlapcnttmp;

void calculate_volumes(double* C_, double* R, double* N_);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs != 6) {
		mexErrMsgIdAndTxt("volumesOfIntersectingCircles3D:Input",
			"Six inputs are required.");
	}

	// centers of circles
	if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]))
	{
		mexErrMsgIdAndTxt("volumesOfIntersectingCircles3D:Input",
			"The matrix of centers of circles must be type double.");
	}
	if (mxGetN(prhs[0]) != 3)
	{
		mexErrMsgIdAndTxt("volumesOfIntersectingCircles3D:Input",
			"The matrix of centers of circles must have three columes.");
	}

	// radii of circles
	if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]))
	{
		mexErrMsgIdAndTxt("volumesOfIntersectingCircles3D:Input",
			"The matrix of radii of circles must be type double.");
	}
	if (mxGetN(prhs[1]) != 1 && mxGetM(prhs[1]) != 1)
	{
		mexErrMsgIdAndTxt("volumesOfIntersectingCircles3D:Input",
			"The matrix of radii of circles must be a vector.");
	}

	// normal vectors of circles
	if (!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]))
	{
		mexErrMsgIdAndTxt("volumesOfIntersectingCircles3D:Input",
			"The matrix of normal vectors of circles must be type double.");
	}
	if (mxGetN(prhs[2]) != 3)
	{
		mexErrMsgIdAndTxt("volumesOfIntersectingCircles3D:Input",
			"The matrix of normal vectors of circles must have three columes.");
	}

	if ( mxGetM(prhs[0]) != mxGetM(prhs[2]) ||
		(mxGetM(prhs[0]) != mxGetM(prhs[1]) && mxGetM(prhs[0]) != mxGetN(prhs[1])))
	{
		mexErrMsgIdAndTxt("volumesOfIntersectingCircles3D:Input",
			"The lengths of all inputs must be same.");
	}

	C_ = mxGetPr(prhs[0]);
	R = mxGetPr(prhs[1]);
	N_ = mxGetPr(prhs[2]);
	nn_ = (int)mxGetM(prhs[0]);

	// pair information of interseting circles
	if (!mxIsClass(prhs[3], "uint8") || mxIsComplex(prhs[3]))
	{
		mexErrMsgIdAndTxt("volumesOfIntersectingCircles3D:Input",
			"The matrix of pair information of interseting circles must be type uint8.");
	}
	if (mxGetN(prhs[3]) != nn_)
	{
		mexErrMsgIdAndTxt("volumesOfIntersectingCircles3D:Input",
			"The number of columns of pair matrix must be the same as the number of circles.");
	}

	bpair_ = (unsigned char*)mxGetData(prhs[3]);

	// file (path) name to triples information
	if (!mxIsChar(prhs[4]) || mxGetM(prhs[4]) != 1) {
		mexErrMsgIdAndTxt("volumesOfIntersectingCircles3D:Input",
			"String(filepath to triples info) is requried.");
	}
	mxGetString(prhs[4], fnametri, mxGetN(prhs[4]) + 1);

	// file (path) name to cross points information
	if (!mxIsChar(prhs[5]) || mxGetM(prhs[5]) != 1) {
		mexErrMsgIdAndTxt("volumesOfIntersectingCircles3D:Input",
			"String(filepath to cross points info) is requried.");
	}
	mxGetString(prhs[5], fnamepts, mxGetN(prhs[5]) + 1);

	time_t t1, t2;

	time(&t1);
	// calculate volumes
	calculate_volumes(C_, R, N_);
	time(&t2);
	printf("Time calculating volumes = %d\n", t2 - t1);

	// output
	mwSize dims[2];
	dims[0] = 1; dims[1] = 1;
	plhs[0] = mxCreateNumericArray(2, dims, mxINT64_CLASS, mxREAL);
	*(__int64*)mxGetData(plhs[0]) = countVolume_;
}

void polyhedron_info_real2temp(bool real2temp)
{
	if (real2temp) {
		vertextristmp = vertextris;
		vertexptstmp = vertexpts;
		patchestmp = patches;
		jcirclePatchtmp = jcirclePatch;
		patchdirstmp = patchdirs;
		cirpatchflagtmp = cirpatchflag;
		vertlapcnttmp = vertlapcnt;
	}
	else {
		vertextris = vertextristmp;
		vertexpts = vertexptstmp;
		patches = patchestmp;
		jcirclePatch = jcirclePatchtmp;
		patchdirs = patchdirstmp;
		cirpatchflag = cirpatchflagtmp;
		vertlapcnt = vertlapcnttmp;
	}
}

bool isnecessarypoint(int j1, int j2, int j3)
{
	int m1 = j1*(nn_ / 8 + 1) + j2 / 8;
	if (!(bpair_[m1] & (1 << (j2 % 8)))) return false;

	int m2 = j1*(nn_ / 8 + 1) + j3 / 8;
	if (!(bpair_[m2] & (1 << (j3 % 8)))) return false;

	int m3 = j2*(nn_ / 8 + 1) + j3 / 8;
	if (!(bpair_[m3] & (1 << (j3 % 8)))) return false;

	return true;
}

int calc_nearest_crossPoint(int j1, int j2, int j3, double *p, 
	int& j31, int& j32, double *p_1, double *p_2)
{
	int i, k;
	unsigned char *bpair1, *bpair2;
	unsigned char ch;
	double u[3], uu[3];
	double c1[3], c2[3], c3[3];
	double n1[3], n2[3], n3[3];
	double pp[3];
	double dist1, dist2, dist;

	j31 = -1; j32 = -1;

	getRow(N_, nn_, n1, j1); getRow(C_, nn_, c1, j1);
	getRow(N_, nn_, n2, j2); getRow(C_, nn_, c2, j2);
	
	cross(n1, n2, u);

	dist1 = 1e13; dist2 = 1e13;

	bpair1 = &bpair_[j1*(nn_ / 8 + 1)];
	bpair2 = &bpair_[j2*(nn_ / 8 + 1)];
	for (i = 0; i < nn_ / 8 + 1; i++) {
		ch = (*(bpair1 + i)) & (*(bpair2 + i));
		if (ch == 0) continue;

		for (k = 0; k < 8; k++) {
			if (!(ch & (1 << k))) continue;
			int jj = i * 8 + k;
			if (jj == j3) continue;

			if (!isnecessarypoint(j1, j2, jj)) continue;
			getRow(N_, nn_, n3, jj); getRow(C_, nn_, c3, jj);
			if (!intersectionPoint(pp, c1, R[j1], n1, c2, R[j2], n2, c3, R[jj], n3)) continue;
			
			if ((float)pp[0] < thrpt[0]) continue;
			else if ((float)pp[0] == thrpt[0]) {
				if ((float)pp[1] < thrpt[1]) continue;
				else if ((float)pp[1] == thrpt[1]) {
					if ((float)pp[2] <= thrpt[2]) continue;
				}
			}

			minus_(pp, p, uu);
			dist = sqrt(dot(uu, uu));
			if (dot(u, uu) > 0) {
				if (dist < dist1) { 
					dist1 = dist;
					memcpy(p_1, pp, sizeof(double) * 3);
					j31 = jj;
				}
			}
			else {
				if (dist < dist2) {
					dist2 = dist;
					memcpy(p_2, pp, sizeof(double) * 3);
					j32 = jj;
				}
			}
		}
	}

	if (j31 > -1 && j32 > -1) return 3;
	if (j31 > -1) return 1;
	if (j32 > -1) return 2;

	return 0;
}

void vertex_tri(int j1, int j2, int j3, int& j1_, int& j2_, int& j3_)
{
	j1_ = min(j1, min(j2, j3));
	j3_ = max(j1, max(j2, j3));
	if (j1 > j1_ && j1 < j3_) j2_ = j1;
	if (j2 > j1_ && j2 < j3_) j2_ = j2;
	if (j3 > j1_ && j3 < j3_) j2_ = j3;
}

int insert_vertex2patch(int ivertex, vector<int>& patch)
{
	int j;
	vector<int> vtri;

	if (patch.size() < 2) return 0;
	if (patch[0] == patch.back()) return 1;

	vtri = vertextris[ivertex];

	int cnt1 = 0, cnt2 = 0;
	for (j = 0; j < 3; j++) {
		if (vertextris[patch[0]][j] == vtri[0] ||
			vertextris[patch[0]][j] == vtri[1] ||
			vertextris[patch[0]][j] == vtri[2]) cnt1++;
		if (vertextris[patch.back()][j] == vtri[0] ||
			vertextris[patch.back()][j] == vtri[1] ||
			vertextris[patch.back()][j] == vtri[2]) cnt2++;
	}

	if (cnt1 == 3) {
		if (cnt2 == 2) {
			patch.push_back(ivertex);
			//vertlapcnt[ivertex]++;
			return 1;
		}
		else { return 0; }
	}
	if (cnt2 == 3) {
		if (cnt1 == 2) {
			patch.insert(patch.begin(), ivertex);
			//vertlapcnt[ivertex]++;
			return 1;
		}
		else { return 0; }
	}

	if (cnt2 == 2) patch.push_back(ivertex); 
	else if (cnt1 == 2) patch.insert(patch.begin(), ivertex);
	else return 0;

	vertlapcnt[ivertex]++;

	return 2;
}

bool initial_vertices(int j1, int j2, int j3, double* p)
{
	double pp1[3], pp2[3], *pp, n[3];
	vector<int> vtrip(3);
	vector<double> vpt(3);
	int jj1, jj2;

	vertextris = vector<vector<int> >(4);
	vertexpts = vector<vector<double> >(4);
	jcirclePatch = vector<int>(3);
	patches = vector<vector<int> >(3);
	patchdirs = vector<vector<double> >(3);
	vertlapcnt = vector<int>(4);

	vtrip[0] = j1; vtrip[1] = j2; vtrip[2] = j3;
	vertextris[0] = vtrip;
	vpt[0] = p[0]; vpt[1] = p[1]; vpt[2] = p[2];
	vertexpts[0] = vpt;

	jcirclePatch[0] = j1; patches[0].push_back(0); cirpatchflag[j1] = 1; vertlapcnt[0]++;
	jcirclePatch[1] = j2; patches[1].push_back(0); cirpatchflag[j2] = 2; vertlapcnt[0]++;
	jcirclePatch[2] = j3; patches[2].push_back(0); cirpatchflag[j3] = 3; vertlapcnt[0]++;
	getRow(N_, nn_, n, j1);
	patchdirs[0] = vector<double>(3); patchdirs[0][0] = n[0]; patchdirs[0][1] = n[1]; patchdirs[0][2] = n[2];
	getRow(N_, nn_, n, j2);
	patchdirs[1] = vector<double>(3); patchdirs[1][0] = n[0]; patchdirs[1][1] = n[1]; patchdirs[1][2] = n[2];
	getRow(N_, nn_, n, j3);
	patchdirs[2] = vector<double>(3); patchdirs[2][0] = n[0]; patchdirs[2][1] = n[1]; patchdirs[2][2] = n[2];

	// j1, j2
	int ret = calc_nearest_crossPoint(j1, j2, j3, p, jj1, jj2, pp1, pp2);
	if (ret == 0 || ret == 3) return false;
	if (jj1 > -1) {
		vtrip[0] = j1; vtrip[1] = j2; vtrip[2] = jj1;
		pp = pp1;
	}
	else if (jj2 > -1) {
		vtrip[0] = j1; vtrip[1] = j2; vtrip[2] = jj2;
		pp = pp2;
	}
	vertextris[1] = vtrip;
	vpt[0] = pp[0]; vpt[1] = pp[1]; vpt[2] = pp[2];
	vertexpts[1] = vpt;
	patches[0].push_back(1); vertlapcnt[1]++;
	patches[1].push_back(1); vertlapcnt[1]++;

	// j2, j3
	ret = calc_nearest_crossPoint(j2, j3, j1, p, jj1, jj2, pp1, pp2);
	if (ret == 0 || ret == 3) return false;
	if (jj1 > -1) {
		vtrip[0] = j2; vtrip[1] = j3; vtrip[2] = jj1;
		pp = pp1;
	}
	else if (jj2 > -1) {
		vtrip[0] = j2; vtrip[1] = j3; vtrip[2] = jj2;
		pp = pp2;
	}
	vertextris[2] = vtrip;
	vpt[0] = pp[0]; vpt[1] = pp[1]; vpt[2] = pp[2];
	vertexpts[2] = vpt;
	patches[1].push_back(2); vertlapcnt[2]++;
	patches[2].push_back(2); vertlapcnt[2]++;

	// j3, j1
	ret = calc_nearest_crossPoint(j3, j1, j2, p, jj1, jj2, pp1, pp2);
	if (ret == 0 || ret == 3) return false;
	if (jj1 > -1) {
		vtrip[0] = j1; vtrip[1] = j3; vtrip[2] = jj1;
		pp = pp1;
	}
	else if (jj2 > -1) {
		vtrip[0] = j1; vtrip[1] = j3; vtrip[2] = jj2;
		pp = pp2;
	}
	vertextris[3] = vtrip;
	vpt[0] = pp[0]; vpt[1] = pp[1]; vpt[2] = pp[2];
	vertexpts[3] = vpt;
	patches[2].push_back(3); vertlapcnt[3]++;
	patches[0].push_back(3); vertlapcnt[3]++;

	return true;
}

void initial_arrange_patches()
{
	int i, j, k;
	vector<double> pt0, n(3);
	double u[3], pt1[3], pt2[3];
	vector<int> vpatch;

	for (i = 0; i < patches.size(); i++) {
		int j1 = jcirclePatch[i];
		vpatch = patches[i];
		for (j = 1; j < vertextris.size(); j++) {
			bool inpatch = false;
			for (k = 0; k < vpatch.size(); k++) {
				if (vpatch[k] == j) {
					inpatch = true; break;
				}
			}
			if (!inpatch) {
				pt0 = vertexpts[j]; break;
			}
		}

		minus_(pt0, vertexpts[vpatch[0]], pt1);
		if (dot(pt1, patchdirs[i]) < 0) {
			patchdirs[i][0] = -patchdirs[i][0];
			patchdirs[i][1] = -patchdirs[i][1];
			patchdirs[i][2] = -patchdirs[i][2];
		}
		minus_(vertexpts[vpatch[0]], vertexpts[vpatch[1]], pt1);
		minus_(vertexpts[vpatch[2]], vertexpts[vpatch[0]], pt2);
		cross(pt1, pt2, u);
		if (dot(u, patchdirs[i]) > 0) {
			patches[i][0] = vpatch[1];
			patches[i][1] = vpatch[0];
			patches[i][2] = vpatch[2];
		}
		else {
			patches[i][0] = vpatch[2];
			patches[i][1] = vpatch[0];
			patches[i][2] = vpatch[1];
		}
	}

	for (i = 1; i < vertextris.size(); i++) {
		j = vertextris[i][2];
		k = abs(cirpatchflag[j]);
		if (k > 0) continue;

		jcirclePatch.push_back(j);
		patches.push_back(vector<int>(1, i)); vertlapcnt[i]++;
		getRow(N_, nn_, n, j);
		patchdirs.push_back(n);
		cirpatchflag[j] = -(int)jcirclePatch.size();
	}

}

bool check_closed_polyhedron(unsigned char* edgeflag)
{
	int i, j;
	vector<vector<int> > unconnectedsides;

	for (i = 0; i < patches.size(); i++) {
		if (patches[i][0] != patches[i].back()) return false;
	}

	int vertsz = vertextris.size();
	memset(edgeflag, 0, vertsz*vertsz);

	for (i = 0; i < patches.size(); i++) {
		for (j = 0; j < patches[i].size()-1; j++) {
			edgeflag[patches[i][j] * vertsz + patches[i][j + 1]]++;
			edgeflag[patches[i][j + 1] * vertsz + patches[i][j]]++;
		}
	}

	for (i = 0; i < patches.size(); i++) {
		for (j = 0; j < patches[i].size() - 1; j++) {
			if (edgeflag[patches[i][j] * vertsz + patches[i][j + 1]] > 1) continue;
			vector<int> edge(2);
			edge[0] = patches[i][j]; edge[1] = patches[i][j + 1];
			unconnectedsides.push_back(edge);
		}
	}
	
	if (unconnectedsides.size() == 0) return true;

	return false;
}

bool determine_one_point(int jcircle, int iv, vector<vector<int> > jjs, vector<vector<double> > pps,
	vector<int>& jj, vector<double>& pp)
{
	int i, k;
	double v1[3], v2[3], u[3];
	double al1, alm;
	vector<double> pt1, pt2;

	pp = vector<double>(3);

	if (cirpatchflag[jcircle] < 0) return false; // the case that direction of patch isn't determined yet

	k = cirpatchflag[jcircle] - 1;
	if (patches[k].back() == iv)
	{
		pt1 = vertexpts[iv];
		pt2 = vertexpts[patches[k][patches[k].size() - 2]];

		minus_(pt2, pt1, v1);
		cross(v1, patchdirs[k], u);

		alm = 2 * mpi;
		for (i = 0; i < jjs.size(); i++) {
			minus_(pps[i], pt1, v2);
			double vcos = dot(v1, v2) / (norm(v1)*norm(v2));
			if (vcos < -1) vcos = -0.999999999999;
			if (vcos > 1) vcos = 0.999999999999;
			al1 = acos(vcos);
			if (dot(u, v2) < 0) al1 = mpi + mpi - al1;
			if (al1 < alm) {
				alm = al1;
				jj = jjs[i];
				pp = pps[i];
			}
		}
	}
	else if (patches[k][0] == iv) 
	{
		pt1 = vertexpts[iv];
		pt2 = vertexpts[patches[k][1]];

		minus_(pt2, pt1, v1);
		cross(v1, patchdirs[k], u);

		alm = 0;
		for (i = 0; i < jjs.size(); i++) {
			minus_(pps[i], pt1, v2);
			double vcos = dot(v1, v2) / (norm(v1)*norm(v2));
			if (vcos < -1) vcos = -0.999999999999;
			if (vcos > 1) vcos = 0.999999999999;
			al1 = acos(vcos);
			if (dot(u, v2) < 0) al1 = mpi + mpi - al1;
			if (al1 > alm) {
				alm = al1;
				jj = jjs[i];
				pp = pps[i];
			}
		}
	}
	else { 
		return false; 
	}

	return true;
}

int insert_polyhedron_info(int jprincir, int jcir1, int jcir2, double *p,
	vector<vector<int> > jjs, vector<vector<double> > pps, int iv)
{
	int i, j, k, m;

	vector<int> vtri(3);
	vector<double> pt1, pt2, pp(3), n(3);

	int j1_, j2_, j3_, jj1, jj2;
	vector<int> jj;

	if (jjs.size() == 0) return 0;
	if (cirpatchflag[jprincir] < 0) return 0; // the case that direction of patch isn't determined yet

	if (jjs.size() > 1) {
		bool bdetermined = false;
		k = cirpatchflag[jprincir] - 1;
		if (patches[k].size() >= 2) {
			bdetermined = determine_one_point(jprincir, iv, jjs, pps, jj, pp);
			jj1 = jj[0]; jj2 = jj[1];
		}

		if (!bdetermined) return 0;
	}
	else {
		jj1 = jjs[0][0]; jj2 = jjs[0][1];
		pp = pps[0];
	}

	vertex_tri(jprincir, jj1, jj2, j1_, j2_, j3_);
	for (i = 0; i < vertextris.size(); i++) {
		vertex_tri(vertextris[i][0], vertextris[i][1], vertextris[i][2],
			vtri[0], vtri[1], vtri[2]);
		if (vtri[0] == j1_ && vtri[1] == j2_ && vtri[2] == j3_) 
		{
			int ninsert = 0;
			m = cirpatchflag[jprincir] - 1;
			for (j = 0; j < patches[m].size(); j++) {
				if (patches[m][j] == i) break;
			}
			if (j > 0 && j < patches[m].size() - 1) {
				patches[m] = vector<int>(patches[m].begin() + j, patches[m].end());
			}
			ninsert = insert_vertex2patch(i, patches[m]);

			return ninsert;
		}
	}

	i = vertextris.size();
	vtri[0] = jprincir; vtri[1] = jj1; vtri[2] = jj2;
	vertextris.push_back(vtri);
	vertexpts.push_back(pp);
	vertlapcnt.push_back(0);

	k = cirpatchflag[jprincir] - 1;
	int ret = insert_vertex2patch(i, patches[k]);
	if (ret == 0) return 0;

	k = cirpatchflag[jj2];
	if (k == 0) {
		patches.push_back(vector<int>(1, i)); vertlapcnt[i]++;
		jcirclePatch.push_back(jj2);
		getRow(N_, nn_, n, jj2);
		patchdirs.push_back(n);
		cirpatchflag[jj2] = -(int)jcirclePatch.size();
	}

	return ret;
}

bool check_3vertices_oneline(vector<int> vert0, vector<int> vert1, vector<int> vert2)
{
	int k, cnt;
	vector<int> comjs(2);

	cnt = 0;
	for (k = 0; k < 3; k++) {
		if (vert0[k] == vert1[0] ||
			vert0[k] == vert1[1] ||
			vert0[k] == vert1[2])
		{
			comjs[cnt] = vert0[k];
			cnt++;
			if (cnt == 2) break;
		}
	}
	if (cnt < 2) return false;

	cnt = 0;
	for (k = 0; k < 2; k++) {
		if (comjs[k] == vert2[0] ||
			comjs[k] == vert2[1] ||
			comjs[k] == vert2[2])
		{
			cnt++;
			if (cnt == 2) break;
		}
	}
	if (cnt < 2) return false;

	return true;
}

void remove_middle_point_of_patch()
{
	int i, j;
	vector<int> patch;
	vector<int> vert0, vert1, vert2;

	for (i = 0; i < patches.size(); i++) {
		patch = patches[i];
		if (patch.size() < 4) continue;
		for (j = 0; j < patch.size() - 1; j++) {
			if (patch.size() <= 4) break;
			vert0 = vertextris[patch[j]];
			if (j == 0) vert1 = vertextris[patch[patch.size() - 2]];
			else vert1 = vertextris[patch[j - 1]];
			vert2 = vertextris[patch[j + 1]];

			if (!check_3vertices_oneline(vert0, vert1, vert2)) continue;
			patch.erase(patch.begin() + j);
			if (j == 0) {
				patch.erase(patch.end() - 1);
				patch.push_back(patch[0]);
			}
			j--;
		}
		patches[i] = patch;
	}
}

bool connectable_two_vertices(int jprincir, int ivert1, int ivert2)
{
	int i;
	int jjoincir;

	vector<int> vertex1 = vertextris[ivert1];
	vector<int> vertex2 = vertextris[ivert2];
	
	for (i = 0; i < 3; i++) {
		if (vertex1[i] == jprincir) break;
	}
	if (i == 3) return false;

	for (i = 0; i < 3; i++) {
		if (vertex2[i] == jprincir) break;
	}
	if (i == 3) return false;

	jjoincir = -1;
	for (i = 0; i < 3; i++) {
		if (vertex1[i] == jprincir) continue;
		if (vertex1[i] == vertex2[0] ||
			vertex1[i] == vertex2[1] ||
			vertex1[i] == vertex2[2])
		{
			jjoincir = vertex1[i]; break;
		}
	}
	if (jjoincir == -1) return false;

	if (cirpatchflag[jjoincir] <= 0) return true;
	int ipat = cirpatchflag[jjoincir] - 1;
	
	vector<int> vpatch = patches[ipat];

	for (i = 0; i < vpatch.size(); i++) {
		if (vpatch[i] == ivert1) break;
	}
	if (i == vpatch.size()) return false;

	if (i == 0) {
		if (vpatch[i + 1] == ivert2) return true;
		if (vpatch[0] == vpatch.back()) {
			if (vpatch[vpatch.size() - 2] == ivert2) return true;
		}
	}
	else if (i == vpatch.size() - 1) {
		if (vpatch[i - 1] == ivert2) return true;
		if (vpatch[0] == vpatch.back()) {
			if (vpatch[1] == ivert2) return true;
		}
	}
	else {
		if (vpatch[i - 1] == ivert2||
			vpatch[i + 1] == ivert2) return true;
	}

	return false;
}

void insert_middle_point_in_patch(vector<int>& patch, int jprincir)
{
	int i, j, k, m, iv, iv1;
	int j1, j2, j3, j31, j32;
	double c1[3], c2[3], c3[3];
	double n1[3], n2[3], n3[3];
	double pp[3], u[3];
	vector<vector<int> > newverttris;
	vector<int> tri(3);
	vector<vector<double> > newvertpts;
	vector<double> vpt(3);
	unsigned char *bpair1, *bpair2;
	unsigned char ch;

	j1 = jprincir;
	bpair1 = &bpair_[j1*(nn_ / 8 + 1)];
	getRow(N_, nn_, n1, j1); getRow(C_, nn_, c1, j1);

	for (m = 0; m < patch.size() - 1; m++) {
		iv = patch[m];
		iv1 = patch[m + 1];
		for (k = 0; k < 3; k++) {
			if (vertextris[iv][k] == j1) continue;
			if (vertextris[iv][k] == vertextris[iv1][0] ||
				vertextris[iv][k] == vertextris[iv1][1] ||
				vertextris[iv][k] == vertextris[iv1][2])
			{
				j2 = vertextris[iv][k]; break;
			}
		}
		for (k = 0; k < 3; k++) {
			if (vertextris[iv][k] == j1 || vertextris[iv][k] == j2) continue;
			j3 = vertextris[iv][k]; break;
		}

		bpair2 = &bpair_[j2*(nn_ / 8 + 1)];
		getRow(N_, nn_, n2, j2); getRow(C_, nn_, c2, j2);
		cross(n1, n2, u);
		
		newverttris.clear();
		newvertpts.clear();

		for (i = 0; i < nn_ / 8 + 1; i++) {
			ch = (*(bpair1 + i)) & (*(bpair2 + i));
			if (ch == 0) continue;

			for (k = 0; k < 8; k++) {
				if (!(ch & (1 << k))) continue;
				int jj = i * 8 + k;
				if (jj == j3) continue;

				if (!isnecessarypoint(j1, j2, jj)) continue;
				getRow(N_, nn_, n3, jj); getRow(C_, nn_, c3, jj);
				if (!intersectionPoint(pp, c1, R[j1], n1, c2, R[j2], n2, c3, R[jj], n3)) continue;

				if (pp[0] < min(vertexpts[iv][0], vertexpts[iv1][0]) || max(vertexpts[iv][0], vertexpts[iv1][0]) > pp[0] ||
					pp[1] < min(vertexpts[iv][1], vertexpts[iv1][1]) || max(vertexpts[iv][1], vertexpts[iv1][1]) > pp[1] ||
					pp[2] < min(vertexpts[iv][2], vertexpts[iv1][2]) || max(vertexpts[iv][2], vertexpts[iv1][0]) > pp[2])
					continue;
				
				tri[0] = j1; tri[1] = j2; tri[2] = jj;
				newverttris.push_back(tri);
				vpt[0] = pp[0]; vpt[1] = pp[1]; vpt[2] = pp[2];
				newvertpts.push_back(vpt);
			}
		}

		if (newverttris.size() == 0) continue;
		while (newverttris.size() > 0) 
		{
			i = 0;
			minus_(vertexpts[iv], newvertpts[0], pp);
			double dmx = sqrt(dot(pp, pp));
			for (k = 1; k < newverttris.size(); k++) {
				minus_(vertexpts[iv], newvertpts[k], pp);
				double d = sqrt(dot(pp, pp));
				if (d >= dmx) continue;
				dmx = d;
				i = k;
			}

			j = vertextris.size();
			vertextris.push_back(newverttris[i]);
			vertexpts.push_back(newvertpts[i]);
			patch.insert(patch.begin() + m + 1, j);
			m++;

			newverttris.erase(newverttris.begin() + i);
			newvertpts.erase(newvertpts.begin() + i);
		}
	}
}

bool determine_patch_direction(int ipatch)
{
	int i, j, k;
	vector<int> patchtmp, patchtmp1;

	int jprincir = jcirclePatch[ipatch];
	int jjoinedcir;
	int ipat1;

	for (i = 0; i < vertextris.size(); i++) {
		if (vertextris[i][0] != jprincir &&
			vertextris[i][1] != jprincir &&
			vertextris[i][2] != jprincir) continue;
		patchtmp.push_back(i);
	}
	if (patchtmp.size() < 3) return false;

	patchtmp1 = vector<int>(1);
	k = 0;
	patchtmp1[0] = patchtmp[0];

	int isel = -1;
	while (true) {
		int cnt = 0;
		for (i = 0; i < patchtmp.size(); i++) {
			if (k == i) continue;
			cnt = 0;
			for (j = 0; j < 3; j++) {
				if (vertextris[patchtmp[i]][j] == vertextris[patchtmp1[0]][0] ||
					vertextris[patchtmp[i]][j] == vertextris[patchtmp1[0]][1] ||
					vertextris[patchtmp[i]][j] == vertextris[patchtmp1[0]][2]) cnt++;
			}
			if (cnt != 2) continue;

			int iv1 = patchtmp1[0];
			int iv2 = patchtmp[i];
			for (j = 0; j < 3; j++) {
				if (vertextris[iv1][j] == jcirclePatch[ipatch]) continue;
				if (vertextris[iv1][j] == vertextris[iv2][0] ||
					vertextris[iv1][j] == vertextris[iv2][1] ||
					vertextris[iv1][j] == vertextris[iv2][2])
				{
					jjoinedcir = vertextris[iv1][j]; break;
				}
			}
			if (cirpatchflag[jjoinedcir] <= 0) continue;

			ipat1 = cirpatchflag[jjoinedcir] - 1;
			for (j = 1; j < patches[ipat1].size() - 1; j++) {
				if (patches[ipat1][j] == iv1) {
					if (patches[ipat1][j - 1] == iv2) {
						patchtmp1.push_back(iv2); break;
					}
					if (patches[ipat1][j + 1] == iv2) {
						patchtmp1.insert(patchtmp1.begin(), iv2); break;
					}
				}
				if (patches[ipat1][j] == iv2) {
					if (patches[ipat1][j - 1] == iv1) {
						patchtmp1.insert(patchtmp1.begin(), iv2); break;
					}
					if (patches[ipat1][j + 1] == iv1) {
						patchtmp1.push_back(iv2); break;
					}
				}
			}
			if (patchtmp1.size() != 2) continue;
			isel = i;
			break;
		}
		if (patchtmp1.size() == 2) break;

		k++;
		if (k == patchtmp.size()) return false;
		patchtmp1[0] = patchtmp[k];
	}

	if (patchtmp1.size() != 2) return false;
	patchtmp.erase(patchtmp.begin() + max(k, isel));
	patchtmp.erase(patchtmp.begin() + min(k, isel));

	while (patchtmp.size() > 0) {
		bool binserted = false;
		for (i = 0; i < patchtmp.size(); i++) {
			int cnt1 = 0, cnt2 = 0;
			for (j = 0; j < 3; j++) {
				if (vertextris[patchtmp[i]][j] == vertextris[patchtmp1[0]][0] ||
					vertextris[patchtmp[i]][j] == vertextris[patchtmp1[0]][1] ||
					vertextris[patchtmp[i]][j] == vertextris[patchtmp1[0]][2]) cnt1++;
				if (vertextris[patchtmp[i]][j] == vertextris[patchtmp1.back()][0] ||
					vertextris[patchtmp[i]][j] == vertextris[patchtmp1.back()][1] ||
					vertextris[patchtmp[i]][j] == vertextris[patchtmp1.back()][2]) cnt2++;
			}
			if (cnt1 != 2 && cnt2 != 2) continue;
			if (cnt1 == 2) {
				if (connectable_two_vertices(jprincir, patchtmp1[0], patchtmp[i])) {
					patchtmp1.insert(patchtmp1.begin(), patchtmp[i]);
					binserted = true;  break;
				}
			}
			if (cnt2 == 2) {
				if (connectable_two_vertices(jprincir, patchtmp1.back(), patchtmp[i])) {
					patchtmp1.push_back(patchtmp[i]);
					binserted = true; break;
				}
			}
		}
		if (!binserted) break;
		patchtmp.erase(patchtmp.begin() + i);
	}

	if (patchtmp1.size() < 3) return false;

	if (patchtmp1[0] != patchtmp1.back()) 
	{
		for (i = 0; i < patchtmp1.size() - 1; i++) {
			int j1, j2;
			j1 = jprincir;
			for (k = 0; k < 3; k++) {
				if (vertextris[patchtmp1[i]][k] == j1) continue;
				if (vertextris[patchtmp1[i]][k] == vertextris[patchtmp1[i + 1]][0] ||
					vertextris[patchtmp1[i]][k] == vertextris[patchtmp1[i + 1]][1] ||
					vertextris[patchtmp1[i]][k] == vertextris[patchtmp1[i + 1]][2])
				{
					j2 = vertextris[patchtmp1[i]][k]; break;
				}
			}
			int cnt = 2;
			vector<int> verts(0);
			verts.push_back(patchtmp1[i]);
			verts.push_back(patchtmp1[i + 1]);
			for (j = 0; j < patchtmp1.size(); j++) {
				if (j == i || j == i + 1) continue;
				for (k = 0; k < 3; k++) {
					if (vertextris[patchtmp1[j]][k] == j2) break;
				}
				if (k == 3) continue;
				if (j > i + 1) {
					if (verts.back() != patchtmp1[j - 1]) cnt++;
					verts.push_back(patchtmp1[j]);
				}
				if (j < i) {
					if (verts[0] != patchtmp1[j + 1]) cnt++;
					verts.insert(verts.begin(), patchtmp1[j]);
				}
			}
			if (cnt < 3) continue;

			patchtmp = patchtmp1;
			patchtmp1.clear();
			int pcnt = 0;
			for (j = 0; j < patchtmp.size(); j++) {
				for (k = 0; k < verts.size(); k++) {
					if (patchtmp[j] == verts[k]) break;
				}
				if (k < verts.size()) {
					if (patchtmp1.size() >= 2) {
						patchtmp1.push_back(patchtmp[j]);
						break;
					}
					else {
						patchtmp1.clear();
					}
					if (pcnt < 1) pcnt++;
					continue;
				}
				patchtmp1.push_back(patchtmp[j]);
				if (pcnt == 1) break;
				else if (pcnt > 1) pcnt = 1;
			}
			if (pcnt == 0) break;

			for (k = j + 1; k < patchtmp.size(); k++) {
				patchtmp1.push_back(patchtmp[k]);
				int k1;
				for (k1 = 0; k1 < verts.size(); k1++) {
					if (patchtmp[k] == verts[k1]) break;
				}
				if (k1 < verts.size()) {
					verts.erase(verts.begin() + k1);
					break;
				}
			}
			for (k = j - 1; k >= 0; k--) {
				patchtmp1.insert(patchtmp1.begin(), patchtmp[k]);
				int k1;
				for (k1 = 0; k1 < verts.size(); k1++) {
					if (patchtmp[k] == verts[k1]) break;
				}
				if (k1 < verts.size()) break;
			}
			break;
		}
	}
	if (patchtmp1.size() < 3) 
		return false;

	for (i = 0; i < patchtmp1.size(); i++) {
		vector<int> vert0, vert1, vert2;
		vector<double> pt1, pt2, pt0;
		if (i == 0) {
			vert0 = vertextris[patchtmp1[i + 1]];
			vert1 = vertextris[patchtmp1[i]];
			vert2 = vertextris[patchtmp1[i + 2]];

			pt0 = vertexpts[patchtmp1[i + 1]];
			pt1 = vertexpts[patchtmp1[i]];
			pt2 = vertexpts[patchtmp1[i + 2]];
		}
		else if (i == patchtmp1.size() - 1) {
			vert0 = vertextris[patchtmp1[patchtmp1.size() - 2]];
			vert1 = vertextris[patchtmp1[patchtmp1.size() - 3]];
			vert2 = vertextris[patchtmp1.back()];

			pt0 = vertexpts[patchtmp1[patchtmp1.size() - 2]];
			pt1 = vertexpts[patchtmp1[patchtmp1.size() - 3]];
			pt2 = vertexpts[patchtmp1.back()];
		}
		else {
			vert0 = vertextris[patchtmp1[i]];
			vert1 = vertextris[patchtmp1[i - 1]];
			vert2 = vertextris[patchtmp1[i + 1]];

			pt0 = vertexpts[patchtmp1[i]];
			pt1 = vertexpts[patchtmp1[i - 1]];
			pt2 = vertexpts[patchtmp1[i + 1]];
		}

		if (check_3vertices_oneline(vert0, vert1, vert2)) continue;

		double u1[3], u2[3], u[3];
		minus_(pt0, pt1, u1);
		minus_(pt2, pt0, u2);
		cross(u1, u2, u);
		if (dot(u, patchdirs[ipatch]) < 0) {
			patchdirs[ipatch][0] = -patchdirs[ipatch][0];
			patchdirs[ipatch][1] = -patchdirs[ipatch][1];
			patchdirs[ipatch][2] = -patchdirs[ipatch][2];
		}

		break;
	}
	if (i == patchtmp1.size())
		return false;

	cirpatchflag[jprincir] = abs(cirpatchflag[jprincir]);

	insert_middle_point_in_patch(patchtmp1, jprincir);

	patches[ipatch] = patchtmp1;
	for (i = 0; i < patchtmp1.size(); i++) {
		vertlapcnt[patchtmp1[i]]++;
	}

	return true;
}

bool check_positive_polygon(vector<int> patch, vector<double> patchdir)
{
	int i;
	vector<double> pt0, pt1, pt2;
	double u1[3], u2[3], n[3], n1[3];

	n[0] = 0; n[1] = 0; n[2] = 0;
	for (i = 0; i < patch.size() - 1; i++) {
		if (i == 0) pt0 = vertexpts[patch[patch.size() - 2]];
		else pt0 = vertexpts[patch[i - 1]];
		pt1 = vertexpts[patch[i]];
		pt2 = vertexpts[patch[i + 1]];

		minus_(pt1, pt0, u1);
		minus_(pt2, pt1, u2);
		cross(u1, u2, n1);

		n[0] += n1[0]; n[1] += n1[1]; n[2] += n1[2];
	}

	if (dot(n, patchdir) < 0) return false;

	return true;
}

bool merge_edge(vector<int>& patch1, vector<int> patch2, int edgept1 = -1, int edgept2 = -1)
{
	int i, j, k, m;
	vector<int> edge(2);
	bool ret = false;

	for (i = 0; i < patch1.size() - 1; i++) {
		if (edgept1 >= 0 && edgept2 >= 0) {
			if ((patch1[i] != edgept1 || patch1[i + 1] != edgept2) &&
				(patch1[i] != edgept2 || patch1[i + 1] != edgept1)) continue;
		}

		for (j = 0; j < patch2.size() - 1; j++) {
			if (patch1[i] == patch2[j + 1] && patch1[i + 1] == patch2[j]) break;
		}
		if (j == patch2.size() - 1) continue;

		m = 1;
		for (k = j + 2; k < patch2.size(); k++) {
			if (patch2[k] == patch1[i + m]) break;
			patch1.insert(patch1.begin() + i + m, patch2[k]);
			m++;
		}
		if (patch2[k] == patch1[i + m]) {
			ret = true;
			break;
		}
		for (k = 1; k < j; k++) {
			if (patch2[k] == patch1[i + m]) break;
			patch1.insert(patch1.begin() + i + m, patch2[k]);
			m++;
		}

		ret = true;
		break;
	}
	if (!ret) return false;

	m--;
	for (k = i + 1; k < i + m; k++) {
		int id = patch1[k];
		for (j = i - 1; j >= 0; j--) {
			if (id == patch1[j]) {
				for (int i1 = i + 1; i1 <= k; i1++) {
					patch1.erase(patch1.begin() + i + 1);
					m--;
				}
				for (int i1 = i; i1 > j; i1--) {
					patch1.erase(patch1.begin() + i1);
				}
				break;
			}
		}
		if (j >= 0) {
			i = j; k = i; 
			continue;
		}

		for (j = patch1.size() - 2; j >= i + m; j--) {
			if (id == patch1[j]) {
				for (int i1 = patch1.size() - 1; i1 > j; i1--) {
					patch1.erase(patch1.begin() + i1);
				}
				for (int i1 = i + 1; i1 < k; i1++) {
					patch1.erase(patch1.begin() + i + 1);
					m--;
				}
				for (int i1 = i; i1 >= 0; i1--) {
					patch1.erase(patch1.begin() + i1);
				}
				break;
			}
		}
		if (j >= i + m) {
			i = 0; k = i;
			continue;
		}
	}

	m++;
	for (k = i + m - 1; k > i; k--) {
		int id = patch1[k];
		for (j = i + m + 1; j < patch1.size(); j++) {
			if (id == patch1[j]) {
				for (int i1 = i + m; i1 < j; i1++) {
					patch1.erase(patch1.begin() + i + m);
				}
				for (int i1 = i + m - 1; i1 >= k; i1--) {
					patch1.erase(patch1.begin() + i1);
					m--;
				}
				break;
			}
		}
		if (j < patch1.size()) {
			k = i + m;
			continue;
		}

		for (j = 1; j < i + 1; j++) {
			if (id == patch1[j]) {
				for (int i1 = i + m - 1; i1 > k; i1--) {
					patch1.erase(patch1.begin() + i + m - 1);
					m--;
				}
				for (int i1 = i + m; i1 < patch1.size(); i1++) {
					patch1.erase(patch1.begin() + i + m);
				}
				for (int i1 = 0; i1 < j; i1++) {
					patch1.erase(patch1.begin());
				}
				break;
			}
		}
		if (j < i + 1) {
			i -= j; k = i + m;
			continue;
		}
	}

	return ret;
}

bool add_patches(unsigned char*& edgeflag)
{
	int i, j, k, m, ipat;
	int iv, iv1;
	vector<int> flag(vertextris.size(), 0);
	vector<int> tri;
	vector<double> n(3);
	double p[3], pp1[3], pp2[3];
	int j1, j2, j3, jj1, jj2;

	vector<vector<int> > patchpts(patches.size());
	vector<int> vpatch;

	int vertsz = vertextris.size();

	for (i = 0; i < vertextris.size(); i++) {
		for (k = 0; k < 3; k++) {
			j = vertextris[i][k];
			if (cirpatchflag[j] == 0) continue;
			ipat = abs(cirpatchflag[j]) - 1;
			patchpts[ipat].push_back(i);
		}
	}

	for (m = 0; m < patchpts.size(); m++) {
		if (patchpts[m].size() == patches[m].size() - 1) continue;

		vector<vector<int> > edgetmps;
		for (j = 0; j < patches[m].size() - 1; j++) {
			if (edgeflag[patches[m][j] * vertsz + patches[m][j + 1]] != 1) continue;
			vector<int> tmp(2);
			tmp[0] = patches[m][j]; tmp[1] = patches[m][j + 1];
			edgetmps.push_back(tmp);
		}

		vector<int> restpts = patchpts[m];
		for (k = 0; k < restpts.size(); k++) {
			int k1;
			for (k1 = 0; k1 < patches[m].size(); k1++) {
				if (restpts[k] == patches[m][k1]) break;
			}
			if (k1 < patches[m].size()) {
				restpts.erase(restpts.begin() + k);
				k--;
				continue;
			}

			vector<int> vert0, vert1, vert2;
			vert0 = vertextris[restpts[k]];
			for (k1 = 0; k1 < patches[m].size() - 1; k1++) {
				vert1 = vertextris[patches[m][k1]];
				vert2 = vertextris[patches[m][k1 + 1]];

				if (!check_3vertices_oneline(vert0, vert1, vert2)) continue;

				minus_(vertexpts[patches[m][k1]], vertexpts[restpts[k]], pp1);
				minus_(vertexpts[patches[m][k1 + 1]], vertexpts[restpts[k]], pp2);
				if (dot(pp1, pp2) < 0) {
					restpts.erase(restpts.begin() + k);
					k--; break;
				}
			}
		}
		if (restpts.size() == 0) continue;

		for (j = 0; j < edgetmps.size(); j++) {
			if (edgeflag[edgetmps[j][0] * vertsz + edgetmps[j][1]] != 1) continue;

			iv = edgetmps[j][0]; iv1 = edgetmps[j][1];
			j1 = jcirclePatch[i];
			j2 = -1;
			for (int i1 = 0; i1 < 3; i1++) {
				if (vertextris[iv][i1] == j1) continue;
				if (vertextris[iv][i1] == vertextris[iv1][0] ||
					vertextris[iv][i1] == vertextris[iv1][1] ||
					vertextris[iv][i1] == vertextris[iv1][2])
				{
					j2 = vertextris[iv][i1]; break;
				}
			}
			if (j2 == -1) continue;

			vpatch = patches[m];
			vector<int> vpatchpt = restpts;

			patches[m].clear();
			patches[m].push_back(edgetmps[j][1]);
			patches[m].push_back(edgetmps[j][0]);

			vector<int> vert0, vert1, vert2;
			vert1 = vertextris[patches[m][0]];
			vert2 = vertextris[patches[m].back()];
			for (k = 0; k < vpatchpt.size(); k++) {
				vert0 = vertextris[vpatchpt[k]];
				int cnt1 = 0, cnt2 = 0;
				for (int k1 = 0; k1 < 3; k1++) {
					if (vert0[k1] == vert1[0] ||
						vert0[k1] == vert1[1] ||
						vert0[k1] == vert1[2]) cnt1++;
					if (vert0[k1] == vert2[0] ||
						vert0[k1] == vert2[1] ||
						vert0[k1] == vert2[2]) cnt2++;
				}
				if (cnt1 != 2 && cnt2 != 2) continue;
				
				if (check_3vertices_oneline(vert0, vert1, vert2)) continue;
				
				if (cnt2 == 2) {
					if (edgeflag[patches[m].back() *vertsz + vpatchpt[k]] >= 2) continue;
					patches[m].push_back(vpatchpt[k]);
				}
				else if (cnt1 == 2) {
					if (edgeflag[patches[m][0] *vertsz + vpatchpt[k]] >= 2) continue;
					patches[m].insert(patches[m].begin(), vpatchpt[k]);
				}
				vpatchpt.erase(vpatchpt.begin() + k);
				k = 0;
			}

			if (patches[m].size() == 2) {
				patches[m] = vpatch;
				continue;
			}

			vector<vector<int> > jjs;
			vector<int> jj(2);
			vector<vector<double> > pps;
			vector<double> ppp(3);
			while (patches[m].back() != patches[m][0]) {
				iv = patches[m].back();
				iv1 = patches[m][patches[m].size() - 2];

				// j1 - principal circle plane
				j1 = jcirclePatch[m];
				// j2 - previous circle plane
				for (int i1 = 0; i1 < 3; i1++) {
					if (vertextris[iv][i1] == j1) continue;
					if (vertextris[iv][i1] == vertextris[iv1][0] ||
						vertextris[iv][i1] == vertextris[iv1][1] ||
						vertextris[iv][i1] == vertextris[iv1][2])
					{
						j2 = vertextris[iv][i1]; break;
					}
				}
				// j3 - new point should be found in the direction (j1,j3)
				for (int i1 = 0; i1 < 3; i1++) {
					if (vertextris[iv][i1] == j1 || vertextris[iv][i1] == j2) continue;
					j3 = vertextris[iv][i1]; break;
				}

				p[0] = vertexpts[iv][0]; p[1] = vertexpts[iv][1]; p[2] = vertexpts[iv][2];

				jjs.clear(); pps.clear();
				int ret = calc_nearest_crossPoint(j1, j3, j2, p, jj1, jj2, pp1, pp2);
				if (jj1 > -1) {
					jj[0] = j3; jj[1] = jj1;
					jjs.push_back(jj);
					ppp[0] = pp1[0]; ppp[1] = pp1[1]; ppp[2] = pp1[2];
					pps.push_back(ppp);
				}
				if (jj2 > -1) {
					jj[0] = j3; jj[1] = jj2;
					jjs.push_back(jj);
					ppp[0] = pp2[0]; ppp[1] = pp2[1]; ppp[2] = pp2[2];
					pps.push_back(ppp);
				}
				ret = calc_nearest_crossPoint(j1, j2, j3, p, jj1, jj2, pp1, pp2);
				int j1t, j2t, j3t, j1s, j2s, j3s;
				vertex_tri(vertextris[iv1][0], vertextris[iv1][1], vertextris[iv1][2], j1t, j2t, j3t);
				if (jj1 > -1) {
					vertex_tri(j1, j2, jj1, j1s, j2s, j3s);
					if (j1t != j1s || j2t != j2s || j3t != j3s) {
						jj[0] = j2; jj[1] = jj1;
						jjs.push_back(jj);
						ppp[0] = pp1[0]; ppp[1] = pp1[1]; ppp[2] = pp1[2];
						pps.push_back(ppp);
					}
				}
				if (jj2 > -1) {
					vertex_tri(j1, j2, jj2, j1s, j2s, j3s);
					if (j1t != j1s || j2t != j2s || j3t != j3s) {
						jj[0] = j2; jj[1] = jj2;
						jjs.push_back(jj);
						ppp[0] = pp2[0]; ppp[1] = pp2[1]; ppp[2] = pp2[2];
						pps.push_back(ppp);
					}
				}

				for (int k1 = 0; k1 < jjs.size(); k1++) {
					vertex_tri(j1, jjs[k1][0], jjs[k1][1], j1s, j2s, j3s);
					for (int k2 = 0; k2 < patchpts[m].size(); k2++) {
						int j1u, j2u, j3u;
						vertex_tri(vertextris[patchpts[m][k2]][0], vertextris[patchpts[m][k2]][1], vertextris[patchpts[m][k2]][2],
							j1u, j2u, j3u);
						if (j1u != j1s || j2u != j2s || j3u != j3s) continue;
						if (edgeflag[iv*vertsz + patchpts[m][k2]] <= 1) break;
						jjs.erase(jjs.begin() + k1);
						pps.erase(pps.begin() + k1);
						k1--;
						break;
					}
				}

				if (jjs.size() == 0) {
					break;
				}

				ret = insert_polyhedron_info(j1, j3, j2, p, jjs, pps, iv);
				if (ret == 0) {
					printf("error : insert_polyhedron_info\n");
					break;
				}
			}

			if (patches[m].back() != patches[m][0]) {
				patches[m] = vpatch;
				continue;
			}

			// merge edge
			if (!merge_edge(vpatch, patches[m], edgetmps[j][0], edgetmps[j][1])) {
				patches[m] = vpatch;
				continue;
			}

			if (vertextris.size() > vertsz) {
				int vertsz1 = vertextris.size();
				unsigned char* edgeflag1 = (unsigned char*)malloc(vertsz1*vertsz1);
				memset(edgeflag1, 0, vertsz1*vertsz1);
				for (k = 0; k < vertsz; k++) {
					memcpy(&edgeflag1[k*vertsz1], &edgeflag[k*vertsz], vertsz);
				}
				free(edgeflag);
				edgeflag = edgeflag1;
				vertsz = vertsz1;
			}

			for (k = 0; k < patches[m].size() - 1; k++) {
				if (edgeflag[patches[m][k] * vertsz + patches[m][k + 1]] > 1) continue;
				edgeflag[patches[m][k] * vertsz + patches[m][k + 1]]++;
				edgeflag[patches[m][k + 1] * vertsz + patches[m][k]]++;
			}

			patches[m] = vpatch;
		}
	}

	for (i = 0; i < patches.size(); i++) {
		for (j = 0; j < patches[i].size(); j++) {
			if (flag[patches[i][j]] > 0) continue;
			tri = vertextris[patches[i][j]];
			for (k = 0; k < 3; k++) {
				if (cirpatchflag[tri[k]] != 0) continue;

				patches.push_back(vector<int>(1, patches[i][j]));
				jcirclePatch.push_back(tri[k]);
				getRow(N_, nn_, n, tri[k]);
				patchdirs.push_back(n);
				cirpatchflag[tri[k]] = -(int)patches.size();
			}
		}
	}

	bool existvolume = true;
	int closedcnt = 0;
	m = 0;
	while (closedcnt < patches.size())
	{
		if (patches[m].size() > 3 &&
			patches[m][0] == patches[m].back()) {
			m++;
			continue;
		}

		if (m == patches.size()) {
			closedcnt = 0;
			for (m = 0; m < patches.size(); m++) {
				if (patches[m].back() != patches[m][0]) break;
				closedcnt++;
			}
			//if (closedcnt == patches.size()) break;
			break;
		}

		if (cirpatchflag[jcirclePatch[m]] <= 0) {
			if (!determine_patch_direction(m))
			{
				cirpatchflag[jcirclePatch[m]] = 0;
				patches.erase(patches.begin() + m);
				jcirclePatch.erase(jcirclePatch.begin() + m);
				patchdirs.erase(patchdirs.begin() + m);
				for (int i1 = m; i1 < jcirclePatch.size(); i1++) {
					if (cirpatchflag[jcirclePatch[i1]] > 0) cirpatchflag[jcirclePatch[i1]]--;
					if (cirpatchflag[jcirclePatch[i1]] < 0) cirpatchflag[jcirclePatch[i1]]++;
				}
				//m++;
				//if (m == patches.size()) break;
				continue;
			}
		}

		vector<vector<int> > vertextri1 = vertextris;
		vector<vector<double> > vertexpt1 = vertexpts;
		int iv, iv1;
		vector<vector<int> > jjs;
		vector<int> jj(2);
		vector<vector<double> > pps;
		vector<double> ppp(3);
		while (patches[m].back() != patches[m][0])
		{
			iv = patches[m].back();
			iv1 = patches[m][patches[m].size() - 2];

			// j1 - principal circle plane
			j1 = jcirclePatch[m];
			// j2 - previous circle plane
			for (int i1 = 0; i1 < 3; i1++) {
				if (vertextris[iv][i1] == j1) continue;
				if (vertextris[iv][i1] == vertextris[iv1][0] ||
					vertextris[iv][i1] == vertextris[iv1][1] ||
					vertextris[iv][i1] == vertextris[iv1][2])
				{
					j2 = vertextris[iv][i1]; break;
				}
			}
			// j3 - new point should be found in the direction (j1,j3)
			for (int i1 = 0; i1 < 3; i1++) {
				if (vertextris[iv][i1] == j1 || vertextris[iv][i1] == j2) continue;
				j3 = vertextris[iv][i1]; break;
			}

			p[0] = vertexpts[iv][0]; p[1] = vertexpts[iv][1]; p[2] = vertexpts[iv][2];

			jjs.clear(); pps.clear();
			int ret = calc_nearest_crossPoint(j1, j3, j2, p, jj1, jj2, pp1, pp2);
			if (jj1 > -1) {
				jj[0] = j3; jj[1] = jj1;
				jjs.push_back(jj);
				ppp[0] = pp1[0]; ppp[1] = pp1[1]; ppp[2] = pp1[2];
				pps.push_back(ppp);
			}
			if (jj2 > -1) {
				jj[0] = j3; jj[1] = jj2;
				jjs.push_back(jj);
				ppp[0] = pp2[0]; ppp[1] = pp2[1]; ppp[2] = pp2[2];
				pps.push_back(ppp);
			}
			ret = calc_nearest_crossPoint(j1, j2, j3, p, jj1, jj2, pp1, pp2);
			int j1t, j2t, j3t, j1s, j2s, j3s;
			vertex_tri(vertextris[iv1][0], vertextris[iv1][1], vertextris[iv1][2], j1t, j2t, j3t);
			if (jj1 > -1) {
				vertex_tri(j1, j2, jj1, j1s, j2s, j3s);
				if (j1t != j1s || j2t != j2s || j3t != j3s) {
					jj[0] = j2; jj[1] = jj1;
					jjs.push_back(jj);
					ppp[0] = pp1[0]; ppp[1] = pp1[1]; ppp[2] = pp1[2];
					pps.push_back(ppp);
				}
			}
			if (jj2 > -1) {
				vertex_tri(j1, j2, jj2, j1s, j2s, j3s);
				if (j1t != j1s || j2t != j2s || j3t != j3s) {
					jj[0] = j2; jj[1] = jj2;
					jjs.push_back(jj);
					ppp[0] = pp2[0]; ppp[1] = pp2[1]; ppp[2] = pp2[2];
					pps.push_back(ppp);
				}
			}

			if (jjs.size() == 0) {
				existvolume = false; break;
			}

			ret = insert_polyhedron_info(j1, j3, j2, p, jjs, pps, iv);
			if (ret == 0) {
				printf("error : insert_polyhedron_info\n");
				existvolume = false; break;
			}

		}

		if (!check_positive_polygon(patches[m], patchdirs[m]))
		{
			cirpatchflag[jcirclePatch[m]] = 0;
			patches.erase(patches.begin() + m);
			jcirclePatch.erase(jcirclePatch.begin() + m);
			patchdirs.erase(patchdirs.begin() + m);
			for (int i1 = m; i1 < jcirclePatch.size(); i1++) {
				if (cirpatchflag[jcirclePatch[i1]] > 0) cirpatchflag[jcirclePatch[i1]]--;
				if (cirpatchflag[jcirclePatch[i1]] < 0) cirpatchflag[jcirclePatch[i1]]++;
			}
			vertextris = vertextri1;
			vertexpts = vertexpt1;
			continue;
		}

		if (!existvolume) break;

		m++;
	}

	if (!existvolume || closedcnt < patches.size()) {
		return false;
	}

	return true;
}

int get_closed_patch(int ipat, unsigned char* edgeflag = NULL, 
	vector<int> pointset = vector<int>(0), bool vertexExpand = true)
{
	int j1, j2, j3, jj1, jj2;
	double p[3], pp1[3], pp2[3], *pp;

	int iv, iv1;
	vector<vector<int> > jjs;
	vector<int> jj(2);
	vector<vector<double> > pps;
	vector<double> ppp(3);
	while (patches[ipat].back() != patches[ipat][0])
	{
		iv = patches[ipat].back();
		iv1 = patches[ipat][patches[ipat].size() - 2];

		// j1 - principal circle plane
		j1 = jcirclePatch[ipat];
		// j2 - previous circle plane
		for (int i1 = 0; i1 < 3; i1++) {
			if (vertextris[iv][i1] == j1) continue;
			if (vertextris[iv][i1] == vertextris[iv1][0] ||
				vertextris[iv][i1] == vertextris[iv1][1] ||
				vertextris[iv][i1] == vertextris[iv1][2])
			{
				j2 = vertextris[iv][i1]; break;
			}
		}
		// j3 - new point should be found in the direction (j1,j3)
		for (int i1 = 0; i1 < 3; i1++) {
			if (vertextris[iv][i1] == j1 || vertextris[iv][i1] == j2) continue;
			j3 = vertextris[iv][i1]; break;
		}

		p[0] = vertexpts[iv][0]; p[1] = vertexpts[iv][1]; p[2] = vertexpts[iv][2];

		jjs.clear(); pps.clear();
		int ret = calc_nearest_crossPoint(j1, j3, j2, p, jj1, jj2, pp1, pp2);
		if (jj1 > -1) {
			jj[0] = j3; jj[1] = jj1;
			jjs.push_back(jj);
			ppp[0] = pp1[0]; ppp[1] = pp1[1]; ppp[2] = pp1[2];
			pps.push_back(ppp);
		}
		if (jj2 > -1) {
			jj[0] = j3; jj[1] = jj2;
			jjs.push_back(jj);
			ppp[0] = pp2[0]; ppp[1] = pp2[1]; ppp[2] = pp2[2];
			pps.push_back(ppp);
		}
		ret = calc_nearest_crossPoint(j1, j2, j3, p, jj1, jj2, pp1, pp2);
		int j1t, j2t, j3t, j1s, j2s, j3s;
		vertex_tri(vertextris[iv1][0], vertextris[iv1][1], vertextris[iv1][2], j1t, j2t, j3t);
		if (jj1 > -1) {
			vertex_tri(j1, j2, jj1, j1s, j2s, j3s);
			if (j1t != j1s || j2t != j2s || j3t != j3s) {
				jj[0] = j2; jj[1] = jj1;
				jjs.push_back(jj);
				ppp[0] = pp1[0]; ppp[1] = pp1[1]; ppp[2] = pp1[2];
				pps.push_back(ppp);
			}
		}
		if (jj2 > -1) {
			vertex_tri(j1, j2, jj2, j1s, j2s, j3s);
			if (j1t != j1s || j2t != j2s || j3t != j3s) {
				jj[0] = j2; jj[1] = jj2;
				jjs.push_back(jj);
				ppp[0] = pp2[0]; ppp[1] = pp2[1]; ppp[2] = pp2[2];
				pps.push_back(ppp);
			}
		}

		if (jjs.size() == 0) return 1;

		ret = insert_polyhedron_info(j1, j3, j2, p, jjs, pps, iv);
		if (ret == 0) {
			printf("error : insert_polyhedron_info\n");
			return 0;
		}

		// 
		if (edgeflag) {
			int vertsz = vertextris.size();
			vector<int> vert = vertextris[patches[ipat].back()];
			vertex_tri(vert[0], vert[1], vert[2], j1s, j2s, j3s);
			bool isbelong = false;
			for (int k2 = 0; k2 < pointset.size(); k2++) {
				int j1u, j2u, j3u;
				vertex_tri(vertextris[pointset[k2]][0], vertextris[pointset[k2]][1], vertextris[pointset[k2]][2],
					j1u, j2u, j3u);
				if (j1u != j1s || j2u != j2s || j3u != j3s) continue;
				isbelong = true;
				if (edgeflag[iv*vertsz + pointset[k2]] <= 1) break;

				return 0;
			}
			if (!vertexExpand && !isbelong) return 0;
		}

	}

	return 2;
}

bool plug_holes(unsigned char*& edgeflag)
{
	double p[3], pp1[3], pp2[3];
	int j1, j2, j3, jj1, jj2;
	int i, j, k, m, ipat;
	int iv, iv1;
	vector<vector<int> > patchpts(patches.size());
	vector<int> vpatch;

	int vertsz = vertextris.size();

	for (i = 0; i < vertextris.size(); i++) {
		for (k = 0; k < 3; k++) {
			j = vertextris[i][k];
			if (cirpatchflag[j] == 0) continue;
			ipat = abs(cirpatchflag[j]) - 1;
			patchpts[ipat].push_back(i);
		}
	}

	for (i = 0; i < patchpts.size(); i++) {
		if (patchpts[i].size() == patches[i].size() - 1) continue;

		vector<vector<int> > edgetmps;
		for (j = 0; j < patches[i].size() - 1; j++) {
			if (edgeflag[patches[i][j] * vertsz + patches[i][j + 1]] != 1) continue;
			vector<int> tmp(2);
			tmp[0] = patches[i][j]; tmp[1] = patches[i][j + 1];
			edgetmps.push_back(tmp);
		}

		for (j = 0; j < edgetmps.size(); j++) {
			if (patchpts[i].size() == patches[i].size() - 1) break;

			if (edgeflag[edgetmps[j][0] * vertsz + edgetmps[j][1]] != 1) continue;

			iv = edgetmps[j][0]; iv1 = edgetmps[j][1];
			j1 = jcirclePatch[i];
			j2 = -1;
			for (int i1 = 0; i1 < 3; i1++) {
				if (vertextris[iv][i1] == j1) continue;
				if (vertextris[iv][i1] == vertextris[iv1][0] ||
					vertextris[iv][i1] == vertextris[iv1][1] ||
					vertextris[iv][i1] == vertextris[iv1][2])
				{
					j2 = vertextris[iv][i1]; break;
				}
			}
			if (j2 == -1) continue;
			m = i;

			polyhedron_info_real2temp(true);
	
			vpatch = patches[m];

			patches[m].clear();
			patches[m].push_back(edgetmps[j][1]);
			patches[m].push_back(edgetmps[j][0]);

			int ret = get_closed_patch(m, edgeflag, patchpts[m], false);

			//vector<vector<int> > jjs;
			//vector<int> jj(2);
			//vector<vector<double> > pps;
			//vector<double> ppp(3);
			//while (patches[m].back() != patches[m][0]) {
			//	iv = patches[m].back();
			//	iv1 = patches[m][patches[m].size() - 2];
			//	// j1 - principal circle plane
			//	j1 = jcirclePatch[m];
			//	// j2 - previous circle plane
			//	for (int i1 = 0; i1 < 3; i1++) {
			//		if (vertextris[iv][i1] == j1) continue;
			//		if (vertextris[iv][i1] == vertextris[iv1][0] ||
			//			vertextris[iv][i1] == vertextris[iv1][1] ||
			//			vertextris[iv][i1] == vertextris[iv1][2])
			//		{
			//			j2 = vertextris[iv][i1]; break;
			//		}
			//	}
			//	// j3 - new point should be found in the direction (j1,j3)
			//	for (int i1 = 0; i1 < 3; i1++) {
			//		if (vertextris[iv][i1] == j1 || vertextris[iv][i1] == j2) continue;
			//		j3 = vertextris[iv][i1]; break;
			//	}
			//	p[0] = vertexpts[iv][0]; p[1] = vertexpts[iv][1]; p[2] = vertexpts[iv][2];
			//
			//	jjs.clear(); pps.clear();
			//	int ret = calc_nearest_crossPoint(j1, j3, j2, p, jj1, jj2, pp1, pp2);
			//	if (jj1 > -1) {
			//		jj[0] = j3; jj[1] = jj1;
			//		jjs.push_back(jj);
			//		ppp[0] = pp1[0]; ppp[1] = pp1[1]; ppp[2] = pp1[2];
			//		pps.push_back(ppp);
			//	}
			//	if (jj2 > -1) {
			//		jj[0] = j3; jj[1] = jj2;
			//		jjs.push_back(jj);
			//		ppp[0] = pp2[0]; ppp[1] = pp2[1]; ppp[2] = pp2[2];
			//		pps.push_back(ppp);
			//	}
			//	ret = calc_nearest_crossPoint(j1, j2, j3, p, jj1, jj2, pp1, pp2);
			//	int j1t, j2t, j3t, j1s, j2s, j3s;
			//	vertex_tri(vertextris[iv1][0], vertextris[iv1][1], vertextris[iv1][2], j1t, j2t, j3t);
			//	if (jj1 > -1) {
			//		vertex_tri(j1, j2, jj1, j1s, j2s, j3s);
			//		if (j1t != j1s || j2t != j2s || j3t != j3s) {
			//			jj[0] = j2; jj[1] = jj1;
			//			jjs.push_back(jj);
			//			ppp[0] = pp1[0]; ppp[1] = pp1[1]; ppp[2] = pp1[2];
			//			pps.push_back(ppp);
			//		}
			//	}
			//	if (jj2 > -1) {
			//		vertex_tri(j1, j2, jj2, j1s, j2s, j3s);
			//		if (j1t != j1s || j2t != j2s || j3t != j3s) {
			//			jj[0] = j2; jj[1] = jj2;
			//			jjs.push_back(jj);
			//			ppp[0] = pp2[0]; ppp[1] = pp2[1]; ppp[2] = pp2[2];
			//			pps.push_back(ppp);
			//		}
			//	}

			//	for (int k1 = 0; k1 < jjs.size(); k1++) {
			//		vertex_tri(j1, jjs[k1][0], jjs[k1][1], j1s, j2s, j3s);
			//		bool isbelong = false;
			//		for (int k2 = 0; k2 < patchpts[m].size(); k2++) {
			//			int j1u, j2u, j3u;
			//			vertex_tri(vertextris[patchpts[m][k2]][0], vertextris[patchpts[m][k2]][1], vertextris[patchpts[m][k2]][2],
			//				j1u, j2u, j3u);
			//			if (j1u != j1s || j2u != j2s || j3u != j3s) continue;
			//			isbelong = true;
			//			if (edgeflag[iv*vertsz + patchpts[m][k2]] <= 1) break;
			//			jjs.erase(jjs.begin() + k1);
			//			pps.erase(pps.begin() + k1);
			//			k1--;
			//			break;
			//		}
			//		if (!isbelong) {
			//			jjs.erase(jjs.begin() + k1);
			//			pps.erase(pps.begin() + k1);
			//			k1--;
			//		}
			//	}

			//	if (jjs.size() == 0) {
			//		break;
			//	}
			//	ret = insert_polyhedron_info(j1, j3, j2, p, jjs, pps, iv);
			//	if (ret == 0) {
			//		printf("error : insert_polyhedron_info\n"); 
			//		break;
			//	}
			//}

			if (ret < 2 || patches[m].back() != patches[m][0])
			{
				polyhedron_info_real2temp(false);
				continue;
			}

			// if direction of the gotten polygon is against the patch direction
			if (!check_positive_polygon(patches[m], patchdirs[m]))
			{
				polyhedron_info_real2temp(false);
				continue;
			}

			// merge edge
			if (!merge_edge(vpatch, patches[m], edgetmps[j][0], edgetmps[j][1])) {
				polyhedron_info_real2temp(false);
				continue;
			}

			//if (vertextris.size() > vertsz) {
			//	int vertsz1 = vertextris.size();
			//	unsigned char* edgeflag1 = (unsigned char*)malloc(vertsz1*vertsz1);
			//	memset(edgeflag1, 0, vertsz1*vertsz1);
			//	for (k = 0; k < vertsz; k++) {
			//		memcpy(&edgeflag1[k*vertsz1], &edgeflag[k*vertsz], vertsz);
			//	}
			//	free(edgeflag);
			//	edgeflag = edgeflag1;
			//	vertsz = vertsz1;
			//}

			for (k = 0; k < patches[m].size() - 1; k++) {
				if (edgeflag[patches[m][k] * vertsz + patches[m][k + 1]] > 1) continue;
				edgeflag[patches[m][k] * vertsz + patches[m][k + 1]]++;
				edgeflag[patches[m][k + 1] * vertsz + patches[m][k]]++;
				if (edgeflag[patches[m][k] * vertsz + patches[m][k + 1]] != 1) continue;
				vector<int> tmp(2);
				tmp[0] = patches[m][k]; tmp[1] = patches[m][k + 1];
				edgetmps.push_back(tmp);
			}

			patches[m] = vpatch;
		}
	}

	return true;
}

void remove_patch(int ipatch)
{
	cirpatchflag[jcirclePatch[ipatch]] = 0;

	for (int i = 0; i < patches[ipatch].size(); i++) {
		vertlapcnt[patches[ipatch][i]]--;
	}

	patches.erase(patches.begin() + ipatch);
	jcirclePatch.erase(jcirclePatch.begin() + ipatch);
	patchdirs.erase(patchdirs.begin() + ipatch);
	for (int i = ipatch; i < jcirclePatch.size(); i++) {
		if (cirpatchflag[jcirclePatch[i]] > 0) cirpatchflag[jcirclePatch[i]]--;
		if (cirpatchflag[jcirclePatch[i]] < 0) cirpatchflag[jcirclePatch[i]]++;
	}
}

double calculate_one_volume()
{
	int i, j;
	double vol = 0;
	double u[3], uu[3], pt1[3], pt2[3];
	vector<double> pt0, pt00;
	vector<int> patch;

	for (i = 0; i < patches.size(); i++) {
		patch = patches[i];
		if (patch.size() < 4) continue;

		while (1) {
			bool bconvex1, bconvex2;
			pt0 = vertexpts[patch[0]];
			minus_(pt0, vertexpts[patch[patch.size() - 2]], pt1);
			minus_(vertexpts[patch[1]], pt0, pt2);
			cross(pt1, pt2, u);
			bconvex1 = dot(u, patchdirs[i]) > 0;
			u[0] = -u[0]; u[1] = -u[1]; u[2] = -u[2];

			bool bdelete = false;
			for (j = 1; j < patch.size() - 1; j++) {
				pt00 = vertexpts[patch[j]];
				minus_(pt00, vertexpts[patch[j - 1]], pt1);
				minus_(vertexpts[patch[j + 1]], pt00, pt2);
				cross(pt1, pt2, uu);
				bconvex2 = dot(uu, patchdirs[i]) > 0;
				uu[0] = -uu[0]; uu[1] = -uu[1]; uu[2] = -uu[2];

				if (bconvex1 == bconvex2) {
					u[0] = uu[0]; u[1] = uu[1]; u[2] = uu[2];
					pt0 = pt00;
					continue;
				}

				if (bconvex1) {
					vol += dot(u, pt0);
					patch.erase(patch.begin() + j - 1);
					if (j == 1) {
						patch.erase(patch.end() - 1);
						patch.push_back(patch[0]);
					}
				}
				if (bconvex2) {
					vol += dot(uu, pt00);
					patch.erase(patch.begin() + j);
				}
				bdelete = true;
				break;
			}

			if (!bdelete) break;
		}

		while (patch.size() > 3) {
			for (j = 0; j < patch.size() - 2; j++) {
				pt0 = vertexpts[patch[j + 1]];
				minus_(pt0, vertexpts[patch[j]], pt1);
				minus_(vertexpts[patch[j + 2]], pt0, pt2);
				cross(pt1, pt2, u);
				u[0] = -u[0]; u[1] = -u[1]; u[2] = -u[2];
				vol += dot(u, pt0);
				patch.erase(patch.begin() + j + 1);
			}
		}
	}

	return (vol / 6);
}

void calculate_volumes(double* C_, double* R, double* N_)
{
	//if (numTriple_ < 4) return;

	int i, k, m;
	int j1, j2, j3, jj1, jj2;
	double p[3], pp1[3], pp2[3], *pp;
	double u1[3], u2[3], u3[3];
	double v1[3], v2[3], v3[3];
	double c1[3], c2[3], c3[3];
	double n1[3], n2[3], n3[3];
	stCrossInfo crsspt;
	int trip[3];
	vector<int> vtrip(3);
	vector<double> vpt(3);

	countVolume_ = 0;

	thrpt[0] = -1.0e13; thrpt[1] = -1.0e13; thrpt[2] = -1.0e13;

	// open triple data to read
	ifstream fcrosspts;
	sprintf(fnamepts, "%s", "crosspts.dat");
	fcrosspts.open(fnamepts, ios_base::in | ios_base::binary);
	if (fcrosspts.bad()) {
		printf("The file %s was not opened.\n", fnamepts);
		return;
	}

	// open triple data to read
	ifstream ftriples;
	ftriples.open(fnametri, ios_base::in | ios_base::binary);
	if (ftriples.bad()) {
		printf("The file %s was not opened.\n", fnametri);
		return;
	}
	ftriples.seekg(0, ios_base::end);
	numTriple_ = ftriples.tellg() / (sizeof(int) * 3);
	ftriples.seekg(0, ios_base::beg);

	int cnttrip = numTriple_ < TRIPLE_NUM_UNIT ? (int)numTriple_ : TRIPLE_NUM_UNIT;
	int cntgroup = (int)((numTriple_ - 1) / cnttrip + 1);
	int cntsub = TRIPLE_NUM_UNIT / cntgroup;
	
#define SET_THRESH_POINT(p)				\
	thrpt[0] = (float)p[0];				\
	thrpt[1] = (float)p[1];				\
	thrpt[2] = (float)p[2];


	for (i = 0; i < cntgroup; i++)
	{
		int cnttrip_r = cnttrip;
		if (i == numTriple_ / cnttrip) {
			cnttrip_r = numTriple_ % cnttrip;
			if (cnttrip_r == 0) break;
		}
	
		fcrosspts.read((char*)crossInfo, sizeof(stCrossInfo) * cnttrip_r);
		for (k = 0; k < cnttrip_r; k++) {
			crsspt = crossInfo[k];
			
			ftriples.seekg(crsspt.nIdx * sizeof(int) * 3);
			ftriples.read((char*)trip, sizeof(int) * 3);

			j1 = trip[0]; j2 = trip[1]; j3 = trip[2];

			if (!isnecessarypoint(j1, j2, j3)) continue;
			getRow(N_, nn_, n1, j1); getRow(C_, nn_, c1, j1);
			getRow(N_, nn_, n2, j2); getRow(C_, nn_, c2, j2);
			getRow(N_, nn_, n3, j3); getRow(C_, nn_, c3, j3);
			if (!intersectionPoint(p, c1, R[j1], n1, c2, R[j2], n2, c3, R[j3], n3)) continue;

			cirpatchflag = vector<int>(nn_, 0);
			// insert initial vertices
			if (!initial_vertices(j1, j2, j3, p)) {
				SET_THRESH_POINT(p);
				continue;
			}

			initial_arrange_patches();
			if (patches.size() <= 3)
				i = i;

			bool existvolume = true;
			int closedcnt = 0;
			m = 0;
			while (closedcnt < patches.size())
			{
				if (m == patches.size()) {
					closedcnt = 0;
					for (m = 0; m < patches.size(); m++) {
						if (patches[m].back() != patches[m][0]) break;
						closedcnt++;
					}
					break;
				}
				
				// if patch direction is not determined
				if (cirpatchflag[jcirclePatch[m]] <= 0) {
					if (!determine_patch_direction(m)){
						remove_patch(m);
						continue;
					}
				}

				polyhedron_info_real2temp(true);

				int ret = get_closed_patch(m);
				
				if (ret == 0) {
					existvolume = false;
					break;
				}

				if (ret == 1) {
					polyhedron_info_real2temp(false);
					remove_patch(m);
					continue;
				}

				if (m < 3 && patches[m][0] != patchestmp[m][0]) {
					existvolume = false;
					break;
				}

				// if direction of the gotten polygon is against the patch direction
				if (!check_positive_polygon(patches[m], patchdirs[m]))
				{
					if (m < 3) {
						existvolume = false;
						break;
					}

					polyhedron_info_real2temp(false);
					remove_patch(m);
					continue;
				}

				m++;
			}

			if (!existvolume || closedcnt < patches.size()) {
				SET_THRESH_POINT(vertexpts[0]);
				continue;
			}
			if (patches.size() <= 3) {
				SET_THRESH_POINT(vertexpts[0]);
				continue;
			}

			//remove_middle_point_of_patch();

			unsigned char* edgeflag = (unsigned char*)malloc(vertextris.size()*vertextris.size());

			bool bclosed = check_closed_polyhedron(edgeflag);
			while (!bclosed) {
				plug_holes(edgeflag);

				//remove_middle_point_of_patch();
			
				bclosed = check_closed_polyhedron(edgeflag);
				if (bclosed) break;

				if (!add_patches(edgeflag))
					break;
				//remove_middle_point_of_patch();
				
				bclosed = check_closed_polyhedron(edgeflag);
				break; // for debug
			}
			if (!bclosed) {
				SET_THRESH_POINT(vertexpts[0]);
				continue;
			}

			free(edgeflag);

			remove_middle_point_of_patch();

			double vol = calculate_one_volume();
			printf("%.4f\n", vol);
			countVolume_++;

			SET_THRESH_POINT(vertexpts[0]);
		}
	}

error_pos:

	ftriples.close();
	fcrosspts.close();

	printf("Volume count = %d\n", countVolume_);

	return;
}

