#pragma once
/*!
*      \class DataAnalysis
*      \brief calculate sub-pixel boundary
*	   \author Wei Feng
*      \date 01/10/2025
*/
#include "ImgProc.h"
#include "FileProc.h"
#define MAXLINESIZE 1024
#define MAGNIFIC 11
#define REGUTRIAREA 0.3031088913

class DataAnalysis
{
public:
	DataAnalysis(ImgProc *proc);
	virtual ~DataAnalysis();
	void GetInnerData(Pt *ps, vector<int> &idx, vector<vector<int>> &triSet);
	void GetOuterData(Pt *ps, vector<int> &idx, vector<vector<int>> &triSet);
	double IdentBdy();
	uchar GetGrayV(Pt &pos, TOPOTYPE type);
	uchar GetGrayV(vector<Pt> &neigh_uvs, Pt &pos);
	double GetBdySubPxPos();
	void GetBdySubPxPos(vector<Pt> &subpxBdy);
	void GetBdySubPxPos1(const char* filePath);
	void OutputIntegralImg(const char *filePath);
private:
	ImgProc *proc;
	vector<Pt> innerPs;
	vector<vector<Pt>> innerTriSet;
	vector<int> innerBdy;

	vector<Pt> outerPs;
	vector<vector<Pt>> outerTriSet;
	vector<int> outerBdy;

	double innerMinY;
	double innerMaxY;
	double outerMinY;
	double outerMaxY;

	double innerMin;
	double innerMax;
	double outerMin;
	double outerMax;
	vector<double> bdyWeigh;
	vector<double> bdyWeighOut;

	double offset4out;

	double matrix[9];
	double anti_matrix[9];
	double y[3];
	double result[3];

	kdtree::KDTree*     inner_ori_tree;
	kdtree::KDTreeArray inner_ori_data;

	kdtree::KDTree*     inner_tri_tree;
	kdtree::KDTreeArray inner_tri_data;

	kdtree::KDTree*     outer_ori_tree;
	kdtree::KDTreeArray outer_ori_data;

	kdtree::KDTree*     outer_tri_tree;
	kdtree::KDTreeArray outer_tri_data;

	int Search(vector<Pt> &neigh_uvs, Pt &uv, TOPOTYPE type);
	void Search(Pt &uv, TOPOTYPE type);
	double P2lineDist(Pt &p, Pt &a, Pt &b);
	double P2lineDist4UV(Pt &p, Pt &a, Pt &b);
	double GetSubpixelPoint(double2 *point_set);
	void GetTriSet(vector<vector<Pt>> &triSet
		, Pt *ps, vector<vector<int>> &tridxSet
		, double deviThres, double maxY, double minY);
	void CalcWeigh();
	double CalcWeigh(Pt &pt, TOPOTYPE type);
	double CalcWeigh(Pt &pt, TOPOTYPE type, int iterCount);
	void FineWeigh(double *wDist, TOPOTYPE type, double thres, int maxIterNum);
	void GetBdbox(vector<Pt> &neigh_uvs, Pt &uv);
	void DelOffset(Pt &p1, Pt &p2, Pt &p3);
	double f(double x, Pt &p, TOPOTYPE type);
	double solve_fx_equals_target(Pt &p, TOPOTYPE type
		, double target, double x0, double step);
};