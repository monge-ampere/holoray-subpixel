#pragma once
/*!
*      \class ImgProc
*      \brief image analysis
*	   \author Wei Feng
*      \date 11/21/2024
*/
#include <opencv2/opencv.hpp>
#include "Algo.h"
#include "KdTree\kdtree.h"
#define MAX_VERTEX_NUM 1048576

typedef double double2[2];
typedef int int3[3];

enum TOPOTYPE
{
	INNER, OUTER
};

class ImgProc
{
	friend class DataAnalysis;
public:
	ImgProc(const char *maskPath, const char *bdyPath, const char *oriPath);
	virtual ~ImgProc();
	void InitKdTree4P();
	void DelKdTree4P();
private:
	Mat ori_mat;
	Mat src_mat;
	Mat bdy_mat;
	Mat sobel_mat;
	Mat canny_mat;
	Mat mask_mat;
	int r;
	int c;
	vector<Point2d> ps;
	vector<uchar> gs;
	vector<int> boundary;
	vector<Point2d> orderedBdy;

	// kd tree of model points
	kdtree::KDTree*     M_tree4P;
	kdtree::KDTreeArray M_data4P;

	double EuclidDist(double2 &v1, double2 &v2);
	void Sobel4Edge(Mat &sMat, Mat &srcMat);
	void DoG4Edge(Mat &dMat, Mat &srcMat);
	void GetBdy();
	bool GetFirstBdyPt(Mat &mat);
};