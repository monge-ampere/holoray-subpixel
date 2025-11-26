#include "stdafx.h"
#include "DataAnalysis.h"
#include "FileProc.h"

typedef double double3[3];

using namespace std;

int main()
{
	printf("beginning----------------------------\n\n");

	Pt *avail_bps = new Pt[MAX_VERTEX_NUM];
	vector<int> avail_idx;
	vector<vector<int>> avail_tris;

	Pt *outer_bps = new Pt[MAX_VERTEX_NUM];
	vector<int> outer_idx;
	vector<vector<int>> outer_tris;

	ImgProc *proc = new ImgProc("data//triangle//mask.bmp"
		, "data//triangle//edge.bmp", "data//triangle//triangle.bmp");

	DataAnalysis analysis(proc);

	readm(avail_bps, avail_idx, avail_tris, "E:\\testPic\\triangle\\inner_integrogram.m");
	readm(outer_bps, outer_idx, outer_tris, "E:\\testPic\\triangle\\outer_integrogram.m");
	printf("reading m is finished\n");
	analysis.GetInnerData(avail_bps, avail_idx, avail_tris);
	analysis.GetOuterData(outer_bps, outer_idx, outer_tris);
	printf("getting data is finished\n");
	double bdy_area = analysis.IdentBdy();
	printf("identifying boundary is finished\n");
	//double bdy_rate = analysis.GetBdySubPxPos();
	//double area = inner_area + bdy_area*bdy_rate;

	analysis.GetBdySubPxPos1("E:\\testPic\\subpixel_edge_show4holoRay_leaf.bmp");
	printf("getting boundary's subpixel position is finished\n");
	analysis.OutputIntegralImg("integrogram.bmp");

	//printf("outputing----------------------------\n\n");
	//printf("sub-pixel area is %f\n", area);

	printf("finishing----------------------------\n\n");
	system("pause");
	delete proc;
	return 0;
}