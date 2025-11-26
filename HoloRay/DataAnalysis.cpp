#include "stdafx.h"
#include "DataAnalysis.h"


DataAnalysis::DataAnalysis(ImgProc *proc)
{
	this->proc = proc;
}


DataAnalysis::~DataAnalysis()
{
}

void DataAnalysis::GetInnerData(Pt *ps, vector<int> &idx, vector<vector<int>> &triSet)
{
	innerMin = MAXV;
	innerMax = -1 * MAXV;
	innerMinY = MAXV;
	innerMaxY = -1 * MAXV;
	int num = idx.size();
	int id;
	Pt pt;

	inner_ori_data.resize(boost::extents[num][2]);
	innerPs.clear();
	for (int i = 0; i != num; ++i)
	{
		id = idx[i];
		pt.x = ps[id].x;
		pt.y = ps[id].y;
		pt.ori = ps[id].ori;
		pt.oci = ps[id].oci;
		pt.val = proc->sobel_mat.at<uchar>(floor(ps[id].ori), floor(ps[id].oci));
		innerPs.push_back(pt);

		inner_ori_data[i][0] = ps[id].oci;
		inner_ori_data[i][1] = ps[id].ori;

		if (pt.x > innerMax)
		{
			innerMax = pt.x;
		}
		if (pt.x < innerMin)
		{
			innerMin = pt.x;
		}

		if (pt.y > innerMax)
		{
			innerMax = pt.y;
		}
		if (pt.y < innerMin)
		{
			innerMin = pt.y;
		}

		if (pt.y > innerMaxY)
		{
			innerMaxY = pt.y;
		}
		if (pt.y < innerMinY)
		{
			innerMinY = pt.y;
		}
	}
	inner_ori_tree = new kdtree::KDTree(inner_ori_data);

	GetTriSet(innerTriSet, ps, triSet, 0.1, innerMaxY, innerMinY);

	num = innerTriSet.size();
	inner_tri_data.resize(boost::extents[num][2]);
	double ctrX, ctrY;
	for (int i = 0; i != num; ++i)
	{
		ctrX = (innerTriSet[i][0].x
			+ innerTriSet[i][1].x
			+ innerTriSet[i][2].x) / 3;

		ctrY = (innerTriSet[i][0].y
			+ innerTriSet[i][1].y
			+ innerTriSet[i][2].y) / 3;
		inner_tri_data[i][0] = ctrX;
		inner_tri_data[i][1] = ctrY;
	}

	inner_tri_tree = new kdtree::KDTree(inner_tri_data);
}

void DataAnalysis::GetOuterData(Pt *ps, vector<int> &idx, vector<vector<int>> &triSet)
{
	outerMin = MAXV;
	outerMax = -1 * MAXV;
	outerMinY = MAXV;
	outerMaxY = -1 * MAXV;
	int num = idx.size();
	int id;
	Pt pt;

	outer_ori_data.resize(boost::extents[num][2]);
	outerPs.clear();
	for (int i = 0; i != num; ++i)
	{
		id = idx[i];
		pt.x = ps[id].x;
		pt.y = ps[id].y;
		pt.ori = ps[id].ori;
		pt.oci = ps[id].oci;
		pt.val = proc->sobel_mat.at<uchar>(floor(ps[id].ori), floor(ps[id].oci));
		outerPs.push_back(pt);

		outer_ori_data[i][0] = ps[id].oci;
		outer_ori_data[i][1] = ps[id].ori;

		if (pt.x > outerMax)
		{
			outerMax = pt.x;
		}
		if (pt.x < outerMin)
		{
			outerMin = pt.x;
		}

		if (pt.y > outerMax)
		{
			outerMax = pt.y;
		}
		if (pt.y < outerMin)
		{
			outerMin = pt.y;
		}

		if (pt.y > outerMaxY)
		{
			outerMaxY = pt.y;
		}
		if (pt.y < outerMinY)
		{
			outerMinY = pt.y;
		}
	}
	outer_ori_tree = new kdtree::KDTree(outer_ori_data);
	GetTriSet(outerTriSet, ps, triSet, 0.1, outerMaxY, outerMinY);

	num = outerTriSet.size();
	outer_tri_data.resize(boost::extents[num][2]);
	double ctrX, ctrY;
	for (int i = 0; i != num; ++i)
	{
		ctrX = (outerTriSet[i][0].x
			+ outerTriSet[i][1].x
			+ outerTriSet[i][2].x) / 3;

		ctrY = (outerTriSet[i][0].y
			+ outerTriSet[i][1].y
			+ outerTriSet[i][2].y) / 3;
		outer_tri_data[i][0] = ctrX;
		outer_tri_data[i][1] = ctrY;
	}

	outer_tri_tree = new kdtree::KDTree(outer_tri_data);
}

double DataAnalysis::IdentBdy()
{
	int n = proc->orderedBdy.size();
	kdtree::KDTreeResultVector result;
	std::vector<float> query(2);
	for (int i = 0; i != n; ++i)
	{
		query[0] = proc->orderedBdy[i].x;
		query[1] = proc->orderedBdy[i].y;
		inner_ori_tree->n_nearest(query, 1, result);
		if (abs(innerPs[result[0].idx].ori - query[1]) < EPSILON
			&& abs(innerPs[result[0].idx].oci - query[0]) < EPSILON)
		{
			innerBdy.push_back(result[0].idx);
		}
		else
		{
			printf("inner No %d boundary segment is not found\n", i);
		}


		outer_ori_tree->n_nearest(query, 1, result);

		if (abs(outerPs[result[0].idx].ori - query[1]) < EPSILON
			&& abs(outerPs[result[0].idx].oci - query[0]) < EPSILON)
		{
			outerBdy.push_back(result[0].idx);
		}
		else
		{
			printf("outer No %d boundary segment is not found\n", i);
		}
	}

	CalcWeigh();

	//double wDist_i[4096] = { 0 };
	//double wDist_o[4096] = { 0 };

	//FineWeigh(wDist_i, INNER, 0.01, 10);
	//FineWeigh(wDist_o, OUTER, 0.01, 10);
	write2f(proc->orderedBdy, "E:\\order_bdy.csv");
	return n;
}

void DataAnalysis::Search(Pt &uv, TOPOTYPE type)
{
	kdtree::KDTreeResultVector result;
	std::vector<float> query(2);
	query[0] = uv.x;
	query[1] = uv.y;
	vector<Pt> Ps;
	vector<vector<Pt>> tris;
	vector<Pt> mTri;
	vector<Point2d> srcPs, qstPs;
	int num;
	double norm = 0;
	if (INNER == type)
	{
		inner_tri_tree->n_nearest(query, 16, result);
		for (int i = 0; i != 16; ++i)
		{
			vector<Pt> tri = innerTriSet[result[i].idx];
			if (pointInTriangle(tri[0], tri[1], tri[2], uv))
			{
				mTri = tri;
			}
			else
			{
				tris.push_back(tri);
			}
		}
	}
	else if (OUTER == type)
	{
		outer_tri_tree->n_nearest(query, 16, result);
		for (int i = 0; i != 16; ++i)
		{
			vector<Pt> tri = outerTriSet[result[i].idx];
			if (pointInTriangle(tri[0], tri[1], tri[2], uv))
			{
				mTri = tri;
			}
			else
			{
				tris.push_back(tri);
			}
		}
	}

	if (0 == mTri.size())
	{
		//printf("-------------mTri's size is zero-------------\n");
		return;
	}

	Ps = extractUniqueVertices(tris);

	double score = triangleQuality(mTri);
	if (score > 0.3)
	{
		tris.clear();

		tris = findTrianglesContainingP(Ps, uv);

		if (0 == tris.size())
		{
			//printf("-------------tris's size is zero-------------\n");
			return;
		}
		mTri = findBestTriangle(tris);
	}
	num = Ps.size();
	double maxP0 = -1 * MAXV, minP0 = MAXV, maxP = -1 * MAXV, minP = MAXV;
	Point2d averP0, averP;
	averP0.x = 0;
	averP0.y = 0;
	averP.x = 0;
	averP.y = 0;

	for (int i = 0; i != num; ++i)
	{
		Point2d p0, p;
		p0.x = Ps[i].x;
		p0.y = Ps[i].y;

		averP0.x += p0.x / num;
		averP0.y += p0.y / num;

		if (p0.x > maxP0)
		{
			maxP0 = p0.x;
		}
		if (p0.x < minP0)
		{
			minP0 = p0.x;
		}

		if (p0.y > maxP0)
		{
			maxP0 = p0.y;
		}
		if (p0.y < minP0)
		{
			minP0 = p0.y;
		}

		p.x = Ps[i].oci;
		p.y = Ps[i].ori;

		averP.x += p.x / num;
		averP.y += p.y / num;

		if (p.x > maxP)
		{
			maxP = p.x;
		}
		if (p.x < minP)
		{
			minP = p.x;
		}

		if (p.y > maxP)
		{
			maxP = p.y;
		}
		if (p.y < minP)
		{
			minP = p.y;
		}

		srcPs.push_back(p0);
		qstPs.push_back(p);
	}

	for (int i = 0; i != num; ++i)
	{
		srcPs[i].x = (srcPs[i].x - averP0.x) / (maxP0 - minP0);
		srcPs[i].y = (srcPs[i].y - averP0.y) / (maxP0 - minP0);

		qstPs[i].x = (qstPs[i].x - averP.x) / (maxP - minP);
		qstPs[i].y = (qstPs[i].y - averP.y) / (maxP - minP);
	}

	for (int i = 0; i != 3; ++i)
	{
		mTri[i].x = (mTri[i].x - averP0.x) / (maxP0 - minP0);
		mTri[i].y = (mTri[i].y - averP0.y) / (maxP0 - minP0);
	}

	Point2d A, B, A_, B_;
	A.x = srcPs[0].x;
	A.y = srcPs[0].y;
	B.x = srcPs[1].x;
	B.y = srcPs[1].y;

	A_.x = qstPs[0].x;
	A_.y = qstPs[0].y;
	B_.x = qstPs[1].x;
	B_.y = qstPs[1].y;
	double angle = estimateTriangleRotationRadian(A, B, A_, B_);
	cv::Matx22d Rt(cos(angle), -1 * sin(angle), sin(angle), cos(angle));
	cv::Point2d t1(0, 0);

	double a1, a2, tranX, tranY;
	Point2d current_pt, uv_copy;
	current_pt.x = (uv.x - averP0.x) / (maxP0 - minP0);
	current_pt.y = (uv.y - averP0.y) / (maxP0 - minP0);
	Point2d res;
	for (int i = 0; i != num; ++i)
	{
		res = applyAffineTransform(Rt, t1, srcPs[i]);
		srcPs[i] = res;
	}
	res = applyAffineTransform(Rt, t1, current_pt);
	current_pt = res;

	for (int i = 0; i != 3; ++i)
	{
		Point2d pt;
		pt.x = mTri[i].x;
		pt.y = mTri[i].y;
		res = applyAffineTransform(Rt, t1, pt);
		mTri[i].x = res.x;
		mTri[i].y = res.y;
	}

	GetTransformParams(a1, a2, tranX, tranY
		, srcPs, qstPs, current_pt);

	cv::Matx22d R(a1, -1 * a2, a2, a1);
	cv::Point2d t(tranX, tranY);

	cv::Point2d P_prime = applyAffineTransform(R, t, current_pt);
	uv_copy.x = P_prime.x*(maxP - minP) + averP.x;
	uv_copy.y = P_prime.y*(maxP - minP) + averP.y;

	current_pt.x = mTri[0].x;
	current_pt.y = mTri[0].y;

	P_prime = applyAffineTransform(R, t, current_pt);

	mTri[0].x = P_prime.x*(maxP - minP) + averP.x;
	mTri[0].y = P_prime.y*(maxP - minP) + averP.y;

	current_pt.x = mTri[1].x;
	current_pt.y = mTri[1].y;

	P_prime = applyAffineTransform(R, t, current_pt);

	mTri[1].x = P_prime.x*(maxP - minP) + averP.x;
	mTri[1].y = P_prime.y*(maxP - minP) + averP.y;

	current_pt.x = mTri[2].x;
	current_pt.y = mTri[2].y;

	P_prime = applyAffineTransform(R, t, current_pt);

	mTri[2].x = P_prime.x*(maxP - minP) + averP.x;
	mTri[2].y = P_prime.y*(maxP - minP) + averP.y;
	srcPs.clear();
	qstPs.clear();
	for (int i = 0; i != 3; ++i)
	{
		Point2d pt;
		pt.x = mTri[i].x;
		pt.y = mTri[i].y;
		srcPs.push_back(pt);

		pt.x = mTri[i].oci;
		pt.y = mTri[i].ori;
		qstPs.push_back(pt);
	}

	Point2d qstP = findCorrespondingPoint(uv_copy, srcPs, qstPs);
	uv.oci = qstP.x;
	uv.ori = qstP.y;
}

int DataAnalysis::Search(vector<Pt> &neigh_uvs, Pt &uv, TOPOTYPE type)
{
	kdtree::KDTreeResultVector result;
	vector<float> query(2);
	query[0] = uv.x;
	query[1] = uv.y;
	vector<Pt> Ps;
	vector<vector<Pt>> tris;
	vector<Pt> mTri;
	vector<Point2d> srcPs, qstPs;
	int num;
	double norm = 0;
	if (INNER == type)
	{
		inner_tri_tree->n_nearest(query, 16, result);
		for (int i = 0; i != 16; ++i)
		{
			vector<Pt> tri = innerTriSet[result[i].idx];
			if (pointInTriangle(tri[0], tri[1], tri[2], uv))
			{
				mTri = tri;
			}
			else
			{
				tris.push_back(tri);
			}
		}
	}
	else if (OUTER == type)
	{
		outer_tri_tree->n_nearest(query, 16, result);
		for (int i = 0; i != 16; ++i)
		{
			vector<Pt> tri = outerTriSet[result[i].idx];
			if (pointInTriangle(tri[0], tri[1], tri[2], uv))
			{
				mTri = tri;
			}
			else
			{
				tris.push_back(tri);
			}
		}
	}

	if (0 == mTri.size())
	{
		printf("-------------mTri's size is zero-------------\n");
		return 0;
	}

	Ps = extractUniqueVertices(tris);

	double score = triangleQuality(mTri);
	if (score > 1.0)
	{
		tris.clear();

		tris = findTrianglesContainingP(Ps, uv);

		if (0 != tris.size())
		{
			mTri = findBestTriangle(tris);
		}
	}
	num = Ps.size();
	double maxP0 = -1 * MAXV, minP0 = MAXV, maxP = -1 * MAXV, minP = MAXV;
	Point2d averP0, averP;
	averP0.x = 0;
	averP0.y = 0;
	averP.x = 0;
	averP.y = 0;

	for (int i = 0; i != num; ++i)
	{
		Point2d p0, p;
		p0.x = Ps[i].x;
		p0.y = Ps[i].y;

		averP0.x += p0.x / num;
		averP0.y += p0.y / num;

		if (p0.x > maxP0)
		{
			maxP0 = p0.x;
		}
		if (p0.x < minP0)
		{
			minP0 = p0.x;
		}

		if (p0.y > maxP0)
		{
			maxP0 = p0.y;
		}
		if (p0.y < minP0)
		{
			minP0 = p0.y;
		}

		p.x = Ps[i].oci;
		p.y = Ps[i].ori;

		averP.x += p.x / num;
		averP.y += p.y / num;

		if (p.x > maxP)
		{
			maxP = p.x;
		}
		if (p.x < minP)
		{
			minP = p.x;
		}

		if (p.y > maxP)
		{
			maxP = p.y;
		}
		if (p.y < minP)
		{
			minP = p.y;
		}

		srcPs.push_back(p0);
		qstPs.push_back(p);
	}

	for (int i = 0; i != num; ++i)
	{
		srcPs[i].x = (srcPs[i].x - averP0.x) / (maxP0 - minP0);
		srcPs[i].y = (srcPs[i].y - averP0.y) / (maxP0 - minP0);

		qstPs[i].x = (qstPs[i].x - averP.x) / (maxP - minP);
		qstPs[i].y = (qstPs[i].y - averP.y) / (maxP - minP);
	}

	for (int i = 0; i != 3; ++i)
	{
		mTri[i].x = (mTri[i].x - averP0.x) / (maxP0 - minP0);
		mTri[i].y = (mTri[i].y - averP0.y) / (maxP0 - minP0);
	}

	Point2d A, B, A_, B_;
	A.x = srcPs[0].x;
	A.y = srcPs[0].y;
	B.x = srcPs[1].x;
	B.y = srcPs[1].y;

	A_.x = qstPs[0].x;
	A_.y = qstPs[0].y;
	B_.x = qstPs[1].x;
	B_.y = qstPs[1].y;
	double angle = estimateTriangleRotationRadian(A, B, A_, B_);
	cv::Matx22d Rt(cos(angle), -1 * sin(angle), sin(angle), cos(angle));
	cv::Point2d t1(0, 0);

	double a1, a2, tranX, tranY;
	Point2d current_pt, uv_copy;
	current_pt.x = (uv.x - averP0.x) / (maxP0 - minP0);
	current_pt.y = (uv.y - averP0.y) / (maxP0 - minP0);
	Point2d res;
	for (int i = 0; i != num; ++i)
	{
		res = applyAffineTransform(Rt, t1, srcPs[i]);
		srcPs[i] = res;
	}
	res = applyAffineTransform(Rt, t1, current_pt);
	current_pt = res;

	for (int i = 0; i != 3; ++i)
	{
		Point2d pt;
		pt.x = mTri[i].x;
		pt.y = mTri[i].y;
		res = applyAffineTransform(Rt, t1, pt);
		mTri[i].x = res.x;
		mTri[i].y = res.y;
	}

	GetTransformParams(a1, a2, tranX, tranY
		, srcPs, qstPs, current_pt);

	cv::Matx22d R(a1, -1 * a2, a2, a1);
	cv::Point2d t(tranX, tranY);

	cv::Point2d P_prime = applyAffineTransform(R, t, current_pt);
	uv_copy.x = P_prime.x*(maxP - minP) + averP.x;
	uv_copy.y = P_prime.y*(maxP - minP) + averP.y;

	srcPs.clear();
	qstPs.clear();

	for (int i = 0; i != 3; ++i)
	{
		current_pt.x = mTri[i].x;
		current_pt.y = mTri[i].y;

		P_prime = applyAffineTransform(R, t, current_pt);

		mTri[i].x = P_prime.x*(maxP - minP) + averP.x;
		mTri[i].y = P_prime.y*(maxP - minP) + averP.y;

		Point2d pt;
		pt.x = mTri[i].x;
		pt.y = mTri[i].y;
		srcPs.push_back(pt);

		pt.x = mTri[i].oci;
		pt.y = mTri[i].ori;
		qstPs.push_back(pt);
	}

	Point2d qstP = findCorrespondingPoint(uv_copy, srcPs, qstPs);
	uv.oci = qstP.x;
	uv.ori = qstP.y;
	GetBdbox(neigh_uvs, uv);
	return neigh_uvs.size();
}

double DataAnalysis::P2lineDist(Pt &p, Pt &a, Pt &b)
{
	double A = b.y - a.y;
	double B = a.x - b.x;
	double C = b.x * a.y - a.x * b.y;
	double denominator = b.euclidDist(a);

	// Check whether a and b coincide to avoid division by zero
	if (denominator == 0) {
		std::cerr << "Error: a and b are the same point!" << std::endl;
		return -1; // Return -1 as an error flag
	}

	return fabs(A * p.x + B * p.y + C) / denominator;
}


double DataAnalysis::P2lineDist4UV(Pt &p, Pt &a, Pt &b)
{
	double A = b.ori - a.ori;
	double B = a.oci - b.oci;
	double C = b.oci * a.ori - a.oci * b.ori;
	double denominator = b.BMPDist(a);

	// Check if a and b are the same point to avoid division by zero
	if (denominator == 0) {
		std::cerr << "Error: a and b are the same point!" << std::endl;
		return -1; // Return -1 as an error flag
	}

	return fabs(A * p.oci + B * p.ori + C) / denominator;
}


uchar DataAnalysis::GetGrayV(vector<Pt> &neigh_uvs, Pt &pos)
{
	uchar g01 = neigh_uvs[1].val;
	uchar g00 = neigh_uvs[0].val;
	uchar g10 = neigh_uvs[3].val;
	uchar g11 = neigh_uvs[2].val;
	double c = neigh_uvs[0].BMPDist(neigh_uvs[1]);
	double r = neigh_uvs[0].BMPDist(neigh_uvs[3]);

	double a = P2lineDist4UV(pos, neigh_uvs[0], neigh_uvs[3]);
	double b = P2lineDist4UV(pos, neigh_uvs[0], neigh_uvs[1]);

	uchar gray = BilineInterpola(g00, g01, g10, g11
		, a / r, b / c);

	return gray;
}

uchar DataAnalysis::GetGrayV(Pt &pos, TOPOTYPE type)
{
	Search(pos, type);
	double i = pos.ori - 0.5;
	double j = pos.oci - 0.5;
	//int gray1 = proc->sobel_mat.at<uchar>(i - 1, j - 1);
	//int gray2 = proc->sobel_mat.at<uchar>(i - 1, j) * 2;
	//int gray3 = proc->sobel_mat.at<uchar>(i - 1, j + 1);
	//int gray4 = proc->sobel_mat.at<uchar>(i, j - 1) * 2;
	//int gray5 = proc->sobel_mat.at<uchar>(i, j) * 4;
	//int gray6 = proc->sobel_mat.at<uchar>(i, j + 1) * 2;
	//int gray7 = proc->sobel_mat.at<uchar>(i + 1, j - 1);
	//int gray8 = proc->sobel_mat.at<uchar>(i + 1, j) * 2;
	//int gray9 = proc->sobel_mat.at<uchar>(i + 1, j + 1);

	Mat gray32; proc->sobel_mat.convertTo(gray32, CV_32F); // 0~255 ∏°µ„

	uchar v_gauss = interp3x3_gauss_LOESS1(gray32, j, i, 0.9, 0.9); // Ω®“È sigma °÷ 0.8~1.1

	return v_gauss;
}

double DataAnalysis::GetSubpixelPoint(double2 *point_set)
{
	memset(anti_matrix, 0, sizeof(double) * 9);
	memset(result, 0, sizeof(double) * 3);

	double norm0 = sqrt(pow(point_set[0][0], 2)
		+ pow(point_set[1][0], 2)
		+ pow(point_set[2][0], 2));

	for (int i = 0; i != 3; ++i)
	{
		matrix[i * 3 + 0] = pow(point_set[i][0] / norm0, 2);
		matrix[i * 3 + 1] = point_set[i][0] / norm0;
		matrix[i * 3 + 2] = 1;
		y[i] = point_set[i][1];
	}

	MatrixAnti(matrix, anti_matrix, 3);

	MultMatrix(anti_matrix, y, result, 3, 3, 1);

	return -0.5*norm0 * result[1] / result[0];
}

void DataAnalysis::GetBdySubPxPos1(const char* filePath)
{
	vector<Pt> subBdy;
	GetBdySubPxPos(subBdy);
	Size newSize(proc->c * MAGNIFIC, proc->r * MAGNIFIC);
	Mat resizedImage;
	resizedImage.create(newSize, CV_8UC3);

	Mat resizedSrcImage;
	resizedSrcImage.create(newSize, CV_8UC3);

	for (int i = 0; i != proc->r; ++i)
	{
		for (int j = 0; j != proc->c; ++j)
		{
			uchar d = proc->sobel_mat.at<uchar>(i, j);
			uchar d_src = proc->src_mat.at<uchar>(i, j);
			for (int k = 0; k != MAGNIFIC; ++k)
			{
				for (int l = 0; l != MAGNIFIC; ++l)
				{
					resizedImage.at<Vec3b>(i * MAGNIFIC + k
						, j * MAGNIFIC + l) = Vec3b(d, d, d);

					resizedSrcImage.at<Vec3b>(i * MAGNIFIC + k
						, j * MAGNIFIC + l) = Vec3b(d_src, d_src, d_src);
				}
			}

		}
	}

	int n = subBdy.size();
	Pt p1, p2;
	int preI;
	for (int i = 0; i != n; ++i)
	{
		preI = (n + i - 1) % n;
		p1 = subBdy[preI];
		p2 = subBdy[i];

		resizedImage.at<Vec3b>(proc->orderedBdy[i].y * MAGNIFIC
			, proc->orderedBdy[i].x * MAGNIFIC) = Vec3b(0, 255, 0);

		cv::Point pt1(p1.oci * MAGNIFIC, p1.ori * MAGNIFIC);
		cv::Point pt2(p2.oci * MAGNIFIC, p2.ori * MAGNIFIC);

		cv::Scalar color(0, 0, 255);
		cv::line(resizedImage, pt1, pt2, color, 1, cv::LINE_8);
	}

	n = proc->orderedBdy.size();

	for (int i = 0; i != n; ++i)
	{
		Point2d p = proc->orderedBdy[i];
		for (int k = 0; k != MAGNIFIC; ++k)
		{
			for (int l = 0; l != MAGNIFIC; ++l)
			{
				resizedSrcImage.at<Vec3b>(p.y * MAGNIFIC + k
					, p.x * MAGNIFIC + l) = Vec3b(0, 255, 0);
			}
		}
	}

	imwrite(filePath, resizedImage);
	imwrite("E:\\testPic\\pixel_edge_show.bmp", resizedSrcImage);
}

void DataAnalysis::GetBdySubPxPos(vector<Pt> &subpxBdy)
{
	int n = innerBdy.size();
	int iI, oI, size;
	Pt ibPt, obPt;
	Pt rPt, lPt;
	double bG, rG, lG, pos;
	double2 ps[3];
	vector<Pt> neigh_uvs;
	double2 *vs4output = new double2[n];
	for (int i = 0; i != n; ++i)
	{
		iI = innerBdy[i];
		ibPt = innerPs[iI];
		oI = outerBdy[i];
		obPt = outerPs[oI];

		lPt = ibPt;
		lPt.x -= bdyWeigh[i];
		rPt = obPt;
		rPt.x += bdyWeighOut[i];

		bG = ibPt.val;
		lG = GetGrayV(lPt, INNER);
		rG = GetGrayV(rPt, OUTER);

		ps[1][0] = 0;
		ps[0][0] = -0.2;
		ps[2][0] = 0.2;

		ps[1][1] = bG / 255.0;
		ps[0][1] = lG / 255.0;
		ps[2][1] = rG / 255.0;

		pos = GetSubpixelPoint(ps);

		if (pos < -0.2 || pos > 0.2)
		{
			printf("pos=%f bG=%f lG=%f rG=%f ibPt.oci=%f ibpt.ori=%f\n bdyWeigh[i]=%f bdyWeighOut[i]=%f\n\n"
				, pos, bG, lG, rG, ibPt.oci, ibPt.ori, bdyWeigh[i], bdyWeighOut[i]);
			//subpxBdy.push_back(ibPt);
			//continue;
		}
		Pt pt;

		neigh_uvs.clear();
		if (0 >= pos)
		{
			pt = ibPt;
			pt.x += pos* bdyWeigh[i];
			size = Search(neigh_uvs, pt, INNER);
		}
		else
		{
			pt = obPt;
			pt.x += pos* bdyWeighOut[i];
			size = Search(neigh_uvs, pt, OUTER);
		}

		vs4output[i][0] = pt.oci;
		vs4output[i][1] = pt.ori;

		if (4 == size)
		{
			pt.val = GetGrayV(neigh_uvs, pt);
			subpxBdy.push_back(pt);
		}

	}

	//Test
	write2f(vs4output, n, "C:\\Users\\1\\Desktop\\subpixel experimental analysis\\HoloRay\\triangle\\sigma1d1.csv");
	delete[] vs4output;
}


double DataAnalysis::GetBdySubPxPos()
{
	int n = innerBdy.size();
	int iI, oI;
	Pt ibPt, obPt;
	Pt rPt, lPt;
	double bG, rG, lG;
	double2 ps[3];
	ps[1][0] = 0;
	ps[0][0] = -1;
	ps[2][0] = 1;

	ps[1][1] = 0;
	ps[0][1] = 0;
	ps[2][1] = 0;
	for (int i = 0; i != n; ++i)
	{
		iI = innerBdy[i];
		ibPt = innerPs[iI];
		oI = outerBdy[i];
		obPt = outerPs[oI];

		lPt = ibPt;
		lPt.x -= bdyWeigh[i];
		rPt = obPt;
		rPt.x += bdyWeigh[i];

		bG = innerPs[iI].val;
		rG = GetGrayV(rPt, OUTER);
		lG = GetGrayV(lPt, INNER);

		ps[1][1] += bG / n;
		ps[0][1] += lG / n;
		ps[2][1] += rG / n;
	}

	double pos = GetSubpixelPoint(ps);
	return 1 + pos;
}



void  DataAnalysis::GetTriSet(vector<vector<Pt>> &triSet
	, Pt *ps, vector<vector<int>> &tridxSet
	, double deviThres, double maxY, double minY)
{
	int n = tridxSet.size();
	Pt p1, p2, p3, p1c, p2c, p3c;
	double ay, deviY;
	for (int i = 0; i != n; ++i)
	{
		p1 = ps[tridxSet[i][0]];
		p2 = ps[tridxSet[i][1]];
		p3 = ps[tridxSet[i][2]];

		ay = (p1.y + p2.y + p3.y) / 3;
		deviY = sqrt(pow(p1.y - ay, 2) / 3
			+ pow(p2.y - ay, 2) / 3
			+ pow(p3.y - ay, 2) / 3);
		if (deviY < deviThres)
		{
			vector<Pt> tri;
			tri.push_back(p1);
			tri.push_back(p2);
			tri.push_back(p3);
			triSet.push_back(tri);

			p1c = p1;
			p2c = p2;
			p3c = p3;

			p1c.y -= 1.0;
			p2c.y -= 1.0;
			p3c.y -= 1.0;
			if (p1c.y>minY || p2c.y > minY || p3c.y > minY)
			{
				vector<Pt> tri1;
				tri1.push_back(p1c);
				tri1.push_back(p2c);
				tri1.push_back(p3c);
				triSet.push_back(tri1);
			}


			p1c = p1;
			p2c = p2;
			p3c = p3;

			p1c.y += 1.0;
			p2c.y += 1.0;
			p3c.y += 1.0;

			if (p1c.y < maxY || p2c.y < maxY || p3c.y < maxY)
			{
				vector<Pt> tri2;
				tri2.push_back(p1c);
				tri2.push_back(p2c);
				tri2.push_back(p3c);
				triSet.push_back(tri2);
			}

		}
		else
		{
			if (p1.y <= ay && p2.y <= ay && p3.y >= ay)
			{
				vector<Pt> tri1;
				p1c = p1;
				p2c = p2;
				p3c = p3;
				p3c.y -= 1.0;
				tri1.push_back(p1c);
				tri1.push_back(p2c);
				tri1.push_back(p3c);
				triSet.push_back(tri1);

				vector<Pt> tri2;
				p1c = p1;
				p2c = p2;
				p3c = p3;
				p1c.y += 1.0;
				p2c.y += 1.0;
				tri2.push_back(p1c);
				tri2.push_back(p2c);
				tri2.push_back(p3c);
				triSet.push_back(tri2);
			}
			else if (p1.y >= ay && p2.y >= ay && p3.y <= ay)
			{
				vector<Pt> tri1;
				p1c = p1;
				p2c = p2;
				p3c = p3;
				p3c.y += 1.0;
				tri1.push_back(p1c);
				tri1.push_back(p2c);
				tri1.push_back(p3c);
				triSet.push_back(tri1);

				vector<Pt> tri2;
				p1c = p1;
				p2c = p2;
				p3c = p3;
				p1c.y -= 1.0;
				p2c.y -= 1.0;
				tri2.push_back(p1c);
				tri2.push_back(p2c);
				tri2.push_back(p3c);
				triSet.push_back(tri2);
			}
			else if (p1.y <= ay && p2.y >= ay && p3.y >= ay)
			{
				vector<Pt> tri1;
				p1c = p1;
				p2c = p2;
				p3c = p3;
				p1c.y += 1.0;
				tri1.push_back(p1c);
				tri1.push_back(p2c);
				tri1.push_back(p3c);
				triSet.push_back(tri1);

				vector<Pt> tri2;
				p1c = p1;
				p2c = p2;
				p3c = p3;
				p2c.y -= 1.0;
				p3c.y -= 1.0;
				tri2.push_back(p1c);
				tri2.push_back(p2c);
				tri2.push_back(p3c);
				triSet.push_back(tri2);
			}
			else if (p1.y >= ay && p2.y <= ay && p3.y <= ay)
			{
				vector<Pt> tri1;
				p1c = p1;
				p2c = p2;
				p3c = p3;
				p1c.y -= 1.0;
				tri1.push_back(p1c);
				tri1.push_back(p2c);
				tri1.push_back(p3c);
				triSet.push_back(tri1);

				vector<Pt> tri2;
				p1c = p1;
				p2c = p2;
				p3c = p3;
				p2c.y += 1.0;
				p3c.y += 1.0;
				tri2.push_back(p1c);
				tri2.push_back(p2c);
				tri2.push_back(p3c);
				triSet.push_back(tri2);
			}
			else if (p1.y <= ay && p2.y >= ay && p3.y <= ay)
			{
				vector<Pt> tri1;
				p1c = p1;
				p2c = p2;
				p3c = p3;
				p2c.y -= 1.0;
				tri1.push_back(p1c);
				tri1.push_back(p2c);
				tri1.push_back(p3c);
				triSet.push_back(tri1);

				vector<Pt> tri2;
				p1c = p1;
				p2c = p2;
				p3c = p3;
				p1c.y += 1.0;
				p3c.y += 1.0;
				tri2.push_back(p1c);
				tri2.push_back(p2c);
				tri2.push_back(p3c);
				triSet.push_back(tri2);
			}
			else if (p1.y >= ay && p2.y <= ay && p3.y >= ay)
			{
				vector<Pt> tri1;
				p1c = p1;
				p2c = p2;
				p3c = p3;
				p2c.y += 1.0;
				tri1.push_back(p1c);
				tri1.push_back(p2c);
				tri1.push_back(p3c);
				triSet.push_back(tri1);

				vector<Pt> tri2;
				p1c = p1;
				p2c = p2;
				p3c = p3;
				p1c.y -= 1.0;
				p3c.y -= 1.0;
				tri2.push_back(p1c);
				tri2.push_back(p2c);
				tri2.push_back(p3c);
				triSet.push_back(tri2);
			}
		}
	}
}

void DataAnalysis::FineWeigh(double *wDist, TOPOTYPE type, double thres, int maxIterNum)
{
	Pt pt[4096], lPt[4096];
	int num, size;
	double w, rate;
	double wAdd[4096][256] = { 0 };
	int iterCount = 16;
	num = innerBdy.size();
	if (INNER == type)
	{
		for (int i = 0; i != num; ++i)
		{
			wAdd[i][0] = bdyWeigh[i];
			pt[i] = innerPs[innerBdy[i]];
		}
	}
	else if (OUTER == type)
	{
		for (int i = 0; i != num; ++i)
		{
			wAdd[i][0] = bdyWeighOut[i];
			pt[i] = outerPs[outerBdy[i]];
		}
	}

	for (int i = 0; i != num; ++i)
	{
		lPt[i] = pt[i];

		w = wAdd[i][0];

		//if (w > 3)
		//{
		//	if (INNER == type)
		//	{
		//		lPt[i].x -= w;
		//	}
		//	else if (OUTER == type)
		//	{
		//		lPt[i].x += w;
		//	}

		//	Search(lPt[i], type);
		//	double d = pt[i].BMPDist(lPt[i]);
		//	wAdd[i][0] /= d;

		//	lPt[i] = pt[i];
		//}

		for (int j = 0; j != maxIterNum; ++j)
		{
			w = wAdd[i][j];

			if (INNER == type)
			{
				lPt[i].x -= w;
			}
			else if (OUTER == type)
			{
				lPt[i].x += w;
			}

			Search(lPt[i], type);
			if (w >= 0)
			{
				wDist[i] += pt[i].BMPDist(lPt[i]);
			}
			else
			{
				wDist[i] -= pt[i].BMPDist(lPt[i]);
			}

			rate = CalcWeigh(lPt[i], type, iterCount) / sqrt(0.7);

			wAdd[i][j + 1] = rate*(1.0 - wDist[i]);
			pt[i] = lPt[i];

			if (INNER == type)
			{
				bdyWeigh[i] += wAdd[i][j + 1];
			}
			else if (OUTER == type)
			{
				bdyWeighOut[i] += wAdd[i][j + 1];
			}

			if (abs((1.0 - wDist[i])) < thres)
			{
				break;
			}
		}
	}

}

double DataAnalysis::f(double x, Pt &p, TOPOTYPE type)
{
	Pt innerP = p;
	if (INNER == type)
	{
		innerP.x -= x;
	}
	else if (OUTER == type)
	{
		innerP.x += x;
	}

	Search(innerP, type);

	return p.BMPDist(innerP);
}

double DataAnalysis::solve_fx_equals_target(Pt &p, TOPOTYPE type,
	double target, double x0, double step)
{
	double low = x0;
	double high = x0 + step;

	// Expand the interval to the right until f(high) > target
	while (f(high, p, type) < target)
	{
		low = high;
		high *= 2; // Exponential growth
		if (high > 1e10)
		{
			std::cerr << "Function may have no solution in this interval or grows too slowly." << std::endl;
			exit(1);
		}
	}

	// Bisection method
	while (high - low > EPSILON) {
		double mid = (low + high) / 2;
		double val = f(mid, p, type);
		if (val < target)
			low = mid;
		else
			high = mid;
	}

	printf("distance to boundary = %lf\n", f((low + high) / 2, p, type));

	return (low + high) / 2;
}


double DataAnalysis::CalcWeigh(Pt &pt, TOPOTYPE type, int iterCount)
{
	kdtree::KDTreeResultVector result;
	std::vector<float> query(2);

	query[0] = pt.x;
	query[1] = pt.y;
	double wMax = 0;
	double area;
	if (INNER == type)
	{
		inner_tri_tree->n_nearest(query, iterCount, result);
		for (int i = 0; i != iterCount; ++i)
		{
			vector<Pt> tri = innerTriSet[result[i].idx];
			area = triangleArea4ori(tri);

			if (abs(area - REGUTRIAREA) < TOLERANCE)
			{
				wMax = averageEdgeLength(tri);
				break;
			}
		}
	}
	else if (OUTER == type)
	{
		outer_tri_tree->n_nearest(query, iterCount, result);
		for (int i = 0; i != iterCount; ++i)
		{
			vector<Pt> tri = outerTriSet[result[i].idx];

			area = triangleArea4ori(tri);

			if (abs(area - REGUTRIAREA) < TOLERANCE)
			{
				wMax = averageEdgeLength(tri);
				break;
			}
		}
	}
	return wMax;
}

double DataAnalysis::CalcWeigh(Pt &pt, TOPOTYPE type)
{
	kdtree::KDTreeResultVector result;
	std::vector<float> query(2);

	query[0] = pt.x;
	query[1] = pt.y;
	Pt p;
	vector<vector<Pt>> tris;
	double wMax = 0;
	double area;
	vector<Pt> tri;
	if (INNER == type)
	{
		inner_tri_tree->n_nearest(query, 1, result);
		tri = innerTriSet[result[0].idx];
	}
	else if (OUTER == type)
	{
		outer_tri_tree->n_nearest(query, 1, result);
		tri = outerTriSet[result[0].idx];
	}
	wMax = averageEdgeLength(tri);
	return wMax;
}

void DataAnalysis::CalcWeigh()
{
	int n = innerBdy.size();
	int preI, postI;
	double w, wOut;
	Pt p;
	for (int i = 0; i != n; ++i)
	{
		p = innerPs[innerBdy[i]];
		w = CalcWeigh(p, INNER);
		w = solve_fx_equals_target(p, INNER
			, 1.0, 0.0, w);

		p = outerPs[outerBdy[i]];
		wOut = CalcWeigh(p, OUTER);
		wOut = solve_fx_equals_target(p, OUTER
			, 1.0, 0.0, wOut);

		printf("w=%.15f wOut=%.15f\n", w, wOut);

		bdyWeigh.push_back(w);
		bdyWeighOut.push_back(wOut);
	}

	/*int n = innerBdy.size();
	int preI, postI;
	double preEdgeLen, postEdgeLen, w, preEdgeLenOut, postEdgeLenOut, wOut;
	Pt preP, postP, p;
	for (int i = 0; i != n; ++i)
	{
	preI = (n + i - 1) % n;
	postI = (n + i + 1) % n;
	preP = innerPs[innerBdy[preI]];
	p = innerPs[innerBdy[i]];
	postP = innerPs[innerBdy[postI]];

	DelOffset(preP, p, postP);

	preEdgeLen = p.euclidDist(preP) / p.BMPDist(preP);
	postEdgeLen = p.euclidDist(postP) / p.BMPDist(postP);
	w = (preEdgeLen + postEdgeLen) / 2;

	preP = outerPs[outerBdy[preI]];
	p = outerPs[outerBdy[i]];
	postP = outerPs[outerBdy[postI]];
	DelOffset(preP, p, postP);

	preEdgeLenOut = p.euclidDist(preP) / p.BMPDist(preP);
	postEdgeLenOut = p.euclidDist(postP) / p.BMPDist(postP);
	wOut = (preEdgeLenOut + postEdgeLenOut) / 2;
	bdyWeigh.push_back(w);
	bdyWeighOut.push_back(wOut);
	}*/
}

void DataAnalysis::GetBdbox(vector<Pt> &neigh_uvs, Pt &uv)
{
	Pt p1, p2, p3, p4;

	p1.oci = floor(uv.oci);
	p1.ori = floor(uv.ori);
	p1.val = proc->sobel_mat.at<uchar>(p1.ori, p1.oci);

	p2.oci = p1.oci + 1;
	p2.ori = p1.ori;
	p2.val = proc->sobel_mat.at<uchar>(p2.ori, p2.oci);

	p3.oci = p1.oci + 1;
	p3.ori = p1.ori + 1;
	p3.val = proc->sobel_mat.at<uchar>(p3.ori, p3.oci);

	p4.oci = p1.oci;
	p4.ori = p1.ori + 1;
	p4.val = proc->sobel_mat.at<uchar>(p4.ori, p4.oci);

	p1.oci += 0.5;
	p1.ori += 0.5;
	neigh_uvs.push_back(p1);

	p2.oci += 0.5;
	p2.ori += 0.5;
	neigh_uvs.push_back(p2);

	p3.oci += 0.5;
	p3.ori += 0.5;
	neigh_uvs.push_back(p3);

	p4.oci += 0.5;
	p4.ori += 0.5;
	neigh_uvs.push_back(p4);
}


void DataAnalysis::DelOffset(Pt &p1, Pt &p2, Pt &p3)
{
	if (p1.y - p2.y > 0.5 && p1.y - p3.y > 0.5)
	{
		p1.y -= 1.0;
	}
	else if (p1.y - p2.y < -0.5 && p1.y - p3.y < -0.5)
	{
		p1.y += 1.0;
	}
	else if (p2.y - p1.y > 0.5 && p2.y - p3.y > 0.5)
	{
		p2.y -= 1.0;
	}
	else if (p2.y - p1.y < -0.5 && p2.y - p3.y < -0.5)
	{
		p2.y += 1.0;
	}
	else if (p3.y - p1.y > 0.5 && p3.y - p2.y > 0.5)
	{
		p3.y -= 1.0;
	}
	else if (p3.y - p1.y < -0.5 && p3.y - p2.y < -0.5)
	{
		p3.y += 1.0;
	}
}

void DataAnalysis::OutputIntegralImg(const char *filePath)
{
	Mat R = Mat(2048, 2048, CV_8UC3);

	for (int i = 0; i != 2048; ++i)
	{
		for (int j = 0; j != 2048; ++j)
		{
			R.at<cv::Vec3b>(i, j) = cv::Vec3b(0, 0, 0);
		}
	}

	for (int i = 0; i != innerPs.size(); ++i)
	{
		Pt p = innerPs[i];

		R.at<Vec3b>((p.y - innerMin) * 1024 / (innerMax - innerMin)
			, (p.x - innerMin) * 1024 / (innerMax - innerMin))
			= Vec3b(p.val, p.val, p.val);
	}

	//for (int i = 0; i != innerPs.size(); ++i)
	//{
	//	Pt p = innerPs[i];

	//	R.at<Vec3b>((p.y - innerMin) * 1024 / (innerMax - innerMin)
	//		, (p.x - innerMin) * 1024 / (innerMax - innerMin))
	//		= Vec3b(p.val, p.val, p.val);
	//}

	//for (int i = 0; i != innerPs.size(); ++i)
	//{
	//	Pt p = innerPs[i];

	//	R.at<Vec3b>((p.y + 1.0 - innerMin) * 1024 / (innerMax - innerMin)
	//		, (p.x - innerMin) * 1024 / (innerMax - innerMin))
	//		= Vec3b(p.val, p.val, p.val);
	//}

	imwrite(filePath, R);
}