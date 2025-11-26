#pragma once
#ifndef ALGORITHM
#define ALGORITHM
/*!
*      \File Algo.h
*      \brief general algorithms
*	   \author Wei Feng
*      \date 11/21/2024
*/

#include <vector>
#include<cassert>
#include<algorithm>
#include <math.h>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#define EIGEN_STRONG_INLINE inline
#include <Eigen/Dense>
#include <set>
using namespace Eigen;
using namespace cv;
using namespace std;
#define PI 3.14159265358979323846264338327950288419716939937510
#define MAXV 1048576
#define TOLERANCE 1e-5
#define EPSILON 1e-11

struct Pt {
	double x, y;//the coordinates on the integral image
	uchar val;
	double ori, oci;

	Pt() : x(0), y(0), val(0), ori(0), oci(0) {}
	Pt(float x, float y, uchar val, float ori, float oci)
		: x(x), y(y), val(val), ori(ori), oci(oci) {}

	double norm()
	{
		return sqrt(pow(x, 2) + pow(y, 2));
	}

	// Vector subtraction
	Pt operator-(const Pt& other) const {
		return Pt(x - other.x
			, y - other.y
			, val - other.val
			, ori - other.ori
			, oci - other.oci);
	}

	// Dot product
	float dot(const Pt& other) const {
		return x * other.x + y * other.y;
	}

	float euclidDist(const Pt& other) const {
		return sqrt(pow(x - other.x, 2) + pow(y - other.y, 2));
	}

	float BMPDist(const Pt& other) const {
		return sqrt(pow(oci - other.oci, 2) + pow(ori - other.ori, 2));
	}
};


struct ComparePoints {
	bool operator()(const Pt& p1, const Pt& p2) const {
		const float epsilon = 1e-6f;
		if (fabs(p1.x - p2.x) > epsilon) {
			return p1.x < p2.x;
		}
		return p1.y < p2.y;
	}
};

struct FootResult {
	cv::Point2d adjustedPoint;
	bool adjusted;
};

template<class T>
void MatrixIdentity(T* matrix, int nrows)
{ //  matrix -  nrows*nrows
	memset(matrix, 0, nrows*nrows * sizeof(T));
	for (int i = 0; i<nrows; ++i)
		matrix[i*nrows + i] = 1;
}

template<class TT1, class TT2, class TT3>//>
void Doolittle(TT1* aa, TT2* bb, TT3* xx, int rows)
{// aa * xx = bb        root - xx[rows][rows]
	int k, i, j, t, ik;
	int* M = new int[rows];
	double  *s, *l, *u, *a, *b;
	double temp, smax = 0, *y, *x;
	s = new double[rows];
	l = new double[rows*rows];
	u = new double[rows*rows];
	a = new double[rows*rows];
	b = new double[rows];
	y = new double[rows];
	x = new double[rows];
	//  QA  =  LU
	for (i = 0; i<rows; ++i)
	{
		M[i] = 0;
		for (j = 0; j<rows; ++j)
		{
			a[i*rows + j] = aa[i*rows + j];
		}
	}
	for (k = 0; k<rows; ++k)
	{
		for (i = k; i<rows; ++i)
		{
			s[i] = a[i*rows + k];
			for (t = 0; t < k; ++t)
			{
				s[i] -= l[i*rows + t] * u[t*rows + k];
			}

			if (i == k)
			{
				smax = s[i];
				ik = i;
			}
			if (fabs(smax)<fabs(s[i]))
			{
				smax = s[i];
				ik = i;
			}
		}
		M[k] = ik;
		if (ik != k)
		{
			for (t = 0; t<k; ++t)
			{
				temp = l[k*rows + t];
				l[k*rows + t] = l[ik*rows + t];
				l[ik*rows + t] = temp;
			}
			for (t = k; t<rows; ++t)
			{
				temp = a[k*rows + t];
				a[k*rows + t] = a[ik*rows + t];
				a[ik*rows + t] = temp;
			}
			temp = s[k];
			s[k] = s[ik];
			s[ik] = temp;
		}
		u[k*rows + k] = s[k];
		if (k<rows - 1)
		{
			for (j = k + 1; j<rows; ++j)
			{
				u[k*rows + j] = a[k*rows + j];
				for (t = 0; t < k; ++t)
				{
					u[k*rows + j] -= l[k*rows + t] * u[t*rows + j];
				}

			}
			for (i = k + 1; i < rows; ++i)
			{
				l[i*rows + k] = s[i] / (u[k*rows + k] + 0.00001);
			}

		}
	}
	//Qb  =  Ly   AND   Ux  =   y
	for (j = 0; j<rows; ++j)
	{
		for (i = 0; i < rows; ++i)
		{
			b[i] = bb[i*rows + j];
		}

		for (k = 0; k<rows - 1; ++k)
		{
			t = M[k];
			temp = b[k];
			b[k] = b[t];
			b[t] = temp;
		}
		y[0] = b[0];
		for (i = 1; i<rows; ++i)
		{
			y[i] = b[i];
			for (t = 0; t < i; ++t)
			{
				y[i] -= l[i*rows + t] * y[t];
			}

		}
		x[rows - 1] = y[rows - 1] / (u[rows*rows - 1] + 0.00001);
		for (i = rows - 2; i>-1; --i)
		{
			x[i] = y[i];
			for (t = i + 1; t < rows; ++t)
			{
				x[i] -= u[i*rows + t] * x[t];
			}

			x[i] /= (u[i*rows + i] + 0.00001);
		}
		for (i = 0; i<rows; ++i)
		{
			xx[i*rows + j] = x[i];
		}
	}
	delete[]M;
	delete[]s;
	delete[]l;
	delete[]u;
	delete[]a;
	delete[]b;
	delete[]y;
	delete[]x;
	M = NULL;
	s = NULL;
	l = NULL;
	u = NULL;
	a = NULL;
	b = NULL;
	y = NULL;
	x = NULL;
}


// invertible matrice
template<class TT1, class TT2>
void MatrixAnti(TT1* Matrix, TT2* MatrixA, int rows)
{//  Matrix * MatrixA = I          I = E
	double* E = new double[rows*rows];
	MatrixIdentity(E, rows);

	//Doolittle solution
	Doolittle(Matrix, E, MatrixA, rows);
	delete[]E;
	E = NULL;
}

template<class TT1, class TT2, class TT3>
void  __declspec(dllexport) MultMatrix(TT1* M, TT2* M1, TT3* M2, int rows, int soms, int cols)
{
	//M - rows*soms    M1 - soms*cols   M2 - rows*cols
	for (int i = 0; i < rows; ++i)
	{
		for (int j = 0; j<cols; ++j)
		{
			M2[i*cols + j] = 0.0;
			for (int k = 0; k < soms; ++k)
			{
				M2[i*cols + j] += M[i*soms + k] * M1[k*cols + j];
			}

		}
	}

}

// Function to check if point P is inside the triangle formed by A, B, C
inline bool isPointInTriangle(const Pt& A, const Pt& B, const Pt& C, const Pt& P) {
	// Calculate vectors
	Pt v0 = C - A;
	Pt v1 = B - A;
	Pt v2 = P - A;

	// Calculate dot products
	float dot00 = v0.dot(v0);
	float dot01 = v0.dot(v1);
	float dot02 = v0.dot(v2);
	float dot11 = v1.dot(v1);
	float dot12 = v1.dot(v2);

	// Calculate the inverse denominator
	float invDenom = 1.0f / (dot00 * dot11 - dot01 * dot01);

	// Calculate barycentric coordinates u and v
	float u = (dot11 * dot02 - dot01 * dot12) * invDenom;
	float v = (dot00 * dot12 - dot01 * dot02) * invDenom;

	// Check if point P is inside the triangle
	return (u >= 0) && (v >= 0) && (u + v <= 1);
}

// Function to check if point p is inside quadrilateral
inline bool isPointInsideQuadrilateral(Pt p, const vector<Pt>& quad)
{
	if (quad.size() != 4)
	{
		std::cerr << "Error: Quadrilateral must have exactly 4 points.\n";
		return false;
	}

	// Step 1: Check bounding box
	double minX = std::min({ quad[0].x, quad[1].x, quad[2].x, quad[3].x });
	double maxX = std::max({ quad[0].x, quad[1].x, quad[2].x, quad[3].x });
	double minY = std::min({ quad[0].y, quad[1].y, quad[2].y, quad[3].y });
	double maxY = std::max({ quad[0].y, quad[1].y, quad[2].y, quad[3].y });

	if (p.x < minX || p.x > maxX || p.y < minY || p.y > maxY)
	{
		return false;
	}

	if (isPointInTriangle(quad[0], quad[1], quad[2], p))
	{
		return true;
	}
	else if (isPointInTriangle(quad[0], quad[2], quad[3], p))
	{
		return true;
	}
	else
	{
		return false;
	}
}

inline int classify2(vector<int> &t1, double &cluster1, vector<int> &t2
	, double &cluster2, vector<Pt> & ps4)
{
	cluster1 = ps4[0].y;
	t1.push_back(0);
	int t_num = 1;
	for (int i = 1; i != 4; ++i)
	{
		if (abs(cluster1 - ps4[i].y) > 0.5)
		{
			cluster2 = ps4[i].y;
			t2.push_back(i);
			t_num = 2;
		}
		else
		{
			t1.push_back(i);
		}
	}

	return t_num;
}

inline uchar BilineInterpola(uchar g00, uchar g01, uchar g10, uchar g11
	, double a, double b)
{
	return b*(a*g11 + (1 - a)*g01) + (1 - b)*(a*g10 + (1 - a)*g00);
}

inline uchar getGray(Mat &mat, Point2d p)
{
	double a = p.y - floor(p.y);
	double b = p.x - floor(p.x);

	uchar g00, g01, g10, g11;
	if (a <= 0.5 && b <= 0.5)
	{
		g00 = mat.at<uchar>(floor(p.y) - 1, floor(p.x) - 1);
		g01 = mat.at<uchar>(floor(p.y) - 1, floor(p.x));
		g10 = mat.at<uchar>(floor(p.y), floor(p.x) - 1);
		g11 = mat.at<uchar>(floor(p.y), floor(p.x));
	}
	else if (a <= 0.5 && b >= 0.5)
	{
		g00 = mat.at<uchar>(floor(p.y) - 1, floor(p.x));
		g01 = mat.at<uchar>(floor(p.y) - 1, floor(p.x) + 1);
		g10 = mat.at<uchar>(floor(p.y), floor(p.x));
		g11 = mat.at<uchar>(floor(p.y), floor(p.x) + 1);
	}
	else if (a >= 0.5 && b <= 0.5)
	{
		g00 = mat.at<uchar>(floor(p.y), floor(p.x) - 1);
		g01 = mat.at<uchar>(floor(p.y), floor(p.x));
		g10 = mat.at<uchar>(floor(p.y) + 1, floor(p.x) - 1);
		g11 = mat.at<uchar>(floor(p.y) + 1, floor(p.x));
	}
	else if (a >= 0.5 && b >= 0.5)
	{
		g00 = mat.at<uchar>(floor(p.y), floor(p.x));
		g01 = mat.at<uchar>(floor(p.y), floor(p.x) + 1);
		g10 = mat.at<uchar>(floor(p.y) + 1, floor(p.x));
		g11 = mat.at<uchar>(floor(p.y) + 1, floor(p.x) + 1);
	}

	return BilineInterpola(g00, g01, g10, g11
		, a, b);
}

template<class T>
double crossProduct(const T& p1, const T& p2, const T& p3)
{
	return (p2.x - p1.x) * (p3.y - p1.y) - (p2.y - p1.y) * (p3.x - p1.x);
}

// Check whether two points are equal (considering floating-point tolerance)
template<class T>
bool __declspec(dllexport) isSamePoint(const T &P1, const T &P2)
{
	return (abs(P1.x - P2.x) < EPSILON) && (abs(P1.y - P2.y) < EPSILON);
}

// Check whether a point lies on a line segment
template<class T>
bool __declspec(dllexport) containsPoint(const T &lineP, const T &lineQ, const T &point)
{
	// First, check if the point coincides with either endpoint
	if (isSamePoint(point, lineP) || isSamePoint(point, lineQ))
	{
		return true; // The point is an endpoint
	}

	// Compute the cross product to test if the point lies on the supporting line of the segment
	double crossPro = crossProduct(lineP, lineQ, point);
	if (std::fabs(crossPro) > EPSILON) {
		return false; // Nonzero cross product: the point is not on the supporting line
	}

	// Compute the dot product to test if the point lies within the segment bounds
	double dotProduct = (point.x - lineP.x) * (lineQ.x - lineP.x)
		+ (point.y - lineP.y) * (lineQ.y - lineP.y);
	if (dotProduct < 0.0f) {
		return false; // The point lies beyond lineP in the opposite direction
	}

	// Compute the squared length of the segment
	double squaredLengthBA = (lineQ.x - lineP.x) * (lineQ.x - lineP.x)
		+ (lineQ.y - lineP.y) * (lineQ.y - lineP.y);
	if (dotProduct > squaredLengthBA) {
		return false; // The point lies beyond lineQ, outside the segment
	}

	return true; // The point lies on the segment
}

template<class T>
int __declspec(dllexport) pointInTriangle(const T& A, const T& B, const T& C, const T& p)
{
	double c1 = crossProduct(A, B, p);
	double c2 = crossProduct(B, C, p);
	double c3 = crossProduct(C, A, p);

	// Check whether the point is strictly inside the triangle (all cross products have the same sign)
	if ((c1 > 0 && c2 > 0 && c3 > 0) || (c1 < 0 && c2 < 0 && c3 < 0))
		return 1; // Inside triangle

				  // Check whether the point lies on an edge
	if (std::abs(c1) < EPSILON && containsPoint(A, B, p)) return 2;  // On edge AB
	if (std::abs(c2) < EPSILON && containsPoint(B, C, p)) return 3;  // On edge BC
	if (std::abs(c3) < EPSILON && containsPoint(C, A, p)) return 4;  // On edge CA

	return 0; // Outside triangle
}


// Compare whether two Point2d objects are approximately equal
inline bool arePointsEqual(const cv::Point2d& p1,
	const cv::Point2d& p2)
{
	return (std::abs(p1.x - p2.x) < EPSILON)
		&& (std::abs(p1.y - p2.y) < EPSILON);
}



inline bool isAlmostEqual(const Point2d& p1, const Point2d& p2)
{
	return (std::fabs(p1.x - p2.x) < TOLERANCE)
		&& (std::fabs(p1.y - p2.y) < TOLERANCE);
}


// Check whether point (px, py) lies on the segment (x1, y1) ¨C (x2, y2)
inline bool isPointOnSegment(Point2d p, Point2d a, Point2d b) {
	double crossProduct = (p.x - a.x) * (b.y - a.y) - (p.y - a.y) * (b.x - a.x);
	if (fabs(crossProduct) > 1e-6) return false;  // Not collinear

	double dotProduct = (p.x - a.x) * (b.x - a.x) + (p.y - a.y) * (b.y - a.y);
	if (dotProduct < 0) return false;  // Beyond endpoint a

	double squaredLengthBA = (b.x - a.x) * (b.x - a.x)
		+ (b.y - a.y) * (b.y - a.y);
	if (dotProduct > squaredLengthBA) return false;  // Beyond endpoint b

	return true;  // Lies on the segment
}


// Normalize a triangle: compute its centroid and scale it to unit size
inline void normalizeTriangle(const vector<Point2d>& T,
	vector<Point2d>& T_norm,
	Point2d& center,
	double& scale)
{
	// Compute centroid
	center.x = (T[0].x + T[1].x + T[2].x) / 3.0;
	center.y = (T[0].y + T[1].y + T[2].y) / 3.0;

	// Compute maximum distance to centroid
	double maxDist = 0;
	for (const auto& pt : T) {
		double dist = norm(pt - center);
		if (dist > maxDist) maxDist = dist;
	}

	// Compute scale factor (normalize the triangle to unit radius)
	scale = (maxDist > 0) ? (1.0 / maxDist) : 1.0;

	// Apply normalization
	T_norm.clear();
	for (const auto& pt : T) {
		T_norm.push_back((pt - center) * scale);
	}
}


// Compute the barycentric coordinates of point A with respect to triangle T1
inline Vec3d computeBarycentricCoords(const Point2d& A,
	const vector<Point2d>& T1)
{
	Matrix3d M;
	M << T1[0].x, T1[1].x, T1[2].x,
		T1[0].y, T1[1].y, T1[2].y,
		1.0, 1.0, 1.0;

	Vector3d rhs(A.x, A.y, 1.0);
	Vector3d lambda = M.colPivHouseholderQr().solve(rhs);

	return Vec3d(lambda(0), lambda(1), lambda(2));
}


// Compute the corresponding point in triangle T2 using barycentric coordinates
inline Point2d computeCorrespondingPoint(const Vec3d& lambda,
	const vector<Point2d>& T2)
{
	return lambda[0] * T2[0] + lambda[1] * T2[1] + lambda[2] * T2[2];
}



// Given triangles T1 and T2, and a query point A in T1, compute the affine-corresponding point in T2
inline Point2d findCorrespondingPoint_affine(const Point2d& A,
	const std::vector<Point2d>& T1,
	const std::vector<Point2d>& T2)
{
	Eigen::Vector2d P0(T1[0].x, T1[0].y);
	Eigen::Vector2d P1(T1[1].x, T1[1].y);
	Eigen::Vector2d P2(T1[2].x, T1[2].y);
	Eigen::Vector2d Q0(T2[0].x, T2[0].y);
	Eigen::Vector2d Q1(T2[1].x, T2[1].y);
	Eigen::Vector2d Q2(T2[2].x, T2[2].y);
	Eigen::Vector2d P(A.x, A.y);

	Eigen::Matrix2d A_mat;
	A_mat.col(0) = P1 - P0;
	A_mat.col(1) = P2 - P0;

	Eigen::Vector2d coeffs = A_mat.inverse() * (P - P0);

	Eigen::Vector2d result = Q0 + coeffs[0] * (Q1 - Q0) + coeffs[1] * (Q2 - Q0);
	return Point2d(result.x(), result.y());
}


// Compute the Euclidean distance between two points
template<class T>
double __declspec(dllexport) dist(const T& P1, const T& P2) {
	return std::sqrt((P1.x - P2.x) * (P1.x - P2.x)
		+ (P1.y - P2.y) * (P1.y - P2.y));
}


// Compute the area of a triangle
template<class PointT>
double __declspec(dllexport) triangleArea(const std::vector<PointT>& tri)
{
	assert(tri.size() == 3);
	return 0.5 * std::fabs(
		tri[0].x * (tri[1].y - tri[2].y) +
		tri[1].x * (tri[2].y - tri[0].y) +
		tri[2].x * (tri[0].y - tri[1].y));
}


inline double triangleArea(const Point2d t1, const Point2d t2, const Point2d t3)
{
	return 0.5 * fabs(
		t1.x * (t2.y - t3.y) +
		t2.x * (t3.y - t1.y) +
		t3.x * (t1.y - t2.y));
}

template<class T>
double __declspec(dllexport) triangleArea4ori(const vector<T>& T) {
	return fabs((T[0].oci * (T[1].ori - T[2].ori) + T[1].oci
		* (T[2].ori - T[0].ori) + T[2].oci * (T[0].ori - T[1].ori)) / 2.0);
}

inline double clamp(double x, double a, double b) {
	return std::max(a, std::min(b, x));
}

inline double angleBetween(const Point2d& u, const Point2d& v)
{
	double dot = u.ddot(v);
	double norm_u = norm(u);
	double norm_v = norm(v);

	if (norm_u == 0 || norm_v == 0) return 0.0;  // Avoid division by zero

	double cos_theta = dot / (norm_u * norm_v);
	cos_theta = clamp(cos_theta, -1.0, 1.0);     // Ensure numerical stability (for VS2013)

	return std::acos(cos_theta);
}


inline double triangleMinAngle(const Point2d& A, const Point2d& B, const Point2d& C)
{
	double angleA = angleBetween(B - A, C - A);
	double angleB = angleBetween(A - B, C - B);
	double angleC = angleBetween(A - C, B - C);
	return std::min(std::min(angleA, angleB), angleC);
}

// Compute equilateral triangle quality score
template<class T>
double __declspec(dllexport) triangleQuality(const vector<T>& Tr) {
	if (Tr.size() != 3) {
		cerr << "Error: Triangle must have exactly 3 points!" << endl;
		return numeric_limits<double>::max();
	}

	double L1 = dist(Tr[0], Tr[1]);
	double L2 = dist(Tr[1], Tr[2]);
	double L3 = dist(Tr[2], Tr[0]);

	// Compute edge length ratio error
	double maxL = max({ L1, L2, L3 });
	double minL = min({ L1, L2, L3 });
	double E_ratio = maxL / minL - 1.0;

	// Compute angle deviation error
	double thetaA = acos((L2 * L2 + L3 * L3 - L1 * L1) / (2 * L2 * L3)) * 180.0 / PI;
	double thetaB = acos((L1 * L1 + L3 * L3 - L2 * L2) / (2 * L1 * L3)) * 180.0 / PI;
	double thetaC = 180.0 - thetaA - thetaB;

	double E_angle = (fabs(thetaA - 60) + fabs(thetaB - 60) + fabs(thetaC - 60)) / 3.0;

	// Compute aspect ratio error
	double s = (L1 + L2 + L3) / 2.0;
	double area = triangleArea(Tr);
	double R = (L1 * L2 * L3) / (4.0 * area);
	double r = area / s;
	double E_aspect = fabs(R / r - 2.0);

	// Final composite score
	double w1 = 1.0, w2 = 1.5, w3 = 1.2;
	return w1 * E_ratio + w2 * (E_angle / 60) + w3 * E_aspect;
}

// Find the triangle closest to an equilateral triangle
template<class T>
vector<T> __declspec(dllexport) findBestTriangle(const vector<vector<T>>& triangles)
{
	vector<T> bestTriangle;
	double minScore = numeric_limits<double>::max();

	for (const auto& S : triangles)
	{
		double score = triangleQuality(S);
		if (score < minScore)
		{
			minScore = score;
			bestTriangle = S;
		}
	}
	return bestTriangle;
}



// Compute the average edge length of a triangle
template<class T>
double __declspec(dllexport) averageEdgeLength(const vector<T>& tri) {
	if (tri.size() != 3) {
		cerr << "Error: Triangle must have exactly 3 points!" << endl;
		return -1.0;  // Return invalid value
	}

	double L1 = dist(tri[0], tri[1]);
	double L2 = dist(tri[1], tri[2]);
	double L3 = dist(tri[2], tri[0]);

	return (L1 + L2 + L3) / 3.0;
}

template<class T>
int __declspec(dllexport) nearestAngle(const vector<T>& tri, T&v) {
	if (tri.size() != 3) {
		cerr << "Error: Triangle must have exactly 3 points!" << endl;
		return -1.0; 
	}

	double L1 = dist(tri[0], v);
	double L2 = dist(tri[1], v);
	double L3 = dist(tri[2], v);

	if (L1 <= L2 && L1 <= L3)
	{
		return 0;
	}
	else if (L2 <= L1 && L2 <= L3)
	{
		return 1;
	}
	else if (L3 <= L1 && L3 <= L2)
	{
		return 2;
	}
	else
	{
		printf("nearestAngle return 3\n");
		return 3;
	}
}

// Extract a unique set of points from a list of triangles
inline vector<Pt> extractUniqueVertices(const vector<vector<Pt>>& triangles)
{
	set<Pt, ComparePoints> uniquePoints;

	// Traverse all triangle vertices and insert them into the set
	for (size_t i = 0; i < triangles.size(); i++) {
		for (size_t j = 0; j < triangles[i].size(); j++) {
			uniquePoints.insert(triangles[i][j]);
		}
	}

	// Convert the set to a vector for output
	return vector<Pt>(uniquePoints.begin(), uniquePoints.end());
}


// Select all triangles from point set M that contain point P
template<class T>
vector<vector<T>> __declspec(dllexport) findTrianglesContainingP(const vector<T>& M, const T& P)
{
	vector<vector<T>> result;

	int N = static_cast<int>(M.size());
	if (N < 3) return result;

	for (int i = 0; i < N - 2; ++i)
	{
		for (int j = i + 1; j < N - 1; ++j)
		{
			for (int k = j + 1; k < N; ++k)
			{
				const T &A = M[i], &B = M[j], &C = M[k];

				if (pointInTriangle(A, B, C, P))
				{
					result.push_back({ A, B, C });
				}
			}
		}
	}
	return result;
}


// Distance from point P to line AB
inline double ptLineDist(const Point2d& P, const Point2d& A, const Point2d& B)
{
	double numerator = fabs((B.x - A.x)*(A.y - P.y) - (A.x - P.x)*(B.y - A.y));
	double denominator = sqrt(pow(B.x - A.x, 2) + pow(B.y - A.y, 2));
	return numerator / denominator;
}


inline double squaredDistance(const Point2d& a, const Point2d& b)
{
	return (a - b).dot(a - b);  // Recommended: more efficient
								// Alternatively: return cv::norm(a - b) * cv::norm(a - b); 
								// Not recommended: redundant square root
}


inline double pointToSegmentDistance(const Point2d& p,
	const Point2d& a, const Point2d& b) {

	double abx = b.x - a.x;
	double aby = b.y - a.y;
	double apx = p.x - a.x;
	double apy = p.y - a.y;

	double abLenSq = abx * abx + aby * aby;
	if (abLenSq == 0.0)
		return std::numeric_limits<double>::max(); // a and b coincide

	double t = (abx * apx + aby * apy) / abLenSq;

	if (t < 0.0 || t > 1.0) {
		return std::numeric_limits<double>::max(); // Projection point is not on the segment
	}

	cv::Point2d projection(a.x + t * abx, a.y + t * aby);
	return std::sqrt(squaredDistance(p, projection));
}



// Coordinate interpolation: fA, fB, fC are coordinates (vectors), not scalar values
inline Vector2d meanValueCoordInterp2D(
	const Vector2d& P,
	const Vector2d& A, const Vector2d& fA,
	const Vector2d& B, const Vector2d& fB,
	const Vector2d& C, const Vector2d& fC)
{
	Vector2d v0 = A - P;
	Vector2d v1 = B - P;
	Vector2d v2 = C - P;

	double d0 = v0.norm();
	double d1 = v1.norm();
	double d2 = v2.norm();

	// Prevent division by zero
	const double eps = 1e-12;
	d0 = std::max(d0, eps);
	d1 = std::max(d1, eps);
	d2 = std::max(d2, eps);

	// Angle computation (use dot product + clamp to prevent precision errors in acos)
	auto safe_acos = [](double x) {
		return std::acos(std::max(-1.0, std::min(1.0, x)));
	};

	double theta0 = safe_acos(v0.dot(v1) / (d0 * d1));
	double theta1 = safe_acos(v1.dot(v2) / (d1 * d2));
	double theta2 = safe_acos(v2.dot(v0) / (d2 * d0));

	double w0 = (std::tan(theta2 / 2.0) + std::tan(theta0 / 2.0)) / d0;
	double w1 = (std::tan(theta0 / 2.0) + std::tan(theta1 / 2.0)) / d1;
	double w2 = (std::tan(theta1 / 2.0) + std::tan(theta2 / 2.0)) / d2;

	double w_sum = w0 + w1 + w2;
	w0 /= w_sum;
	w1 /= w_sum;
	w2 /= w_sum;

	// Coordinate interpolation: f(P) = w0 * fA + w1 * fB + w2 * fC
	return w0 * fA + w1 * fB + w2 * fC;
}


// Compute the barycentric coordinates of point P with respect to triangle ABC
inline Vec3d computeBarycentricCoords(const Point2d& P,
	const Point2d& A,
	const Point2d& B,
	const Point2d& C) {
	// Can use either area-based or vector-based method
	Vector2d v0(B.x - A.x, B.y - A.y);
	Vector2d v1(C.x - A.x, C.y - A.y);
	Vector2d v2(P.x - A.x, P.y - A.y);

	double d00 = v0.dot(v0);
	double d01 = v0.dot(v1);
	double d11 = v1.dot(v1);
	double d20 = v2.dot(v0);
	double d21 = v2.dot(v1);

	double denom = d00 * d11 - d01 * d01;
	if (std::abs(denom) < 1e-10)
		return Vec3d(1.0 / 3, 1.0 / 3, 1.0 / 3); // Prevent degenerate case

	double v = (d11 * d20 - d01 * d21) / denom;
	double w = (d00 * d21 - d01 * d20) / denom;
	double u = 1.0 - v - w;

	return Vec3d(u, v, w);
}


// Barycentric interpolation version for triangles
inline Point2d findCorrespondingPoint_Barycentric(const Point2d& A,
	const vector<Point2d>& T1,
	const vector<Point2d>& T2) {
	// Compute barycentric coordinates of A in triangle T1
	Vec3d lambda = computeBarycentricCoords(A, T1[0], T1[1], T1[2]);

	// Use barycentric coordinates to interpolate the corresponding point in T2
	Point2d B = lambda[0] * T2[0] + lambda[1] * T2[1] + lambda[2] * T2[2];
	return B;
}


// Find the corresponding point between two similar triangles
inline Point2d findCorrespondingPoint(const Point2d& A, const vector<Point2d>& T1,
	const vector<Point2d>& T2) {
	vector<Point2d> T1_norm, T2_norm;
	Point2d center1, center2;
	double scale1, scale2;

	normalizeTriangle(T1, T1_norm, center1, scale1);
	normalizeTriangle(T2, T2_norm, center2, scale2);

	Point2d A_norm = (A - center1) * scale1;

	Vector2d A_norm_v, T1_a, T1_b, T1_c, T2_a, T2_b, T2_c;
	A_norm_v(0) = A_norm.x;
	A_norm_v(1) = A_norm.y;
	T1_a(0) = T1_norm[0].x;
	T1_a(1) = T1_norm[0].y;

	T1_b(0) = T1_norm[1].x;
	T1_b(1) = T1_norm[1].y;

	T1_c(0) = T1_norm[2].x;
	T1_c(1) = T1_norm[2].y;

	T2_a(0) = T2_norm[0].x;
	T2_a(1) = T2_norm[0].y;

	T2_b(0) = T2_norm[1].x;
	T2_b(1) = T2_norm[1].y;

	T2_c(0) = T2_norm[2].x;
	T2_c(1) = T2_norm[2].y;

	Vector2d B_norm = meanValueCoordInterp2D(A_norm_v, T1_a,
		T2_a, T1_b, T2_b, T1_c, T2_c);

	return Point2d(B_norm(0) / scale2 + center2.x, B_norm(1) / scale2 + center2.y);
}

template<class T>
void __declspec(dllexport) GetWu1(T Wu1[], const T *matrixB
	, const T *P, const T *l, int N, int U)
{
	T *transposedB = new T[N*U];
	T *transposedBP = new T[N*U];
	MatrixTranspose(matrixB, transposedB, N, U);
	MultMatrix(transposedB, P, transposedBP, U, N, N);
	MultMatrix(transposedBP, l, Wu1, U, N, 1);

	delete[] transposedB;
	delete[] transposedBP;
}

//matrix[3][4] change to  matrix0[4][3]
template<class T, class T0>
void __declspec(dllexport) MatrixTranspose(T* matrix, T0* matrix0, int nrows, int ncols)
{// matrix - nrows*ncols    matrix0 - ncols*nrows
	for (int i = 0; i < nrows; ++i)
	{
		for (int j = 0; j < ncols; ++j)
		{
			matrix0[j*nrows + i] = matrix[i*ncols + j];
		}
	}
}


template<class TT1>
void __declspec(dllexport) GetNBB(TT1 nbb[], const TT1 *matrixB, const TT1 *P, int N, int U)
{
	TT1 *transposedB = new TT1[N*U];
	TT1 *transposedBP = new TT1[N*U];
	MatrixTranspose(matrixB, transposedB, N, U);
	MultMatrix(transposedB, P, transposedBP, U, N, N);
	MultMatrix(transposedBP, matrixB, nbb, U, N, U);

	delete[] transposedB;
	delete[] transposedBP;
}


/*!
weighted indirect adjustment model
\param correction
\param matrixB
\param l
\param P
\param N
\param U
*/
template<class T>
void __declspec(dllexport) GetCorrection(T correction[], const T *matrixB
	, const T *l, const T *P, int N, int U)
{
	T *nbb = new T[U*U];
	T *inverseForNbb = new T[U*U];
	T *Wu1 = new T[U];

	GetNBB<T>(nbb, matrixB, P, N, U);
	MatrixAnti(nbb, inverseForNbb, U);
	GetWu1<T>(Wu1, matrixB, P, l, N, U);

	MultMatrix(inverseForNbb, Wu1, correction, U, U, 1);

	delete[] nbb;
	delete[] inverseForNbb;
	delete[] Wu1;
}

template<class T>
void getL(T l[], const double a1, const double a2
	, vector<Point2d> &p0, vector<Point2d> &p, T delta_x, T delta_y)
{
	assert(p0.size() == p.size());
	int N = p0.size();
	memset(l, 0, sizeof(T)*N * 2);

	for (int i = 0; i != N; ++i)
	{
		l[i * 2 + 0] = -1 * a1*p0[i].x + a2*p0[i].y - delta_x + p[i].x;
		l[i * 2 + 1] = -1 * a2*p0[i].x - a1*p0[i].y - delta_y + p[i].y;
	}
}

template<class T>
void getMatrixB(T matrixB[], vector<Point2d> &ps)
{
	int N = ps.size();
	memset(matrixB, 0, sizeof(T)*N * 2 * 4);

	for (int i = 0; i != N; ++i)
	{
		matrixB[i * 8 + 0] = ps[i].x;
		matrixB[i * 8 + 1] = -1 * ps[i].y;
		matrixB[i * 8 + 2] = 1;
		matrixB[i * 8 + 3] = 0;

		matrixB[i * 8 + 4 + 0] = ps[i].y;
		matrixB[i * 8 + 4 + 1] = ps[i].x;
		matrixB[i * 8 + 4 + 2] = 0;
		matrixB[i * 8 + 4 + 3] = 1;
	}
}

// Compute the squared distances between all sample pairs and return the median
template<class T>
T __declspec(dllexport) compute_median_distance(const vector<Point2d>& data) {
	std::vector<T> distances;

	for (size_t i = 0; i < data.size(); ++i) {
		for (size_t j = i + 1; j < data.size(); ++j) {
			distances.push_back(norm(data[i] - data[j]));
		}
	}

	if (distances.empty())
		return 0.0;

	size_t n = distances.size();
	std::nth_element(distances.begin(), distances.begin() + n / 2, distances.end());
	T median;

	if (n % 2 == 0) {
		T a = *std::max_element(distances.begin(), distances.begin() + n / 2);
		T b = *std::min_element(distances.begin() + n / 2, distances.end());
		median = (a + b) / 2.0;
	}
	else {
		median = distances[n / 2];
	}

	return std::sqrt(median / 2.0);  // sigma
}


// Gaussian kernel function: compute the kernel weight of x relative to center c
template<class T>
T __declspec(dllexport) gaussian_kernel(const Point2d& x, const Point2d& c, T sigma) {
	T dist2 = norm(x - c);
	return std::exp(-dist2 / (2.0 * sigma * sigma));
}



template<class T>
void __declspec(dllexport) getP4Curv(T *P, vector<Point2d> &ps, Point2d curr_point)
{
	T dist = 0;
	T w = 0;

	int num = ps.size();

	T sigma = compute_median_distance<T>(ps);

	for (int i = 0; i != num * 2; ++i)
	{
		for (int j = 0; j != num * 2; ++j)
		{
			if (i != j)
			{
				P[i*num * 2 + j] = 0;
			}
			else
			{
				w = gaussian_kernel(ps[i / 2], curr_point, sigma);
				P[i*num * 2 + j] = w;
			}
		}
	}
}

template<class T>
void __declspec(dllexport) getDeltaVector(T DeltaV[], const T a1, const T a2
	, vector<Point2d> &p0, vector<Point2d> &p, Point2d curr_point, T delta_x, T delta_y)
{
	assert(p0.size() == p.size());
	int N = p0.size();
	T *matrixB = new T[N * 2 * 4];
	T *l = new T[N * 2];
	T P[4096] = { 0 };
	getP4Curv(P, p0, curr_point);
	getMatrixB(matrixB, p0);
	getL(l, a1, a2, p0, p, delta_x, delta_y);
	GetCorrection(DeltaV, matrixB, l, P, N * 2, 4);

	delete[] matrixB;
	delete[] l;
}


inline void GetTransformParams(double &a1, double &a2, double &tran_x, double &tran_y
	, vector<Point2d> &p0, vector<Point2d> &p, Point2d current_pt)
{
	a1 = 0;
	a2 = 0;
	tran_x = 0;
	tran_y = 0;
	double deltaV[4] = { 0 };
	assert(p0.size() == p.size());

	do
	{
		getDeltaVector(deltaV, a1, a2, p0, p, current_pt, tran_x, tran_y);

		a1 += deltaV[0];
		a2 += deltaV[1];

		tran_x += deltaV[2];
		tran_y += deltaV[3];
	} while (sqrt(pow(deltaV[0], 2) + pow(deltaV[1], 2)
		+ pow(deltaV[2], 2) + pow(deltaV[3], 2)) < 0.00001);
}

inline Point2d applyAffineTransform(
	const Matx22d& A,
	const Point2d& t,
	const Point2d& P)
{
	double x_new = A(0, 0) * P.x + A(0, 1) * P.y + t.x;
	double y_new = A(1, 0) * P.x + A(1, 1) * P.y + t.y;
	return Point2d(x_new, y_new);
}

// Given vectors A¡úB and A¡ä¡úB¡ä from two triangles, return the rotation angle (in radians, range: 0 to 2¦Ð)
inline double estimateTriangleRotationRadian(
	const cv::Point2d& A, const cv::Point2d& B,
	const cv::Point2d& A_, const cv::Point2d& B_)
{
	cv::Point2d v1 = B - A;
	cv::Point2d v2 = B_ - A_;

	double angle1 = std::atan2(v1.y, v1.x);
	double angle2 = std::atan2(v2.y, v2.x);

	double angle_rad = angle2 - angle1;

	// Map to [0, 2¦Ð)
	if (angle_rad < 0)
		angle_rad += 2 * CV_PI;

	return angle_rad; // Unit: radians
}


inline Mat gaussian1D(float sigma, int ksize)
{
	// 1D Gaussian kernel (column vector), L1-normalized
	CV_Assert(sigma > 0 && (ksize % 2 == 1));
	int half = ksize / 2;
	cv::Mat k(ksize, 1, CV_32F);
	float s2 = sigma * sigma;
	float sum = 0.f;

	for (int i = -half; i <= half; ++i) {
		float val = std::exp(-0.5f * (i * i) / s2);
		k.at<float>(i + half, 0) = val;
		sum += val;
	}
	k /= sum; // L1 normalize
	return k;
}

inline Mat dgaussian1D(float sigma, int ksize)
{
	// 1D derivative-of-Gaussian kernel (column vector), zero-sum
	// G'(x) = - (x / sigma^2) * G(x)
	CV_Assert(sigma > 0 && (ksize % 2 == 1));
	int half = ksize / 2;
	cv::Mat g = gaussian1D(sigma, ksize);
	cv::Mat dg(ksize, 1, CV_32F);

	float s2 = sigma * sigma;
	for (int i = -half; i <= half; ++i) {
		float val = -(static_cast<float>(i) / s2) * g.at<float>(i + half, 0);
		dg.at<float>(i + half, 0) = val;
	}

	// Optional: L2 normalize (keeps scale consistent across ¦Ò)
	double nrm = cv::norm(dg, cv::NORM_L2);
	if (nrm > 0) dg /= static_cast<float>(nrm);

	return dg;
}

inline int ksizeForSigma(float sigma)
{
	// Standard rule: 6*sigma rounded up to odd
	int k = static_cast<int>(std::ceil(6.f * sigma)) | 1;
	// Clamp to a reasonable minimum/maximum if you like:
	k = std::max(k, 3);
	return k;
}

inline void computeDoGGradients(const Mat& srcGray8U,
	float sigma,
	Mat& gradX32F,
	Mat& gradY32F,
	int borderType = cv::BORDER_REPLICATE)
{
	CV_Assert(srcGray8U.type() == CV_8UC1);
	CV_Assert(sigma > 0);

	// Convert to float [0,1]
	cv::Mat f32;
	srcGray8U.convertTo(f32, CV_32F, 1.0 / 255.0);

	// Build separable kernels
	int ksize = ksizeForSigma(sigma);
	cv::Mat G = gaussian1D(sigma, ksize);   // column vector
	cv::Mat dG = dgaussian1D(sigma, ksize);  // column vector (zero mean)

											 // Ix = dG(x) * G(y) * I
	cv::sepFilter2D(f32, gradX32F, CV_32F, dG, G, cv::Point(-1, -1), 0.0, borderType);

	// Iy = G(x) * dG(y) * I
	cv::sepFilter2D(f32, gradY32F, CV_32F, G, dG, cv::Point(-1, -1), 0.0, borderType);
}

// ---- Utility: boundary-safe pixel access & clamping ----
inline float pix32f(const Mat& img32f, int y, int x) {
	if (x < 0) x = 0; else if (x >= img32f.cols) x = img32f.cols - 1;
	if (y < 0) y = 0; else if (y >= img32f.rows) y = img32f.rows - 1;
	return img32f.at<float>(y, x);
}

inline double clamp255(double v) {
	return (v < 0.0) ? 0.0 : (v > 255.0 ? 255.0 : v);
}


// ---- 3¡Á3 Gaussian weights: supports anisotropy (sigmaX/Y adjustable; typically around 0.9) ----
inline void gaussianWeights3x3(double fx, double fy, double W[3][3],
	double sigmaX, double sigmaY)
{
	// Use the fractional part (fx, fy) of the target point (x, y) to define 3¡Á3 relative coordinates:
	// sample locations (i, j), where i, j ¡Ê {?1, 0, 1}
	const double inv2sx2 = 1.0 / (2.0 * sigmaX * sigmaX);
	const double inv2sy2 = 1.0 / (2.0 * sigmaY * sigmaY);

	for (int j = -1; j <= 1; ++j) {
		for (int i = -1; i <= 1; ++i) {
			double dx = (i - fx);
			double dy = (j - fy);
			W[j + 1][i + 1] = std::exp(-(dx * dx * inv2sx2 + dy * dy * inv2sy2));
		}
	}
}


// ---- Core: 3¡Á3 Gaussian-weighted local linear regression (LOESS-1) ----
// Fit a local linear model z = a + b*X + c*Y over a 3¡Á3 window using Gaussian weights,
// where X, Y are relative coordinates centered at the target pixel.
// Since we evaluate at X=Y=0, the output is simply a.
inline double interp3x3_gauss_LOESS1(const Mat& img32f, double x, double y,
	double sigmaX = 0.9, double sigmaY = 0.9)
{
	CV_Assert(img32f.type() == CV_32F && img32f.channels() == 1);

	int ix = (int)std::floor(x), iy = (int)std::floor(y);
	double fx = x - ix;  // [0,1)
	double fy = y - iy;

	// 3¡Á3 Gaussian weights
	double W[3][3];
	gaussianWeights3x3(fx, fy, W, sigmaX, sigmaY);

	// Construct weighted normal equations (3¡Á3): A^T W A c = A^T W b
	// Each row of A is [1, X, Y], where X = (ix+i)-x and Y = (iy+j)-y
	// Since we're evaluating at (0,0), the result is simply coefficient a
	double S00 = 0, S10 = 0, S01 = 0, S20 = 0, S11 = 0, S02 = 0; // Corresponding to ¡Æw, ¡ÆwX, ¡ÆwY, ¡ÆwX^2, ¡ÆwXY, ¡ÆwY^2
	double T0 = 0, T1 = 0, T2 = 0;                               // Corresponding to ¡Æw z, ¡Æw X z, ¡Æw Y z

	for (int j = -1; j <= 1; ++j) {
		for (int i = -1; i <= 1; ++i) {
			double w = W[j + 1][i + 1];
			double X = (ix + i) - x;
			double Y = (iy + j) - y;
			double z = (double)pix32f(img32f, iy + j, ix + i);

			S00 += w;
			S10 += w * X;
			S01 += w * Y;
			S20 += w * X * X;
			S11 += w * X * Y;
			S02 += w * Y * Y;

			T0 += w * z;
			T1 += w * X * z;
			T2 += w * Y * z;
		}
	}

	// Solve the 3¡Á3 linear system:
	// [S00 S10 S01][a] = [T0]
	// [S10 S20 S11][b]   [T1]
	// [S01 S11 S02][c]   [T2]
	cv::Matx33d A(S00, S10, S01,
		S10, S20, S11,
		S01, S11, S02);
	cv::Vec3d   B(T0, T1, T2), C;
	bool ok = cv::solve(A, B, C, cv::DECOMP_CHOLESKY);
	if (!ok) ok = cv::solve(A, B, C, cv::DECOMP_QR);

	double val;
	if (ok) {
		val = C[0]; // a (value at X=Y=0)
	}
	else {
		// Rare degenerate case: fall back to 3¡Á3 Gaussian-weighted average (zero-order kernel regression)
		double num = 0.0, den = 0.0;
		for (int j = -1; j <= 1; ++j) for (int i = -1; i <= 1; ++i) {
			double w = W[j + 1][i + 1];
			double z = (double)pix32f(img32f, iy + j, ix + i);
			num += w * z; den += w;
		}
		val = (den > 0 ? num / den : (double)pix32f(img32f, iy, ix));
	}
	return clamp255(val);
}


#endif