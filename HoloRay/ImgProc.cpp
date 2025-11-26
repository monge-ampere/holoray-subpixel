#include "stdafx.h"
#include "ImgProc.h"


ImgProc::ImgProc(const char *maskPath, const char *bdyPath, const char *oriPath)
{
	src_mat = imread(oriPath, IMREAD_GRAYSCALE);
	bdy_mat = imread(bdyPath, IMREAD_GRAYSCALE);
	mask_mat = imread(maskPath, IMREAD_GRAYSCALE);
	assert(src_mat.rows == mask_mat.rows);
	assert(src_mat.cols == mask_mat.cols);
	r = src_mat.rows;
	c = src_mat.cols;
	canny_mat.create(r, c, CV_8UC1);

	//ori_mat.create(r, c, CV_8UC3);

	//for (int i = 0; i != r; ++i)
	//{
	//	for (int j = 0; j != c; ++j)
	//	{
	//		uchar d = src_mat.at<uchar>(i, j);
	//		ori_mat.at<Vec3b>(i, j) = Vec3b(d, d, d);
	//	}
	//}
	ori_mat = imread(oriPath, IMREAD_COLOR);
	imwrite("E:\\gray.bmp", src_mat);
	//Sobel4Edge(sobel_mat, src_mat);
	DoG4Edge(sobel_mat, src_mat);
	imwrite("E:\\sobel.bmp", sobel_mat);
	double lowThreshold = 80;
	double highThreshold = 120;
	Canny(src_mat, canny_mat, lowThreshold, highThreshold);

	imwrite("E:\\canny.bmp", canny_mat);

	GetBdy();
}


ImgProc::~ImgProc()
{
}

bool ImgProc::GetFirstBdyPt(Mat &mat)
{
	for (int i = 0; i != r; ++i)
	{
		for (int j = 0; j != c; ++j)
		{
			uchar data = mat.at<uchar>(i, j);

			if (255 == data)
			{
				boundary.push_back(i*c + j);
				Point2d p;
				p.x = j + 0.5;
				p.y = i + 0.5;
				orderedBdy.push_back(p);
				mat.at<uchar>(i, j) = 0;
				return true;
			}
		}
	}

	return false;
}

void ImgProc::GetBdy()
{
	CV_Assert(bdy_mat.type() == CV_8UC1);

	orderedBdy.clear();
	boundary.clear();

	// 1) Work on a deep copy to avoid modifying the source mask.
	cv::Mat bw = bdy_mat.clone();

	// 2) Pad by 1 pixel to make neighbor checks safe (no bounds branches).
	cv::Mat img;
	cv::copyMakeBorder(bw, img, 1, 1, 1, 1, cv::BORDER_CONSTANT, 0);

	const int H = img.rows;
	const int W = img.cols;
	const int origCols = bdy_mat.cols;  // for linear indexing like (y * cols + x)

										// 3) Find the first foreground pixel as the starting point.
	cv::Point cur(-1, -1);
	for (int y = 1; y < H - 1 && cur.x < 0; ++y)
		for (int x = 1; x < W - 1; ++x)
			if (img.at<uchar>(y, x) == 255) { cur = { x, y }; break; }

	if (cur.x < 0) return; // No boundary found.

						   // 4) Moore-neighborhood directions (clockwise):
						   //    E, SE, S, SW, W, NW, N, NE
	static const int dx[8] = { +1, +1,  0, -1, -1, -1,  0, +1 };
	static const int dy[8] = { 0, +1, +1, +1,  0, -1, -1, -1 };

	// 5) Backtrack direction (the neighbor behind current). Initialize to "W".
	//    This means our search will begin from the neighbor after W (i.e., NW)
	//    which corresponds to a consistent ¡°hug the boundary¡± rule (left-hand).
	int prevIdx = 4; // W
	const cv::Point start = cur;
	const int startPrevIdx = prevIdx;

	// Safety bound to avoid infinite loops in malformed inputs.
	const int maxIter = (H - 2) * (W - 2) * 8;
	int iter = 0;

	do {
		// 6) Record current boundary pixel:
		//    - Convert from padded coords to original image coords by subtracting 1.
		//    - Store pixel center as (x + 0.5, y + 0.5).
		const int ox = cur.x - 1;
		const int oy = cur.y - 1;

		orderedBdy.emplace_back(ox + 0.5, oy + 0.5);
		boundary.push_back(oy * origCols + ox);

		// 7) Scan neighbors starting from the one after prevIdx, clockwise.
		//    This ensures we ¡°follow¡± the boundary and don¡¯t zig-zag or jump.
		int foundIdx = -1;
		for (int k = 1; k <= 8; ++k) {
			const int idx = (prevIdx + k) & 7; // modulo 8
			const int nx = cur.x + dx[idx];
			const int ny = cur.y + dy[idx];

			// Treat any non-zero as boundary; tighten to 255 if your mask is binary.
			if (img.at<uchar>(ny, nx) != 0) {
				foundIdx = idx;
				break;
			}
		}

		if (foundIdx < 0) {
			// Isolated point or broken mask; stop gracefully.
			break;
		}

		// 8) Optional: mark current as visited to prevent re-visiting in noisy masks.
		//    Moore tracing doesn¡¯t require this strictly, but it helps robustness.
		img.at<uchar>(cur.y, cur.x) = 0;

		// 9) Step to the found neighbor and update the backtrack direction
		//    to the opposite of foundIdx.
		cur.x += dx[foundIdx];
		cur.y += dy[foundIdx];
		prevIdx = (foundIdx + 4) & 7; // opposite direction

									  // 10) Stop when we have returned to the start *and* the backtrack
									  //     direction is also the original one. This guarantees a closed loop.
		if (++iter > maxIter) break;
	} while (!(cur == start && prevIdx == startPrevIdx));

	// If your downstream expects the first point repeated at the end,
	// uncomment the line below to explicitly close the polyline:
	// if (!orderedBdy.empty()) orderedBdy.push_back(orderedBdy.front());
}

double ImgProc::EuclidDist(double2 &v1, double2 &v2)
{
	return sqrt(pow(v1[0] - v2[0], 2)
		+ pow(v1[1] - v2[1], 2));
}

void ImgProc::DoG4Edge(Mat &dMat, Mat &srcMat)
{
	float sigma = 1.1f;

	cv::Mat gx, gy;
	computeDoGGradients(srcMat, sigma, gx, gy);

	cv::Mat mag, ang;
	cv::cartToPolar(gx, gy, mag, ang, false); // angle in radians

	cv::normalize(mag, mag, 0, 255, cv::NORM_MINMAX);
	mag.convertTo(dMat, CV_8U);

}

void ImgProc::Sobel4Edge(Mat &sMat, Mat &srcMat)
{
	cv::Mat blurred_image;
	cv::GaussianBlur(srcMat, blurred_image, cv::Size(7, 7), 2.5);
	cv::Mat sobel_x, sobel_y;
	cv::Sobel(blurred_image, sobel_x, CV_64F, 1, 0, 3);
	cv::Sobel(blurred_image, sobel_y, CV_64F, 0, 1, 3);

	cv::Mat sobel_magnitude;
	cv::magnitude(sobel_x, sobel_y, sobel_magnitude);
	cv::normalize(sobel_magnitude, sMat, 0, 255, cv::NORM_MINMAX);
	sMat.convertTo(sMat, CV_8U);
}

void ImgProc::InitKdTree4P()
{
	int M_num = ps.size();
	// copy model points to M_data
	M_data4P.resize(boost::extents[M_num][2]);
	for (int m = 0; m < M_num; m++)
	{
		M_data4P[m][0] = ps[m].x;
		M_data4P[m][1] = ps[m].y;
	}

	// build a kd tree from the model point cloud
	M_tree4P = new kdtree::KDTree(M_data4P);
}

void ImgProc::DelKdTree4P()
{
	delete M_tree4P;
}