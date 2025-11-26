#pragma once
/*!
*      \file FileProc.h
*      \brief read and write files
*	   \author Wei Feng
*      \date 11/24/2024
*/
#include <direct.h>
#include <iostream>
#include <fstream>
#include <assert.h>
#include <list>
#include <vector>
#include "Parser\strutil.h"

using namespace std;

inline void write2m(double(*taget)[2],
	int pN, uchar *gray, int gN, int(*tri_set)[3], int tN, string path)
{
	assert(pN == gN);
	FILE *fp = fopen(path.c_str(), "w+");
	if (fp == NULL) {
		fprintf(stderr, "fopen() failed.\n");
		exit(EXIT_FAILURE);
	}

	for (int i = 0; i != pN; ++i)
	{
		double g01 = gray[i] / 255.0;
		fprintf(fp, "Vertex %d %f %f %f {rgb=(%f %f %f) normal=(0 0 1.0) father=(%d)}\n"
			, i, taget[i][0], taget[i][1], 0.0, g01, g01, g01, i);
	}

	for (int i = 0; i != tN; ++i)
	{
		fprintf(fp, "Face %d %d %d %d\n"
			, i + 1, tri_set[i][0], tri_set[i][1], tri_set[i][2]);
	}


	fclose(fp);
}

inline void write2f(double(*taget)[2], int pN, string path)
{
	FILE *fp = fopen(path.c_str(), "w+");
	if (fp == NULL) {
		fprintf(stderr, "fopen() failed.\n");
		exit(EXIT_FAILURE);
	}

	for (int i = 0; i != pN; ++i)
	{
		fprintf(fp, "%f,%f\n", taget[i][0], taget[i][1]);
	}

	fclose(fp);
}

inline void write2f(vector<Point2d> target, string path)
{
	int pN = target.size();
	FILE *fp = fopen(path.c_str(), "w+");
	if (fp == NULL) {
		fprintf(stderr, "fopen() failed.\n");
		exit(EXIT_FAILURE);
	}

	for (int i = 0; i != pN; ++i)
	{
		fprintf(fp, "%f,%f\n", target[i].x, target[i].y);
	}

	fclose(fp);
}

inline void writeBdy2f(int *taget, int bN, string path)
{
	FILE *fp = fopen(path.c_str(), "w+");
	if (fp == NULL) {
		fprintf(stderr, "fopen() failed.\n");
		exit(EXIT_FAILURE);
	}
	int pre_i;
	for (int i = 0; i != bN; ++i)
	{
		pre_i = (bN + i - 1) % bN;
		fprintf(fp, "%d,%d\n"
			, taget[pre_i], taget[i]);
	}

	fclose(fp);
}

inline void readm(Pt *pts, vector<int> &ids, vector<vector<int>> &tris, const char* path)
{
	std::fstream is(path, std::fstream::in);

	if (is.fail())
	{
		fprintf(stderr, "Error in opening file %s\n", path);
		return;
	}

	char buffer[_MAX_PATH];
	std::string token;
	int id = 0;
	while (is.getline(buffer, _MAX_PATH))
	{
		std::string line(buffer);
		line = strutil::trim(line);

		strutil::Tokenizer stokenizer(line, " \r\n");

		stokenizer.nextToken();
		std::string token = stokenizer.getToken();

		if (token == "Vertex")
		{
			Pt p;
			stokenizer.nextToken();
			token = stokenizer.getToken();
			id = strutil::parseString<int>(token);

			stokenizer.nextToken();
			token = stokenizer.getToken();
			p.oci = strutil::parseString<double>(token);

			stokenizer.nextToken();
			token = stokenizer.getToken();
			p.ori = strutil::parseString<double>(token);

			stokenizer.nextToken();
			token = stokenizer.getToken();
			p.x = strutil::parseString<double>(token);

			stokenizer.nextToken();
			token = stokenizer.getToken();
			p.y = strutil::parseString<double>(token);

			pts[id] = p;
			ids.push_back(id);
			continue;
		}
		else if (token == "Face")
		{
			vector<int> tri;
			stokenizer.nextToken();
			token = stokenizer.getToken();
			id = strutil::parseString<int>(token);

			stokenizer.nextToken();
			token = stokenizer.getToken();
			int i1 = strutil::parseString<int>(token);
			tri.push_back(i1);

			stokenizer.nextToken();
			token = stokenizer.getToken();
			int i2 = strutil::parseString<int>(token);
			tri.push_back(i2);

			stokenizer.nextToken();
			token = stokenizer.getToken();
			int i3 = strutil::parseString<int>(token);
			tri.push_back(i3);

			tris.push_back(tri);
		}
	}
}