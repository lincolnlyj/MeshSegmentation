//函数定义
#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <queue>
#include <functional>
#include "StructsAndConsts.h"
using namespace std;

//分割字符串的函数
void splitString(const std::string& s, std::vector<std::string>& v, const std::string& c);
//读取obj文件
bool readObjFile(string& FileName, vector<Vertice>& Vertices, vector<Face>& Faces);
//输出obj文件
void outputObjFile(string NewFileName, const vector<Vertice>& Vertices, const vector<Face>& Faces);
//计算测地距离和角距离
pair<double, double> calculateDist(vector<Vertice>& Vertices, vector<Face>& Faces);
//计算任意两个面片间的距离
Edge dijkstra(vector<Vertice>& Vertices, vector<Face>& Faces, double**& ppEdges);
//清除所有标记
void clearState(vector<Face>& Faces);
//聚类
void cluster(vector<Face>& Faces, int& Repa, int& Repb, double**& ppEdges);
//模糊划分
void fuzzyDiv(vector<Vertice>& Vertices, vector<Face>& Faces, vector<int>& FuzzyFacesIdx);
//精细划分，利用最大流进行划分
void accurateDiv(vector<Vertice>& Vertices, vector<Face>& Faces, vector<int>& FuzzyFacesIdx, double AverageAngleDist, double**& ppEdges, int Repa, int Repb);
//清除父亲面片编号
void clearParent(vector<Face>& Faces, vector<int>& FacesIdx);
//求绝对值
inline double absd(double a)
{
	if (a > 0)
		return a;
	else
		return -a;
}


#endif
