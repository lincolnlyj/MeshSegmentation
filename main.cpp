#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <queue>
#include <functional>
#include <algorithm>
#include "StructsAndConsts.h"
#include "Functions.h"
using namespace std;


int main(int argc, char* argv[])
{
	string FileName;//Obj文件的名称
	vector<Vertice> Vertices;//所有顶点
	vector<Face> Faces;//所有面
	double** ppEdges;//任意两点间距离的矩阵，精确分割时为流量
	pair<double, double> AverageDist;//第一个是AverageGeod，第二个是AverageAngleDist
	vector<int> Reps;//聚类中心
	//读取文件
	if (!readObjFile(FileName, Vertices, Faces))
		return 0;

	unsigned int FacesCnt = Faces.size();//
	ppEdges = new double*[FacesCnt];
	for (unsigned int i = 0; i < FacesCnt; i++)
	{
		ppEdges[i] = new double[FacesCnt];
	}

	AverageDist = calculateDist(Vertices, Faces);
	dijkstra(Vertices, Faces, ppEdges);//得到初始聚类种子

	decideReps(Faces, ppEdges, Reps);

	cluster(Faces, Reps, ppEdges);//进行聚类迭代
	Div(Vertices, Faces, Reps, AverageDist.second, ppEdges, FileName);

	outputObjFile("New" + FileName, Vertices, Faces);

	for (unsigned int i = 0; i < FacesCnt; i++)//清除边的矩阵
	{
		delete[] ppEdges[i];
	}
	delete ppEdges;
	return 0;
}
