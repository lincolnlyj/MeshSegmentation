//结构体、常量、枚举类型的定义
#ifndef STRUCTSANDCONSTS_H
#define STRUCTSANDCONSTS_H

#include <vector>
using namespace std;

typedef double Coordinate[3];

const double DELTA = 0.5;//初始化种子的加权参数
const double CONVEXETA = 0.1;//凹面的Eta值
const double INF = 999999;//无穷
const double EPSILON = 0.036;//模糊分割的阈值

enum FaceState
{
	UNDISCOVERED, DISCOVERED, VISITIED
};

enum FaceKind
{
	FUZZY, DETERMINED
};


struct Edge
{
	int Start;
	int End;
	double Distance;
	friend bool operator < (const Edge& A, const Edge& B);
};

struct Neighbor
{
	int EdgeV[2];
	int Face;
	double Geod;//测地距离
	double AngDist;//角距离
	double Weight;//权重
};

struct Face
{
	int V[3];//组成面的三个顶点
	Coordinate Center;
	vector<Neighbor> Neighbors;
	double NormalVector[3];
	FaceState State;
	vector<double> Pl;//落在不同个区域中的概率
	FaceKind Kind;//面片的类型
	pair<int, int> PossibleRep;//落在的聚类种子
	int Parent;//父亲面片编号，用于最大流问题广度优先搜索
	double CurMaxCap;//到达当前节点的最大流
	Face();
};

struct Color
{
	int R;
	int G;
	int B;
	Color(int r = 0, int g = 0, int b = 0);
};

struct Vertice
{
	Coordinate Position;
	Color ColorData;
};

#endif