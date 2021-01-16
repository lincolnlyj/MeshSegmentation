//结构体、常量、枚举类型的定义
#ifndef STRUCTSANDCONSTS_H
#define STRUCTSANDCONSTS_H

#include <vector>
using namespace std;

typedef double Coordinate[3];

const double DELTA = 0.6;//初始化种子的加权参数
const double CONVEXETA = 0.1;//凹面的Eta值
const double INF = 999999;//无穷
const double EPSILON = 0.1;//模糊分割的阈值

enum FaceState
{
	UNDISCOVERED, DISCOVERED, VISITIED
};

enum FaceClass
{
	A, B, FUZZY
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
	double Pa;//落在两个区域中的概率
	double Pb;
	FaceClass Kind;//面片的类型
	int Parent;//父亲面片编号，用于最大流问题广度优先搜索
	double CurMaxCap;//到达当前节点的最大流
	Face();
};

struct Color
{
	int R;
	int G;
	int B;
};

struct Vertice
{
	Coordinate Position;
	Color ColorData;
	Vertice();
};


#endif
