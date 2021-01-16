#include "Functions.h"

void splitString(const std::string& s, std::vector<std::string>& v, const std::string& c)
{
	std::string::size_type pos1, pos2;
	pos2 = s.find(c);
	pos1 = 0;
	while (std::string::npos != pos2)
	{
		v.push_back(s.substr(pos1, pos2 - pos1));

		pos1 = pos2 + c.size();
		pos2 = s.find(c, pos1);
	}
	if (pos1 != s.length())
		v.push_back(s.substr(pos1));
}


bool readObjFile(string& FileName, vector<Vertice>& Vertices, vector<Face>& Faces)
{
	cin >> FileName;
	ifstream ObjFile(FileName.c_str());//打开Obj文件
	if (!ObjFile)
	{
		cout << "打开文件失败" << endl;
		return false;
	}
	
	string Temp;
	vector<string> VecTemp;
	map<pair<int, int>, int> Edges;//用于得到相邻的边

	while (getline(ObjFile, Temp))
	{
		VecTemp.clear();
		splitString(Temp, VecTemp, " ");
		if (!VecTemp.empty())
		{

			if (VecTemp[0] == "v")//如果输入的是顶点
			{
				Vertice TempV;
				if (VecTemp.size() == 4)
				{
					TempV.Position[0] = stod(VecTemp[1]);
					TempV.Position[1] = stod(VecTemp[2]);
					TempV.Position[2] = stod(VecTemp[3]);
					Vertices.push_back(TempV);
				}
				else if (VecTemp.size() == 7)
				{
					TempV.Position[0] = stod(VecTemp[1]);
					TempV.Position[1] = stod(VecTemp[2]);
					TempV.Position[2] = stod(VecTemp[3]);
					TempV.ColorData.R = stoi(VecTemp[4]);
					TempV.ColorData.G = stoi(VecTemp[5]);
					TempV.ColorData.B = stoi(VecTemp[6]);
					Vertices.push_back(TempV);
				}
				else
				{
					cout << "文件格式有误" << endl;
					ObjFile.close();
					return false;
				}
			}
			else if (VecTemp[0] == "f")//如果输入的是面片
			{
				Face TempFace;
				TempFace.V[0] = stoi(VecTemp[1]) - 1;
				TempFace.V[1] = stoi(VecTemp[2]) - 1;
				TempFace.V[2] = stoi(VecTemp[3]) - 1;
				//计算质心
				TempFace.Center[0] =
					(Vertices[TempFace.V[0]].Position[0]
						+ Vertices[TempFace.V[1]].Position[0]
						+ Vertices[TempFace.V[2]].Position[0]) / 3;
				TempFace.Center[1] =
					(Vertices[TempFace.V[0]].Position[1]
						+ Vertices[TempFace.V[1]].Position[1]
						+ Vertices[TempFace.V[2]].Position[1]) / 3;
				TempFace.Center[2] =
					(Vertices[TempFace.V[0]].Position[2]
						+ Vertices[TempFace.V[1]].Position[2]
						+ Vertices[TempFace.V[2]].Position[2]) / 3;
				//计算法向量
				TempFace.NormalVector[0] = (Vertices[TempFace.V[1]].Position[1] - Vertices[TempFace.V[0]].Position[1]) * (Vertices[TempFace.V[2]].Position[2] - Vertices[TempFace.V[1]].Position[2])
					- (Vertices[TempFace.V[2]].Position[1] - Vertices[TempFace.V[1]].Position[1]) * (Vertices[TempFace.V[1]].Position[2] - Vertices[TempFace.V[0]].Position[2]);
				TempFace.NormalVector[1] = (Vertices[TempFace.V[1]].Position[2] - Vertices[TempFace.V[0]].Position[2]) * (Vertices[TempFace.V[2]].Position[0] - Vertices[TempFace.V[1]].Position[0])
					- (Vertices[TempFace.V[2]].Position[2] - Vertices[TempFace.V[1]].Position[2]) * (Vertices[TempFace.V[1]].Position[0] - Vertices[TempFace.V[0]].Position[0]);
				TempFace.NormalVector[2] = (Vertices[TempFace.V[1]].Position[0] - Vertices[TempFace.V[0]].Position[0]) * (Vertices[TempFace.V[2]].Position[1] - Vertices[TempFace.V[1]].Position[1])
					- (Vertices[TempFace.V[2]].Position[0] - Vertices[TempFace.V[1]].Position[0]) * (Vertices[TempFace.V[1]].Position[1] - Vertices[TempFace.V[0]].Position[1]);
				double VectorLength = sqrt(pow(TempFace.NormalVector[0], 2) + pow(TempFace.NormalVector[1], 2) + pow(TempFace.NormalVector[2], 2));
				//向量单位化
				TempFace.NormalVector[0] = TempFace.NormalVector[0] / VectorLength;
				TempFace.NormalVector[1] = TempFace.NormalVector[1] / VectorLength;
				TempFace.NormalVector[2] = TempFace.NormalVector[2] / VectorLength;

				//寻找相邻的面
				if (TempFace.V[0] < TempFace.V[1])
				{
					int& Temp = Edges[pair<int, int>(TempFace.V[0], TempFace.V[1])];
					if (Temp)//如果该边存储过面的信息，则记录为相邻边
					{
						Neighbor TempNeighbor;
						TempNeighbor.EdgeV[0] = TempFace.V[0];
						TempNeighbor.EdgeV[1] = TempFace.V[1];
						TempNeighbor.Face = Temp - 1;
						TempFace.Neighbors.push_back(TempNeighbor);
						TempNeighbor.Face = Faces.size();
						Faces[Temp - 1].Neighbors.push_back(TempNeighbor);
					}
					else
					{
						Temp = Faces.size() + 1;
					}
				}
				else
				{
					int& Temp = Edges[pair<int, int>(TempFace.V[1], TempFace.V[0])];
					if (Temp)
					{
						Neighbor TempNeighbor;
						TempNeighbor.EdgeV[0] = TempFace.V[0];
						TempNeighbor.EdgeV[1] = TempFace.V[1];
						TempNeighbor.Face = Temp - 1;
						TempFace.Neighbors.push_back(TempNeighbor);
						TempNeighbor.Face = Faces.size();
						Faces[Temp - 1].Neighbors.push_back(TempNeighbor);
					}
					else
					{
						Temp = Faces.size() + 1;
					}
				}
				if (TempFace.V[0] < TempFace.V[2])
				{
					int& Temp = Edges[pair<int, int>(TempFace.V[0], TempFace.V[2])];
					if (Temp)//如果该边存储过面的信息，则记录为相邻边
					{
						Neighbor TempNeighbor;
						TempNeighbor.EdgeV[0] = TempFace.V[0];
						TempNeighbor.EdgeV[1] = TempFace.V[2];
						TempNeighbor.Face = Temp - 1;
						TempFace.Neighbors.push_back(TempNeighbor);
						TempNeighbor.Face = Faces.size();
						Faces[Temp - 1].Neighbors.push_back(TempNeighbor);
					}
					else
					{
						Temp = Faces.size() + 1;
					}
				}
				else
				{
					int& Temp = Edges[pair<int, int>(TempFace.V[2], TempFace.V[0])];
					if (Temp)
					{
						Neighbor TempNeighbor;
						TempNeighbor.EdgeV[0] = TempFace.V[0];
						TempNeighbor.EdgeV[1] = TempFace.V[2];
						TempNeighbor.Face = Temp - 1;
						TempFace.Neighbors.push_back(TempNeighbor);
						TempNeighbor.Face = Faces.size();
						Faces[Temp - 1].Neighbors.push_back(TempNeighbor);
					}
					else
					{
						Temp = Faces.size() + 1;
					}
				}
				if (TempFace.V[1] < TempFace.V[2])
				{
					int& Temp = Edges[pair<int, int>(TempFace.V[1], TempFace.V[2])];
					if (Temp)//如果该边存储过面的信息，则记录为相邻边
					{
						Neighbor TempNeighbor;
						TempNeighbor.EdgeV[0] = TempFace.V[1];
						TempNeighbor.EdgeV[1] = TempFace.V[2];
						TempNeighbor.Face = Temp - 1;
						TempFace.Neighbors.push_back(TempNeighbor);
						TempNeighbor.Face = Faces.size();
						Faces[Temp - 1].Neighbors.push_back(TempNeighbor);
					}
					else
					{
						Temp = Faces.size() + 1;
					}
				}
				else
				{
					int& Temp = Edges[pair<int, int>(TempFace.V[2], TempFace.V[1])];
					if (Temp)
					{
						Neighbor TempNeighbor;
						TempNeighbor.EdgeV[0] = TempFace.V[1];
						TempNeighbor.EdgeV[1] = TempFace.V[2];
						TempNeighbor.Face = Temp - 1;
						TempFace.Neighbors.push_back(TempNeighbor);
						TempNeighbor.Face = Faces.size();
						Faces[Temp - 1].Neighbors.push_back(TempNeighbor);
					}
					else
					{
						Temp = Faces.size() + 1;
					}
				}

				Faces.push_back(TempFace);
			}
		}
	}

	ObjFile.close();
	return true;
}

void outputObjFile(string NewFileName, const vector<Vertice>& Vertices, const vector<Face>& Faces)
{
	ofstream OutputObj(NewFileName.c_str());
	OutputObj << "# " << Vertices.size() << " vertices, " << Faces.size() << " faces" << endl;

	unsigned int VerticesCnt = Vertices.size();
	unsigned int FacesCnt = Faces.size();

	//输出点的信息
	for (unsigned int i = 0; i < VerticesCnt; i++)
	{
		OutputObj << "v " << Vertices[i].Position[0] << ' '
			<< Vertices[i].Position[1] << ' '
			<< Vertices[i].Position[2] << ' '
			<< Vertices[i].ColorData.R << ' '
			<< Vertices[i].ColorData.G << ' '
			<< Vertices[i].ColorData.B << endl;
	}

	//输出面的信息
	for (unsigned int i = 0; i < FacesCnt; i++)
	{
		OutputObj << "f " << Faces[i].V[0] + 1 << ' '
			<< Faces[i].V[1] + 1 << ' '
			<< Faces[i].V[2] + 1 << endl;
	}
	OutputObj.close();
}

pair<double, double> calculateDist(vector<Vertice>& Vertices, vector<Face>& Faces)
{
	unsigned int FacesCnt = Faces.size();
	double AverageGeod = 0;
	double AverageAngleDist = 0;
	unsigned int Cnt = 0;
	for (unsigned int i = 0; i < FacesCnt; i++)
	{
		for (unsigned int j = 0; j < Faces[i].Neighbors.size(); j++)
		{
			++Cnt;
			//计算测地距离
			int NeiFace;//相邻边的序号
			NeiFace = Faces[i].Neighbors[j].Face;
			Coordinate Direction;//公共边方向向量
			Direction[0] = Vertices[Faces[i].Neighbors[j].EdgeV[1]].Position[0] - Vertices[Faces[i].Neighbors[j].EdgeV[0]].Position[0];
			Direction[1] = Vertices[Faces[i].Neighbors[j].EdgeV[1]].Position[1] - Vertices[Faces[i].Neighbors[j].EdgeV[0]].Position[1];
			Direction[2] = Vertices[Faces[i].Neighbors[j].EdgeV[1]].Position[2] - Vertices[Faces[i].Neighbors[j].EdgeV[0]].Position[2];
			double DirLength;//方向向量长度
			DirLength = sqrt(pow(Direction[0], 2) + pow(Direction[1], 2) + pow(Direction[2], 2));
			Direction[0] /= DirLength;//向量单位化
			Direction[1] /= DirLength;
			Direction[2] /= DirLength;
			Coordinate L1;//顶点1到质心的方向向量
			L1[0] = Faces[i].Center[0] - Vertices[Faces[i].Neighbors[j].EdgeV[0]].Position[0];
			L1[1] = Faces[i].Center[1] - Vertices[Faces[i].Neighbors[j].EdgeV[0]].Position[1];
			L1[2] = Faces[i].Center[2] - Vertices[Faces[i].Neighbors[j].EdgeV[0]].Position[2];
			double L1Length;//L1的长度
			L1Length = sqrt(pow(L1[0], 2) + pow(L1[1], 2) + pow(L1[2], 2));
			L1[0] /= L1Length;//向量单位化
			L1[1] /= L1Length;
			L1[2] /= L1Length;
			Coordinate L2;//顶点2到质心的方向向量
			L2[0] = Faces[NeiFace].Center[0] - Vertices[Faces[i].Neighbors[j].EdgeV[0]].Position[0];
			L2[1] = Faces[NeiFace].Center[1] - Vertices[Faces[i].Neighbors[j].EdgeV[0]].Position[1];
			L2[2] = Faces[NeiFace].Center[2] - Vertices[Faces[i].Neighbors[j].EdgeV[0]].Position[2];
			double L2Length;//L2的长度
			L2Length = sqrt(pow(L2[0], 2) + pow(L2[1], 2) + pow(L2[2], 2));
			L2[0] /= L2Length;//向量单位化
			L2[1] /= L2Length;
			L2[2] /= L2Length;
			double Alpha = acos(Direction[0] * L1[0] + Direction[1] * L1[1] + Direction[2] * L1[2]);//计算两个方向向量的夹角
			double Beta = acos(Direction[0] * L2[0] + Direction[1] * L2[1] + Direction[2] * L2[2]);
			//展平后计算测地距离
			Faces[i].Neighbors[j].Geod = sqrt(pow(L1Length, 2) + pow(L2Length, 2) - 2 * L1Length * L2Length * cos(Alpha + Beta));
			AverageGeod += Faces[i].Neighbors[j].Geod;

			//计算角距离
			double CosAngle = Faces[i].NormalVector[0] * Faces[NeiFace].NormalVector[0] + Faces[i].NormalVector[1] * Faces[NeiFace].NormalVector[1] + Faces[i].NormalVector[2] * Faces[NeiFace].NormalVector[2];
			double CenterVec[3];//两个质心之间的向量
			CenterVec[0] = Faces[NeiFace].Center[0] - Faces[i].Center[0];
			CenterVec[1] = Faces[NeiFace].Center[1] - Faces[i].Center[1];
			CenterVec[2] = Faces[NeiFace].Center[2] - Faces[i].Center[2];
			double VecLength = sqrt(pow(CenterVec[0], 2) + pow(CenterVec[1], 2) + pow(CenterVec[2], 2));
			CenterVec[0] /= VecLength;//向量单位化
			CenterVec[1] /= VecLength;
			CenterVec[2] /= VecLength;
			double TempCos;//用于判断凹凸的角，大于零为凹面，小于零为凸面
			TempCos = Faces[i].NormalVector[0] * CenterVec[0] + Faces[i].NormalVector[1] * CenterVec[1] + Faces[i].NormalVector[2] * CenterVec[2];
			double Eta;
			if (TempCos < 0)
			{
				Eta = CONVEXETA;
			}
			else
			{
				Eta = 1;
			}
			Faces[i].Neighbors[j].AngDist = Eta * (1 - CosAngle);
			if (Faces[i].Neighbors[j].AngDist < 0)
				Faces[i].Neighbors[j].AngDist = 0;
			AverageAngleDist += Faces[i].Neighbors[j].AngDist;
			//cout << Angle << ' ' << AverageAngleDist << ' ' << i << endl;
		}
	}
	AverageGeod /= Cnt;
	AverageAngleDist /= Cnt;
	for (unsigned int i = 0; i < FacesCnt; i++)
	{
		for (unsigned int j = 0; j < Faces[i].Neighbors.size(); j++)
		{
			Faces[i].Neighbors[j].Weight = DELTA * (Faces[i].Neighbors[j].Geod / AverageGeod) + (1 - DELTA) * (Faces[i].Neighbors[j].AngDist / AverageAngleDist);
		}
	}
	return pair<double, double>(AverageGeod, AverageAngleDist);
}

Edge dijkstra(vector<Vertice>& Vertices, vector<Face>& Faces, double**& ppEdges)
{
	unsigned int FacesCnt = Faces.size();
	Edge LargestEdge;//最长的边
	LargestEdge.Distance = 0;
	priority_queue<Edge> EdgeQ;
	for (unsigned int i = 0; i < FacesCnt; i++)
	{
		Edge E;
		E.Start = i;
		E.End = i;
		E.Distance = 0;
		EdgeQ.push(E);
		Faces[i].State = DISCOVERED;
		while (!EdgeQ.empty())
		{
			E = EdgeQ.top();
			EdgeQ.pop();
			if (Faces[E.End].State == VISITIED)
				continue;
			Faces[E.End].State = VISITIED;
			ppEdges[E.Start][E.End] = E.Distance;
			if (E.Distance > LargestEdge.Distance)
				LargestEdge = E;
			int EdgeCnt = Faces[E.End].Neighbors.size();
			for (int j = 0; j < EdgeCnt; j++)
			{
				Edge Temp;
				Temp.Start = i;
				Temp.End = Faces[E.End].Neighbors[j].Face;
				Temp.Distance = Faces[E.End].Neighbors[j].Weight + E.Distance;
				if (Faces[Temp.End].State != VISITIED)
				{
					if (!ppEdges[Temp.Start][Temp.End])//如果没有存过临时距离，就保存临时距离
					{
						EdgeQ.push(Temp);
						ppEdges[Temp.Start][Temp.End] = Temp.Distance;
					}
					else if (Temp.Distance < ppEdges[Temp.Start][Temp.End])//如果存过临时距离，就判断新距离是否小于临时距离，如果小于再入队
					{
						EdgeQ.push(Temp);
						ppEdges[Temp.Start][Temp.End] = Temp.Distance;
					}
				}
			}
		}
		clearState(Faces);
	}
	return LargestEdge;
}

void clearState(vector<Face>& Faces)
{
	unsigned int FacesCnt = Faces.size();
	for (unsigned int i = 0; i < FacesCnt; i++)
	{
		Faces[i].State = UNDISCOVERED;
	}
}

void cluster(vector<Face>& Faces, int & Repa, int & Repb, double**& ppEdges)
{
	unsigned int FacesCnt = Faces.size();
	int LastRepa = Repa;//上一次的聚类种子
	int LastRepb = Repb;
	cout << "初始聚类种子: " << Repa << ' ' << Repb << endl;

	pair<int, double> SmallestRepa(Repa, 0);//新的聚类种子
	pair<int, double> SmallestRepb(Repb, 0);
	//初始化新聚类种子
	for (unsigned int i = 0; i < FacesCnt; i++)
	{
		Faces[i].Pa = ppEdges[i][Repb] / (ppEdges[i][Repa] + ppEdges[i][Repb]);
		Faces[i].Pb = ppEdges[i][Repa] / (ppEdges[i][Repa] + ppEdges[i][Repb]);
	}
	for (unsigned int j = 0; j < FacesCnt; j++)
	{
		SmallestRepa.second += Faces[j].Pa * ppEdges[Repa][j];
		SmallestRepb.second += Faces[j].Pb * ppEdges[Repb][j];
	}
	cout << "初始加权距离: " << SmallestRepa.second << ' ' << SmallestRepb.second << endl;

	int Cnt = 0;
	do
	{
		Cnt++;
		LastRepa = Repa;
		LastRepb = Repb;

		for (unsigned int i = 0; i < FacesCnt; i++)
		{
			double TempSumA = 0;//临时的加权和
			double TempSumB = 0;
			for (unsigned int j = 0; j < FacesCnt; j++)
			{
				TempSumA += Faces[j].Pa * ppEdges[i][j];
				if (TempSumA > SmallestRepa.second)
					break;
			}
			//更新聚类种子
			if (TempSumA < SmallestRepa.second)
			{
				SmallestRepa.first = i;
				SmallestRepa.second = TempSumA;
			}
			for (unsigned int j = 0; j < FacesCnt; j++)
			{
				TempSumB += Faces[j].Pb * ppEdges[i][j];
				if (TempSumB > SmallestRepb.second)
					break;
			}
			if (TempSumB < SmallestRepb.second)
			{
				SmallestRepb.first = i;
				SmallestRepb.second = TempSumB;
			}
		}
		Repa = SmallestRepa.first;
		Repb = SmallestRepb.first;
		if (Repa == LastRepa && Repb == LastRepb)//如果聚类种子没变就直接跳出循环
			break;

		cout << "第" << Cnt << "迭代: " << endl;
		cout << "聚类种子: " << Repa << ' ' << Repb << endl;
		cout << "加权距离: " << SmallestRepa.second << ' ' << SmallestRepb.second << endl;
		for (unsigned int i = 0; i < FacesCnt; i++)
		{
			Faces[i].Pa = ppEdges[i][Repb] / (ppEdges[i][Repa] + ppEdges[i][Repb]);
			Faces[i].Pb = ppEdges[i][Repa] / (ppEdges[i][Repa] + ppEdges[i][Repb]);
		}

	} while (Cnt < 20);
}

void fuzzyDiv(vector<Vertice>& Vertices, vector<Face>& Faces, vector<int>& FuzzyFacesIdx)
{
	unsigned int FacesCnt = Faces.size();
	for (unsigned int i = 0; i < FacesCnt; i++)
	{
		if (Faces[i].Pa > 0.5 + EPSILON)
		{
			Faces[i].Kind = A;
			for (unsigned int j = 0; j < 3; j++)
			{
				Vertices[Faces[i].V[j]].ColorData.R = 255;
				Vertices[Faces[i].V[j]].ColorData.G = 0;
				Vertices[Faces[i].V[j]].ColorData.B = 0;
			}
		}
		else if (Faces[i].Pb > 0.5 + EPSILON)
		{
			Faces[i].Kind = B;
			for (unsigned int j = 0; j < 3; j++)
			{
				Vertices[Faces[i].V[j]].ColorData.R = 0;
				Vertices[Faces[i].V[j]].ColorData.G = 255;
				Vertices[Faces[i].V[j]].ColorData.B = 0;
			}
		}
		else if (Faces[i].Pb <= 0.5 + EPSILON && Faces[i].Pb >= 0.5 - EPSILON)
		{
			Faces[i].Kind = FUZZY;
			for (unsigned int j = 0; j < 3; j++)
			{
				Vertices[Faces[i].V[j]].ColorData.R = 0;
				Vertices[Faces[i].V[j]].ColorData.G = 0;
				Vertices[Faces[i].V[j]].ColorData.B = 255;
			}
			FuzzyFacesIdx.push_back(i);
		}
	}
	unsigned int FuzzyFacesCnt = FuzzyFacesIdx.size();
	for (unsigned int i = 0; i < FuzzyFacesCnt; i++)
	{
		for (unsigned int j = 0; j < Faces[FuzzyFacesIdx[i]].Neighbors.size(); j++)
		{
			if (Faces[Faces[FuzzyFacesIdx[i]].Neighbors[j].Face].Kind != FUZZY)
			{
				Faces[Faces[FuzzyFacesIdx[i]].Neighbors[j].Face].Kind = FUZZY;
				FuzzyFacesIdx.push_back(Faces[FuzzyFacesIdx[i]].Neighbors[j].Face);//将边界上与模糊区域相交的面片也变成边界面
			}
		}
	}
}

void accurateDiv(vector<Vertice>& Vertices, vector<Face>& Faces, vector<int>& FuzzyFacesIdx, double AverageAngleDist, double**& ppEdges, int Repa, int Repb)
{
	vector<int> ANeighbors;
	//计算每条边的流量
	unsigned int FuzzyFacesCnt = FuzzyFacesIdx.size();
	for (unsigned int i = 0; i < FuzzyFacesCnt; i++)
	{
		for (unsigned int j = 0; j < Faces[FuzzyFacesIdx[i]].Neighbors.size(); j++)
		{
			if (Faces[Faces[FuzzyFacesIdx[i]].Neighbors[j].Face].Kind != FUZZY)//如果边上是A或B区域的面片，则流量为无穷
			{
				if (Faces[Faces[FuzzyFacesIdx[i]].Neighbors[j].Face].Kind == A)//如果边上是A区域面片，就把面片放入A区域的邻域中
				{
					ppEdges[FuzzyFacesIdx[i]][Repa] = INF;
					ppEdges[Repa][FuzzyFacesIdx[i]] = INF;
					ANeighbors.push_back(FuzzyFacesIdx[i]);
				}
				else if (Faces[Faces[FuzzyFacesIdx[i]].Neighbors[j].Face].Kind == B)//如果边上是B区域面片
				{
					ppEdges[FuzzyFacesIdx[i]][Repb] = INF;
					ppEdges[Repb][FuzzyFacesIdx[i]] = INF;
				}
			}
			else//其他情况的流量根据公式计算
			{
				ppEdges[FuzzyFacesIdx[i]][Faces[FuzzyFacesIdx[i]].Neighbors[j].Face] = 1 / (1 + Faces[FuzzyFacesIdx[i]].Neighbors[j].AngDist / AverageAngleDist);
			}
		}
	}

	bool PathExist = false;
	queue<int> NeighborQ;

	do
	{
		PathExist = false;
		clearState(Faces);
		clearParent(Faces, FuzzyFacesIdx);

		//广度优先搜索
		double CapMax = 0;//最大流量
		int EndFace = 0;
		for (unsigned int i = 0; i < ANeighbors.size(); i++)
		{
			if (ppEdges[Repa][ANeighbors[i]] <= 1E-6)//如果流量为0就不入队
				continue;
			Faces[ANeighbors[i]].State = DISCOVERED;
			Faces[ANeighbors[i]].Parent = Repa;
			NeighborQ.push(ANeighbors[i]);//从A区域开始进行广度优先遍历
			Faces[ANeighbors[i]].CurMaxCap = ppEdges[Repa][ANeighbors[i]];
		}
		while (!NeighborQ.empty())
		{
			int CurFace = NeighborQ.front();
			NeighborQ.pop();
			for (unsigned int j = 0; j < Faces[CurFace].Neighbors.size(); j++)
			{
				if (Faces[Faces[CurFace].Neighbors[j].Face].State == DISCOVERED)
					continue;
				if (absd(ppEdges[CurFace][Faces[CurFace].Neighbors[j].Face]) <= 1E-6)//如果到下一个面片的流量为0，下一个面片不入队
					continue;
				if (Faces[Faces[CurFace].Neighbors[j].Face].Kind == B)//如果下一个面片是在B区域内，找到通路
				{
					if (Faces[CurFace].CurMaxCap > CapMax)
					{
						PathExist = true;
						CapMax = Faces[CurFace].CurMaxCap;//设置最大流为当前面片中存储的最大流
						EndFace = CurFace;//当前面片为最后结束的面片
					}
				}
				if (Faces[Faces[CurFace].Neighbors[j].Face].Kind == FUZZY)//如果没有到达B区域，就入队，并记录父节点位置
				{
					Faces[Faces[CurFace].Neighbors[j].Face].Parent = CurFace;
					Faces[Faces[CurFace].Neighbors[j].Face].State = DISCOVERED;
					NeighborQ.push(Faces[CurFace].Neighbors[j].Face);
					if (ppEdges[CurFace][Faces[CurFace].Neighbors[j].Face] < Faces[CurFace].CurMaxCap)//如果到下一个面片的流量小于当前面片的最大流，就设置到下一个面片的流量为最大流
					{
						Faces[Faces[CurFace].Neighbors[j].Face].CurMaxCap = ppEdges[CurFace][Faces[CurFace].Neighbors[j].Face];
					}
					else
					{
						Faces[Faces[CurFace].Neighbors[j].Face].CurMaxCap = Faces[CurFace].CurMaxCap;
					}
				}
			}
		}

		if (PathExist)
		{
			do
			{
				ppEdges[Faces[EndFace].Parent][EndFace] -= CapMax;//减去当前通路上的最大流
				//ppEdges[EndFace][Faces[EndFace].Parent] -= CapMax;
				EndFace = Faces[EndFace].Parent;//继续回溯
			} while (EndFace != Repa);//如果没有回溯到起点，就继续回溯
		}
	} while (PathExist);

	clearState(Faces);
	//精细分割上色
	for (unsigned int i = 0; i < ANeighbors.size(); i++)
	{
		NeighborQ.push(ANeighbors[i]);//从A区域开始进行广度优先遍历
		Faces[ANeighbors[i]].Kind = A;
		Faces[ANeighbors[i]].State = DISCOVERED;
	}
	while (!NeighborQ.empty())
	{
		int CurFace = NeighborQ.front();
		NeighborQ.pop();
		for (unsigned int j = 0; j < 3; j++)
		{
			Vertices[Faces[CurFace].V[j]].ColorData.R = 255;
			Vertices[Faces[CurFace].V[j]].ColorData.G = 0;
			Vertices[Faces[CurFace].V[j]].ColorData.B = 0;
		}
		for (unsigned int j = 0; j < Faces[CurFace].Neighbors.size(); j++)
		{
			if (Faces[Faces[CurFace].Neighbors[j].Face].State == DISCOVERED)
				continue;
			if (Faces[Faces[CurFace].Neighbors[j].Face].Kind != FUZZY)//如果下一个面片是在不在模糊区域内，不入队
				continue;
			if (absd(ppEdges[CurFace][Faces[CurFace].Neighbors[j].Face]) <= 1E-6)//如果到下一个面片的流量为0，下一个面片不入队
				continue;
			if (Faces[Faces[CurFace].Neighbors[j].Face].Kind == FUZZY)//如果没有到达分解区域，就入队，并更改颜色
			{
				Faces[Faces[CurFace].Neighbors[j].Face].Kind = A;
				Faces[Faces[CurFace].Neighbors[j].Face].State = DISCOVERED;
				NeighborQ.push(Faces[CurFace].Neighbors[j].Face);
			}
		}
	}
	for (unsigned int i = 0; i < FuzzyFacesIdx.size(); i++)
	{
		if (Faces[FuzzyFacesIdx[i]].Kind == FUZZY)
		{
			Faces[FuzzyFacesIdx[i]].Kind = B;
			for (unsigned int j = 0; j < 3; j++)
			{
				Vertices[Faces[FuzzyFacesIdx[i]].V[j]].ColorData.R = 0;
				Vertices[Faces[FuzzyFacesIdx[i]].V[j]].ColorData.G = 255;
				Vertices[Faces[FuzzyFacesIdx[i]].V[j]].ColorData.B = 0;
			}
		}
	}
}

void clearParent(vector<Face>& Faces, vector<int>& FacesIdx)
{
	unsigned int FacesCnt = FacesIdx.size();
	for (unsigned int i = 0; i < FacesCnt; i++)
	{
		Faces[FacesIdx[i]].Parent = -1;
		Faces[FacesIdx[i]].CurMaxCap = 0;//初始的最大流为0
	}
}
