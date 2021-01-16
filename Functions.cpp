#include "Functions.h"

void splitString(const std::string & s, std::vector<std::string>& v, const std::string & c)
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

void dijkstra(vector<Vertice>& Vertices, vector<Face>& Faces, double**& ppEdges)
{
	unsigned int FacesCnt = Faces.size();
	priority_queue<Edge> EdgeQ;
	for (unsigned int i = 0; i < FacesCnt; i++)
	{
		Edge E, Temp;
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
			int EdgeCnt = Faces[E.End].Neighbors.size();
			for (int j = 0; j < EdgeCnt; j++)
			{
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
}

void clearState(vector<Face>& Faces)
{
	unsigned int FacesCnt = Faces.size();
	for (unsigned int i = 0; i < FacesCnt; i++)
	{
		Faces[i].State = UNDISCOVERED;
	}
}

void decideReps(vector<Face>& Faces, double **& ppEdges, vector<int>& Reps)
{
	unsigned int FacesCnt = Faces.size();
	double LastMinRepDist;//上一次聚类中心之间的最小距离
	double MaxRepDistChange;//聚类中心之间最小距离变化的最大值
	unsigned int RepsCnt = 0;//变化前聚类中心数量
	unsigned int MaxChangeStayCnt;//变化最大值保持次数的长短


	double MinDist = INF;//求解最初的聚类中心
	int InitialRep;
	for (unsigned int i = 0; i < FacesCnt; i++)
	{
		double DistSum = 0;
		for (unsigned int j = 0; j < FacesCnt; j++)
		{
			DistSum += ppEdges[i][j];
			if (DistSum > MinDist)
				break;
		}
		if (DistSum < MinDist)//如果和小于当前的最小值，更新初始聚类中心
		{
			MinDist = DistSum;
			InitialRep = i;
		}
	}
	Reps.push_back(InitialRep);

	while (true)
	{
		double MaxDist = 0;//最大的距离
		int NewRep;//新的聚类中心
		for (unsigned int i = 0; i < FacesCnt; i++)
		{
			double MinDist = INF;//到聚类中心的最小距离
			for (unsigned int j = 0; j < Reps.size(); j++)
			{
				if (ppEdges[i][Reps[j]] < MinDist)
					MinDist = ppEdges[i][Reps[j]];
			}
			if (MinDist > MaxDist)//如果当前面片到已有聚类中心的最小距离大于之前的最小距离，就用当前的面片代替之前的面片
			{
				MaxDist = MinDist;
				NewRep = i;
			}
		}

		Reps.push_back(NewRep);

		double CurRepsMinDist = INF;//当前不同聚类中心间的最小距离

		for (unsigned int i = 0; i < Reps.size() - 1; i++)
		{
			if (ppEdges[Reps[i]][NewRep] < CurRepsMinDist)//如果当前两个聚类中心间的距离更小，就更新聚类中心之间的最小距离
			{
				CurRepsMinDist = ppEdges[Reps[i]][NewRep];
			}
		}

		if (!RepsCnt)//如果之前只有一个聚类中心，就不计算差值
		{
			RepsCnt = 1;
			LastMinRepDist = CurRepsMinDist;
			MaxRepDistChange = 0;//变化量设为0
			MaxChangeStayCnt = 0;//维持次数设为0
		}
		else
		{
			if (LastMinRepDist - CurRepsMinDist > MaxRepDistChange)//如果变化更快，就更新聚类中心最大变化斜率
			{
				MaxRepDistChange = LastMinRepDist - CurRepsMinDist;
				RepsCnt = Reps.size() - 1;
				MaxChangeStayCnt = 0;
			}
			else//如果变化更慢，最大变化保持
			{
				MaxChangeStayCnt++;
			}
			LastMinRepDist = CurRepsMinDist;
		}

		if (MaxChangeStayCnt >= 15)
		{
			while (Reps.size() > RepsCnt)
			{
				Reps.pop_back();//将后续聚类中心清除
			}
			return;
		}
	}
}

void cluster(vector<Face>& Faces, vector<int>& Reps, double**& ppEdges)
{
	unsigned int FacesCnt = Faces.size();
	vector<int> LastReps = Reps;
	vector<pair<int, double> > SmallestReps;//新聚类种子
	cout << "初始聚类种子: ";
	for (unsigned int i = 0; i < Reps.size(); i++)
	{
		cout << Reps[i] + 1 << ' ';
		pair<int, double> Temp(Reps[i], 0);
		SmallestReps.push_back(Temp);
	}
	cout << endl;


	//初始化新聚类种子
	for (unsigned int i = 0; i < FacesCnt; i++)
	{
		double SumDist = 0;;//到各个聚类中心距离倒数和
		for (unsigned int j = 0; j < Reps.size(); j++)
		{
			if (ppEdges[i][Reps[j]] <= 1E-6)//避免出现除0的情况
				continue;
			SumDist += 1 / ppEdges[i][Reps[j]];
		}
		for (unsigned int j = 0; j < Reps.size(); j++)
		{
			if (ppEdges[i][Reps[j]] <= 1E-6)//避免出现除0的情况
			{
				Faces[i].Pl.push_back(1);
				continue;
			}
			double Temp = (1 / ppEdges[i][Reps[j]]) / SumDist;
			Faces[i].Pl.push_back(Temp);
		}
	}
	for (unsigned int j = 0; j < FacesCnt; j++)
	{
		for (unsigned int k = 0; k < Reps.size(); k++)
		{
			SmallestReps[k].second += Faces[j].Pl[k] * ppEdges[Reps[k]][j];
		}
	}
	cout << "初始加权距离和: ";
	for (unsigned int k = 0; k < Reps.size(); k++)
	{
		cout << SmallestReps[k].second << ' ';
	}
	cout << endl;

	int Cnt = 0;
	do
	{
		Cnt++;
		LastReps = Reps;

		for (unsigned int i = 0; i < FacesCnt; i++)
		{
			for (unsigned int k = 0; k < Reps.size(); k++)
			{
				double TempSum = 0;//临时的加权和
				for (unsigned int j = 0; j < FacesCnt; j++)
				{
					TempSum += Faces[j].Pl[k] * ppEdges[i][j];
					if (TempSum > SmallestReps[k].second)
						break;
				}
				//更新聚类种子
				if (TempSum < SmallestReps[k].second)
				{
					SmallestReps[k].first = i;
					SmallestReps[k].second = TempSum;
				}
			}
		}
		for (unsigned int j = 0; j < Reps.size(); j++)
			Reps[j] = SmallestReps[j].first;
		if (Reps == LastReps)//如果聚类种子没变就直接跳出循环
			break;

		cout << "第" << Cnt << "次迭代:" << endl <<
			"聚类种子: ";
		for (unsigned int i = 0; i < Reps.size(); i++)
		{
			cout << Reps[i] << ' ';
		}
		cout << endl;
		cout << "加权距离和: ";
		for (unsigned int k = 0; k < Reps.size(); k++)
		{
			cout << SmallestReps[k].second << ' ';
		}
		cout << endl;
		for (unsigned int i = 0; i < FacesCnt; i++)
		{
			double SumDist = 0;;//到各个聚类中心距离倒数和
			for (unsigned int j = 0; j < Reps.size(); j++)
			{
				if (ppEdges[i][Reps[j]] <= 1E-6)//避免出现除0的情况
					continue;
				SumDist += 1 / ppEdges[i][Reps[j]];
			}
			for (unsigned int j = 0; j < Reps.size(); j++)
			{
				if (ppEdges[i][Reps[j]] <= 1E-6)//避免出现除0的情况
				{
					Faces[i].Pl[j] = 1;
					continue;
				}
				double Temp = (1 / ppEdges[i][Reps[j]]) / SumDist;
				Faces[i].Pl[j] = Temp;
			}
		}
	} while (Cnt < 20);
}

void Div(vector<Vertice>& Vertices, vector<Face>& Faces, vector<int>& Reps, double AverageAngleDist, double**& ppEdges, string& FileName)
{
	map<pair<int, int>, vector<int> > FuzzyIdx;
	//模糊分割
	unsigned int FacesCnt = Faces.size();
	for (unsigned int i = 0; i < FacesCnt; i++)
	{
		vector<double> TempPl = Faces[i].Pl;
		sort(TempPl.rbegin(), TempPl.rend());
		double MaxP = TempPl[0];//最大的两个概率
		double MaxSecP = TempPl[1];
		int R[2];
		unsigned int k = 0;
		for (unsigned int j = 0; j < Faces[i].Pl.size(); j++)//找到最大的两个概率对应的聚类中心
		{
			if ((Faces[i].Pl[j] <= MaxP + 1E-6 && Faces[i].Pl[j] >= MaxP - 1E-6)
				|| (Faces[i].Pl[j] <= MaxSecP + 1E-6 && Faces[i].Pl[j] >= MaxSecP - 1E-6))
			{
				R[k] = j;
				if (k == 1)
					break;
				k++;
			}
		}

		double NormFactor = Faces[i].Pl[R[0]] + Faces[i].Pl[R[1]];
		Faces[i].Pl[R[0]] /= NormFactor;//归一化
		Faces[i].Pl[R[1]] /= NormFactor;

		if (Faces[i].Pl[R[0]] > 0.5 + EPSILON)//如果最大的概率比第二大的大一个阈值，就认为落在第一大的聚类中心内部，反之同理
		{
			Faces[i].Kind = DETERMINED;
			Faces[i].PossibleRep.first = Reps[R[0]];
			Faces[i].PossibleRep.second = Reps[R[1]];
		}
		else if (Faces[i].Pl[R[1]] > 0.5 + EPSILON)
		{
			Faces[i].Kind = DETERMINED;
			Faces[i].PossibleRep.first = Reps[R[1]];
			Faces[i].PossibleRep.second = Reps[R[0]];
		}
		else
		{
			Faces[i].Kind = FUZZY;
			Faces[i].PossibleRep.first = Reps[R[0]];
			Faces[i].PossibleRep.second = Reps[R[1]];

			FuzzyIdx[Faces[i].PossibleRep].push_back(i);
		}
	}

	//上色
	unordered_map<int, Color> NewColor;
	int ColorDiff;
	if ((Reps.size() + FuzzyIdx.size()) % 6 == 0)
		ColorDiff = 255 / (maxNum(Reps.size(), FuzzyIdx.size()) / 6);
	else
		ColorDiff = 255 / (maxNum(Reps.size(), FuzzyIdx.size()) / 6 + 1);
	for (unsigned int i = 0; i < Reps.size(); i++)
	{
		if (i % 6 == 0)
		{
			NewColor[Reps[i]] = Color((i / 6 + 1) * ColorDiff, 0, 0);
		}
		else if (i % 6 == 1)
		{
			NewColor[Reps[i]] = Color(0, (i / 6 + 1) * ColorDiff, 0);
		}
		else if(i % 6 == 2)
		{
			NewColor[Reps[i]] = Color(0, 0, (i / 6 + 1) * ColorDiff);
		}
		else if (i % 6 == 3)
		{
			NewColor[Reps[i]] = Color((i / 6 + 1) * ColorDiff, (i / 6 + 1) * ColorDiff, 0);
		}
		else if (i % 6 == 4)
		{
			NewColor[Reps[i]] = Color((i / 6 + 1) * ColorDiff, 0, (i / 6 + 1) * ColorDiff);
		}
		else
		{
			NewColor[Reps[i]] = Color(0, (i / 6 + 1) * ColorDiff, (i / 6 + 1) * ColorDiff);
		}
	}

	//模糊区域上色以及模糊区域的重新划分
	unsigned int FuzzyColorCnt = 0;
	map<pair<int, int>, vector<int> >::iterator It;
	for (It = FuzzyIdx.begin(); It != FuzzyIdx.end(); It++)
	{
		unsigned int FuzzyCnt = (*It).second.size();
		//上色
		Color Temp;
		if (FuzzyColorCnt % 6 == 0)
			Temp = Color((FuzzyColorCnt / 6 + 1) * ColorDiff, 0, 0);
		else if (FuzzyColorCnt % 3 == 1)
			Temp = Color(0, (FuzzyColorCnt / 6 + 1) * ColorDiff, 0);
		else if (FuzzyColorCnt % 3 == 2)
			Temp = Color(0, 0, (FuzzyColorCnt / 6 + 1) * ColorDiff);
		else if (FuzzyColorCnt % 3 == 3)
			Temp = Color((FuzzyColorCnt / 6 + 1) * ColorDiff, (FuzzyColorCnt / 6 + 1) * ColorDiff, 0);
		else if (FuzzyColorCnt % 3 == 4)
			Temp = Color((FuzzyColorCnt / 6 + 1) * ColorDiff, 0, (FuzzyColorCnt / 6 + 1) * ColorDiff);
		else
			Temp = Color(0, (FuzzyColorCnt / 6 + 1) * ColorDiff, (FuzzyColorCnt / 6 + 1) * ColorDiff);

		FuzzyColorCnt++;

		for (unsigned int i = 0; i < FuzzyCnt; i++)
		{
			for (unsigned int j = 0; j < 3; j++)
				Vertices[Faces[(*It).second[i]].V[j]].ColorData = Temp;

			//重新划分模糊区域
			for (unsigned int j = 0; j < Faces[(*It).second[i]].Neighbors.size(); j++)
			{
				if (Faces[Faces[(*It).second[i]].Neighbors[j].Face].Kind != FUZZY)
				{
					if (Faces[Faces[(*It).second[i]].Neighbors[j].Face].PossibleRep == Faces[(*It).second[i]].PossibleRep
						|| (Faces[Faces[(*It).second[i]].Neighbors[j].Face].PossibleRep.first == Faces[(*It).second[i]].PossibleRep.second
							&& Faces[Faces[(*It).second[i]].Neighbors[j].Face].PossibleRep.second == Faces[(*It).second[i]].PossibleRep.first))
					{
						Faces[Faces[(*It).second[i]].Neighbors[j].Face].Kind = FUZZY;
						(*It).second.push_back(Faces[(*It).second[i]].Neighbors[j].Face);//将边界上与模糊区域相交的面片也变成模糊面
					}
				}
			}
		}

	}

	outputObjFile("Mid" + FileName, Vertices, Faces);//不输出黑色

	for (unsigned int i = 0; i < FacesCnt; i++)
	{
		if (Faces[i].Kind == DETERMINED)
		{
			for (unsigned int j = 0; j < 3; j++)
				Vertices[Faces[i].V[j]].ColorData = NewColor[Faces[i].PossibleRep.first];
		}
	}

	for (It = FuzzyIdx.begin(); It != FuzzyIdx.end(); It++)
	{
		accurateDiv(Vertices, Faces, (*It).second, AverageAngleDist, ppEdges, (*It).first.first, (*It).first.second, NewColor);
	}
}

void accurateDiv(vector<Vertice>& Vertices, vector<Face>& Faces, vector<int>& FuzzyFacesIdx, double AverageAngleDist, double**& ppEdges, int Repa, int Repb, unordered_map<int, Color>& NewColor)
{
	vector<int> ANeighbors;
	pair<int, int> CurReps(Repa, Repb);
	//计算每条边的流量
	unsigned int FuzzyFacesCnt = FuzzyFacesIdx.size();
	for (unsigned int i = 0; i < FuzzyFacesCnt; i++)
	{
		for (unsigned int j = 0; j < Faces[FuzzyFacesIdx[i]].Neighbors.size(); j++)
		{
			if (Faces[Faces[FuzzyFacesIdx[i]].Neighbors[j].Face].Kind != FUZZY)//如果边上是A或B区域的面片，则流量为无穷
			{
				if (Faces[Faces[FuzzyFacesIdx[i]].Neighbors[j].Face].PossibleRep.first == Repa)//如果边上是A区域面片，就把面片放入A区域的邻域中
				{
					ppEdges[FuzzyFacesIdx[i]][Repa] = INF;
					ppEdges[Repa][FuzzyFacesIdx[i]] = INF;
					ANeighbors.push_back(FuzzyFacesIdx[i]);
				}
				else if (Faces[Faces[FuzzyFacesIdx[i]].Neighbors[j].Face].PossibleRep.second == Repb)//如果边上是B区域面片
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
				if (ppEdges[CurFace][Faces[CurFace].Neighbors[j].Face] <= 1E-6)//如果到下一个面片的流量为0，下一个面片不入队
					continue;
				if (Faces[Faces[CurFace].Neighbors[j].Face].Kind == DETERMINED && Faces[Faces[CurFace].Neighbors[j].Face].PossibleRep.first == Repb)//如果下一个面片是在B区域内，找到通路
				{
					if (Faces[CurFace].CurMaxCap > CapMax)
					{
						PathExist = true;
						CapMax = Faces[CurFace].CurMaxCap;//设置最大流为当前面片中存储的最大流
						EndFace = CurFace;//当前面片为最后结束的面片
					}
				}
				if (Faces[Faces[CurFace].Neighbors[j].Face].Kind == FUZZY &&
					(Faces[Faces[CurFace].Neighbors[j].Face].PossibleRep == CurReps ||
					(Faces[Faces[CurFace].Neighbors[j].Face].PossibleRep.first == CurReps.second &&
						Faces[Faces[CurFace].Neighbors[j].Face].PossibleRep.second == CurReps.first)))//如果没有到达B区域，就入队，并记录父节点位置
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
		Faces[ANeighbors[i]].Kind = DETERMINED;
		Faces[ANeighbors[i]].State = DISCOVERED;
		Faces[ANeighbors[i]].PossibleRep.first = Repa;
		Faces[ANeighbors[i]].PossibleRep.second = Repb;
	}
	while (!NeighborQ.empty())
	{
		int CurFace = NeighborQ.front();
		NeighborQ.pop();
		//上色
		for (unsigned int k = 0; k < 3; k++)
			Vertices[Faces[CurFace].V[k]].ColorData = NewColor[Repa];

		for (unsigned int j = 0; j < Faces[CurFace].Neighbors.size(); j++)
		{
			if (Faces[Faces[CurFace].Neighbors[j].Face].State == DISCOVERED)
				continue;
			if (Faces[Faces[CurFace].Neighbors[j].Face].Kind != FUZZY)//如果下一个面片是在不在模糊区域内，不入队
				continue;
			if (ppEdges[CurFace][Faces[CurFace].Neighbors[j].Face] <= 1E-6)//如果到下一个面片的流量为0，下一个面片不入队
				continue;
			if (Faces[Faces[CurFace].Neighbors[j].Face].Kind == FUZZY &&
				(Faces[Faces[CurFace].Neighbors[j].Face].PossibleRep == CurReps ||
				(Faces[Faces[CurFace].Neighbors[j].Face].PossibleRep.first == CurReps.second &&
					Faces[Faces[CurFace].Neighbors[j].Face].PossibleRep.second == CurReps.first)))//如果没有到达分解区域，就入队
			{
				Faces[Faces[CurFace].Neighbors[j].Face].Kind = DETERMINED;
				Faces[Faces[CurFace].Neighbors[j].Face].State = DISCOVERED;
				Faces[Faces[CurFace].Neighbors[j].Face].PossibleRep.first = Repa;
				Faces[Faces[CurFace].Neighbors[j].Face].PossibleRep.second = Repb;
				NeighborQ.push(Faces[CurFace].Neighbors[j].Face);
			}
		}
	}
	for (unsigned int i = 0; i < FuzzyFacesIdx.size(); i++)
	{
		if (Faces[FuzzyFacesIdx[i]].Kind == FUZZY)
		{
			Faces[FuzzyFacesIdx[i]].Kind = DETERMINED;
			Faces[FuzzyFacesIdx[i]].PossibleRep.first = Repb;
			Faces[FuzzyFacesIdx[i]].PossibleRep.second = Repa;
			//上色
			for (unsigned int k = 0; k < 3; k++)
				Vertices[Faces[FuzzyFacesIdx[i]].V[k]].ColorData = NewColor[Repb];
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
