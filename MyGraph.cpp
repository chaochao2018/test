#include "MyGraph.h"


double MyGraph::errorDeviation = eps*100;//不能太大
MyGraph *MyGraph::ginstance = NULL;

MyGraph::MyGraph()
{
	cout << "ERROR:please input a filename" << endl;
}
MyGraph::~MyGraph()
{
}
MyGraph::MyGraph(string filename) :MySvgFile(filename)
{
	//cout << "MyGraph(char *filename)" << endl;
	 file_width=0;
	 file_height=0;
	 setGraphInstance();
}

void MyGraph::constuctEdgeGraph(double w, double h,string file_v)
{
	//cout<<this->cubicBezier_vec_origin.size()<<endl;

	file_width = w;
	file_height = h;
	file_view = file_v;
	vector<PolyLine>::iterator iter_pl = polyline_vec.begin();
	double eps_dis = eps_isOnsegment * 2;
	int dele_count = 0;
	while (iter_pl != polyline_vec.end())//删除掉非常短（不必要的）的线，
	{
		//if (LineSegmentLength(iter_pl->pointlist.front(), iter_pl->pointlist.back()) < eps_dis)//eps_isOnsegment
		if (isPointRoughEqual(iter_pl->pointlist.front(), iter_pl->pointlist.back()))//
		{
			iter_pl = polyline_vec.erase(iter_pl);
			dele_count++;
		}
		else
			++iter_pl;
	}

	int indexcount = 0;
	myEdgeGraph.edgeNumber = 0;
	if (!cubicBezier_vec.empty())
	{
		for (unsigned int i = 0; i < cubicBezier_vec.size(); ++i)
		{
			EdgeStruct anedge;
			anedge.index = indexcount;
			anedge.curveType = 'C';//anedge.curveType = 'c';
			anedge.curveIndex = i;
			anedge.T1.p_position = cubicBezier_vec[i].p0;
			anedge.T2.p_position = cubicBezier_vec[i].p3;
			myEdgeGraph.edgeNumber++;
			myEdgeGraph.edgesOfGraph.push_back(anedge);
			indexcount++;
		}
	}
	if (!quadraticBezier_vec.empty())
	{
		for (unsigned int i = 0; i < quadraticBezier_vec.size(); ++i)
		{
			EdgeStruct anedge;
			anedge.index = indexcount;
			anedge.curveType = 'Q';//anedge.curveType = 'q';
			anedge.curveIndex = i;
			anedge.T1.p_position = quadraticBezier_vec[i].p0;
			anedge.T2.p_position = quadraticBezier_vec[i].p2;
			myEdgeGraph.edgeNumber++;
			myEdgeGraph.edgesOfGraph.push_back(anedge);
			indexcount++;
		}
	}
	if (!polyline_vec.empty())
	{
		//int tmp_num = indexcount;	
		//vector<PolyLine> tmp_poly;
		//vector<EdgeStruct> tmp_edge_vec;
		//tbb::parallel_for(tbb::blocked_range<size_t>(0, polyline_vec.size()), [&tmp_num, &polyline_vec](tbb::blocked_range<size_t>& r, vector<PolyLine> &tmp_poly, vector<EdgeStruct> &tmp_edge_vec){for (size_t ii = r.begin(); ii != r.end(); ++ii){
		//	//tbb::mutex mu;
		//	//mu.lock();
		//	EdgeStruct anedge;
		//	anedge.index = tmp_num;
		//	anedge.curveType = 'L';//'l'
		//	anedge.curveIndex = ii;
		//	anedge.T1.p_position = tmp_poly[ii].pointlist.front();
		//	anedge.T2.p_position = tmp_poly[ii].pointlist.back();
		//	tmp_edge_vec.push_back(anedge);
		//	tmp_num++;
		//	//mu.unlock();
		//} });

		for (unsigned int i = 0; i < polyline_vec.size(); ++i)
		{
			EdgeStruct anedge;
			anedge.index = indexcount;
			anedge.curveType = 'L';//'l'
			anedge.curveIndex = i;
			anedge.T1.p_position = polyline_vec[i].pointlist.front();
			anedge.T2.p_position = polyline_vec[i].pointlist.back();
			myEdgeGraph.edgeNumber++;
			myEdgeGraph.edgesOfGraph.push_back(anedge);
			indexcount++;
		}
	}
	for (int i = 0; i < myEdgeGraph.edgeNumber; ++i)
	{
		if (myEdgeGraph.edgesOfGraph.at(i).index >= 0)//index=-1的话说明不可用，（太短，标记已被删除）
		{
			vector<int> combineFlag;// 记录下短线这种特殊情况，记录j的值
			for (int j = 0; j < myEdgeGraph.edgeNumber; ++j)
			{
				if (i != j && myEdgeGraph.edgesOfGraph.at(j).index >= 0)
				{
					string strerr;
					if (isPointRoughEqual(myEdgeGraph.edgesOfGraph.at(i).T1.p_position, myEdgeGraph.edgesOfGraph.at(j).T1.p_position))
					{
						strerr.push_back('a');
						myEdgeGraph.edgesOfGraph.at(i).T1.adjEdgeIndex_vec.push_back(myEdgeGraph.edgesOfGraph.at(j).index);//myEdgeGraph.edgesOfGraph.at(j)!=j!!!
					}
					else if (isPointRoughEqual(myEdgeGraph.edgesOfGraph.at(i).T1.p_position, myEdgeGraph.edgesOfGraph.at(j).T2.p_position))
					{
						strerr.push_back('b');
						myEdgeGraph.edgesOfGraph.at(i).T1.adjEdgeIndex_vec.push_back(myEdgeGraph.edgesOfGraph.at(j).index);
					}
					if (isPointRoughEqual(myEdgeGraph.edgesOfGraph.at(i).T2.p_position, myEdgeGraph.edgesOfGraph.at(j).T1.p_position))
					{
						strerr.push_back('c');
						myEdgeGraph.edgesOfGraph.at(i).T2.adjEdgeIndex_vec.push_back(myEdgeGraph.edgesOfGraph.at(j).index);
					}
					else if (isPointRoughEqual(myEdgeGraph.edgesOfGraph.at(i).T2.p_position, myEdgeGraph.edgesOfGraph.at(j).T2.p_position))
					{
						strerr.push_back('d');
						myEdgeGraph.edgesOfGraph.at(i).T2.adjEdgeIndex_vec.push_back(myEdgeGraph.edgesOfGraph.at(j).index);
					}
					if (strerr.find("ac") != string::npos)//说明第i条线很短，其T1，T2两个端点都等于 第j条线的T1端
					{
						combineFlag.push_back(j);
					}
					if (strerr.find("bd") != string::npos)//同理，第i条线很短，其T1，T2两个端点都等于 第j条线的T2端
					{
						combineFlag.push_back(j);
					}
				}
			}
			if (!combineFlag.empty())
			{
				TerminalVertex_E t1, t2, t3, new_t4;
				t1 = myEdgeGraph.edgesOfGraph.at(i).T1;
				t2 = myEdgeGraph.edgesOfGraph.at(i).T2;
				if (isPointRoughEqual(myEdgeGraph.edgesOfGraph.at(combineFlag[0]).T1.p_position, myEdgeGraph.edgesOfGraph.at(i).T1.p_position))
				{
					t3 = myEdgeGraph.edgesOfGraph.at(combineFlag[0]).T1;
				}
				else if (isPointRoughEqual(myEdgeGraph.edgesOfGraph.at(combineFlag[0]).T2.p_position, myEdgeGraph.edgesOfGraph.at(i).T1.p_position))
				{
					t3 = myEdgeGraph.edgesOfGraph.at(combineFlag[0]).T2;
				}
				mypoint2f cmposition = (t1.p_position + t2.p_position + t3.p_position);
				new_t4.p_position.x = cmposition.x / 3;
				new_t4.p_position.y = cmposition.y / 3;
				vector<int> allAdjEdges;  //记录第i条线两端相邻的所有edge
				allAdjEdges.insert(allAdjEdges.end(), myEdgeGraph.edgesOfGraph.at(i).T1.adjEdgeIndex_vec.begin(), myEdgeGraph.edgesOfGraph.at(i).T1.adjEdgeIndex_vec.end());
				allAdjEdges.insert(allAdjEdges.end(), myEdgeGraph.edgesOfGraph.at(i).T2.adjEdgeIndex_vec.begin(), myEdgeGraph.edgesOfGraph.at(i).T2.adjEdgeIndex_vec.end());
				deleRepeatPoint(&allAdjEdges);
				//遍历allAdjEdges，查看其中有没有很短的edge并删除，存入shortedge
				vector<int> shortedge;
				vector<int>::iterator it_adjedge = allAdjEdges.begin();
				while (it_adjedge != allAdjEdges.end())
				{
					int cmEdge = *it_adjedge;
					TerminalVertex_E tv1 = myEdgeGraph.edgesOfGraph.at(cmEdge).T1;
					TerminalVertex_E tv2 = myEdgeGraph.edgesOfGraph.at(cmEdge).T2;//cmedge的
					if (isPointRoughEqual(myEdgeGraph.edgesOfGraph.at(i).T1.p_position, tv1.p_position) && isPointRoughEqual(myEdgeGraph.edgesOfGraph.at(i).T1.p_position, tv2.p_position))
					{
						myEdgeGraph.edgesOfGraph.at(cmEdge).index = -1;
						shortedge.push_back(cmEdge);
						it_adjedge = allAdjEdges.erase(it_adjedge);
					}
					else if (isPointRoughEqual(myEdgeGraph.edgesOfGraph.at(i).T2.p_position, tv1.p_position) && isPointRoughEqual(myEdgeGraph.edgesOfGraph.at(i).T2.p_position, tv2.p_position))
					{
						myEdgeGraph.edgesOfGraph.at(cmEdge).index = -1;
						shortedge.push_back(cmEdge);
						it_adjedge = allAdjEdges.erase(it_adjedge);
					}
					else
						it_adjedge++;
				}
				for (unsigned int k = 0; k < allAdjEdges.size(); ++k)//遍历 线段i两端所连接的每条线。(1)改变端点坐标(2)删除其中T1/T2端点中的adjEdgeIndex中的i、shortedge(3)修改T1/T2的adjEdgeIndex的值（增添）(4) 修改对应curve的pointlist
				{
					int cmEdge = allAdjEdges[k];
					TerminalVertex_E tv1 = myEdgeGraph.edgesOfGraph.at(cmEdge).T1;
					TerminalVertex_E tv2 = myEdgeGraph.edgesOfGraph.at(cmEdge).T2;
					if (isPointRoughEqual(myEdgeGraph.edgesOfGraph.at(i).T1.p_position, tv1.p_position))//cmEdge的T1端在中心
					{
						myEdgeGraph.edgesOfGraph.at(cmEdge).T1.p_position = new_t4.p_position;//1、修改端点的坐标
						if (!myEdgeGraph.edgesOfGraph.at(cmEdge).T1.adjEdgeIndex_vec.empty())
						{
							vector<int>::iterator iter_adjedge = myEdgeGraph.edgesOfGraph.at(cmEdge).T1.adjEdgeIndex_vec.begin();
							while (iter_adjedge != myEdgeGraph.edgesOfGraph.at(cmEdge).T1.adjEdgeIndex_vec.end())//2、删除端点中adjEdgeIndex_ver含有i的值
							{
								if (*iter_adjedge == i)
									iter_adjedge = myEdgeGraph.edgesOfGraph.at(cmEdge).T1.adjEdgeIndex_vec.erase(iter_adjedge);
								else
									iter_adjedge++;
							}
							for (unsigned int si = 0; si < shortedge.size(); ++si)
							{
								iter_adjedge = myEdgeGraph.edgesOfGraph.at(cmEdge).T1.adjEdgeIndex_vec.begin();
								while (iter_adjedge != myEdgeGraph.edgesOfGraph.at(cmEdge).T1.adjEdgeIndex_vec.end())//2、删除端点中adjEdgeIndex_ver含有shortedge的值
								{
									if (*iter_adjedge == shortedge[si])
										iter_adjedge = myEdgeGraph.edgesOfGraph.at(cmEdge).T1.adjEdgeIndex_vec.erase(iter_adjedge);
									else
										iter_adjedge++;
								}
							}
							//3、（增添or更新）端点中adjEdgeIndex_ver的值
							myEdgeGraph.edgesOfGraph.at(cmEdge).T1.adjEdgeIndex_vec.clear();
							myEdgeGraph.edgesOfGraph.at(cmEdge).T1.adjEdgeIndex_vec.swap(vector<int>());
							for (unsigned int k_i = 0; k_i < allAdjEdges.size(); ++k_i)
							{
								if (allAdjEdges[k_i] != cmEdge)
									myEdgeGraph.edgesOfGraph.at(cmEdge).T1.adjEdgeIndex_vec.push_back(allAdjEdges[k_i]);
							}
						}						
						//4、修改curve的长度
						int curindex = myEdgeGraph.edgesOfGraph.at(cmEdge).curveIndex;
						if (myEdgeGraph.edgesOfGraph.at(cmEdge).curveType == 'L')// 'l' 若是polyline直接修改polyline_vec中对应的线条
						{
							if (isPointRoughEqual(polyline_vec[curindex].pointlist.front(), tv1.p_position))
								polyline_vec[curindex].pointlist.front()=new_t4.p_position;
							else
								polyline_vec[curindex].pointlist.back()=new_t4.p_position;
						}
						else if (myEdgeGraph.edgesOfGraph.at(cmEdge).curveType == 'C')// 'c'   
						{
							//原方法一。若是贝塞尔曲线，则离散化成折线，
							PolyLine pl;
							getPListfromCubic(0, cubicBezier_vec[curindex], 0,&pl.pointlist);
							if (isPointRoughEqual(pl.pointlist.front(), tv1.p_position))
								pl.pointlist.insert(pl.pointlist.begin(), new_t4.p_position);
							else
								pl.pointlist.push_back(new_t4.p_position);
							pl.origi_type = 'C'; pl.origi_index = myEdgeGraph.edgesOfGraph.at(cmEdge).curveIndex;
							polyline_vec.push_back(pl);
							myEdgeGraph.edgesOfGraph.at(cmEdge).curveType = 'L';//'l'
							myEdgeGraph.edgesOfGraph.at(cmEdge).curveIndex = polyline_vec.size() - 1;

						}
						else if (myEdgeGraph.edgesOfGraph.at(cmEdge).curveType == 'Q')//q
						{
							PolyLine pl;
							getPListfromQuadratic(0, quadraticBezier_vec[curindex], 0,&pl.pointlist);
							if (isPointRoughEqual(pl.pointlist.front(), tv1.p_position))
								pl.pointlist.insert(pl.pointlist.begin(), new_t4.p_position);
							else
								pl.pointlist.push_back(new_t4.p_position);
							pl.origi_type = 'Q'; pl.origi_index = myEdgeGraph.edgesOfGraph.at(cmEdge).curveIndex;
							polyline_vec.push_back(pl);
							myEdgeGraph.edgesOfGraph.at(cmEdge).curveType = 'L';//'l'
							myEdgeGraph.edgesOfGraph.at(cmEdge).curveIndex = polyline_vec.size() - 1;
						}
					}
					else if (isPointRoughEqual(myEdgeGraph.edgesOfGraph.at(i).T1.p_position, tv2.p_position))//T2端在中心,同上
					{
						myEdgeGraph.edgesOfGraph.at(cmEdge).T2.p_position = new_t4.p_position;//修改端点的坐标
						if (!myEdgeGraph.edgesOfGraph.at(cmEdge).T2.adjEdgeIndex_vec.empty())
						{
							vector<int>::iterator iter_adjedge = myEdgeGraph.edgesOfGraph.at(cmEdge).T2.adjEdgeIndex_vec.begin();
							while (iter_adjedge != myEdgeGraph.edgesOfGraph.at(cmEdge).T2.adjEdgeIndex_vec.end())//2、删除端点中adjEdgeIndex_ver含有i的值
							{
								if (*iter_adjedge == i)
									iter_adjedge = myEdgeGraph.edgesOfGraph.at(cmEdge).T2.adjEdgeIndex_vec.erase(iter_adjedge);
								else
									iter_adjedge++;
							}
							for (unsigned int si = 0; si < shortedge.size(); ++si)
							{
								iter_adjedge = myEdgeGraph.edgesOfGraph.at(cmEdge).T2.adjEdgeIndex_vec.begin();
								while (iter_adjedge != myEdgeGraph.edgesOfGraph.at(cmEdge).T2.adjEdgeIndex_vec.end())//2、删除端点中adjEdgeIndex_ver含有shortedge的值
								{
									if (*iter_adjedge == shortedge[si])
										iter_adjedge = myEdgeGraph.edgesOfGraph.at(cmEdge).T2.adjEdgeIndex_vec.erase(iter_adjedge);
									else
										iter_adjedge++;
								}
							}
							//3、（增添or更新）端点中adjEdgeIndex_ver的值
							myEdgeGraph.edgesOfGraph.at(cmEdge).T2.adjEdgeIndex_vec.clear();
							myEdgeGraph.edgesOfGraph.at(cmEdge).T2.adjEdgeIndex_vec.swap(vector<int>());
							for (unsigned int k_i = 0; k_i < allAdjEdges.size(); ++k_i)
							{
								if (allAdjEdges[k_i] != cmEdge)
									myEdgeGraph.edgesOfGraph.at(cmEdge).T2.adjEdgeIndex_vec.push_back(allAdjEdges[k_i]);
							}
						}					
						//4、修改线的长度
						int curindex = myEdgeGraph.edgesOfGraph.at(cmEdge).curveIndex;
						if (myEdgeGraph.edgesOfGraph.at(cmEdge).curveType == 'L')//'l'若是polyline直接修改polyline_vec中对应的线条
						{
							if (isPointRoughEqual(polyline_vec[curindex].pointlist.front(), tv2.p_position))
								polyline_vec[curindex].pointlist.front()= new_t4.p_position;
							else
								polyline_vec[curindex].pointlist.back()=new_t4.p_position;
						}
						else if (myEdgeGraph.edgesOfGraph.at(cmEdge).curveType == 'C')//若是贝塞尔曲线，则离散化成折线，
						{
							PolyLine pl;
							getPListfromCubic(0, cubicBezier_vec[curindex],0, &pl.pointlist);
							if (isPointRoughEqual(pl.pointlist.front(), tv2.p_position))
								pl.pointlist.insert(pl.pointlist.begin(), new_t4.p_position);
							else
								pl.pointlist.push_back(new_t4.p_position);
							pl.origi_type = 'C'; pl.origi_index = myEdgeGraph.edgesOfGraph.at(cmEdge).curveIndex;
							polyline_vec.push_back(pl);
							myEdgeGraph.edgesOfGraph.at(cmEdge).curveType = 'L';//'l'
							myEdgeGraph.edgesOfGraph.at(cmEdge).curveIndex = polyline_vec.size() - 1;
						}
						else if (myEdgeGraph.edgesOfGraph.at(cmEdge).curveType == 'Q')
						{
							PolyLine pl;
							getPListfromQuadratic(0, quadraticBezier_vec[curindex], 0,&pl.pointlist);
							if (isPointRoughEqual(pl.pointlist.front(), tv2.p_position))
								pl.pointlist.insert(pl.pointlist.begin(), new_t4.p_position);
							else
								pl.pointlist.push_back(new_t4.p_position);
							pl.origi_type = 'Q'; pl.origi_index = myEdgeGraph.edgesOfGraph.at(cmEdge).curveIndex;
							polyline_vec.push_back(pl);
							myEdgeGraph.edgesOfGraph.at(cmEdge).curveType = 'L';//'l'
							myEdgeGraph.edgesOfGraph.at(cmEdge).curveIndex = polyline_vec.size() - 1;
						}
					}
				}
				myEdgeGraph.edgesOfGraph.at(i).T1.adjEdgeIndex_vec.clear();
				myEdgeGraph.edgesOfGraph.at(i).T2.adjEdgeIndex_vec.clear();
				myEdgeGraph.edgesOfGraph.at(i).index = -1;//-1表示不可用
			}
		}
	}
	cout << "cubic curve size=" << cubicBezier_vec.size() << ",quadratic curve size=" << quadraticBezier_vec.size() << ",polyline size=" << polyline_vec.size() << endl;
	cout << "Edge Number=" << myEdgeGraph.edgeNumber << endl;
}
void MyGraph::constructPointGraph()
{
	int *visited = new int[myEdgeGraph.edgesOfGraph.size()];
	for (unsigned int i = 0; i < myEdgeGraph.edgesOfGraph.size(); ++i)
	{
		visited[i] = 0;//0:unvisited,  1:visit T1,  2:visit T2,  3:visit both T1&T2
	}
	for (unsigned int i = 0; i < myEdgeGraph.edgesOfGraph.size(); ++i)
	{
		if (myEdgeGraph.edgesOfGraph.at(i).index >= 0)
		{
			if (visited[i] == 0)
			{
				visited[i] = 3;
				mypoint2f cp1 = myEdgeGraph.edgesOfGraph[i].T1.p_position;
				vector<TerminalPoint> terpoints_cp1;
				for (unsigned int j = 0; j < myEdgeGraph.edgesOfGraph[i].T1.adjEdgeIndex_vec.size(); ++j)
				{
					int adjeg = myEdgeGraph.edgesOfGraph[i].T1.adjEdgeIndex_vec[j];
					TerminalPoint tp;
					tp.endpoint = getOtherTPofEdge(myEdgeGraph.edgesOfGraph[adjeg], cp1)->p_position;
					tp.curveindex = myEdgeGraph.edgesOfGraph[adjeg].curveIndex;
					tp.curvetype = myEdgeGraph.edgesOfGraph[adjeg].curveType;
					terpoints_cp1.push_back(tp);
					if (visited[adjeg] == 0)
					{
						if (isPointRoughEqual(myEdgeGraph.edgesOfGraph[adjeg].T1.p_position, cp1))
						{
							visited[adjeg] = 1;
						}
						else
						{
							visited[adjeg] = 2;
						}
					}
					else
						visited[adjeg] = 3;
				}
				TerminalPoint tp_origi;
				tp_origi.endpoint = myEdgeGraph.edgesOfGraph[i].T2.p_position;
				tp_origi.curveindex = myEdgeGraph.edgesOfGraph[i].curveIndex;
				tp_origi.curvetype = myEdgeGraph.edgesOfGraph[i].curveType;
				terpoints_cp1.push_back(tp_origi);
				pointTopl.push_back(make_pair(cp1, terpoints_cp1));

				mypoint2f cp2 = myEdgeGraph.edgesOfGraph[i].T2.p_position;
				vector<TerminalPoint> terpoints_cp2;
				for (unsigned int j = 0; j < myEdgeGraph.edgesOfGraph[i].T2.adjEdgeIndex_vec.size(); ++j)
				{
					int adjeg = myEdgeGraph.edgesOfGraph[i].T2.adjEdgeIndex_vec[j];
					TerminalPoint tp;
					tp.endpoint = getOtherTPofEdge(myEdgeGraph.edgesOfGraph[adjeg], cp2)->p_position;
					tp.curveindex = myEdgeGraph.edgesOfGraph[adjeg].curveIndex;
					tp.curvetype = myEdgeGraph.edgesOfGraph[adjeg].curveType;
					terpoints_cp2.push_back(tp);
					if (visited[adjeg] == 0)
					{
						if (isPointRoughEqual(myEdgeGraph.edgesOfGraph[adjeg].T1.p_position, cp2))
						{
							visited[adjeg] = 1;
						}
						else
						{
							visited[adjeg] = 2;
						}
					}
					else
						visited[adjeg] = 3;
				}
				tp_origi.endpoint = myEdgeGraph.edgesOfGraph[i].T1.p_position;
				terpoints_cp2.push_back(tp_origi);
				pointTopl.push_back(make_pair(cp2, terpoints_cp2));
			}
			else if (visited[i] == 1)
			{
				visited[i] = 3;
				mypoint2f cp2 = myEdgeGraph.edgesOfGraph[i].T2.p_position;
				vector<TerminalPoint> terpoints_cp2;
				for (unsigned int j = 0; j < myEdgeGraph.edgesOfGraph[i].T2.adjEdgeIndex_vec.size(); ++j)
				{
					int adjeg = myEdgeGraph.edgesOfGraph[i].T2.adjEdgeIndex_vec[j];
					TerminalPoint tp;
					tp.endpoint = getOtherTPofEdge(myEdgeGraph.edgesOfGraph[adjeg], cp2)->p_position;
					tp.curveindex = myEdgeGraph.edgesOfGraph[adjeg].curveIndex;
					tp.curvetype = myEdgeGraph.edgesOfGraph[adjeg].curveType;
					terpoints_cp2.push_back(tp);
					if (visited[adjeg] == 0)
					{
						if (isPointRoughEqual(myEdgeGraph.edgesOfGraph[adjeg].T1.p_position, cp2))
						{
							visited[adjeg] = 1;
						}
						else
						{
							visited[adjeg] = 2;
						}
					}
					else
						visited[adjeg] = 3;
				}
				TerminalPoint tp_origi;
				tp_origi.endpoint = myEdgeGraph.edgesOfGraph[i].T1.p_position;
				tp_origi.curveindex = myEdgeGraph.edgesOfGraph[i].curveIndex;
				tp_origi.curvetype = myEdgeGraph.edgesOfGraph[i].curveType;
				terpoints_cp2.push_back(tp_origi);
				pointTopl.push_back(make_pair(cp2, terpoints_cp2));
			}
			else if (visited[i] == 2)
			{
				visited[i] = 3;
				mypoint2f cp1 = myEdgeGraph.edgesOfGraph[i].T1.p_position;
				vector<TerminalPoint> terpoints_cp1;
				for (unsigned int j = 0; j < myEdgeGraph.edgesOfGraph[i].T1.adjEdgeIndex_vec.size(); ++j)
				{
					int adjeg = myEdgeGraph.edgesOfGraph[i].T1.adjEdgeIndex_vec[j];
					TerminalPoint tp;
					tp.endpoint = getOtherTPofEdge(myEdgeGraph.edgesOfGraph[adjeg], cp1)->p_position;
					tp.curveindex = myEdgeGraph.edgesOfGraph[adjeg].curveIndex;
					tp.curvetype = myEdgeGraph.edgesOfGraph[adjeg].curveType;
					terpoints_cp1.push_back(tp);
					if (visited[adjeg] == 0)
					{
						if (isPointRoughEqual(myEdgeGraph.edgesOfGraph[adjeg].T1.p_position, cp1))
						{
							visited[adjeg] = 1;
						}
						else
						{
							visited[adjeg] = 2;
						}
					}
					else
						visited[adjeg] = 3;
				}
				TerminalPoint tp_origi;
				tp_origi.endpoint = myEdgeGraph.edgesOfGraph[i].T2.p_position;
				tp_origi.curveindex = myEdgeGraph.edgesOfGraph[i].curveIndex;
				tp_origi.curvetype = myEdgeGraph.edgesOfGraph[i].curveType;
				terpoints_cp1.push_back(tp_origi);
				pointTopl.push_back(make_pair(cp1, terpoints_cp1));
			}
		}
	}
	delete[] visited;
	visited = NULL;
	// combine two lines (delete point which degree is 2) & check and add aided midpoint 
	//for (unsigned int i = 0; i < pointTopl.size(); ++i)
	//{
	//	//vector<TerminalPoint> tpoints_of_i = pointTopl[i].second;
	//	if (pointTopl[i].second.size() == 2 )
	//	{
	//		if (!isPointRoughEqual(pointTopl[i].second[0].endpoint, pointTopl[i].second[1].endpoint))
	//		{
	//			combineTwoCurves(i, 0, 1);
	//		}
	//	}		
	//}
	for (unsigned int i = 0; i < pointTopl.size(); ++i)
	{
		if (pointTopl[i].second.size() > 1)
		{
			for (unsigned int j = 0; j < pointTopl[i].second.size(); ++j)
			{
				for (unsigned int k = j + 1; k < pointTopl[i].second.size(); ++k)
				{
					if (isPointRoughEqual(pointTopl[i].second[j].endpoint, pointTopl[i].second[k].endpoint))
					{
						addAssistantPoint(i, k, j);
					}
				}
			}
		}
	}
}
void MyGraph::combineTwoCurves(int i_ptopl, int tp_k, int tp_j)
{
	mypoint2f startp = pointTopl[i_ptopl].first;
	TerminalPoint tp1 = pointTopl[i_ptopl].second[tp_k];
	TerminalPoint tp2 = pointTopl[i_ptopl].second[tp_j];
	vector<mypoint2f> plist_vec1;
	vector<mypoint2f> plist_vec2;
	int stopflag = 0;
	int tp1_i = 0;
	int tp2_i = 0;
	for (unsigned int i = 0; i < pointTopl.size(); ++i)
	{
		tp1.curveindex;
		if (isPointRoughEqual(tp1.endpoint, pointTopl[i].first))
		{
			getPListOfaCurve(tp1.endpoint, startp, i, -1, &plist_vec1,0);
			if (!plist_vec1.empty())
			{
				tp1_i = i;
				stopflag++;
			}
			else
				cout << "error in combineTwoCurves";
		}
		else if (isPointRoughEqual(tp2.endpoint, pointTopl[i].first))
		{
			getPListOfaCurve(startp, tp2.endpoint, i_ptopl, -1, &plist_vec2,0);
			if (!plist_vec2.empty())
			{
				tp2_i = i;
				stopflag++;
			}
			else
				cout << "error in combineTwoCurves";
		}
		if (stopflag == 2)
			break;
	}
	//判断是不是能合并
	stopflag = 1;
	if (pointTopl[tp1_i].second.size()>2 && pointTopl[tp2_i].second.size() > 2)
	{
		for (unsigned int j = 0; j < pointTopl[tp1_i].second.size(); ++j)
		{
			if (isPointRoughEqual(pointTopl[tp1_i].second[j].endpoint, pointTopl[tp2_i].first))
			{
				stopflag = 0;//不能合并
				break;
			}
		}
	}
	if (stopflag == 1)
	{
		plist_vec1.pop_back();
		PolyLine apl;
		apl.pointlist.insert(apl.pointlist.end(), plist_vec1.begin(), plist_vec1.end());
		apl.pointlist.insert(apl.pointlist.end(), plist_vec2.begin(), plist_vec2.end());
		apl.origi_index = polyline_vec.size();
		apl.origi_type = 'L';
		polyline_vec.push_back(apl);
		for (unsigned int j = 0; j < pointTopl[tp1_i].second.size(); ++j)
		{
			if (isPointRoughEqual(pointTopl[tp1_i].second[j].endpoint, startp))
			{
				pointTopl[tp1_i].second[j].endpoint = tp2.endpoint;
				pointTopl[tp1_i].second[j].curvetype = 'L'; //'l' ?
				pointTopl[tp1_i].second[j].curveindex = polyline_vec.size() - 1;
				break;
			}
		}
		for (unsigned int j = 0; j < pointTopl[tp2_i].second.size(); ++j)
		{
			if (isPointRoughEqual(pointTopl[tp2_i].second[j].endpoint, startp))
			{
				pointTopl[tp2_i].second[j].endpoint = tp1.endpoint;
				pointTopl[tp2_i].second[j].curvetype = 'L';//'l' ?
				pointTopl[tp2_i].second[j].curveindex = polyline_vec.size() - 1;
				break;
			}
		}
		pointTopl[i_ptopl].second.clear();
		pointTopl[i_ptopl].second.swap(vector<TerminalPoint>());
	}
}
void MyGraph::addAssistantPoint(int i_ptopl, int tp_k, int tp_j)//vector<PolyLine> &polylineVec, vector<pair<mypoint2f, vector<TerminalPoint>>> &ptp_vec
{
	TerminalPoint tp;
	int k;
	if (pointTopl[i_ptopl].second.at(tp_k).curvetype != 'L')//'l'
	{
		tp = pointTopl[i_ptopl].second.at(tp_k);
		k = tp_k;
	}	
	else if (pointTopl[i_ptopl].second.at(tp_j).curvetype != 'L')//'l'
	{
		tp = pointTopl[i_ptopl].second.at(tp_j);
		k = tp_j;
	}
	else
	{
		int plist1 = polyline_vec[pointTopl[i_ptopl].second.at(tp_k).curveindex].pointlist.size();
		int plist2 = polyline_vec[pointTopl[i_ptopl].second.at(tp_j).curveindex].pointlist.size();
		if (plist2 > plist1)
		{
			tp = pointTopl[i_ptopl].second.at(tp_j);
			k = tp_j;
		}
		else
		{
			tp = pointTopl[i_ptopl].second.at(tp_k);
			k = tp_k;
		}
	}	
	if (tp.curvetype == 'C')//'c'
	{
		cubicBezier cbcurve = cubicBezier_vec[tp.curveindex];
		double p1[2] = { cbcurve.p0.x, cbcurve.p0.y };
		double p2[2] = { cbcurve.p1.x, cbcurve.p1.y };
		double p3[2] = { cbcurve.p2.x, cbcurve.p2.y };
		double p4[2] = { cbcurve.p3.x, cbcurve.p3.y };
		vector<double> split_t;
		split_t.push_back(0.5);
		split_t.push_back(1.0);
		vector<mypoint2f> seg_curve1;
		double a1, a2, a3, a4;
		double start_t = 0;
		mypoint2f apoint, midpoint;
		//三次公式 B(t)=p0*(1-t)^3 + p1*3t(1-t)^2 + p2*3t^2(1-t) + p3*t^3
		for (unsigned int j = 0; j < split_t.size(); ++j)
		{
			for (double t = start_t; t < split_t[j]; t += 0.02)
			{
				a1 = pow((1 - t), 3);
				a2 = pow((1 - t), 2) * 3 * t;
				a3 = 3 * t*t*(1 - t);
				a4 = t*t*t;
				apoint.x = a1*p1[0] + a2*p2[0] + a3*p3[0] + a4*p4[0];
				apoint.y = a1*p1[1] + a2*p2[1] + a3*p3[1] + a4*p4[1];
				seg_curve1.push_back(apoint);
			}
			a1 = pow((1 - split_t[j]), 3);
			a2 = pow((1 - split_t[j]), 2) * 3 * split_t[j];
			a3 = 3 * pow(split_t[j], 2)*(1 - split_t[j]);
			a4 = pow(split_t[j], 3);
			apoint.x = a1*p1[0] + a2*p2[0] + a3*p3[0] + a4*p4[0];
			apoint.y = a1*p1[1] + a2*p2[1] + a3*p3[1] + a4*p4[1];
			seg_curve1.push_back(apoint);

			PolyLine apl;
			apl.pointlist.insert(apl.pointlist.end(), seg_curve1.begin(), seg_curve1.end());
			apl.origi_type = 'C'; apl.origi_index = tp.curveindex;
			polyline_vec.push_back(apl);//polylineVec.push_back(seg_curve1);
			seg_curve1.swap(vector<mypoint2f>());
			start_t = split_t[j];
		}
		midpoint = polyline_vec.at(polyline_vec.size()-2).pointlist.back();
		//mycircle acir;///////////////////////////添加点-》circle_vec
		//acir.cx = midpoint.x;
		//acir.cy = midpoint.y;
		//acir.radius = 1;
		//circle_vec.push_back(acir);//pointTopl.push_back()
		vector<TerminalPoint> terpoints;
		TerminalPoint addtp;
		addtp.endpoint = cbcurve.p0;
		addtp.curveindex = polyline_vec.size() - 2;//polylineVec.size() - 2;
		addtp.curvetype = 'L';// 'l';
		terpoints.push_back(addtp);
		addtp.endpoint = cbcurve.p3;
		addtp.curveindex = polyline_vec.size() - 1;//polylineVec.size() - 1;
		addtp.curvetype = 'L';// 'l';
		terpoints.push_back(addtp);
		pointTopl.push_back(make_pair(midpoint, terpoints));//ptp_vec.push_back(make_pair(midpoint, terpoints));添加新的辅助点
		//modify pointTopl[i].second[k]... 修改被添加辅助点的那条线
		mypoint2f origistartp = pointTopl[i_ptopl].first;//ptp_vec[i].first;
		mypoint2f origiendp = pointTopl[i_ptopl].second[k].endpoint;//ptp_vec[i].second[k].endpoint;
		int origiIndex = pointTopl[i_ptopl].second[k].curveindex;//ptp_vec[i].second[k].curveindex;
		pointTopl[i_ptopl].second[k].curvetype = 'L';// 'l';  //ptp_vec[i].second[k].curvetype = 'l';
		if (isPointRoughEqual(origistartp, cbcurve.p0))//if (origistartp == cbcurve.p0)
			pointTopl[i_ptopl].second[k].curveindex = polyline_vec.size() - 2;//ptp_vec[i].second[k].curveindex = polylineVec.size() - 2;
		else
			pointTopl[i_ptopl].second[k].curveindex = polyline_vec.size() - 1;//ptp_vec[i].second[k].curveindex = polylineVec.size() - 1;
		pointTopl[i_ptopl].second[k].endpoint = midpoint;//ptp_vec[i].second[k].endpoint = midpoint;
		for (unsigned int v = 0; v < pointTopl.size(); ++v)
		{
			if (isPointRoughEqual(pointTopl[v].first, origiendp))//if (pointTopl[v].first == origiendp)// 以pointTopl[i_ptopl].second[k].endpoint为端点的线也要修改
			{
				for (unsigned int u = 0; u < pointTopl[v].second.size(); ++u)
				{
					if (isPointRoughEqual(pointTopl[v].second.at(u).endpoint, origistartp) && origiIndex == pointTopl[v].second.at(u).curveindex)
					//if (pointTopl[v].second.at(u).endpoint == origistartp && origiIndex == pointTopl[v].second.at(u).curveindex)
					{
						pointTopl[v].second[u].curvetype = 'L';// 'l';
						if (isPointRoughEqual(pointTopl[v].first , cbcurve.p0))//if (pointTopl[v].first == cbcurve.p0)
							pointTopl[v].second[u].curveindex = polyline_vec.size() - 2;
						else
							pointTopl[v].second[u].curveindex = polyline_vec.size() - 1;
						pointTopl[v].second[u].endpoint = midpoint;
						break;
					}
				}
				break;
			}
		}
	}
	else if (tp.curvetype == 'Q')//'q'
	{
		quaBezier qbcurve = quadraticBezier_vec[tp.curveindex];
		double p1[2] = { qbcurve.p0.x, qbcurve.p0.y };
		double p2[2] = { qbcurve.p1.x, qbcurve.p1.y };
		double p3[2] = { qbcurve.p2.x, qbcurve.p2.y };
		vector<double> split_t;
		split_t.push_back(0.5);
		split_t.push_back(1.0);
		vector<mypoint2f> seg_curve1;
		double a1, a2, a3;
		double start_t = 0;
		mypoint2f apoint, midpoint;
		for (unsigned int i = 0; i < split_t.size(); ++i)
		{
			for (double t = start_t; t < split_t[i]; t += 0.02)
			{
				a1 = pow((1 - t), 2);
				a2 = 2 * t*(1 - t);
				a3 = t*t;
				apoint.x = a1*p1[0] + a2*p2[0] + a3*p3[0];
				apoint.y = a1*p1[1] + a2*p2[1] + a3*p3[1];
				seg_curve1.push_back(apoint);
			}
			a1 = pow((1 - split_t[i]), 2);
			a2 = 2 * split_t[i] * (1 - split_t[i]);
			a3 = pow(split_t[i], 2);
			apoint.x = a1*p1[0] + a2*p2[0] + a3*p3[0];
			apoint.y = a1*p1[1] + a2*p2[1] + a3*p3[1];
			seg_curve1.push_back(apoint);
			
			PolyLine apl;
			apl.pointlist.insert(apl.pointlist.end(), seg_curve1.begin(), seg_curve1.end());
			apl.origi_type = 'Q'; apl.origi_index = tp.curveindex;
			polyline_vec.push_back(apl);
			seg_curve1.swap(vector<mypoint2f>());
			start_t = split_t[i];
		}
		midpoint = polyline_vec.at(polyline_vec.size() - 2).pointlist.back();
		//mycircle acir;///////////////////////////添加点-》circle_vec
		//acir.cx = midpoint.x;
		//acir.cy = midpoint.y;
		//acir.radius = 1;
		//circle_vec.push_back(acir);	//pointTopl.push_back()		
		vector<TerminalPoint> terpoints;
		TerminalPoint addtp;
		addtp.endpoint = qbcurve.p0;
		addtp.curveindex = polyline_vec.size() - 2;
		addtp.curvetype = 'L';
		terpoints.push_back(addtp);
		addtp.endpoint = qbcurve.p2;
		addtp.curveindex = polyline_vec.size() - 1;
		addtp.curvetype = 'L';
		terpoints.push_back(addtp);
		pointTopl.push_back(make_pair(midpoint, terpoints));
		//modify pointTopl[i].second[k]...
		mypoint2f origistartp = pointTopl[i_ptopl].first;
		mypoint2f origiendp = pointTopl[i_ptopl].second[k].endpoint;
		int origiIndex = pointTopl[i_ptopl].second[k].curveindex;
		pointTopl[i_ptopl].second[k].curvetype = 'L';
		if (isPointRoughEqual(origistartp , qbcurve.p0))//if (origistartp == qbcurve.p0)
			pointTopl[i_ptopl].second[k].curveindex = polyline_vec.size() - 2;
		else
			pointTopl[i_ptopl].second[k].curveindex = polyline_vec.size() - 1;
		pointTopl[i_ptopl].second[k].endpoint = midpoint;
		for (unsigned int v = 0; v < pointTopl.size(); ++v)
		{
			if (isPointRoughEqual(pointTopl[v].first, origiendp))//if (pointTopl[v].first == origiendp)
			{
				for (unsigned int u = 0; u < pointTopl[v].second.size(); ++u)
				{
					//if (pointTopl[v].second.at(u).endpoint == origistartp && origiIndex == pointTopl[v].second.at(u).curveindex)
					if (isPointRoughEqual(pointTopl[v].second.at(u).endpoint, origistartp) && origiIndex == pointTopl[v].second.at(u).curveindex)
					{
						pointTopl[v].second[u].curvetype = 'L';
						if (isPointRoughEqual(pointTopl[v].first , qbcurve.p0))//if (pointTopl[v].first == qbcurve.p0)
							pointTopl[v].second[u].curveindex = polyline_vec.size() - 2;
						else
							pointTopl[v].second[u].curveindex = polyline_vec.size() - 1;
						pointTopl[v].second[u].endpoint = midpoint;
						break;
					}
				}
				break;
			}
		}
	}
	else if (tp.curvetype == 'L')//'l'
	{
		int origi_index = polyline_vec[tp.curveindex].origi_index;
		char origi_tp = polyline_vec[tp.curveindex].origi_type;
		vector<mypoint2f> tp_polyline = polyline_vec[tp.curveindex].pointlist;	
		PolyLine anewpl;
		//vector<mypoint2f> seg_curve;
		mypoint2f midpoint;
		vector<TerminalPoint> terpoints;
		if (tp_polyline.size() <= 2)
		{
			mypoint2f origistartp = pointTopl[i_ptopl].first;
			mypoint2f origiendp = pointTopl[i_ptopl].second[k].endpoint;
			int origiIndex = pointTopl[i_ptopl].second[k].curveindex;
			//vector<vector<mypoint2f>>::iterator iter_poly = polylineVec.begin()+origiIndex;
			//polylineVec.erase(iter_poly);
			vector<TerminalPoint>::iterator iter_ptp = pointTopl[i_ptopl].second.begin() + k;
			pointTopl[i_ptopl].second.erase(iter_ptp);
			for (unsigned int v = 0; v < pointTopl.size(); ++v)
			{
				if (isPointRoughEqual(pointTopl[v].first, origiendp))//if (pointTopl[v].first == origiendp)
				{
					for (unsigned int u = 0; u < pointTopl[v].second.size(); ++u)
					{
						if (isPointRoughEqual(pointTopl[v].second.at(u).endpoint , origistartp) && origiIndex == pointTopl[v].second.at(u).curveindex)
						{
							iter_ptp = pointTopl[v].second.begin() + u;
							pointTopl[v].second.erase(iter_ptp);
							break;
						}
					}
					break;
				}
			}
			return;
		}
		for (unsigned int j = 0; j <= tp_polyline.size() / 2; ++j)
		{
			anewpl.pointlist.push_back(tp_polyline[j]);//seg_curve.push_back(tp_polyline[j]);
		}
		midpoint = anewpl.pointlist.back();//seg_curve.back();
		anewpl.origi_index = origi_index;
		anewpl.origi_type = origi_tp;
		polyline_vec.push_back(anewpl);
		//mycircle acir;///////////////////////////添加点-》circle_vec
		//acir.cx = midpoint.x;
		//acir.cy = midpoint.y;
		//acir.radius = 1;
		//circle_vec.push_back(acir);
		TerminalPoint addtp1;
		addtp1.endpoint = anewpl.pointlist.front();//seg_curve.front();
		addtp1.curveindex = polyline_vec.size() - 1;
		addtp1.curvetype = 'L';
		terpoints.push_back(addtp1);
		anewpl.pointlist.swap(vector<mypoint2f>());//seg_curve.swap(vector<mypoint2f>());
		for (unsigned int j = tp_polyline.size() / 2; j <tp_polyline.size(); ++j)
		{
			anewpl.pointlist.push_back(tp_polyline[j]);//seg_curve.push_back(tp_polyline[j]);
		}
		anewpl.origi_index = origi_index;
		anewpl.origi_type = origi_tp;
		polyline_vec.push_back(anewpl);
		TerminalPoint addtp2;
		addtp2.endpoint = anewpl.pointlist.back();//seg_curve.back();
		addtp2.curveindex = polyline_vec.size() - 1;
		addtp2.curvetype = 'L';
		terpoints.push_back(addtp2);
		pointTopl.push_back(make_pair(midpoint, terpoints));
		//modify pointTopl[i].second[k]...
		mypoint2f origistartp = pointTopl[i_ptopl].first;
		mypoint2f origiendp = pointTopl[i_ptopl].second[k].endpoint;
		int origiIndex = pointTopl[i_ptopl].second[k].curveindex;
		pointTopl[i_ptopl].second[k].curvetype = 'L';
		if (origiendp == polyline_vec.back().pointlist.back())
			pointTopl[i_ptopl].second[k].curveindex = polyline_vec.size() - 2;
		else
			pointTopl[i_ptopl].second[k].curveindex = polyline_vec.size() - 1;
		pointTopl[i_ptopl].second[k].endpoint = midpoint;
		for (unsigned int v = 0; v < pointTopl.size(); ++v)
		{
			if (isPointRoughEqual(pointTopl[v].first, origiendp))//if (pointTopl[v].first == origiendp)
			{
				for (unsigned int u = 0; u < pointTopl[v].second.size(); ++u)
				{
					if (isPointRoughEqual(pointTopl[v].second.at(u).endpoint, origistartp) && origiIndex == pointTopl[v].second.at(u).curveindex)
					{
						pointTopl[v].second[u].curvetype = 'L';
						if(isPointRoughEqual(pointTopl[v].first, polyline_vec.back().pointlist.back()))//if (pointTopl[v].first == tp_polyline.front())
							pointTopl[v].second[u].curveindex = polyline_vec.size() - 1;//pointTopl[v].second[u].curveindex = polyline_vec.size() - 2;
						else
							pointTopl[v].second[u].curveindex = polyline_vec.size() - 2;//pointTopl[v].second[u].curveindex = polyline_vec.size() - 1;
						pointTopl[v].second[u].endpoint = midpoint;
						break;
					}
				}
				break;
			}
		}
	}
}

bool MyGraph::isPointRoughEqual(mypoint2f p1, mypoint2f p2)
{
	if (LineSegmentLength(p1, p2)<eps_isOnsegment)
	//if (fabs(p1.x - p2.x) < errorDeviation && fabs(p1.y - p2.y) < errorDeviation)
		return true;
	else
	{
		return false;
	}
}
void MyGraph::createGraph()//undirected graphs
{
	sgraph.arcNum = 0;
	sgraph.vertexNum = pointTopl.size();
	//sgraph.nodeList = new VertexNode[sgraph.vertexNum];
	TerminalPoint tp;
	for (int i = 0; i < sgraph.vertexNum; ++i)
	{
		VertexNode vnode;
		vnode.position = pointTopl[i].first;
		vnode.verIndex = i;
		vnode.firstEdge = NULL;
		sgraph.nodeList.push_back(vnode);
		/*sgraph.nodeList[i].position = pointTopl[i].first;
		sgraph.nodeList[i].verIndex = i;
		sgraph.nodeList[i].firstEdge = NULL;*/
		if (pointTopl[i].second.size()>0)//有的=0是因为去掉了度为2的点
		{
			for (unsigned int j = 0; j < pointTopl[i].second.size(); ++j)
			{
				tp = pointTopl[i].second[j];
				int pindex = getAdjVerIndex(tp.endpoint);
				if (pindex >= 0)
				{
					ArcNode *anode = new ArcNode;
					anode->adjVertex = pindex;
					anode->curveIndex = tp.curveindex;
					anode->curveType = tp.curvetype;
					anode->lineFlag = 0;  ////0:初始值/围城patch的边 1：innerline  2:outline  3:isolatedline
					anode->next = sgraph.nodeList[i].firstEdge;
					sgraph.nodeList[i].firstEdge = anode;
				}
			}
		}
	}
	cout << "create Graph done" << endl;
}
int MyGraph::getAdjVerIndex(mypoint2f point)
{
	int pointindex = -1;
	for (unsigned int i = 0; i < pointTopl.size(); ++i)
	{
		if (isPointRoughEqual(pointTopl[i].first, point))//if (LineSegmentLength(pointTopl[i].first, point)<eps)//
		{
			pointindex = (int)i;
			break;
		}
	}
	return pointindex;
}

void MyGraph::findLines()
{
	//fond lines from leaf nodes,remove leaf on by one 先找outline
	vector<int> alinepath;//vector<mypoint2f> apointlist;
	bool *visited = new bool[sgraph.vertexNum];
	for (int i = 0; i < sgraph.vertexNum; ++i)
	{
		visited[i] = false;
	}
	//原方法一，找到的outline/attachlines是不连续的,arc的lineflag被大体分为两类，标记为2，0，在findpatch之后需要细分
	for (int i = 0; i < sgraph.vertexNum; ++i)
	{
		if (countEdgeNum(i, 1) == 1 && visited[i] ==false) //第二个参数counttype=0表示计算所有的arcnode ，=1表示除了innerline outline isolateline其他的arc个数
		{
			stack<int> treestack;
			treestack.push(i);
			visited[i] = true;
			while (!treestack.empty())
			{
				int verindex = treestack.top();
				treestack.pop();
				int edgecount = countEdgeNum(verindex, 1); //= 1表示除了innerline outline isolateline其他的arc个数,即只计算linflag=0的arc!!
				if (edgecount == 1)
				{
					int delVerIndex = 0;
					ArcNode *edgenode = sgraph.nodeList[verindex].firstEdge;
					while (edgenode)
					{
						if (visited[edgenode->adjVertex] == false && edgenode->lineFlag == 0)//edgenode->innerLine==0!!!
						{
							treestack.push(edgenode->adjVertex);
							visited[edgenode->adjVertex] = true;
							delVerIndex = edgenode->adjVertex;
						}
						edgenode = edgenode->next;
					}
					alinepath.push_back(verindex);
					//delVerIndex = sgraph.nodeList[verindex].firstEdge->adjVertex;
					removeArcNode(verindex, delVerIndex,2); //2表示只修改标志位为2(attachline)，实际不删除，相应的标志位：0:初始值 1：innerline  2:outline  3:isolatedline 4:不是一条线，是围成patch的一条边
					removeArcNode(delVerIndex, verindex,2);//attachline里包含了 isolateline，后期还需要分开！
				}
				else
				{
					alinepath.push_back(verindex);
					visited[verindex] = false;//!!
					break;
				}
			}
			ICurve plrcurve;
			plrcurve.pt_vec = alinepath;
			plrcurve.type = 0;
			plrGraph.attachlines.push_back(plrcurve);//plrGraph.outlines[i]中记录的是sgraph.nodelist[]中的下标，（从点到点可以找到对应的curvetype curveindex）
			alinepath.swap(vector<int>());
		}
	}

	////check repeat line
	//if (!plrGraph.outlines.empty())
	//{
	//	for (unsigned int i = 0; i < plrGraph.outlines.size() - 1; ++i)
	//	{
	//		vector<int> l1 = plrGraph.outlines[i];
	//		for (unsigned int j = i + 1; j < plrGraph.outlines.size(); ++j)
	//		{
	//			vector<int> l2 = plrGraph.outlines[j];
	//			for (unsigned int k = 1; k < l1.size(); ++k)
	//			{
	//				for (unsigned int h = 1; h < l2.size(); ++h)
	//				{
	//					if (l1[k - 1] == l2[h - 1] && l1[k] == l2[h])
	//					{
	//						cout << "same segment" << endl;
	//					}
	//				}
	//			}
	//			reverse(l2.begin(), l2.end());
	//			for (unsigned int k = 1; k < l1.size(); ++k)
	//			{
	//				for (unsigned int h = 1; h < l2.size(); ++h)
	//				{
	//					if (l1[k - 1] == l2[h - 1] && l1[k] == l2[h])
	//					{
	//						cout << "same segment" << endl;
	//					}
	//				}
	//			}
	//		}
	//	}
	//}

	delete[] visited;
	visited = NULL;
	cout << "find line done" << endl;
}
int MyGraph::countEdgeNum(int startVerIndex, int flag)//第二个参数counttype=0表示计算所有的arcnode ，=1表示除了innerline outline isolateline其他的arc个数
{
	int count = 0;
	if (flag == 0)//计算所有arc个数
	{
		ArcNode *edgenode = sgraph.nodeList[startVerIndex].firstEdge;
		while (edgenode)
		{
			//if (edgenode->loopEdgeflag==false)
			count = count + 1;
			edgenode = edgenode->next;
		}
	}
	else
	{
		ArcNode *edgenode = sgraph.nodeList[startVerIndex].firstEdge;
		while (edgenode)
		{
			if (edgenode->lineFlag == 0)  //0:初始值 1：innerline  2:outline  3:isolatedline 
				count = count + 1;
			edgenode = edgenode->next;
		}
	}
	return count;
}
ArcNode *MyGraph::getArcNode(int startVerIndex, int endVerIndex)
{
	ArcNode *arcnode = sgraph.nodeList[startVerIndex].firstEdge;
	while (arcnode)
	{
		if (arcnode->adjVertex == endVerIndex)
		{
			break;
		}
		arcnode = arcnode->next;
	}
	return arcnode;
}
void MyGraph::removeArcNode(int startVerIndex, int endVerIndex, int flag)
{
	//int cc = countEdgeNum(startVerIndex);
	//int bb = countEdgeNum(endVerIndex);
	if (flag == 0)//0表示真正的删除
	{
		ArcNode *arcnode = sgraph.nodeList[startVerIndex].firstEdge;
		ArcNode *prior_arc = NULL;
		while (arcnode)
		{
			if (arcnode->adjVertex == endVerIndex)
			{
				//prior_arc->next = arcnode->next;
				if (prior_arc == NULL)
				{
					sgraph.nodeList[startVerIndex].firstEdge = arcnode->next;
				}
				else
				{
					prior_arc->next = arcnode->next;
				}
				delete arcnode;
				arcnode = NULL;
				break;
			}
			prior_arc = arcnode;
			arcnode = arcnode->next;
		}
		//delete orderarc
		if (!sgraph.nodeList[startVerIndex].orderedArcs.empty())
		{
			vector<ArcNode>::iterator iter_arc = sgraph.nodeList[startVerIndex].orderedArcs.begin();
			while (iter_arc != sgraph.nodeList[startVerIndex].orderedArcs.end())
			{
				if (iter_arc->adjVertex == endVerIndex)
				{
					iter_arc = sgraph.nodeList[startVerIndex].orderedArcs.erase(iter_arc);
				}
				else
					iter_arc++;
			}
		}
	}
	else//用flag修改标志位
	{
		ArcNode *arcnode = sgraph.nodeList[startVerIndex].firstEdge;
		while (arcnode)
		{
			if (arcnode->adjVertex == endVerIndex)
			{
				arcnode->lineFlag = flag;
				break;
			}
			arcnode = arcnode->next;
		}
		//delete orderarc
		if (!sgraph.nodeList[startVerIndex].orderedArcs.empty())
		{
			vector<ArcNode>::iterator iter_arc = sgraph.nodeList[startVerIndex].orderedArcs.begin();
			while (iter_arc != sgraph.nodeList[startVerIndex].orderedArcs.end())
			{
				if (iter_arc->adjVertex == endVerIndex)
				{
					iter_arc->lineFlag = flag;
					break;
				}
				iter_arc++;
			}
		}
	}
}

void MyGraph::divideSubGraph()
{
	int *visitStatus = new int[sgraph.vertexNum];
	//0:unvisited,1:push in queue,2: out of queue
	for (int i = 0; i < sgraph.vertexNum; ++i)
	{
		visitStatus[i] = 0;
	}
	for (int i = 0; i < sgraph.vertexNum; ++i)//breadth-first traversal  BFT
	{
		if (visitStatus[i] == 0)//if (visitStatus[i] == 0 && countEdgeNum(i,1)>0)//countEdgeNum只计算arc.innerline=0的arc
		{		
			int arcCount = 0;
			int verCount = 0;
			int tmp_ver = 0;
			queue<int> verQueue;
			verQueue.push(i);
			visitStatus[i] = 1;
			while (!verQueue.empty())
			{
				tmp_ver = verQueue.front();
				verQueue.pop();
				verCount = verCount + 1;
				visitStatus[tmp_ver] = 2;
				ArcNode *tmp_arc = sgraph.nodeList[tmp_ver].firstEdge;
				while (tmp_arc)
				{
					//if (tmp_arc->innerLine == 0)
					//{
						if (visitStatus[tmp_arc->adjVertex] == 0)
						{
							verQueue.push(tmp_arc->adjVertex);
							visitStatus[tmp_arc->adjVertex] = 1;
							arcCount = arcCount + 1;
						}
						else if (visitStatus[tmp_arc->adjVertex] == 1)
						{
							arcCount = arcCount + 1;
						}
					//}							
					tmp_arc = tmp_arc->next;
				}
			}
			cout << "this subgraph's arc number=" << arcCount << ",vertex number=" << verCount << endl;
			//if (loopflag==false)
			//	subgraph_map.insert(map<int, string>::value_type(i, "tree"));
			//else
			//	subgraph_map.insert(map<int, string>::value_type(i, "loop"));
			if (verCount - 1 == arcCount)
			{
				subgraph_map.insert(map<int, string>::value_type(i, "tree"));
			}
			else
			{
				subgraph_map.insert(map<int, string>::value_type(i, "loop"));
			}
		}
	}
	delete[]visitStatus;
	visitStatus = NULL;
	cout << "subgraph_map=" << subgraph_map.size() << endl;
	map<int, string>::iterator iter = subgraph_map.begin();
	for (iter; iter != subgraph_map.end(); ++iter)
	{
		cout <<"start index:"<< iter->first << "," <<"type:"<< iter->second.c_str() << endl;
	}
	cout << "divideSubGraph done" << endl;
}
void MyGraph::MakeOrderedArcs()
{
	cout << "MakeOrderOfArcs" << endl;
	//store ordered arcs(edge) of every points
	for (int i = 0; i < sgraph.vertexNum; ++i)//breadth-first traversal  BFT
	{
		int edgenum = countEdgeNum(i,1);
		if (edgenum>2)
		{
			mypoint2f startp = sgraph.nodeList[i].position;
			vector<pair<int,mypoint2f>> unitVector;  //int 表示arcnode->adjVertex
			ArcNode *anode_p = sgraph.nodeList[i].firstEdge;
			while (anode_p)//unitVector记录所有相邻arc的单位方向
			{
				if (anode_p->lineFlag == 0)
				{
					mypoint2f endp = sgraph.nodeList[anode_p->adjVertex].position;
					vector<mypoint2f> cur_plist;
					getPListOfaCurve(startp, endp, i, anode_p->adjVertex, &cur_plist,0);//从点i到 anode_p->adjVertex
					if (cur_plist.size() > 2)
					{
						endp = cur_plist[1];//endp = cur_plist[2];若两条线爱的很近，取[2]会出错...
					}
					mypoint2f norvec = normalization(startp, endp);//方向向量 归一化
					unitVector.push_back(make_pair(anode_p->adjVertex, norvec));
				}
				anode_p = anode_p->next;
			}
			vector<pair<int, double>> uv_left;
			vector<pair<int, double>> uv_right;
			for (unsigned int j = 1; j < unitVector.size(); ++j)//以unitVector[0] 为中心线，分左右两组
			{			
				if (cross(unitVector[0].second, unitVector[j].second)>0)//unitVector[j].second在左边
				{
					double costheta=getCosineAngle(mypoint2f(0, 0), unitVector[0].second, mypoint2f(0, 0), unitVector[j].second);
					if (uv_left.empty())
						uv_left.push_back(make_pair(unitVector[j].first,costheta));
					else
					{
						vector<pair<int, double>>::iterator iter_l= uv_left.begin();
						bool pushflag = false;
						while (iter_l != uv_left.end())//theta从小到大排序，即costheta从大到小
						{
							if (costheta > iter_l->second)
							{
								uv_left.insert(iter_l, make_pair(unitVector[j].first, costheta));
								pushflag = true;
								break;
							}						
							iter_l++;
						}
						if (pushflag == false)
							uv_left.push_back(make_pair(unitVector[j].first, costheta));
					}
				}				
				else//unitVector[j].second在右边，操作相似，但costheta从小到大排序
				{
					double costheta = getCosineAngle(mypoint2f(0, 0), unitVector[0].second, mypoint2f(0, 0), unitVector[j].second);
					if (uv_right.empty())
						uv_right.push_back(make_pair(unitVector[j].first, costheta));
					else
					{
						vector<pair<int, double>>::iterator iter_r = uv_right.begin();
						bool pushflag = false;
						while (iter_r != uv_right.end())//theta从小到大排序，即costheta从大到小
						{
							if (costheta < iter_r->second)
							{
								uv_right.insert(iter_r, make_pair(unitVector[j].first, costheta));
								pushflag = true;
								break;
							}
							iter_r++;
						}
						if (pushflag == false)
							uv_right.push_back(make_pair(unitVector[j].first, costheta));
					}
				}
			}			
			ArcNode *arcnd_p = getArcNode(i, unitVector[0].first);//sgraph.nodeList[i].firstEdge;//先将中心线unitVector[0]存入sgraph.nodeList[i].orderedArcs
			ArcNode arcnd;
			arcnd.adjVertex = arcnd_p->adjVertex;
			arcnd.curveIndex = arcnd_p->curveIndex;
			arcnd.curveType = arcnd_p->curveType;
			arcnd.lineFlag = 0;
			arcnd.next = NULL;
			sgraph.nodeList[i].orderedArcs.push_back(arcnd);
			arcnd_p = arcnd_p->next;
			for (unsigned int k = 0; k < uv_left.size(); ++k)
			{
				arcnd_p = sgraph.nodeList[i].firstEdge;
				while (arcnd_p)
				{
					if (arcnd_p->adjVertex == uv_left[k].first)
					{
						ArcNode tmparc;
						tmparc.adjVertex = arcnd_p->adjVertex;
						tmparc.curveIndex = arcnd_p->curveIndex;
						tmparc.curveType = arcnd_p->curveType;
						tmparc.lineFlag = 0;
						tmparc.next = NULL;
						sgraph.nodeList[i].orderedArcs.push_back(tmparc);
						break;
					}
					arcnd_p = arcnd_p->next;
				}				
			}
			//arcnd_p = sgraph.nodeList[i].firstEdge->next;
			for (unsigned int k = 0; k < uv_right.size(); ++k)
			{
				arcnd_p = sgraph.nodeList[i].firstEdge;
				while (arcnd_p)
				{
					if (arcnd_p->adjVertex == uv_right[k].first)
					{
						ArcNode tmparc;
						tmparc.adjVertex = arcnd_p->adjVertex;
						tmparc.curveIndex = arcnd_p->curveIndex;
						tmparc.curveType = arcnd_p->curveType;
						tmparc.lineFlag = 0;
						tmparc.next = NULL;
						sgraph.nodeList[i].orderedArcs.push_back(tmparc);
						break;
					}
					arcnd_p = arcnd_p->next;
				}
			}
			if (sgraph.nodeList[i].orderedArcs.size() != edgenum)
				cout << "error" << endl;
		}
		else if (edgenum == 2)
		{
			ArcNode *anode_p = sgraph.nodeList[i].firstEdge;
			while (anode_p)
			{
				if (anode_p->lineFlag == 0)
				{
					ArcNode tmpnode;
					tmpnode.adjVertex = anode_p->adjVertex;
					tmpnode.curveIndex = anode_p->curveIndex;
					tmpnode.curveType = anode_p->curveType;
					tmpnode.lineFlag = 0;
					tmpnode.next = NULL;
					sgraph.nodeList[i].orderedArcs.push_back(tmpnode);				
				}
				anode_p = anode_p->next;
			}
		}
	}
	cout << "MakeOrderOfArcs done" << endl;
}

void MyGraph::findPatches_new()
{
	cout << endl;
	cout << "findPatches_new" << endl;
	int innerline_count = 0;
	int overlap_count = 0;
	int patchCount = 0;
	int lastpNum = 0;//累加记录这一个子图中找到的patch的数量
	vector<int> patch_segment;
	int *visitStatus = new int[sgraph.vertexNum];	//0:unvisited,1:push in queue,2: out of queue
	for (int i = 0; i < sgraph.vertexNum; ++i)
	{
		visitStatus[i] = 0;
	}
	map<int, string>::iterator iter = subgraph_map.begin();
	/*1. 一个while循环是对一个子图的处理，找完patch后，接着删除掉覆盖小patch的那个大patch*/
	while (iter != subgraph_map.end())
	{
		if (iter->second == "tree")
		{
			iter++;
			continue;
		}
		int graph_p = iter->first;
		bool ft_flag = false;//标记是否有innerline，并不记录innerline，因为method()中找到的innerline可能不完整
		vector<int> innerPoints;//记录该子图中的innerPoints(注意删除重复的) ,,, 没有用可以删除
		int tmp_ver = 0;
		queue<int> verQueue;
		verQueue.push(graph_p);
		visitStatus[graph_p] = 1;
		while (!verQueue.empty())
		{
			tmp_ver = verQueue.front();
			verQueue.pop();
			visitStatus[tmp_ver] = 2;			
			ArcNode *tmp_arc = sgraph.nodeList[tmp_ver].firstEdge;
			while (tmp_arc)
			{
				if (tmp_arc->lineFlag == 0)
				{
					if (visitStatus[tmp_arc->adjVertex] != 2)
					{
						if (visitStatus[tmp_arc->adjVertex] == 0)
						{
							verQueue.push(tmp_arc->adjVertex);
							visitStatus[tmp_arc->adjVertex] = 1;
						}
						method(tmp_ver, tmp_arc->adjVertex, patchCount, ft_flag, innerPoints);//ft_flag=true返回有innerline存在
						method(tmp_arc->adjVertex, tmp_ver, patchCount, ft_flag, innerPoints);
					}
				}
				else
				{
					//tmp_ver--tmp_arc这条边的类型=2，是一条outline
					if (visitStatus[tmp_arc->adjVertex] != 2)
					{
						if (visitStatus[tmp_arc->adjVertex] == 0)
						{
							verQueue.push(tmp_arc->adjVertex);
							visitStatus[tmp_arc->adjVertex] = 1;
						}
					}
				}
				tmp_arc = tmp_arc->next;
			}
		}
		//在这个子图中,计算patch的周长，面积，中心，n内接圆半径
		vector<pair<int, double>> ordered_patch;//patchIndex vs. perimeter
		for (int i = lastpNum; i < patchCount; i++)
		{
			vector<mypoint2f> total_plist;
			double pmt = 0; double areas = 0;
			for (unsigned int j = 1; j < plrGraph.patches[i].pointIndexOfPatch.size(); ++j)
			{
				int verindex1 = plrGraph.patches[i].pointIndexOfPatch.at(j - 1);
				int verindex2 = plrGraph.patches[i].pointIndexOfPatch.at(j);
				mypoint2f point1 = sgraph.nodeList[verindex1].position;//=lineWithCommonPoints[verindex1].first
				mypoint2f point2 = sgraph.nodeList[verindex2].position;
				if (!total_plist.empty())
					total_plist.pop_back();
				getPListOfaCurve(point1, point2, verindex1, verindex2, &total_plist,0);
			}
			double aver_x = 0;
			double aver_y = 0;
			for (unsigned int j = 0; j < total_plist.size(); ++j)
			{
				aver_x = aver_x + total_plist[j].x;
				aver_y = aver_y + total_plist[j].y;
			}
			aver_x = aver_x / (double)total_plist.size();
			aver_y = aver_y / (double)total_plist.size();
			plrGraph.patches[i].centre_position = mypoint2f(aver_x, aver_y);//计算中心

			pmt = getPerimeterofVec(&total_plist);//周长
			plrGraph.patches[i].perimeter = pmt;
			areas = GetPolygonArea(&total_plist);//面积
			if (areas >= 0)//若面积>=0 说明是正常patch，应为patch都是沿同一个顺序的（逆时针）
				plrGraph.patches[i].area = fabs(areas);
			else//否则是错误的patch，即repeat patch或overlay patch
			{
				plrGraph.patches[i].type = -1;
				deleOverlayPatches(i, 0);
				overlap_count++;
			}
			/*计算内接圆 圆心,有错误*/
				//if (areas >= 0)
				//{
				//	total_plist.pop_back();//删除最后一个点，现在total_plist  首！=尾
				//	mypoint2f incir_center(0, 0);
				//	double incir_r = 0;
				//	maxInscribedCircle(&total_plist, incir_center, incir_r);//total_plist可能改变了，数据点增多
				//	plrGraph.patches[i].centre_position = incir_center;
				//	plrGraph.patches[i].incircle_r = incir_r;
				//	//添加组成patch的节点， 
				//	mycircle acir;
				//	acir.cx = incir_center.x; acir.cy = incir_center.y;
				//	acir.radius = incir_r;
				//	circle_vec.push_back(acir);
				//}
			
			/*另一种计算中心的方法 （cx,cy）公式： cx=1/(6A)* 累加（Xi+Xi+1）(XiYi+1-Xi+1Yi)， cy=1/(6A)* 累加（Yi+Yi+1）(XiYi+1-Xi+1Yi)，A是面积（正负都可以）*/
			//注意区别，计算时需要先把svg坐标转换成opengl坐标！！算得的中心更准一些？
			//double add_x = 0;
			//double add_y = 0;
			//for (unsigned int j = 1; j < total_plist.size(); ++j)
			//{
			//	mypoint2f tmp_p1(total_plist[j - 1].x, file_height - total_plist[j - 1].y);
			//	mypoint2f tmp_p2(total_plist[j].x, file_height - total_plist[j].y);
			//	double tmp_x = (tmp_p1.x + tmp_p2.x)*cross(tmp_p1, tmp_p2);
			//	double tmp_y = (tmp_p1.y + tmp_p2.y)*cross(tmp_p1, tmp_p2);
			//	add_x = add_x + tmp_x;
			//	add_y = add_y + tmp_y;
			//	/*double tmp_x = (total_plist[j - 1].x + total_plist[j].x)*cross(total_plist[j - 1], total_plist[j]);
			//	double tmp_y = (total_plist[j - 1].y + total_plist[j].y)*cross(total_plist[j - 1], total_plist[j]);
			//	add_x = add_x + tmp_x;
			//	add_y = add_y + tmp_y;*/
			//}
			//double fact = 1 / (6 * areas);
			//add_x = add_x*fact;
			//add_y = add_y*fact;
			//plrGraph.patches[i].centre_position = mypoint2f(add_x, add_y);

			////以上 是计算周长和面积，下面是对patch的周长进行从大到小排序（可以删除）
			//bool inflag = false;
			//vector<pair<int, double>>::iterator iter_op = ordered_patch.begin();//int：patch的编号，double:patch的周长
			//while (iter_op != ordered_patch.end())
			//{
			//	if (pmt > iter_op->second)//从大到小排序
			//	{
			//		ordered_patch.insert(iter_op, make_pair(i, pmt));
			//		inflag = true;
			//		break;
			//	}
			//	iter_op++;
			//}
			//if (inflag == false)
			//	ordered_patch.push_back(make_pair(i, pmt));
		}
		if (ft_flag == true)//说明有patch-line-patch情况，findInnerLine()是找到这种line,设置为1（原来为0）
		{
			//在这个子图中，找patch-line-patch情况中的line，返回inner line的个数
			int il_count = findInnerLine(graph_p);
			innerline_count = innerline_count + il_count;
		}
		lastpNum = patchCount;//!!!!
		iter++;
	}
	//2.删除重复的patch(),应该没有重复的，该函数可以删除
	int delepatch = deleRepeatPatch();
	//3.然后把有向的patch补全, 这里主要是指边界处的patch！！不是内部的。rever_node->patchesIndex.push_back(pch_index)+sgraph.nodeList[anode->adjVertex].orderedArcs[j].patchesIndex.push_back(pch_index);  
	int ndsize = sgraph.nodeList.size();
	for (int i = 0; i < ndsize; ++i)
	{
		ArcNode* anode = sgraph.nodeList[i].firstEdge;
		while (anode)
		{
			if (anode->patchesIndex.size() == 1)//先找到i中patchsize=1的arcnode
			{
				ArcNode* rever_node = getArcNode(anode->adjVertex, i);
				if (rever_node->patchesIndex.size() == 0)
				{
					int pch_index = anode->patchesIndex[0];
					//plrGraph.patches[pch_index].type = 1;
					rever_node->patchesIndex.push_back(pch_index);//!!!!
					for (unsigned int j = 0; j < sgraph.nodeList[anode->adjVertex].orderedArcs.size(); ++j)
					{
						if (sgraph.nodeList[anode->adjVertex].orderedArcs[j].adjVertex == i && sgraph.nodeList[anode->adjVertex].orderedArcs[j].patchesIndex.empty())
						{
							sgraph.nodeList[anode->adjVertex].orderedArcs[j].patchesIndex.push_back(pch_index);
						}
					}
				}
			}
			anode = anode->next;
		}
	}
	//4.找完patch后 才可以划分isolatedline和attachline
	//划分outline 和isolateline，在findline()中都被标记为2，现在将两者分开。plrGraph.outlines[i]中记录的是sgraph.nodelist[]中的下标，（从点到点可以找到对应的curvetype curveindex）
	//（1）isolateline arc标记为3，直接存入plrGraph.isolatedlines.push_back （2）outline需要重新trace	重新放入一个新的变量里vector<vector<int>> tmp_outlines 
	//(3)把patch 和line 结合起来patch 结构体中的<pair> point_line
	divideAttachLines();

	delete[]visitStatus;
	visitStatus = NULL;

	////check
	//for (unsigned int i = 0; i < plrGraph.patches.size(); ++i)
	//{
	//	/*if (plrGraph.patches[i].patch_index != i)
	//		cout << "0" << endl;
	//	for (unsigned int j = 1; j < plrGraph.patches[i].pointIndexOfPatch.size(); ++j)
	//	{
	//		if (plrGraph.patches[i].edges_vec[j - 1].startp_index != plrGraph.patches[i].pointIndexOfPatch[j - 1])
	//			cout << "1" << endl;
	//		if (plrGraph.patches[i].edges_vec[j - 1].endp_index != plrGraph.patches[i].pointIndexOfPatch[j])
	//			cout << "2" << endl;		
	//	}*/
	//	for (unsigned int j = 0; j < plrGraph.patches[i].edges_vec.size(); ++j)
	//	{
	//		if (plrGraph.patches[i].type >= 0)
	//		{
	//			if (plrGraph.patches[i].edges_vec[j].adjPatchIndex[0] != i)
	//				cout << "3第一个不是自己的index" << endl;
	//			for (unsigned int k = 0; k < plrGraph.patches[i].edges_vec[j].adjPatchIndex.size(); ++k)
	//			{
	//			int pchindex = plrGraph.patches[i].edges_vec[j].adjPatchIndex[k];
	//			int ss = plrGraph.patches[i].edges_vec[j].startp_index;
	//			int ee = plrGraph.patches[i].edges_vec[j].endp_index;
	//			bool flag = false;
	//			ArcNode* anode = getArcNode(ss, ee);
	//			if (anode == NULL)
	//			cout<<"NULL" << endl;
	//			else
	//			{
	//			for (unsigned int g = 0; g < anode->patchesIndex.size(); ++g)
	//			{
	//			if (anode->patchesIndex[g] == pchindex)
	//			flag = true;
	//			}
	//			if (flag == false)
	//			cout << "4" << endl;
	//			}
	//			flag = false;
	//			anode = getArcNode(ee, ss);
	//			if (anode == NULL)
	//			cout << "NULL" << endl;
	//			else
	//			{
	//			for (unsigned int g = 0; g < anode->patchesIndex.size(); ++g)
	//			{
	//			if (anode->patchesIndex[g] == pchindex)
	//			flag = true;
	//			}
	//			if (flag == false)
	//			cout << "5" << endl;
	//			}
	//			}
	//		}		
	//	}
	//}

	//printpathOfPatches();
	cout << endl;
	cout << "inner line count=" << innerline_count << endl;
	cout << "total lineCount =" << plrGraph.attachlines.size() + plrGraph.innerlines.size() + plrGraph.isolatedlines.size() << endl << endl;
	cout << "total patches:" << plrGraph.patches.size() << endl;
	cout << "overlap patch=" << overlap_count << endl;//closure patch是指最大的覆盖其他小patch的 patch
	cout << "deleted repeat patches=" << delepatch << endl;
	cout << "finally patchCount=" << plrGraph.patches.size() - delepatch - overlap_count << endl;
   	cout << "findPatches_new done" << endl << endl;
}
void MyGraph::method(int pindex1, int pindex2, int &pcount, bool &cflag, vector<int> &inn_pt)
{
	//cflag在method调用之前已被赋值false
	int p1 = pindex1;
	int p2 = pindex2;
	int edgenum = countEdgeNum(p1,1);
	ArcNode *tmp_arc = getArcNode(pindex1, pindex2);//沿着p1->p2这条边的方向，寻找一个patch
	if (tmp_arc==NULL)
		cout << "error" << endl;
	if (edgenum > 1)
	{
		if (tmp_arc->patchesIndex.empty())//???if (tmp_arc->patchesIndex.size()<2)//
		{
			int fromVerindex = p1;
			int nextVerindex = p2;
			vector<int> apatch;
			apatch.push_back(fromVerindex);
			apatch.push_back(nextVerindex);
			while (nextVerindex != p1)
			{
				if (sgraph.nodeList[nextVerindex].orderedArcs.size() >= 2)
				{
					bool ffalg = false;
					ArcNode h_arc = sgraph.nodeList[nextVerindex].orderedArcs.front();
					sgraph.nodeList[nextVerindex].orderedArcs.push_back(h_arc);//把头push进去，构成循环，下面再pop出来
					for (unsigned int j = 0; j < sgraph.nodeList[nextVerindex].orderedArcs.size(); ++j)
					{
						if (sgraph.nodeList[nextVerindex].orderedArcs.at(j).adjVertex == fromVerindex)
						{
							if (sgraph.nodeList[nextVerindex].orderedArcs.at(j + 1).patchesIndex.size() < 2)
							{
								/*if (sgraph.nodeList[nextVerindex].orderedArcs.at(j + 1).innerLine == false)
								{
								int adjV = sgraph.nodeList[nextVerindex].orderedArcs.at(j + 1).adjVertex;
								apatch.push_back(adjV);
								fromVerindex = nextVerindex;
								nextVerindex = adjV;
								sgraph.nodeList[fromVerindex].orderedArcs.pop_back();
								ffalg = true;
								}*/
								int adjV = sgraph.nodeList[nextVerindex].orderedArcs.at(j + 1).adjVertex;
								ArcNode *adjNode = getArcNode(nextVerindex, adjV);
								if (adjNode != NULL)
								{
									apatch.push_back(adjV);
									fromVerindex = nextVerindex;
									nextVerindex = adjV;
									sgraph.nodeList[fromVerindex].orderedArcs.pop_back();
									ffalg = true;
								}
								else
									break;
							}
							break;
						}
					}
					if (ffalg == false)
					{
						sgraph.nodeList[nextVerindex].orderedArcs.pop_back();
						break;//!!
					}
				}
				else
					break;
			}
			if (apatch.size() > 2 && apatch.back() == apatch.front())//判断是不是一个正确的patch
			{
				bool inn_pt_flag = false;
				if (apatch.at(1) == apatch.at(apatch.size() - 2))//line-patch 的情况
				{		
					cflag = true;
					inn_pt_flag = true;
				}
				else//在判断patch-line-patch或者patch-point-patch的情况,并且记录innerpoint
				{
					vector<int> alinepath;
					for (unsigned int g = 1; g < apatch.size() - 1; ++g)
					{
						for (unsigned int g_inn = g + 1; g_inn < apatch.size() - 1; ++g_inn)
						{
							if (apatch[g] == apatch[g_inn])
							{
								inn_pt_flag = true;
								int k = 0;
								while (apatch[g + k] == apatch[g_inn - k] && (g + k) <= (g_inn - k))
								{
									alinepath.push_back(apatch[g + k]);
									++k;
								}
								if (alinepath.size() > 0)
								{
									if (alinepath.size() == 1)//=1说明是一个点，>1是一条线
									{
										inn_pt.push_back(alinepath.back());
									}
									else
									{
										g = g + alinepath.size() - 1;
										cflag = true;
									}
									alinepath.clear();
									break;//跳出循环
								}
							}
						}
					}
				}
				if (inn_pt_flag == false)//没有innerline innerpoint,正确的patch
				{
					//有方向的patch
					Patch ap;
					ap.type = 0;
					ap.pointIndexOfPatch.insert(ap.pointIndexOfPatch.begin(), apatch.begin(), apatch.end());
					ap.patch_index = pcount;
					ap.perimeter = 0;
					ap.area = 0;
					ap.centre_position = mypoint2f(0,0);//ap.incircle_r = 0;
					for (unsigned int i_p = 1; i_p < apatch.size(); ++i_p)
					{
						EdgeOfPatch anedge;
						anedge.startp_index = apatch[i_p - 1];
						anedge.endp_index = apatch[i_p];
						anedge.adjPatchIndex.push_back(pcount);
						anedge.edgeflag = 0;
						ArcNode *tpnode1 = getArcNode(apatch[i_p - 1], apatch[i_p]);//patch的边是有方向的
						tpnode1->patchesIndex.push_back(pcount);
						for (unsigned int i_oa = 0; i_oa < sgraph.nodeList[apatch[i_p - 1]].orderedArcs.size(); ++i_oa)
						{
							if (sgraph.nodeList[apatch[i_p - 1]].orderedArcs[i_oa].adjVertex == apatch[i_p])
							{
								sgraph.nodeList[apatch[i_p - 1]].orderedArcs[i_oa].patchesIndex.push_back(pcount);
								break;
							}
						}
						ArcNode *tpnode2 = getArcNode(apatch[i_p], apatch[i_p - 1]);
						if (tpnode2->patchesIndex.size() == 1)
						{
							anedge.adjPatchIndex.push_back(tpnode2->patchesIndex.front());
							Patch *pth2 = getAPatch(tpnode2->patchesIndex.front());
							for (unsigned int e_i = 0; e_i < pth2->edges_vec.size(); ++e_i)
							{
								if (pth2->edges_vec[e_i].startp_index == apatch[i_p - 1] && pth2->edges_vec[e_i].endp_index == apatch[i_p])
									pth2->edges_vec[e_i].adjPatchIndex.push_back(pcount);
								else if (pth2->edges_vec[e_i].endp_index == apatch[i_p - 1] && pth2->edges_vec[e_i].startp_index == apatch[i_p])
									pth2->edges_vec[e_i].adjPatchIndex.push_back(pcount);
							}
							tpnode2->patchesIndex.push_back(pcount);
							for (unsigned int i_oa = 0; i_oa < sgraph.nodeList[apatch[i_p]].orderedArcs.size(); ++i_oa)
							{
								if (sgraph.nodeList[apatch[i_p]].orderedArcs[i_oa].adjVertex == apatch[i_p - 1])
								{
									sgraph.nodeList[apatch[i_p]].orderedArcs[i_oa].patchesIndex.push_back(pcount);
									break;
								}
							}
							tpnode1->patchesIndex.push_back(tpnode2->patchesIndex.front());
							for (unsigned int i_oa = 0; i_oa < sgraph.nodeList[apatch[i_p - 1]].orderedArcs.size(); ++i_oa)
							{
								if (sgraph.nodeList[apatch[i_p - 1]].orderedArcs[i_oa].adjVertex == apatch[i_p])
								{
									sgraph.nodeList[apatch[i_p - 1]].orderedArcs[i_oa].patchesIndex.push_back(tpnode2->patchesIndex.front());
									break;
								}
							}
						}
						ap.edges_vec.push_back(anedge);
					}
					pcount++;
					plrGraph.patches.push_back(ap);

					////无方向的
					//Patch ap;
					//ap.type = 0;
					//ap.pointIndexOfPatch.insert(ap.pointIndexOfPatch.begin(), apatch.begin(), apatch.end());
					//ap.patch_index = pcount;
					//ap.perimeter = 0;
					//ap.area = 0;
					//ap.centre_position = mypoint2f(0, 0);
					//for (unsigned int i_p = 1; i_p < apatch.size(); ++i_p)
					//{
					//	EdgeOfPatch anedge;
					//	anedge.startp_index = apatch[i_p - 1];
					//	anedge.endp_index = apatch[i_p];
					//	anedge.adjPatchIndex.push_back(pcount);
					//	anedge.edgeflag = 0;
					//	ArcNode *tpnode1 = getArcNode(apatch[i_p - 1], apatch[i_p]);//方向1
					//	tpnode1->patchesIndex.push_back(pcount);
					//	for (unsigned int i_oa = 0; i_oa < sgraph.nodeList[apatch[i_p - 1]].orderedArcs.size(); ++i_oa)
					//	{
					//		if (sgraph.nodeList[apatch[i_p - 1]].orderedArcs[i_oa].adjVertex == apatch[i_p])
					//		{
					//			sgraph.nodeList[apatch[i_p - 1]].orderedArcs[i_oa].patchesIndex.push_back(pcount);
					//			break;
					//		}
					//	}
					//	ArcNode *tpnode2 = getArcNode(apatch[i_p], apatch[i_p - 1]);//方向2
					//	if (tpnode2->patchesIndex.size() == 1)
					//	{
					//		anedge.adjPatchIndex.push_back(tpnode2->patchesIndex.front());
					//		Patch *pth2 = getAPatch(tpnode2->patchesIndex.front());
					//		for (unsigned int e_i = 0; e_i < pth2->edges_vec.size(); ++e_i)
					//		{
					//			if (pth2->edges_vec[e_i].startp_index == apatch[i_p - 1] && pth2->edges_vec[e_i].endp_index == apatch[i_p])
					//				pth2->edges_vec[e_i].adjPatchIndex.push_back(pcount);
					//			else if (pth2->edges_vec[e_i].endp_index == apatch[i_p - 1] && pth2->edges_vec[e_i].startp_index == apatch[i_p])
					//				pth2->edges_vec[e_i].adjPatchIndex.push_back(pcount);
					//		}
					//		tpnode2->patchesIndex.push_back(pcount);
					//		for (unsigned int i_oa = 0; i_oa < sgraph.nodeList[apatch[i_p]].orderedArcs.size(); ++i_oa)
					//		{
					//			if (sgraph.nodeList[apatch[i_p]].orderedArcs[i_oa].adjVertex == apatch[i_p - 1])
					//			{
					//				sgraph.nodeList[apatch[i_p]].orderedArcs[i_oa].patchesIndex.push_back(pcount);
					//				break;
					//			}
					//		}
					//		tpnode1->patchesIndex.push_back(tpnode2->patchesIndex.front());
					//		for (unsigned int i_oa = 0; i_oa < sgraph.nodeList[apatch[i_p - 1]].orderedArcs.size(); ++i_oa)
					//		{
					//			if (sgraph.nodeList[apatch[i_p - 1]].orderedArcs[i_oa].adjVertex == apatch[i_p])
					//			{
					//				sgraph.nodeList[apatch[i_p - 1]].orderedArcs[i_oa].patchesIndex.push_back(tpnode2->patchesIndex.front());
					//				break;
					//			}
					//		}
					//	}
					//	else if (tpnode2->patchesIndex.size() == 0)
					//	{
					//		tpnode2->patchesIndex.push_back(pcount);
					//		for (unsigned int i_oa = 0; i_oa < sgraph.nodeList[apatch[i_p]].orderedArcs.size(); ++i_oa)
					//		{
					//			if (sgraph.nodeList[apatch[i_p]].orderedArcs[i_oa].adjVertex == apatch[i_p - 1])
					//			{
					//				sgraph.nodeList[apatch[i_p]].orderedArcs[i_oa].patchesIndex.push_back(pcount);
					//				break;
					//			}
					//		}
					//	}
					//	ap.edges_vec.push_back(anedge);
					//}
					//pcount++;
					//plrGraph.patches.push_back(ap);
				}
			}
		}
	}	
}
void MyGraph::deleOverlayPatches(int dele_pch, int store_pch)//store_pch没用,注意和deleRepeatPatch()的区别
{
	int dele_pchindex = dele_pch;
	int store_pchindex = store_pch;
	plrGraph.patches[dele_pchindex].type = -1;
	for (unsigned int vec = 0; vec < plrGraph.patches[dele_pchindex].edges_vec.size(); vec++)//要删除的那个patch的edge
	{
		int sp = plrGraph.patches[dele_pchindex].edges_vec[vec].startp_index;
		int ep = plrGraph.patches[dele_pchindex].edges_vec[vec].endp_index;
		if (plrGraph.patches[dele_pchindex].edges_vec[vec].adjPatchIndex.size()>1)//plrGraph.patches[dele_pchindex].的每一个edge
		{
			for (unsigned int i = 0; i < plrGraph.patches[dele_pchindex].edges_vec[vec].adjPatchIndex.size(); ++i)
			{
				if (plrGraph.patches[dele_pchindex].edges_vec[vec].adjPatchIndex[i] != dele_pchindex)
				{
					int another_pch = plrGraph.patches[dele_pchindex].edges_vec[vec].adjPatchIndex[i];
					EdgeOfPatch *e_pch = getBorderofPatch(another_pch, sp, ep);
					if (e_pch->adjPatchIndex.back() == dele_pchindex)
						e_pch->adjPatchIndex.pop_back();
					else if (e_pch->adjPatchIndex.front() == dele_pchindex)
						e_pch->adjPatchIndex.erase(e_pch->adjPatchIndex.begin());//一般不会执行这步
				}
			}
			/*
			//原方法，为了确保删除正确，用了上方法
			int another_pch = plrGraph.patches[dele_pchindex].edges_vec[vec].adjPatchIndex[1];
			EdgeOfPatch *e_pch = getBorderofPatch(another_pch, sp, ep);
			if (e_pch->adjPatchIndex.back()==dele_pchindex)
				e_pch->adjPatchIndex.pop_back();*/
		}
		plrGraph.patches[dele_pchindex].edges_vec[vec].adjPatchIndex.clear();
		for (unsigned int order_i = 0; order_i < sgraph.nodeList[sp].orderedArcs.size(); ++order_i)//修改sgraph.nodelist[i].orderedArcs的patchIndex
		{
			//方法一，找到对应的sp-->ep
			if (sgraph.nodeList[sp].orderedArcs[order_i].adjVertex == ep)
			{
				//正向修改sp-->ep
				vector<int>::iterator iter_pi = sgraph.nodeList[sp].orderedArcs[order_i].patchesIndex.begin();
				while (iter_pi != sgraph.nodeList[sp].orderedArcs[order_i].patchesIndex.end())
				{
					if (*iter_pi == dele_pchindex)
						iter_pi = sgraph.nodeList[sp].orderedArcs[order_i].patchesIndex.erase(iter_pi);
					else
						iter_pi++;
				}
				//逆向修改ep-->sp
				for (unsigned int arc_i = 0; arc_i < sgraph.nodeList[ep].orderedArcs.size(); arc_i++)
				{
					if (sgraph.nodeList[ep].orderedArcs[arc_i].adjVertex == sp)
					{
						iter_pi = sgraph.nodeList[ep].orderedArcs[arc_i].patchesIndex.begin();
						while (iter_pi != sgraph.nodeList[ep].orderedArcs[arc_i].patchesIndex.end())
						{
							if (*iter_pi == dele_pchindex)
								iter_pi = sgraph.nodeList[ep].orderedArcs[arc_i].patchesIndex.erase(iter_pi);
							else
								iter_pi++;
						}
						break; 
					}
				}
				//原方法，为了确保删除正确 用上方法
				//if (sgraph.nodeList[sp].orderedArcs[order_i].patchesIndex.size()>0)
				//{
				//	//逆向修改,从点endp开始
				//	for (unsigned int arc_i = 0; arc_i < sgraph.nodeList[ep].orderedArcs.size(); arc_i++)
				//	{
				//		if (sgraph.nodeList[ep].orderedArcs[arc_i].adjVertex == sp)
				//		{
				//			sgraph.nodeList[ep].orderedArcs[arc_i].patchesIndex.pop_back();
				//			break; 
				//		}
				//	}
				//}
				//sgraph.nodeList[sp].orderedArcs[order_i].patchesIndex.clear();//sp-ep这条edge的相邻patch删除掉，因为这个dele_petch不存在，其边上的信息也要删掉（有可能ep-sp 边上的patch是存在的）
				break;
			}
			////或者方法三，
			//vector<int>::iterator iter_pi = sgraph.nodeList[sp].orderedArcs[order_i].patchesIndex.begin();
			//while (iter_pi != sgraph.nodeList[sp].orderedArcs[order_i].patchesIndex.end())
			//{
			//	if (*iter_pi == dele_pchindex)
			//		iter_pi = sgraph.nodeList[sp].orderedArcs[order_i].patchesIndex.erase(iter_pi);
			//	else
			//		iter_pi++;
			//}			
		}
		ArcNode *anode = getArcNode(sp, ep);
		vector<int>::iterator iter_i = anode->patchesIndex.begin();
		while (iter_i != anode->patchesIndex.end())
		{
			if (*iter_i == dele_pchindex)
				iter_i = anode->patchesIndex.erase(iter_i);
			else
				iter_i++;
		}
		anode = getArcNode(ep, sp);
		iter_i = anode->patchesIndex.begin();
		while (iter_i != anode->patchesIndex.end())
		{
			if (*iter_i == dele_pchindex)
				iter_i = anode->patchesIndex.erase(iter_i);
			else
				iter_i++;
		}
	}
}
void MyGraph::deleRepeatPoint(vector<int> *plist)
{
	vector<int>::iterator iter_tp = plist->begin();
	while (iter_tp != plist->end())
	{
		vector<int>::iterator iter_tp_inn = iter_tp + 1;
		while (iter_tp_inn != plist->end())
		{
			if (*iter_tp == *iter_tp_inn)
				iter_tp_inn = plist->erase(iter_tp_inn);
			else
				iter_tp_inn++;
		}
		iter_tp++;
	}
}
int MyGraph::findInnerLine(int graph_p)//graph_p是指从第几个subgraph开始
{
	int count = 0;
	int *visitStatus = new int[sgraph.vertexNum];
	//0:unvisited,1:push in queue,2: out of queue
	for (int i = 0; i < sgraph.vertexNum; ++i)
	{
		visitStatus[i] = 0;
	}
	if (countEdgeNum(graph_p, 1) == 2)//从度>=3  /=1的点开始
	{
		vector<int> storevisit; //存  visitstatus改变了的变量
		bool fflag = false;
		stack<int> myst;//深度遍历
		myst.push(graph_p);
		storevisit.push_back(graph_p);
		visitStatus[graph_p] = 1;
		while (!myst.empty() && fflag == false)
		{
			int ver = myst.top();
			myst.pop();
			visitStatus[ver] = 2;
			ArcNode *tmp_arc = sgraph.nodeList[ver].firstEdge;
			while (tmp_arc)
			{
				if (visitStatus[tmp_arc->adjVertex] == 0)
				{
					myst.push(tmp_arc->adjVertex);
					storevisit.push_back(tmp_arc->adjVertex);
					visitStatus[tmp_arc->adjVertex] = 1;
					if (countEdgeNum(tmp_arc->adjVertex, 1) != 2)
					{
						graph_p = tmp_arc->adjVertex;
						fflag = true;
						break;
					}
				}
				tmp_arc = tmp_arc->next;
			}
		}
		for (size_t i = 0; i < storevisit.size(); ++i)
		{
			visitStatus[storevisit[i]] = 0;
		}
	}
	int tmp_ver = 0;
	stack<int> verQueue;//深度遍历
	verQueue.push(graph_p);
	visitStatus[graph_p] = 1;
	while (!verQueue.empty())
	{
		tmp_ver = verQueue.top();
		verQueue.pop();
		visitStatus[tmp_ver] = 2;
		ArcNode *tmp_arc = sgraph.nodeList[tmp_ver].firstEdge;
		while (tmp_arc)
		{
				if (visitStatus[tmp_arc->adjVertex] != 2)
				{
					if (visitStatus[tmp_arc->adjVertex] == 0)
					{
						verQueue.push(tmp_arc->adjVertex);
						visitStatus[tmp_arc->adjVertex] = 1;
					}
					if (tmp_arc->lineFlag == 0)
					{
						ArcNode *arc_reverse = getArcNode(tmp_arc->adjVertex, tmp_ver);//patch+line+patch
						if (arc_reverse == NULL)
							cout << "error in findInnerLine()";
						if (tmp_arc->patchesIndex.empty() && arc_reverse->patchesIndex.empty())//!!!说明是一条line//if (tmp_arc->patchesIndex.empty() && tt->patchesIndex.empty() && tmp_arc->innerLine == false && tt->innerLine == false)//
						{
							vector<int> tmpinnerline;
							tmpinnerline.push_back(tmp_ver);
							tmpinnerline.push_back(tmp_arc->adjVertex);
							int p1 = tmp_ver;
							int p2 = tmp_arc->adjVertex;
							int edge_count = countEdgeNum(p2, 1);
							while (edge_count == 2)
							{
								ArcNode *nextArc = sgraph.nodeList[p2].firstEdge;
								while (nextArc != NULL)
								{
									if (nextArc->lineFlag == 0)
									{
										if (nextArc->adjVertex != p1)
										{
											break;
										}
									}
									nextArc = nextArc->next;
								}
								tmpinnerline.push_back(nextArc->adjVertex);
								p1 = p2;
								p2 = nextArc->adjVertex;
								edge_count = countEdgeNum(p2, 1);
							}
							ICurve t_cv;
							t_cv.type = 0;
							t_cv.pt_vec = tmpinnerline;
							plrGraph.innerlines.push_back(t_cv);//plrGraph.innerlines.push_back(tmpinnerline);
							count++;
							for (unsigned int line_i = 1; line_i < tmpinnerline.size(); ++line_i)
							{
								getArcNode(tmpinnerline[line_i - 1], tmpinnerline[line_i])->lineFlag = 1;
								for (unsigned int i_oa = 0; i_oa < sgraph.nodeList[tmpinnerline[line_i - 1]].orderedArcs.size(); ++i_oa)
								{
									if (sgraph.nodeList[tmpinnerline[line_i - 1]].orderedArcs.at(i_oa).adjVertex == tmpinnerline[line_i])
										sgraph.nodeList[tmpinnerline[line_i - 1]].orderedArcs.at(i_oa).lineFlag = 1;
								}
								getArcNode(tmpinnerline[line_i], tmpinnerline[line_i - 1])->lineFlag = 1;
								for (unsigned int i_oa = 0; i_oa < sgraph.nodeList[tmpinnerline[line_i]].orderedArcs.size(); ++i_oa)
								{
									if (sgraph.nodeList[tmpinnerline[line_i]].orderedArcs.at(i_oa).adjVertex == tmpinnerline[line_i - 1])
										sgraph.nodeList[tmpinnerline[line_i]].orderedArcs.at(i_oa).lineFlag = 1;
								}
							}
						}
					}		
				}
			
			tmp_arc = tmp_arc->next;
		}
	}
	delete[]visitStatus;
	visitStatus = NULL;
	return count;
}

Patch *MyGraph::getAPatch(int patchindex)
{
	Patch *pp = &plrGraph.patches[patchindex];
	if (pp->patch_index == patchindex)
		return pp;
	else
		return NULL;
}
int MyGraph::deleRepeatPatch()
{
	//delete repeat patches删除patch，要维护好两个数据结构：plrGraph.patch 和 sgraph.nodelist[i].orderedArcs / firstEdge
	int delecount = 0;
	vector<int> vec1, vec2, vec1_add;
	vector<Patch>::iterator iter = plrGraph.patches.begin();
	while (iter != plrGraph.patches.end())
	{
		if (iter->type != -1)
		{
			vec1 = iter->pointIndexOfPatch;
			vector<Patch>::iterator iter_inn = iter + 1;
			while (iter_inn != plrGraph.patches.end())
			{
				vec2 = iter_inn->pointIndexOfPatch;
				if (vec1.size() == vec2.size() && iter_inn->type != -1)
				{
					bool dflag = isDifferentPath(vec1, vec2);
					if (dflag == false)//false-> vec1=vec2
					{
						delecount++;
						iter_inn->type = -1;//-1  表示patch不可用
						for (unsigned int vec = 0; vec < iter->edges_vec.size(); vec++)//iter指向要保留的那个patch，遍历该patch的所有edge，
						{
							int sp = iter->edges_vec[vec].startp_index;
							int ep = iter->edges_vec[vec].endp_index;
							if (iter->edges_vec[vec].adjPatchIndex.size()>1)//该判断可以删去
							{
								iter->edges_vec[vec].adjPatchIndex.pop_back();// 修改plrGraph.patch中的edge的patchIndex            getBorderofPatch(iter_inn->edges_vec[vec].adjPatchIndex[1], sp, ep)->adjPatchIndex.pop_back();
							}
							for (unsigned int order_i = 0; order_i < sgraph.nodeList[sp].orderedArcs.size(); ++order_i)//修改sgraph.nodelist[i].orderedArcs的patchIndex
							{
								if (sgraph.nodeList[sp].orderedArcs[order_i].adjVertex == ep)
								{
									sgraph.nodeList[sp].orderedArcs[order_i].patchesIndex.pop_back();
									break;
								}
							}
							ArcNode *anode = getArcNode(sp, ep);//修改sgraph.nodelist[i]。firstEdge...的patchIndex
							if (anode->patchesIndex.size() > 1)  //edge是按顺序排放的,边sp-ep相邻的patchIndex===iter->edges_vec[vec].adjPatchIndex					
								anode->patchesIndex.pop_back();
						}
						for (unsigned int vec = 0; vec < iter_inn->edges_vec.size(); vec++)//iter_inn指向要删除的那个patch,删除掉每条边相邻的patch，但是这个patch的edge，point还都有
						{
							int sp = iter_inn->edges_vec[vec].startp_index;
							int ep = iter_inn->edges_vec[vec].endp_index;
							iter_inn->edges_vec[vec].adjPatchIndex.clear();
							for (unsigned int order_i = 0; order_i < sgraph.nodeList[sp].orderedArcs.size(); ++order_i)//修改sgraph.nodelist[i].orderedArcs的patchIndex
							{
								if (sgraph.nodeList[sp].orderedArcs[order_i].adjVertex == ep)
								{
									sgraph.nodeList[sp].orderedArcs[order_i].patchesIndex.clear();
									break;
								}
							}
							ArcNode *anode = getArcNode(sp, ep);
							anode->patchesIndex.clear();
						}
						iter_inn++;//break;//相同的patch只有一对儿，找到了就break掉，而不是 iter_inn++;(iter_inn++;保险一点)
					}
					else
						iter_inn++;
				}
				else
					iter_inn++;
			}
		}
		iter++;
	}
	return delecount;
}
bool MyGraph::isDifferentPath(vector<int> p1, vector<int> p2)
{
	vector<int> vec1, vec2, vec1_add;
	vec1 = p1;
	vec1_add = p1;
	int k = 1;
	while (vec1[k] != vec1[0])//eg vec1=1,2,3,4,1  vec1_add=1,2,3,4,1 + 2,3,4 =1,2,3,4,1,2,3,4
	{
		vec1_add.push_back(vec1[k]);
		k = k + 1;
	}
	bool sflag = false;
	bool findflag = false;
	vec2 = p2;
	for (unsigned int h = 0; h < vec1_add.size(); ++h)
	{
		if (vec1_add[h] == vec2[0])
		{
			findflag = true;
			for (unsigned int u = 0; u < vec2.size(); ++u)
			{
				if (vec2[u] != vec1_add[h + u])//vec2!=vec1
				{
					sflag = true;
					break;
				}
			}
			break;
		}
	}
	if (findflag == false)//第一个相等的数字都没有找到，vec2!=vec1
	{
		sflag = true;
	}
	else//倒序情况
	{
		if (sflag == true)
		{
			sflag = false;
			reverse(vec2.begin(), vec2.end());
			for (unsigned int h = 0; h < vec1_add.size(); ++h)
			{
				if (vec1_add[h] == vec2[0])
				{
					for (unsigned int u = 0; u < vec2.size(); ++u)
					{
						if (vec2[u] != vec1_add[h + u])
						{
							sflag = true;
							break;
						}
					}
					break;
				}
			}
		}
	}
	return sflag;
}
void MyGraph::printpathOfPatches()
{
	cout << "path Of Patches:" << plrGraph.patches.size() << endl;
	for (unsigned int i = 0; i < plrGraph.patches.size(); ++i)
	{
		if (plrGraph.patches[i].type >= 0)
		{
			for (unsigned int j = 0; j < plrGraph.patches[i].pointIndexOfPatch.size(); ++j)
			{
				cout << plrGraph.patches[i].pointIndexOfPatch.at(j) << ",";
			}
			cout << endl;
		}
	}
}
/*findPatch()之后执行该函数，（1）把isolatedline和attachline分开，(2)otherline=attach+isolated （3）对patch中的vector<pair<int,int>>point_line赋值*/
void MyGraph::divideAttachLines()
{
	/*1、（在找完patch之后执行该函数）把isolatedline和attachline分开(在findline中统一标记为2)。*/
	int* visitStatus = new int[sgraph.nodeList.size()];
	for (int i = 0; i < sgraph.vertexNum; ++i)
	{
		visitStatus[i] = 0;
	}
	vector<ICurve> tmp_outlines;
	vector<ICurve>::iterator iter_l = plrGraph.attachlines.begin();//vector<vector<int>>::iterator iter_l = plrGraph.outlines.begin();
	while (iter_l != plrGraph.attachlines.end())
	{
		vector<int> tmpline = iter_l->pt_vec;
		int p1 = tmpline.front();
		int p2 = tmpline[1];
		int p3 = tmpline.back();
		int ct1 = countEdgeNum(p1, 0);
		int ct2 = countEdgeNum(p3, 0);
		if (visitStatus[p1] != 0)//说明这条线已经trace过，
		{
			iter_l++;
			continue;
		}
		if (ct1 == 1 && ct2 == 1)//（1.isolated line 3的情况)
		{
			bool iso_flag = true;
			for (unsigned int j = 1; j < tmpline.size() - 1; ++j)//check
			{
				int ct_tmp = countEdgeNum(tmpline[j], 0);
				if (ct_tmp != 2)
				{
					iso_flag = false;
					break;
				}
			}
			if (iso_flag == true)
			{
				for (unsigned int j = 1; j < tmpline.size(); ++j)
				{
					getArcNode(tmpline[j - 1], tmpline[j])->lineFlag = 3;
					getArcNode(tmpline[j], tmpline[j - 1])->lineFlag = 3;
					visitStatus[tmpline[j - 1]] = 2;
				}
				visitStatus[tmpline.back()] = 2;
				plrGraph.isolatedlines.push_back(*iter_l);//！！！
				iter_l = plrGraph.attachlines.erase(iter_l);
				//iter_l++;
				continue;
			}//若中间有点的度>2,则执行else中的操作
		}
		//else//outline 的情况，重新trace，存入tmp_outlines
		//{
		//（2.剩下的是attachline的情况，attachline的第一个点（p1）可能是度为1的 末端点，也有可能是度>2的点
		ArcNode* i_arc = sgraph.nodeList[p1].firstEdge;
		while (i_arc != NULL)
		{
			if (i_arc->adjVertex == p2)
			{
				break;
			}
			i_arc = i_arc->next;
		}
		char arc_ctype = i_arc->curveType;//找到该条line所对应的原curve的index和type
		int arc_cindex = i_arc->curveIndex;
		char arc_otype; int arc_oindex; vector<mypoint2f> arc_plist;
		getOrigiCurveList(arc_ctype, arc_cindex, arc_otype, arc_oindex, arc_plist, 0);
		vector<int> a_line;
		a_line.push_back(p1);
		vector<int> record_degree;//记录点的度数>2的点，
		stack<int> treestack;
		treestack.push(p1);
		visitStatus[p1] = 1;
		while (!treestack.empty())
		{
			int verindex = treestack.top();
			treestack.pop();
			visitStatus[verindex] = 2;
			ArcNode* tmp_arc = sgraph.nodeList[verindex].firstEdge;
			int neighbor_ct = countEdgeNum(verindex, 0);
			while (tmp_arc != NULL)
			{
				if (visitStatus[tmp_arc->adjVertex] != 2)
				{
					if (tmp_arc->lineFlag == 2)
					{
						if (neighbor_ct > 2)
						{
							char tmp_otype; int tmp_oindex; vector<mypoint2f> tmparc_plist;
							getOrigiCurveList(tmp_arc->curveType, tmp_arc->curveIndex, tmp_otype, tmp_oindex, tmparc_plist, 0);
							//if (tmp_arc->curveIndex == arc_oindex && tmp_arc->curveType == arc_otype && tmp_arc->patchesIndex.empty())
							if (tmp_otype == arc_otype && tmp_oindex == arc_oindex &&tmp_arc->patchesIndex.empty())
							{
								if (countEdgeNum(tmp_arc->adjVertex, 0)>2)
									record_degree.push_back(tmp_arc->adjVertex);
								treestack.push(tmp_arc->adjVertex);
								visitStatus[tmp_arc->adjVertex] = 1;
								a_line.push_back(tmp_arc->adjVertex);
								//getOrigiCurveList(tmp_arc->curveType, tmp_arc->curveIndex, arc_otype, arc_oindex, arc_plist, 0);//找到该条line所对应的原curve的index和type
								arc_otype = tmp_otype;
								arc_oindex = tmp_oindex;
								break;
							}
						}
						else
						{
							int ttt_c = countEdgeNum(tmp_arc->adjVertex, 0);
							if (ttt_c > 2)
							{
								record_degree.push_back(tmp_arc->adjVertex);
								getOrigiCurveList(tmp_arc->curveType, tmp_arc->curveIndex, arc_otype, arc_oindex, arc_plist, 0);//找到该条line所对应的原curve的index和type
							}
							treestack.push(tmp_arc->adjVertex);
							visitStatus[tmp_arc->adjVertex] = 1;
							a_line.push_back(tmp_arc->adjVertex);
							break;
						}
					}
				}
				tmp_arc = tmp_arc->next;
			}
		}
		//(3.attachline 分段)
		vector<int> seg_pindex;
		for (unsigned int j = 1; j < a_line.size()-1; ++j)
		{
			ArcNode* anode = sgraph.nodeList[a_line[j]].firstEdge;//第一种情况 line被patch 的边截成两半
			bool bflag = false;
			while (anode != NULL)
			{
				if (!anode->patchesIndex.empty())
				{
					bflag = true;
					seg_pindex.push_back(j);
					break;
				}
				anode = anode->next;
			}
			if (bflag == false)//第二种 拐角<150
			{
				mypoint2f p1 = sgraph.nodeList[a_line[j - 1]].position;
				mypoint2f p2 = sgraph.nodeList[a_line[j]].position;
				mypoint2f p3 = sgraph.nodeList[a_line[j + 1]].position;
				vector<mypoint2f> p2p1_plist;
				vector<mypoint2f> p2p3_plist;
				getPListOfaCurve(p2, p1, a_line[j], a_line[j - 1], &p2p1_plist, 0);
				getPListOfaCurve(p2, p3, a_line[j], a_line[j + 1], &p2p3_plist, 0);

				double costheta = getCosineAngle(p2p1_plist[0], p2p1_plist[1], p2p3_plist[0], p2p3_plist[1]);
				double agl = acos(costheta) * 180 / 3.1415926;//转化成角度 //
				if (agl < 150)
				{
					seg_pindex.push_back(j);
				}
			}		
		}
		if (seg_pindex.empty())
		{
			ICurve icurve;
			icurve.type = 0;
			icurve.pt_vec = a_line;
			tmp_outlines.push_back(icurve);//tmp_outlines.push_back(a_line);
		}
		else
		{
			seg_pindex.insert(seg_pindex.begin(),0);
			seg_pindex.push_back(a_line.size()-1);
			for (unsigned int seg_i = 1; seg_i < seg_pindex.size(); ++seg_i)
			{
				int bb = seg_pindex[seg_i - 1];
				int ee = seg_pindex[seg_i];
				vector<int> a_seg;
				for (int si = bb; si <= ee; ++si)
				{
					a_seg.push_back(a_line[si]);
				}
				ICurve icurve;
				icurve.type = 0;
				icurve.pt_vec = a_seg;
				tmp_outlines.push_back(icurve);//tmp_outlines.push_back(a_line);
			}
		}
		/*原来未分割，直接存入tmp_outline*/
		//ICurve icurve;
		//icurve.type = 0;
		//icurve.pt_vec = a_line;
		//tmp_outlines.push_back(icurve);//tmp_outlines.push_back(a_line);
		/*度>=3的点 visit重新设为0*/
		for (unsigned int j = 0; j < record_degree.size(); ++j)
		{
			visitStatus[record_degree[j]] = 0;
		}
		iter_l++;
		//}
	}

	plrGraph.attachlines.clear(); plrGraph.attachlines.swap(vector<ICurve>());
	plrGraph.attachlines.insert(plrGraph.attachlines.begin(), tmp_outlines.begin(), tmp_outlines.end());
	delete[]visitStatus;
	visitStatus = NULL;
	/*2、otherlines = outline+isolateline，用于mixtopology中查找line，或者patch结构体中point_line中的line的index*/
	if(!plrGraph.attachlines.empty())
		plrGraph.otherlines.insert(plrGraph.otherlines.end(), plrGraph.attachlines.begin(), plrGraph.attachlines.end());
	if (!plrGraph.isolatedlines.empty())
		plrGraph.otherlines.insert(plrGraph.otherlines.end(), plrGraph.isolatedlines.begin(), plrGraph.isolatedlines.end());
	/*3、找到attachline对应的那个patch，添加vector<pair<int,int>> point_line */
	for (unsigned int i = 0; i < plrGraph.attachlines.size(); ++i)
	{
		int p_front = plrGraph.attachlines[i].pt_vec.front();
		int p_back = plrGraph.attachlines[i].pt_vec.back();
		//先查看front点是否连接有patch，
		vector<int> pt_patch;
		ArcNode* anode = sgraph.nodeList[p_front].firstEdge;
		while (anode != NULL)
		{
			if (!anode->patchesIndex.empty())//说明这条arc有patch，或者用anode->lineFlag == 0来判断也可以
			{
				if (anode->lineFlag != 0)
					cout << "error" << endl;//说明有错误，其中一个条件满足即可，这里是为了检测
				for (unsigned int j = 0; j < anode->patchesIndex.size(); ++j)
				{
					pt_patch.push_back(anode->patchesIndex[j]);
				}
			}
			anode = anode->next;
		}
		if (!pt_patch.empty())
		{
			deleRepeatPoint(&pt_patch);
			for (unsigned int j = 0; j < pt_patch.size(); ++j)
			{
				int pch = pt_patch[j];   //patch index
				plrGraph.patches[pch].point_line.push_back(make_pair(p_front, i));//在pch中的p_front点连接有线i (i是plrGraph.attachlines的下标！！)
			}
		}
		else
		{
			//front点没有连接的patch，再查看back点是否有patch
			ArcNode* anode_back = sgraph.nodeList[p_back].firstEdge;
			while (anode_back != NULL)
			{
				if (!anode_back->patchesIndex.empty())//说明这条arc有patch，或者用anode_back->lineFlag == 0来判断也可以
				{
					if (anode_back->lineFlag != 0)
						cout << "error" << endl;//说明有错误，其中一个条件满足即可，这里是为了检测
					for (unsigned int j = 0; j < anode_back->patchesIndex.size(); ++j)
					{
						pt_patch.push_back(anode_back->patchesIndex[j]);
					}
				}
				anode_back = anode_back->next;
			}
			if (!pt_patch.empty())
			{
				deleRepeatPoint(&pt_patch);
				for (unsigned int j = 0; j < pt_patch.size(); ++j)
				{
					int pch = pt_patch[j];   //patch index
					plrGraph.patches[pch].point_line.push_back(make_pair(p_back, i));//在pch中的p_front点连接有线i (i是plrGraph.attachlines的下标！！)
				}
			}
			//else//说明这条attachline 没有patch相邻
			//{
			//	cout << "!!" << endl;
			//}
		}
		
	}
}


void MyGraph::jointPolylines()
{
	int verindex1 = 0;
	int verindex2 = 0;
	
	for (unsigned int i = 0; i < plrGraph.innerlines.size(); ++i)
	{
		if (!plrGraph.innerlines[i].pt_vec.empty() && plrGraph.innerlines[i].type >= 0)
		{
			PolyLine apl;//vector<mypoint2f> pvec;
			for (unsigned int j = 1; j < plrGraph.innerlines[i].pt_vec.size(); ++j)
			{
				verindex1 = plrGraph.innerlines[i].pt_vec.at(j - 1);
				verindex2 = plrGraph.innerlines[i].pt_vec.at(j);
				mypoint2f point1 = sgraph.nodeList[verindex1].position;//=lineWithCommonPoints[verindex1].first
				mypoint2f point2 = sgraph.nodeList[verindex2].position;
				if (!apl.pointlist.empty())
					apl.pointlist.pop_back();
				getPListOfaCurve(point1, point2, verindex1, verindex2, &apl.pointlist,0);
				/*mycircle acir;
				acir.cx = point1.x; acir.cy = point1.y; acir.radius = 0.5;
				circle_vec.push_back(acir);
				acir.cx = point2.x; acir.cy = point2.y; acir.radius = 0.5;
				circle_vec.push_back(acir);*/
			}
			pointListOfLines.push_back(apl);//plrGraph.patches.push_back(apatch);
		}
	}
	cout << "line size after merging:" << pointListOfLines.size() << endl;

	for (unsigned int i = 0; i < plrGraph.otherlines.size(); ++i)//plrGraph.outlines.size()
	{
		if (!plrGraph.otherlines[i].pt_vec.empty() && plrGraph.otherlines[i].type >= 0)
		{
			PolyLine apl;//vector<mypoint2f> pvec;
			for (unsigned int j = 1; j < plrGraph.otherlines[i].pt_vec.size(); ++j)
			{
				verindex1 = plrGraph.otherlines[i].pt_vec.at(j - 1);
				verindex2 = plrGraph.otherlines[i].pt_vec.at(j);
				mypoint2f point1 = sgraph.nodeList[verindex1].position;//=lineWithCommonPoints[verindex1].first
				mypoint2f point2 = sgraph.nodeList[verindex2].position;
				if (!apl.pointlist.empty())
					apl.pointlist.pop_back();
				getPListOfaCurve(point1, point2, verindex1, verindex2, &apl.pointlist, 0);
			}
			pointListOfLines.push_back(apl);//plrGraph.patches.push_back(apatch);
		}
	}
	cout << "line size after merging:" << pointListOfLines.size() << endl;


	//for (unsigned int i = 0; i < plrGraph.attachlines.size(); ++i)//plrGraph.outlines.size()
	//{
	//	if (!plrGraph.attachlines[i].pt_vec.empty() && plrGraph.attachlines[i].type >= 0)
	//	{
	//		PolyLine apl;//vector<mypoint2f> pvec;
	//		for (unsigned int j = 1; j < plrGraph.attachlines[i].pt_vec.size(); ++j)
	//		{
	//			verindex1 = plrGraph.attachlines[i].pt_vec.at(j - 1);
	//			verindex2 = plrGraph.attachlines[i].pt_vec.at(j);
	//			mypoint2f point1 = sgraph.nodeList[verindex1].position;//=lineWithCommonPoints[verindex1].first
	//			mypoint2f point2 = sgraph.nodeList[verindex2].position;
	//			if (!apl.pointlist.empty())
	//				apl.pointlist.pop_back();
	//			getPListOfaCurve(point1, point2, verindex1, verindex2, &apl.pointlist, 0);
	//		}
	//		pointListOfLines.push_back(apl);//plrGraph.patches.push_back(apatch);
	//	}
	//}
	//for (unsigned int i = 0; i < plrGraph.isolatedlines.size(); ++i)
	//{
	//	if (!plrGraph.isolatedlines[i].pt_vec.empty() && plrGraph.isolatedlines[i].type >= 0)
	//	{
	//		PolyLine apl;//vector<mypoint2f> pvec;
	//		for (unsigned int j = 1; j < plrGraph.isolatedlines[i].pt_vec.size(); ++j)
	//		{
	//			verindex1 = plrGraph.isolatedlines[i].pt_vec.at(j - 1);
	//			verindex2 = plrGraph.isolatedlines[i].pt_vec.at(j);
	//			mypoint2f point1 = sgraph.nodeList[verindex1].position;//=lineWithCommonPoints[verindex1].first
	//			mypoint2f point2 = sgraph.nodeList[verindex2].position;
	//			if (!apl.pointlist.empty())
	//				apl.pointlist.pop_back();
	//			getPListOfaCurve(point1, point2, verindex1, verindex2, &apl.pointlist,0);
	//		}
	//		pointListOfLines.push_back(apl);//plrGraph.patches.push_back(apatch);
	//	}
	//}
}
void MyGraph::jointCurveTPatch()
{
	for (unsigned int i = 0; i < plrGraph.patches.size(); ++i)
	{		
		if (!plrGraph.patches[i].pointIndexOfPatch.empty() && plrGraph.patches[i].type >= 0)//不能=-1
		{
			vector<mypoint2f> apatch;
			jointMethod(&plrGraph.patches[i].pointIndexOfPatch, &apatch, 0);/*0是density 以下注释是已经封装成这个函数*/
			//for (unsigned int j = 1; j < plrGraph.patches[i].pointIndexOfPatch.size(); ++j)
			//{
			//	verindex1 = plrGraph.patches[i].pointIndexOfPatch.at(j - 1);
			//	verindex2 = plrGraph.patches[i].pointIndexOfPatch.at(j);
			//	mypoint2f point1 = sgraph.nodeList[verindex1].position;//=lineWithCommonPoints[verindex1].first
			//	mypoint2f point2 = sgraph.nodeList[verindex2].position;
			//	if (!apatch.empty())
			//		apatch.pop_back();
			//	getPListOfaCurve(point1, point2, verindex1, verindex2, &apatch,0);
			//	////添加组成patch的节点， 
			//		//	mycircle acir;
			//		//	acir.cx = point1.x; acir.cy = point1.y; acir.radius = 0.5;
			//		//	circle_vec.push_back(acir);
			//}
			if (isPointRoughEqual(apatch.front(), apatch.back()))
				apatch.pop_back();
			pointListOfPatches.push_back(make_pair(plrGraph.patches[i].type, apatch));//plrGraph.patches.push_back(apatch);	
		}
	}
	cout << "patch size after merging:" << pointListOfPatches.size() << endl;
	
	//ofstream outfile("E:/out.txt");
	//if (outfile.is_open())
	//{
	//	for (unsigned int i = 0; i < pointListOfPatches.size(); ++i)
	//	{
	//		vector<mypoint2f> tmpvec = pointListOfPatches[i];
	//		mypoint2f spoint = tmpvec.front();
	//		tmpvec.insert(tmpvec.begin(), tmpvec.back());
	//		tmpvec.push_back(spoint);
	//		for (unsigned int j = 2; j < tmpvec.size(); ++j)
	//		{
	//			mypoint2f v1 = tmpvec[j] - tmpvec[j - 1];
	//			mypoint2f v2 = tmpvec[j - 1] - tmpvec[j - 2];
	//			double theta = dot(v1, v2) / (LineSegmentLength(v1, mypoint2f(0, 0)) * LineSegmentLength(v2, mypoint2f(0, 0)));
	//			double degree = acos(theta);
	//			outfile << degree;
	//			outfile << " ";
	//		}
	//		outfile << "\n";
	//	}
	//}
	//outfile.close();
	////排序，patch周长从大到小
	//for (unsigned int i = 0; i < pointListOfPatches.size(); ++i)
	//{
	//	double len1 = getPatchPerimeter(pointListOfPatches[i]);
	//	for (unsigned int j = i ; j < plrGraph.patches.size(); ++j)
	//	{
	//		double len2 = getPatchPerimeter(pointListOfPatches[j]);
	//		if (len1 < len2)
	//		{
	//			vector<mypoint2f> tmp_pvec = pointListOfPatches[i];
	//			pointListOfPatches[i].swap(vector<mypoint2f>());
	//			pointListOfPatches[i] = pointListOfPatches[j];
	//			pointListOfPatches[j].swap(vector<mypoint2f>());
	//			pointListOfPatches[j] = tmp_pvec;
	//		}
	//	}
	//}
}
void MyGraph::jointMethod(vector<int>* pindex_vec, vector<mypoint2f>* plist,int density)
{
	for (unsigned int j = 1; j < pindex_vec->size(); ++j)
	{
		int verindex1 = pindex_vec->at(j - 1);
		int verindex2 = pindex_vec->at(j);
		mypoint2f point1 = sgraph.nodeList[verindex1].position;//=lineWithCommonPoints[verindex1].first
		mypoint2f point2 = sgraph.nodeList[verindex2].position;
		if (!plist->empty())
			plist->pop_back();
		getPListOfaCurve(point1, point2, verindex1, verindex2, plist, density);
		////添加组成patch的节点， 
		//	mycircle acir;
		//	acir.cx = point1.x; acir.cy = point1.y; acir.radius = 0.5;
		//	circle_vec.push_back(acir);
	}
}
void MyGraph::joinCurveToSpecialPatch(vector<int> pch)//为了测试观察某几个patch用
{
	int verindex1 = 0;
	int verindex2 = 0;
	for (unsigned int i = 0; i < pch.size(); ++i)
	{		
		vector<mypoint2f> apatch;
		for (unsigned int j = 1; j < plrGraph.patches[pch[i]].pointIndexOfPatch.size(); ++j)
		{
			verindex1 = plrGraph.patches[pch[i]].pointIndexOfPatch.at(j - 1);
			verindex2 = plrGraph.patches[pch[i]].pointIndexOfPatch.at(j);
			mypoint2f point1 = sgraph.nodeList[verindex1].position;//=lineWithCommonPoints[verindex1].first
			mypoint2f point2 = sgraph.nodeList[verindex2].position;
			if (!apatch.empty())
				apatch.pop_back();
			getPListOfaCurve(point1, point2, verindex1, verindex2, &apatch,0);
			////添加组成patch的节点， 
			//	mycircle acir;
			//	acir.cx = point1.x; acir.cy = point1.y; acir.radius = 0.5;
			//	circle_vec.push_back(acir);
		}
		if (isPointRoughEqual(apatch.front(), apatch.back()))
			apatch.pop_back();
		pointListOfPatches.push_back(make_pair(plrGraph.patches[pch[i]].type, apatch));//plrGraph.patches.push_back(apatch);	
	}
}
void MyGraph::getPListOfaCurve(mypoint2f point1, mypoint2f point2, int verindex1, int verindex2, vector<mypoint2f> *vec, int density)
{
	//density 0表示最密集，1 表示隔一个点取一个，与MySvgFile.cpp的splitACurve()不一样，splitACurve()中0表示每条取50个点，1表示平均取点比较密集
	char ctype = ' ';
	int cindex = 0;
	int order = -1;//0:positive order  1:reversed order
	if (verindex2 == -1)//说明 verindex2不知道，只能通过point2 --> pointTopl
	{
		vector<TerminalPoint>::iterator iter = pointTopl[verindex1].second.begin();
		for (iter; iter != pointTopl[verindex1].second.end(); ++iter)
		{
			if (isPointRoughEqual(iter->endpoint, point2))//if (LineSegmentLength(iter->endpoint, point2) < eps) 
			{
				if (iter->curvetype == 'C')
				{
					ctype = 'C';
					cindex = iter->curveindex;
					if (isPointRoughEqual(cubicBezier_vec[iter->curveindex].p0, point2))//if (LineSegmentLength(cubicBezier_vec[iter->curveindex].p0, point2)<eps)
					{			
						order = 1;
					}
					else 
					{
						order = 0;
					}
				}
				else if (iter->curvetype == 'Q')
				{
					ctype = 'Q';
					cindex = iter->curveindex;
					if (isPointRoughEqual(quadraticBezier_vec[iter->curveindex].p0, point2))//if (LineSegmentLength(quadraticBezier_vec[iter->curveindex].p0, point2)<eps)
					{
						order = 1;
					}
					else
					{
						order = 0;
					}
				}
				else if (iter->curvetype == 'L')
				{
					ctype = 'L';
					cindex = iter->curveindex;
					if (isPointRoughEqual(polyline_vec[iter->curveindex].pointlist.front(), point2))//if (LineSegmentLength(polyline_vec[iter->curveindex].front(), point2)<eps)
					{		
						order = 1;
					}
					else
					{
						order = 0;
					}
				}
				break;
			}
		}
		if (ctype == ' ')
			cout << "?" << endl;
	}
	else
	{
		ArcNode *edgenode = sgraph.nodeList[verindex1].firstEdge;
		while (edgenode)
		{
			if (edgenode->adjVertex == verindex2)
			{
				ctype = edgenode->curveType;
				cindex = edgenode->curveIndex;
				if (edgenode->curveType == 'C')
				{
					if (isPointRoughEqual(cubicBezier_vec[edgenode->curveIndex].p0, point2))
						order = 1;
					else
						order = 0;
				}
				else if (edgenode->curveType == 'Q')
				{
					if (isPointRoughEqual(quadraticBezier_vec[edgenode->curveIndex].p0, point2))
						order = 1;
					else
						order = 0;
				}
				else if (edgenode->curveType == 'L')
				{
					if (isPointRoughEqual(polyline_vec[edgenode->curveIndex].pointlist.front(), point2))
						order = 1;
					else
						order = 0;
				}
				break;
			}
			edgenode = edgenode->next;
		}
		if (order == -1)
		{
			cout << "?? getPListOfaCurve() error" << endl;
			return;
		}
	}
	switch (ctype)
	{
		case 'C':
		{
			getPListfromCubic(order, cubicBezier_vec[cindex], density, vec);//getPListfromCubic(order, cindex, vec);
			break;
		}
		case 'Q':
		{
			getPListfromQuadratic(order, quadraticBezier_vec[cindex], density, vec);
			break;
		}
		case 'L':
		{
			getPListfromPolyline(order, polyline_vec[cindex],density, vec);
			break;
		}
	}
}
void MyGraph::getPListfromCubic(int order, cubicBezier aCBcurve, int density,vector<mypoint2f> *vec)
{
	GLfloat p1[2], p2[2], p3[2], p4[2];
	if (order == 0) 
	{
		p1[0] = aCBcurve.p0.x; p1[1] = aCBcurve.p0.y;
		p2[0] = aCBcurve.p1.x; p2[1] = aCBcurve.p1.y;
		p3[0] = aCBcurve.p2.x; p3[1] = aCBcurve.p2.y;
		p4[0] = aCBcurve.p3.x; p4[1] = aCBcurve.p3.y;
	}
	else
	{
		p1[0] = aCBcurve.p3.x; p1[1] = aCBcurve.p3.y;
		p2[0] = aCBcurve.p2.x; p2[1] = aCBcurve.p2.y;
		p3[0] = aCBcurve.p1.x; p3[1] = aCBcurve.p1.y;
		p4[0] = aCBcurve.p0.x; p4[1] = aCBcurve.p0.y;
	}	
	double totallong = 0;
	totallong = LineSegmentLength(aCBcurve.p0, aCBcurve.p1) + LineSegmentLength(aCBcurve.p1, aCBcurve.p2) + LineSegmentLength(aCBcurve.p2, aCBcurve.p3);
	int diviNum = totallong / 0.8;//0.8
	double step = 1 / (double)diviNum;
	if (density == 0)///密集
	{
		//三次公式 B(t)=p0*(1-t)^3 + p1*3t(1-t)^2 + p2*3t^2(1-t) + p3*t^3
		step = step*(density + 1);
		if (step >= 1)
			cout << endl;
		mypoint2f apoint(0, 0);
		for (double t = 0; t < 1; t += step)
		{
			double a1 = pow((1 - t), 3);
			double a2 = pow((1 - t), 2) * 3 * t;
			double a3 = 3 * t*t*(1 - t);
			double a4 = t*t*t;
			apoint.x = a1*p1[0] + a2*p2[0] + a3*p3[0] + a4*p4[0];
			apoint.y = a1*p1[1] + a2*p2[1] + a3*p3[1] + a4*p4[1];
			vec->push_back(apoint);
		}
		apoint = mypoint2f(p4[0], p4[1]);
		if (!isPointRoughEqual(vec->back(),apoint))
			vec->push_back(apoint);
	}
	else //if (density == 1)
	{
		//GLint i = 0;
		//GLfloat pcurve[51][2];
		////GLfloat step = 1 / 50;//t[i]均匀的
		////三次公式 B(t)=p0*(1-t)^3 + p1*3t(1-t)^2 + p2*3t^2(1-t) + p3*t^3
		//for (GLfloat t = 0.0; t <= 1.0; t += 0.02)//t!=1, do not put the last point into vector,because the next curve's first point is its last point
		//{
		//	GLfloat a1 = pow((1 - t), 3);
		//	GLfloat a2 = pow((1 - t), 2) * 3 * t;
		//	GLfloat a3 = 3 * t*t*(1 - t);
		//	GLfloat a4 = t*t*t;
		//	pcurve[i][0] = a1*p1[0] + a2*p2[0] + a3*p3[0] + a4*p4[0];
		//	pcurve[i][1] = a1*p1[1] + a2*p2[1] + a3*p3[1] + a4*p4[1];
		//	i = i + 1;
		//}
		//for (int j = 0; j < 51; j++)
		//{
		//	mypoint2f tp(pcurve[j][0], pcurve[j][1]);
		//	vec->push_back(tp);
		//}
		step = step*(density + 1);
		if (step >=0.25)
			step = 0.02;
		mypoint2f apoint(0, 0);
		for (double t = 0; t < 1; t += step)
		{
			double a1 = pow((1 - t), 3);
			double a2 = pow((1 - t), 2) * 3 * t;
			double a3 = 3 * t*t*(1 - t);
			double a4 = t*t*t;
			apoint.x = a1*p1[0] + a2*p2[0] + a3*p3[0] + a4*p4[0];
			apoint.y = a1*p1[1] + a2*p2[1] + a3*p3[1] + a4*p4[1];
			vec->push_back(apoint);
		}
		apoint = mypoint2f(p4[0], p4[1]);
		if (vec->back() != apoint)
			vec->push_back(apoint);
	}
}
void MyGraph::getPListfromQuadratic(int order, quaBezier aQBcurve, int density, vector<mypoint2f> *vec)
{
	GLfloat p1[2], p2[2], p3[2];
	if (order == 0)
	{
		p1[0] = aQBcurve.p0.x; p1[1] = aQBcurve.p0.y;
		p2[0] = aQBcurve.p1.x; p2[1] = aQBcurve.p1.y;
		p3[0] = aQBcurve.p2.x; p3[1] = aQBcurve.p2.y;
	}
	else
	{
		p1[0] = aQBcurve.p2.x; p1[1] = aQBcurve.p2.y;
		p2[0] = aQBcurve.p1.x; p2[1] = aQBcurve.p1.y;
		p3[0] = aQBcurve.p0.x; p3[1] = aQBcurve.p0.y;
	}
	double totallong = 0;
	totallong = LineSegmentLength(aQBcurve.p0, aQBcurve.p1) + LineSegmentLength(aQBcurve.p1, aQBcurve.p2);
	int diviNum = totallong / 0.8;//0.8
	double step = 1 / (double)diviNum;
	if (density == 0)//密集
	{
		step = step*(density + 1);
		mypoint2f apoint(0, 0);
		for (double t = 0; t < 1; t += step)
		{
			double a1 = pow((1 - t), 2);
			double a2 = 2 * t*(1 - t);
			double a3 = t*t;
			apoint.x = a1*p1[0] + a2*p2[0] + a3*p3[0];
			apoint.y = a1*p1[1] + a2*p2[1] + a3*p3[1];
			vec->push_back(apoint);
		}
		apoint = mypoint2f(p3[0], p3[1]);
		if (vec->back() != apoint)
			vec->push_back(apoint);
	}
	else// if (density==1)
	{
		//GLfloat pcurve[51][2];//GLfloat pcurve[51][2];
		//GLint i = 0;
		//for (GLfloat t = 0.0; t <= 1.0; t += 0.02)//for (GLfloat t = 0.0; t <=1.0; t += 0.02)
		//{
		//	GLfloat a1 = pow((1 - t), 2);
		//	GLfloat a2 = 2 * t*(1 - t);
		//	GLfloat a3 = t*t;
		//	pcurve[i][0] = a1*p1[0] + a2*p2[0] + a3*p3[0];
		//	pcurve[i][1] = a1*p1[1] + a2*p2[1] + a3*p3[1];
		//	i = i + 1;
		//}
		//for (int j = 0; j < 51; j++)//for (int i = 0; i < 51; i++)
		//{
		//	mypoint2f tp(pcurve[j][0], pcurve[j][1]);
		//	vec->push_back(tp);
		//}
		step = step*(density + 1);
		if (step >= 0.25)
			step = 0.02;
		mypoint2f apoint(0, 0);
		for (double t = 0; t < 1; t += step)
		{
			double a1 = pow((1 - t), 2);
			double a2 = 2 * t*(1 - t);
			double a3 = t*t;
			apoint.x = a1*p1[0] + a2*p2[0] + a3*p3[0];
			apoint.y = a1*p1[1] + a2*p2[1] + a3*p3[1];
			vec->push_back(apoint);
		}
		apoint = mypoint2f(p3[0], p3[1]);
		if (vec->back() != apoint)
			vec->push_back(apoint);
	}
}
void MyGraph::getPListfromPolyline(int order, PolyLine aPLline, int density, vector<mypoint2f> *vec)
{
	int step = density + 1;
	if (aPLline.pointlist.size()/step<1)
	{
		step = 1;
	}
	if (order == 0)
	{	
		for (unsigned int i = 0; i < aPLline.pointlist.size(); i+=step)
		{
			vec->push_back(aPLline.pointlist.at(i));
		}
		if (vec->back() != aPLline.pointlist.back())
			vec->push_back(aPLline.pointlist.back());
	}
	else
	{
		int psize = aPLline.pointlist.size()-1;
		for (int i = psize; i >= 0; i -= step)
		{
			vec->push_back(aPLline.pointlist.at(i));
		}
		if (vec->back() != aPLline.pointlist.front())
			vec->push_back(aPLline.pointlist.front());
	/*	vector<mypoint2f>::iterator iter = aPLline.pointlist.end() - 1;
		for (iter; iter != aPLline.pointlist.begin(); iter-=step)
		{
			vec->push_back(*iter);
		}
		vec->push_back(aPLline.pointlist.front());*/
	}
}


//最一开始的对patch合并相关操作   getBorderofPatch()获取patch的一条边border，  getNeighborPatches()得到与某patch相邻的所有patch
void MyGraph::detectShorterEdge()
{
	cout << "enter detect Shorter Edge and merging" << endl;
	vector<pair<int, double>> ordered_patch;
	int psize = plrGraph.patches.size();
	for (int i = 0; i < psize; i++)
	{
		if (plrGraph.patches[i].type >= 0)
		{
			vector<double> edge_per;
			double mini_edge = 655350;//0
			int record_j = -1;
			for (unsigned int j = 1; j < plrGraph.patches[i].pointIndexOfPatch.size(); ++j)
			{
				//if (plrGraph.patches[i].type >= 0)
				//{
					vector<mypoint2f> edgeplist;
					int verindex1 = plrGraph.patches[i].pointIndexOfPatch.at(j - 1);
					int verindex2 = plrGraph.patches[i].pointIndexOfPatch.at(j);
					mypoint2f point1 = sgraph.nodeList[verindex1].position;//=lineWithCommonPoints[verindex1].first
					mypoint2f point2 = sgraph.nodeList[verindex2].position;
					getPListOfaCurve(point1, point2, verindex1, verindex2, &edgeplist,0);
					double anEdge_p = getPerimeterofVec(&edgeplist);
					int pchcount = getBorderofPatch(i, verindex1, verindex2)->adjPatchIndex.size();
					if (anEdge_p < mini_edge && pchcount>1) //最duan 且这条边有两个patch
					{
						mini_edge = anEdge_p;
						record_j = j;
					}
					edge_per.push_back(anEdge_p);
				//}				
			}
			if (record_j > 0)
			{
				plrGraph.patches[i].edges_vec[record_j - 1].edgeflag = 1;//edge标记为1，表示最短/最长边		
				double perimeter = 0;  //对所有patch的周长进行从大到小排序，插入法
				for (unsigned int k = 0; k < edge_per.size(); ++k)
				{
					perimeter += edge_per[k];
				}
				plrGraph.patches[i].perimeter = perimeter;//可能有更新
				bool inflag = false;
				vector<pair<int, double>>::iterator iter = ordered_patch.begin();
				while (iter != ordered_patch.end())
				{
					if (perimeter >iter->second)//从大到小排序
					{
						ordered_patch.insert(iter, make_pair(i, perimeter));//i是patch的index,perimeter是周长
						inflag = true;
						break;
					}
					iter++;
				}
				if (inflag == false)
					ordered_patch.push_back(make_pair(i, perimeter));
			}
			//int verindex1 = plrGraph.patches[i].edges_vec[record_j - 1].endp_index;
			//int verindex2 = plrGraph.patches[i].edges_vec[record_j - 1].startp_index;
			//if (sgraph.nodeList[verindex1].orderedArcs.size() == 2 || sgraph.nodeList[verindex2].orderedArcs.size() == 2)
			//	cout << "error" << endl;
		}
	}
	//对plrGraph.patches[...].type幅值，平均分为5个级别
	int addnum = ordered_patch.size() / 5;
	int start_i = 0;
	int segnum = addnum;
	for (int j = 0; j < 5; ++j)
	{
		for (int i = start_i; i < segnum; ++i)
		{
			int pindex = ordered_patch[i].first;
			if (plrGraph.patches[pindex].type>=0)//-1说明是repeat/overlay patch
				plrGraph.patches[pindex].type = j + 1;//patch等级1,2,3,4,5
		}
		start_i = segnum;
		if (j == 3)
			segnum = ordered_patch.size();
		else
			segnum = segnum + addnum;
	}
	mergeEdge();
}
EdgeOfPatch *MyGraph::getBorderofPatch(int patchindex, int startp, int endp)
{
	EdgeOfPatch* ep = NULL;
	for (unsigned int i = 0; i < plrGraph.patches[patchindex].edges_vec.size(); ++i)
	{
		if (plrGraph.patches[patchindex].edges_vec[i].startp_index == startp &&plrGraph.patches[patchindex].edges_vec[i].endp_index == endp)
		{
			ep = &plrGraph.patches[patchindex].edges_vec[i];
			break;
		}
		else if (plrGraph.patches[patchindex].edges_vec[i].endp_index == startp &&plrGraph.patches[patchindex].edges_vec[i].startp_index == endp)
		{
			ep = &plrGraph.patches[patchindex].edges_vec[i];
			break;
		}
	}
	return ep;
}
double MyGraph::getPerimeterofVec(vector<mypoint2f>* pvec)
{
	double len = 0;
	for (unsigned int i = 1; i < pvec->size(); ++i)
	{
		len = len + LineSegmentLength(pvec->at(i), pvec->at(i-1));
	}
	return len;
}
void MyGraph::mergeEdge()//大的吞并小的，大的是origipatch 小的是adjpatch,
{
	for (unsigned int i = 0; i < plrGraph.patches.size(); ++i)
	{
		if (plrGraph.patches[i].type < 5 && plrGraph.patches[i].type>0)//1234等级
		{
			int temp_i = i;//记录下当前i的值，
			int fromp = 0;// != plrGraph.patches[i].edges_vec[j].startp_index;//记录第i个patch第j条边的两个端点
			int endp = 0;// != plrGraph.patches[i].edges_vec[j].endp_index;
			int adgpatch_index = -1;//该条边相邻的另一个面的index
			for (unsigned int j = 0; j < plrGraph.patches[i].edges_vec.size(); ++j)
			{		
				if (plrGraph.patches[i].edges_vec[j].edgeflag == 1)//且第j条边的标记edgeflag=1
				{
					if (plrGraph.patches[i].edges_vec[j].adjPatchIndex.size()>1)//第j 条边的adjPatchIndex.size 有2个
					{
							if (plrGraph.patches[i].edges_vec[j].adjPatchIndex.at(1) != i)
							{
								adgpatch_index = plrGraph.patches[i].edges_vec[j].adjPatchIndex.at(1);
								if (plrGraph.patches[adgpatch_index].type < plrGraph.patches[i].type || plrGraph.patches[adgpatch_index].perimeter>plrGraph.patches[i].perimeter)
								{
									i = adgpatch_index;
									adgpatch_index = temp_i;
								}
								if (plrGraph.patches[adgpatch_index].type == -1)
								{
									cout << "adjpatch type =-1 ?" << endl;
									adgpatch_index = -1;
								}						
								break;
							}
							else
								cout << "adjPatchIndex error" << endl;
					}
				}	
			}
			if (adgpatch_index != -1)
			{
				plrGraph.patches[adgpatch_index].type = -1;//type, patch_index, edges_vec, pointIndexOfPatch。相邻面的type变为-1（表示已merge掉）
				vector<int> record_addp;//记录该相邻面的除共有点外其余顶点的index
				vector<int> remove_p;
				//找到两个patch公共的边
				vector<int> pindex_origi = plrGraph.patches[i].pointIndexOfPatch;
				pindex_origi.insert(pindex_origi.end(), plrGraph.patches[i].pointIndexOfPatch.begin() + 1, plrGraph.patches[i].pointIndexOfPatch.end() - 1);
				vector<int> pindex_adj = plrGraph.patches[adgpatch_index].pointIndexOfPatch;
				pindex_adj.insert(pindex_adj.end(), plrGraph.patches[adgpatch_index].pointIndexOfPatch.begin() + 1, plrGraph.patches[adgpatch_index].pointIndexOfPatch.end() - 1);
				reverse(pindex_adj.begin(), pindex_adj.end());
				getAddpFromAdjpatch(pindex_origi, pindex_adj, record_addp, fromp, endp);
				//未包含关系的patch
				//将这些点(record_addp)加入到第i个patch中的pointIndexOfPatch,采用方法：插入到fromp之后					
				if (!record_addp.empty())
				{
					reverse(record_addp.begin(), record_addp.end());//getAddpFromAdjpatch()之前reverse过，再翻转回来	
					vector<int>::iterator iter_pi = plrGraph.patches[i].pointIndexOfPatch.begin();
					if (*iter_pi != fromp)//对pointIndexOfPatch调整顺序,调整到以fromp开头的顺序
					{
						plrGraph.patches[i].pointIndexOfPatch.pop_back();
						while (*iter_pi != fromp)
						{
							int tmp = *iter_pi;
							iter_pi = plrGraph.patches[i].pointIndexOfPatch.erase(iter_pi);
							plrGraph.patches[i].pointIndexOfPatch.push_back(tmp);
						}
						plrGraph.patches[i].pointIndexOfPatch.push_back(fromp);
						iter_pi = plrGraph.patches[i].pointIndexOfPatch.begin();
					}
					//此时iter_pi指向pointIndexOfPatch的begin(),然后插入需要添加的点
					iter_pi = plrGraph.patches[i].pointIndexOfPatch.insert(iter_pi + 1, record_addp.begin(), record_addp.end());//insert(..begin()，val)是在第一个元素之前插入val，并返回val迭代器的位置
					iter_pi = iter_pi + record_addp.size();
					while (iter_pi != plrGraph.patches[i].pointIndexOfPatch.end())
					{
						if (*iter_pi != endp)
						{
							remove_p.push_back(*iter_pi);//记录被删除的点，公共点，除了两头
							iter_pi = plrGraph.patches[i].pointIndexOfPatch.erase(iter_pi);//执行这一步，说明公共边不止一个
						}
						else
							break;
					}
				}

				plrGraph.patches[adgpatch_index].patch_index = plrGraph.patches[i].patch_index;//!!!!						
				//将adgpatch的其余edge加入到第i个patch的edges_vec中。并且修改其相邻path的edge的信息。
				record_addp.push_back(endp);
				record_addp.insert(record_addp.begin(), fromp);
				vector<EdgeOfPatch>::iterator iter_oe = plrGraph.patches[i].edges_vec.begin();
				if (iter_oe->startp_index != fromp)//对edges_vec调整顺序
				{
					while (iter_oe->startp_index != fromp)
					{
						EdgeOfPatch tmp_e = *iter_oe;
						iter_oe = plrGraph.patches[i].edges_vec.erase(iter_oe);
						plrGraph.patches[i].edges_vec.push_back(tmp_e);
					}
				}
				iter_oe = plrGraph.patches[i].edges_vec.begin();
				while (iter_oe != plrGraph.patches[i].edges_vec.end())
				{
					if (iter_oe->startp_index == fromp)
					{
						for (int edge_i = 1; edge_i < record_addp.size(); edge_i++)
						{
							EdgeOfPatch *ep = getBorderofPatch(adgpatch_index, record_addp[edge_i - 1], record_addp[edge_i]);
							EdgeOfPatch newedge = *ep;
							if (ep->startp_index == record_addp[edge_i - 1])
							{
								newedge.adjPatchIndex[0] = plrGraph.patches[i].patch_index;//newedge.adjPatchIndex[0]
								if (newedge.adjPatchIndex.size()>1)
								{
									//修改其相邻patch的edge的信息
									int ano_pindex = newedge.adjPatchIndex[1];//newedge.adjPatchIndex[1]
									EdgeOfPatch* epch=getBorderofPatch(ano_pindex, newedge.endp_index, newedge.startp_index);
									if (epch->adjPatchIndex.size() > 1)//epch->adjPatchIndex
										epch->adjPatchIndex[1] = i;
									else
										epch->adjPatchIndex.push_back(i);
								}
								//修改相关node
								ArcNode *c_arc = getArcNode(ep->startp_index, ep->endp_index);
								c_arc->patchesIndex.push_back(i);
								vector<int>::iterator pch_iter = c_arc->patchesIndex.begin();
								while (pch_iter != c_arc->patchesIndex.end())
								{
									if (*pch_iter == adgpatch_index)
									{
										c_arc->patchesIndex.erase(pch_iter);
										break;
									}
									pch_iter++;
								}
								//反向修改
								c_arc = getArcNode(ep->endp_index, ep->startp_index);
								c_arc->patchesIndex.push_back(i);
								pch_iter = c_arc->patchesIndex.begin();
								while (pch_iter != c_arc->patchesIndex.end())
								{
									if (*pch_iter == adgpatch_index)
									{
										c_arc->patchesIndex.erase(pch_iter);
										break;
									}
									pch_iter++;
								}
								iter_oe = plrGraph.patches[i].edges_vec.insert(iter_oe, newedge);//依次累加，最后成了addedge+deleedge+origiedge
								iter_oe = iter_oe + 1;
							}
							else
								cout << "error" << endl;
						}
						//删除掉plrGraph.patches[i].edges_vec中原有的dele edge(公共边)
						while (iter_oe != plrGraph.patches[i].edges_vec.end())
						{
							if (iter_oe->endp_index != endp)
								iter_oe = plrGraph.patches[i].edges_vec.erase(iter_oe);
							else
								break;
						}
						if (iter_oe != plrGraph.patches[i].edges_vec.end())
							plrGraph.patches[i].edges_vec.erase(iter_oe);//删除掉iter_oe->endp_index = endp的那条边
						//更新周长
						vector<double> edge_per;
 						vector<mypoint2f> edgeplist;
						for (unsigned int iEdge = 1; iEdge < plrGraph.patches[i].pointIndexOfPatch.size(); ++iEdge)
						{
							int verindex1 = plrGraph.patches[i].pointIndexOfPatch.at(iEdge - 1);
							int verindex2 = plrGraph.patches[i].pointIndexOfPatch.at(iEdge);
							mypoint2f point1 = sgraph.nodeList[verindex1].position;//=lineWithCommonPoints[verindex1].first
							mypoint2f point2 = sgraph.nodeList[verindex2].position;
							getPListOfaCurve(point1, point2, verindex1, verindex2, &edgeplist,0);
							double anEdge_p = getPerimeterofVec(&edgeplist);
							edge_per.push_back(anEdge_p);
						}
						double perimeter = 0;  //对所有patch的周长进行从大到小排序，插入法
						for (unsigned int k = 0; k < edge_per.size(); ++k)
						{
							perimeter += edge_per[k];
						}
						plrGraph.patches[i].perimeter = perimeter;//可能有更新

						break;
					}
					else
						cout << "error" << endl;
					iter_oe++;
				}
			}
			//最后i一定变回原来的值!!!!!
			i = temp_i;
		}
	}

	////取每个patch较长的边来merge
	//for (unsigned int i = 0; i < plrGraph.patches.size(); ++i)
	//{
	//	if (plrGraph.patches[i].type <5 && plrGraph.patches[i].type>0)//1234等级
	//	{
	//		for (unsigned int j = 0; j < plrGraph.patches[i].edges_vec.size(); ++j)
	//		{
	//			if (plrGraph.patches[i].edges_vec[j].edgeflag == 1)
	//			{
	//				int temp_i = i;//记录下当前i的值，
	//				int fromp = 0;//= plrGraph.patches[i].edges_vec[j].startp_index;//记录第i个patch第j条边的两个端点
	//				int endp = 0;//= plrGraph.patches[i].edges_vec[j].endp_index;
	//				int adgpatch_index = 0;//该条边相邻的另一个面的index
	//				if (plrGraph.patches[i].edges_vec[j].adjPatchIndex.size()>1)
	//				{
	//					if (plrGraph.patches[i].edges_vec[j].adjPatchIndex.at(0) != i)
	//					{
	//						adgpatch_index = plrGraph.patches[i].edges_vec[j].adjPatchIndex.at(0);
	//					}
	//					else if (plrGraph.patches[i].edges_vec[j].adjPatchIndex.at(1) != i)
	//					{
	//						adgpatch_index = plrGraph.patches[i].edges_vec[j].adjPatchIndex.at(1);
	//					}
	//					if (plrGraph.patches[adgpatch_index].type == -1)
	//					{
	//						break;//应该不会有这种情况？？
	//						//adgpatch_index = plrGraph.patches[adgpatch_index].patch_index;
	//					}
	//						plrGraph.patches[adgpatch_index].type = -1;//type, patch_index, edges_vec, pointIndexOfPatch。相邻面的type变为-1（表示已merge掉）
	//					vector<int> record_addp;//记录该相邻面的除共有点外其余顶点的index
	//					vector<int> remove_p;
	//					//找到两个patch公共的边
	//					vector<int> pindex_origi = plrGraph.patches[i].pointIndexOfPatch;
	//					pindex_origi.insert(pindex_origi.end(), plrGraph.patches[i].pointIndexOfPatch.begin() + 1, plrGraph.patches[i].pointIndexOfPatch.end() - 1);
	//					vector<int> pindex_adj = plrGraph.patches[adgpatch_index].pointIndexOfPatch;
	//					pindex_adj.insert(pindex_adj.end(), plrGraph.patches[adgpatch_index].pointIndexOfPatch.begin() + 1, plrGraph.patches[adgpatch_index].pointIndexOfPatch.end() - 1);
	//					reverse(pindex_adj.begin(), pindex_adj.end());
	//					getAddpFromAdjpatch(pindex_origi, pindex_adj, record_addp,fromp,endp);
	//					
	//					////判断是否有包含关系
	//					//bool closeflag = false;
	//					//if (!record_addp.empty())
	//					//	reverse(record_addp.begin(), record_addp.end());//getAddpFromAdjpatch()之前reverse过，再翻转回来	
	//					//record_addp.push_back(endp);
	//					//record_addp.insert(record_addp.begin(), fromp);
	//					//for (unsigned int edge_i = 1; edge_i < record_addp.size(); edge_i++)
	//					//{
	//					//	EdgeOfPatch *ep = getBorderofPatch(adgpatch_index, record_addp[edge_i - 1], record_addp[edge_i]);
	//					//	if (ep->adjPatchIndex.size() == 1)
	//					//		break;
	//					//	else if (ep->adjPatchIndex.size() >1)
	//					//	{
	//					//		if (ep->adjPatchIndex[1] == plrGraph.patches[i].patch_index)//错了，不会出现这种情况，出现的话说明公共边没有找全！！
	//					//		{
	//					//			closeflag = true;
	//					//			plrGraph.patches[adgpatch_index].patch_index = plrGraph.patches[i].patch_index;
	//					//			//对edges_vec调整顺序
	//					//			vector<EdgeOfPatch>::iterator iter_adje = plrGraph.patches[adgpatch_index].edges_vec.begin();
	//					//			if (iter_adje->startp_index != endp)
	//					//			{
	//					//				while (iter_adje->startp_index != endp)
	//					//				{
	//					//					EdgeOfPatch tmp_e = *iter_adje;
	//					//					iter_adje = plrGraph.patches[adgpatch_index].edges_vec.erase(iter_adje);
	//					//					plrGraph.patches[adgpatch_index].edges_vec.push_back(tmp_e);
	//					//				}
	//					//			}
	//					//			iter_adje = plrGraph.patches[adgpatch_index].edges_vec.begin();
	//					//			while (iter_adje != plrGraph.patches[adgpatch_index].edges_vec.end())
	//					//			{
	//					//				if (iter_adje->endp_index != fromp)
	//					//				{
	//					//					EdgeOfPatch *tte = getBorderofPatch(iter_adje->adjPatchIndex[1], iter_adje->startp_index, iter_adje->endp_index);
	//					//					tte->adjPatchIndex[1] = plrGraph.patches[i].patch_index;
	//					//					iter_adje++;
	//					//				}
	//					//				else
	//					//					break;										
	//					//			}
	//					//			if (iter_adje != plrGraph.patches[adgpatch_index].edges_vec.end())
	//					//				getBorderofPatch(iter_adje->adjPatchIndex[1], iter_adje->startp_index, iter_adje->endp_index)->adjPatchIndex[1] = plrGraph.patches[i].patch_index;
	//					//			break;//处理完包围关系
	//					//		}
	//					//	}
	//					//}
	//					//if (closeflag == true)
	//					//	break;//不执行以下操作

	//					//未包含关系的patch
	//					//将这些点(record_addp)加入到第i个patch中的pointIndexOfPatch,采用方法：插入到fromp之后					
	//					if (!record_addp.empty())
	//					{				
	//						reverse(record_addp.begin(), record_addp.end());//getAddpFromAdjpatch()之前reverse过，再翻转回来	
	//						vector<int>::iterator iter_pi = plrGraph.patches[i].pointIndexOfPatch.begin();
	//						if (*iter_pi != fromp)//对pointIndexOfPatch调整顺序,调整到以fromp开头的顺序
	//						{
	//							plrGraph.patches[i].pointIndexOfPatch.pop_back();
	//							while (*iter_pi != fromp)
	//							{
	//								int tmp = *iter_pi;
	//								iter_pi = plrGraph.patches[i].pointIndexOfPatch.erase(iter_pi);
	//								plrGraph.patches[i].pointIndexOfPatch.push_back(tmp);
	//							}
	//							plrGraph.patches[i].pointIndexOfPatch.push_back(fromp);
	//							iter_pi = plrGraph.patches[i].pointIndexOfPatch.begin();
	//						}
	//						//此时iter_pi指向pointIndexOfPatch的begin()
	//						iter_pi = plrGraph.patches[i].pointIndexOfPatch.insert(iter_pi + 1, record_addp.begin(), record_addp.end());//insert(..begin()，val)是在第一个元素之前插入val，并返回val迭代器的位置
	//						iter_pi = iter_pi + record_addp.size();
	//						while (iter_pi != plrGraph.patches[i].pointIndexOfPatch.end())
	//						{
	//							if (*iter_pi != endp)
	//							{
	//								remove_p.push_back(*iter_pi);//记录被删除的点，公共点，除了两头
	//								iter_pi = plrGraph.patches[i].pointIndexOfPatch.erase(iter_pi);//执行这一步，说明公共边不止一个
	//							}
	//							else
	//								break;
	//						}
	//					}	

	//						plrGraph.patches[adgpatch_index].patch_index = plrGraph.patches[i].patch_index;//!!!!						
	//						//将adgpatch的其余edge加入到第i个patch的edges_vec中。并且修改其相邻path的edge的信息。
	//						record_addp.push_back(endp);
	//						record_addp.insert(record_addp.begin(), fromp);
	//						vector<EdgeOfPatch>::iterator iter_oe = plrGraph.patches[i].edges_vec.begin();
	//						if (iter_oe->startp_index != fromp)//对edges_vec调整顺序
	//						{
	//							while (iter_oe->startp_index != fromp)
	//							{
	//								EdgeOfPatch tmp_e = *iter_oe;
	//								iter_oe = plrGraph.patches[i].edges_vec.erase(iter_oe);
	//								plrGraph.patches[i].edges_vec.push_back(tmp_e);
	//							}
	//						}
	//						iter_oe = plrGraph.patches[i].edges_vec.begin();
	//						while (iter_oe != plrGraph.patches[i].edges_vec.end())
	//						{
	//							if (iter_oe->startp_index == fromp)
	//							{
	//								for (int edge_i = 1; edge_i < record_addp.size(); edge_i++)
	//								{
	//									EdgeOfPatch *ep = getBorderofPatch(adgpatch_index, record_addp[edge_i - 1], record_addp[edge_i]);
	//									EdgeOfPatch newedge = *ep;
	//									if (ep->startp_index == record_addp[edge_i - 1])
	//									{
	//										newedge.adjPatchIndex[0] = plrGraph.patches[i].patch_index;//newedge.adjPatchIndex[0]
	//										if (newedge.adjPatchIndex.size()>1)
	//										{
	//											//修改其相邻patch的edge的信息
	//											int ano_pindex = newedge.adjPatchIndex[1];//newedge.adjPatchIndex[1]
	//											getBorderofPatch(ano_pindex, newedge.endp_index, newedge.startp_index)->adjPatchIndex[1] = plrGraph.patches[i].patch_index;
	//										}
	//										iter_oe = plrGraph.patches[i].edges_vec.insert(iter_oe, newedge);//依次累加，最后成了addedge+deleedge+origiedge
	//										iter_oe = iter_oe + 1;
	//									}
	//									else
	//										cout << "error" << endl;
	//								}
	//								while (iter_oe != plrGraph.patches[i].edges_vec.end())
	//								{
	//									if (iter_oe->endp_index != endp)
	//										iter_oe = plrGraph.patches[i].edges_vec.erase(iter_oe);
	//									else
	//										break;
	//								}
	//								if (iter_oe != plrGraph.patches[i].edges_vec.end())
	//									plrGraph.patches[i].edges_vec.erase(iter_oe);//删除掉iter_oe->endp_index = endp的那条边
	//								break;
	//							}
	//							else
	//								cout << "error" << endl;
	//							iter_oe++;
	//						}
	//					break;//plrGraph.patches[i].edges_vec已经改变
	//				}
	//				else//adjpatch size()=1 说明这个patch的这条长边没有临街patch，不需要merge，但是type不能变，edgeflag也不能变，有可能其他patch与他merge时，需要先判断一下type谁大谁小
	//				{
	//					//plrGraph.patches[i].type = 0;
	//					//plrGraph.patches[i].edges_vec[j].edgeflag = 0;
	//				}
	//			}
	//		}
	//	}
	//}
	cout << "merge done" << endl;
}
bool MyGraph::getAddpFromAdjpatch(vector<int> pvec_origi, vector<int> pvec_adj, vector<int> &addp, int &startp, int &top)
{//参数pvec_origi和pvec_adj都是加长的，且pvec_adj是倒序的
	bool sucflag = false;
	//vector<int> commp,commp2;
	bool firstflag = false;
	vector<int> delepoint;//delepoint即两个patch之间的common point，包括了起始端点//参数addp是originalpatch应该添加的adjpatch的那部分边（是从倒序的adjptch中提取的，但对于origipatch来说是也倒序的，此函数之外还需要reverse）
	for (unsigned int k_adj = 0; k_adj < pvec_adj.size(); k_adj++)
	{
		for (unsigned int k_origi = 0; k_origi < pvec_origi.size(); ++k_origi)
		{
			if (pvec_adj[k_adj] == pvec_origi[k_origi])
			{
				startp = pvec_adj[k_adj];//!!!startp是指沿origipatch的方向，与adjpatch公共边的起始点。同理top是指公共边的终点。这两个变量是对origipatch进行添加和删除操作的标记位
				int add_i = 0;//add_i从第一个相等的起始点开始计数，到最后一个相等点的下一个结束。add_i等于commonpoint(delepoint)的个数
				bool sflag = false;
				while ((k_adj + add_i) < pvec_adj.size() && (k_origi + add_i) < pvec_origi.size() && sflag == false)
				{
					if (pvec_adj[k_adj + add_i] == pvec_origi[k_origi + add_i])
					{
						delepoint.push_back(pvec_adj[k_adj + add_i]);//注意delepoint包含了startp 和 top
						add_i++;
					}
					else
						sflag = true;
				}
				if (add_i > 1)//找到common point(delepoint),//若add_i=1则跳出循环，选取第下一个adjpatch的点并依次遍历origipatch，在下次循环中还会找到相同的点并记录delepoint
				{
					if (k_adj != 0)//if (k_adj != 0 && firstflag==false)
					{
						top = pvec_adj[k_adj + add_i - 1];//或 pindex_origi[k_origi + add_i - 1];//!!!!!!
						while ((k_adj + add_i) < pvec_adj.size())
						{
							if (pvec_adj[k_adj + add_i] != startp)
							{
								addp.push_back(pvec_adj[k_adj + add_i]);//注意得到的addp，对于adjpatch是逆序的，若要加到origipatch中还需要reverse一下
								add_i++;
							}
							else
								break;
						}
						sucflag = true;//record_addp可能是空的,光靠判断addp是否不为空来结束该函数是不行的，所以设立了标志位sucflag
						break;
					}
					else//如果k_adj=0，得到的recor_addp可能大于实际值。即可能有一部分要删除的commonpoint没有找到，delepoint不全
					{
						top = pvec_adj[k_adj + add_i - 1];//!!!!!!
						int tmp = 1;
						bool equalflag = false;
						while (equalflag == false)
						{
							if (((int)k_origi - tmp) < 0)//pvec_origi之前插入一个值，即循环移动pvec_origi
							{
								int aback = pvec_origi.back();
								pvec_origi.pop_back();
								pvec_origi.insert(pvec_origi.begin(), aback);
								k_origi = (unsigned int)tmp;
							}
							if (((int)k_adj - tmp) <0)//pvec_adj之后
							{
								int aback = pvec_adj.back();
								pvec_adj.pop_back();
								pvec_adj.insert(pvec_adj.begin(), aback);
								k_adj = (unsigned int)tmp;
							}
							if (pvec_adj[(int)k_adj - tmp] == pvec_origi[(int)k_origi - tmp])//一步一步的向delepoint添加遗漏的，并且更新startp起始点等等
							{
								startp = pvec_adj[(int)k_adj - tmp];
								add_i++;
								delepoint.insert(delepoint.begin(), startp);
								tmp++;
							}
							else
								equalflag = true;	//找到完整的commonpoint(delepoint),结束					
						}
						k_adj = (unsigned int)((int)k_adj - tmp + 1);////因为上面while循环改变了pvec_adj，所以k_adj也要变
						while ((k_adj + add_i) < pvec_adj.size())//有的delepoint,所以可以找到adjpatch中其余的点，即需要添加到origipatch中的点列
						{
							if (pvec_adj[(int)k_adj + add_i] != startp)
							{
								addp.push_back(pvec_adj[(int)k_adj + add_i]);
								add_i++;
							}
							else
								break;
						}
						sucflag = true;//record_addp可能是空的		
						break;

						//if (firstflag == false)
						//{
						//	firstflag = true;
						//	first_addsize = add_i;
						//	fir_k_posi = k_adj;
						//	int k = 0;
						//	while (k < add_i)
						//	{
						//		commp.push_back(pvec_adj[k_adj + k]);
						//		k++;
						//	}
						//	break;
						//}
						//else
						//{		
						//	int k = 0;
						//	while (k < add_i)
						//	{
						//		commp2.push_back(pvec_adj[k_adj + k]);
						//		k++;
						//	}
						//	bool valiflag = false;
						//	for (unsigned int c_1 = 0; c_1 < commp.size(); ++c_1)////
						//	{
						//		for (unsigned int c_2 = 0; c_2 < commp2.size(); ++c_2)
						//		{
						//			if (commp[c_1] == commp2[c_2])
						//			{
						//				valiflag = true;
						//				break;
						//			}
						//		}
						//	}
						//	if (valiflag == true)
						//	{
						//		sucflag = true;
						//		if (add_i >= first_addsize)
						//		{									
						//			top = pvec_adj[k_adj + add_i - 1];//或 pindex_origi[k_origi + add_i - 1];//!!!!!!
						//			while ((k_adj + add_i) < pvec_adj.size())
						//			{
						//				if (pvec_adj[k_adj + add_i] != startp)
						//				{
						//					addp.push_back(pvec_adj[k_adj + add_i]);
						//					add_i++;
						//				}
						//				else
						//					break;
						//			}
						//		}
						//		else
						//		{
						//			add_i = first_addsize;
						//			k_adj = fir_k_posi;
						//			startp = pvec_adj[k_adj];
						//			top = pvec_adj[k_adj + add_i - 1];//或 pindex_origi[k_origi + add_i - 1];//!!!!!!
						//			while ((k_adj + add_i) < pvec_adj.size())
						//			{
						//				if (pvec_adj[k_adj + add_i] != startp)
						//				{
						//					addp.push_back(pvec_adj[k_adj + add_i]);
						//					add_i++;
						//				}
						//				else
						//					break;
						//			}
						//		}
						//	}
						//	
						//	break;
						//}
					}
				}
				else
					break;
			}
		}//
		if (!addp.empty() || sucflag == true)
		{
			break;
		}
	}
	/*cout << "common point:";
	for (unsigned int i = 0; i < delepoint.size(); ++i)
	{
		cout << delepoint[i] << ",";
	}
	cout << endl;*/

	//check addp中是否有其他点和pvec_origi中的点重复	
	bool result = false;
	if (sucflag == true)
	{
		if (!addp.empty())
		{
			bool check_flag = false;
			int halfsize = pvec_origi.size() / 2;
			for (unsigned int i = 0; i < halfsize; i++)
			{
				for (unsigned int j = 0; j < addp.size(); ++j)
				{
					if (addp[j] == pvec_origi[i])
					{
						check_flag = true;
						break;
					}
				}
				if (check_flag == true)
					break;
			}
			if (check_flag == true)//又找到相同的点，说明combine有问题，返回false
				result = false;
			else
				result = true;
		}
		else
			result = true;
	}
	else
		result = false;
	return result;
}
vector<int> MyGraph::getNeighborPatches(int patchindex, int flag)
{
	vector<int> neipth;
	if (flag == 0)
	{
		//求只有公共边的patch
		for (unsigned int i = 0; i<plrGraph.patches[patchindex].edges_vec.size(); ++i)
		{
			if (plrGraph.patches[patchindex].edges_vec[i].adjPatchIndex.size()>1)
			{
				int adj_p = plrGraph.patches[patchindex].edges_vec[i].adjPatchIndex[1];
				if (adj_p != patchindex && plrGraph.patches[adj_p].type >= 0)
					neipth.push_back(adj_p);
			}
		}
	}
	else if (flag==1)
	{
		//公共边patch + 对角线方向的patch
		for (unsigned int i = 0; i < plrGraph.patches[patchindex].pointIndexOfPatch.size() - 1; ++i)
		{
			int nodeindex = plrGraph.patches[patchindex].pointIndexOfPatch[i];
			ArcNode *anode = sgraph.nodeList[nodeindex].firstEdge;
			while (anode)
			{
				if (!anode->patchesIndex.empty())
				{
					for (unsigned int j = 0; j < anode->patchesIndex.size(); ++j)
					{
						if (anode->patchesIndex[j] != patchindex && plrGraph.patches[anode->patchesIndex[j]].type >= 0)
							neipth.push_back(anode->patchesIndex[j]);
					}
				}
				anode = anode->next;
			}
			/*for (unsigned int j = 0; j < sgraph.nodeList[nodeindex].orderedArcs.size(); ++j)
			{
			if (!sgraph.nodeList[nodeindex].orderedArcs[j].patchesIndex.empty())
			neipth.insert(neipth.begin(), sgraph.nodeList[nodeindex].orderedArcs[j].patchesIndex.begin(), sgraph.nodeList[nodeindex].orderedArcs[j].patchesIndex.end());
			}*/
		}
	}
	deleRepeatPoint(&neipth);
	return neipth;
}

//patch topology的相关操作
void MyGraph::ConstructPatchTopology()
{
	int pchsize=plrGraph.patches.size();
	patchTopology.patch_number = pchsize;
	for (int i = 0; i < pchsize; ++i)
	{
		//不管patch的type是不是-1 ，都存入patch topology(vec)中。若是type=-1，这就是一个孤立的点，可以直接用patchindex找到
		PatchNode apnode;
		apnode.patchType = plrGraph.patches[i].type;
		apnode.patch_index = plrGraph.patches[i].patch_index;
		apnode.firstAdjPatch = NULL;
		if (apnode.patchType >= 0)
		{
			vector<int> adjpatch;
			for (unsigned int j = 0; j < plrGraph.patches[i].edges_vec.size(); ++j)
			{
				adjpatch.insert(adjpatch.end(), plrGraph.patches[i].edges_vec[j].adjPatchIndex.begin(), plrGraph.patches[i].edges_vec[j].adjPatchIndex.end());
			}
			deleRepeatPoint(&adjpatch);
			for (unsigned int j = 0; j < adjpatch.size(); ++j)
			{
				if (adjpatch[j] != i)
				{
					PatchArc *anarc = new PatchArc;
					anarc->nextPatch = NULL;
					anarc->adj_patchindex = adjpatch[j];
					anarc->weight = 0;
					anarc->nextPatch = apnode.firstAdjPatch;
					apnode.firstAdjPatch = anarc;
				}
			}
		}
		//计算每个patch的重心，先得到plist,再求均值
		//vector<mypoint2f> pch_plist;
		//for (unsigned int j = 1; j < plrGraph.patches[i].pointIndexOfPatch.size(); ++j)
		//{
		//	vector<mypoint2f> edgeplist;
		//	int verindex1 = plrGraph.patches[i].pointIndexOfPatch.at(j - 1);
		//	int verindex2 = plrGraph.patches[i].pointIndexOfPatch.at(j);
		//	mypoint2f point1 = sgraph.nodeList[verindex1].position;//=lineWithCommonPoints[verindex1].first
		//	mypoint2f point2 = sgraph.nodeList[verindex2].position;
		//	getPListOfaCurve(point1, point2, verindex1, 0, &edgeplist);//only know the startpoint's and endpoint's coordination and startpoint'sindex.need to find which kind of curves it belongs to,like cubicBezier,quadraticBezier...and store curve's pointList in 'apatch'
		//	pch_plist.insert(pch_plist.end(), edgeplist.begin(), edgeplist.end());
		//}
		//double aver_x = 0;
		//double aver_y = 0;
		//for (unsigned int j = 0; j < pch_plist.size(); ++j)
		//{
		//	aver_x = aver_x + pch_plist[j].x;
		//	aver_y = aver_y + pch_plist[j].y;
		//}
		//aver_x = aver_x / (double)pch_plist.size();
		//aver_y = aver_y / (double)pch_plist.size();
		//apnode.centre_position = mypoint2f(aver_x, aver_y);
		patchTopology.allpatches.push_back(apnode);
	}
	//划分patch的子图
	int *visitStatus = new int[pchsize];//0:unvisited,1:push in queue,2: out of queue
	for (int i = 0; i < pchsize; ++i)
	{
		visitStatus[i] = 0;
	}
	for (int i = 0; i < pchsize; ++i)//breadth-first traversal  BFT
	{
		if (visitStatus[i] == 0 && patchTopology.allpatches[i].patchType>=0)
		{
			if (patchTopology.allpatches[i].firstAdjPatch != NULL)
			{
				int arcCount = 0;
				int verCount = 0;
				int tmp_ver = 0;
				queue<int> verQueue;
				verQueue.push(i);
				visitStatus[i] = 1;
				while (!verQueue.empty())
				{
					tmp_ver = verQueue.front();
					verQueue.pop();
					verCount = verCount + 1;
					visitStatus[tmp_ver] = 2;
					PatchArc *tmp_arc = patchTopology.allpatches[tmp_ver].firstAdjPatch;
					while (tmp_arc)
					{
						if (visitStatus[tmp_arc->adj_patchindex] == 0)
						{
							verQueue.push(tmp_arc->adj_patchindex);
							visitStatus[tmp_arc->adj_patchindex] = 1;
							arcCount = arcCount + 1;
						}
						else if (visitStatus[tmp_arc->adj_patchindex] == 1)
						{
							arcCount = arcCount + 1;
						}
						tmp_arc = tmp_arc->nextPatch;
					}
				}
				if (verCount - 1 == arcCount)
				{
					patch_subgraphs.push_back(make_pair(i, "tree"));
				}
				else
				{
					patch_subgraphs.push_back(make_pair(i, "loop"));
				}
			}
			else
				patch_subgraphs.push_back(make_pair(i, "node"));
		}
	}
	delete[] visitStatus;
	visitStatus = NULL;
}
int MyGraph::countPatchEdgeNum(int pchindex)
{
	int count = 0;
	PatchArc *arcnode = patchTopology.allpatches[pchindex].firstAdjPatch;
	while (arcnode)
	{
		count = count + 1;
		arcnode = arcnode->nextPatch;
	}
	return count;
}
void MyGraph::removePatchNode(int hostpch, int delpch)
{
	//int cc = countEdgeNum(startVerIndex);
	//int bb = countEdgeNum(endVerIndex);
	PatchArc *arcnode = patchTopology.allpatches[hostpch].firstAdjPatch;
	PatchArc *prior_arc = NULL;
	while (arcnode)
	{
		if (arcnode->adj_patchindex == delpch)
		{
			if (prior_arc == NULL)
			{
				patchTopology.allpatches[hostpch].firstAdjPatch = arcnode->nextPatch;
			}
			else
			{
				prior_arc->nextPatch = arcnode->nextPatch;
			}
			delete arcnode;
			arcnode = NULL;
			break;
		}
		prior_arc = arcnode;
		arcnode = arcnode->nextPatch;
	}
}
PatchArc* MyGraph::getPatchArc(int pindex, int nei_pindex)
{
	PatchArc *tmparc = NULL;
	PatchArc* arc=patchTopology.allpatches[pindex].firstAdjPatch;
	while (arc != NULL)
	{
		if (arc->adj_patchindex == nei_pindex)
		{
			tmparc = arc;
			break;
		}
		arc = arc->nextPatch;
	}
	return tmparc;
}
void MyGraph::deleUndesirablePatch()//在下面的KmeansClustering()函数之后执行
{
	/*int subNum = patch_subgraphs.size();
	for (int i = 0; i < subNum; ++i)
	{
		int startpch = patch_subgraphs[i].first;
	}*/
	for (unsigned int j = 0; j < plrGraph.patches.size(); ++j)
	{
		//一开始patch.type的初始值是0，某些重复的patch是-1。执行KmeansClustering()之后，patch被分为1类和2类。之后令噪声点和其他不需要的patch的type又被赋值为0，所以现在只要是0的patch就delete掉！！
		if (plrGraph.patches[j].type == 0)
		{
			int pch = plrGraph.patches[j].patch_index;
			int pch_dele = 0;
			PatchArc *tmpnode = patchTopology.allpatches[pch].firstAdjPatch;
			while (tmpnode!=NULL)
			{
				pch_dele = tmpnode->adj_patchindex;
				tmpnode = tmpnode->nextPatch;
				removePatchNode(pch, pch_dele);
				removePatchNode(pch_dele, pch);
			}
		}
	}
	//原来的graph被破坏了，重新构建子图
	int pchsize = plrGraph.patches.size();
	for (int i = 0; i < pchsize; ++i)
	{
		//修改patch的type
		patchTopology.allpatches[i].patchType = plrGraph.patches[i].type;
	}
	//重新划分patch的子图，patch_subgraphs--->patch_subgraphs_new
	int *visitStatus = new int[pchsize];
	//0:unvisited,1:push in queue,2: out of queue
	for (int i = 0; i < pchsize; ++i)
	{
		visitStatus[i] = 0;
	}
	for (int i = 0; i < pchsize; ++i)//breadth-first traversal  BFT
	{
		if (visitStatus[i] == 0 && patchTopology.allpatches[i].patchType > 0)
		{
			if (patchTopology.allpatches[i].firstAdjPatch == NULL)//一个孤立的点
			{
				patch_subgraphs_group.push_back(make_pair(i, "node"));
			}
			else
			{
				int arcCount = 0;
				int verCount = 0;
				int tmp_ver = 0;
				queue<int> verQueue;
				verQueue.push(patchTopology.allpatches[i].patch_index);
				visitStatus[patchTopology.allpatches[i].patch_index] = 1;
				while (!verQueue.empty())
				{
					tmp_ver = verQueue.front();
					verQueue.pop();
					verCount = verCount + 1;
					visitStatus[tmp_ver] = 2;
					PatchArc *tmp_arc = patchTopology.allpatches[tmp_ver].firstAdjPatch;
					while (tmp_arc)
					{
						if (visitStatus[tmp_arc->adj_patchindex] == 0)
						{
							verQueue.push(tmp_arc->adj_patchindex);
							visitStatus[tmp_arc->adj_patchindex] = 1;
							arcCount = arcCount + 1;
						}
						else if (visitStatus[tmp_arc->adj_patchindex] == 1)
						{
							arcCount = arcCount + 1;
						}
						tmp_arc = tmp_arc->nextPatch;
					}
				}
				if (verCount - 1 == arcCount)
				{
					patch_subgraphs_group.push_back(make_pair(patchTopology.allpatches[i].patch_index, "tree"));
				}
				else
				{
					patch_subgraphs_group.push_back(make_pair(patchTopology.allpatches[i].patch_index, "loop"));
				}
			}	
		}
	}
	delete[] visitStatus;
	visitStatus = NULL;
}
void MyGraph::jointCurveTPatch_subgraph()
{
	//另一种方法通过patch_subgraphs输出patch ，一个subgraph是一个颜色。前提是patchTopology已经被初始化，patch_subgraphs也被初始化
	int pchsize = plrGraph.patches.size();
	int *visitStatus = new int[pchsize];
	for (int i = 0; i < pchsize; ++i)	//0:unvisited,1:push in queue,2: out of queue
	{
		visitStatus[i] = 0;
	}
	int typenum = 1;
	for (unsigned int i = 0; i < patch_subgraphs_group.size(); ++i)
	{
		if (patch_subgraphs_group[i].second != "node")
		{
			typenum++;
			int startnum = patch_subgraphs_group[i].first;
			int tmp_ver = 0;
			queue<int> verQueue;
			verQueue.push(startnum);
			visitStatus[startnum] = 1;

			vector<mypoint2f> apatch_start;
			for (unsigned int j = 1; j < plrGraph.patches[startnum].pointIndexOfPatch.size(); ++j)
			{
				int verindex1 = plrGraph.patches[startnum].pointIndexOfPatch.at(j - 1);
				int verindex2 = plrGraph.patches[startnum].pointIndexOfPatch.at(j);
				mypoint2f point1 = sgraph.nodeList[verindex1].position;//=lineWithCommonPoints[verindex1].first
				mypoint2f point2 = sgraph.nodeList[verindex2].position;
				if (!apatch_start.empty())
					apatch_start.pop_back();
				getPListOfaCurve(point1, point2, verindex1, verindex2, &apatch_start,0);
				////添加组成patch的节点， 
				//mycircle acir;
				//acir.cx = point1.x; acir.cy = point1.y; acir.radius = 0.5;
				//circle_vec.push_back(acir);
			}
			if (isPointRoughEqual(apatch_start.front(), apatch_start.back()))
				apatch_start.pop_back();
			pointListOfPatches.push_back(make_pair(typenum, apatch_start));//plrGraph.patches.push_back(apatch);

			while (!verQueue.empty())
			{
				tmp_ver = verQueue.front();
				verQueue.pop();
				visitStatus[tmp_ver] = 2;
				PatchArc *tmp_arc = patchTopology.allpatches[tmp_ver].firstAdjPatch;
				while (tmp_arc)
				{
					if (visitStatus[tmp_arc->adj_patchindex] != 2)
					{
						if (visitStatus[tmp_arc->adj_patchindex] == 0)
						{
							verQueue.push(tmp_arc->adj_patchindex);
							visitStatus[tmp_arc->adj_patchindex] = 1;
							vector<mypoint2f> apatch;
							for (unsigned int j = 1; j < plrGraph.patches[tmp_arc->adj_patchindex].pointIndexOfPatch.size(); ++j)
							{
								int verindex1 = plrGraph.patches[tmp_arc->adj_patchindex].pointIndexOfPatch.at(j - 1);
								int verindex2 = plrGraph.patches[tmp_arc->adj_patchindex].pointIndexOfPatch.at(j);
								mypoint2f point1 = sgraph.nodeList[verindex1].position;//=lineWithCommonPoints[verindex1].first
								mypoint2f point2 = sgraph.nodeList[verindex2].position;
								if (!apatch.empty())
									apatch.pop_back();
								getPListOfaCurve(point1, point2, verindex1, verindex2, &apatch,0);
								////添加组成patch的节点， 
								//mycircle acir;
								//acir.cx = point1.x; acir.cy = point1.y; acir.radius = 0.5;
								//circle_vec.push_back(acir);
							}
							if (isPointRoughEqual(apatch.front(), apatch.back()))
								apatch.pop_back();
							pointListOfPatches.push_back(make_pair(typenum, apatch));//plrGraph.patches.push_back(apatch);
						}	
					}
					tmp_arc = tmp_arc->nextPatch;
				}
			}
		}
	}
	delete[] visitStatus;
	visitStatus = NULL;
}

//用聚类的方法吧patch进行分类，先getfeature()-->kmeans/其他
void MyGraph::getPatchFeatures()
{
	//1、圆度，公式4*pi*area/permeter^2. 如果是圆，改特征值=1
	//2、长度(elongation)、(离心率eccentricity)、length of major axis/length of minor axis
	//周长or面积算吗？

	//vector<pair<int,vector<double>>> patchfeatures;//int 对应patchindex，vector<double>对应特征向量
	cout << "getPatchFeatures()" << endl;
	vector<double> cir_vec;
	vector<double> long_vec;
	vector<double> area_vec;
	int pchsize = plrGraph.patches.size();
	//计算这个图形所有patch的总面积
	double total_area = 0; double min_area = INT_MAX; double max_area = 0;
	for (int i = 0; i < pchsize; ++i)
	{
		//if (plrGraph.patches[i].type >= 0)
		//{
		////	total_area = total_area + plrGraph.patches[i].area;
		//	if (plrGraph.patches[i].area>max_area)
		//		max_area = plrGraph.patches[i].area;
		//	if (plrGraph.patches[i].area < min_area)
		//		min_area = plrGraph.patches[i].area;
		//}
		if (plrGraph.patches[i].type >= 0)
		{
			if (plrGraph.patches[i].area>max_area)
				max_area = plrGraph.patches[i].area;
			if (plrGraph.patches[i].area < min_area)
				min_area = plrGraph.patches[i].area;
		}
	}
	for (int i = 0; i < pchsize; ++i)
	{
		if (plrGraph.patches[i].type >= 0)
		{
			vector<mypoint2f> total_plist;
			jointMethod(&plrGraph.patches[i].pointIndexOfPatch, &total_plist, 0);// 0是density，返回的结果total_plist首=尾
			total_plist.pop_back();
			/*1、圆度*/
				double circularity = 4 * 3.1415926*plrGraph.patches[i].area / pow(plrGraph.patches[i].perimeter, 2);
			/*2、长度,方法一（公式来自论文，直接计算出比值）*/
				//mypoint2f centrol = plrGraph.patches[i].centre_position;
				//double cxx = 0; double cyy = 0; double cxy = 0; double cyx = 0;
				//for (unsigned int j = 0; j < total_plist.size(); ++j)
				//{
				//	cxx = cxx + pow(total_plist[j].x - centrol.x, 2);
				//	cyy = cyy + pow(total_plist[j].y - centrol.y, 2);
				//	cxy = cxy + (total_plist[j].x - centrol.x)*(total_plist[j].y - centrol.y);
				//}
				//cxx = cxx / total_plist.size();
				//cyy = cyy / total_plist.size();
				//cxy = cxy / total_plist.size();
				//double tmp1 = cxx + cyy;
				//double tmp2 = cxx*cyy;
				//double tmp3 = sqrt((pow(tmp1, 2) - 4 * (tmp2 - pow(cxy, 2))));
				//double elongation = (tmp1 - tmp3) / (tmp1 + tmp3);
				
			/*2、长度，方法二，计算出最大矩形bounding box*/
				//PCA 主成分分析（ 直角坐标系） ，先计算协方差矩阵A ,初始化 数据点矩阵MatrixXd B（2，n） 数据按列排列！！
				int plist_size = total_plist.size();
				vector<double> vec_x;	//vec_x，vec_y，B 都是直角坐标系下的，
				vector<double> vec_y;
				MatrixXd B = MatrixXd::Zero(2, plist_size);
				MatrixXd B_revolve = MatrixXd::Zero(2, plist_size);
				for (int j = 0; j < plist_size; ++j)
				{
					vec_x.push_back(total_plist[j].x);
					vec_y.push_back(total_plist[j].y);
					B(0, j) = total_plist[j].x;
					B(1, j) = total_plist[j].y;
				}
				//cout << B[0] << endl;
				double aver_x = GetAverange(vec_x);
				double aver_y = GetAverange(vec_y);
				double var_cxx = GetVariance(vec_x, aver_x, vec_x, aver_x);
				double var_cyy = GetVariance(vec_y, aver_y, vec_y, aver_y);
				double var_cxy = GetVariance(vec_x, aver_x, vec_y, aver_y);
				Matrix2d A;
				A << var_cxx, var_cxy, var_cxy, var_cyy;  //cout << "Here is a covariance 2*2 matrix, A:" << endl << A << endl << endl;
				//计算协方差矩阵的特征值 特征向量
				EigenSolver<MatrixXd> es(A);
				cout << "The eigenvalues of A are:" << endl << es.eigenvalues() << endl;//注意特征向量和特征值都是 复数类型！！complex
				cout << "The eigenvectors is:" << endl << es.eigenvectors() << endl;//特征向量是按列排列的,且已经是单位化的
				//现将特征向量变成double类型的matrix，
				VectorXd an_vec;
				Matrix2d eigenvec_double = MatrixXd::Zero(2, 2);
				for (int m_i = 0; m_i<2; m_i++)
				{
					an_vec = es.eigenvectors().col(m_i).real();
					eigenvec_double.col(m_i) = an_vec;
				}
				// B中的数据点乘以特征向量矩阵的转置，进行坐标变换
				B_revolve = eigenvec_double.transpose()*B;  //(2*2) * (2*n)， eigenvec_double.transpose()返回转置后的矩阵，eigenvec_double本身没有变
				//求在新坐标系下的 上下左右 边界 
				Vector2d apt_vec;//将某一点进行坐标变换
				apt_vec << total_plist[0].x, total_plist[0].y;
				Vector2d apt_revolve = eigenvec_double.transpose()*apt_vec;
				double posi_right = apt_revolve[0]; double posi_left = apt_revolve[0];
				double posi_top = apt_revolve[1]; double posi_bottom = apt_revolve[1];
				for (int j = 0; j < plist_size; ++j)
				{
					if (B_revolve(0, j) < posi_left)
						posi_left = B_revolve(0, j);
					if (B_revolve(0, j) > posi_right)
						posi_right = B_revolve(0, j);
					if (B_revolve(1, j) < posi_bottom)
						posi_bottom = B_revolve(1, j);
					if (B_revolve(1, j) > posi_top)
						posi_top = B_revolve(1, j);					
				}
				//计算新坐标系下的 原点,
				mypoint2f revol_center_pt(0, 0);
				revol_center_pt.x = (posi_left + posi_right) / 2;
				revol_center_pt.y = (posi_bottom + posi_top) / 2;
				//计算 长宽 和 比例
				double box_w = fabs(posi_right - posi_left);
				double box_h = fabs(posi_top - posi_bottom);
				double elongation = 0;
				if (box_w < box_h)
					elongation = box_w / box_h;
				else
					elongation = box_h / box_w;
				/*showBoundingBox(),计算bounding box四个点*/
				//vector<mypoint2f> plist_box;
				//mypoint2f bp(revol_center_pt.x - box_w / 2, revol_center_pt.y + box_h / 2);//左上角
				//plist_box.push_back(bp);
				//bp.x = revol_center_pt.x + box_w / 2;  bp.y = revol_center_pt.y + box_h / 2;//右上角
				//plist_box.push_back(bp);
				//bp.x = revol_center_pt.x + box_w / 2;  bp.y = revol_center_pt.y - box_h / 2;//右下角
				//plist_box.push_back(bp);
				//bp.x = revol_center_pt.x - box_w / 2;  bp.y = revol_center_pt.y - box_h / 2;//左下角
				//plist_box.push_back(bp); 
				//bp.x = revol_center_pt.x - box_w / 2;  bp.y = revol_center_pt.y + box_h / 2;
				//plist_box.push_back(bp);
				////4点坐标变换到 原来坐标系，
				//MatrixXd box_matrix = MatrixXd::Zero(2, 5);//先存到矩阵中
				//for (int j = 0; j < 5; j++)
				//{
				//	box_matrix(0, j) = plist_box[j].x;
				//	box_matrix(1, j) = plist_box[j].y;
				//}
				//MatrixXd box_matrix_o = eigenvec_double*box_matrix;//再左乘 特征向量矩阵
				//for (int j = 0; j < 5; j++)
				//{
				//	plist_box[j].x = box_matrix_o(0, j);
				//	plist_box[j].y = box_matrix_o(1, j);
				//}
				//polygon_vec.push_back(plist_box);
				////原点转换到 原来的坐标系
				//Vector2d revol_center_vec;
				//revol_center_vec << revol_center_pt.x, revol_center_pt.y;
				//Vector2d center_vec1 = eigenvec_double*revol_center_vec;	
				//mycircle acir;
				//acir.cx = center_vec1[0];
				//acir.cy = center_vec1[1];
				//acir.radius = 1;
				//circle_vec.push_back(acir);
						
			/*3、面积大小（归一化的）*/
				double are_ratio = plrGraph.patches[i].area;//(plrGraph.patches[i].area - min_area) / (max_area - min_area);

			vector<double> fvec;//特征向量
			fvec.push_back(circularity);//圆度
			fvec.push_back(elongation);//长度
			fvec.push_back(are_ratio);//面积比
			//fvec.push_back(terp_num);
			cir_vec.push_back(circularity);
			long_vec.push_back(elongation);
			area_vec.push_back(are_ratio);
			patchfeatures.push_back(make_pair(i, fvec));
		}
	}
	//输出特征值
	ofstream outfile1("E:/circularity_output.txt");
	ofstream outfile2("E:/elongation_output.txt");
	ofstream outfile3("E:/area_output.txt");
	if (outfile1.is_open() && outfile2.is_open()&&outfile3.is_open())
	{
		for (unsigned int i = 0; i < area_vec.size(); ++i)
		{
			//if (fabs(cir_vec[i] - long_vec[i])<=0.4)
			//	plrGraph.patches[patchfeatures[i].first].type = 1;//red
			//else
			//	plrGraph.patches[patchfeatures[i].first].type = 2;//blue
			outfile1 << cir_vec[i];
			outfile1 << " ";
			outfile1 << "\n";
			outfile2 << long_vec[i];
			outfile2 << " ";
			outfile2 << "\n";
			outfile3 << area_vec[i];
			outfile3 << " ";
			outfile3 << "\n";
		}
	}
	outfile1.close(); outfile2.close(); outfile3.close();
	cout << "get feature done" << endl;
}
void MyGraph::KmeansClustering(int classnum)//指定分类数量
{
	cout << "KmeansClustering" << endl;
	int ft_dimention = patchfeatures[0].second.size();//ft_dimension 是特征向量的维度
	double maxlong = -1; double minlong = 1; double maxArea = 0;
	 //int init_center1;
	 //int init_center2;
	 //int init_center3;
	//for (unsigned int i = 0; i < patchfeatures.size(); ++i)
	//{
	//	if (patchfeatures[i].second.at(1)>maxlong)
	//	{
	//		init_center1 =i;
	//		maxlong = patchfeatures[i].second.at(1);
	//	}
	//	if (patchfeatures[i].second[1] < minlong)
	//	{
	//		init_center2 = i;
	//		minlong = patchfeatures[i].second.at(1);
	//	}
	//	
	//}	

	//产生patch index的随机数，即随机产生类中心
	vector<vector<double>> centroids;
	srand((unsigned)time(NULL));
	for (int i = 0; i < classnum; i++)
	{
		int rnum = rand() % (patchfeatures.size() - 1);
		cout << rnum << ' ';
		centroids.push_back(patchfeatures[rnum].second);
	}
	/*centroids.push_back(patchfeatures[init_center1].second);
	centroids.push_back(patchfeatures[init_center2].second);
	centroids.push_back(patchfeatures[init_center3].second);*/
	bool clusterChanged = true;
	int patchsize = patchfeatures.size();
	//初始化一个距离向量，在聚类完之后便于观察分析距离
	vector<pair<int, double>> type_dis;//int 对应第几类（不是patchindex）， ddouble 对应到中心点的距离
	for (int i = 0; i < patchsize; ++i)
	{
		type_dis.push_back(make_pair(0, 0));
	}
	//开始k-means聚类
	while (clusterChanged)
	{
		clusterChanged = false;
		//step one : find the nearest centroid of each point  
		for (int i = 0; i<patchsize; i++)
		{
			int typeIndex = 0;
			double minDist = INT_MAX;
			for (int j = 0; j< classnum; j++)
			{
				double dist = getEuclidDistance(&centroids[j], &patchfeatures[i].second);
				if (dist < minDist )
				{
					minDist = dist;
					typeIndex = j + 1;//类别数字是1，2,3（不是0，1,2）
				}
			}
			if (plrGraph.patches[patchfeatures[i].first].type != typeIndex)//update the cluster which the dataSet[i] belongs to... 
			{
				clusterChanged = true;
				plrGraph.patches[patchfeatures[i].first].type = typeIndex;
				type_dis[i].first = typeIndex; type_dis[i].second = minDist;//	clusterAssment[i].minDist = minDist;
			}
		}
		//step two : update the centroids  
		for (int type_i = 0; type_i<classnum; type_i++)
		{
			vector<double> a_centroid(ft_dimention, 0);
			int cnt = 0;
			int typeindex = type_i + 1;
			for (int i = 0; i<patchsize; i++)
			{
				if (plrGraph.patches[patchfeatures[i].first].type == typeindex)//if (type_dis[i].first == typeindex)
				{
					++cnt;
					//sum of two vectors  
					for (int j = 0; j<ft_dimention; j++)//ft_dimention是特征向量的维度
					{
						a_centroid[j] += patchfeatures[i].second.at(j);
					}
				}
			}
			cout << "class " << typeindex << "'s count=" << cnt << endl;
			//mean of the vector and update the centroids[i]  
			for (int i = 0; i<ft_dimention; i++)
			{
				if (cnt != 0)
					a_centroid[i] = a_centroid[i] / cnt;
				centroids[type_i].at(i) = a_centroid[i];
			}
		}
		cout << "cycle---" << endl;
	}
	//计算两个类 的类内距离 的平均值，标准差，
	double mean_t1 = 0;	double mean_t2 = 0;double mean_area=0; double var_t1, var_t2,var_area;
	vector<double> vec_type1, vec_type2,vec_area;
	for (unsigned int i = 0; i < type_dis.size(); ++i)
	{
		if (type_dis[i].first == 1)
		{
			vec_type1.push_back(type_dis[i].second);
		}
		else if (type_dis[i].first == 2)
		{
			vec_type2.push_back(type_dis[i].second);
		}
		//vec_area.push_back(patchfeatures[i].second[2]);
	}
	mean_t1 = GetAverange(vec_type1);//均值
	mean_t2 = GetAverange(vec_type2);
	mean_area = GetAverange(vec_area);
	var_t1 = sqrt(GetVariance(vec_type1, mean_t1, vec_type1, mean_t1));//标准差
	var_t2 = sqrt(GetVariance(vec_type2, mean_t2, vec_type2, mean_t2));
	var_area = sqrt(GetVariance(vec_area, mean_area, vec_area, mean_area));//通过area来删除噪声点？？
	//删除均值较大的那一组patch（比较分散的）
	int dele_type = 0;//比较松散的那一类
	if (mean_t1 > mean_t2)
		dele_type = 1;
	else
		dele_type = 2;
	//for (unsigned int i = 0; i < type_dis.size(); ++i)
	//{
	//	if (type_dis[i].first != dele_type)
	//	{
	//		vec_area.push_back(patchfeatures[i].second[2]);
	//	}
	//}
	//mean_area = GetAverange(vec_area);
	//var_area = sqrt(GetVariance(vec_area, mean_area, vec_area, mean_area));//通过area来删除噪声点？？
	for ( int i = 0; i < patchsize; ++i)
	{
		////超出 均值+-标准差 的范围，标记为噪声点
		//if (type_dis[i].first == 1)
		//{
		//	if ( type_dis[i].second>(mean_t1 + var_t1))
		//	{
		//		plrGraph.patches[patchfeatures[i].first].type = 0;//超出 均值+-标准差 的范围，标记为噪声点
		//	}
		//}
		//else if (type_dis[i].first == 2)
		//{
		//	if (type_dis[i].second>(mean_t2 + var_t2))
		//	{
		//		plrGraph.patches[patchfeatures[i].first].type = 0;//超出 均值+-标准差 的范围，标记为噪声点
		//	}
		//}

		//删除均值较大的那一组patch（比较分散的）
		if (type_dis[i].first == dele_type)
		{
			plrGraph.patches[patchfeatures[i].first].type = 0;
		}
		else//剩下要保留的patch里再删除面积较大的（超出 均值+标准差的）
		{
			if (patchfeatures[i].second[2] > (mean_area + var_area))
			{
				plrGraph.patches[patchfeatures[i].first].type = 0;
			}
		}
	}
	////二次分类，对cross patch 进行二次分类，里面还惨杂着些大的patch,结果不理想！！------------------------------------------
	//int recluster_type = 0;//需要重新聚类的 类别
	//if (mean_t1 > mean_t2)
	//	recluster_type = 1;
	//else
	//	recluster_type = 2;
	//centroids.clear(); //随机类中心 vector<vector<double>> centroids
	////样本数量比较少，随机选取中心可能结果不好，选择面积最大和最小的两个样本作为中心  double maxlong = -1; double minlong = 1; double maxArea = 0;
	//int init_center1;
	//int init_center2;
	//for (unsigned int i = 0; i < patchfeatures.size(); ++i)
	//{
	//	if (type_dis[i].first == recluster_type)
	//	{
	//		if (patchfeatures[i].second.at(2)>maxlong)
	//		{
	//			init_center1 = i;
	//			maxlong = patchfeatures[i].second.at(2);
	//		}
	//		if (patchfeatures[i].second[2] < minlong)
	//		{
	//			init_center2 = i;
	//			minlong = patchfeatures[i].second.at(2);
	//		}
	//	}	
	//}	
	//centroids.push_back(patchfeatures[init_center1].second);
	//centroids.push_back(patchfeatures[init_center2].second);
	//clusterChanged = true; //k-means循环的标志位
	//vector<pair<int, double>> type_dis_recluster;//int 对应第几类（不是patchindex）， double 对应到中心点的距离
	//for (int i = 0; i < patchsize; ++i)
	//{
	//	type_dis_recluster.push_back(make_pair(0, 0));
	//}
	//
	//while (clusterChanged)
	//{
	//	clusterChanged = false;
	//	//step one : find the nearest centroid of each point  
	//	for (int i = 0; i<patchsize; i++)
	//	{
	//		if (type_dis[i].first == recluster_type)
	//		{
	//			int typeIndex = 0;
	//			double minDist = INT_MAX;
	//			for (int j = 0; j< classnum; j++)
	//			{
	//				double dist = getEuclidDistance(centroids[j], patchfeatures[i].second);
	//				if (dist < minDist)
	//				{
	//					minDist = dist;
	//					typeIndex = j + 3;//0,1,2已经用过了，从3，4开始
	//				}
	//			}
	//			if (plrGraph.patches[patchfeatures[i].first].type != typeIndex)//update the cluster which the dataSet[i] belongs to... 
	//			{
	//				clusterChanged = true;
	//				plrGraph.patches[patchfeatures[i].first].type = typeIndex;
	//				type_dis_recluster[i].first = typeIndex; type_dis_recluster[i].second = minDist;//	clusterAssment[i].minDist = minDist;
	//			}
	//		}		
	//	}
	//	//step two : update the centroids  
	//	for (int type_i = 0; type_i<classnum; type_i++)
	//	{
	//		vector<double> a_centroid(ft_dimention, 0);
	//		int cnt = 0;
	//		int typeindex = type_i + 3;//0,1,2已经用过了，从3，4开始
	//		for (int i = 0; i<patchsize; i++)
	//		{
	//			if (plrGraph.patches[patchfeatures[i].first].type == typeindex)//if (type_dis[i].first == typeindex)
	//			{
	//				++cnt;
	//				//sum of two vectors  
	//				for (int j = 0; j<ft_dimention; j++)//ft_dimention是特征向量的维度
	//				{
	//					a_centroid[j] += patchfeatures[i].second.at(j);
	//				}
	//			}
	//		}
	//		cout << "class " << typeindex << "'s count=" << cnt << endl;
	//		//mean of the vector and update the centroids[i]  
	//		for (int i = 0; i<ft_dimention; i++)
	//		{
	//			if (cnt != 0)
	//				a_centroid[i] = a_centroid[i] / cnt;
	//			centroids[type_i].at(i) = a_centroid[i];
	//		}
	//	}
	//	cout << "cycle---" << endl;
	//}

	//输出k-mean中每个类的类间距离----------------------------------------------------------------
	ofstream outfile1("E:/distance_type1.txt");
	ofstream outfile2("E:/distance_type2.txt");
	if (outfile1.is_open() && outfile2.is_open())
	{
		for (unsigned int i = 0; i < type_dis.size(); ++i)
		{
			if (type_dis[i].first == 1)
			{
				outfile1 << type_dis[i].second;
				outfile1 << " ";
				outfile1 << "\n";
			}
			else if(type_dis[i].first == 2)
			{
				outfile2 << type_dis[i].second;
				outfile2 << " ";
				outfile2 << "\n";
			}
		}
	}
	outfile1.close();
	outfile2.close();
	
	/*
	K-MEDODIS的具体流程如下：
1）任意选取K个对象作为medoids（O1,O2,…Oi…Ok）。　　
2）将余下的对象分到各个类中去（根据与medoid最相近的原则）；　　
3）对于每个类（Oi）中，顺序选取一个Or，计算用Or代替Oi后的消耗—E（Or）。选择E最小的那个Or来代替Oi。这样K个medoids就改变了。*/
 	cout << "k-mean done" << endl;
}
void MyGraph::KmedoidClustering(int classnum)
{
	cout << "KmeansClustering" << endl;
	int ft_dimention = patchfeatures[0].second.size();//ft_dimension 是特征向量的维度
	double maxlong = -1; double minlong = 1; double maxArea = 0;
	//产生patch index的随机数，即随机产生类中心
	vector<vector<double>> centroids;
	srand((unsigned)time(NULL));
	for (int i = 0; i < classnum; i++)
	{
		int rnum = rand() % (patchfeatures.size() - 1);
		cout << rnum << ' ';
		centroids.push_back(patchfeatures[rnum].second);
	}
	bool clusterChanged = true;
	int patchsize = patchfeatures.size();
	//初始化一个距离向量，在聚类完之后便于观察分析距离
	vector<pair<int, double>> type_dis;//int 对应第几类（不是patchindex）， ddouble 对应到中心点的距离
	for (int i = 0; i < patchsize; ++i)
	{
		type_dis.push_back(make_pair(0, 0));
	}
	//开始k-means聚类
	while (clusterChanged)
	{
		clusterChanged = false;
		//step one : find the nearest centroid of each point  
		for (int i = 0; i<patchsize; i++)
		{
			int typeIndex = 0;
			double minDist = INT_MAX;
			for (int j = 0; j< classnum; j++)
			{
				double dist = getEuclidDistance(&centroids[j], &patchfeatures[i].second);
				if (dist < minDist)
				{
					minDist = dist;
					typeIndex = j + 1;//类别数字是1，2,3（不是0，1,2）
				}
			}
			if (plrGraph.patches[patchfeatures[i].first].type != typeIndex)//update the cluster which the dataSet[i] belongs to... 
			{
				clusterChanged = true;
				plrGraph.patches[patchfeatures[i].first].type = typeIndex;
				type_dis[i].first = typeIndex; type_dis[i].second = minDist;//	clusterAssment[i].minDist = minDist;
			}
		}
		//step two : update the centroids
		vector<double> nei_count(patchsize,0);
		for (int type_i = 0; type_i<classnum; type_i++)
		{
			int typeindex = type_i + 1;
			for (int i = 0; i<patchsize; i++)
			{
				if (plrGraph.patches[patchfeatures[i].first].type == typeindex)//if (type_dis[i].first == typeindex)
				{
					for (int j = i+1; j < patchsize; ++j)
					{
						if (plrGraph.patches[patchfeatures[j].first].type == typeindex)
						{
							double dist = getEuclidDistance(&patchfeatures[i].second, &patchfeatures[j].second);
							nei_count[i] = nei_count[i] + dist;
							nei_count[j] = nei_count[j] + dist;
						}
					}	
				}
			}			
			// update the centroids[i]  
			vector<double> a_centroid(ft_dimention, 0);
			double mindis = INT_MAX;
			for (int i = 0; i<patchsize; i++)
			{
				if (plrGraph.patches[patchfeatures[i].first].type == typeindex)
				{
					if (nei_count[i] < mindis)
					{
						a_centroid = patchfeatures[i].second;
						mindis = nei_count[i];
					}
				}
			}
			centroids[type_i] = a_centroid;
		}
		cout << "cycle---" << endl;
	}
	//计算两个类 的类内距离 的平均值，标准差，删除超出 均值+-标准差 的范围的patch
	double mean_t1 = 0;	double mean_t2 = 0; double var_t1, var_t2;
	vector<double> vec_type1, vec_type2;
	for (unsigned int i = 0; i < type_dis.size(); ++i)
	{
		if (type_dis[i].first == 1)
		{
			vec_type1.push_back(type_dis[i].second);
		}
		else if (type_dis[i].first == 2)
		{
			vec_type2.push_back(type_dis[i].second);
		}
	}
	mean_t1 = GetAverange(vec_type1);//均值
	mean_t2 = GetAverange(vec_type2);
	var_t1 = sqrt(GetVariance(vec_type1, mean_t1, vec_type1, mean_t1));//标准差
	var_t2 = sqrt(GetVariance(vec_type2, mean_t2, vec_type2, mean_t2));
	for (int i = 0; i < patchsize; ++i)
	{
		if (type_dis[i].first == 1)
		{
			if (type_dis[i].second<(mean_t1 - var_t1) || type_dis[i].second>(mean_t1 + var_t1))
			{
				plrGraph.patches[patchfeatures[i].first].type = 0;//超出 均值+-标准差 的范围，标记为噪声点
			}
		}
		else if (type_dis[i].first == 2)
		{
			if (type_dis[i].second<(mean_t2 - var_t2) || type_dis[i].second>(mean_t2 + var_t2))
			{
				plrGraph.patches[patchfeatures[i].first].type = 0;//超出 均值+-标准差 的范围，标记为噪声点
			}
		}
	}
	//输出k-mean中每个类的类间距离----------------------------------------------------------------
	ofstream outfile1("E:/distance_type1.txt");
	ofstream outfile2("E:/distance_type2.txt");
	if (outfile1.is_open() && outfile2.is_open())
	{
		for (unsigned int i = 0; i < type_dis.size(); ++i)
		{
			if (type_dis[i].first == 1)
			{
				outfile1 << type_dis[i].second;
				outfile1 << " ";
				outfile1 << "\n";
			}
			else if (type_dis[i].first == 2)
			{
				outfile2 << type_dis[i].second;
				outfile2 << " ";
				outfile2 << "\n";
			}
		}
	}
	outfile1.close();
	outfile2.close();

	/*
	K-MEDODIS的具体流程如下：
	1）任意选取K个对象作为medoids（O1,O2,…Oi…Ok）。　　
	2）将余下的对象分到各个类中去（根据与medoid最相近的原则）；　　
	3）对于每个类（Oi）中，顺序选取一个Or，计算用Or代替Oi后的消耗—E（Or）。选择E最小的那个Or来代替Oi。这样K个medoids就改变了。*/
	cout << "k-mean done" << endl;
}
void MyGraph::DBSCANClusering()//不可行，在交叉点的patch数量本来就少，规定radius半径之内的点的个数，若radius大了找不到crosspatch
{
	//初始化数据集dataset 
	int pchsize = patchfeatures.size();
	vector<DataPoint> patchDataSet;//patchDataSet的大小和patchfeature的大小相同（小于plrgraph.patch）
	DataPoint dd;
	for (int i = 0; i < pchsize; ++i)
	{
		dd.dataID = i;// dataID不等于patchfeatures[i].first，即不等于patchindex
		dd.pch_ft = patchfeatures[i].second;
		dd.clusterId = 0;
		dd.datatype = 1;//1 noise 2 border 3 core 
		dd.visited = false;
		//dd.neigbor_ct = 0;
		patchDataSet.push_back(dd);
	}
	int cluster_count=dbscanclustering(patchDataSet, 0.35, 3);

	cout << "cluster number=" << cluster_count << endl;
	int noisy_cnt = 0;
	int type1_count = 0;
	int type2_count = 0;
	for (int i = 0; i < pchsize; ++i)
	{
		if (patchDataSet[i].datatype == 1)//噪音点
		{
			noisy_cnt++;
			int pch_index = patchfeatures[i].first;
			plrGraph.patches[pch_index].type = 3;//绿色
		}
		else
		{
			if (patchDataSet[i].clusterId == 1)
			{
				int pch_index = patchfeatures[i].first;
				plrGraph.patches[pch_index].type = 1;
				type1_count++;
			}
			else if (patchDataSet[i].clusterId == 2)
			{
				int pch_index = patchfeatures[i].first;
				plrGraph.patches[pch_index].type = 2;
				type2_count++;
			}
		}		
	}
	cout << "noisy count=" << noisy_cnt << ",type1_count=" << type1_count << ",type2_count=" << type2_count << endl;
	cout << "dbscan done" << endl;
}
void MyGraph::FCMClustering(int classnum, double para_m)
{
	int dimension = patchfeatures[0].second.size();
	vector<vector<double>> U;//隶属度矩阵，设有100个数数据，2个分类，则是2*100的矩阵
	vector<vector<double>> u_tmp;//U的每个元素的para_m次方的矩阵；
	vector<vector<double>> centers;//类中心
	//初始化隶属度矩阵，列数是待分析数据的个数，行数是要分的类的个数。需满足每一列的和是1
	srand((unsigned)time(NULL));
	for (int i = 0; i < classnum; ++i)//随机产生小数
	{
		vector<double> tmp_row;
		for (int j = 0; j < patchfeatures.size(); j++)
		{
			double rnum = (double)((rand() % 999) + 1) / 1000;
			tmp_row.push_back(rnum);
		}
		U.push_back(tmp_row);
		u_tmp.push_back(tmp_row);
	}
	//使每一列和为1
	for (int i = 0; i < patchfeatures.size(); i++)
	{
		double sum = 0;
		for (int j = 0; j < classnum; ++j)
		{
			sum = sum + U[j][i];
		}
		for (int j = 0; j < classnum; ++j)
		{
			U[j][i] = U[j][i] / sum;
		}
	}	
	//开始迭代
	for (int i = 0; i < 10; ++i)//循环次数20
	{
		//更新centers
		for (unsigned int j = 0; j < classnum; ++j)
		{
			vector<double> a_ct(dimension, 0);//中心点（是个向量），初始化为维度是dimension的零向量
			double u_sum_j = 0;//第j行的和
			for (unsigned int h = 0; h < patchfeatures.size(); ++h)
			{
				double tmp = pow(U[j][h], para_m);
				u_tmp[j][h] = tmp;//对u_tmp幅值(更新)
				u_sum_j = u_sum_j + tmp;
				for (int k = 0; k < dimension; k++)
				{
					a_ct[k] = a_ct[k] + tmp*patchfeatures[h].second[k];
				}
			}
			for (int k = 0; k < dimension; k++)
			{
				a_ct[k] = a_ct[k] / u_sum_j;
			}
			centers.push_back(a_ct);
		}
		//计算目标函数J
		double J_result = 0;
		for (int j = 0; j < classnum; ++j)
		{
			for (unsigned int h = 0; h < patchfeatures.size(); ++h)
			{
				J_result = J_result + u_tmp[j][h] * pow(getEuclidDistance(&patchfeatures[h].second, &centers[j]),2);
			}
		}
		cout << "J_result=" << J_result << endl;
		//更新U
		for (int j = 0; j < classnum; ++j)
		{
			for (unsigned int h = 0; h <patchfeatures.size(); ++h)
			{
				double sum_uh=0;//计算U 第h列的和
				double dis1 = getEuclidDistance(&patchfeatures[h].second, &centers[j]);
				for (int k = 0; k < classnum; ++k)
				{
					double dis2 = getEuclidDistance(&patchfeatures[h].second, &centers[k]);
					sum_uh = sum_uh + pow((dis1 / dis2), (2 / (para_m - 1)));
				}
				U[j][h] = 1/sum_uh;
			}
		}
	}
	//选择隶属度较大的那个
	for (unsigned int i = 0; i < patchfeatures.size(); ++i)
	{
		double max_u = 0;
		int type_index = 0;
		for (int j = 0; j < classnum; ++j)
		{
			if (U[j][i]>max_u)
			{
				max_u = U[j][i];
				type_index = j;
			}
		}
		plrGraph.patches[patchfeatures[i].first].type = type_index+1;

			//double max_u = max(U[0][i], U[1][i]);
			//int type_index = 0;
			//if (fabs(U[0][i] - U[1][i])>0.09)//隶属度之差>0.1
			//{
			//	if (U[0][i]==max_u)
			//	{
			//		type_index = 0;
			//	}
			//	else
			//	{
			//		type_index = 1;
			//	}
			//	plrGraph.patches[patchfeatures[i].first].type = type_index + 1;
			//}
			//else
			//	plrGraph.patches[patchfeatures[i].first].type = 0;//隶属度之差太小，可能是噪声点？？？
	}	
	//输出隶属度
	ofstream outfile1("E:/u_i.txt");
	ofstream outfile2("E:/u_j.txt");
	if (outfile1.is_open() && outfile2.is_open())
	{
		for (unsigned int i = 0; i < patchfeatures.size(); ++i)
		{
			if (U[0][i]>U[1][i])
			{
				outfile1 << U[0].at(i);
				outfile1 << " ";
				outfile1 << "\n";
			}
			else
			{
				outfile2 << U[1].at(i);
				outfile2 << " ";
				outfile2 << "\n";
			}		
		}
	}
	outfile1.close();
	outfile2.close();
	////输出隶属度
	//ofstream outfile1("E:/u_i.txt");
	//ofstream outfile2("E:/u_j.txt");
	//if (outfile1.is_open() && outfile2.is_open())
	//{
	//	for (unsigned int i = 0; i < patchfeatures.size(); ++i)
	//	{
	//			outfile1 << U[0].at(i);
	//			outfile1 << " ";
	//			outfile1 << "\n";
	//			outfile2 << U[1].at(i);
	//			outfile2 << " ";
	//			outfile2 << "\n";
	//	}
	//}
	//outfile1.close();
	//outfile2.close();

	cout << "FCMClustering done" << endl;
}

//得到相邻的线（角度，距离接近）
void MyGraph::getOrigiCurveList(char currt_type, int currt_index, char& origin_type, int& origin_index, vector<mypoint2f> &plist, int flag)
{
	if (currt_type == 'L')
	{
		origin_type = polyline_vec[currt_index].origi_type;
		origin_index = polyline_vec[currt_index].origi_index;
	}
	else if (currt_type == 'C')
	{
		origin_type = 'C';
		origin_index = currt_index;
	}
	else if (currt_type == 'Q')
	{
		origin_type = 'Q';
		origin_index = currt_index;
	}
	if (flag == 0)
		return;
	char ttt = origin_type;
	switch (ttt)
	{
	case 'C':
		//cubicBezier cb1=cubicBezier_vec[origin_index];
		getPListfromCubic(0, cubicBezier_vec[origin_index], 0, &plist);
		break;
	case 'c':
		//cubicBezier cb2=cubicBezier_vec_origin[origin_index];
		getPListfromCubic(0, cubicBezier_vec_origin[origin_index], 0, &plist);
		break;
	case 'Q':
		//quaBezier qb1=quadraticBezier_vec[origin_index];
		getPListfromQuadratic(0, quadraticBezier_vec[origin_index], 0, &plist);
		break;
	case 'q':
		//quaBezier qb2 = quadraticBezier_vec_origin[origin_index];
		getPListfromQuadratic(0, quadraticBezier_vec_origin[origin_index], 0, &plist);
		break;
	case 'L':
		plist = polyline_vec[origin_index].pointlist;
		break;
	case 'l':
		plist = polyline_vec_origin[origin_index].pointlist;
		break;
	}
}
vector<EdgeStruct> MyGraph::getNeiborLines(int nodeindex1, int nodeindex2, int scope, double width_threshold, double angle_threshold)
{
	////方法一，找line segment ，
	//vector<EdgeStruct> nei_diffCurve;
	////起点终点坐标， pointlist
	//mypoint2f node_position1 = sgraph.nodeList[nodeindex1].position;
	//mypoint2f node_position2 = sgraph.nodeList[nodeindex2].position;
	//vector<mypoint2f> edge_plist;
	//getPListOfaCurve(sgraph.nodeList[nodeindex1].position, sgraph.nodeList[nodeindex2].position, nodeindex1, 0,&edge_plist);
	////平均斜率  nodeindex1-->nodeindex2 线条的斜率
	//vector<double> aver_tan_vec;
	//for (unsigned int i = 1; i < edge_plist.size()-1; ++i)
	//{
	//	mypoint2f p1 = edge_plist[i] - edge_plist[i - 1];
	//	mypoint2f p2 = edge_plist[i + 1] - edge_plist[i];
	//	double tmp_atan1 = getAtan(p1);
	//	double tmp_atan2 = getAtan(p2);
	//	double atn = (tmp_atan1 + tmp_atan2) / 2;
	//	aver_tan_vec.push_back(atn);
	//}
	//double aver_atan=GetAverange(aver_tan_vec);
	////从头至尾 总体斜率
	//double rough_atan = getAtan(sgraph.nodeList[nodeindex2].position - sgraph.nodeList[nodeindex1].position);
	////广度遍历  ,距离起始点最远距离为scope
	//int* visited = new int[sgraph.vertexNum];
	//for (int i = 0; i < sgraph.vertexNum; ++i)
	//{
	//	visited[i] = 0;
	//}
	//queue<int> st;
	//st.push(nodeindex1);
	//visited[nodeindex1] = 1;
	//while (!st.empty())
	//{
	//	int tmpnode = st.front();
	//	st.pop();
	//	visited[tmpnode] = 2;
	//	ArcNode *arc = sgraph.nodeList[tmpnode].firstEdge;
	//	while (arc)
	//	{
	//		if (visited[arc->adjVertex] != 2)
	//		{
	//			if (visited[arc->adjVertex] == 0)
	//			{
	//				if (LineSegmentLength(node_position1, sgraph.nodeList[arc->adjVertex].position)<=scope)
	//				{
	//					st.push(arc->adjVertex);
	//					visited[arc->adjVertex] = 1;
	//				}
	//				else if (LineSegmentLength(node_position2, sgraph.nodeList[arc->adjVertex].position) <= scope)
	//				{
	//					st.push(arc->adjVertex);
	//					visited[arc->adjVertex] = 1;
	//				}
	//			}
	//			char ctype = arc->curveType;
	//			int cindex = arc->curveIndex;
	//			vector<mypoint2f> tmp_plist;
	//			getPListOfaCurve(sgraph.nodeList[tmpnode].position, sgraph.nodeList[arc->adjVertex].position, tmpnode, 0, &tmp_plist);
	//			int pindex1, pindex2;
	//			double proxi = getProximity(edge_plist, tmp_plist, pindex1, pindex2);
	//			double conti = getContinuity(edge_plist, tmp_plist, pindex1, pindex2);
	//			if (proxi < 10 && conti < 0.1)// 5cm以内  25.85度以下
	//			{
	//				EdgeStruct est;
	//				est.curveType = arc->curveType;
	//				est.curveIndex = arc->curveIndex;
	//				nei_diffCurve.push_back(est);
	//				PolyLine apl;
	//				apl.pointlist = tmp_plist;
	//				plist_group.push_back(apl);  //plist_group .h文件中  测试用
	//			}
	//		}
	//		arc = arc->next;
	//	}
	//}
	//delete[] visited;
	//visited = NULL;
	//return nei_diffCurve;

	//方法二 找原来输入的线origin line
	vector<EdgeStruct> nei_diffCurve;//存放该短边对应的邻域边 的origi_curve
	char curve_type = ' '; char origi_type = ' ';
	int curve_index = 0;  int origi_index = 0;
	ArcNode* anode = sgraph.nodeList[nodeindex1].firstEdge;
	while (anode)
	{
		if (anode->adjVertex == nodeindex2)
		{
			curve_type = anode->curveType;
			curve_index = anode->curveIndex;
			break;
		}
		anode = anode->next;
	}
	vector<mypoint2f> edge_plist;
	getOrigiCurveList(curve_type, curve_index, origi_type, origi_index, edge_plist,1);//edge_plist是origincurve的点列
	PolyLine apl;
	apl.pointlist = edge_plist;
	neighbor_plist_group.push_back(apl);  //在.h文件中 用于测试显示 
	EdgeStruct anedge;
	anedge.curveIndex = origi_index;
	anedge.curveType = origi_type;
	nei_diffCurve.push_back(anedge);//先把node1-->node2的origicurve放进去

	mypoint2f node_position1 = sgraph.nodeList[nodeindex1].position;
	mypoint2f node_position2 = sgraph.nodeList[nodeindex2].position;
	//getPListOfaCurve(sgraph.nodeList[nodeindex1].position, sgraph.nodeList[nodeindex2].position, nodeindex1, 0, &edge_plist);
	//广度遍历  ,距离起始点最远距离为scope
	int* visited = new int[sgraph.vertexNum];
	for (int i = 0; i < sgraph.vertexNum; ++i)
	{
		visited[i] = 0;
	}
	queue<int> st;
	st.push(nodeindex1);
	visited[nodeindex1] = 1;
	while (!st.empty())
	{
		int tmpnode = st.front();
		st.pop();
		visited[tmpnode] = 2;
		ArcNode *arc = sgraph.nodeList[tmpnode].firstEdge;
		while (arc)
		{
			if (visited[arc->adjVertex] != 2)
			{
				if (visited[arc->adjVertex] == 0)
				{
					if (LineSegmentLength(node_position1, sgraph.nodeList[arc->adjVertex].position) <= scope)
					{
						st.push(arc->adjVertex);
						visited[arc->adjVertex] = 1;
					}
					else if (LineSegmentLength(node_position2, sgraph.nodeList[arc->adjVertex].position) <= scope)
					{
						st.push(arc->adjVertex);
						visited[arc->adjVertex] = 1;
					}
				}
				char ctype = arc->curveType;   char otype = ' ';
				int cindex = arc->curveIndex;  int oindex = 0;
				vector<mypoint2f> tmp_plist;
				getOrigiCurveList(ctype, cindex, otype, oindex, tmp_plist,1);
				//getPListOfaCurve(sgraph.nodeList[tmpnode].position, sgraph.nodeList[arc->adjVertex].position, tmpnode, 0, &tmp_plist);
				int pindex1, pindex2;
				double proxi = getProximity(&edge_plist, &tmp_plist, pindex1, pindex2);//原origicurve 和 neighbor的origicurve 相比较：最近距离
				double conti = getContinuity(&edge_plist, &tmp_plist, pindex1, pindex2);//最近距离处的夹角
				if (proxi < width_threshold && conti < angle_threshold)//最近距离可以设置的较大一些 if (proxi < 5 && conti < 0.1)// 5cm以内  25.85度以下
				{
					EdgeStruct est;
					est.curveType = otype;
					est.curveIndex = oindex;
					nei_diffCurve.push_back(est);

					PolyLine apl;
					apl.pointlist = tmp_plist;
					neighbor_plist_group.push_back(apl);  //neighbor_plist_group  在.h文件中  测试用
				}
			}
			arc = arc->next;
		}
	}
	delete[] visited;
	visited = NULL;
	//删除重复的curve
	vector<EdgeStruct>::iterator iter_edge = nei_diffCurve.begin();
	vector<PolyLine>::iterator iter_polyline = neighbor_plist_group.begin();
	int index_out = 0;
	int index_inner = 0;
	while (iter_edge != nei_diffCurve.end())
	{
		vector<EdgeStruct>::iterator edge_inner = iter_edge + 1;
		index_inner = index_out + 1;
		while (edge_inner != nei_diffCurve.end())
		{
			if (iter_edge->curveType == edge_inner->curveType && iter_edge->curveIndex == edge_inner->curveIndex)
			{
				edge_inner = nei_diffCurve.erase(edge_inner);
				neighbor_plist_group.erase(iter_polyline + index_inner);
				iter_polyline = neighbor_plist_group.begin();
			}
			else
			{
				edge_inner++;
				index_inner++;
			}
		}
		iter_edge++;
		index_out++;
	}
	return nei_diffCurve;
}

//腐蚀
void MyGraph::Erosion_test()
{
	//对边界patch进行腐蚀，注意判断拐点数量2个
	//EdgeOfPatch * getBorderofPatch()获取patch的一条边border，  vector<int> getNeighborPatches()得到与某patch相邻的所有patch
	int pchsize = plrGraph.patches.size();

	//vector<int> edgepatch;
	//for (int i = 0; i < pchsize; ++i)
	//{
	//	bool eflag = false;
	//	for (unsigned int k = 0; k < plrGraph.patches[i].edges_vec.size(); ++k)
	//	{
	//		if (plrGraph.patches[i].edges_vec[k].adjPatchIndex.size() == 1)
	//		{
	//			eflag = true; break;
	//		}
	//	}
	//	if (plrGraph.patches[i].type >= 0 && eflag == true)//可用 且是边界patch
	//	{
	//		edgepatch.push_back(i);
	//	}
	//}
	//for (int i = 0; i <edgepatch.size(); ++i)
	//{
	//	if (plrGraph.patches[edgepatch[i]].type >= 0)//可用 且是边界patch
	//	{
	//		//getpatch pointlist,判断拐点
	//		vector<mypoint2f> apatch;
	//		for (unsigned int j = 1; j < plrGraph.patches[edgepatch[i]].pointIndexOfPatch.size(); ++j)
	//		{
	//			int verindex1 = plrGraph.patches[edgepatch[i]].pointIndexOfPatch.at(j - 1);
	//			int verindex2 = plrGraph.patches[edgepatch[i]].pointIndexOfPatch.at(j);
	//			mypoint2f point1 = sgraph.nodeList[verindex1].position;//=lineWithCommonPoints[verindex1].first
	//			mypoint2f point2 = sgraph.nodeList[verindex2].position;
	//			if (!apatch.empty())
	//				apatch.pop_back();
	//			getPListOfaCurve(point1, point2, verindex1, verindex2, &apatch);
	//		}
	//		if (isPointRoughEqual(apatch.front(), apatch.back()))
	//			apatch.pop_back();
	//		mypoint2f pfront = apatch.front();
	//		mypoint2f pback = apatch.back();
	//		apatch.insert(apatch.begin(), pback);
	//		apatch.push_back(pfront);
	//		int turn_count = 0;
	//		for (unsigned int j = 1; j < apatch.size() - 1; ++j)
	//		{
	//			mypoint2f vec1 = apatch[j] - apatch[j - 1];
	//			mypoint2f vec2 = apatch[j + 1] - apatch[j];
	//			double cos_theta = getCosineAngle(apatch[j - 1], apatch[j], apatch[j], apatch[j + 1]);
	//			if (cos_theta >= -1 && cos_theta <= 0.707106)//   135-->-0.7071  150--> -0.866
	//			{
	//				turn_count++; //cout << cos_theta << endl;
	//			}
	//		}
	//		if (turn_count < 4)//拐点小于4个才可以 腐蚀成一条线？？
	//		{
	//			//找最远的两个点,
	//			int p1, p2; double fdis = 0;
	//			vector<int> tp_vec = plrGraph.patches[edgepatch[i]].pointIndexOfPatch;
	//			tp_vec.pop_back();
	//			for (unsigned int j = 0; j < tp_vec.size() - 1; ++j)
	//			{
	//				int node1 = tp_vec[j];
	//				for (unsigned int k = j + 1; k < tp_vec.size(); ++k)
	//				{
	//					int node2 = tp_vec[k];
	//					double tmp_d = LineSegmentLength(sgraph.nodeList[node1].position, sgraph.nodeList[node2].position);
	//					if (tmp_d > fdis)
	//					{
	//						fdis = tmp_d;
	//						p1 = node1;
	//						p2 = node2;
	//					}
	//				}
	//			}
	//			//两条路径plist1 plist2
	//			vector<mypoint2f> plist1, plist2;//plist1 和plist2都是从p1到p2的，有序的
	//			for (unsigned int j = 0; j < plrGraph.patches[edgepatch[i]].pointIndexOfPatch.size() - 1; ++j)//12341+234
	//			{
	//				tp_vec.push_back(plrGraph.patches[edgepatch[i]].pointIndexOfPatch[j]);
	//			}
	//			for (unsigned int j = 0; j < tp_vec.size(); ++j)
	//			{
	//				if (tp_vec[j] == p1)
	//				{
	//					int k = j;
	//					while (tp_vec[k] != p2)
	//					{
	//						if (!plist1.empty())
	//							plist1.pop_back();
	//						getPListOfaCurve(sgraph.nodeList[tp_vec[k]].position, sgraph.nodeList[tp_vec[k + 1]].position, tp_vec[k], tp_vec[k + 1], &plist1);
	//						k++;
	//						if (k >= tp_vec.size())
	//							break;
	//					}
	//					break;
	//				}
	//			}
	//			reverse(tp_vec.begin(), tp_vec.end());
	//			for (unsigned int j = 0; j < tp_vec.size(); ++j)
	//			{
	//				if (tp_vec[j] == p1)
	//				{
	//					int k = j;
	//					while (tp_vec[k] != p2)
	//					{
	//						if (!plist2.empty())
	//							plist2.pop_back();
	//						getPListOfaCurve(sgraph.nodeList[tp_vec[k]].position, sgraph.nodeList[tp_vec[k + 1]].position, tp_vec[k], tp_vec[k + 1], &plist2);
	//						k++;
	//						if (k >= tp_vec.size())
	//							break;
	//					}
	//					break;
	//				}
	//			}
	//			//计算中间路径mid_plist, 存入临时变量 addplist; 如果patch非常小，路径可能就4，5个点，需注意一下
	//			vector<vector<mypoint2f>> addplist;
	//			vector<mypoint2f> mid_plist;
	//			getMidPlist(plist1, plist2, &mid_plist); //p1-->p2
	//			addplist.push_back(mid_plist);
	//			//删除patch，并且删除patch边界处de边sgraph.nodelist.firstedge->delete   sgraph.nodelist.orderarc-->delete
	//			plrGraph.patches[edgepatch[i]].type = -1;
	//			deleOverlayPatches(plrGraph.patches[edgepatch[i]].patch_index, 0);
	//			for (unsigned int j = 0; j < plrGraph.patches[edgepatch[i]].edges_vec.size(); ++j)
	//			{
	//				int edge_s = plrGraph.patches[edgepatch[i]].edges_vec[j].startp_index;
	//				int edge_e = plrGraph.patches[edgepatch[i]].edges_vec[j].endp_index;
	//				ArcNode* anode = sgraph.nodeList[edge_s].firstEdge;
	//				while (anode)
	//				{
	//					if (anode->adjVertex == edge_e)
	//					{
	//						if (anode->patchesIndex.empty() && anode->innerLine == 0)
	//						{
	//							removeArcNode(edge_s, edge_e, 0);//0表示真删除
	//							removeArcNode(edge_e, edge_s, 0);
	//						}
	//						break;
	//					}
	//					anode = anode->next;
	//				}
	//			}
	//			//修改与已删除patch相邻的patch的线.先找到非p1p2的端点，这些端点相连的线又分为3种情况处理：有patch相连，无patch相连（outline or innerline）
	//			vector<int> leftpoints;//leftpoints 现在是无序的
	//			for (unsigned int j = 0; j < plrGraph.patches[edgepatch[i]].pointIndexOfPatch.size() - 1; ++j)
	//			{
	//				if (plrGraph.patches[edgepatch[i]].pointIndexOfPatch[j] != p1 &&plrGraph.patches[edgepatch[i]].pointIndexOfPatch[j] != p2)
	//				{
	//					leftpoints.push_back(plrGraph.patches[edgepatch[i]].pointIndexOfPatch[j]);
	//				}
	//			}
	//			if (!leftpoints.empty())
	//			{
	//				//计算leftpoint到midplist的投影点(即距离最近的点)
	//				vector<mypoint2f> project_p;
	//				for (unsigned int j = 0; j < leftpoints.size(); ++j)
	//				{
	//					double mini_dis = 6553500;
	//					mypoint2f pro_p(0, 0);
	//					for (unsigned int k = 0; k < mid_plist.size(); ++k)
	//					{
	//						double tmp_dis = LineSegmentLength(mid_plist[k], sgraph.nodeList[leftpoints[j]].position);
	//						if (tmp_dis< mini_dis)
	//						{
	//							mini_dis = tmp_dis;
	//							pro_p = mid_plist[k];
	//						}
	//					}
	//					project_p.push_back(pro_p);
	//				}
	//				//根据投影点， 对leftpoint排序: p1--> leftpoint -->p2,存储到leftpoint_order
	//				vector<int> leftpoint_order;
	//				vector<mypoint2f> project_p_order;
	//				if (leftpoints.size() == 1)
	//				{
	//					leftpoint_order.push_back(p1); leftpoint_order.push_back(leftpoints[0]); leftpoint_order.push_back(p2);
	//					project_p_order.push_back(sgraph.nodeList[p1].position);
	//					project_p_order.push_back(project_p[0]);
	//					project_p_order.push_back(sgraph.nodeList[p2].position);
	//				}
	//				else
	//				{
	//					leftpoint_order.push_back(p1);
	//					project_p_order.push_back(sgraph.nodeList[p1].position);
	//					vector<pair<int, double>> ptp_dis;// ptp_dis.push_back(0);
	//					mypoint2f p1_position = sgraph.nodeList[p1].position;
	//					mypoint2f p2_position = sgraph.nodeList[p2].position;
	//					for (unsigned int j = 0; j < project_p.size(); ++j)
	//					{
	//						double tmp = LineSegmentLength(p1_position, project_p[j]);
	//						bool inflag = false;
	//						vector<pair<int, double>>::iterator iter_op = ptp_dis.begin();//int：patch的编号，double:patch的周长
	//						while (iter_op != ptp_dis.end())
	//						{
	//							if (tmp < iter_op->second)//从xiao到da排序
	//							{
	//								ptp_dis.insert(iter_op, make_pair(leftpoints[j], tmp));
	//								inflag = true;
	//								break;
	//							}
	//							iter_op++;
	//						}
	//						if (inflag == false)
	//							ptp_dis.push_back(make_pair(leftpoints[j], tmp));
	//					}
	//					for (unsigned int j = 0; j < ptp_dis.size(); ++j)
	//					{
	//						leftpoint_order.push_back(ptp_dis[j].first);
	//						for (unsigned int k = 0; k < leftpoints.size(); ++k)
	//						{
	//							if (leftpoints[k] == ptp_dis[j].first)
	//							{
	//								project_p_order.push_back(project_p[k]);
	//								break;
	//							}
	//						}
	//					}
	//					leftpoint_order.push_back(p2);
	//					project_p_order.push_back(sgraph.nodeList[p2].position);
	//				}
	//				//一个一个的处理（改变node顶点坐标，替换polyline等）
	//				//sgraph.nodeList[leftpoints[j]].position 和 sgraph.nodeList[anode->adjVertex].position 在原来的patch上，先不做变换，在for之后变换。若不在原patch上例如outline或另一个patch上度为2的点则改变坐标？
	//				vector<int> usednode;//该变量，存储countedge=2的点，表示已经用过，用在specialTransformEdge()函数里
	//				for (unsigned int j = 0; j < leftpoints.size(); ++j)
	//				{
	//					int alp = leftpoints[j];
	//					ArcNode* anode = sgraph.nodeList[alp].firstEdge;
	//					mypoint2f trans_p = project_p[j] - sgraph.nodeList[leftpoints[j]].position;//增量
	//					while (anode)
	//					{
	//						if (anode->innerLine == 2)//condition1:无patch相连的线，即outline ，整体向内平移
	//						{
	//							sgraph.nodeList[anode->adjVertex].position = trans_p + sgraph.nodeList[anode->adjVertex].position;
	//							translateEdge(anode->curveType, anode->curveIndex, trans_p);
	//							//void translateOutline(leftpoints[j],anode->adjVertex,trans_p);
	//						}
	//						else if (anode->innerLine == 1)
	//						{
	//							//还没有考虑这种情况
	//							cout << endl;
	//						}
	//						else if (!anode->patchesIndex.empty())
	//						{
	//							if (isInVector(tp_vec, anode->adjVertex))//是不是在原patch周边上， tp_vec = plrGraph.patches[i].pointIndexOfPatch;在上面定义过
	//							{
	//								//贴着原来patch的edge，截取midplist的一段当做新的edge
	//								//截取+替换
	//								vector<mypoint2f> cutout;//pstart pend是midplist中的坐标位置
	//								mypoint2f pstart = project_p[j];//sgraph.nodeList[leftpoints[j]].position;
	//								mypoint2f pend = sgraph.nodeList[anode->adjVertex].position;
	//								if (isInVector(leftpoints, anode->adjVertex))
	//								{
	//									for (unsigned int k = 0; k < project_p.size(); ++k)
	//									{
	//										if (leftpoints[k] == anode->adjVertex)
	//										{
	//											pend = project_p[k];
	//											break;
	//										}
	//									}
	//								}
	//								truncateAnEdge(mid_plist, &cutout, pstart, pend);
	//								substituteEdge(anode->curveType, anode->curveIndex, cutout);
	//								//shansgraph.nodeList[leftpoints[j]].position = cutout.front();//shan
	//								//sgraph.nodeList[anode->adjVertex].position = cutout.back();//先不做变换，会影响会面的操作
	//							}
	//							else
	//							{
	//								//不在原patch上，且是另一个patch上的一点。
	//								//1)若anode->adjVertex的countEdgeNum(，，0)!=2  发生形变（膨胀？）增量*百分比+原坐标=new edge
	//								//2)countedge()=2，找到结束点，判断结束点是否在原patch上
	//								int ecount = countEdgeNum(anode->adjVertex, 0);
	//								if (ecount != 2)
	//								{
	//									mypoint2f pstart = sgraph.nodeList[leftpoints[j]].position;//增量100%
	//									mypoint2f pend = sgraph.nodeList[anode->adjVertex].position;//增量0%
	//									vector<mypoint2f> ori_plist;
	//									vector<mypoint2f> new_plist;
	//									getPListOfaCurve(pstart, pend, leftpoints[j], anode->adjVertex, &ori_plist);
	//									transformEdge(&new_plist, trans_p, &ori_plist);  //增量在while循环外最上面定义了	
	//									//sgraph.nodeList[leftpoints[j]].position = new_plist.front();//shan
	//									//sgraph.nodeList[anode->adjVertex].position = new_plist.back();//不用改，坐标没有发生变化
	//									substituteEdge(anode->curveType, anode->curveIndex, new_plist);
	//								}
	//								else
	//								{
	//									//已知leftpoints[j]和anode->ver, 再往下寻找，一直找到度！=2的结束点。然后判断结束点是在原patch上海市其他patch上
	//									//sgraph.nodeList[anode->adjVertex].position = trans_p + sgraph.nodeList[anode->adjVertex].position;
	//									//translateEdge(anode->curveType, anode->curveIndex, trans_p);							
	//									specialTransformEdge(leftpoints[j], anode, i, &leftpoints, &project_p, &usednode);//i式原来被删掉的patch的index									
	//								}
	//							}
	//						}
	//						anode = anode->next;
	//					}
	//					//sgraph.nodeList[leftpoints[j]].position = project_p[j];
	//				}
	//				for (unsigned int j = 0; j < leftpoints.size(); ++j)
	//				{
	//					sgraph.nodeList[leftpoints[j]].position = project_p[j];
	//				}
	//				//根据leftpoint_order & project_p_order对相关点进行连接，有的连接之后变成了patch的一部分有的变成了outline，需要添加到plrgraph.outline中
	//				vector<pair<int, int>> add_arcnode;
	//				for (unsigned int j = 1; j < leftpoint_order.size(); ++j)
	//				{
	//					//先连接0-1，1-2，2-3，，，
	//					int tp1 = leftpoint_order[j - 1];
	//					int tp2 = leftpoint_order[j];
	//					ArcNode* anode = getArcNode(tp1, tp2);
	//					if (anode == NULL)
	//					{
	//						ArcNode* newnode = new ArcNode;
	//						newnode->adjVertex = tp2;
	//						newnode->curveType = 'l';
	//						vector<mypoint2f> cutout;
	//						mypoint2f pstart = project_p_order[j - 1];  //tp1的投影点
	//						mypoint2f pend = project_p_order[j];  //tp2的投影点
	//						truncateAnEdge(mid_plist, &cutout, pstart, pend);
	//						PolyLine apl; apl.pointlist = cutout;
	//						polyline_vec.push_back(apl);
	//						newnode->curveIndex = polyline_vec.size() - 1;
	//						newnode->innerLine = 0;
	//						newnode->next = sgraph.nodeList[tp1].firstEdge;
	//						sgraph.nodeList[tp1].firstEdge = newnode;
	//						ArcNode* newnode2 = new ArcNode;
	//						newnode2->adjVertex = tp1;
	//						newnode2->curveType = 'l';
	//						newnode2->curveIndex = polyline_vec.size() - 1;
	//						newnode2->innerLine = 0;
	//						newnode2->next = sgraph.nodeList[tp2].firstEdge;
	//						sgraph.nodeList[tp2].firstEdge = newnode2;
	//						sgraph.nodeList[tp2].orderedArcs;
	//					}
	//					//leftpoint_order[j]->leftpoint_order[j-2/j-3.../0]向后排查有没有 应该断开的edge，
	//					//从这条长arc得到对应的patchindex，去修改patch->plist和patch->edge_vec。同时修改 新加入的arcnode的patchindex
	//					int k = j - 2;
	//					while (k >= 0)
	//					{
	//						ArcNode* enode = getArcNode(leftpoint_order[j], leftpoint_order[k]);
	//						if (enode != NULL)
	//						{
	//							int pchindex = enode->patchesIndex[0];//或者for循环，找到patchindex不是i的那个
	//							vector<int> pch_verlist = plrGraph.patches[pchindex].pointIndexOfPatch;															
	//							vector<int> origi_piece;
	//							origi_piece.push_back(leftpoint_order[j]);
	//							origi_piece.push_back(leftpoint_order[k]);
	//							vector<int> target_piece;
	//							int g = j;
	//							while (g != k)
	//							{
	//								target_piece.push_back(leftpoint_order[g]);
	//								g--;
	//							}
	//							target_piece.push_back(leftpoint_order[k]);
	//							bool tar_flag = false;
	//							for (unsigned int ti = 1; ti < target_piece.size() - 1; ++ti)
	//							{
	//								if (isInVector(pch_verlist, target_piece[ti]))//需要添加的点在patch中本来就有，则不做修改
	//								{
	//									tar_flag = true;
	//									break;
	//								}
	//							}
	//							if (tar_flag == false)
	//							{
	//								Patch *modify_pch = getAPatch(pchindex); //modify_pch->pointIndexOfPatch; modify_pch->edges_vec;
	//								modifyPatch_plist(modify_pch, origi_piece, target_piece);
	//								removeArcNode(leftpoint_order[j], leftpoint_order[k], 0);//删除已经拆分掉的arc
	//								removeArcNode(leftpoint_order[k], leftpoint_order[j], 0);
	//							}							
	//						}
	//						k--;
	//					}
	//					//新添加的arc（tp1,tp2）,还不确定是outline还是patch的一部分，先记录下来,for循环外进行判断
	//					if (anode == NULL)
	//					{
	//						add_arcnode.push_back(make_pair(tp1, tp2));
	//					}
	//				}
	//				if (!add_arcnode.empty())
	//				{
	//					for (unsigned int ai = 0; ai < add_arcnode.size(); ++ai)
	//					{
	//						int n1 = add_arcnode[ai].first;
	//						int n2 = add_arcnode[ai].second;
	//						ArcNode* nnarc = getArcNode(n1, n2); //new arcnode 的innerlin = 0且patchindex.empty()说明是一条新的线
	//						if (nnarc->innerLine == 0 && nnarc->patchesIndex.empty())
	//						{
	//							nnarc->innerLine = 2;
	//							nnarc = getArcNode(n2, n1);//反向设置一遍
	//							nnarc->innerLine = 2;
	//							vector<int> tmp_vec;
	//							tmp_vec.push_back(n1);
	//							tmp_vec.push_back(n2);
	//							plrGraph.outlines.push_back(tmp_vec);
	//						}
	//					}
	//				}
	//			}
	//		}
	//	}
	//}

/////////////////////////////////////////////////////////

	for (int i = 0; i <pchsize; ++i)
	{
		bool eflag = false;
		for (unsigned int k = 0; k < plrGraph.patches[i].edges_vec.size(); ++k)
		{
			if (plrGraph.patches[i].edges_vec[k].adjPatchIndex.size() == 1)
			{
				eflag = true; break;
			}
		}
		if (plrGraph.patches[i].type >= 0 && eflag == true)//可用 且是边界patch
		{
			//getpatch pointlist,判断拐点
			vector<mypoint2f> apatch;
			for (unsigned int j = 1; j < plrGraph.patches[i].pointIndexOfPatch.size(); ++j)
			{
				int verindex1 = plrGraph.patches[i].pointIndexOfPatch.at(j - 1);
				int verindex2 = plrGraph.patches[i].pointIndexOfPatch.at(j);
				mypoint2f point1 = sgraph.nodeList[verindex1].position;//=lineWithCommonPoints[verindex1].first
				mypoint2f point2 = sgraph.nodeList[verindex2].position;
				if (!apatch.empty())
					apatch.pop_back();
				getPListOfaCurve(point1, point2, verindex1, verindex2, &apatch,0);//最终得到的apatch 首=尾
			}
			if (isPointRoughEqual(apatch.front(), apatch.back()))
				apatch.pop_back();
			mypoint2f pfront = apatch.front();
			mypoint2f pback = apatch.back();
			apatch.insert(apatch.begin(), pback);
			apatch.push_back(pfront);
			int turn_count = 0;
			for (unsigned int j = 1; j < apatch.size() - 1; ++j)
			{
				mypoint2f vec1 = apatch[j] - apatch[j - 1];
				mypoint2f vec2 = apatch[j + 1] - apatch[j];
				double cos_theta = getCosineAngle(apatch[j - 1], apatch[j], apatch[j], apatch[j + 1]);
				if (cos_theta >= -1 &&cos_theta <= 0.707106)//转角在45到180认为是一个拐点   135-->-0.7071  150--> -0.866
				{
					turn_count++; //cout << cos_theta << endl;
				}
			}
			if (turn_count < 4)//拐点小于4个才可以 腐蚀成一条线？？
			{
				//找最远的两个点,
				int p1, p2; double fdis=0;
				vector<int> tp_vec = plrGraph.patches[i].pointIndexOfPatch;
				tp_vec.pop_back();
				for (unsigned int j = 0; j < tp_vec.size()-1; ++j)
				{
					int node1 = tp_vec[j];
					for (unsigned int k =j+1; k < tp_vec.size(); ++k)
					{
						int node2 = tp_vec[k];
						double tmp_d = LineSegmentLength(sgraph.nodeList[node1].position, sgraph.nodeList[node2].position);
						if (tmp_d > fdis)
						{
							fdis = tmp_d;
							p1 = node1;
							p2 = node2;
						}
					}
				}
				//两条路径plist1 plist2
				vector<mypoint2f> plist1, plist2;//plist1 和plist2都是从p1到p2的，有序的
				for (unsigned int j = 0; j < plrGraph.patches[i].pointIndexOfPatch.size() - 1; ++j)//12341+234
				{
					tp_vec.push_back(plrGraph.patches[i].pointIndexOfPatch[j]);
				}
				for (unsigned int j = 0; j < tp_vec.size(); ++j)
				{
					if (tp_vec[j] == p1)
					{
						int k = j;
						while (tp_vec[k] != p2)
						{				
							if (!plist1.empty())
								plist1.pop_back();
							getPListOfaCurve(sgraph.nodeList[tp_vec[k]].position, sgraph.nodeList[tp_vec[k + 1]].position, tp_vec[k], tp_vec[k + 1], &plist1,0);
							k++;
							if (k >= tp_vec.size())
								break;
						}
						break;
					}
				}
				reverse(tp_vec.begin(), tp_vec.end());
				for (unsigned int j = 0; j < tp_vec.size(); ++j)
				{
					if (tp_vec[j] == p1)
					{
						int k = j;
						while (tp_vec[k] != p2)
						{
							if (!plist2.empty())
								plist2.pop_back();
							getPListOfaCurve(sgraph.nodeList[tp_vec[k]].position, sgraph.nodeList[tp_vec[k + 1]].position, tp_vec[k], tp_vec[k + 1], &plist2,0);
							k++;
							if (k >= tp_vec.size())
								break;
						}
						break;
					}
				}
				//计算中间路径mid_plist, 存入临时变量 addplist; 如果patch非常小，路径可能就4，5个点，需注意一下
				vector<vector<mypoint2f>> addplist;
				vector<mypoint2f> mid_plist;
				getMidPlist(plist1, plist2, &mid_plist); //p1-->p2
				addplist.push_back(mid_plist);
				//删除patch，并且删除patch边界处de边sgraph.nodelist.firstedge->delete   sgraph.nodelist.orderarc-->delete
				plrGraph.patches[i].type = -1;
				deleOverlayPatches(plrGraph.patches[i].patch_index, 0);
				for (unsigned int j = 0; j < plrGraph.patches[i].edges_vec.size(); ++j)
				{
					int edge_s = plrGraph.patches[i].edges_vec[j].startp_index;
					int edge_e = plrGraph.patches[i].edges_vec[j].endp_index;
					ArcNode* anode = sgraph.nodeList[edge_s].firstEdge;
					while (anode)
					{
						if (anode->adjVertex == edge_e)
						{
							if (anode->patchesIndex.empty() && anode->lineFlag == 0)
							{				
								removeArcNode(edge_s, edge_e, 0);//0表示真删除
								removeArcNode(edge_e, edge_s, 0);
							}
							break;
						}
						anode = anode->next;
					}
				}		
				//修改与已删除patch相邻的patch的线.先找到非p1p2的端点，这些端点相连的线又分为3种情况处理：有patch相连，无patch相连（outline or innerline）
				vector<int> leftpoints;//leftpoints 现在是无序的
				for (unsigned int j = 0; j < plrGraph.patches[i].pointIndexOfPatch.size()-1; ++j)
				{
					if (plrGraph.patches[i].pointIndexOfPatch[j] != p1 &&plrGraph.patches[i].pointIndexOfPatch[j] != p2)
					{
						leftpoints.push_back(plrGraph.patches[i].pointIndexOfPatch[j]);
					}
				}
				if (!leftpoints.empty())
				{
					//计算leftpoint到midplist的投影点(即距离最近的点)
					vector<mypoint2f> project_p;
					for (unsigned int j = 0; j < leftpoints.size(); ++j)
					{
						double mini_dis = 6553500;
						mypoint2f pro_p(0, 0);
						for (unsigned int k = 0; k < mid_plist.size(); ++k)
						{
							double tmp_dis = LineSegmentLength(mid_plist[k], sgraph.nodeList[leftpoints[j]].position);
							if (tmp_dis< mini_dis)
							{
								mini_dis = tmp_dis;
								pro_p = mid_plist[k];
							}
						}
						project_p.push_back(pro_p);
					}
					//根据投影点， 对leftpoint排序: p1--> leftpoint -->p2,存储到leftpoint_order
					vector<int> leftpoint_order;
					vector<mypoint2f> project_p_order;
					if (leftpoints.size() == 1)
					{
						leftpoint_order.push_back(p1); leftpoint_order.push_back(leftpoints[0]); leftpoint_order.push_back(p2);
						project_p_order.push_back(sgraph.nodeList[p1].position);
						project_p_order.push_back(project_p[0]);
						project_p_order.push_back(sgraph.nodeList[p2].position);
					}
					else
					{
						leftpoint_order.push_back(p1);
						project_p_order.push_back(sgraph.nodeList[p1].position);
						vector<pair<int,double>> ptp_dis;// ptp_dis.push_back(0);
						mypoint2f p1_position = sgraph.nodeList[p1].position;
						mypoint2f p2_position = sgraph.nodeList[p2].position;		
						for (unsigned int j = 0; j < project_p.size(); ++j)
						{
							double tmp = LineSegmentLength(p1_position, project_p[j]);
							bool inflag = false;
							vector<pair<int, double>>::iterator iter_op = ptp_dis.begin();//int：patch的编号，double:patch的周长
							while (iter_op != ptp_dis.end())
							{
								if (tmp < iter_op->second)//从xiao到da排序
								{
									ptp_dis.insert(iter_op, make_pair(leftpoints[j],tmp));
									inflag = true;
									break;
								}
								iter_op++;
							}
							if (inflag == false)
								ptp_dis.push_back(make_pair(leftpoints[j],tmp));
						}
						for (unsigned int j = 0; j < ptp_dis.size(); ++j)
						{
							leftpoint_order.push_back(ptp_dis[j].first);
							for (unsigned int k = 0; k < leftpoints.size(); ++k)
							{
								if (leftpoints[k] == ptp_dis[j].first)
								{
									project_p_order.push_back(project_p[k]);
									break;
								}
							}
						}
						leftpoint_order.push_back(p2);
						project_p_order.push_back(sgraph.nodeList[p2].position);
					}			
					//一个一个的处理（改变node顶点坐标，替换polyline等）
					//sgraph.nodeList[leftpoints[j]].position 和 sgraph.nodeList[anode->adjVertex].position 在原来的patch上，先不做变换，在for之后变换。若不在原patch上例如outline或另一个patch上度为2的点则改变坐标？
					vector<int> usednode;//该变量，存储countedge=2的点，表示已经用过，用在specialTransformEdge()函数里
					for (unsigned int j = 0; j < leftpoints.size(); ++j)
					{
						int alp = leftpoints[j];
						ArcNode* anode = sgraph.nodeList[alp].firstEdge;
						mypoint2f trans_p = project_p[j] - sgraph.nodeList[leftpoints[j]].position;//增量
						while (anode)
						{					
							if (anode->lineFlag == 2)//condition1:无patch相连的线，即outline ，整体向内平移
							{
								//sgraph.nodeList[anode->adjVertex].position = trans_p + sgraph.nodeList[anode->adjVertex].position;
								//translateEdge(anode->curveType, anode->curveIndex, trans_p);
								translateOutline(leftpoints[j],anode,trans_p);
							}
							else if (anode->lineFlag == 1)
							{
								//还没有考虑这种情况
								cout << endl;
							}
							else if (!anode->patchesIndex.empty())
							{
								if (isInVector(tp_vec,anode->adjVertex))//是不是在原patch周边上， tp_vec = plrGraph.patches[i].pointIndexOfPatch;在上面定义过
								{
									//贴着原来patch的edge，截取midplist的一段当做新的edge
									//截取+替换
									vector<mypoint2f> cutout;//pstart pend是midplist中的坐标位置
									mypoint2f pstart = project_p[j];//sgraph.nodeList[leftpoints[j]].position;
									mypoint2f pend = sgraph.nodeList[anode->adjVertex].position;
									if (isInVector(leftpoints, anode->adjVertex))
									{
										for (unsigned int k = 0; k < project_p.size(); ++k)
										{
											if (leftpoints[k] == anode->adjVertex)
											{
												pend = project_p[k];
												break;
											}
										}
									}
									truncateAnEdge(mid_plist, &cutout, pstart, pend);
									substituteEdge(anode->curveType, anode->curveIndex, cutout);
									//shansgraph.nodeList[leftpoints[j]].position = cutout.front();//shan
									//sgraph.nodeList[anode->adjVertex].position = cutout.back();//先不做变换，会影响会面的操作
								}
								else
								{
									//不在原patch上，且是另一个patch上的一点。
									//1)若anode->adjVertex的countEdgeNum(，，0)!=2  发生形变（膨胀？）增量*百分比+原坐标=new edge
									//2)countedge()=2，找到结束点，判断结束点是否在原patch上
									int ecount = countEdgeNum(anode->adjVertex, 0);
									if (ecount != 2)
									{		
										mypoint2f pstart = sgraph.nodeList[leftpoints[j]].position;//增量100%
										mypoint2f pend = sgraph.nodeList[anode->adjVertex].position;//增量0%
										vector<mypoint2f> ori_plist;
										vector<mypoint2f> new_plist;
										getPListOfaCurve(pstart, pend, leftpoints[j], anode->adjVertex, &ori_plist,0);
										transformEdge(&new_plist, trans_p, &ori_plist);  //增量在while循环外最上面定义了	
										//sgraph.nodeList[leftpoints[j]].position = new_plist.front();//shan
										//sgraph.nodeList[anode->adjVertex].position = new_plist.back();//不用改，坐标没有发生变化
										substituteEdge(anode->curveType, anode->curveIndex, new_plist);
									}
									else
									{
										//已知leftpoints[j]和anode->ver, 再往下寻找，一直找到度！=2的结束点。然后判断结束点是在原patch上海市其他patch上
										//sgraph.nodeList[anode->adjVertex].position = trans_p + sgraph.nodeList[anode->adjVertex].position;
										//translateEdge(anode->curveType, anode->curveIndex, trans_p);							
										specialTransformEdge(leftpoints[j], anode, i, &leftpoints, &project_p, &usednode);//i式原来被删掉的patch的index									
									}						
								}								
							}
							anode = anode->next;
						}
						//sgraph.nodeList[leftpoints[j]].position = project_p[j];
					}
					for (unsigned int j = 0; j < leftpoints.size(); ++j)
					{
						sgraph.nodeList[leftpoints[j]].position = project_p[j];
					}
					//根据leftpoint_order & project_p_order对相关点进行连接，有的连接之后变成了patch的一部分有的变成了outline，需要添加到plrgraph.outline中
					vector<pair<int, int>> add_arcnode;
					for (unsigned int j = 1; j < leftpoint_order.size(); ++j)
					{
						//先连接0-1，1-2，2-3，，，
						int tp1 = leftpoint_order[j - 1];
						int tp2 = leftpoint_order[j];
						ArcNode* anode = getArcNode(tp1, tp2);
						if (anode == NULL)
						{
							ArcNode* newnode = new ArcNode;
							newnode->adjVertex = tp2;
							newnode->curveType = 'l';
							vector<mypoint2f> cutout;
							mypoint2f pstart = project_p_order[j - 1];  //tp1的投影点
							mypoint2f pend = project_p_order[j];  //tp2的投影点
							truncateAnEdge(mid_plist, &cutout, pstart, pend);
							PolyLine apl; apl.pointlist = cutout;
							polyline_vec.push_back(apl);
							newnode->curveIndex = polyline_vec.size() - 1;
							newnode->lineFlag = 0;
							newnode->next = sgraph.nodeList[tp1].firstEdge;
							sgraph.nodeList[tp1].firstEdge = newnode;

							ArcNode* newnode2 = new ArcNode;
							newnode2->adjVertex = tp1;
							newnode2->curveType = 'l';
							newnode2->curveIndex = polyline_vec.size() - 1;
							newnode2->lineFlag = 0;
							newnode2->next = sgraph.nodeList[tp2].firstEdge;
							sgraph.nodeList[tp2].firstEdge = newnode2;
							sgraph.nodeList[tp2].orderedArcs;
						}
						//leftpoint_order[j]->leftpoint_order[j-2/j-3.../0]向后排查有没有 应该断开的edge，
						//从这条长arc得到对应的patchindex，去修改patch->plist和patch->edge_vec。同时修改 新加入的arcnode的patchindex
						int k = j - 2;
						while (k >= 0)
						{
							ArcNode* enode = getArcNode(leftpoint_order[j], leftpoint_order[k]);
							if (enode != NULL)
							{
								int pchindex = enode->patchesIndex[0];//或者for循环，找到patchindex不是i的那个
								vector<int> pch_verlist = plrGraph.patches[pchindex].pointIndexOfPatch;
								vector<int> origi_piece;
								origi_piece.push_back(leftpoint_order[j]);
								origi_piece.push_back(leftpoint_order[k]);
								vector<int> target_piece;
								int g = j;
								while (g != k)
								{
									target_piece.push_back(leftpoint_order[g]);
									g--;
								}
								target_piece.push_back(leftpoint_order[k]);
								bool tar_flag = false;
								for (unsigned int ti = 1; ti < target_piece.size() - 1; ++ti)
								{
									if (isInVector(pch_verlist, target_piece[ti]))//需要添加的点在patch中本来就有，则不做修改
									{
										tar_flag = true;
										break;
									}
								}
								if (tar_flag == false)
								{
									Patch *modify_pch = getAPatch(pchindex); //modify_pch->pointIndexOfPatch; modify_pch->edges_vec;
									modifyPatch_plist(modify_pch, origi_piece, target_piece);
									removeArcNode(leftpoint_order[j], leftpoint_order[k], 0);//删除已经拆分掉的arc
									removeArcNode(leftpoint_order[k], leftpoint_order[j], 0);
								}
							}
							k--;
						}
						//新添加的arc（tp1,tp2）,还不确定是outline还是patch的一部分，先记录下来,for循环外进行判断
						if (anode == NULL)
						{
							add_arcnode.push_back(make_pair(tp1, tp2));
						}
					}
					if (!add_arcnode.empty())
					{
						for (unsigned int ai = 0; ai < add_arcnode.size(); ++ai)
						{
							int n1 = add_arcnode[ai].first;
							int n2 = add_arcnode[ai].second;
							ArcNode* nnarc = getArcNode(n1, n2); //new arcnode 的innerlin = 0且patchindex.empty()说明是一条新的线
							if (nnarc->lineFlag == 0 && nnarc->patchesIndex.empty())
							{
								nnarc->lineFlag = 2;
								nnarc = getArcNode(n2, n1);//反向设置一遍
								nnarc->lineFlag = 2;
								vector<int> tmp_vec;
								tmp_vec.push_back(n1);
								tmp_vec.push_back(n2);
								ICurve icv;
								icv.type = 0;
								icv.pt_vec = tmp_vec;
								plrGraph.attachlines.push_back(icv);
							}
						}
					}
		
				}
			}
		}
	}
}
void MyGraph::getMidPlist(vector<mypoint2f> plist1, vector<mypoint2f> plist2, vector<mypoint2f> *midplist)
{
	midplist->push_back(plist1.front());
	if (plist1.size() > plist2.size())
	{
		if (plist2.size() > 2)
		{
			int j_start = 1;
			for (unsigned int i = 1; i < plist1.size() - 1; ++i)
			{
				mypoint2f tmp_p1 = plist1[i];
				double minidis = 6553500;
				int index_plist2 = 0;
				for (unsigned int j = j_start; j < plist2.size(); ++j)
				{
					double dis = LineSegmentLength(tmp_p1, plist2[j]);
					if (dis < minidis)
					{
						minidis = dis;
						index_plist2 = j;
					}
					else
						break;
				}
				j_start = index_plist2;
				mypoint2f midp((tmp_p1.x + plist2[index_plist2].x) / 2, (tmp_p1.y + plist2[index_plist2].y) / 2);
				midplist->push_back(midp);
			}
		}
		else
		{
			//plist2.size() = 2
			mypoint2f midp((plist2[0].x + plist2[1].x) / 2, (plist2[0].y + plist2[1].y) / 2);
			midplist->push_back(midp);
		}
	}
	else
	{
		if (plist1.size() > 2)
		{
			/*int j_start = 1;
			for (unsigned int i = 1; i < plist2.size() - 1; ++i)
			{
				mypoint2f tmp_p1 = plist2[i];
				double minidis = 6553500;
				int index_plist2 = 0;
				for (unsigned int j = j_start; j < plist1.size(); ++j)
				{
					double dis = LineSegmentLength(tmp_p1, plist1[j]);
					if (dis < minidis)
					{
						minidis = dis;
						index_plist2 = j;
					}
					else
						break;
				}
				j_start = index_plist2;
				mypoint2f midp((tmp_p1.x + plist1[index_plist2].x) / 2, (tmp_p1.y + plist1[index_plist2].y) / 2);
				midplist->push_back(midp);
			}*/

			for (unsigned int i = 1; i < plist2.size() - 1; ++i)
			{
				mypoint2f tmp_p1 = plist2[i];
				double minidis = 6553500;
				int index_plist2 = 0;
				for (unsigned int j = 0; j < plist1.size(); ++j)
				{
					double dis = LineSegmentLength(tmp_p1, plist1[j]);
					if (dis < minidis)
					{
						minidis = dis;
						index_plist2 = j;
					}
					else
						break;
				}
				mypoint2f midp((tmp_p1.x + plist1[index_plist2].x) / 2, (tmp_p1.y + plist1[index_plist2].y) / 2);
				midplist->push_back(midp);
			}
		}
		else
		{
			//plist2.size() = 2
			mypoint2f midp((plist1[0].x + plist1[1].x) / 2, (plist1[0].y + plist1[1].y) / 2);
			midplist->push_back(midp);
		}
	}
	midplist->push_back(plist1.back());
}
bool MyGraph::isInVector(vector<int> nodevec, int node)
{
	bool ff = false;
	if (!nodevec.empty())
	{
		for (unsigned int i = 0; i < nodevec.size(); ++i)
		{
			if (nodevec[i] == node)
			{
				ff = true;
				break;
			}
		}
	}
	return ff;
}
void MyGraph::translateOutline(int startp, ArcNode* anodep, mypoint2f trans_dis)
{
	sgraph.nodeList[anodep->adjVertex].position = trans_dis + sgraph.nodeList[anodep->adjVertex].position;
	translateEdge(anodep->curveType, anodep->curveIndex, trans_dis);
	int nsize = sgraph.nodeList.size();
	int* visit = new int[nsize];
	for (int i = 0; i < nsize; ++i)
	{
		visit[i] = 0;
	}
	visit[anodep->adjVertex] = 2;//!!!
	stack<int> ss;
	ArcNode *rv_node = sgraph.nodeList[anodep->adjVertex].firstEdge;
	while (rv_node)
	{
		if (rv_node->adjVertex != startp &&rv_node->lineFlag == 2)
		{
			ss.push(rv_node->adjVertex);
			visit[rv_node->adjVertex] = 1;
			sgraph.nodeList[rv_node->adjVertex].position = trans_dis + sgraph.nodeList[rv_node->adjVertex].position;
			translateEdge(rv_node->curveType, rv_node->curveIndex, trans_dis);
		}
		rv_node = rv_node->next;
	}
	if (!ss.empty())
	{
		while (!ss.empty())
		{
			int tmpver = ss.top();
			ss.pop();
			visit[tmpver] = 2;
			ArcNode* tmpnode = sgraph.nodeList[tmpver].firstEdge;
			while (tmpnode)
			{
				if (visit[tmpnode->adjVertex] != 2)
				{
					if (visit[tmpnode->adjVertex] == 0 && tmpnode->lineFlag == 2)
					{
						ss.push(tmpnode->adjVertex);
						visit[tmpnode->adjVertex] = 1;
						sgraph.nodeList[tmpnode->adjVertex].position = trans_dis + sgraph.nodeList[tmpnode->adjVertex].position;
						translateEdge(tmpnode->curveType, tmpnode->curveIndex, trans_dis);
					}
				}
				tmpnode = tmpnode->next;
			}
		}
	}	
}
void MyGraph::translateEdge(char ctype, int cindex, mypoint2f trans_dis)
{
	if (ctype == 'l')
	{
		for (unsigned int i = 0; i < polyline_vec[cindex].pointlist.size(); ++i)
		{
			polyline_vec[cindex].pointlist[i] = polyline_vec[cindex].pointlist[i] +trans_dis;
		}
	}
	else if (ctype == 'c')
	{

	}
	else if (ctype == 'q')
	{

	}
}
void MyGraph::truncateAnEdge(vector<mypoint2f> origi_plist, vector<mypoint2f> *cutout, mypoint2f fromp, mypoint2f endp)
{
	for (unsigned int k = 0; k < origi_plist.size(); ++k)
	{
		if (origi_plist[k] == fromp )
		{
			cutout->push_back(origi_plist[k]);
			for (unsigned int g = k + 1; g < origi_plist.size(); ++g)
			{
				cutout->push_back(origi_plist[g]);
				if (origi_plist[g] == endp)
					break;
			}
			break;
		}
		else if (origi_plist[k] == endp)
		{
			cutout->push_back(origi_plist[k]);
			for (unsigned int g = k + 1; g < origi_plist.size(); ++g)
			{
				cutout->push_back(origi_plist[g]);
				if (origi_plist[g] == fromp)
					break;
			}
			reverse(cutout->begin(), cutout->end());
			break;
		}
	}
}
void MyGraph::substituteEdge(char ctype, int cindex, vector<mypoint2f> newplist)
{
	if (ctype == 'l')
	{	
		polyline_vec[cindex].pointlist.clear(); polyline_vec[cindex].pointlist.swap(vector<mypoint2f>());
		polyline_vec[cindex].pointlist = newplist;
	}
	else if (ctype == 'c')
	{

	}
	else if (ctype == 'q')
	{

	}
}
void MyGraph::transformEdge(vector<mypoint2f> *new_plist, mypoint2f addp, vector<mypoint2f> *origi_plist)
{
	//sp的增量是1，ep的增量是0，中间递减。即保持ep的坐标不变，
		int psize = origi_plist->size();
		double step = 1.0 / (psize-1);
		double factor = 1;
		for (int j = 0; j < psize;++j)
		{
			factor = 1.0 - pow(j*step,2);
			mypoint2f tmp1(addp.x*factor, addp.y*factor);
			mypoint2f tmp2 = tmp1 + origi_plist->at(j);
			new_plist->push_back(tmp2);
		}
		if (factor != 0)
		{
			new_plist->push_back(origi_plist->back());
		}
}
void MyGraph::modifyPatch_plist(Patch *pch, vector<int> origip, vector<int> targetp)
{
	//eg  origip=1,2 targetp=1,2,3;
	//先判断origip的顺序是不是和patch的一样,origip一定是两个
	vector<int> plist = pch->pointIndexOfPatch;
	bool order = false;//false 顺序不一致 true 顺序一致
	for (unsigned int i = 0; i < plist.size(); ++i)
	{
		if (plist[i] == origip[0])
		{
			if (i + 1 < plist.size())
			{
				if (plist[i + 1] == origip[1])
				{
					order = true;
				}
			}
			break;
		}
	}
	if (order == false)
	{
		reverse(targetp.begin(), targetp.end());
		reverse(origip.begin(), origip.end());
	}
	//step1 修改pointIndexOfPatch
	vector<int>::iterator iter = pch->pointIndexOfPatch.begin();
	while (iter != pch->pointIndexOfPatch.end())
	{
		if (*iter == origip[0])
		{
			iter++;
			for (unsigned int j = 1; j < targetp.size() - 1; ++j)
			{
				iter = pch->pointIndexOfPatch.insert(iter, targetp[j]);
				iter++;
			}
			break;
		}
		else
			iter++;
	}
	//step2 修改edge_ved
	vector<EdgeOfPatch>::iterator iter_edge = pch->edges_vec.begin();
	while (iter_edge != pch->edges_vec.end())
	{
		if (iter_edge->startp_index == origip[0] && iter_edge->endp_index == origip[1])
		{
			iter_edge = pch->edges_vec.erase(iter_edge);//先删除旧的edge
			for (unsigned int j = 1; j < targetp.size(); ++j)
			{
				EdgeOfPatch ee;
				ee.startp_index = targetp[j - 1];
				ee.endp_index = targetp[j];
				ee.edgeflag = 0;
				ee.adjPatchIndex.push_back(pch->patch_index);//存入本patch的index
				ArcNode* tmp_arc = getArcNode(targetp[j - 1], targetp[j]);			
				if (!tmp_arc->patchesIndex.empty())
				{
					bool havepch = false;
					int ccc = 0;
					for (unsigned int k = 0; k < tmp_arc->patchesIndex.size(); ++k)
					{
						if (tmp_arc->patchesIndex[k] != pch->patch_index)
						{
							ccc++;
							ee.adjPatchIndex.push_back(tmp_arc->patchesIndex[k]);//存入相邻patch的index,tmp_arc->patchesIndex[k]和pch变成相邻的patch了
							int another_pch = tmp_arc->patchesIndex[k];
							for (unsigned int ei = 0; ei < plrGraph.patches[another_pch].edges_vec.size(); ++ei)
							{
								if (plrGraph.patches[another_pch].edges_vec[ei].startp_index == ee.startp_index &&plrGraph.patches[another_pch].edges_vec[ei].endp_index == ee.endp_index)
								{
									plrGraph.patches[another_pch].edges_vec[ei].adjPatchIndex.push_back(pch->patch_index); break;
								}
								else if (plrGraph.patches[another_pch].edges_vec[ei].startp_index == ee.endp_index &&plrGraph.patches[another_pch].edges_vec[ei].endp_index == ee.startp_index)
								{
									plrGraph.patches[another_pch].edges_vec[ei].adjPatchIndex.push_back(pch->patch_index); break;
								}
							}
						}
						else
							havepch = true;
					}
					if (ccc > 1)
						cout << endl;
					if (havepch == false)
					{
						tmp_arc->patchesIndex.push_back(pch->patch_index);
						tmp_arc = getArcNode(targetp[j], targetp[j-1]);
						tmp_arc->patchesIndex.push_back(pch->patch_index);
					}
				}
				else//给新添的的arcnode 的patchindex赋值，双向的
				{
					tmp_arc->patchesIndex.push_back(pch->patch_index);
					tmp_arc = getArcNode(targetp[j], targetp[j - 1]);
					tmp_arc->patchesIndex.push_back(pch->patch_index);
				}
				iter_edge = pch->edges_vec.insert(iter_edge, ee);//最后把新edge加入
				iter_edge++;//!!!!!!!!!!!!
			}
			break;
		}
		else
			iter_edge++;
	}
}
void MyGraph::specialTransformEdge(int startp, ArcNode *anode, int origi_pch, vector<int> *left_p, vector<mypoint2f> *leftp_project , vector<int> *used_node)
{
	vector<int> apath;
	apath.push_back(startp);
	apath.push_back(anode->adjVertex);
	int lastindex = startp;
	int ecount = countEdgeNum(anode->adjVertex, 0);
	while (ecount == 2)
	{
		ArcNode *tmpnode = sgraph.nodeList[apath.back()].firstEdge;
		while (tmpnode)
		{
			if (tmpnode->adjVertex != lastindex)
			{
				lastindex = apath.back();
				apath.push_back(tmpnode->adjVertex);
				ecount = countEdgeNum(tmpnode->adjVertex, 0);
				break;
			}
			tmpnode = tmpnode->next;
		}
	}
	for (unsigned int i = 1; i < apath.size() - 1; ++i)
	{
		if (isInVector(*used_node, apath[i]) == true)
			return;
		else
		{
			used_node->push_back(apath[i]);
		}
	}
	//先得到这几个edge的origi_plist , 并记录分段点的下标segment_index
	vector<mypoint2f> origi_plist;
	vector<mypoint2f> newplist;

	vector<int> segment_index;
	segment_index.push_back(0);
	for (unsigned int j = 1; j < apath.size(); ++j)
	{
		int verindex1 = apath[j - 1];
		int verindex2 = apath[j];
		mypoint2f point1 = sgraph.nodeList[verindex1].position;//=lineWithCommonPoints[verindex1].first
		mypoint2f point2 = sgraph.nodeList[verindex2].position;
		if (!origi_plist.empty())
			origi_plist.pop_back();
		getPListOfaCurve(point1, point2, verindex1, verindex2, &origi_plist,0);
		segment_index.push_back(origi_plist.size() - 1);
	}
	bool inflag = isInVector(plrGraph.patches[origi_pch].pointIndexOfPatch, apath.back());
	if (inflag == false)
	{
		//若不在原patch里面  ,
		mypoint2f addp(0, 0);//计算增量
		for (unsigned int i = 0; i < left_p->size(); ++i)
		{
			if (left_p->at(i) == startp)
			{
				mypoint2f origip = sgraph.nodeList[startp].position;
				addp = leftp_project->at(i) - origip;
				break;
			}
		}
		transformEdge(&newplist, addp, &origi_plist);  // new edge=增量*百分比+原坐标
		//分段更新polyline  segment_index.size()==apathy.size() ?
		for (unsigned int i = 1; i < segment_index.size(); ++i)
		{
			//得到arcnode
			ArcNode* tmp_node = getArcNode(apath[i - 1], apath[i]);
			//截取新的polyline
			vector<mypoint2f> tmp_pl;
			for (int j = segment_index[i - 1]; j <= segment_index[i]; j++)
			{
				tmp_pl.push_back(newplist[j]);
			}
			//更新节点的坐标, 第0个不用更新
			sgraph.nodeList[apath[i]].position = newplist[segment_index[i]];
			//替换polyline
			substituteEdge(tmp_node->curveType, tmp_node->curveIndex, tmp_pl);
		}
	}
	else
	{
		mypoint2f addp1(0, 0);//计算增量
		mypoint2f addp2(0, 0);
		for (unsigned int i = 0; i < left_p->size(); ++i)
		{
			if (left_p->at(i) == startp)
			{
				mypoint2f origip = sgraph.nodeList[startp].position;
				addp1 = leftp_project->at(i) - origip;
			}
			else if (left_p->at(i) == apath.back())
			{
				mypoint2f origip = sgraph.nodeList[apath.back()].position;
				addp2 = leftp_project->at(i) - origip;
			}
		}
		//线性递减，x=0,y=addp1  x=1.y=addp2; y=(d2-d1)x+d1;
		mypoint2f a = addp2 - addp1;
		mypoint2f b = addp1;
		double step = 1.0 / (origi_plist.size() - 1);
		for (unsigned int i = 0; i < origi_plist.size(); ++i)
		{
			double tmp = step*(double)i;
			mypoint2f tmp_add(a.x*tmp + b.x, a.y*tmp + b.y);
			mypoint2f tmp_p = tmp_add + origi_plist[i];
			newplist.push_back(tmp_p);
		}
		//分段更新polyline  segment_index.size()==apathy.size() ?
		for (unsigned int i = 1; i < segment_index.size(); ++i)
		{
			//得到arcnode
			ArcNode* tmp_node = getArcNode(apath[i - 1], apath[i]);
			//截取新的polyline
			vector<mypoint2f> tmp_pl;
			for (int j = segment_index[i - 1]; j <= segment_index[i]; j++)
			{
				tmp_pl.push_back(newplist[j]);
			}
			//更新节点的坐标, 第0个和最后一个不用更新
			if (apath[i] != startp && apath[i]!=apath.back())
				sgraph.nodeList[apath[i]].position = newplist[segment_index[i]];
			//替换polyline
			substituteEdge(tmp_node->curveType, tmp_node->curveIndex, tmp_pl);
		}
	}
	
}

//尝试 局部操作，先计算每个patch的weight（跟neighbor有关），然后根据weight从小到大，一个一个处理
void MyGraph::localProcess() //localProcess
{
	//1、圆度，公式4*pi*area/permeter^2. 如果是圆，改特征值=1
	//2、长度(elongation)、(离心率eccentricity)、length of major axis/length of minor axis	
	cout << "localProcess()" << endl;
	vector<double> cir_vec;
	vector<double> long_vec;
	vector<double> area_vec;
	int pchsize = plrGraph.patches.size();
	for (int i = 0; i < pchsize; ++i)//付初值 vector<pair<int,vector<double>>> patchfeatures;在.h文件中//int 对应patchindex，vector<double>对应特征向量
	{
		vector<double> tmp;
		patchfeatures.push_back(make_pair(0, tmp));
	}
	for (int i = 0; i < pchsize; ++i)
	{
		if (plrGraph.patches[i].type >= 0)
		{
			vector<mypoint2f> total_plist;
			jointMethod(&plrGraph.patches[i].pointIndexOfPatch, &total_plist, 0);//	total_plist.pop_back();

			vector<double> fvec;//特征向量
			getThreeFeature(&fvec, &total_plist, plrGraph.patches[i].area, plrGraph.patches[i].perimeter);
			patchfeatures[i].first = i;
			patchfeatures[i].second = fvec;
		}
	}
	vector<pair<int, double>> tt_w;//所有patch的权值,各特征在局部所占比例的和
	for (int i = 0; i < pchsize; ++i)
	{
		tt_w.push_back(make_pair(i, 0));
		if (plrGraph.patches[i].type >= 0)
		{
			vector<int> adjpch = getNeighborPatches(i,1);
			if (!adjpch.empty())
			{
				//计算局部总
				double local_cir = 0;
				double local_long = 0;
				double local_area = 0;	
				for (unsigned int j = 0; j < adjpch.size(); ++j)
				{
					local_cir = local_cir + patchfeatures[adjpch[j]].second[0];
					local_long = local_long + patchfeatures[adjpch[j]].second[1];
					local_area = local_area + patchfeatures[adjpch[j]].second[2];
				}
				//邻域patch不包括自己，把本patch的特征加进去
				local_cir = local_cir + patchfeatures[i].second[0];
				local_long = local_long + patchfeatures[i].second[1];
				local_area = local_area + patchfeatures[i].second[2];
				//计算各权值(圆度+长度+面积之比)		
				double weight = patchfeatures[i].second[0] / local_cir + patchfeatures[i].second[1] / local_long + patchfeatures[i].second[2] / local_area;
				tt_w[i].second = weight;//tt_w.push_back(make_pair(i, weight));
			}
			else
			{
				//单个patch，没有neighbor,权值=各特征相加
				tt_w[i].second = patchfeatures[i].second[0] + patchfeatures[i].second[1] + patchfeatures[i].second[2];
			}
		}
	}
	sort(tt_w.begin(), tt_w.end(), mycompare);//从小到大排序
	for (int i = 0; i < tt_w.size(); ++i)// tt_w.size()  pchsize 31
	{
		int pindex = tt_w[i].first;
		if (plrGraph.patches[pindex].type >= 0)
		{
			/*if (tt_w[i].second > 0 && tt_w[i].second < 0.2)
			{
				plrGraph.patches[pindex].type = 1;
			}
			else if (tt_w[i].second >= 0.2 && tt_w[i].second < 0.4)
			{
				plrGraph.patches[pindex].type = 2;
			}
			else if (tt_w[i].second >= 0.4 && tt_w[i].second < 0.6)
			{
				plrGraph.patches[pindex].type = 3;
			}
			else if (tt_w[i].second >= 0.6 && tt_w[i].second < 0.8)
			{
				plrGraph.patches[pindex].type =4;
			}
			else if (tt_w[i].second >= 0.8 && tt_w[i].second <1.0)
			{
				plrGraph.patches[pindex].type = 5;
			}
			else if (tt_w[i].second >= 1.0 && tt_w[i].second <2)
			{
				plrGraph.patches[pindex].type = 6;
			}
			else
			{
				plrGraph.patches[pindex].type = 0;
			}*/
			
			//判断拐点个数
			vector<int> turn_point;//记录拐点 拐点一定在point index of patch上吗????????????
			vector<mypoint2f> apatch;
			for (unsigned int j = 1; j < plrGraph.patches[pindex].pointIndexOfPatch.size(); ++j)
			{
				int verindex1 = plrGraph.patches[pindex].pointIndexOfPatch.at(j - 1);
				int verindex2 = plrGraph.patches[pindex].pointIndexOfPatch.at(j);
				mypoint2f point1 = sgraph.nodeList[verindex1].position;//=lineWithCommonPoints[verindex1].first
				mypoint2f point2 = sgraph.nodeList[verindex2].position;
				if (!apatch.empty())
					apatch.pop_back();
				turn_point.push_back((int)apatch.size());  //pointindex对应plist中的位置
				getPListOfaCurve(point1, point2, verindex1, verindex2, &apatch,0);//最终得到的apatch 首=尾
			}
			mypoint2f pfront = apatch.front();
			mypoint2f pback = apatch[apatch.size()-2];
			apatch.insert(apatch.begin(), pback);
		
			int turn_count = 0;
			for (unsigned int j = 1; j < apatch.size() - 1; ++j)
			{
				mypoint2f vec1 = apatch[j] - apatch[j - 1];
				mypoint2f vec2 = apatch[j + 1] - apatch[j];
				double cos_theta = getCosineAngle(apatch[j - 1], apatch[j], apatch[j], apatch[j + 1]);
				if (cos_theta >= -1 && cos_theta <= 0.707106)//转角在45到180认为是一个拐点   135-->-0.7071  150--> -0.866
				{
					turn_count++; //cout << cos_theta << endl;
					for (unsigned int k = 0; k < turn_point.size(); ++k)
					{
						if (turn_point[k] == j - 1)
						{
							cout << k << "," << plrGraph.patches[pindex] .pointIndexOfPatch[k]<< endl;
							break;
						}
					}
				}
			}
			if (turn_count < 3)//有两个拐点的 -->erosion
			{
				//查看是不是neighbor中最大weight
				vector<int> adjpch = getNeighborPatches(pindex,1);
				bool maxornot = true;
				if (!adjpch.empty())
				{
					for (unsigned int j = 0; j < adjpch.size(); ++j)
					{
						if (tt_w[pindex].second < tt_w[adjpch[j]].second)
						{
							maxornot = false;
							break;
						}
					}
					if (maxornot == false)
					{
						erosionApatch(pindex);
					}
					
				}
				//else
					//return;
			}
			else if (turn_count<=4)//>2拐点的combine
			{
				//找到最短边，
				vector<double> edge_per;
				double mini_edge = 655350;//0
				int record_j = -1;
				for (unsigned int j = 1; j < plrGraph.patches[pindex].pointIndexOfPatch.size(); ++j)
				{
					//if (plrGraph.patches[i].type >= 0)
					//{
					vector<mypoint2f> edgeplist;
					int verindex1 = plrGraph.patches[pindex].pointIndexOfPatch.at(j - 1);
					int verindex2 = plrGraph.patches[pindex].pointIndexOfPatch.at(j);
					mypoint2f point1 = sgraph.nodeList[verindex1].position;//=lineWithCommonPoints[verindex1].first
					mypoint2f point2 = sgraph.nodeList[verindex2].position;
					getPListOfaCurve(point1, point2, verindex1, verindex2, &edgeplist,0);
					double anEdge_p = getPerimeterofVec(&edgeplist);
					int pchcount = getBorderofPatch(pindex, verindex1, verindex2)->adjPatchIndex.size();
					if (anEdge_p < mini_edge && pchcount>1) //最duan 且这条边有两个patch
					{
						mini_edge = anEdge_p;
						record_j = j;
					}
					edge_per.push_back(anEdge_p);
					//}				
				}
				if (record_j > 0)
				{
					//plrGraph.patches[i].edges_vec[record_j - 1].edgeflag = 1;//edge标记为1，表示最短/最长边		
					int adj_pindex = plrGraph.patches[pindex].edges_vec[record_j - 1].adjPatchIndex[1];
					combineTwoPatch(adj_pindex, pindex);
				}
			}
		}
	}

	//vector<double> tt_w(plrGraph.patches.size(),0);//所有patch的权值,各特征在局部所占比例的和
	//for (int i = 0; i < pchsize; ++i)
	//{
	//	if (plrGraph.patches[i].type >= 0)
	//	{
	//		vector<int> adjpch = getNeighborPatches(i);
	//		if (!adjpch.empty())
	//		{
	//			//计算局部总面
	//			double local_cir = 0;
	//			double local_long = 0;
	//			double local_area = 0;
	//			for (unsigned int j = 0; j < adjpch.size(); ++j)
	//			{
	//				local_cir = local_cir + patchfeatures[adjpch[j]].second[0];
	//				local_long = local_long + patchfeatures[adjpch[j]].second[1];
	//				local_area = local_area + patchfeatures[adjpch[j]].second[2];
	//			}
	//			//邻域patch不包括自己，把本patch的特征加进去
	//			local_cir = local_cir + patchfeatures[i].second[0];
	//			local_long = local_long + patchfeatures[i].second[1];
	//			local_area = local_area + patchfeatures[i].second[2];
	//			//计算各权值(圆度+长度+面积之比)		
	//			double weight = patchfeatures[i].second[0] / local_cir + patchfeatures[i].second[1] / local_long + patchfeatures[i].second[2] / local_area;
	//			tt_w[i]=weight;				
	//		}
	//	}
	//}
	//sort(tt_w.begin(), tt_w.end());//从小到大排序
	//for (int i = 0; i < 5; ++i)//pchsize 31
	//{
	//	if (plrGraph.patches[i].type >= 0)
	//	{
	//		//判断拐点个数
	//		vector<mypoint2f> apatch;
	//		for (unsigned int j = 1; j < plrGraph.patches[i].pointIndexOfPatch.size(); ++j)
	//		{
	//			int verindex1 = plrGraph.patches[i].pointIndexOfPatch.at(j - 1);
	//			int verindex2 = plrGraph.patches[i].pointIndexOfPatch.at(j);
	//			mypoint2f point1 = sgraph.nodeList[verindex1].position;//=lineWithCommonPoints[verindex1].first
	//			mypoint2f point2 = sgraph.nodeList[verindex2].position;
	//			if (!apatch.empty())
	//				apatch.pop_back();
	//			getPListOfaCurve(point1, point2, verindex1, verindex2, &apatch);//最终得到的apatch 首=尾
	//		}
	//		if (isPointRoughEqual(apatch.front(), apatch.back()))
	//			apatch.pop_back();
	//		mypoint2f pfront = apatch.front();
	//		mypoint2f pback = apatch.back();
	//		apatch.insert(apatch.begin(), pback);
	//		apatch.push_back(pfront);
	//		int turn_count = 0;
	//		for (unsigned int j = 1; j < apatch.size() - 1; ++j)
	//		{
	//			mypoint2f vec1 = apatch[j] - apatch[j - 1];
	//			mypoint2f vec2 = apatch[j + 1] - apatch[j];
	//			double cos_theta = getCosineAngle(apatch[j - 1], apatch[j], apatch[j], apatch[j + 1]);
	//			if (cos_theta >= -1 && cos_theta <= 0.707106)//转角在45到180认为是一个拐点   135-->-0.7071  150--> -0.866
	//			{
	//				turn_count++; //cout << cos_theta << endl;
	//			}
	//		}
	//		if (turn_count < 3)
	//		{
	//			//查看是不是neighbor中最大weight
	//			vector<int> adjpch = getNeighborPatches(i);
	//			bool maxornot = true;
	//			if (!adjpch.empty())
	//			{
	//				for (unsigned int j = 0; j < adjpch.size(); ++j)
	//				{
	//					if (tt_w[i] < tt_w[adjpch[j]])
	//					{
	//						maxornot = false;
	//						break;
	//					}
	//				}
	//				if (maxornot == false)
	//					erosionApatch(i);
	//			}
	//			else
	//				return;
	//		}
	//		else
	//		{
	//			//找到最短边，
	//			vector<double> edge_per;
	//			double mini_edge = 655350;//0
	//			int record_j = -1;
	//			for (unsigned int j = 1; j < plrGraph.patches[i].pointIndexOfPatch.size(); ++j)
	//			{
	//				//if (plrGraph.patches[i].type >= 0)
	//				//{
	//				vector<mypoint2f> edgeplist;
	//				int verindex1 = plrGraph.patches[i].pointIndexOfPatch.at(j - 1);
	//				int verindex2 = plrGraph.patches[i].pointIndexOfPatch.at(j);
	//				mypoint2f point1 = sgraph.nodeList[verindex1].position;//=lineWithCommonPoints[verindex1].first
	//				mypoint2f point2 = sgraph.nodeList[verindex2].position;
	//				getPListOfaCurve(point1, point2, verindex1, verindex2, &edgeplist);
	//				double anEdge_p = getPerimeterofVec(edgeplist);
	//				int pchcount = getBorderofPatch(i, verindex1, verindex2)->adjPatchIndex.size();
	//				if (anEdge_p < mini_edge && pchcount>1) //最duan 且这条边有两个patch
	//				{
	//					mini_edge = anEdge_p;
	//					record_j = j;
	//				}
	//				edge_per.push_back(anEdge_p);
	//				//}				
	//			}
	//			if (record_j > 0)
	//			{
	//				//plrGraph.patches[i].edges_vec[record_j - 1].edgeflag = 1;//edge标记为1，表示最短/最长边		
	//				int adj_pindex=plrGraph.patches[i].edges_vec[record_j - 1].adjPatchIndex[1];
	//				combineTwoPatch(adj_pindex, i);
	//			}
	//		}
	//	}
	//}

	/*vector<int> adjpch = getNeighborPatches(0);
	combineTwoPatch(0,adjpch[0]);
	erosionApatch(3);*/
}



//计算相似度(尝试)
void MyGraph::patchSimilarity()//shan
{
	int simply_num = 64;
	int pchsize = plrGraph.patches.size();
	vector<vector<double>> descriptor_An;

	for (int i = 0; i < pchsize; ++i)
	{
		if (plrGraph.patches[i].type >= 0)
		{
			vector<mypoint2f> total_plist;
			jointMethod(&plrGraph.patches[i].pointIndexOfPatch, &total_plist, 0);// total_plist首=尾
			total_plist.pop_back();
			vector<mypoint2f> evol_plist;
			if (total_plist.size()>simply_num)// >128个点，则需要化简
				dataReduction(&total_plist, &plrGraph.patches[i].pointIndexOfPatch, i, simply_num, &evol_plist);//evol_plist是简化后的point list,front()!=back()
			else//否则就添加点，再化简
			{
				total_plist.push_back(total_plist.front());
				dataIncreasing(&total_plist, simply_num);
				total_plist.pop_back();
				dataReduction(&total_plist, &plrGraph.patches[i].pointIndexOfPatch, i, simply_num, &evol_plist);
			}	
			pointListOfPatches.push_back(make_pair(plrGraph.patches[i].type, evol_plist));//用于patch的显示
			/*方法一,*/
			////计算bend angle function
			//mypoint2f frontp = evol_plist.front();//添加一个头，使evol_plist 首=尾
			//evol_plist.push_back(frontp);
			//vector<vector<double>> bendAngle_vec;//[x,y] x累计边长， y旋转角
			//getbendAngleFunction(&evol_plist, &bendAngle_vec, simply_num * 2);//先在简化的contour的基础上计算，再等距采样的 
			////计算fourier coefficient
			//vector<double> An_pch;
			//for (unsigned int j = 0; j < bendAngle_vec.size(); ++j)
			//{
			//	double an_sum = 0;
			//	double bn_sum = 0;
			//	for (unsigned int k = 0; k < bendAngle_vec.size(); ++k)
			//	{
			//		double t_ff = 2 * 3.1415926*(j + 1)*bendAngle_vec[k][0];
			//		double ttmp=bendAngle_vec[k][1] * sin(t_ff);
			//		an_sum = an_sum + ttmp;
			//		ttmp = bendAngle_vec[k][1] * cos(t_ff);
			//		bn_sum = bn_sum + ttmp;
			//	}
			//	an_sum = an_sum / ((j + 1)*(-3.1415926));
			//	bn_sum = bn_sum / ((j + 1)*3.1415926);
			//	double tt = sqrt(pow(an_sum, 2) + pow(bn_sum, 2));
			//	An_pch.push_back(tt);
			//}
			//descriptor_An.push_back(An_pch);
			//ofstream outfile1("E:/aaa.txt");
			//ofstream outfile2("E:/bbb.txt");
			//if (outfile1.is_open() && outfile2.is_open())
			//{
			//	for (unsigned int j = 0; j < bendAngle_vec.size(); ++j)//单位化
			//	{
			//		bendAngle_vec[j][0] = bendAngle_vec[j][0] / cumulate_len;
			//		outfile1 << bendAngle_vec[j][0];
			//		outfile1 << " ";
			//		outfile2 << bendAngle_vec[j][1];
			//		outfile2 << " ";
			//	}
			//}
			//outfile1.close();
			//outfile2.close();
			/*方法二*/
			vector<complex<double>> FD;
			mypoint2f frontp = evol_plist.front();//添加一个头，使evol_plist 首=尾
			evol_plist.push_back(frontp);
			getFourierDescriptor(&evol_plist, &FD, i);//用radial destance计算
			vector<double> norm_fd;
			double norm_len = sqrt(pow(FD[0].real(), 2) + pow(FD[0].imag(), 2));
			int half_size = FD.size() / 2+1;
			for (int j = 1; j < half_size; ++j)//从1开始到N/2！！！
			{
				double fd_len = sqrt(pow(FD[j].real(), 2) + pow(FD[j].imag(), 2)) / norm_len;
				norm_fd.push_back(fd_len);
			}
			//getFDFeature(&total_plist, &plrGraph.patches[i].pointIndexOfPatch, i, simply_num, 0, &norm_fd);
			descriptor_An.push_back(norm_fd);
			///*The lower frequency descriptors contain information about the general shape,
			//and the higher frequency descriptors contain information about
			//smaller details. For classification purposes, a subset of the Fourier
			//descriptors is often enough to discriminate different shapes.*/
		}
	}
	if (descriptor_An.size() > 1)//至少有2个patch才能比较
	{
		for (unsigned int i = 0; i < descriptor_An.size(); ++i)
		{	
			vector<double> d_1 = descriptor_An[i];
			for (unsigned int j = i + 1; j < descriptor_An.size(); ++j)
			{
				double diff = 0;
				vector<double> d_2 = descriptor_An[j];
				for (unsigned int k = 0; k < d_1.size(); ++k)
				{
					diff = diff + pow((d_1[k] - d_2[k]), 2);
				}
				diff = sqrt(diff);
				cout << i << "," << j << " similarity is:" << diff << endl;
			}
		}
	}
	cout << endl;
}
//数据集的化简，使其有固定的fixed point number？ 并且保留特征点（拐点）
void MyGraph::dataReduction(vector<mypoint2f> *plist, vector<int>* p_vec, int pchindex, int fixednum, vector<mypoint2f> *adjust_plist)
{
	//plist是一定有的，p_vec为了间隔取点（为NULL，说明是一条newline），pchindex为了获得patch周长（<0要重心计算周长）
	//plist的front() != back()，化简得evolu_plist, pchindex有可能<0,有了pchindex可以知道perimeter，没有pchindex,周长通过contour计算
	if (plist->front()==plist->back())
	{
		plist->pop_back();
	}
	//1.若没有perimeter，先计算一下
	double perimeter = 0;
	if (pchindex < 0 || pchindex>=plrGraph.patches.size())
	{
		mypoint2f frontp = plist->front();
		plist->push_back(frontp);
		perimeter = getPerimeterofVec(plist);
		plist->pop_back();
	}
	else
	{
		perimeter = plrGraph.patches[pchindex].perimeter;
	}
	//2.若contour plist的数量等于fixednum，
	if (plist->size() <= fixednum)
	{
		for (unsigned int i = 0; i < plist->size(); ++i)
		{
			adjust_plist->push_back(plist->at(i));
		}
	}
	else
	{
		if (plist->size() > fixednum * 4)
		{	//3.若数量很大(大于fixednum的三倍)，先间隔取点，再用堆删除
			if (p_vec->empty())//为NULL，说明是一条newline
			{
				int first_seg = plist->size() / 2;//这是一条线且front!=back。eg:12345432 (5+3)/2=4,index=4个点是一个拐点。
				int tmp = fixednum ;//每条边化简后的点数*2
				int interval = (first_seg - 1) / tmp;
				vector<mypoint2f> mid_plist;
				for (int i = 0; i < first_seg; i += interval)
				{
					mid_plist.push_back(plist->at(i));
				}
				if (mid_plist.back() == plist->at(first_seg))
					mid_plist.pop_back();
				for (int i = first_seg; i < plist->size(); i += interval)
				{
					mid_plist.push_back(plist->at(i));
				}
				dataReduction_byheap(&mid_plist, fixednum, perimeter, adjust_plist);
			}
			else//patch or line
			{
				int ptsize = p_vec->size();
				vector<mypoint2f> mid_plist;
				for (int i = 1; i < ptsize; ++i)
				{
					vector<int> anedge;
					vector<mypoint2f> edgeplist;
					anedge.push_back(p_vec->at(i - 1));
					anedge.push_back(p_vec->at(i));
					jointMethod(&anedge, &edgeplist, 0);
					double edge_peri = getPerimeterofVec(&edgeplist);
					int tmp = (edge_peri / perimeter)* (fixednum * 2);//每条边化简后的点数*2
					if (tmp <=2)
						tmp = edgeplist.size() - 1;
					int	interval = (edgeplist.size() - 1) / tmp;//每隔interval个取一个点，两段之间共取tmp个点
					if (interval <= 0)
						interval = 1;
					for (unsigned int j = 0; j < edgeplist.size(); j += interval)//取样点会比实际的多,interval=2.3、2.8-->2
					{
						mid_plist.push_back(edgeplist[j]);
					}
					if (mid_plist.back() == edgeplist.back())
						mid_plist.pop_back();
				}
				dataReduction_byheap(&mid_plist, fixednum, perimeter, adjust_plist);
			}		
		}
		else
		{
			//4.数量不多，直接用堆删除。
			dataReduction_byheap(plist, fixednum, perimeter, adjust_plist);
		}	
	}
	/*原来直接用堆删除的方法*/
	//double perimeter = 0;
	//if (pchindex < 0)
	//{
	//	mypoint2f frontp = plist->front();
	//	plist->push_back(frontp);
	//	perimeter = getPerimeterofVec(plist);
	//	plist->pop_back();
	//}
	//else
	//{
	//	perimeter = plrGraph.patches[pchindex].perimeter;
	//}
	//vector<double> kk;// relevance measure
	//double costheta = getCosineAngle(plist->at(0), plist->back(), plist->at(0), plist->at(1));//plist->back() != plist->at(0)
	//if (costheta > 1)
	//	costheta = 1;
	//double agl = acos(costheta) * 180 / 3.1415926;
	//double l1 = LineSegmentLength(plist->back(), plist->at(0)) / perimeter;
	//double l2 = LineSegmentLength(plist->at(0), plist->at(1)) / perimeter;
	//double relevance = abs(agl - 180)*l1*l2 / (l1 + l2);
	//kk.push_back(relevance);
	//for (unsigned int i = 1; i < plist->size() - 1; ++i)
	//{
	//	costheta = getCosineAngle(plist->at(i), plist->at(i - 1), plist->at(i), plist->at(i + 1));
	//	if (costheta > 1)
	//		costheta = 1;
	//	agl = acos(costheta) * 180 / 3.1415926;
	//	l1 = LineSegmentLength(plist->at(i - 1), plist->at(i)) / perimeter;
	//	l2 = LineSegmentLength(plist->at(i), plist->at(i + 1)) / perimeter;
	//	relevance = abs(agl - 180)*l1*l2 / (l1 + l2);
	//	kk.push_back(relevance);
	//}
	//costheta = getCosineAngle(plist->back(), plist->at(plist->size() - 2), plist->back(), plist->at(0));
	//if (costheta > 1)
	//	costheta = 1;
	//agl = acos(costheta) * 180 / 3.1415926;
	//l1 = LineSegmentLength(plist->at(plist->size() - 2), plist->back()) / perimeter;
	//l2 = LineSegmentLength(plist->back(), plist->at(0)) / perimeter;
	//relevance = abs(agl - 180)*l1*l2 / (l1 + l2);
	//kk.push_back(relevance);
	////建堆，kk保存了从第0个点到底n个点relevance值，堆heap_vec保存relevance在kk中的下标
	//vector<int> heap_vec;
	//vector<int> ele_index;
	//for (unsigned int i = 0; i < kk.size(); ++i)
	//{
	//	heap_vec.push_back(i);
	//	ele_index.push_back(0);
	//}
	////建立左右neighbor关系，
	//vector<vector<int>> l_nei;
	//vector<int> vec_tmp;  vec_tmp.push_back(kk.size() - 1); vec_tmp.push_back(1);
	//l_nei.push_back(vec_tmp);
	//for (unsigned int i = 0; i < kk.size() - 2; ++i)
	//{
	//	vector<int> vv; vv.push_back(i); vv.push_back(i + 2);
	//	l_nei.push_back(vv);
	//}
	//vec_tmp[0] = kk.size() - 2; vec_tmp[1] = 0;
	//l_nei.push_back(vec_tmp);
	////建堆
	//int h_size = heap_vec.size();
	//for (int i = h_size / 2; i >= 1; i--)
	//{
	//	heap_siftdown(i, &heap_vec, &kk, &ele_index);
	//}
	////建立索引，知道某个relevanc的下标，可以直接找到在heap的位置，位置是从1开始，注意更新
	//for (int i = 0; i < h_size; ++i)
	//{
	//	ele_index[heap_vec[i]] = i;// ele_index
	//}
	////每次删除最小值，并且更新左右邻域的relevance值 并且调整其在堆中的位置
	//while (h_size > fixednum)
	//{
	//	int tmp_back = heap_vec.back();//删除最小的,即最后一个点与第一个交换位置，然后进行下沉调整
	//	int tmp_mini = heap_vec[0];
	//	kk[tmp_mini] = -1;
	//	ele_index[tmp_mini] = -1;
	//	heap_vec[0] = tmp_back;
	//	ele_index[tmp_back] = 0;//!!
	//	heap_siftdown(1, &heap_vec, &kk, &ele_index);//下沉调整
	//	heap_vec.pop_back();
	//	h_size--;
	//	//修改左右节点，并更新调整
	//	int n_left = l_nei[tmp_mini][0];
	//	int n_right = l_nei[tmp_mini][1];
	//	l_nei[n_left][1] = n_right;
	//	l_nei[n_right][0] = n_left;
	//	//更新kk
	//	costheta = getCosineAngle(plist->at(n_left), plist->at(l_nei[n_left][0]), plist->at(n_left), plist->at(l_nei[n_left][1]));
	//	if (costheta > 1)
	//		costheta = 1;
	//	agl = acos(costheta) * 180 / 3.1415926;
	//	l1 = LineSegmentLength(plist->at(l_nei[n_left][0]), plist->at(n_left)) / perimeter;
	//	l2 = LineSegmentLength(plist->at(n_left), plist->at(l_nei[n_left][1])) / perimeter;
	//	relevance = abs(agl - 180)*l1*l2 / (l1 + l2);
	//	//kk[n_left] = relevance;
	//	heap_update(ele_index[n_left] + 1, relevance, &heap_vec, &kk, &ele_index);//注意+1
	//	costheta = getCosineAngle(plist->at(n_right), plist->at(l_nei[n_right][0]), plist->at(n_right), plist->at(l_nei[n_right][1]));
	//	if (costheta > 1)
	//		costheta = 1;
	//	agl = acos(costheta) * 180 / 3.1415926;
	//	l1 = LineSegmentLength(plist->at(l_nei[n_right][0]), plist->at(n_right)) / perimeter;
	//	l2 = LineSegmentLength(plist->at(n_right), plist->at(l_nei[n_right][1])) / perimeter;
	//	relevance = abs(agl - 180)*l1*l2 / (l1 + l2);
	//	heap_update(ele_index[n_right] + 1, relevance, &heap_vec, &kk, &ele_index);//kk[n_right] = relevance;
	//}
	////kk中标记为-1的表示已删除，把剩下的存入 参数adjust_plist中，
	//for (unsigned int i = 0; i < kk.size(); ++i)
	//{
	//	if (kk[i] != -1)
	//	{
	//		adjust_plist->push_back(plist->at(i));
	//	}
	//}
}
void MyGraph::dataReduction_byheap(vector<mypoint2f>* plist, int fixednum, double perimeter, vector<mypoint2f> *adjust_plist)
{
	/*原来的方法，直接用堆删除的方法*/
	vector<double> kk;// relevance measure
	double costheta = getCosineAngle(plist->at(0), plist->back(), plist->at(0), plist->at(1));//plist->back() != plist->at(0)
	if (costheta > 1)
		costheta = 1;
	double agl = acos(costheta) * 180 / 3.1415926;
	double l1 = LineSegmentLength(plist->back(), plist->at(0)) / perimeter;
	double l2 = LineSegmentLength(plist->at(0), plist->at(1)) / perimeter;
	double relevance = abs(agl - 180)*l1*l2 / (l1 + l2);
	kk.push_back(relevance);
	for (unsigned int i = 1; i < plist->size() - 1; ++i)
	{
		costheta = getCosineAngle(plist->at(i), plist->at(i - 1), plist->at(i), plist->at(i + 1));
		agl = acos(costheta) * 180 / 3.1415926;
		l1 = LineSegmentLength(plist->at(i - 1), plist->at(i)) / perimeter;
		l2 = LineSegmentLength(plist->at(i), plist->at(i + 1)) / perimeter;
		relevance = abs(agl - 180)*l1*l2 / (l1 + l2);
		kk.push_back(relevance);
	}
	costheta = getCosineAngle(plist->back(), plist->at(plist->size() - 2), plist->back(), plist->at(0));
	agl = acos(costheta) * 180 / 3.1415926;
	l1 = LineSegmentLength(plist->at(plist->size() - 2), plist->back()) / perimeter;
	l2 = LineSegmentLength(plist->back(), plist->at(0)) / perimeter;
	relevance = abs(agl - 180)*l1*l2 / (l1 + l2);
	kk.push_back(relevance);
	//建堆，kk保存了从第0个点到底n个点relevance值，堆heap_vec保存relevance在kk中的下标
	vector<int> heap_vec;
	vector<int> ele_index;
	for (unsigned int i = 0; i < kk.size(); ++i)
	{
		heap_vec.push_back(i);
		ele_index.push_back(0);
	}
	//建立左右neighbor关系，
	vector<vector<int>> l_nei;
	vector<int> vec_tmp;  vec_tmp.push_back(kk.size() - 1); vec_tmp.push_back(1);
	l_nei.push_back(vec_tmp);
	for (unsigned int i = 0; i < kk.size() - 2; ++i)
	{
		vector<int> vv; vv.push_back(i); vv.push_back(i + 2);
		l_nei.push_back(vv);
	}
	vec_tmp[0] = kk.size() - 2; vec_tmp[1] = 0;
	l_nei.push_back(vec_tmp);
	//建堆
	int h_size = heap_vec.size();
	for (int i = h_size / 2; i >= 1; i--)
	{
		heap_siftdown(i, &heap_vec, &kk, &ele_index);
	}
	//建立索引，知道某个relevanc的下标，可以直接找到在heap的位置，位置是从1开始，注意更新
	for (int i = 0; i < h_size; ++i)
	{
		ele_index[heap_vec[i]] = i;// ele_index
	}
	//每次删除最小值，并且更新左右邻域的relevance值 并且调整其在堆中的位置
	while (h_size > fixednum)
	{
		int tmp_back = heap_vec.back();//删除最小的,即最后一个点与第一个交换位置，然后进行下沉调整
		int tmp_mini = heap_vec[0];
		kk[tmp_mini] = -1;
		ele_index[tmp_mini] = -1;
		heap_vec[0] = tmp_back;
		ele_index[tmp_back] = 0;//!!
		heap_siftdown(1, &heap_vec, &kk, &ele_index);//下沉调整
		heap_vec.pop_back();
		h_size--;
		//修改左右节点，并更新调整
		int n_left = l_nei[tmp_mini][0];
		int n_right = l_nei[tmp_mini][1];
		l_nei[n_left][1] = n_right;
		l_nei[n_right][0] = n_left;
		//更新kk
		costheta = getCosineAngle(plist->at(n_left), plist->at(l_nei[n_left][0]), plist->at(n_left), plist->at(l_nei[n_left][1]));
		agl = acos(costheta) * 180 / 3.1415926;
		l1 = LineSegmentLength(plist->at(l_nei[n_left][0]), plist->at(n_left)) / perimeter;
		l2 = LineSegmentLength(plist->at(n_left), plist->at(l_nei[n_left][1])) / perimeter;
		relevance = abs(agl - 180)*l1*l2 / (l1 + l2);
		//kk[n_left] = relevance;
		heap_update(ele_index[n_left] + 1, relevance, &heap_vec, &kk, &ele_index);//注意+1
		costheta = getCosineAngle(plist->at(n_right), plist->at(l_nei[n_right][0]), plist->at(n_right), plist->at(l_nei[n_right][1]));
		agl = acos(costheta) * 180 / 3.1415926;
		l1 = LineSegmentLength(plist->at(l_nei[n_right][0]), plist->at(n_right)) / perimeter;
		l2 = LineSegmentLength(plist->at(n_right), plist->at(l_nei[n_right][1])) / perimeter;
		relevance = abs(agl - 180)*l1*l2 / (l1 + l2);
		heap_update(ele_index[n_right] + 1, relevance, &heap_vec, &kk, &ele_index);//kk[n_right] = relevance;
	}
	//kk中标记为-1的表示已删除，把剩下的存入 参数adjust_plist中，
	for (unsigned int i = 0; i < kk.size(); ++i)
	{
		if (kk[i] != -1)
		{
			adjust_plist->push_back(plist->at(i));
		}
	}
}
void MyGraph::equalSampling(vector<mypoint2f> *plist, vector<int>* p_vec, int pchindex, int fixednum, vector<mypoint2f> *adjust_plist)
{
	//plist是一定有的，p_vec为了间隔取点（为NULL，说明是一条newline），pchindex为了获得patch周长（<0要重心计算周长）
	//plist的front() != back()，化简得evolu_plist, pchindex有可能<0,有了pchindex可以知道perimeter，没有pchindex,周长通过contour计算
	if (plist->front() == plist->back())
	{
		plist->pop_back();
	}
	//1.若没有perimeter，先计算一下
	double perimeter = 0;
	if (pchindex < 0 || pchindex>=plrGraph.patches.size())
	{
		mypoint2f frontp = plist->front();
		plist->push_back(frontp);
		perimeter = getPerimeterofVec(plist);
		plist->pop_back();
	}
	else
	{ 
		if (p_vec->empty())
		{
			mypoint2f frontp = plist->front();
			plist->push_back(frontp);
			perimeter = getPerimeterofVec(plist);
			plist->pop_back();
		}
		else
		{
			perimeter = plrGraph.patches[pchindex].perimeter;
		}
	}
	//2.若contour plist的数量等于fixednum，
	if (plist->size() == fixednum)
	{
		for (unsigned int i = 0; i < plist->size(); ++i)
		{
			adjust_plist->push_back(plist->at(i));
		}
	}
	else
	{
		vector<int> edge_count;   //存储edge简化后的点数
		vector<double> edge_peri; //存储每条edge 的长度
		vector<vector<mypoint2f>> edge_plist;  //存储每条edge 中实际的点plist
		int pvecsize =0;
		if (p_vec->empty())//为NULL，说明是一条newline 或者 变形patch（patch+line）
		{
			if (pchindex < 0)// l+l newline
			{
				//vector<int> edge_count;   //存储edge简化后的点数
				//vector<double> edge_peri; //存储每条edge 的长度
				//vector<vector<mypoint2f>> edge_plist;
				pvecsize = 3;
				int first_seg = plist->size() / 2;//这是一条线且front!=back。eg:12345432 (5+3)/2=4,index=4个点是一个拐点。
				int tmp = fixednum * 0.5;//每条边化简后的点数(大概的)（每边包括两个端点）
				edge_count.push_back(tmp);	edge_count.push_back(tmp);
				edge_peri.push_back(perimeter / 2);  edge_peri.push_back(perimeter / 2);

				vector<mypoint2f> mid_plist;
				for (int i = 0; i <= first_seg; ++i)
				{
					mid_plist.push_back(plist->at(i));
				}
				edge_plist.push_back(mid_plist);
				mid_plist.swap(vector<mypoint2f>());
				for (unsigned int i = first_seg; i < plist->size(); ++i)
				{
					mid_plist.push_back(plist->at(i));
				}
				mid_plist.push_back(plist->front());
				edge_plist.push_back(mid_plist);

				int cha = fixednum - tmp * 2 + 2;
				if (cha != 0)
				{
					edge_count[0] = edge_count[0] + cha;
				}
			}
			else// p+l -->变形的transform_patch
			{
				//没有p_vec ，但是有pchindex, plist.front!=back
				//诸葛店判断拐角，先对edge_peri，mid_plist 赋值
				plist->push_back(plist->front());  //front=back !
				vector<mypoint2f> mid_plist;
				for (unsigned int i = 1; i < plist->size()-1; ++i)
				{
					mid_plist.push_back(plist->at(i-1));
					double bend_cos = getCosineAngle(plist->at(i - 1), plist->at(i), plist->at(i), plist->at(i + 1));
					double bend_angle = acos(bend_cos) * 180 / 3.1415926;
					if (bend_angle > 15)
					{
						mid_plist.push_back(plist->at(i));
						edge_plist.push_back(mid_plist);
						double ss_peri= getPerimeterofVec(&mid_plist);

						edge_peri.push_back(ss_peri);
						mid_plist.swap(vector<mypoint2f>());
					}
				}
				mid_plist.push_back(plist->at(plist->size() - 2));
				mid_plist.push_back(plist->back());
				plist->pop_back();  //front!=back
				edge_plist.push_back(mid_plist);
				double ss_peri = getPerimeterofVec(&mid_plist);
				edge_peri.push_back(ss_peri);
				//然后计算edge_count,同时调整edge_plist中的点的个数 ,跟下面重复的代码
				int sim_sum = 0;
				int max_index = 0; int max_c = 0;  //最大plist的edge 的index
				for (unsigned int i = 0; i < edge_plist.size(); ++i)
				{
					//vector<mypoint2f> tmp_plist=edge_plist[i];
					double peri = edge_peri[i];
					int tmp = fixednum * peri / perimeter;//每条边化简后的点数(大概的)（每边包括两个端点）
					if (tmp < 2)
					{
						tmp = 2;
					}
					int ensure_ecount = tmp + 15;  //!!!!!!!!!!!!!
					if (ensure_ecount > edge_plist[i].size())//化简后的point >  原本的point个数
					{
						int csize = edge_plist[i].size();
						int accumulate = csize - 1;
						int j = 1;
						int addnum = 0;
						while ((addnum + csize - 1) < ensure_ecount)//不把末尾点算上
						{
							addnum = accumulate *j;
							j++;
						}
						//每两个点之间插入i-1个点(平均分成i份，即除以i)，共addnum个
						vector<mypoint2f>::iterator iter = edge_plist[i].begin() + 1;
						while (iter != edge_plist[i].end())
						{
							mypoint2f p1 = *(iter - 1);
							mypoint2f p2 = *iter;
							double x_seg = (p2.x - p1.x) / (double)j;
							double y_seg = (p2.y - p1.y) / (double)j;
							for (int k = 1; k < j; ++k)
							{
								mypoint2f midp(p2.x - x_seg*k, p2.y - y_seg*k);
								iter = edge_plist[i].insert(iter, midp);
							}
							iter = iter + j;
						}
					}
					if (tmp > max_c)
					{
						max_c = tmp;
						max_index = i;
					}
					edge_count.push_back(tmp);
					sim_sum = sim_sum + tmp;
				}
				/*把少添加的点平均分到edge中*/
				int cha = fixednum - sim_sum + (int)edge_plist.size();
				if (cha != 0)
				{
					int ttt = cha / (int)edge_count.size();
					if (cha > 0)
					{
						if (ttt <1)
							ttt = 1;
						int sum_ttt = 0;
						int k = 0;
						int k_count = edge_count.size();
						while (sum_ttt < cha && k < k_count)
						{
							if (edge_count[k] != 2)
							{
								if ((edge_count[k] + ttt) <= edge_plist[k].size())
								{
									edge_count[k] = edge_count[k] + ttt;
									sum_ttt = sum_ttt + ttt;
								}
							}
							k++;
						}
						if (sum_ttt < cha)
							edge_count[max_index] = edge_count[max_index] + (cha - sum_ttt);
					}
					else
					{
						if (ttt >-1)
							ttt = -1;
						int sum_ttt = 0;
						int k = 0;
						int k_count = edge_count.size();
						while (sum_ttt > cha && k < k_count)
						{
							if (edge_count[k] != 2)
							{
								if ((edge_count[k] + ttt) >2)
								{
									edge_count[k] = edge_count[k] + ttt;
									sum_ttt = sum_ttt + ttt;
								}
							}
							k++;
						}
						if (sum_ttt > cha)
							edge_count[max_index] = edge_count[max_index] + (cha - sum_ttt);
					}
					//edge_count[max_index] = edge_count[max_index] + cha;
				}

			}	// P+L done	
		}
		else//patch , line，newpatch(p+p)即p_vec不为空
		{
			//vector<int> edge_count;   //存储edge简化后的点数
			//vector<double> edge_peri; //存储每条edge 的长度
			//vector<vector<mypoint2f>> edge_plist;
			/*现在  1.先通过p_vec 把可以连接的edge连接起来，while{}之后得到了edge_peri ,edge_plist。2.之后在通过遍历edge_plist得到edge_count*/
			pvecsize = p_vec->size();			
			int pvec_i = 1;
			while (pvec_i < pvecsize)//根据pointvec分段，但是有的段是可以连起来的
			{
				mypoint2f last_direct(0, 0);
				double bend_angle = 0;
				vector<mypoint2f> connect_plist;
				while (bend_angle < 15 && pvec_i<pvecsize)
				{
					vector<int> anedge;
					vector<mypoint2f> tmp_plist1;
					anedge.push_back(p_vec->at(pvec_i - 1));
					anedge.push_back(p_vec->at(pvec_i));
					jointMethod(&anedge, &tmp_plist1, 0);
					
					if (last_direct == mypoint2f(0, 0))  //第1条edge
					{
						last_direct = tmp_plist1.back() - tmp_plist1[tmp_plist1.size() - 2];
						connect_plist.insert(connect_plist.end(), tmp_plist1.begin(), tmp_plist1.end());
						pvec_i++;
					}
					else
					{
						//mypoint2f curr_direct = tmp_plist1[1] - tmp_plist1[0];
						double bend_cos = getCosineAngle(tmp_plist1[0], tmp_plist1[1], mypoint2f(0, 0), last_direct);
						bend_angle = acos(bend_cos) * 180 / 3.1415926;
						last_direct = tmp_plist1.back() - tmp_plist1[tmp_plist1.size() - 2];
						if (bend_angle < 15)
						{
							connect_plist.pop_back();
							connect_plist.insert(connect_plist.end(), tmp_plist1.begin(), tmp_plist1.end());
							pvec_i++;
						}
					}
				}
				if (!connect_plist.empty())
				{
					double peri = getPerimeterofVec(&connect_plist);
					edge_peri.push_back(peri);
					edge_plist.push_back(connect_plist);
				}
			}
			int sim_sum = 0;
			int max_index = 0; int max_c = 0;  //最大plist的edge 的index
			for (unsigned int i = 0; i < edge_plist.size(); ++i)
			{
				//vector<mypoint2f> tmp_plist=edge_plist[i];
				double peri = edge_peri[i];
				int tmp = fixednum * peri / perimeter;//每条边化简后的点数(大概的)（每边包括两个端点）
				if (tmp < 2)
				{
					tmp = 2;
				}
				int ensure_ecount = tmp + 15;  //!!!!!!!!!!!!!
				if (ensure_ecount > edge_plist[i].size())//化简后的point >  原本的point个数
				{
					int csize = edge_plist[i].size();
					int accumulate = csize - 1;
					int j = 1;
					int addnum = 0;
					while ((addnum + csize - 1) < ensure_ecount)//不把末尾点算上
					{
						addnum = accumulate *j;
						j++;
					}
					//每两个点之间插入i-1个点(平均分成i份，即除以i)，共addnum个
					vector<mypoint2f>::iterator iter = edge_plist[i].begin() + 1;
					while (iter != edge_plist[i].end())
					{
						mypoint2f p1 = *(iter - 1);
						mypoint2f p2 = *iter;
						double x_seg = (p2.x - p1.x) / (double)j;
						double y_seg = (p2.y - p1.y) / (double)j;
						for (int k = 1; k < j; ++k)
						{
							mypoint2f midp(p2.x - x_seg*k, p2.y - y_seg*k);
							iter = edge_plist[i].insert(iter, midp);
						}
						iter = iter + j;
					}
				}
				if (tmp > max_c)
				{
					max_c = tmp;
					max_index = i;
				}
				edge_count.push_back(tmp);	
				sim_sum = sim_sum + tmp;
			}
			/*原来:通过p_vec ，一段段 ，对edge_peri ，edge_plist  ，edge_count赋值*/
			//pvecsize = p_vec->size();
			//int sim_sum = 0;
			//int max_index = 0; int max_c = 0;  //最大plist的edge 的index
			//for (int i = 1; i < pvecsize; ++i)
			//{
			//	vector<int> anedge;
			//	vector<mypoint2f> tmp_plist;
			//	anedge.push_back(p_vec->at(i - 1));
			//	anedge.push_back(p_vec->at(i));
			//	jointMethod(&anedge, &tmp_plist, 0);
			//	double peri = getPerimeterofVec(&tmp_plist);
			//	edge_peri.push_back(peri);
			//	int tmp = fixednum * peri / perimeter;//每条边化简后的点数(大概的)（每边包括两个端点）
			//	if (tmp < 2)
			//	{
			//		tmp = 2;
			//	}
			//	int ensure_ecount = tmp + 10;  //?????????????
			//	if (ensure_ecount > tmp_plist.size())//化简后的point >  原本的point个数
			//	{
			//		int csize = tmp_plist.size();
			//		int accumulate = csize - 1;
			//		int j = 1;
			//		int addnum = 0;
			//		while ((addnum + csize - 1) < ensure_ecount)//不把末尾点算上
			//		{
			//			addnum = accumulate *j;
			//			j++;
			//		}
			//		//每两个点之间插入i-1个点(平均分成i份，即除以i)，共addnum个
			//		vector<mypoint2f>::iterator iter = tmp_plist.begin() + 1;
			//		while (iter != tmp_plist.end())
			//		{
			//			mypoint2f p1 = *(iter - 1);
			//			mypoint2f p2 = *iter;
			//			double x_seg = (p2.x - p1.x) / (double)j;
			//			double y_seg = (p2.y - p1.y) / (double)j;
			//			for (int k = 1; k < j; ++k)
			//			{
			//				mypoint2f midp(p2.x - x_seg*k, p2.y - y_seg*k);
			//				iter = tmp_plist.insert(iter, midp);
			//			}
			//			iter = iter + j;
			//		}
			//	}
			//	if (tmp > max_c)
			//	{
			//		max_c = tmp;
			//		max_index = i - 1;
			//	}
			//	edge_plist.push_back(tmp_plist);
			//	edge_count.push_back(tmp);	
			//	sim_sum = sim_sum + tmp;
			//}
			/*把少添加的点平均分到edge中*/
			int cha = fixednum - sim_sum + (int)edge_plist.size();
			if (cha != 0)  
			{
				int ttt = cha / (int)edge_count.size();
				if (cha > 0)
				{
					if (ttt <1)
						ttt = 1;
					int sum_ttt = 0;
					int k = 0;
					int k_count = edge_count.size();
					while (sum_ttt < cha && k < k_count)
					{
						if (edge_count[k] != 2)
						{
							if ((edge_count[k] + ttt) <= edge_plist[k].size())
							{
								edge_count[k] = edge_count[k] + ttt;
								sum_ttt = sum_ttt + abs(ttt);
							}
						}
						k++;
					}
					if (sum_ttt < cha)
						edge_count[max_index] = edge_count[max_index] + (cha - sum_ttt);
				}
				else
				{
					if (ttt >-1)
						ttt = -1;
					int sum_ttt = 0;
					int k = 0;
					int k_count = edge_count.size();
					while (sum_ttt > cha && k < k_count)
					{
						if (edge_count[k] != 2)
						{
							if ((edge_count[k] + ttt) >2)
							{
								edge_count[k] = edge_count[k] + ttt;
								sum_ttt = sum_ttt + ttt;
							}
						}
						k++;
					}
					if (sum_ttt > cha)
						edge_count[max_index] = edge_count[max_index] + (cha - sum_ttt);
				}		
				//edge_count[max_index] = edge_count[max_index] + cha;
			}
		}
		//edge_count[]存储了每条边应简化的个数，包括边的起始端点。现在从edge_plist中等距采样
		for (unsigned int i = 0; i < edge_count.size(); ++i)
		{
			int object_edge = edge_count[i];
			vector<mypoint2f> edgeplist = edge_plist[i];
			if (edgeplist.size() < object_edge)
				cout << "error" << endl;
			if (edgeplist.size() == object_edge)//若目标点数>实际的点数，直接存入adjust_plist
			{
				for (unsigned int j = 0; j < edgeplist.size()-1; ++j)
				{
					adjust_plist->push_back(edgeplist[j]);
				}
			}
			else if (object_edge==2)
			{
				adjust_plist->push_back(edgeplist.front());
			}
			else
			{
				double interval = edge_peri[i] / (double)(object_edge - 1);
				double accu_len = 0;
				double accu_interval = interval;
				adjust_plist->push_back(edgeplist.front());
				int add_number = 1;
				for (unsigned int j = 1; j < edgeplist.size(); ++j)
				{
					double tmp_len = LineSegmentLength(edgeplist[j], edgeplist[j - 1]);
					accu_len = accu_len + tmp_len;
					while (accu_len >= accu_interval)
					{
						double cutoff = 1 - (accu_len - accu_interval) / tmp_len;  //比例
						mypoint2f dic_vec = edgeplist[j] - edgeplist[j - 1];
						dic_vec.x = dic_vec.x*cutoff;
						dic_vec.y = dic_vec.y*cutoff;
						mypoint2f insert_p = dic_vec + edgeplist[j - 1];
						adjust_plist->push_back(insert_p);
						add_number++;
						accu_interval = accu_interval + interval;  //!!!
					}
				}
				if (adjust_plist->back() == edgeplist.back() || add_number == object_edge)
					adjust_plist->pop_back();
			}
		}//for

	}
}
//数据集的增加，，，，
void MyGraph::dataIncreasing(vector<mypoint2f> *plist, int fixednum)
{
	//plist  front() != back()，以下注释掉的 是原来的代码。原来的代码必须保证plist是一个front!=back的patch，
	/*if (plist->front() == plist->back())
		plist->pop_back();
	mypoint2f frontp = plist->front();
	plist->push_back(frontp);*/
	//现在，plist可以是一个front=back的patch ，也可以是一个line
	int csize = plist->size();
	int accumulate = csize - 1;
	int i = 1;
	int addnum = 0;
	while ((addnum + csize-1) < fixednum)//不把加进去的点算上
	{
		addnum = accumulate *i;
		i++;
	}
	//每两个点之间插入i-1个点(平均分成i份，即除以i)，共addnum个
	vector<mypoint2f>::iterator iter = plist->begin()+1;
	while (iter != plist->end())
	{
		mypoint2f p1 = *(iter-1);
		mypoint2f p2 = *iter;
		double x_seg = (p2.x - p1.x) / (double)i;
		double y_seg = (p2.y - p1.y) / (double)i;
		for (int k = 1; k < i; ++k)
		{
			mypoint2f midp(p2.x - x_seg*k, p2.y - y_seg*k);
			iter = plist->insert(iter, midp);
		}
		iter = iter + i;
	}
	if (plist->size() != (addnum + csize))
		cout << "error" << endl;
	//plist->pop_back();//最后再pop()掉，还是front()!=back()

	//while (plist->size() <= fixednum)
	//{
	//	vector<mypoint2f>::iterator iter_p = plist->begin()+1;
	//	while (iter_p != plist->end())
	//	{
	//		mypoint2f tmp_p((iter_p->x + (iter_p - 1)->x) / 2, (iter_p->y + (iter_p - 1)->y) / 2);
	//		iter_p = plist->insert(iter_p, tmp_p);
	//		iter_p = iter_p + 2;
	//	}
	//}
	//plist->pop_back();//最后再pop()掉，还是front()!=back()
}
//方法一:计算bend angle function，没有转换成fourier的形式，adjusted_plist的front()=back()
void MyGraph::getbendAngleFunction(vector<mypoint2f> *adjusted_plist,vector<vector<double>>* bendAf,int sample_num)
{
	//先在简化的contour的基础上计算bendAngleFunction，再等距采样 ,采样点数为sample_num
	vector<vector<double>> tmp_bd;
	double cumulate_len = 0;
	//double cumulate_angle = 0;
	for (unsigned int j = 1; j < adjusted_plist->size() - 1; ++j)//front=back
	{
		cumulate_len = cumulate_len + LineSegmentLength(adjusted_plist->at(j-1), adjusted_plist->at(j));
		double costheta = getCosineAngle(adjusted_plist->at(j - 1), adjusted_plist->at(j), adjusted_plist->at(j), adjusted_plist->at(j + 1));
		double cross_tmp = cross((adjusted_plist->at(j) - adjusted_plist->at(j - 1)), (adjusted_plist->at(j + 1) - adjusted_plist->at(j)));//cross(a,b) aXb<0, b在a的右边（顺时针方向） 
		double agl = costheta;   //acos(costheta) * 180 / 3.1415926;//转化成角度 //
		if (cross_tmp<0 && agl!=0)
			agl = agl*(-1);
		//cumulate_angle = cumulate_angle + agl;
		vector<double> baf;
		baf.push_back(cumulate_len);
		baf.push_back(agl);
		tmp_bd.push_back(baf);
	}
	cumulate_len = cumulate_len + LineSegmentLength(adjusted_plist->at(adjusted_plist->size() - 2), adjusted_plist->back());//back()=at(0)
	double costheta = getCosineAngle(adjusted_plist->at(adjusted_plist->size() - 2), adjusted_plist->at(0), adjusted_plist->at(0), adjusted_plist->at(1));
	double cross_tmp = cross((adjusted_plist->at(0) - adjusted_plist->at(adjusted_plist->size() - 2)), (adjusted_plist->at(1) - adjusted_plist->at(0)));//cross(a,b) aXb<0, b在a的右边（顺时针方向） 
	double agl = costheta; // acos(costheta) * 180 / 3.1415926;  //转化成角度 //
	if (cross_tmp<0 && agl != 0)
		agl = agl*(-1);
	//cumulate_angle = cumulate_angle + agl;
	vector<double> baf;
	baf.push_back(cumulate_len);
	baf.push_back(agl);
	tmp_bd.push_back(baf);
	for (unsigned int j = 0; j < tmp_bd.size(); ++j)//单位化
	{
		tmp_bd.at(j)[0] = tmp_bd.at(j)[0] / cumulate_len;
		//bendAf->at(j)[0] = bendAf->at(j)[0] / cumulate_len;
	}
	//等距采样
	double inter_len = 1 / double(sample_num);
	double cumu_interval = inter_len;
	int index_origi = 0;
	for (int i = 0; i < sample_num-1; ++i)
	{
		if(cumu_interval > tmp_bd[index_origi][0])
		{
			index_origi++;
		}
		vector<double> dis_angle = tmp_bd[index_origi];
		dis_angle[0] = cumu_interval;
		bendAf->push_back(dis_angle);
		cumu_interval = cumu_interval + inter_len;
	}
	bendAf->push_back(tmp_bd.back());
}
//方法二：计算傅里叶描述子，输入contour point，输出一个一维的向量FD， 用了radial destance(shape signature), evol_plist的front()=back()
void MyGraph::getFourierDescriptor(vector<mypoint2f>* contour, vector<complex<double>>* FD, int pch)//最一开始计算两种FD
{
	////计算cumulative bendAngleFunction，， contour front()=back()
	//vector<double> tmp_bd;
	//vector<double> cumu_len;
	//double cumulate_len = 0;
	//double cumulate_angle = 0;
	//int c_size = contour->size() - 1;
	//for (int j = 1; j < c_size; ++j)//front=back
	//{
	//	cumulate_len = cumulate_len + LineSegmentLength(contour->at(j - 1), contour->at(j));
	//	double costheta = getCosineAngle(contour->at(j - 1), contour->at(j), contour->at(j), contour->at(j + 1));
	//	if (costheta>1)
	//		costheta=1;
	//	double cross_tmp = cross((contour->at(j) - contour->at(j - 1)), (contour->at(j + 1) - contour->at(j)));//cross(a,b) aXb<0, b在a的右边（顺时针方向） 
	//	double agl = acos(costheta) * 180 / 3.1415926;
	//	if (cross_tmp<0 && agl!=0)
	//		agl = agl*(-1);
	//	cumulate_angle = cumulate_angle + agl;
	//	cumu_len.push_back(cumulate_len);
	//	tmp_bd.push_back(cumulate_angle);//  
	//}
	//cumulate_len = cumulate_len + LineSegmentLength(contour->at(contour->size() - 2), contour->back());//back()=at(0)
	//double costheta = getCosineAngle(contour->at(contour->size() - 2), contour->at(0), contour->at(0), contour->at(1));
	//if (costheta>1)
	//	costheta = 1;
	//double cross_tmp = cross((contour->at(0) - contour->at(contour->size() - 2)), (contour->at(1) - contour->at(0)));//cross(a,b) aXb<0, b在a的右边（顺时针方向） 
	//double agl = acos(costheta) * 180 / 3.1415926;
	//if (cross_tmp<0 && agl!=0)
	//	agl = agl*(-1);
	//cumulate_angle = cumulate_angle + agl;
	//cumu_len.push_back(cumulate_len);
	//tmp_bd.push_back(cumulate_angle);//  
	//int plistsize = tmp_bd.size();
	//complex<double> ii(0, 1);//复数 j
	//for (int i = 0; i < plistsize; ++i)
	//{
	//	complex<double> sum_cmp = 0;
	//	for (int j = 0; j < plistsize; ++j)
	//	{
	//		double tmp = (-1) * 2 * 3.1415926*i*j / plistsize;
	//		complex<double> cmp(tmp, 0);
	//		complex<double> cmp_2 = exp(ii*cmp);
	//		complex<double> cmp_rt(tmp_bd[j], 0);
	//		complex<double> cmp_3 = cmp_rt*cmp_2;
	//		sum_cmp = sum_cmp + cmp_3;
	//	}
	//	complex<double> cmpsize(plistsize, 0);
	//	sum_cmp = sum_cmp / cmpsize;
	//	FD->push_back(sum_cmp);
	//}

	//2. contour front()=back() radial function 
	mypoint2f ct_p = plrGraph.patches[pch].centre_position;
	int plistsize = contour->size();
	//计算radial destance，存入radial_vec
	vector<double> radial_vec;
	for (int i = 0; i < plistsize; ++i)
	{
		double tmp_l = LineSegmentLength(contour->at(i), ct_p);// contour->front() != back()
		radial_vec.push_back(tmp_l);
	}
	//Fourier transform 求出a0-an
	complex<double> ii(0, 1);//复数 j
	for (int i = 0; i < plistsize; ++i)
	{
		complex<double> sum_cmp = 0;
		for (int j = 0; j < plistsize; ++j)
		{		
			double tmp = (-1) * 2 * 3.1415926*i*j/plistsize;
			complex<double> cmp(tmp,0);
			complex<double> cmp_2 = exp(ii*cmp);
			complex<double> cmp_rt(radial_vec[j], 0);
			complex<double> cmp_3 = cmp_rt*cmp_2;
			sum_cmp = sum_cmp + cmp_3;
		}
		complex<double> cmpsize(plistsize, 0);
		sum_cmp = sum_cmp / cmpsize;
		FD->push_back(sum_cmp);
	}
}


void MyGraph::localProcess_two(int simply_num)
{
	//1、圆度，公式4*pi*area/permeter^2. 如果是圆，改特征值=1
	//2、长度(elongation)、(离心率eccentricity)、length of major axis/length of minor axis	
	cout << "localProcess_two()" << endl;
	vector<pair<int, vector<double>>> patch_tmpweight;
	vector<double> cir_vec;
	vector<double> long_vec;
	vector<double> area_vec;
	int pchsize = plrGraph.patches.size();
	for (int i = 0; i < pchsize; ++i)//付初值 vector<pair<int,vector<double>>> patchfeatures;//int 对应patchindex，vector<double>对应特征向量
	{
		vector<double> tmp;
		patch_tmpweight.push_back(make_pair(i, tmp));
	}
	for (int i = 0; i < pchsize; ++i)
	{
		if (plrGraph.patches[i].type >= 0)
		{
			vector<mypoint2f> total_plist;
			jointMethod(&plrGraph.patches[i].pointIndexOfPatch, &total_plist, 0); //total_plist首 = 尾

			vector<double> fvec;//特征向量
			getThreeFeature(&fvec, &total_plist, plrGraph.patches[i].area, plrGraph.patches[i].perimeter);

			patch_tmpweight[i].first = i;
			patch_tmpweight[i].second = fvec;
		}
	}
	vector<pair<int, double>> tt_w;//patch的权值=三个ft的和
	for (int i = 0; i < pchsize; ++i)
	{
		tt_w.push_back(make_pair(i, 0));
		if (plrGraph.patches[i].type >= 0)
		{
			double feature_sum = 0;
			for (unsigned int j = 0; j < patch_tmpweight[i].second.size(); ++j)
			{
				feature_sum = feature_sum + patch_tmpweight[i].second[j];
			}
			tt_w[i].second = feature_sum;
			//vector<int> adjpch = getNeighborPatches(i);
			//if (!adjpch.empty())
			//{
			//	//计算局部总
			//	double local_cir = 0;
			//	double local_long = 0;
			//	double local_area = 0;
			//	for (unsigned int j = 0; j < adjpch.size(); ++j)
			//	{
			//		local_cir = local_cir + patchfeatures[adjpch[j]].second[0];
			//		local_long = local_long + patchfeatures[adjpch[j]].second[1];
			//		local_area = local_area + patchfeatures[adjpch[j]].second[2];
			//	}
			//	//邻域patch不包括自己，把本patch的特征加进去
			//	local_cir = local_cir + patchfeatures[i].second[0];
			//	local_long = local_long + patchfeatures[i].second[1];
			//	local_area = local_area + patchfeatures[i].second[2];
			//	//计算各权值(圆度+长度+面积之比)		
			//	double weight = patchfeatures[i].second[0] / local_cir + patchfeatures[i].second[1] / local_long + patchfeatures[i].second[2] / local_area;
			//	tt_w[i].second = weight;//tt_w.push_back(make_pair(i, weight));
			//}
			//else
			//{
			//	//单个patch，没有neighbor,权值=各特征相加
			//	tt_w[i].second = patchfeatures[i].second[0] + patchfeatures[i].second[1] + patchfeatures[i].second[2];
			//}
		}
	}
	sort(tt_w.begin(), tt_w.end(), mycompare);//从小到大排序

	/*权值平均分成5份，用不同颜色表示，查看分布*/
	//int segment = tt_w.size() / 5;
	//int i_seg = 0;
	//int step = segment;
	//for (int j = 0; j < 5; ++j)
	//{
	//	for (i_seg; i_seg < step; ++i_seg)// tt_w.size() 
	//	{
	//		int pindex = tt_w[i_seg].first;
	//		if (plrGraph.patches[pindex].type >= 0)
	//		{
	//			plrGraph.patches[pindex].type = j+1;
	//		}		
	//	}
	//	if (j == 3)
	//	{
	//		step = tt_w.size();
	//	}
	//	else
	//	{
	//		step = step + segment;
	//	}
	//}

	vector<vector<double>> similar_dis;
	for (int i = 0; i < tt_w.size(); ++i)// tt_w.size() 
	{
		int pindex = tt_w[i].first;
		if (plrGraph.patches[pindex].type >= 0)
		{
			vector<int> adjpch = getNeighborPatches(pindex, 0);//pindex 与 adjpch[]分别合并，并计算合并前后的相似度大小，选择相似度最大的进行合并
			if (!adjpch.empty())
			{
				//先计算pindex与邻域patch的pointlist大小，若有非常小的（如<10个点）则设置化简个数为10，否则是128
				adjpch.push_back(pindex);//把中心patch加进去
				vector<vector<mypoint2f>> nei_plist;
				int min_count = 65535; //simply_num
				for (unsigned int j = 0; j < adjpch.size(); ++j)
				{
					vector<mypoint2f> contour_plist;
					jointMethod(&plrGraph.patches[adjpch[j]].pointIndexOfPatch, &contour_plist, 0);
					contour_plist.pop_back();
				
					if (contour_plist.size() < min_count)
					{
						min_count = contour_plist.size();
					}
					nei_plist.push_back(contour_plist); //front()!=back()
				}
				if (min_count >0 && min_count<128)
				{
					simply_num = min_count;
				}
				//计算adjpch中的(neighbor patch+pindex)的 fourier descriptor
				vector<vector<double>> four_des;
				for (unsigned int j = 0; j < adjpch.size(); ++j)//adjpch最后一个是中心patch: pindex
				{
					vector<mypoint2f> evolu_plist;
					dataReduction_byheap(&nei_plist[j], simply_num, plrGraph.patches[adjpch[j]].perimeter, &evolu_plist);
					//dataReduction(&nei_plist[j], &evolu_plist, plrGraph.patches[adjpch[j]].perimeter, simply_num);//front()!=back()
					mypoint2f frontp = evolu_plist.front();//添加一个头，使evol_plist 首=尾
					evolu_plist.push_back(frontp);
					vector<double> pch_FD;
					getFourierDescriptor(&evolu_plist, &pch_FD,1);
					four_des.push_back(pch_FD);
				}
				//分别合并patch,并计算新合并的FD,之后比较pindex & combinepatch , neighpatch & combinepatch,选其中较小的similarity值，记录下neighborpatch 的index
				double min_similarity = 65535;
				int combine_pindex = -1;
				for (unsigned int j = 0; j < adjpch.size() - 1; ++j)//adjpch最后一个是中心patch: pindex
				{
					vector<int> combine_pvec;
					combineTwoPatch_getpvec(pindex, adjpch[j], &combine_pvec);//combine_plist  front()！= back()
					if (!combine_pvec.empty())
					{
						vector<mypoint2f> combine_plist;
						jointMethod(&combine_pvec, &combine_plist, 0);
						combine_plist.pop_back();
						vector<mypoint2f> evolu_plist;
						dataReduction(&combine_plist, &combine_pvec, -1, simply_num, &evolu_plist); //front()!=back()
						mypoint2f frontp = evolu_plist.front();//添加一个头，使evol_plist的front()=back()
						evolu_plist.push_back(frontp);
						vector<double> pch_FD;
						getFourierDescriptor(&evolu_plist, &pch_FD,1);
						//double dis1 = similarityMeasure(&four_des[j], &pch_FD);//和第j个neighborpatch融合的，与four_des中的第j个patch的FD比较
						//double dis2 = similarityMeasure(&four_des.back(), &pch_FD);//与pindex的FD比较
						//if (dis1 < min_similarity)
						//{
						//	min_similarity = dis1;
						//	combine_pindex = adjpch[j];
						//}
						//else if (dis2 < min_similarity)
						//{
						//	min_similarity = dis2;
						//	combine_pindex = adjpch[j];
						//}
						double dis1 =similarityMeasure(&four_des.back(), &pch_FD);//pindex & combinepatch
						if (dis1 < min_similarity)
						{
							min_similarity = dis1;
							combine_pindex = adjpch[j];
						}
					}
				}
				if (combine_pindex != -1)
				{
					combineTwoPatch(pindex, combine_pindex);
				}
			}
		}
	}
}

//Method 1  结果在文件夹merge
double MyGraph::localProcess_heap(int final_pnum, double similarity_t, int simply_num, int iteration)
{
	//用heap， 每次删除一个patch，例如把A与B合并，B的weight需要重新计算，并且更新他在堆中的位置
	//1、圆度，公式4*pi*area/permeter^2. 如果是圆，改特征值=1
	//2、长度(elongation)、(离心率eccentricity)、length of major axis/length of minor axis	
	cout << "localProcess_two()" << endl;
	vector<pair<int, vector<double>>> patch_tmpweight;
	vector<double> cir_vec;
	//vector<double> long_vec;
	vector<double> area_vec;

	ofstream outfile_width("E:/svgfile/m2/patch_width.txt");//"E:/record_similarity.txt"
	outfile_width.is_open();
	vector<double> aver_w;
	vector<double> long_vec;

	int pchsize = plrGraph.patches.size();
	int pch_count = pchsize;//patch 从heap中删除时减1，最后通过比较pch_count < stopNum来结束简化
	int stop_num = pchsize *((double)final_pnum / 100.0);

	for (int i = 0; i < pchsize; ++i)
	{
		vector<double> tmp;
		patch_tmpweight.push_back(make_pair(i, tmp)); //付初值 vector<pair<int,vector<double>>> patchfeatures;//int 对应patchindex，vector<double>对应特征向量
		if (plrGraph.patches[i].type >= 0)
		{
			vector<mypoint2f> total_plist;
			jointMethod(&plrGraph.patches[i].pointIndexOfPatch, &total_plist, 0);
			total_plist.pop_back();
			//1、圆度
			double circularity = 4 * 3.1415926*plrGraph.patches[i].area / pow(plrGraph.patches[i].perimeter, 2);
			//2、长度，方法二，计算出最大矩形bounding box
			double elongation = 0;
			vector<mypoint2f> plist_box;
			mypoint2f rec_center(0, 0);
			getBoundRectangle(&total_plist, rec_center, &plist_box);
			double box_w = fabs(plist_box[0].x - plist_box[1].x);
			double box_h = fabs(plist_box[1].y - plist_box[2].y);
			if (box_w < box_h)
			{
				elongation = box_w / box_h;
				outfile_width << box_w << ' ';
				aver_w.push_back(box_w);
			}
			else
			{
				elongation = box_h / box_w;
				outfile_width << box_h << ' ';
				aver_w.push_back(box_h);
			}
			long_vec.push_back(elongation);
			//3、面积比，object_area/min-bounding rectangle_area,  convex or concave?
			double are_rect = box_w*box_h;
			double are_ratio = plrGraph.patches[i].area / are_rect;
			vector<double> fvec;//特征向量
			fvec.push_back(circularity);//圆度
			fvec.push_back(elongation);//长度
			fvec.push_back(are_ratio);//面积

			cir_vec.push_back(fvec[0]);
			long_vec.push_back(fvec[1]);
			area_vec.push_back(fvec[2]);

			patch_tmpweight[i].second = fvec;
		}
	}
	outfile_width.close();
	double retn_num = 0;
	retn_num = GetAverange(aver_w);
	double aver_long = GetAverange(long_vec);//patch的长宽比小于这个数，才可以拉伸成一条线，，

	vector<double> tt_w;//patch的权值（三个ft的和）
	for (int i = 0; i < pchsize; ++i)
	{
		tt_w.push_back(-1);
		if (plrGraph.patches[i].type >= 0)
		{
			double feature_sum = 0;
			for (unsigned int j = 0; j < patch_tmpweight[i].second.size(); ++j)
			{
				feature_sum = feature_sum + patch_tmpweight[i].second[j];
			}
			tt_w[i]= feature_sum;
		}
	}
	//sort(tt_w.begin(), tt_w.end(), mycompare);//从小到大排序 ,被删除的patch的weight=-1
	bool stopflag = false;
	while (pch_count >= stop_num && stopflag == false)  //stop_num=pchsize *((double)final_pnum / 100.0);
	{
		//建堆，tt_w保存了从第0个点到底n个点weight值，堆heap_vec保存每个patch的权值在tt_w中的下标，也是patch 的index
		vector<int> heap_vec;
		vector<int> ele_index;
		for (unsigned int i = 0; i < tt_w.size(); ++i)
		{
			ele_index.push_back(-1);
			if (plrGraph.patches[i].type >= 0)
			{
				heap_vec.push_back(i);
			}
			//heap_vec.push_back(i);
			//ele_index.push_back(0);
		}
		//建堆
		int h_size = heap_vec.size();
		for (int i = h_size / 2; i >= 1; i--)
		{
			heap_siftdown(i, &heap_vec, &tt_w, &ele_index);
		}
		//建立索引，知道某个权值的下标，可以直接找到在heap的位置，位置是从1开始，注意更新
		for (int i = 0; i < h_size; ++i)
		{
			ele_index[heap_vec[i]] = i;// ele_index
		}
		//记录每个patch的化简后的pointlist ，在每次合并之后记得更新！！
		vector<vector<mypoint2f>> patch_plistnum;
		for (unsigned int i = 0; i < plrGraph.patches.size(); ++i)
		{
			patch_plistnum.push_back(vector<mypoint2f>());
			if (plrGraph.patches[i].type >= 0)
			{
				vector<mypoint2f> contour_plist;
				vector<mypoint2f> simplified_plist;
				jointMethod(&plrGraph.patches[i].pointIndexOfPatch, &contour_plist, 0); //front()=back()
				contour_plist.pop_back();
				
				if (contour_plist.size() < simply_num)
				{
					contour_plist.push_back(contour_plist.front());
					dataIncreasing(&contour_plist, simply_num);
					contour_plist.pop_back();
					dataReduction_byheap(&contour_plist, simply_num, plrGraph.patches[i].perimeter, &simplified_plist);
					//dataReduction(&contour_plist, &simplified_plist, plrGraph.patches[i].perimeter, simply_num);
				}
				else
				{
					dataReduction_byheap(&contour_plist, simply_num, plrGraph.patches[i].perimeter, &simplified_plist);
				}
				patch_plistnum[i] = simplified_plist;//front() != back()
			}
		}
		string file_str = "E:/mini_similarity";
		file_str.append(to_string(iteration));
		file_str.append(".txt");
		ofstream outfile1(file_str);//"E:/record_similarity.txt"
		outfile1.is_open();

		vector<double> record_similarity;//输出相似度
		int reduction_num = simply_num;//规定计算相似度的点的个数（即化间个数），若数据点过小，该值可以减少一点
		//	int stop_num = final_pnum;
		stopflag = true;//若以下while循环没有执行合并操作，即patch没有减少，则stopflag=true，整体结束
		while (h_size >0 && pch_count >= stop_num)//while (!heap_vec.empty())//
		{
			int pindex = heap_vec[0];
			//if (plrGraph.patches[pindex].type >= 0) //if (tt_w[pindex] != 0)
			vector<int> adjpch = getNeighborPatches(pindex, 0);//pindex 与 adjpch[]分别合并，并计算合并前后的相似度大小，选择相似度最大的进行合并
			if (!adjpch.empty())
			{
				reduction_num = simply_num;

				/*计算pindex的 fourier descriptor*/
				//vector<vector<double>> four_des;
				//for (unsigned int j = 0; j < adjpch.size(); ++j)//adjpch最后一个是中心patch: pindex
				//{
				//	vector<mypoint2f> evolu_plist;
				//	if (nei_plist[j].size() < reduction_num)
				//	{
				//		dataIncreasing(&nei_plist[j], reduction_num);
				//		dataReduction(&nei_plist[j], &evolu_plist, plrGraph.patches[adjpch[j]].perimeter, reduction_num);
				//	}
				//	else
				//	{
				//		dataReduction(&nei_plist[j], &evolu_plist, plrGraph.patches[adjpch[j]].perimeter, reduction_num);
				//	}
				//	//dataReduction(&nei_plist[j], &evolu_plist, plrGraph.patches[adjpch[j]].perimeter, reduction_num);//front()!=back()
				//	mypoint2f frontp = evolu_plist.front();//添加一个头，使evol_plist 首=尾
				//	evolu_plist.push_back(frontp);
				//	vector<double> pch_FD;
				//	getFourierDescriptor(&evolu_plist, &pch_FD);
				//	four_des.push_back(pch_FD);
				//}
				vector<double> center_FD;
				vector<mypoint2f> center_evoplist = patch_plistnum[pindex];
				mypoint2f ftp = center_evoplist.front();//添加一个头，使evol_plist 首=尾
				center_evoplist.push_back(ftp);
				getFourierDescriptor(&center_evoplist, &center_FD, 0);//计算傅里叶descriptor
				//分别合并patch,并计算新合并的FD,之后比较pindex & combinepatch , neighpatch & combinepatch,选其中较小的similarity值，记录下neighborpatch 的index
				double min_similarity = 65535;
				int combine_pindex = -1;
				vector<mypoint2f> record_complist;
				for (unsigned int j = 0; j < adjpch.size(); ++j)//adjpch最后一个是中心patch: pindex
				{
					vector<int> combine_pvec;
					combineTwoPatch_getpvec(pindex, adjpch[j], &combine_pvec);  //新patch的pointlist，front()！= back()，间隔1个取一个点
					if (!combine_pvec.empty())
					{
						vector<mypoint2f> combine_plist;//没化简之前的pointlist
						jointMethod(&combine_pvec, &combine_plist, 0); //front()= back()
						double peri_comb = getPerimeterofVec(&combine_plist);//周长
						combine_plist.pop_back();  //front()！= back()
						vector<double> pch_FD;  //新patch的FD
						//getFDFeature(&combine_plist,&combine_pvec,-1, reduction_num, 0, &pch_FD);
						vector<mypoint2f> evolu_plist;
						if (combine_plist.size() < reduction_num)
						{
							combine_plist.push_back(combine_plist.front());
							dataIncreasing(&combine_plist, reduction_num);
							combine_plist.pop_back();
							dataReduction(&combine_plist, &combine_pvec, -1, reduction_num, &evolu_plist); //front() != back()化简得evolu_plist
						}
						else
						{
							dataReduction(&combine_plist, &combine_pvec, -1, reduction_num, &evolu_plist);//front()!=back()化简得evolu_plist
						}
						mypoint2f frontp = evolu_plist.front();//添加一个头，使evol_plist的front()=back()
						evolu_plist.push_back(frontp);
						getFourierDescriptor(&evolu_plist, &pch_FD, 0);//计算傅里叶descriptor		
						double dis1 = similarityMeasure(&center_FD, &pch_FD);//similarityMeasure(&four_des.back(), &pch_FD);//pindex & combinepatch
						//outfile1 << dis1;
						//outfile1 << " ";

						if (dis1 < min_similarity)
						{
							min_similarity = dis1;
							combine_pindex = adjpch[j];
							evolu_plist.pop_back();
							record_complist = evolu_plist;// combine_plist;//记录下来这个合并的且化简的patch的plist,
						}
						/*
						vector<mypoint2f> plist_adj = patch_plistnum[adjpch[j]];
						frontp = plist_adj.front();//使front()=back()
						plist_adj.push_back(frontp);
						vector<double> FD_adjpch;
						getFourierDescriptor(&plist_adj, &FD_adjpch,0);//计算傅里叶descriptor
						double dis2 = similarityMeasure(&FD_adjpch, &pch_FD);
						if (dis2 < min_similarity )
						{
						if (dis1 > min_similarity)
						{
						min_similarity = dis2;
						combine_pindex = adjpch[j];
						evolu_plist.pop_back();
						record_complist = evolu_plist;
						}
						else
						min_similarity = dis2;
						}*/
					}
				}
				//outfile1 << "\n";
				outfile1 << min_similarity;
				outfile1 << " ";
				if (combine_pindex != -1 && min_similarity <= similarity_t)
				{
					stopflag = false;//只要有patch删除，就不会停止，若while循环不执行这句话说明没有可以合并的patch
					//attachLineProcess(int pchindex)
					record_similarity.push_back(min_similarity);
					//combineTwoPatch(pindex, combine_pindex);//留pindex，或删除combine_pindex
					combineTwoPatch(combine_pindex, pindex);//留combine_pindex ，删除pindex //从堆顶删除的patch要和合并掉的patch（typr变为-1）同时！！
					//更新新合并的patch的pointlist的数据点！！
					patch_plistnum[combine_pindex].swap(vector<mypoint2f>());
					patch_plistnum[combine_pindex] = record_complist; //front()!=back()
					//1.先删除堆顶的pindex
					int tmp_back = heap_vec.back();
					int tmp_mini = heap_vec[0];
					tt_w[tmp_mini] = -1;
					ele_index[tmp_mini] = -1;
					ele_index[tmp_back] = 0;//!!!
					heap_vec[0] = tmp_back;
					//heap_vec.pop_back();
					heap_siftdown(1, &heap_vec, &tt_w, &ele_index);//下沉调整
					heap_vec.pop_back();
					h_size--;
					pch_count--;
					//2.更新新合并的patch的weight,和在堆中的位置
					double feature_sum = 0;
					//	//1、圆度	
					//double circularity = 4 * 3.1415926*plrGraph.patches[combine_pindex].area / pow(plrGraph.patches[combine_pindex].perimeter, 2);
					//feature_sum = feature_sum + circularity;
					//	//2、长度，方法二，计算出最大矩形bounding box
					//double elongation = 0;
					//vector<mypoint2f> plist_box;
					//mypoint2f rec_center(0, 0);
					//getBoundRectangle(&record_complist, rec_center, &plist_box);//用化简之hou的plist计算bounding box
					//double box_w = fabs(plist_box[0].x - plist_box[1].x);
					//double box_h = fabs(plist_box[1].y - plist_box[2].y);
					//if (box_w < box_h)
					//	elongation = box_w / box_h;
					//else
					//	elongation = box_h / box_w;
					//feature_sum = feature_sum + elongation;
					//	//3、面积比，object_area/min-bounding rectangle_area,  convex or concave?
					//double are_rect = box_w*box_h;
					//double are_ratio = plrGraph.patches[combine_pindex].area / are_rect;
					//feature_sum = feature_sum + are_ratio;			
					vector<double> new_ft;
					getThreeFeature(&new_ft, &record_complist, plrGraph.patches[combine_pindex].area, plrGraph.patches[combine_pindex].perimeter);//用这个函数替代上面的代码
					for (unsigned int i_ft = 0; i_ft < new_ft.size(); ++i_ft)
					{
						feature_sum = feature_sum + new_ft[i_ft];
					}
					patch_tmpweight[combine_pindex].second = new_ft;
					//tt_w[combine_pindex] = feature_sum;
					if (ele_index[combine_pindex] == -1)//说明combine_pindex在之前 从堆顶删除 但是没有合并（变-1），再把他加进heap中
					{
						tt_w[combine_pindex] = feature_sum;
						heap_vec.push_back(combine_pindex);
						ele_index[combine_pindex] = heap_vec.size() - 1;
						heap_siftup(heap_vec.size(), &heap_vec, &tt_w, &ele_index);
					}
					else
						heap_update(ele_index[combine_pindex] + 1, feature_sum, &heap_vec, &tt_w, &ele_index);//注意+1
				}
				else
				{
					//similarity 大于阈值，不需要合并，从堆顶端删除	
					int tmp_back = heap_vec.back();
					int tmp_mini = heap_vec[0];
					tt_w[tmp_mini] = -1;
					ele_index[tmp_mini] = -1;
					ele_index[tmp_back] = 0;//!!!
					heap_vec[0] = tmp_back;
					heap_siftdown(1, &heap_vec, &tt_w, &ele_index);//下沉调整
					heap_vec.pop_back();
					h_size--;
				}
			}
			else
			{
				//没有neighbor,不需要合并。先用patchToCurve()把一个patch变成一条线，然后从堆顶端删除	
				if (patch_tmpweight[pindex].second[1] < aver_long)
				{
					vector<mypoint2f> nwcurve;
					patchToCurve(pindex, &nwcurve);//!!!
				}
				int tmp_back = heap_vec.back();
				int tmp_mini = heap_vec[0];
				tt_w[tmp_mini] = -1;
				ele_index[tmp_mini] = -1;
				ele_index[tmp_back] = 0;//!!!
				heap_vec[0] = tmp_back;
				//heap_vec.pop_back();
				heap_siftdown(1, &heap_vec, &tt_w, &ele_index);//下沉调整
				heap_vec.pop_back();
				h_size--;
			}
		}
		cout << "combines patch size=" << record_similarity.size() << endl;
		outfile1.close();
	}
	////建堆，tt_w保存了从第0个点到底n个点weight值，堆heap_vec保存每个patch的权值在tt_w中的下标，也是patch 的index
	//vector<int> heap_vec;
	//vector<int> ele_index;
	//for (unsigned int i = 0; i < tt_w.size(); ++i)
	//{
	//	ele_index.push_back(-1);
	//	if (plrGraph.patches[i].type >= 0)
	//	{
	//		heap_vec.push_back(i);		
	//	}
	//	//heap_vec.push_back(i);
	//	//ele_index.push_back(0);
	//}
	////建堆
	//int h_size = heap_vec.size();
	//for (int i = h_size / 2; i >= 1; i--)
	//{
	//	heap_siftdown(i, &heap_vec, &tt_w, &ele_index);
	//}
	////建立索引，知道某个权值的下标，可以直接找到在heap的位置，位置是从1开始，注意更新
	//for (int i = 0; i < h_size; ++i)
	//{
	//	ele_index[heap_vec[i]] = i;// ele_index
	//}
	////记录每个patch的化简后的pointlist ，在每次合并之后记得更新！！
	//vector<vector<mypoint2f>> patch_plistnum;
	//for (unsigned int i = 0; i < plrGraph.patches.size(); ++i)
	//{
	//	patch_plistnum.push_back(vector<mypoint2f>());
	//	if (plrGraph.patches[i].type >= 0)
	//	{
	//		vector<mypoint2f> contour_plist;
	//		vector<mypoint2f> simplified_plist;
	//		for (unsigned int k = 1; k < plrGraph.patches[i].pointIndexOfPatch.size(); ++k)
	//		{
	//			int verindex1 = plrGraph.patches[i].pointIndexOfPatch.at(k - 1);
	//			int verindex2 = plrGraph.patches[i].pointIndexOfPatch.at(k);
	//			mypoint2f point1 = sgraph.nodeList[verindex1].position;//=lineWithCommonPoints[verindex1].first
	//			mypoint2f point2 = sgraph.nodeList[verindex2].position;
	//			getPListOfaCurve(point1, point2, verindex1, verindex2, &contour_plist,0);
	//			contour_plist.pop_back();//删除一个重复的点
	//		}
	//		//patch_plistnum[i] = contour_plist;
	//		if (contour_plist.size() < simply_num)
	//		{
	//			dataIncreasing(&contour_plist, simply_num);
	//			dataReduction(&contour_plist, &simplified_plist, plrGraph.patches[i].perimeter, simply_num);
	//		}
	//		else
	//		{
	//			dataReduction(&contour_plist, &simplified_plist, plrGraph.patches[i].perimeter, simply_num);
	//		}
	//		patch_plistnum[i] = simplified_plist;//front() != back()
	//	}
	//}
	//string file_str = "E:/mini_similarity";
	//file_str.append(to_string(iteration));
	//file_str.append(".txt");
	//ofstream outfile1(file_str);//"E:/record_similarity.txt"
	//outfile1.is_open();
	//vector<double> record_similarity;//输出相似度
	//int reduction_num = simply_num;//规定计算相似度的点的个数（即化间个数），若数据点过小，该值可以减少一点
	////	int stop_num = final_pnum;
	//while (h_size >0 && pch_count >= stop_num)//while (!heap_vec.empty())//
	//{
	//	int pindex = heap_vec[0];
	//	//if (plrGraph.patches[pindex].type >= 0) //if (tt_w[pindex] != 0)
	//	vector<int> adjpch = getNeighborPatches(pindex, 0);//pindex 与 adjpch[]分别合并，并计算合并前后的相似度大小，选择相似度最大的进行合并
	//	if (!adjpch.empty())
	//	{
	//		adjpch.push_back(pindex);//把中心patch加进去
	//		reduction_num = simply_num;
	//			
	//		/*计算pindex的 fourier descriptor*/
	//		//vector<vector<double>> four_des;
	//			//for (unsigned int j = 0; j < adjpch.size(); ++j)//adjpch最后一个是中心patch: pindex
	//			//{
	//			//	vector<mypoint2f> evolu_plist;
	//			//	if (nei_plist[j].size() < reduction_num)
	//			//	{
	//			//		dataIncreasing(&nei_plist[j], reduction_num);
	//			//		dataReduction(&nei_plist[j], &evolu_plist, plrGraph.patches[adjpch[j]].perimeter, reduction_num);
	//			//	}
	//			//	else
	//			//	{
	//			//		dataReduction(&nei_plist[j], &evolu_plist, plrGraph.patches[adjpch[j]].perimeter, reduction_num);
	//			//	}
	//			//	//dataReduction(&nei_plist[j], &evolu_plist, plrGraph.patches[adjpch[j]].perimeter, reduction_num);//front()!=back()
	//			//	mypoint2f frontp = evolu_plist.front();//添加一个头，使evol_plist 首=尾
	//			//	evolu_plist.push_back(frontp);
	//			//	vector<double> pch_FD;
	//			//	getFourierDescriptor(&evolu_plist, &pch_FD);
	//			//	four_des.push_back(pch_FD);
	//			//}
	//		vector<double> center_FD;
	//		vector<mypoint2f> center_evoplist = patch_plistnum[pindex];
	//		mypoint2f ftp = center_evoplist.front();//添加一个头，使evol_plist 首=尾
	//		center_evoplist.push_back(ftp);
	//		getFourierDescriptor(&center_evoplist, &center_FD,0);//计算傅里叶descriptor
	//		//分别合并patch,并计算新合并的FD,之后比较pindex & combinepatch , neighpatch & combinepatch,选其中较小的similarity值，记录下neighborpatch 的index
	//		double min_similarity = 65535;
	//		int combine_pindex = -1;
	//		vector<mypoint2f> record_complist;
	//		for (unsigned int j = 0; j < adjpch.size() - 1; ++j)//adjpch最后一个是中心patch: pindex
	//		{
	//			vector<mypoint2f> combine_plist;//没化简之前的pointlist
	//			combineTwoPatch_getplist(pindex, adjpch[j], &combine_plist, 1);  //新patch的pointlist，front()！= back()，间隔1个取一个点
	//			if (!combine_plist.empty())
	//			{
	//				mypoint2f tmp_front = combine_plist.front();
	//				combine_plist.push_back(tmp_front);   //front()= back()
	//				double peri_comb = getPerimeterofVec(&combine_plist);//周长
	//				combine_plist.pop_back();  //front()！= back()
	//				vector<double> pch_FD;  //新patch的FD
	//				//getFDFeature(&combine_plist, peri_comb, reduction_num, 0, &pch_FD);
	//				vector<mypoint2f> evolu_plist;
	//				if (combine_plist.size() < reduction_num)
	//				{
	//					dataIncreasing(&combine_plist, reduction_num);
	//					dataReduction(&combine_plist, &evolu_plist, peri_comb, reduction_num); //front() != back()化简得evolu_plist
	//				}
	//				else
	//				{
	//					dataReduction(&combine_plist, &evolu_plist, peri_comb, reduction_num);//front()!=back()化简得evolu_plist
	//				}
	//				mypoint2f frontp = evolu_plist.front();//添加一个头，使evol_plist的front()=back()
	//				evolu_plist.push_back(frontp);	
	//				getFourierDescriptor(&evolu_plist, &pch_FD, 0);//计算傅里叶descriptor		
	//				double dis1 = similarityMeasure(&center_FD, &pch_FD);//similarityMeasure(&four_des.back(), &pch_FD);//pindex & combinepatch
	//				//outfile1 << dis1;
	//				//outfile1 << " ";
	//				
	//				if (dis1 < min_similarity)
	//				{
	//					min_similarity = dis1;
	//					combine_pindex = adjpch[j];
	//					evolu_plist.pop_back();
	//					record_complist = evolu_plist;// combine_plist;//记录下来这个合并的且化简的patch的plist,
	//				}
	//				/*
	//					vector<mypoint2f> plist_adj = patch_plistnum[adjpch[j]];
	//					frontp = plist_adj.front();//使front()=back()
	//					plist_adj.push_back(frontp);
	//					vector<double> FD_adjpch;
	//					getFourierDescriptor(&plist_adj, &FD_adjpch,0);//计算傅里叶descriptor	
	//					double dis2 = similarityMeasure(&FD_adjpch, &pch_FD);
	//					if (dis2 < min_similarity )
	//					{
	//						if (dis1 > min_similarity)
	//						{
	//							min_similarity = dis2;
	//							combine_pindex = adjpch[j];
	//							evolu_plist.pop_back();
	//							record_complist = evolu_plist;
	//						}
	//						else
	//							min_similarity = dis2;
	//					}*/
	//			}
	//		}
	//		//outfile1 << "\n";
	//		outfile1 << min_similarity;
	//		outfile1 << " ";
	//		if (combine_pindex != -1 && min_similarity <= similarity_t)
	//		{
	//			//attachLineProcess(int pchindex)
	//			record_similarity.push_back(min_similarity);
	//			//combineTwoPatch(pindex, combine_pindex);//留pindex，或删除combine_pindex
	//			combineTwoPatch(combine_pindex, pindex);//留combine_pindex ，删除pindex //从堆顶删除的patch要和合并掉的patch（typr变为-1）同时！！
	//			//更新新合并的patch的pointlist的数据点！！
	//			patch_plistnum[combine_pindex].swap(vector<mypoint2f>());
	//			patch_plistnum[combine_pindex] = record_complist; //front()!=back()
	//			//1.先删除堆顶的pindex
	//			int tmp_back = heap_vec.back();
	//			int tmp_mini = heap_vec[0];
	//			tt_w[tmp_mini] = -1;
	//			ele_index[tmp_mini] = -1;
	//			ele_index[tmp_back] = 0;//!!!
	//			heap_vec[0] = tmp_back;
	//			//heap_vec.pop_back();
	//			heap_siftdown(1, &heap_vec, &tt_w, &ele_index);//下沉调整
	//			heap_vec.pop_back();
	//			h_size--;
	//			//2.更新新合并的patch的weight,和在堆中的位置
	//			double feature_sum = 0;
	//			//	//1、圆度	
	//			//double circularity = 4 * 3.1415926*plrGraph.patches[combine_pindex].area / pow(plrGraph.patches[combine_pindex].perimeter, 2);
	//			//feature_sum = feature_sum + circularity;
	//			//	//2、长度，方法二，计算出最大矩形bounding box
	//			//double elongation = 0;
	//			//vector<mypoint2f> plist_box;
	//			//mypoint2f rec_center(0, 0);
	//			//getBoundRectangle(&record_complist, rec_center, &plist_box);//用化简之hou的plist计算bounding box
	//			//double box_w = fabs(plist_box[0].x - plist_box[1].x);
	//			//double box_h = fabs(plist_box[1].y - plist_box[2].y);
	//			//if (box_w < box_h)
	//			//	elongation = box_w / box_h;
	//			//else
	//			//	elongation = box_h / box_w;
	//			//feature_sum = feature_sum + elongation;
	//			//	//3、面积比，object_area/min-bounding rectangle_area,  convex or concave?
	//			//double are_rect = box_w*box_h;
	//			//double are_ratio = plrGraph.patches[combine_pindex].area / are_rect;
	//			//feature_sum = feature_sum + are_ratio;			
	//			vector<double> new_ft;
	//			getThreeFeature(&new_ft, &record_complist);//用这个函数替代上面的代码
	//			for (unsigned int i_ft = 0; i_ft < new_ft.size(); ++i_ft)
	//			{
	//				feature_sum = feature_sum + new_ft[i_ft];
	//			}
	//			patch_tmpweight[combine_pindex].second = new_ft;
	//			//tt_w[combine_pindex] = feature_sum;
	//			if (ele_index[combine_pindex] == -1)//说明combine_pindex在之前 从堆顶删除 但是没有合并（变-1），再把他加进heap中
	//			{
	//				tt_w[combine_pindex] = feature_sum;
	//				heap_vec.push_back(combine_pindex);
	//				ele_index[combine_pindex] = heap_vec.size()-1;
	//				heap_siftup(heap_vec.size(), &heap_vec, &tt_w, &ele_index);
	//			}
	//			else
	//				heap_update(ele_index[combine_pindex] + 1, feature_sum, &heap_vec, &tt_w, &ele_index);//注意+1
	//		}
	//		else
	//		{
	//			//similarity 大于阈值，不需要合并，从堆顶端删除	
	//			int tmp_back = heap_vec.back();
	//			int tmp_mini = heap_vec[0];
	//			tt_w[tmp_mini] = -1;
	//			ele_index[tmp_mini] = -1;
	//			ele_index[tmp_back] = 0;//!!!
	//			heap_vec[0] = tmp_back;
	//			heap_siftdown(1, &heap_vec, &tt_w, &ele_index);//下沉调整
	//			heap_vec.pop_back();
	//			h_size--;
	//		}
	//	}
	//	else
	//	{
	//		//没有neighbor,不需要合并，从堆顶端删除	
	//		int tmp_back = heap_vec.back();
	//		int tmp_mini = heap_vec[0];
	//		tt_w[tmp_mini] = -1;
	//		ele_index[tmp_mini] = -1;
	//		ele_index[tmp_back] = 0;//!!!
	//		heap_vec[0] = tmp_back;
	//		//heap_vec.pop_back();
	//		heap_siftdown(1, &heap_vec, &tt_w, &ele_index);//下沉调整
	//		heap_vec.pop_back();
	//		h_size--;
	//	}
	//}
	//cout <<"combines patch size="<< record_similarity.size() << endl;
	//outfile1.close();

	return retn_num;
}
//计算两个patch合并后的plist
void MyGraph::combineTwoPatch_getpvec(int pch1, int pch2, vector<int>* totalpvec)
{
	vector<int> record_addp;//记录该相邻面的除共有点外其余顶点的index
	vector<int> remove_p;
	//找到两个patch公共的边
	vector<int> pindex_origi = plrGraph.patches[pch1].pointIndexOfPatch;
	pindex_origi.insert(pindex_origi.end(), plrGraph.patches[pch1].pointIndexOfPatch.begin() + 1, plrGraph.patches[pch1].pointIndexOfPatch.end() - 1);
	vector<int> pindex_adj = plrGraph.patches[pch2].pointIndexOfPatch;
	pindex_adj.insert(pindex_adj.end(), plrGraph.patches[pch2].pointIndexOfPatch.begin() + 1, plrGraph.patches[pch2].pointIndexOfPatch.end() - 1);
	reverse(pindex_adj.begin(), pindex_adj.end());
	int fromp, endp;
	bool succflag = getAddpFromAdjpatch(pindex_origi, pindex_adj, record_addp, fromp, endp);//record_addp不包括起始端点,返回true/false是否可以合并，有些虽然可以找到commonpoint，但是。。。
	if (succflag == false)
		return;
	vector<int> pindex_pch1 = plrGraph.patches[pch1].pointIndexOfPatch;
	//vector<int> pindex_pch2 = plrGraph.patches[pch2].pointIndexOfPatch;
	vector<int>::iterator iter_pi = pindex_pch1.begin();// plrGraph.patches[pch1].pointIndexOfPatch.begin();
	if (*iter_pi != fromp)//若不是fromp开头，则对pointIndexOfPatch调整顺序,调整到以fromp开头的顺序
	{
		pindex_pch1.pop_back();//plrGraph.patches[pch1].pointIndexOfPatch.pop_back();
		while (*iter_pi != fromp)
		{
			int tmp = *iter_pi;
			iter_pi = pindex_pch1.erase(iter_pi); //plrGraph.patches[pch1].pointIndexOfPatch.erase(iter_pi);
			pindex_pch1.push_back(tmp); // plrGraph.patches[pch1].pointIndexOfPatch.push_back(tmp);
		}
		pindex_pch1.push_back(fromp); //plrGraph.patches[pch1].pointIndexOfPatch.push_back(fromp);
		iter_pi = pindex_pch1.begin(); //plrGraph.patches[pch1].pointIndexOfPatch.begin();//注意最后修改iter_pi还是指向begin!!
	}
	//此时iter_pi指向pointIndexOfPatch的begin(),然后将这些点(record_addp)加入到第i个(pch1)patch中的pointIndexOfPatch的fromp（即begin之后）之后！！		
	if (!record_addp.empty())
	{
		reverse(record_addp.begin(), record_addp.end());//getAddpFromAdjpatch()之前reverse过，再翻转回来	
		iter_pi = pindex_pch1.insert(iter_pi + 1, record_addp.begin(), record_addp.end());//insert(..begin()，val)是在第一个元素之前插入val，并返回val的第一个位置
		iter_pi = iter_pi + record_addp.size();
	}
	else
	{
		//没有要添加的点，则iter_pi从begin向前移一个位置
		iter_pi++;
	}
	//添加完了点，之后就是删除公共点
	while (iter_pi != pindex_pch1.end())//while (iter_pi != plrGraph.patches[pch1].pointIndexOfPatch.end())
	{
		if (*iter_pi != endp)
		{
			remove_p.push_back(*iter_pi);//记录被删除的点，公共点，除了两头
			iter_pi = pindex_pch1.erase(iter_pi);//plrGraph.patches[pch1].pointIndexOfPatch.erase(iter_pi);//执行这一步，说明公共边不止一个
		}
		else
			break;
	}
	for (unsigned int i = 0; i < pindex_pch1.size(); ++i)
	{
		totalpvec->push_back(pindex_pch1[i]);
	}
}
//输入pointlist（化简之后的same datapoint），输出相应的fourier descriptor (fd_type=0计算的是getbendAngleFunction，=1计算的FD)
void MyGraph::getFourierDescriptor(vector<mypoint2f>* contour, vector<double>* FD,int fd_type)
{
	/*fd_type = 0, 用bendangle function然后等采样计算FD， = 1用radial function计算FD. 注意contour是化简之后的plist,且 front()=back() */
	if (fd_type == 0)//5.0 
	{
		/*方法一 ,计算bend angle function，时间bendAngle_vec=contour*2,  (contour*2)^2*/
		// contour是化简之后的plist,且 front()=back()
		vector<vector<double>> bendAngle_vec;//[x,y] x累计边长， y旋转角
		getbendAngleFunction(contour, &bendAngle_vec, contour->size() * 2);//先在简化的contour的基础上计算，再等距采样的得 bendAngle_vec
		//计算fourier coefficient
		for (unsigned int j = 0; j < bendAngle_vec.size(); ++j)
		{
			double an_sum = 0;
			double bn_sum = 0;
			for (unsigned int k = 0; k < bendAngle_vec.size(); ++k)
			{
				double t_ff = 2 * 3.1415926*(j + 1)*bendAngle_vec[k][0];
				double ttmp = bendAngle_vec[k][1] * sin(t_ff);
				an_sum = an_sum + ttmp;
				ttmp = bendAngle_vec[k][1] * cos(t_ff);
				bn_sum = bn_sum + ttmp;
			}
			an_sum = an_sum / ((j + 1)*(-3.1415926));
			bn_sum = bn_sum / ((j + 1)*3.1415926);
			double tt = sqrt(pow(an_sum, 2) + pow(bn_sum, 2));
			FD->push_back(tt);// An_pch.push_back(tt);
		}
	}
	else//0.06 ,0.07 ,0.07
	{
		/*方法二，front()=back().... radial function,时间 contour^2 */
		int c_size = contour->size();
		double c_x = 0;
		double c_y = 0;
		for (int i = 0; i < c_size; ++i)
		{
			c_x = c_x + contour->at(i).x;
			c_y = c_y + contour->at(i).y;
		}
		mypoint2f ct_p(c_x / c_size, c_y / c_size);
		vector<double> radial_vec;
		/*（现）contour是简化后的，但可能不是均匀分布的，所以类似上面的方法，等距重采样，得到两倍的point(contour->size() * 2)*/
		//double contour_per = getPerimeterofVec(contour);
		//double inter_len = contour_per * 1 / ((double)contour->size() * 2);//两点之间插入contour->size() * 2-2个点
		//double accumu_interval = inter_len;
		////int resample_num = contour->size() * 2 - 1;
		//vector<mypoint2f> equa_contour;
		//for (unsigned int i = 1; i < contour->size(); ++i)
		//{
		//	double tmp_l = LineSegmentLength(contour->at(i - 1), contour->at(i));
		//	mypoint2f dic_vec=contour->at(i) - contour->at(i - 1);
		//}
		/*（原）计算radial destance，存入radial_vec*/
		for (int i = 0; i < c_size; ++i)
		{
			double tmp_l = LineSegmentLength(contour->at(i), ct_p);
			radial_vec.push_back(tmp_l);
		}
		/*Fourier transform 求出a0-an，算出来是复数*/
		vector<complex<double>> fd_comp;
		complex<double> ii(0, 1);//复数 j
		for (int i = 0; i < c_size; ++i)
		{
			complex<double> sum_cmp = 0;
			for (int j = 0; j < c_size; ++j)
			{
				double tmp = (-1) * 2 * 3.1415926*i*j / c_size;
				complex<double> cmp(tmp, 0);
				complex<double> cmp_2 = exp(ii*cmp);
				complex<double> cmp_rt(radial_vec[j], 0);
				complex<double> cmp_3 = cmp_rt*cmp_2;
				sum_cmp = sum_cmp + cmp_3;
			}
			complex<double> cmpsize(c_size, 0);
			sum_cmp = sum_cmp / cmpsize;
			fd_comp.push_back(sum_cmp);
		}
		//准换成实数，并标准化,1 - N/2
		double norm_len = sqrt(pow(fd_comp[0].real(), 2) + pow(fd_comp[0].imag(), 2));
		int half_size = fd_comp.size() / 2+1;
		for (int j = 1; j < half_size; ++j)//从1开始到N/2！！！
		{
			double fd_len = sqrt(pow(fd_comp[j].real(), 2) + pow(fd_comp[j].imag(), 2)) / norm_len;
			FD->push_back(fd_len);
		}
	}
}
//计算两个向量之间的欧氏距离，可以用GeometryHeader.h中的 getEuclidDistance() 代替
double MyGraph::similarityMeasure(vector<double>* fd1, vector<double>* fd2)
{
	if (fd1->size() != fd2->size())
		return 0;
	double diff = 0;
	for (unsigned int i = 0; i < fd1->size(); ++i)
	{
		diff = diff + pow((fd1->at(i)- fd2->at(i)), 2);
	}
	diff = sqrt(diff);
	return diff;
}

void MyGraph::patchToCurve(int pchindex, vector<mypoint2f>* newcurve)
{
	vector<int> curve1 = plrGraph.patches[pchindex].pointIndexOfPatch;
	vector<mypoint2f> total_plist;
	jointMethod(&curve1, &total_plist, 1);
	//for (unsigned int j = 1; j < curve1.size(); ++j)
	//{
	//	int verindex1 = curve1.at(j - 1);
	//	int verindex2 = curve1.at(j);
	//	mypoint2f point1 = sgraph.nodeList[verindex1].position;//=lineWithCommonPoints[verindex1].first
	//	mypoint2f point2 = sgraph.nodeList[verindex2].position;
	//	if (!total_plist.empty())
	//		total_plist.pop_back();
	//	getPListOfaCurve(point1, point2, verindex1, verindex2, &total_plist, 1);
	//}
	cubicBezier bz;
	fitCubBezier(&total_plist, newcurve, &bz, 1);
	
	//新增newcurve
	PolyLine pl;
	pl.origi_index = polyline_vec.size();
	pl.origi_type = 'L';
	pl.pointlist = *newcurve;
	polyline_vec.push_back(pl);
	VertexNode nnode1;
	nnode1.verIndex = sgraph.nodeList.size();
	nnode1.position = newcurve->front();
	nnode1.firstEdge = NULL;
	sgraph.nodeList.push_back(nnode1);
	VertexNode nnode2;
	nnode2.verIndex = sgraph.nodeList.size();
	nnode2.position = newcurve->back();
	nnode2.firstEdge = NULL;
	sgraph.nodeList.push_back(nnode2);

	ArcNode* arc1 = new ArcNode;
	arc1->adjVertex = nnode2.verIndex;
	arc1->curveType = 'L';
	arc1->curveIndex = polyline_vec.size() - 1;
	arc1->lineFlag = 2;
	arc1->next = sgraph.nodeList[nnode1.verIndex].firstEdge;
	sgraph.nodeList[nnode1.verIndex].firstEdge = arc1;
	ArcNode* arc2 = new ArcNode;
	arc2->adjVertex = nnode1.verIndex;
	arc2->curveType = 'L';
	arc2->curveIndex = polyline_vec.size() - 1;
	arc2->lineFlag = 2;
	arc2->next = sgraph.nodeList[nnode2.verIndex].firstEdge;
	sgraph.nodeList[nnode2.verIndex].firstEdge = arc2;
	//new curve添加到outline中
	ICurve ncurve;
	ncurve.type = 0;
	vector<int> pp;
	pp.push_back(nnode1.verIndex);
	pp.push_back(nnode2.verIndex);
	ncurve.pt_vec = pp;
	plrGraph.attachlines.push_back(ncurve);
}


//Method 2  结果在文件夹m2
double MyGraph::similarityRange(int final_pnum, double similarity_t, int simply_num, int file_num)//这个方法不适合用FD
{
	//1、圆度，公式4*pi*area/permeter^2. 如果是圆，改特征值=1
	//2、长度(elongation)、(离心率eccentricity)、length of major axis/length of minor axis	
	cout << "similarityRange()" << endl;
	vector<double> cir_vec;
//	vector<double> long_vec;
	vector<double> area_vec;

	ofstream outfile_width("E:/svgfile/m2/patch_width.txt");//"E:/record_similarity.txt"
	outfile_width.is_open();
	vector<double> aver_w;
	vector<double> long_vec;

	int pchsize = plrGraph.patches.size();
	vector<pair<int, vector<double>>> patch_tmpfeature;//int 对应patchindex，vector<double>对应特征向量
	vector<vector<double>> patch_FDs;
	for (int i = 0; i < pchsize; ++i)
	{
		vector<double> tmp;
		patch_tmpfeature.push_back(make_pair(i, tmp));//赋初值 三个feature
		patch_FDs.push_back(tmp);                    //赋初值 FD
		if (plrGraph.patches[i].type >= 0)
		{
			vector<mypoint2f> total_plist;
			jointMethod(&plrGraph.patches[i].pointIndexOfPatch, &total_plist, 0); //total_plist首 = 尾
			total_plist.pop_back();

			//1、圆度
			double circularity = 4 * 3.1415926*plrGraph.patches[i].area / pow(plrGraph.patches[i].perimeter, 2);
			//2、长度，方法二，计算出最大矩形bounding box
			double elongation = 0;
			vector<mypoint2f> plist_box;
			mypoint2f rec_center(0, 0);
			getBoundRectangle(&total_plist, rec_center, &plist_box);
			double box_w = fabs(plist_box[0].x - plist_box[1].x);
			double box_h = fabs(plist_box[1].y - plist_box[2].y);
			if (box_w < box_h)
			{
				elongation = box_w / box_h;
				outfile_width << box_w << ' ';
				aver_w.push_back(box_w); 
			}
			else
			{
				elongation = box_h / box_w;
				outfile_width << box_h << ' ';
				aver_w.push_back(box_h);
			}
			long_vec.push_back(elongation);
			//3、面积比，object_area/min-bounding rectangle_area,  convex or concave?
			double are_rect = box_w*box_h;
			double are_ratio = plrGraph.patches[i].area / are_rect;
			vector<double> fvec;//特征向量
			fvec.push_back(circularity);//圆度
			fvec.push_back(elongation);//长度
			fvec.push_back(are_ratio);//面积
			patch_tmpfeature[i].second = fvec;
			
			//vector<double> pchFD;
			////vector<mypoint2f> evolu_plist;
			////if (total_plist.size() < simply_num)
			////{
			////	dataIncreasing(&total_plist, simply_num);
			////	dataReduction(&total_plist, &evolu_plist, plrGraph.patches[i].perimeter, simply_num); //front() != back()化简得evolu_plist
			////}
			////else
			////{
			////	dataReduction(&total_plist, &evolu_plist, plrGraph.patches[i].perimeter, simply_num);//front()!=back()化简得evolu_plist
			////}
			////mypoint2f frontp = evolu_plist.front();//添加一个头，使evol_plist的front()=back()
			////evolu_plist.push_back(frontp);
			////getFourierDescriptor(&evolu_plist, &pchFD,0);//计算傅里叶descriptor		
			//getFDFeature(&total_plist, plrGraph.patches[i].perimeter, simply_num, 0, &pchFD);
			//patch_FDs[i] = pchFD;
		}
	}
	outfile_width.close();
	double retn_num = GetAverange(aver_w);//计算平均宽度
	double aver_long = GetAverange(long_vec);
	//计算每个patch与neighbor合并后的特征值, dis(center_ft,new_ft)->patchTopology.allpatches[i].arc.weight中 (//patch_subgraphs.push_back(make_pair(i, "tree"));  )
	vector<double> tt_w;//记录权值
	vector<int> weight_pindex;//记录某个patch的最小weight对应的patchindex
	for (int i = 0; i < pchsize; ++i)
	{
		tt_w.push_back(-1);
		weight_pindex.push_back(-1);
		if (patchTopology.allpatches[i].patchType >= 0)	//if (plrGraph.patches[i].type >= 0)
		{
			PatchArc* parc = patchTopology.allpatches[i].firstAdjPatch;
			double min_dis = 6553500;
			int min_pindex = -1;
			while (parc != NULL)//每一个arc 是一个neighbor patch, 中心patch与neighborpatch合并，然后计算合并后的feature，然后计算该newfeature与patch_tmpfeature[i].second的二范数，存入arc 的weight
			{
				if (parc->weight > 0)//
				{
					if (parc->weight < min_dis)
					{
						min_dis = parc->weight;
						min_pindex = parc->adj_patchindex;
					}
					parc = parc->nextPatch;
				}
				else
				{
					vector<int> combine_pvec;//合并的pointlist
					combineTwoPatch_getpvec(i,parc->adj_patchindex, &combine_pvec);//combine_plist  front()！= back()
					if (!combine_pvec.empty())//是NULL的话说明不能合并
					{
						////计算新patch的周长+FD
						//vector<double> new_ft;											
						//mypoint2f tmp_front = combine_plist.front();
						//combine_plist.push_back(tmp_front);   //front()= back()
						//double peri_comb = getPerimeterofVec(&combine_plist);//周长
						//combine_plist.pop_back();  //front()！= back()	
						//getFDFeature(&combine_plist, peri_comb, simply_num, 0, &new_ft);
						//double dis1 = getEuclidDistance(&patch_FDs[i], &new_ft);	
						//parc->weight = dis1;//parc->weight初始化是赋值为0，在选取最小值的注意把0排除!			
						//double dis2 = getEuclidDistance(&patch_FDs[parc->adj_patchindex], &new_ft);
						//计算三个ft

						vector<mypoint2f> combine_plist;//合并的pointlist
						jointMethod(&combine_pvec, &combine_plist, 0);
						combine_plist.pop_back();
						vector<double> new_ft;	
						getThreeFeature(&new_ft, &combine_plist,0,0);  //注意combine_plist front()!=back()
						double dis1 = getEuclidDistance(&patch_tmpfeature[i].second, &new_ft);
						double dis2 = getEuclidDistance(&patch_tmpfeature[parc->adj_patchindex].second, &new_ft);
						if (dis1<dis2)
							parc->weight = dis1;   //dis1较小
						else
							parc->weight = dis2;   //dis2较小

						PatchArc* parc_reve = getPatchArc(parc->adj_patchindex, i);//patchTopology.allpatches[parc->adj_patchindex].firstAdjPatch;
						if (parc_reve != NULL)
						{
							parc_reve->weight = parc->weight;
						}
						if (parc->weight < min_dis)
						{
							min_dis = parc->weight;
							min_pindex = parc->adj_patchindex;
						}
						parc = parc->nextPatch;
					}
					else
					{
						int tmp_dele = parc->adj_patchindex; 
						parc = parc->nextPatch;
						removePatchNode(tmp_dele, i);
						removePatchNode(i, tmp_dele);
					}
				}
				//parc = parc->nextPatch;
			}
			if (min_dis== 6553500)//patch 没有neighbor，
			{			
			}
			tt_w[i]=min_dis;
			weight_pindex[i] = min_pindex;
		}
	}	
	//排序
	//sort(tt_w.begin(), tt_w.end(), mycompare);//从小到大排序
	//建堆，tt_w保存了从第0个点到底n个点weight值，堆heap_vec保存每个patch的权值在tt_w中的下标，也是patch 的index
	vector<int> heap_vec;
	vector<int> ele_index;
	for (unsigned int i = 0; i < tt_w.size(); ++i)
	{
		ele_index.push_back(-1);
		if (plrGraph.patches[i].type >= 0)
		{
			heap_vec.push_back(i);
		}
		//heap_vec.push_back(i);
		//ele_index.push_back(0);
	}
	//建堆
	int h_size = heap_vec.size();
	for (int i = h_size / 2; i >= 1; i--)
	{
		heap_siftdown(i, &heap_vec, &tt_w, &ele_index);
	}
	//建立索引，知道某个权值的下标，可以直接找到在heap的位置，位置是从1开始，注意更新
	for (int i = 0; i < h_size; ++i)
	{
		ele_index[heap_vec[i]] = i;// ele_index
	}

	string file_str = "E:/svgfile/m2/similarity_t";
	file_str.append(to_string(file_num));
	file_str.append(".txt");
	ofstream outfile1(file_str);//"E:/record_similarity.txt"
	outfile1.is_open();

	//h_size是目前堆中的patch的数量，每次合并一个patch,h_size--，这是以数量为停止条件
	int stop_num = h_size *((double)final_pnum / 100.0);
	while (h_size >stop_num)//while (!heap_vec.empty())  final_pnum
	{
		int pindex = heap_vec[0];
		if (tt_w[pindex] == 6553500)
		{
			cout << "left:" << h_size << " ,stop_num " << stop_num << endl;
			break;//若final_pnum设置的过小，使weight=65535，说明剩下的patch都没有neighbor,可以直接跳出循环，此时最小weight对应的patchindex， weight_pindex[pindex]=-1
		}
		else if (tt_w[pindex] > similarity_t)
		{
			//可以直接break 跳出循环，因为剩下的都是大于similarity_t的patch
			cout << "left:" << h_size << " ,stop_num " << stop_num << endl;
			break;
		}
		int combine_pindex=-1;
		PatchArc* parc = patchTopology.allpatches[pindex].firstAdjPatch;
		//while (parc != NULL)
		//{
		//	if (parc->weight == tt_w[pindex])
		//	{
		//		combine_pindex = parc->adj_patchindex;				
		//		break;
		//	}
		//	parc = parc->nextPatch;
		//}
		//int ttt = weight_pindex[pindex];
		//if (ttt != combine_pindex)
		//	cout << "error";
		combine_pindex = weight_pindex[pindex];//第pindex个patch最小权值对应的那个patch的index
		if (combine_pindex != -1)
		{
			outfile1 << tt_w[pindex];
			outfile1 << " ";
			bool succ=combineTwoPatch(combine_pindex, pindex);//combine_pindex的周长，面积，center也会更新
			if (succ == false)
			{
				int tmp_back = heap_vec.back();
				int tmp_mini = heap_vec[0];
				tt_w[tmp_mini] = -1;
				ele_index[tmp_mini] = -1;
				ele_index[tmp_back] = 0;//!!!
				heap_vec[0] = tmp_back;
				//heap_vec.pop_back();
				heap_siftdown(1, &heap_vec, &tt_w, &ele_index);//下沉调整
				heap_vec.pop_back();
				h_size--;
				continue;
			}
		}
		else
			cout << "error" << endl;
		//1.先删除堆顶的pindex
		int tmp_back = heap_vec.back();
		int tmp_mini = heap_vec[0];
		tt_w[tmp_mini] = -1;
		ele_index[tmp_mini] = -1;
		ele_index[tmp_back] = 0;//!!!
		heap_vec[0] = tmp_back;
		heap_siftdown(1, &heap_vec, &tt_w, &ele_index);//下沉调整
		heap_vec.pop_back();
		h_size--;
		//2.更新combine_patch（combine_pindex）的三个特征
		vector<mypoint2f> total_plist;
		jointMethod(&plrGraph.patches[combine_pindex].pointIndexOfPatch, &total_plist, 0); //total_plist首 = 尾
		total_plist.pop_back();
			//1、圆度
		double circularity = 4 * 3.1415926*plrGraph.patches[combine_pindex].area / pow(plrGraph.patches[combine_pindex].perimeter, 2);
			//2、长度，计算出最大矩形bounding box
		double elongation = 0;
		vector<mypoint2f> plist_box;
		mypoint2f rec_center(0, 0);
		getBoundRectangle(&total_plist, rec_center, &plist_box);
		double box_w = fabs(plist_box[0].x - plist_box[1].x);
		double box_h = fabs(plist_box[1].y - plist_box[2].y);
		if (box_w < box_h)
			elongation = box_w / box_h;
		else
			elongation = box_h / box_w;
			//3、面积比，object_area/min-bounding rectangle_area,  convex or concave?
		double are_rect = box_w*box_h;
		double are_ratio = plrGraph.patches[combine_pindex].area / are_rect;
		vector<double> fvec;//特征向量
		fvec.push_back(circularity);//圆度
		fvec.push_back(elongation);//长度
		fvec.push_back(are_ratio);//面积
		patch_tmpfeature[combine_pindex].second = fvec;

	/*	vector<double> tmppchFD;
		getFDFeature(&total_plist, plrGraph.patches[combine_pindex].perimeter, simply_num, 0, &tmppchFD);
		patch_FDs[combine_pindex].clear(); patch_FDs[combine_pindex].swap(vector<double>());
		patch_FDs[combine_pindex] = tmppchFD;*/

		//3. 把neighbor中的pindex换成combine_pindex，除了combine_pindex
		PatchArc* pindex_nei = patchTopology.allpatches[pindex].firstAdjPatch;
		while (pindex_nei != NULL)
		{
			PatchArc* tmp_arc = getPatchArc(pindex_nei->adj_patchindex, combine_pindex);
			if (tmp_arc != NULL)
			{
				removePatchNode(pindex_nei->adj_patchindex, pindex);
			}
			else
			{
				tmp_arc = getPatchArc(pindex_nei->adj_patchindex, pindex);
				tmp_arc->adj_patchindex = combine_pindex;
			}		
			pindex_nei = pindex_nei->nextPatch;
		}
		patchTopology.allpatches[pindex].patchType = -1;//!!!
		patchTopology.allpatches[combine_pindex].firstAdjPatch = NULL;//清空neighbor，重新添加
		vector<int> compch_nei = getNeighborPatches(combine_pindex, 0);
		for (unsigned int j = 0; j < compch_nei.size(); ++j)
		{
			PatchArc* neiarc = new PatchArc;
			neiarc->adj_patchindex = compch_nei[j];
			neiarc->weight = 0;
			neiarc->nextPatch = patchTopology.allpatches[combine_pindex].firstAdjPatch;
			patchTopology.allpatches[combine_pindex].firstAdjPatch = neiarc;
		}
		//4.计算combine_patch 与所有neighbor的相似度tt_w[],更新在堆中的位置
		parc = patchTopology.allpatches[combine_pindex].firstAdjPatch;
		double mini_simi = 6553500;
		int simi_pindex = -1;
		while (parc != NULL)
		{
			vector<int> combine_pvec;//合并的pointlist
			combineTwoPatch_getpvec(combine_pindex, parc->adj_patchindex, &combine_pvec);//combine_plist  front()！= back()
			if (!combine_pvec.empty())
			{
				//vector<double> newpchFD;//FD
				//mypoint2f tmp_front = combine_plist.front();
				//combine_plist.push_back(tmp_front);   //front()= back()
				//double peri_comb = getPerimeterofVec(&combine_plist);//周长
				//combine_plist.pop_back();  //front()！= back()	
				//getFDFeature(&combine_plist, peri_comb, simply_num, 0, &newpchFD);
				//double dis1 = getEuclidDistance(&newpchFD, &patch_FDs[combine_pindex]);
				//parc->weight = dis1;
				//if (dis1 < mini_simi)
				//{
				//	mini_simi = dis1;
				//	simi_pindex = parc->adj_patchindex;//
				//}
				//double dis2 = getEuclidDistance(&newpchFD, &patch_FDs[parc->adj_patchindex]);
				//PatchArc* rever_ptcharc = getPatchArc(parc->adj_patchindex, combine_pindex);
				//rever_ptcharc->weight = dis2;

				vector<mypoint2f> combine_plist;//合并的pointlist
				jointMethod(&combine_pvec, &combine_plist, 0);
				vector<double> ft_vec;//3ge特征
				getThreeFeature(&ft_vec, &combine_plist, 0, 0);
				double dis1 = getEuclidDistance(&ft_vec, &patch_tmpfeature[combine_pindex].second);
				double dis2 = getEuclidDistance(&ft_vec, &patch_tmpfeature[parc->adj_patchindex].second);
				if (dis1<dis2)
					parc->weight = dis1;
				else
					parc->weight = dis2;
				if (parc->weight < mini_simi)
				{
					mini_simi = parc->weight;
					simi_pindex = parc->adj_patchindex;//
				}
				PatchArc* rever_ptcharc = getPatchArc(parc->adj_patchindex, combine_pindex);
				rever_ptcharc->weight = parc->weight;
					
				//重新选parc->adj_patchindex的最小weight
				double arc_minweight = 6553500;
				int arc_minipindex = -1;
				PatchArc* arc_neighbor = patchTopology.allpatches[parc->adj_patchindex].firstAdjPatch;
				while (arc_neighbor != NULL)
				{
					if (arc_neighbor->weight < arc_minweight && arc_neighbor->weight!=0)
					{
						arc_minweight = arc_neighbor->weight;
						arc_minipindex = arc_neighbor->adj_patchindex;
					}
					arc_neighbor = arc_neighbor->nextPatch;
				}
				//跟新
				heap_update(ele_index[parc->adj_patchindex] + 1, arc_minweight, &heap_vec, &tt_w, &ele_index);
				weight_pindex[parc->adj_patchindex] = arc_minipindex;
					parc = parc->nextPatch;
			}
			else
			{
				int tmp_dele = parc->adj_patchindex;
				parc = parc->nextPatch;
				removePatchNode(tmp_dele, combine_pindex);
				removePatchNode(combine_pindex, tmp_dele);
				double arc_minweight = 6553500;
				int arc_minipindex = -1;
				PatchArc* arc_neighbor = patchTopology.allpatches[tmp_dele].firstAdjPatch;
				while (arc_neighbor != NULL)
				{
					if (arc_neighbor->weight < arc_minweight && arc_neighbor->weight != 0)
					{
						arc_minweight = arc_neighbor->weight;
						arc_minipindex = arc_neighbor->adj_patchindex;
					}
					arc_neighbor = arc_neighbor->nextPatch;
				}
				//跟新
				heap_update(ele_index[tmp_dele] + 1, arc_minweight, &heap_vec, &tt_w, &ele_index);
				weight_pindex[tmp_dele] = arc_minipindex;
			}
			//parc = parc->nextPatch;
		}		
		heap_update(ele_index[combine_pindex] + 1, mini_simi, &heap_vec, &tt_w, &ele_index);//注意+1	tt_w[combine_pindex] = dis1;	
		weight_pindex[combine_pindex] = simi_pindex;
	}
	//cout << "combines patch size=" << record_similarity.size() << endl;
	outfile1.close();
	return retn_num;
}
void MyGraph::getThreeFeature(vector<double>* t_ft, vector<mypoint2f> *plist, double p_area, double p_perimeter)
{
	double pcharea = p_area;
	double pmt = p_perimeter;
	if (pcharea <= 0 || pmt <= 0)//若不知道周长和面，就重新计算
	{
		if (!isPointRoughEqual(plist->front(), plist->back()))//如果前后两个点不相等，则添加一个点使首=尾（为了计算周长和面积用）
		{
			mypoint2f startp_p = plist->front();
			plist->push_back(startp_p);
		}
		pmt = getPerimeterofVec(plist);
		pcharea = GetPolygonArea(plist);		
	}
	//1、圆度
	double circularity = 4 * 3.1415926*pcharea / pow(pmt, 2);
	//2、长度，方法二，计算出最大矩形bounding box
	double elongation = 0;
	vector<mypoint2f> plist_box;
	mypoint2f rec_center(0, 0);
	getBoundRectangle(plist, rec_center, &plist_box);
	double box_w = fabs(plist_box[0].x - plist_box[1].x);
	double box_h = fabs(plist_box[1].y - plist_box[2].y);
	if (box_w < box_h)
		elongation = box_w / box_h;
	else
		elongation = box_h / box_w;
	//3、面积比，object_area/min-bounding rectangle_area,  convex or concave?
	double are_rect = box_w*box_h;
	double are_ratio = pcharea / are_rect;
	//存入 向量
	t_ft->push_back(circularity);
	t_ft->push_back(elongation);
	t_ft->push_back(are_ratio);
}
void MyGraph::getThreeFeature_forline(vector<double>* t_ft, vector<mypoint2f> *plist)
{
	double pmt = getPerimeterofVec(plist);
	double pcharea =pmt;//area=pmt*1
	pmt = pmt * 2 + 2;
	
    //1、圆度
	double circularity = 4 * 3.1415926*pcharea / pow(pmt, 2);
	//2、长度，方法二，计算出最大矩形bounding box
	double elongation = 0;
	vector<mypoint2f> plist_box;
	mypoint2f rec_center(0, 0);
	getBoundRectangle(plist, rec_center, &plist_box);
	double box_w = fabs(plist_box[0].x - plist_box[1].x);
	double box_h = fabs(plist_box[1].y - plist_box[2].y);
	if (box_w < box_h)
		elongation = box_w / box_h;
	else
		elongation = box_h / box_w;
	//3、面积比，object_area/min-bounding rectangle_area,  convex or concave?
	double are_rect = box_w*box_h;
	double are_ratio = pcharea / are_rect;
	//存入 向量
	t_ft->push_back(circularity);
	t_ft->push_back(elongation);
	t_ft->push_back(are_ratio);
}
void MyGraph::getFDFeature(vector<mypoint2f>* contour, vector<int>* p_vec, int pchindex, int fixednum, int fd_type, vector<double>* fd_feature)//有可能是一个有index的patch，也有可能是两个patch临时合并的contour
{
	//1.contour 化简  contour  front()!=back()，且contour是全的!!
	if (contour->front() == contour->back())
		contour->pop_back();
	vector<mypoint2f> evolu_plist;
	if (fd_type == 0)//band angle function 方法
	{
		if (contour->size() < fixednum)
		{
			if (contour->front() != contour->back())
				contour->push_back(contour->front());
			dataIncreasing(contour, fixednum);
			contour->pop_back();
			dataReduction(contour, p_vec, pchindex, fixednum, &evolu_plist); //contour的front() != back()，化简得evolu_plist,pchindex有可能<0,有了pchindex可以知道perimeter，没有pchindex,周长通过contour计算
		}
		else
		{
			dataReduction(contour, p_vec, pchindex, fixednum, &evolu_plist);//front()!=back()化简得evolu_plist
		}
	}
	else //fd_type=1，用另一种data Reduction 方法
	{
		//等距采点
		if (contour->size() < fixednum)
		{
			if (contour->front() != contour->back())
				contour->push_back(contour->front());
			dataIncreasing(contour, fixednum);
			contour->pop_back();
			equalSampling(contour, p_vec, pchindex, fixednum, &evolu_plist);
		}
		else
		{
			equalSampling(contour, p_vec, pchindex, fixednum, &evolu_plist);
		}
	}
	
	mypoint2f frontp = evolu_plist.front();//添加一个头，使evol_plist的front()=back()
	evolu_plist.push_back(frontp);
	patch_simplified_plist.push_back(make_pair(pchindex, evolu_plist));////////////////////////////////////////
	vector<double> pchFD;
	//2.根据fd_type 获取FD,fd_feature
	getFourierDescriptor(&evolu_plist, fd_feature, fd_type);//计算傅里叶descriptor	
	double half_fixed = fixednum / 2;
	if (fd_type==1 && fd_feature->size() != half_fixed)
		cout << "error" << endl;
}


//Method 3 结果在文件夹m3
void MyGraph::shapeSimilarityRange(int final_pnum, double similarity_t, int simply_num, int file_num)//按照，自己的形状与neighbor的形状的相似度进行排
{
	//1、圆度，公式4*pi*area/permeter^2. 如果是圆，改特征值=1
	//2、长度(elongation)、(离心率eccentricity)、length of major axis/length of minor axis	
	cout << "shapeSimilarityRange()" << endl;
	vector<double> cir_vec;
	vector<double> long_vec;
	vector<double> area_vec;
	int pchsize = plrGraph.patches.size();
	vector<pair<int, vector<double>>> patch_tmpfeature;//int 对应patchindex，vector<double>对应特征向量
	vector<vector<double>> patch_FDs;
	for (int i = 0; i < pchsize; ++i)
	{
		vector<double> tmp;
		patch_tmpfeature.push_back(make_pair(i, tmp));//付初值
		patch_FDs.push_back(tmp);
		if (plrGraph.patches[i].type >= 0)
		{
			vector<mypoint2f> total_plist;
			for (unsigned int j = 1; j < plrGraph.patches[i].pointIndexOfPatch.size(); ++j)
			{
				int verindex1 = plrGraph.patches[i].pointIndexOfPatch.at(j - 1);
				int verindex2 = plrGraph.patches[i].pointIndexOfPatch.at(j);
				mypoint2f point1 = sgraph.nodeList[verindex1].position;//=lineWithCommonPoints[verindex1].first
				mypoint2f point2 = sgraph.nodeList[verindex2].position;
				getPListOfaCurve(point1, point2, verindex1, verindex2, &total_plist, 1);
				total_plist.pop_back();//删除一个重复的点 front()!= back()
			}
			////1、圆度
			//double circularity = 4 * 3.1415926*plrGraph.patches[i].area / pow(plrGraph.patches[i].perimeter, 2);
			////2、长度，方法二，计算出最大矩形bounding box
			//double elongation = 0;
			//vector<mypoint2f> plist_box;
			//mypoint2f rec_center(0, 0);
			//getBoundRectangle(&total_plist, rec_center, &plist_box);
			//double box_w = fabs(plist_box[0].x - plist_box[1].x);
			//double box_h = fabs(plist_box[1].y - plist_box[2].y);
			//if (box_w < box_h)
			//	elongation = box_w / box_h;
			//else
			//	elongation = box_h / box_w;
			////3、面积比，object_area/min-bounding rectangle_area,  convex or concave?
			//double are_rect = box_w*box_h;
			//double are_ratio = plrGraph.patches[i].area / are_rect;
			//vector<double> fvec;//特征向量
			//fvec.push_back(circularity);//圆度
			//fvec.push_back(elongation);//长度
			//fvec.push_back(are_ratio);//面积
			////cir_vec.push_back(circularity);
			////long_vec.push_back(elongation);
			////area_vec.push_back(are_ratio);
			//patch_tmpfeature[i].second = fvec;

			vector<mypoint2f> evolu_plist;
			if (total_plist.size() < simply_num)
			{
				total_plist.push_back(total_plist.front());
				dataIncreasing(&total_plist, simply_num);
				total_plist.pop_back();
				//dataReduction(&total_plist, &evolu_plist, plrGraph.patches[i].perimeter, simply_num); //front() != back()化简得evolu_plist
				dataReduction_byheap(&total_plist, simply_num, plrGraph.patches[i].perimeter, &evolu_plist);
			}
			else
			{
				dataReduction_byheap(&total_plist, simply_num, plrGraph.patches[i].perimeter, &evolu_plist);
			}
			mypoint2f frontp = evolu_plist.front();//添加一个头，使evol_plist的front()=back()
			evolu_plist.push_back(frontp);
			vector<double> pchFD;
			getFourierDescriptor(&evolu_plist, &pchFD,0);//计算傅里叶descriptor		
			patch_FDs[i]=pchFD;
		}
	}
	//计算每个patch与neighbor的distance, mini_dis(center_ft,neighbor_ft)->patchTopology.allpatches[i].arc.weight中 (//patch_subgraphs.push_back(make_pair(i, "tree"));  )
	vector<double> tt_w;//记录权值
	vector<int> weight_pindex;//记录最小权值对应的patchindex
	for (int i = 0; i < pchsize; ++i)
	{
		tt_w.push_back(-1);
		weight_pindex.push_back(-1);
		if (patchTopology.allpatches[i].patchType >= 0)	//if (plrGraph.patches[i].type >= 0)
		{
			PatchArc* parc = patchTopology.allpatches[i].firstAdjPatch;
			double min_dis = 6553500;
			int min_pindex = -1;
			while (parc != NULL)//每一个arc 是一个neighbor patch, 中心patch与neighborpatch合并，然后计算合并后的feature，然后计算该newfeature与patch_tmpfeature[i].second的二范数，存入arc 的weight
			{
				if (parc->weight > 0)
				{
					if (parc->weight < min_dis)
					{
						min_dis = parc->weight;
						min_pindex = parc->adj_patchindex;
					}
					parc = parc->nextPatch;
				}
				else
				{
					bool mgable=checkMergeable(i, parc->adj_patchindex);
					if (mgable == true)
					{
						double dis = getEuclidDistance(&patch_FDs[i], &patch_FDs[parc->adj_patchindex]);	//计算距离
						//double dis = getEuclidDistance(patch_tmpfeature[i].second, patch_tmpfeature[parc->adj_patchindex].second);	//计算距离
						parc->weight = dis;//parc->weight初始化是赋值为0，在选取最小值的注意把0排除!
						//e.g.0的neighbor5赋值weight，反过来 5的neighbor0也赋值weight
						PatchArc* parc_reve = getPatchArc(parc->adj_patchindex, i);//patchTopology.allpatches[parc->adj_patchindex].firstAdjPatch;
						if (parc_reve != NULL)
						{
							parc_reve->weight = dis;
						}
						else
							cout << "patch arc=NULL" << endl;
						if (parc->weight < min_dis)
						{
							min_dis = parc->weight;
							min_pindex = parc->adj_patchindex;
						}
						parc = parc->nextPatch;
					}
					else
					{
						int arc_index = parc->adj_patchindex;
						parc = parc->nextPatch;
						removePatchNode(i, arc_index);
						removePatchNode(arc_index, i);
					}
				}
				//parc = parc->nextPatch;
			}
			if (min_dis == 6553500)//patch 没有neighbor，
			{
			}
			tt_w[i] = min_dis;
			weight_pindex[i] = min_pindex;
		}
	}
	//建堆，tt_w保存了从第0个点到底n个点weight值，堆heap_vec保存每个patch的权值在tt_w中的下标，也是patch 的index
	vector<int> heap_vec;
	vector<int> ele_index;
	for (unsigned int i = 0; i < tt_w.size(); ++i)
	{
		ele_index.push_back(-1);
		if (plrGraph.patches[i].type >= 0)
		{
			heap_vec.push_back(i);
		}
		//heap_vec.push_back(i);
		//ele_index.push_back(0);
	}
	//建堆
	int h_size = heap_vec.size();
	for (int i = h_size / 2; i >= 1; i--)
	{
		heap_siftdown(i, &heap_vec, &tt_w, &ele_index);
	}
	//建立索引，知道某个权值的下标，可以直接找到在heap的位置，位置是从1开始，注意更新
	for (int i = 0; i < h_size; ++i)
	{
		ele_index[heap_vec[i]] = i;// ele_index
	}

	string file_str = "E:/mini_similarity";
	file_str.append(to_string(file_num));
	file_str.append(".txt");
	ofstream outfile1(file_str);//"E:/record_similarity.txt"
	outfile1.is_open();
	//vector<double> record_similarity;//输出相似度
	//int reduction_num = simply_num;//规定计算相似度的点的个数（即化间个数），若数据点过小，该值可以减少一点
	//int stop_num = h_size / 5;

	while (h_size >final_pnum)//while (!heap_vec.empty())//
	{
		int pindex = heap_vec[0];
		if (tt_w[pindex] == 65535)
		{
			break;//可以直接跳出循环，说明剩下的patch都没有neighbor  weight_pindex[pindex]=-1
		}
		else if (tt_w[pindex] > similarity_t)//0.6  10
		{
			//可以直接break 跳出循环，因为剩下的都是大于similarity_t的patch
			break;
		}
		int combine_pindex = -1;
		PatchArc* parc = patchTopology.allpatches[pindex].firstAdjPatch;
		//while (parc != NULL)
		//{
		//	if (parc->weight == tt_w[pindex])
		//	{
		//		combine_pindex = parc->adj_patchindex;				
		//		break;
		//	}
		//	parc = parc->nextPatch;
		//}
		//int ttt = weight_pindex[pindex];
		//if (ttt != combine_pindex)
		//	cout << "error";
		combine_pindex = weight_pindex[pindex];
		if (combine_pindex != -1)
		{
			outfile1 << tt_w[pindex];
			outfile1 << " ";
			bool succ = combineTwoPatch(combine_pindex, pindex);//combine_pindex的周长，面积，center也会更新
			if (succ == false)
			{
				int tmp_back = heap_vec.back();
				int tmp_mini = heap_vec[0];
				tt_w[tmp_mini] = -1;
				ele_index[tmp_mini] = -1;
				ele_index[tmp_back] = 0;//!!!
				heap_vec[0] = tmp_back;
				//heap_vec.pop_back();
				heap_siftdown(1, &heap_vec, &tt_w, &ele_index);//下沉调整
				heap_vec.pop_back();
				h_size--;
				continue;
			}
		}
		else
		{
			cout << "error" << endl;
		}
		//1.先删除堆顶的pindex
		int tmp_back = heap_vec.back();
		int tmp_mini = heap_vec[0];
		tt_w[tmp_mini] = -1;
		ele_index[tmp_mini] = -1;
		ele_index[tmp_back] = 0;//!!!
		heap_vec[0] = tmp_back;
		heap_siftdown(1, &heap_vec, &tt_w, &ele_index);//下沉调整
		heap_vec.pop_back();
		h_size--;
		//2.更新combinr_patch（combine_pindex）的三个特征
		vector<mypoint2f> total_plist;
		for (unsigned int j = 1; j < plrGraph.patches[combine_pindex].pointIndexOfPatch.size(); ++j)
		{
			int verindex1 = plrGraph.patches[combine_pindex].pointIndexOfPatch.at(j - 1);
			int verindex2 = plrGraph.patches[combine_pindex].pointIndexOfPatch.at(j);
			mypoint2f point1 = sgraph.nodeList[verindex1].position;//=lineWithCommonPoints[verindex1].first
			mypoint2f point2 = sgraph.nodeList[verindex2].position;
			getPListOfaCurve(point1, point2, verindex1, verindex2, &total_plist, 1);
			total_plist.pop_back();//删除一个重复的点 front!= back
		}
		////1、圆度
		//double circularity = 4 * 3.1415926*plrGraph.patches[combine_pindex].area / pow(plrGraph.patches[combine_pindex].perimeter, 2);
		////2、长度，计算出最大矩形bounding box
		//double elongation = 0;
		//vector<mypoint2f> plist_box;
		//mypoint2f rec_center(0, 0);
		//getBoundRectangle(&total_plist, rec_center, &plist_box);
		//double box_w = fabs(plist_box[0].x - plist_box[1].x);
		//double box_h = fabs(plist_box[1].y - plist_box[2].y);
		//if (box_w < box_h)
		//	elongation = box_w / box_h;
		//else
		//	elongation = box_h / box_w;
		////3、面积比，object_area/min-bounding rectangle_area,  convex or concave?
		//double are_rect = box_w*box_h;
		//double are_ratio = plrGraph.patches[combine_pindex].area / are_rect;
		//vector<double> fvec;//特征向量
		//fvec.push_back(circularity);//圆度
		//fvec.push_back(elongation);//长度
		//fvec.push_back(are_ratio);//面积
		//patch_tmpfeature[combine_pindex].second = fvec;

		vector<mypoint2f> evolu_plist;
		if (total_plist.size() < simply_num)
		{
			total_plist.push_back(total_plist.front());
			dataIncreasing(&total_plist, simply_num);
			total_plist.pop_back();
			//dataReduction(&total_plist, &evolu_plist, plrGraph.patches[combine_pindex].perimeter, simply_num); //front() != back()化简得evolu_plist
			dataReduction_byheap(&total_plist, simply_num, plrGraph.patches[combine_pindex].perimeter, &evolu_plist);
		}
		else
		{
			dataReduction_byheap(&total_plist, simply_num, plrGraph.patches[combine_pindex].perimeter, &evolu_plist);
		}
		mypoint2f frontp = evolu_plist.front();//添加一个头，使evol_plist的front()=back()
		evolu_plist.push_back(frontp);
		vector<double> tmppchFD;
		getFourierDescriptor(&evolu_plist, &tmppchFD,0);
		patch_FDs[combine_pindex].clear(); patch_FDs[combine_pindex].swap(vector<double>());
		patch_FDs[combine_pindex]=tmppchFD;

		//3. 把neighbor中的pindex换成combine_pindex，除了combine_pindex
		PatchArc* pindex_nei = patchTopology.allpatches[pindex].firstAdjPatch;
		while (pindex_nei != NULL)
		{
			PatchArc* tmp_arc = getPatchArc(pindex_nei->adj_patchindex, combine_pindex);
			if (tmp_arc != NULL)
			{
				removePatchNode(pindex_nei->adj_patchindex, pindex);
			}
			else
			{
				tmp_arc = getPatchArc(pindex_nei->adj_patchindex, pindex);
				tmp_arc->adj_patchindex = combine_pindex;
			}
			pindex_nei = pindex_nei->nextPatch;
		}
		patchTopology.allpatches[pindex].patchType = -1;
		patchTopology.allpatches[combine_pindex].firstAdjPatch = NULL;//清空neighbor，重新添加
		vector<int> compch_nei = getNeighborPatches(combine_pindex, 0);
		for (unsigned int j = 0; j < compch_nei.size(); ++j)
		{	
			PatchArc* neiarc = new PatchArc;
			neiarc->adj_patchindex = compch_nei[j];
			neiarc->weight = 0;
			neiarc->nextPatch = patchTopology.allpatches[combine_pindex].firstAdjPatch;
			patchTopology.allpatches[combine_pindex].firstAdjPatch = neiarc;
			//if (checkMergeable(compch_nei[j], combine_pindex) == true)
			//{
			//	PatchArc* neiarc = new PatchArc;
			//	neiarc->adj_patchindex = compch_nei[j];
			//	neiarc->weight = 0;
			//	neiarc->nextPatch = patchTopology.allpatches[combine_pindex].firstAdjPatch;
			//	patchTopology.allpatches[combine_pindex].firstAdjPatch = neiarc;
			//}
			//else
			//{
			//	PatchArc* chee = getPatchArc(compch_nei[j], combine_pindex);
			//	if (chee != NULL)
			//	{
			//		removePatchNode(compch_nei[j], combine_pindex);
			//		double arc_minweight = 65535;
			//		int arc_minipindex = -1;
			//		PatchArc* arc_neighbor = patchTopology.allpatches[compch_nei[j]].firstAdjPatch;
			//		while (arc_neighbor != NULL)
			//		{
			//			if (arc_neighbor->weight < arc_minweight && arc_neighbor->weight != 0)
			//			{
			//				arc_minweight = arc_neighbor->weight;
			//				arc_minipindex = arc_neighbor->adj_patchindex;
			//			}
			//			arc_neighbor = arc_neighbor->nextPatch;
			//		}
			//		//跟新
			//		heap_update(ele_index[compch_nei[j]] + 1, arc_minweight, &heap_vec, &tt_w, &ele_index);
			//		weight_pindex[compch_nei[j]] = arc_minipindex;
			//	}
			//}
		}
		//4.计算cumbine_patch 与所有neighbor的相似度tt_w[],即更新在堆中的位置
		parc = patchTopology.allpatches[combine_pindex].firstAdjPatch;
		double mini_simi = 6553500;
		int simi_pindex = -1;
		while (parc != NULL)
		{
			if (checkMergeable(parc->adj_patchindex, combine_pindex) == true)
			{
				double dis1 = getEuclidDistance(&patch_FDs[parc->adj_patchindex], &patch_FDs[combine_pindex]);	//计算距离
				//double dis1 = getEuclidDistance(patch_tmpfeature[parc->adj_patchindex].second, patch_tmpfeature[combine_pindex].second);
				parc->weight = dis1;
				if (dis1 < mini_simi)
				{
					mini_simi = dis1;
					simi_pindex = parc->adj_patchindex;//
				}
				PatchArc* rever_ptcharc = getPatchArc(parc->adj_patchindex, combine_pindex);
				if (rever_ptcharc == NULL)
				{
					PatchArc* neiarc = new PatchArc;
					neiarc->adj_patchindex = combine_pindex;
					neiarc->weight = dis1;
					neiarc->nextPatch = patchTopology.allpatches[combine_pindex].firstAdjPatch;
					patchTopology.allpatches[combine_pindex].firstAdjPatch = neiarc;
				}
				else
				{
					rever_ptcharc->weight = dis1;
				}
				//重新选最小weight
				double arc_minweight = 65535;
				int arc_minipindex = -1;
				PatchArc* arc_neighbor = patchTopology.allpatches[parc->adj_patchindex].firstAdjPatch;
				while (arc_neighbor != NULL)
				{
					if (arc_neighbor->weight < arc_minweight && arc_neighbor->weight != 0)
					{
						arc_minweight = arc_neighbor->weight;
						arc_minipindex = arc_neighbor->adj_patchindex;
					}
					arc_neighbor = arc_neighbor->nextPatch;
				}
				heap_update(ele_index[parc->adj_patchindex] + 1, arc_minweight, &heap_vec, &tt_w, &ele_index);
				weight_pindex[parc->adj_patchindex] = arc_minipindex;
				parc = parc->nextPatch;
			}
			else
			{
				int tmp_dele = parc->adj_patchindex;
				parc = parc->nextPatch;
				removePatchNode(tmp_dele, combine_pindex);
				removePatchNode(combine_pindex, tmp_dele);
				double arc_minweight = 65535;
				int arc_minipindex = -1;
				PatchArc* arc_neighbor = patchTopology.allpatches[tmp_dele].firstAdjPatch;
				while (arc_neighbor != NULL)
				{
					if (arc_neighbor->weight < arc_minweight && arc_neighbor->weight != 0)
					{
						arc_minweight = arc_neighbor->weight;
						arc_minipindex = arc_neighbor->adj_patchindex;
					}
					arc_neighbor = arc_neighbor->nextPatch;
				}
				//跟新
				heap_update(ele_index[tmp_dele] + 1, arc_minweight, &heap_vec, &tt_w, &ele_index);
				weight_pindex[tmp_dele] = arc_minipindex;
			}
			////parc = parc->nextPatch;
		}
		heap_update(ele_index[combine_pindex] + 1, mini_simi, &heap_vec, &tt_w, &ele_index);//注意+1	tt_w[combine_pindex] = dis1;	
		weight_pindex[combine_pindex] = simi_pindex;
	}
	outfile1.close();
}
bool MyGraph::checkMergeable(int pch1, int pch2)
{
	vector<int> record_addp;//记录该相邻面的除共有点外其余顶点的index
	vector<int> remove_p;
	//找到两个patch公共的边
	vector<int> pindex_origi = plrGraph.patches[pch1].pointIndexOfPatch;
	pindex_origi.insert(pindex_origi.end(), plrGraph.patches[pch1].pointIndexOfPatch.begin() + 1, plrGraph.patches[pch1].pointIndexOfPatch.end() - 1);
	vector<int> pindex_adj = plrGraph.patches[pch2].pointIndexOfPatch;
	pindex_adj.insert(pindex_adj.end(), plrGraph.patches[pch2].pointIndexOfPatch.begin() + 1, plrGraph.patches[pch2].pointIndexOfPatch.end() - 1);
	reverse(pindex_adj.begin(), pindex_adj.end());
	int fromp, endp;
	bool succflag = getAddpFromAdjpatch(pindex_origi, pindex_adj, record_addp, fromp, endp);//record_addp不包括起始端点,返回true/false是否可以合并，有些虽然可以找到commonpoint，但是。。。
	if (succflag == false)
		return false;
	else
		return true;
}

//对patch的单一操作
bool MyGraph::combineTwoPatch(int pch1, int pch2)//保留pch1并修改edge_vec等等，删除pch2,添加了一点处理：若消失的边上有attachline，直接删除掉？
{
	if (plrGraph.patches[pch1].type == -1)
		cout << "error" << endl;
	vector<int> record_addp;//记录该相邻面的除共有点外其余顶点的index
	vector<int> remove_p;  //公共边上的点，用这些点查看是否有attachline在上面
	//找到两个patch公共的边
	vector<int> pindex_origi = plrGraph.patches[pch1].pointIndexOfPatch;
	pindex_origi.insert(pindex_origi.end(), plrGraph.patches[pch1].pointIndexOfPatch.begin() + 1, plrGraph.patches[pch1].pointIndexOfPatch.end() - 1);
	vector<int> pindex_adj = plrGraph.patches[pch2].pointIndexOfPatch;
	pindex_adj.insert(pindex_adj.end(), plrGraph.patches[pch2].pointIndexOfPatch.begin() + 1, plrGraph.patches[pch2].pointIndexOfPatch.end() - 1);
	reverse(pindex_adj.begin(), pindex_adj.end());
	int fromp, endp;
	bool succflag = getAddpFromAdjpatch(pindex_origi, pindex_adj, record_addp, fromp, endp);//record_addp不包括起始端点
	if (succflag == false)
	{
		return false;
	}
	plrGraph.patches[pch2].type = -1;//type, patch_index, edges_vec, pointIndexOfPatch。相邻面的type变为-1（表示已merge掉）
	//plrGraph.patches[pch2].patch_index = plrGraph.patches[pch1].patch_index;//!!!!						

	vector<int>::iterator iter_pi = plrGraph.patches[pch1].pointIndexOfPatch.begin();
		if (*iter_pi != fromp)//若不是fromp开头，则对pointIndexOfPatch调整顺序,调整到以fromp开头的顺序
		{
			plrGraph.patches[pch1].pointIndexOfPatch.pop_back();
			while (*iter_pi != fromp)
			{
				int tmp = *iter_pi;
				iter_pi = plrGraph.patches[pch1].pointIndexOfPatch.erase(iter_pi);
				plrGraph.patches[pch1].pointIndexOfPatch.push_back(tmp);
			}
			plrGraph.patches[pch1].pointIndexOfPatch.push_back(fromp);
			iter_pi = plrGraph.patches[pch1].pointIndexOfPatch.begin();//注意最后修改iter_pi还是指向begin!!
		}
		//此时iter_pi指向pointIndexOfPatch的begin(),然后将这些点(record_addp)加入到第i个patch中的pointIndexOfPatch的fromp（即begin之后）之后！！		
		if (!record_addp.empty())
		{
			reverse(record_addp.begin(), record_addp.end());//getAddpFromAdjpatch()之前reverse过，再翻转回来	
			iter_pi = plrGraph.patches[pch1].pointIndexOfPatch.insert(iter_pi + 1, record_addp.begin(), record_addp.end());//insert(..begin()，val)是在第一个元素之前插入val，并返回val的第一个位置
			iter_pi = iter_pi + record_addp.size();
		}
		else
		{
			//没有要添加的点，则iter_pi从begin向前移一个位置
			iter_pi++;
		}
		//添加完了点，之后就是删除公共点。
		while (iter_pi != plrGraph.patches[pch1].pointIndexOfPatch.end())
		{
			if (*iter_pi != endp)
			{
				remove_p.push_back(*iter_pi);//记录被删除的点，公共点，除了两头
				iter_pi = plrGraph.patches[pch1].pointIndexOfPatch.erase(iter_pi);//执行这一步，说明公共边不止一个
			}
			else
				break;
		}
		//新添加的代码，remove_p记录下这些公共点，看是否有attachline相连，若有，则令plrGraph.otherlines[rmove_li].type = -1; 且删除sgraph.nodelist
		if (!remove_p.empty())
		{
			if (!plrGraph.patches[pch1].point_line.empty())
			{
				for (unsigned int ri = 0; ri < remove_p.size(); ++ri)
				{
					vector<pair<int, int>>::iterator iter_pl = plrGraph.patches[pch1].point_line.begin();
					while (iter_pl != plrGraph.patches[pch1].point_line.end())
					{
						if (iter_pl->first == remove_p[ri])
						{
							int line_index = iter_pl->second;
							plrGraph.otherlines[line_index].type = -1;
							plrGraph.attachlines[line_index].type = -1;  //mixedTopology.allnode[].status还是0， 没有被修改，，
							//vector<int> rm_line=plrGraph.otherlines[line_index].pt_vec;
							//for (size_t ti = 1; ti < rm_line.size(); ++ti)
							//{
							//	removeArcNode(rm_line[ti - 1], rm_line[ti], 0);
							//	removeArcNode(rm_line[ti], rm_line[ti - 1], 0);
							//}//removeArcNode
							iter_pl = plrGraph.patches[pch1].point_line.erase(iter_pl);
						}
						else
							iter_pl++;
					}
				}
			}
		}
		//把pch2中的attachline添加到pch1中，即从record_addp中的找attachline 添加到pch1的point_line里
		if (!plrGraph.patches[pch2].point_line.empty())
		{
			if (!record_addp.empty())
			{
				for (unsigned int i = 0; i < plrGraph.patches[pch2].point_line.size(); ++i)
				{
					int tmp_pindex = plrGraph.patches[pch2].point_line[i].first;
					int tmp_lindex = plrGraph.patches[pch2].point_line[i].second;
					//bool fflag = false;
					//for (unsigned int j = 0; j < record_addp.size(); ++j)
					//{
					//	if (record_addp[j] == tmp_pindex)
					//	{
					//		fflag = true;
					//		break;
					//	}
					//}
					//if (fflag == true)
					//{
					//	if (plrGraph.otherlines[tmp_lindex].type >= 0)
					//		plrGraph.patches[pch1].point_line.push_back(plrGraph.patches[pch2].point_line[i]);
					//}
					if (tmp_pindex != fromp && tmp_pindex != endp)
					{
						if (plrGraph.otherlines[tmp_lindex].type >= 0)//!!
							plrGraph.patches[pch1].point_line.push_back(plrGraph.patches[pch2].point_line[i]);
					}
				}
				/*检查是否有重复的添加进去
				vector<pair<int, int>>::iterator iter_pl = plrGraph.patches[pch1].point_line.begin();
				while (iter_pl != plrGraph.patches[pch1].point_line.end())
				{
				vector<pair<int, int>>::iterator iter_inn = iter_pl + 1;
				while (iter_inn != plrGraph.patches[pch1].point_line.end())
				{
				if (iter_pl->first == iter_inn->first && iter_pl->second == iter_inn->second)
				iter_inn = plrGraph.patches[pch1].point_line.erase(iter_inn);
				else
				iter_inn++;
				}
				iter_pl++;
				}*/
			}
		}
	//将adgpatch的其余edge加入到第i个patch的edges_vec中。并且修改其相邻path的edge的信息。
	record_addp.push_back(endp);
	record_addp.insert(record_addp.begin(), fromp);
	vector<EdgeOfPatch>::iterator iter_oe = plrGraph.patches[pch1].edges_vec.begin();
	if (iter_oe->startp_index != fromp)//对edges_vec调整顺序
	{
		while (iter_oe->startp_index != fromp)
		{
			EdgeOfPatch tmp_e = *iter_oe;
			iter_oe = plrGraph.patches[pch1].edges_vec.erase(iter_oe);
			plrGraph.patches[pch1].edges_vec.push_back(tmp_e);
		}
	}
	iter_oe = plrGraph.patches[pch1].edges_vec.begin();
	while (iter_oe != plrGraph.patches[pch1].edges_vec.end())
	{
		if (iter_oe->startp_index == fromp)
		{
			for (int edge_i = 1; edge_i < record_addp.size(); edge_i++)
			{
				EdgeOfPatch *ep = getBorderofPatch(pch2, record_addp[edge_i - 1], record_addp[edge_i]);
				EdgeOfPatch newedge;	
				if (ep->startp_index == record_addp[edge_i - 1])
				{
					newedge.startp_index = ep->startp_index;
					newedge.endp_index = ep->endp_index;
					newedge.edgeflag = 0;
					newedge.adjPatchIndex.push_back(plrGraph.patches[pch1].patch_index);//newedge.adjPatchIndex[0]
					if (ep->adjPatchIndex.size()>1)
					{
						//修改其相邻patch的edge的信息
						newedge.adjPatchIndex.push_back(ep->adjPatchIndex[1]);
						int ano_pindex = ep->adjPatchIndex[1];//newedge.adjPatchIndex[1]
						EdgeOfPatch* edp = getBorderofPatch(ano_pindex, newedge.endp_index, newedge.startp_index);
						if (edp->adjPatchIndex.size() > 1)
							edp->adjPatchIndex[1] = pch1;
						else
							edp->adjPatchIndex.push_back(pch1);
					}
					//修改相关node
					ArcNode *c_arc = getArcNode(ep->startp_index, ep->endp_index);
					c_arc->patchesIndex.push_back(pch1);
					vector<int>::iterator pch_iter = c_arc->patchesIndex.begin();
					while (pch_iter != c_arc->patchesIndex.end())
					{
						if (*pch_iter == pch2)
						{
							c_arc->patchesIndex.erase(pch_iter);
							break;
						}
						pch_iter++;
					}
					//反向修改
					c_arc = getArcNode(ep->endp_index, ep->startp_index);
					c_arc->patchesIndex.push_back(pch1);
					pch_iter = c_arc->patchesIndex.begin();
					while (pch_iter != c_arc->patchesIndex.end())
					{
						if (*pch_iter == pch2)
						{
							c_arc->patchesIndex.erase(pch_iter);
							break;
						}
						pch_iter++;
					}
					iter_oe = plrGraph.patches[pch1].edges_vec.insert(iter_oe, newedge);//依次累加，最后成了addedge+deleedge+origiedge
					iter_oe = iter_oe + 1;
				}
				else
					cout << "error" << endl;
			}
			//删除掉plrGraph.patches[i].edges_vec中原有的dele edge(公共边),并且删除sgraph.nodelist中的连接关系
			while (iter_oe != plrGraph.patches[pch1].edges_vec.end())
			{
				if (iter_oe->endp_index != endp)
				{
					ArcNode *arc_node = getArcNode(iter_oe->startp_index, iter_oe->endp_index);
					removeArcNode(iter_oe->startp_index, iter_oe->endp_index, 0);
					removeArcNode(iter_oe->endp_index, iter_oe->startp_index, 0);

					iter_oe = plrGraph.patches[pch1].edges_vec.erase(iter_oe);
				}
				else
					break;
			}
			if (iter_oe != plrGraph.patches[pch1].edges_vec.end())
			{
				removeArcNode(iter_oe->startp_index, iter_oe->endp_index, 0);
				removeArcNode(iter_oe->endp_index, iter_oe->startp_index, 0);

				plrGraph.patches[pch1].edges_vec.erase(iter_oe);//删除掉iter_oe->endp_index = endp的那条边
			}
			//更新周长,面积，centerpt
			vector<mypoint2f> edgeplist;
			for (unsigned int iEdge = 1; iEdge < plrGraph.patches[pch1].pointIndexOfPatch.size(); ++iEdge)
			{
				int verindex1 = plrGraph.patches[pch1].pointIndexOfPatch.at(iEdge - 1);
				int verindex2 = plrGraph.patches[pch1].pointIndexOfPatch.at(iEdge);
				mypoint2f point1 = sgraph.nodeList[verindex1].position;//=lineWithCommonPoints[verindex1].first
				mypoint2f point2 = sgraph.nodeList[verindex2].position;
				if (!edgeplist.empty())
					edgeplist.pop_back();
				getPListOfaCurve(point1, point2, verindex1, verindex2, &edgeplist,0);
			}
			double pch_perimeter = getPerimeterofVec(&edgeplist);
			plrGraph.patches[pch1].perimeter = pch_perimeter;//更新周长type
			double aver_x = 0;
			double aver_y = 0;
			for (unsigned int j = 0; j < edgeplist.size(); ++j)
			{
				aver_x = aver_x + edgeplist[j].x;
				aver_y = aver_y + edgeplist[j].y;
			}
			aver_x = aver_x / (double)edgeplist.size();
			aver_y = aver_y / (double)edgeplist.size();
			plrGraph.patches[pch1].centre_position = mypoint2f(aver_x, aver_y);//更新中心
			double pch_areas = GetPolygonArea(&edgeplist);//更新面积
			if (pch_areas >= 0)//若面积>=0 说明是正常patch，应为patch都是沿同一个顺序的（逆时针）
				plrGraph.patches[pch1].area = fabs(pch_areas);

			break;
		}
		else
			cout << "error" << endl;
		iter_oe++;
	}
	return true;
}
void MyGraph::erosionApatch(int pch)//单个操作 combine(merge) erosion(变成一条线)
{
	//getpatch pointlist,已知有两个拐点，判断两个最远点
	vector<mypoint2f> apatch;
	for (unsigned int j = 1; j < plrGraph.patches[pch].pointIndexOfPatch.size(); ++j)
	{
		int verindex1 = plrGraph.patches[pch].pointIndexOfPatch.at(j - 1);
		int verindex2 = plrGraph.patches[pch].pointIndexOfPatch.at(j);
		mypoint2f point1 = sgraph.nodeList[verindex1].position;//=lineWithCommonPoints[verindex1].first
		mypoint2f point2 = sgraph.nodeList[verindex2].position;
		if (!apatch.empty())
			apatch.pop_back();
		getPListOfaCurve(point1, point2, verindex1, verindex2, &apatch,0);
	}
	if (isPointRoughEqual(apatch.front(), apatch.back()))
		apatch.pop_back();
	mypoint2f pfront = apatch.front();
	mypoint2f pback = apatch.back();
	apatch.insert(apatch.begin(), pback);
	apatch.push_back(pfront);
	//找最远的两个点,
	int p1, p2; double fdis = 0;
	vector<int> tp_vec = plrGraph.patches[pch].pointIndexOfPatch;
	tp_vec.pop_back();
	for (unsigned int j = 0; j < tp_vec.size() - 1; ++j)
	{
		int node1 = tp_vec[j];
		for (unsigned int k = j + 1; k < tp_vec.size(); ++k)
		{
			int node2 = tp_vec[k];
			double tmp_d = LineSegmentLength(sgraph.nodeList[node1].position, sgraph.nodeList[node2].position);
			if (tmp_d > fdis)
			{
				fdis = tmp_d;
				p1 = node1;
				p2 = node2;
			}
		}
	}
	//判断是否可以腐蚀 erosion
	int ano_pch = 0;
	bool eflag = checkErosion(pch, p1, p2,ano_pch);
	if (eflag == false)
	{
		combineTwoPatch(pch, ano_pch);
	}
	else
	{
		//两条路径plist1 plist2
		vector<mypoint2f> plist1, plist2;//plist1 和plist2都是从p1到p2的，有序的
		for (unsigned int j = 0; j < plrGraph.patches[pch].pointIndexOfPatch.size() - 1; ++j)//12341+234
		{
			tp_vec.push_back(plrGraph.patches[pch].pointIndexOfPatch[j]);
		}
		for (unsigned int j = 0; j < tp_vec.size(); ++j)
		{
			if (tp_vec[j] == p1)
			{
				int k = j;
				while (tp_vec[k] != p2)
				{
					if (!plist1.empty())
						plist1.pop_back();
					getPListOfaCurve(sgraph.nodeList[tp_vec[k]].position, sgraph.nodeList[tp_vec[k + 1]].position, tp_vec[k], tp_vec[k + 1], &plist1,0);
					k++;
					if (k >= tp_vec.size())
						break;
				}
				break;
			}
		}
		reverse(tp_vec.begin(), tp_vec.end());
		for (unsigned int j = 0; j < tp_vec.size(); ++j)
		{
			if (tp_vec[j] == p1)
			{
				int k = j;
				while (tp_vec[k] != p2)
				{
					if (!plist2.empty())
						plist2.pop_back();
					getPListOfaCurve(sgraph.nodeList[tp_vec[k]].position, sgraph.nodeList[tp_vec[k + 1]].position, tp_vec[k], tp_vec[k + 1], &plist2,0);
					k++;
					if (k >= tp_vec.size())
						break;
				}
				break;
			}
		}
		//计算中间路径mid_plist, 存入临时变量 addplist; 如果patch非常小，路径可能就4，5个点，需注意一下
		vector<vector<mypoint2f>> addplist;
		vector<mypoint2f> mid_plist;
		getMidPlist(plist1, plist2, &mid_plist); //p1-->p2
		addplist.push_back(mid_plist);
		//删除patch，并且删除patch边界处de边sgraph.nodelist.firstedge->delete   sgraph.nodelist.orderarc-->delete
		plrGraph.patches[pch].type = -1;
		deleOverlayPatches(plrGraph.patches[pch].patch_index, 0);
		for (unsigned int j = 0; j < plrGraph.patches[pch].edges_vec.size(); ++j)
		{
			int edge_s = plrGraph.patches[pch].edges_vec[j].startp_index;
			int edge_e = plrGraph.patches[pch].edges_vec[j].endp_index;
			ArcNode* anode = sgraph.nodeList[edge_s].firstEdge;
			while (anode)
			{
				if (anode->adjVertex == edge_e)
				{
					if (anode->patchesIndex.empty() && anode->lineFlag == 0)
					{
						removeArcNode(edge_s, edge_e, 0);//0表示真删除
						removeArcNode(edge_e, edge_s, 0);
					}
					break;
				}
				anode = anode->next;
			}
		}
		//修改与已删除patch相邻的patch的线.先找到非p1p2的端点，这些端点相连的线又分为3种情况处理：有patch相连，无patch相连（outline or innerline）
		vector<int> leftpoints;//leftpoints 现在是无序的
		for (unsigned int j = 0; j < plrGraph.patches[pch].pointIndexOfPatch.size() - 1; ++j)
		{
			if (plrGraph.patches[pch].pointIndexOfPatch[j] != p1 &&plrGraph.patches[pch].pointIndexOfPatch[j] != p2)
			{
				leftpoints.push_back(plrGraph.patches[pch].pointIndexOfPatch[j]);
			}
		}
		if (!leftpoints.empty())
		{
			//计算leftpoint到midplist的投影点(即距离最近的点)
			vector<mypoint2f> project_p;
			for (unsigned int j = 0; j < leftpoints.size(); ++j)
			{
				double mini_dis = 6553500;
				mypoint2f pro_p(0, 0);
				for (unsigned int k = 0; k < mid_plist.size(); ++k)
				{
					double tmp_dis = LineSegmentLength(mid_plist[k], sgraph.nodeList[leftpoints[j]].position);
					if (tmp_dis< mini_dis)
					{
						mini_dis = tmp_dis;
						pro_p = mid_plist[k];
					}
				}
				project_p.push_back(pro_p);
			}
			//根据投影点， 对leftpoint排序: p1--> leftpoint -->p2,存储到leftpoint_order
			vector<int> leftpoint_order;
			vector<mypoint2f> project_p_order;
			if (leftpoints.size() == 1)
			{
				leftpoint_order.push_back(p1); leftpoint_order.push_back(leftpoints[0]); leftpoint_order.push_back(p2);
				project_p_order.push_back(sgraph.nodeList[p1].position);
				project_p_order.push_back(project_p[0]);
				project_p_order.push_back(sgraph.nodeList[p2].position);
			}
			else
			{
				leftpoint_order.push_back(p1);
				project_p_order.push_back(sgraph.nodeList[p1].position);
				vector<pair<int, double>> ptp_dis;// ptp_dis.push_back(0);
				mypoint2f p1_position = sgraph.nodeList[p1].position;
				mypoint2f p2_position = sgraph.nodeList[p2].position;
				for (unsigned int j = 0; j < project_p.size(); ++j)
				{
					double tmp = LineSegmentLength(p1_position, project_p[j]);
					bool inflag = false;
					vector<pair<int, double>>::iterator iter_op = ptp_dis.begin();//int：patch的编号，double:patch的周长
					while (iter_op != ptp_dis.end())
					{
						if (tmp < iter_op->second)//从xiao到da排序
						{
							ptp_dis.insert(iter_op, make_pair(leftpoints[j], tmp));
							inflag = true;
							break;
						}
						iter_op++;
					}
					if (inflag == false)
						ptp_dis.push_back(make_pair(leftpoints[j], tmp));
				}
				for (unsigned int j = 0; j < ptp_dis.size(); ++j)
				{
					leftpoint_order.push_back(ptp_dis[j].first);
					for (unsigned int k = 0; k < leftpoints.size(); ++k)
					{
						if (leftpoints[k] == ptp_dis[j].first)
						{
							project_p_order.push_back(project_p[k]);
							break;
						}
					}
				}
				leftpoint_order.push_back(p2);
				project_p_order.push_back(sgraph.nodeList[p2].position);
			}
			//一个一个的处理（改变node顶点坐标，替换polyline等）
			//sgraph.nodeList[leftpoints[j]].position 和 sgraph.nodeList[anode->adjVertex].position 在原来的patch上，先不做变换，在for之后变换。若不在原patch上例如outline或另一个patch上度为2的点则改变坐标？
			vector<int> usednode;//该变量，存储countedge=2的点，表示已经用过，用在specialTransformEdge()函数里
			for (unsigned int j = 0; j < leftpoints.size(); ++j)
			{
				int alp = leftpoints[j];
				ArcNode* anode = sgraph.nodeList[alp].firstEdge;
				mypoint2f trans_p = project_p[j] - sgraph.nodeList[leftpoints[j]].position;//增量
				while (anode)
				{
					if (anode->lineFlag == 2)//condition1:无patch相连的线，即outline ，整体向内平移
					{
						translateOutline(leftpoints[j], anode, trans_p);
					}
					else if (anode->lineFlag == 1)
					{
						//还没有考虑这种情况
						cout << endl;
					}
					else if (!anode->patchesIndex.empty())
					{
						if (isInVector(tp_vec, anode->adjVertex))//是不是在原patch周边上， tp_vec = plrGraph.patches[i].pointIndexOfPatch;在上面定义过
						{
							//贴着原来patch的edge，截取midplist的一段当做新的edge
							//截取+替换
							vector<mypoint2f> cutout;//pstart pend是midplist中的坐标位置
							mypoint2f pstart = project_p[j];//sgraph.nodeList[leftpoints[j]].position;
							mypoint2f pend = sgraph.nodeList[anode->adjVertex].position;
							if (isInVector(leftpoints, anode->adjVertex))
							{
								for (unsigned int k = 0; k < project_p.size(); ++k)
								{
									if (leftpoints[k] == anode->adjVertex)
									{
										pend = project_p[k];
										break;
									}
								}
							}
							truncateAnEdge(mid_plist, &cutout, pstart, pend);
							substituteEdge(anode->curveType, anode->curveIndex, cutout);
						}
						else
						{
							//不在原patch上，且是另一个patch上的一点。
							//1)若anode->adjVertex的countEdgeNum(，，0)!=2  发生形变（膨胀？）增量*百分比+原坐标=new edge
							//2)countedge()=2，找到结束点，判断结束点是否在原patch上
							int ecount = countEdgeNum(anode->adjVertex, 0);
							if (ecount != 2)
							{
								mypoint2f pstart = sgraph.nodeList[leftpoints[j]].position;//增量100%
								mypoint2f pend = sgraph.nodeList[anode->adjVertex].position;//增量0%
								vector<mypoint2f> ori_plist;
								vector<mypoint2f> new_plist;
								getPListOfaCurve(pstart, pend, leftpoints[j], anode->adjVertex, &ori_plist,0);
								transformEdge(&new_plist, trans_p, &ori_plist);  //增量在while循环外最上面定义了	
								//sgraph.nodeList[leftpoints[j]].position = new_plist.front();//shan
								//sgraph.nodeList[anode->adjVertex].position = new_plist.back();//不用改，坐标没有发生变化
								substituteEdge(anode->curveType, anode->curveIndex, new_plist);
							}
							else
							{
								//已知leftpoints[j]和anode->ver, 再往下寻找，一直找到度！=2的结束点。然后判断结束点是在原patch上海市其他patch上					
								specialTransformEdge(leftpoints[j], anode, pch, &leftpoints, &project_p, &usednode);//i式原来被删掉的patch的index									
							}
						}
					}
					anode = anode->next;
				}
			}
			for (unsigned int j = 0; j < leftpoints.size(); ++j)
			{
				sgraph.nodeList[leftpoints[j]].position = project_p[j];
			}
			//根据leftpoint_order & project_p_order对相关点进行连接，有的连接之后变成了patch的一部分有的变成了outline，需要添加到plrgraph.outline中
			vector<pair<int, int>> add_arcnode;
			for (unsigned int j = 1; j < leftpoint_order.size(); ++j)
			{
				//先连接0-1，1-2，2-3，，，
				int tp1 = leftpoint_order[j - 1];
				int tp2 = leftpoint_order[j];
				ArcNode* anode = getArcNode(tp1, tp2);
				if (anode == NULL)
				{
					ArcNode* newnode = new ArcNode;
					newnode->adjVertex = tp2;
					newnode->curveType = 'l';
					vector<mypoint2f> cutout;
					mypoint2f pstart = project_p_order[j - 1];  //tp1的投影点
					mypoint2f pend = project_p_order[j];  //tp2的投影点
					truncateAnEdge(mid_plist, &cutout, pstart, pend);
					PolyLine apl; apl.pointlist = cutout;
					polyline_vec.push_back(apl);
					newnode->curveIndex = polyline_vec.size() - 1;
					newnode->lineFlag = 0;
					newnode->next = sgraph.nodeList[tp1].firstEdge;
					sgraph.nodeList[tp1].firstEdge = newnode;

					ArcNode* newnode2 = new ArcNode;
					newnode2->adjVertex = tp1;
					newnode2->curveType = 'l';
					newnode2->curveIndex = polyline_vec.size() - 1;
					newnode2->lineFlag = 0;
					newnode2->next = sgraph.nodeList[tp2].firstEdge;
					sgraph.nodeList[tp2].firstEdge = newnode2;
					sgraph.nodeList[tp2].orderedArcs;
				}
				//leftpoint_order[j]->leftpoint_order[j-2/j-3.../0]向后排查有没有 应该断开的edge，
				//从这条长arc得到对应的patchindex，去修改patch->plist和patch->edge_vec。同时修改 新加入的arcnode的patchindex
				int k = j - 2;
				while (k >= 0)
				{
					ArcNode* enode = getArcNode(leftpoint_order[j], leftpoint_order[k]);
					if (enode != NULL)
					{
						int pchindex = enode->patchesIndex[0];//或者for循环，找到patchindex不是i的那个
						vector<int> pch_verlist = plrGraph.patches[pchindex].pointIndexOfPatch;
						vector<int> origi_piece;
						origi_piece.push_back(leftpoint_order[j]);
						origi_piece.push_back(leftpoint_order[k]);
						vector<int> target_piece;
						int g = j;
						while (g != k)
						{
							target_piece.push_back(leftpoint_order[g]);
							g--;
						}
						target_piece.push_back(leftpoint_order[k]);
						bool tar_flag = false;
						for (unsigned int ti = 1; ti < target_piece.size() - 1; ++ti)
						{
							if (isInVector(pch_verlist, target_piece[ti]))//需要添加的点在patch中本来就有，则不做修改
							{
								tar_flag = true;
								break;
							}
						}
						if (tar_flag == false)
						{
							Patch *modify_pch = getAPatch(pchindex); //modify_pch->pointIndexOfPatch; modify_pch->edges_vec;
							modifyPatch_plist(modify_pch, origi_piece, target_piece);
							removeArcNode(leftpoint_order[j], leftpoint_order[k], 0);//删除已经拆分掉的arc
							removeArcNode(leftpoint_order[k], leftpoint_order[j], 0);
						}
					}
					k--;
				}
				//新添加的arc（tp1,tp2）,还不确定是outline还是patch的一部分，先记录下来,for循环外进行判断
				if (anode == NULL)
				{
					add_arcnode.push_back(make_pair(tp1, tp2));
				}
			}
			if (!add_arcnode.empty())
			{
				for (unsigned int ai = 0; ai < add_arcnode.size(); ++ai)
				{
					int n1 = add_arcnode[ai].first;
					int n2 = add_arcnode[ai].second;
					ArcNode* nnarc = getArcNode(n1, n2); //new arcnode 的innerlin = 0且patchindex.empty()说明是一条新的线
					if (nnarc->lineFlag == 0 && nnarc->patchesIndex.empty())
					{
						nnarc->lineFlag = 2;
						nnarc = getArcNode(n2, n1);//反向设置一遍
						nnarc->lineFlag = 2;
						vector<int> tmp_vec;
						tmp_vec.push_back(n1);
						tmp_vec.push_back(n2);
						ICurve icv;
						icv.type = 0;
						icv.pt_vec = tmp_vec;
						plrGraph.attachlines.push_back(icv);
					}
				}
			}
		}

	}
}//单个操作 combine(merge) erosion(变成一条线)
bool MyGraph::checkErosion(int pch, int p1, int p2, int& anotherpch)
{
	bool flag = false;
	vector<int> pch_index = plrGraph.patches[pch].pointIndexOfPatch;
	pch_index.pop_back();
	EdgeOfPatch *edgep1 = NULL;
	EdgeOfPatch *edgep2 = NULL;
	//先判断p1点
	if (pch_index.front() == p1)
	{
		edgep1 = getBorderofPatch(pch, p1, pch_index.back());
		edgep2 = getBorderofPatch(pch, p1, pch_index[1]);
	}
	else if (pch_index.back() == p1)
	{
		edgep1 = getBorderofPatch(pch, p1, pch_index.front());
		edgep2 = getBorderofPatch(pch, p1, pch_index[pch_index.size()-2]);
	}
	else
	{
		for (unsigned int i = 0; i < pch_index.size(); ++i)
		{
			if (pch_index[i] == p1)
			{
				edgep1 = getBorderofPatch(pch, p1, pch_index[i-1]);
				edgep2 = getBorderofPatch(pch, p1, pch_index[i+1]);
				break;
			}
		}
	}
	if (edgep1->adjPatchIndex.size() > 1 && edgep2->adjPatchIndex.size() > 1)
	{
		if (edgep1->adjPatchIndex[1] != edgep2->adjPatchIndex[1])
			flag = true;
		else
		{
			flag = false;
			anotherpch = edgep1->adjPatchIndex[1];
		}
	}
	else
	{
		flag = true;
	}
	//再判断p2点
	if (flag == true)
	{
		if (pch_index.front() == p2)
		{
			edgep1 = getBorderofPatch(pch, p2, pch_index.back());
			edgep2 = getBorderofPatch(pch, p2, pch_index[1]);
		}
		else if (pch_index.back() == p2)
		{
			edgep1 = getBorderofPatch(pch, p2, pch_index.front());
			edgep2 = getBorderofPatch(pch, p2, pch_index[pch_index.size() - 2]);
		}
		else
		{
			for (unsigned int i = 0; i < pch_index.size(); ++i)
			{
				if (pch_index[i] == p2)
				{
					edgep1 = getBorderofPatch(pch, p2, pch_index[i - 1]);
					edgep2 = getBorderofPatch(pch, p2, pch_index[i + 1]);
					break;
				}
			}
		}
		if (edgep1->adjPatchIndex.size() > 1 && edgep2->adjPatchIndex.size() > 1)
		{
			if (edgep1->adjPatchIndex[1] != edgep2->adjPatchIndex[1])
				flag = true;
			else
			{
				flag = false;
				anotherpch = edgep1->adjPatchIndex[1];
			}
		}
		else
		{
			flag = true;
		}
	}
	return flag;
}

/* 线 ，method 1/2执行完了，再处理线 */
void MyGraph::attachLineProcess(double width_threshold, double angle_threshold)
{
	//vector<int> patch_points = plrGraph.patches[pchindex].pointIndexOfPatch;
	//for (unsigned int i = 0; i < patch_points.size(); ++i)
	//{
	//	bool attachflag = false;
	//	ArcNode* arc=sgraph.nodeList[patch_points[i]].firstEdge;
	//	while (arc != NULL)
	//	{
	//		if (arc->lineFlag != 0)
	//		{
	//			attachflag = true;
	//			break;
	//		}
	//		arc = arc->next;
	//	}
	//	if (attachflag == true)//patch_points[i]该点有悬着的线
	//	{
	//	}
	//}

	if (width_threshold == 0)
		width_threshold = 10;
	cout << "attachLineProcess" << endl;
	cout << "width_threshold=" << width_threshold << endl;
	//在divideAttachLines()里，otherline= attach + isolated （isolatedlines和 outline 先放到一块处理）
	int linesize = plrGraph.otherlines.size();  //outline和isolatedline 存到了一起
	vector<double> line_weight;//每条线都有一个weight，weight=proximity+continuity
	vector<int> minwei_index;
	vector<double> proxi_vec;//输出用
	vector<double> conti_vec;//输出用
	vector<vector<mypoint2f>> record_plist;
	for (int i = 0; i < linesize; ++i)
	{
		line_weight.push_back(6553500);
		minwei_index.push_back(-1);
		proxi_vec.push_back(6553500);//输出用
		conti_vec.push_back(6553500);//输出用
		vector<mypoint2f> plist;
		if (plrGraph.otherlines[i].type >= 0)
		{
			vector<int> cv_pindex = plrGraph.otherlines[i].pt_vec;
			for (unsigned int j = 1; j < cv_pindex.size(); ++j)
			{
				int verindex1 = cv_pindex.at(j - 1);
				int verindex2 = cv_pindex.at(j);
				mypoint2f point1 = sgraph.nodeList[verindex1].position;//=lineWithCommonPoints[verindex1].first
				mypoint2f point2 = sgraph.nodeList[verindex2].position;
				if (!plist.empty())
					plist.pop_back();
				getPListOfaCurve(point1, point2, verindex1, verindex2, &plist, 0);
			}
		}		
		record_plist.push_back(plist);
	}
	//计算每条curve的最近邻neighborcurve
	for (int i = 0; i < linesize - 1; ++i)//注意linesize-1
	{
		vector<mypoint2f> plist1 = record_plist[i];
		if (plist1.empty())
			continue;
		for (int j = i + 1; j < linesize; ++j)
		{
			vector<mypoint2f> plist2 = record_plist[j];
			if (plist2.empty())
				continue;
			int n_pindex1 = -1; int n_pindex2=-1;
			double proxi = getProximity(&plist1, &plist2, n_pindex1, n_pindex2);//找两条curve距离最近的点，并返回他们的位置pindex
			double conti = getContinuity(&plist1, &plist2, n_pindex1, n_pindex2);//最近距离处的夹角
			double w_sum = proxi;  // +conti;
			if (w_sum < line_weight[i])
			{
				if (conti < angle_threshold)//若不满足这个条件，可能是交叉但垂直
				{
					minwei_index[i] = j;
					line_weight[i] = w_sum;
				}	
			}
			if (w_sum < line_weight[j])
			{
				if (conti < angle_threshold)
				{
					minwei_index[j] = i;
					line_weight[j] = w_sum;
				}			
			}
			if (proxi < proxi_vec[i])//输出用
				proxi_vec[i] = proxi;
			if (conti < conti_vec[i])
				conti_vec[i] = conti;
			if (proxi < proxi_vec[j])
				proxi_vec[j] = proxi;
			if (conti < conti_vec[j])
				conti_vec[j] = conti;
		}
	}
	string file_str = "E:/proxi_weight.txt";//	file_str.append(to_string(iteration));
	ofstream outfile1(file_str);//"E:/record_similarity.txt"
	outfile1.is_open();
	ofstream outfile2("E:/conti_weight.txt");
	outfile2.is_open();
	for (int i = 0; i < linesize; ++i)//注意linesize-1
	{
		outfile1 << proxi_vec[i] << ' ';
		outfile2 << conti_vec[i] << ' ';
	}
	outfile1.close();
	outfile2.close();
	//heap
	//建堆，line_weight保存了从第0个点到底n个点weight值，堆heap_vec保存每个patch的权值在tt_w中的下标，也是patch 的index
	vector<int> heap_vec;
	vector<int> ele_index;
	for (unsigned int i = 0; i < line_weight.size(); ++i)
	{
		ele_index.push_back(-1);
		if (plrGraph.otherlines[i].type >= 0)
		{
			heap_vec.push_back(i);
		}
		//heap_vec.push_back(i);
		//ele_index.push_back(0);
	}
	//建堆
	int h_size = heap_vec.size();
	for (int i = h_size / 2; i >= 1; i--)
	{
		heap_siftdown(i, &heap_vec, &line_weight, &ele_index);
	}
	//建立索引，知道某个权值的下标，可以直接找到在heap的位置，位置是从1开始，注意更新
	for (int i = 0; i < h_size; ++i)
	{
		ele_index[heap_vec[i]] = i;// ele_index
	}
	
	//int stop_num = final_pnum;//h_size / 5;
	double t_sum = width_threshold;// +angle_threshold;
	while (h_size > 0)//while (!heap_vec.empty())//
	{
		int pindex = heap_vec[0];
		if (line_weight[pindex] > t_sum)//剩下的weight都大于t_sum
		{
			cout << "left line number:" << h_size << endl;
			break;
		}
		else if (line_weight[pindex] == 6553500)
		{
			cout << "left line number:" << h_size << endl;
			break;//若final_pnum设置的过小，使weight=65535，说明剩下的patch都没有neighbor,可以直接跳出循环，此时最小weight对应的patchindex， weight_pindex[pindex]=-1
		}
		int combine_line = minwei_index[pindex];
		outfile1 << line_weight[pindex];
		outfile1 << " ";
		vector<mypoint2f> new_plist;
		combineTwoCurves(combine_line, pindex, &new_plist);//process!!//删除pindex,保留combine_line

		//1.先删除堆顶的pindex
		int tmp_back = heap_vec.back();
		int tmp_mini = heap_vec[0];
		line_weight[tmp_mini] = -1;
		minwei_index[tmp_mini] = -1;
		ele_index[tmp_mini] = -1;
		ele_index[tmp_back] = 0;//!!!
		heap_vec[0] = tmp_back;
		heap_siftdown(1, &heap_vec, &line_weight, &ele_index);//下沉调整
		heap_vec.pop_back();
		h_size--;
		record_plist[pindex].swap(vector<mypoint2f>());
		//2.更新newCurve即combinr_line的weight,及对应的lineindex
		record_plist[combine_line].swap(vector<mypoint2f>());
		record_plist[combine_line] = new_plist;
		int n_lindex = -1;
		double n_weight = 6553500;
		for (int i = 0; i < linesize; ++i)
		{
			vector<mypoint2f> cur_plist = record_plist[i];
			if (cur_plist.empty() || i == combine_line)
				continue;
			int n_pindex1, n_pindex2;
			double proxi = getProximity(&new_plist, &cur_plist, n_pindex1, n_pindex2);//找两条curve距离最近的点，并返回他们的位置pindex
			double conti = getContinuity(&new_plist, &cur_plist, n_pindex1, n_pindex2);//最近距离处的夹角
			double w_sum = proxi;// +conti;
			if (w_sum < n_weight)
			{
				if (conti < angle_threshold)//若不满足这个条件，可能是交叉但垂直
				{
					n_lindex = i;
					n_weight = w_sum;
				}
			}
		}
		heap_update(ele_index[combine_line] + 1, n_weight, &heap_vec, &line_weight, &ele_index);
		minwei_index[combine_line] = n_lindex;
		//改变neighbor是pindex combine_line的曲线的新neighbor，
		for (int i = 0; i < linesize; ++i)
		{
			if (minwei_index[i] == pindex || minwei_index[i] == combine_line)
			{
				int chageindex = i;
				vector<mypoint2f> change_plist = record_plist[i];
				int mini_lindex = -1;
				double mini_weight = 6553500;
				for (int j = 0; j < linesize; ++j)
				{
					vector<mypoint2f> cur_plist = record_plist[j];
					if (cur_plist.empty() || j == chageindex)
						continue;
					int n_pindex1, n_pindex2;
					double proxi = getProximity(&change_plist, &cur_plist, n_pindex1, n_pindex2);//找两条curve距离最近的点，并返回他们的位置pindex
					double conti = getContinuity(&change_plist, &cur_plist, n_pindex1, n_pindex2);//最近距离处的夹角
					double w_sum = proxi;// +conti;
					if (w_sum < mini_weight)
					{
						if (conti < angle_threshold)//若不满足这个条件，可能是交叉但垂直
						{
							mini_lindex = j;
							mini_weight = w_sum;
						}
					}
				}
				heap_update(ele_index[chageindex] + 1, mini_weight, &heap_vec, &line_weight, &ele_index);
				minwei_index[chageindex] = mini_lindex;
			}
		}
	
	}
}
void MyGraph::combineTwoCurves(int lindex1, int lindex2, vector<mypoint2f>* newcurve)//lindex1  lindex2是plrgraph.otherline中的下标！！！
{
	//保留lindex1 ，删除lindex2,
	vector<int> curve1 = plrGraph.otherlines[lindex1].pt_vec;
	vector<mypoint2f> plist1;
	jointMethod(&curve1, &plist1, 1);///////////////
	vector<int> curve2 = plrGraph.otherlines[lindex2].pt_vec;
	vector<mypoint2f> plist2;
	jointMethod(&curve2, &plist2, 1);/////////////////

	vector<mypoint2f> total_plist;
	total_plist.insert(total_plist.end(), plist1.begin(), plist1.end());
	total_plist.insert(total_plist.end(), plist2.begin(), plist2.end());
	
	cubicBezier bz;
	int fit_flag = 1;
	fitCubBezier(&total_plist, newcurve, &bz, fit_flag);//拟合新的Beziercurve
	if (_isnan(bz.p1.x) || !_finite(bz.p1.x))
	{
		newcurve->swap(vector<mypoint2f>());
		return;
	}
	//删除curve2 ，
	plrGraph.otherlines[lindex2].type = -1;
	vector<int> rm_line = plrGraph.otherlines[lindex2].pt_vec;
	for (size_t ti = 1; ti < rm_line.size(); ++ti)
	{
		removeArcNode(rm_line[ti - 1], rm_line[ti], 0);
		removeArcNode(rm_line[ti], rm_line[ti - 1], 0);
	}//removeArcNode
	//新增newcurve
	PolyLine pl;
	pl.origi_index = polyline_vec.size();
	pl.origi_type = 'L';
	pl.pointlist = *newcurve;
	polyline_vec.push_back(pl);
	//若lindex1/lindex2 有连接的patch，把patch 中的point_line：删除掉second中带有lindex2的,先计算所有的adjpatch。。nei_pch
	vector<int> terminal_p;
	terminal_p.push_back(curve1.front());
	terminal_p.push_back(curve1.back());
	terminal_p.push_back(curve2.front());
	terminal_p.push_back(curve2.back());
	vector<int> nei_pch;
	for (int i = 0; i < 4; ++i)
	{
		ArcNode* anode = sgraph.nodeList[terminal_p[i]].firstEdge;
		while (anode != NULL)
		{
			if (!anode->patchesIndex.empty())
			{
				nei_pch.insert(nei_pch.end(), anode->patchesIndex.begin(), anode->patchesIndex.end());
			}
			anode = anode->next;
		}
		anode = sgraph.nodeList[terminal_p[i]].firstEdge;
		while (anode != NULL)
		{
			if (!anode->patchesIndex.empty())
			{
				nei_pch.insert(nei_pch.end(), anode->patchesIndex.begin(), anode->patchesIndex.end());
			}
			anode = anode->next;
		}
	}
	//更改sgraph相关
	if (fit_flag == 0)
	{
		//添加点
		VertexNode nnode1;
		VertexNode nnode2;
		nnode1.position = newcurve->front();
		nnode1.verIndex = sgraph.nodeList.size();
		nnode1.firstEdge = NULL;
		sgraph.nodeList.push_back(nnode1);

		nnode2.verIndex = sgraph.nodeList.size();
		nnode2.position = newcurve->back();
		nnode2.firstEdge = NULL;
		sgraph.nodeList.push_back(nnode2);

		ArcNode* arc1 = new ArcNode;
		arc1->adjVertex = nnode2.verIndex;
		arc1->curveType = 'L';
		arc1->curveIndex = polyline_vec.size() - 1;
		arc1->lineFlag = 2;
		arc1->next = sgraph.nodeList[nnode1.verIndex].firstEdge;
		sgraph.nodeList[nnode1.verIndex].firstEdge = arc1;
		ArcNode* arc2 = new ArcNode;
		arc2->adjVertex = nnode1.verIndex;
		arc2->curveType = 'L';
		arc2->curveIndex = polyline_vec.size() - 1;
		arc2->lineFlag = 2;
		arc2->next = sgraph.nodeList[nnode2.verIndex].firstEdge;
		sgraph.nodeList[nnode2.verIndex].firstEdge = arc2;
		////判断new curve的前后端点 是否与原来的两个curve的前后端点相同，若相同则不用向sgraph中添加新点
		//int ssflag_1 = -1; int ssflag_2 = -1;
		//vector<mypoint2f> four_cor;
		//four_cor.push_back(plist1.front()); four_cor.push_back(plist1.back()); four_cor.push_back(plist2.front()); four_cor.push_back(plist2.back());
		//for (int i = 0; i < 4; ++i)
		//{
		//	if (four_cor[i] == newcurve->front())
		//	{
		//		ssflag_1 = i;
		//	}
		//	else if (four_cor[i] == newcurve->back())
		//	{
		//		ssflag_2 = i;
		//	}
		//}
		//if (ssflag_1 < 0)
		//{
		//	if (ssflag_2 < 0)
		//	{
		//		nnode1.position = newcurve->front();
		//		nnode1.verIndex = sgraph.nodeList.size();
		//		nnode1.firstEdge = NULL;
		//		sgraph.nodeList.push_back(nnode1);
		//		
		//		nnode2.verIndex = sgraph.nodeList.size();
		//		nnode2.position = newcurve->back();
		//		nnode2.firstEdge = NULL;
		//		sgraph.nodeList.push_back(nnode2);
		//		ArcNode* arc1 = new ArcNode;
		//		arc1->adjVertex = nnode2.verIndex;
		//		arc1->curveType = 'L';
		//		arc1->curveIndex = polyline_vec.size() - 1;
		//		arc1->lineFlag = 2;
		//		arc1->next = sgraph.nodeList[nnode1.verIndex].firstEdge;
		//		sgraph.nodeList[nnode1.verIndex].firstEdge = arc1;
		//		ArcNode* arc2 = new ArcNode;
		//		arc2->adjVertex = nnode1.verIndex;
		//		arc2->curveType = 'L';
		//		arc2->curveIndex = polyline_vec.size() - 1;
		//		arc2->lineFlag = 2;
		//		arc2->next = sgraph.nodeList[nnode2.verIndex].firstEdge;
		//		sgraph.nodeList[nnode2.verIndex].firstEdge = arc2;
		//	}
		//	else
		//	{
		//		//某点与 newcurve->back()相等
		//	}
		//}
		//else
		//{
		//	if (ssflag_2 < 0)
		//	{
		//		//某点与 newcurve->front()相等
		//	}
		//	else
		//	{
		//		//new curve 的两个点都存在
		//	}
		//}

		//并且变更curve1。没有对sgraph中的点
		plrGraph.otherlines[lindex1].pt_vec.swap(vector<int>());
		vector<int> ncurve;
		ncurve.push_back(nnode1.verIndex);
		ncurve.push_back(nnode2.verIndex);
		plrGraph.otherlines[lindex1].pt_vec = ncurve;
		//若lindex1/lindex2 有连接的patch，把patch 中的point_line：删除掉second中带有lindex2的
		if (!nei_pch.empty())
		{
			deleRepeatPoint(&nei_pch);
			for (unsigned int i = 0; i < nei_pch.size(); ++i)
			{
				int tmp_pindex = nei_pch[i];
				vector<pair<int, int>>::iterator iter_pl = plrGraph.patches[tmp_pindex].point_line.begin();
				while (iter_pl != plrGraph.patches[tmp_pindex].point_line.end())
				{
					if (iter_pl->second == lindex2)
					{
						iter_pl = plrGraph.patches[tmp_pindex].point_line.erase(iter_pl);
					}
					else
						iter_pl++;
				}
			}
		}
	}
	else
	{
		//找到newBeziercurve 两端点对应的index
		int p0_index = 0;
		int p3_index = 0;
		if (bz.p0 == plist1.front() || bz.p0 == plist1.back())
		{
			if (bz.p0 == plist1.front())
			{
				p0_index = curve1[0];
			}
			else
			{
				p0_index = curve1.back();
			}
		}
		else if (bz.p0 == plist2.front() || bz.p0 == plist2.back())
		{
			if (bz.p0 == plist2.front())
			{
				p0_index = curve2[0];
			}
			else
			{
				p0_index = curve2.back();
			}
		}
		if (bz.p3 == plist1.front() || bz.p3 == plist1.back())
		{
			if (bz.p3 == plist1.front())
			{
				p3_index = curve1[0];
			}
			else
			{
				p3_index = curve1.back();
			}
		}
		else if (bz.p3 == plist2.front() || bz.p3 == plist2.back())
		{
			if (bz.p3 == plist2.front())
			{
				p3_index = curve2[0];
			}
			else
			{
				p3_index = curve2.back();
			}
		}
		//修改arc
		ArcNode* arc1 = getArcNode(p0_index, p3_index);
		if (arc1 != NULL)
		{
			arc1->curveType = 'L';
			arc1->curveIndex = polyline_vec.size() - 1;
			ArcNode* arc1_rever = getArcNode(p3_index, p0_index);
			arc1_rever->curveType = 'L';
			arc1_rever->curveIndex = polyline_vec.size() - 1;
		}
		else
		{
			ArcNode* arc1 = new ArcNode;
			arc1->adjVertex = p3_index;
			arc1->curveType = 'L';
			arc1->curveIndex = polyline_vec.size() - 1;
			arc1->lineFlag = 2;
			arc1->next = sgraph.nodeList[p0_index].firstEdge;
			sgraph.nodeList[p0_index].firstEdge = arc1;
			ArcNode* arc2 = new ArcNode;
			arc2->adjVertex = p0_index;
			arc2->curveType = 'L';
			arc2->curveIndex = polyline_vec.size() - 1;
			arc2->lineFlag = 2;
			arc2->next = sgraph.nodeList[p3_index].firstEdge;
			sgraph.nodeList[p3_index].firstEdge = arc2;
		}
		//变更curve1
		plrGraph.otherlines[lindex1].pt_vec.swap(vector<int>());
		vector<int> ncurve;
		ncurve.push_back(p0_index);
		ncurve.push_back(p3_index);
		plrGraph.otherlines[lindex1].pt_vec = ncurve;
		//若lindex1/lindex2 有连接的patch，把patch 中的point_line：删除掉second中带有lindex2的
		if (!nei_pch.empty())
		{
			deleRepeatPoint(&nei_pch);
			for (unsigned int i = 0; i < nei_pch.size(); ++i)
			{
				int tmp_pindex = nei_pch[i];
				vector<pair<int, int>>::iterator iter_pl = plrGraph.patches[tmp_pindex].point_line.begin();
				while (iter_pl != plrGraph.patches[tmp_pindex].point_line.end())
				{
					if (iter_pl->first == p0_index || iter_pl->first == p3_index)
					{
						if (iter_pl->second == lindex2)
						{
							iter_pl->second = lindex1;
						}
						iter_pl++;
					}
					else if (iter_pl->second == lindex2)
					{
						iter_pl = plrGraph.patches[tmp_pindex].point_line.erase(iter_pl);
					}
					else if (iter_pl->second == lindex1)
					{
						if (iter_pl->first != p0_index && iter_pl->first != p3_index)
						{
							iter_pl = plrGraph.patches[tmp_pindex].point_line.erase(iter_pl);
						}
					}
					else
						iter_pl++;
				}
			}
		}
	}
	
}
void MyGraph::fitCubBezier(vector<mypoint2f>* allplist, vector<mypoint2f>* cbcurve, cubicBezier* cubBezier, int flag)
{
	//flag=0表示不固定两端点，=1固定两个端点
	//删除重复的数据点
	vector<mypoint2f>::iterator iter_p = allplist->begin();
	while (iter_p != allplist->end())
	{
		vector<mypoint2f>::iterator iter_p_inner = iter_p + 1;
		while (iter_p_inner != allplist->end())
		{
			if (LineSegmentLength(*iter_p, *iter_p_inner) < eps)
			{
				iter_p_inner = allplist->erase(iter_p_inner);
			}
			else
				iter_p_inner++;
		}
		iter_p++;
	}
	//PCA 降维  直角坐标系
	//先计算协方差矩阵A ,初始化 数据点矩阵MatrixXd B（2，n） 数据按列排列！！
	int plist_size = allplist->size();
	//vec_x，vec_y，B 都是直角坐标系下的，
	vector<double> vec_x;
	vector<double> vec_y;
	MatrixXd B = MatrixXd::Zero(2, plist_size);
	for (int i = 0; i < plist_size; ++i)
	{
		vec_x.push_back(allplist->at(i).x);
		vec_y.push_back(allplist->at(i).y);
		B(0, i) = allplist->at(i).x;
		B(1, i) = allplist->at(i).y;
	}
	double aver_x = GetAverange(vec_x);
	double aver_y = GetAverange(vec_y);
	double var_cxx = GetVariance(vec_x, aver_x, vec_x, aver_x);
	double var_cyy = GetVariance(vec_y, aver_y, vec_y, aver_y);
	double var_cxy = GetVariance(vec_x, aver_x, vec_y, aver_y);
	Matrix2d A;
	A << var_cxx, var_cxy, var_cxy, var_cyy;//cout << "Here is a covariance 2*2 matrix, A:" << endl << A << endl << endl;
	//计算协方差矩阵的特征值 特征向量
	EigenSolver<MatrixXd> es(A);
	//cout << "The eigenvalues of A are:" << endl << es.eigenvalues() << endl;//注意特征向量和特征值都是 复数类型！！complex
	//cout << "The eigenvectors is:" << endl << es.eigenvectors() << endl;
	//选择较大的那个特征值和其对应的特征向量，作为PCA的主方向
	double max_evalue = 0;
	VectorXd max_vec;
	if (es.eigenvalues()[0].real()> es.eigenvalues()[1].real())
	{
		max_evalue = es.eigenvalues()[0].real();
		max_vec = es.eigenvectors().col(0).real();
	}
	else
	{
		max_evalue = es.eigenvalues()[1].real();
		max_vec = es.eigenvectors().col(1).real();
	}//cout << "max_evalue=" << max_evalue << endl;//cout << "max_v=" << max_vec << endl;
	//projection_m是所有数据点 投影到主特征向量上的结果（投影点），二维降到一维
	VectorXd projection_m = max_vec.transpose()*B;
	//把vectorXd 类型转换成 vector<double>, 并sort()从小到大排序，得到有序点
	vector<pair<int, double>> proMatrix_order; //int对应元数据点在total_plist中的下标，double 是projection_m中的值
	for (int i = 0; i < plist_size; ++i)
	{
		proMatrix_order.push_back(make_pair(i, projection_m[i]));
	}
	sort(proMatrix_order.begin(), proMatrix_order.end(), mycompare);
	//最小二乘法 拟合cubic Bezier curve ，变量的命名参考自http://jimherold.com/2012/04/20/least-squares-bezier-fit/
	//初始化有序坐标点vec_x_ordered , vec_y_ordered
	VectorXd vec_x_ordered = VectorXd::Zero(plist_size, 1);//即x坐标值 列向量
	VectorXd vec_y_ordered = VectorXd::Zero(plist_size, 1);
	for (int i = 0; i < plist_size; ++i)
	{
		vec_x_ordered[i] = allplist->at(proMatrix_order[i].first).x;
		vec_y_ordered[i] = allplist->at(proMatrix_order[i].first).y;
		//ordered_plist.push_back(mypoint2f(vec_x_ordered[i], vec_y_ordered[i]));
	}
	//初始化参数矩阵T，弦长参数化？proMatrix_order
	double startp = proMatrix_order.front().second;
	double endp = proMatrix_order.back().second;
	double arclen = endp - startp;
	vector<double> para_t;
	for (int i = 0; i < plist_size; ++i)
	{
		para_t.push_back((proMatrix_order[i].second - startp) / arclen);
	}

	if (flag == 0)
	{
		//初始化矩阵M, 赋值符号<<是一行一行的赋值 ，变量的命名参考自http://jimherold.com/2012/04/20/least-squares-bezier-fit/
		Matrix4d cub_m = Matrix4d::Zero(4, 4);
		cub_m << -1, 3, -3, 1, 3, -6, 3, 0, -3, 3, 0, 0, 1, 0, 0, 0;
		Matrix4d cub_m_inverse = cub_m.inverse();
		MatrixXd T = MatrixXd::Zero(plist_size, 4);
		for (int i = 0; i < plist_size; ++i)
		{
			double tmp_t = para_t[i];
			T(i, 0) = pow(tmp_t, 3);
			T(i, 1) = pow(tmp_t, 2);
			T(i, 2) = tmp_t;
			T(i, 3) = 1;
		}
		MatrixXd T_tran = T.transpose();
		MatrixXd TT_inverse = (T_tran*T).inverse();
		//cout << cub_m_inverse.rows() << "," << cub_m_inverse.cols() << endl;
		//cout << TT_inverse.rows() << "," << TT_inverse.cols() << endl;
		//cout << T_tran.rows() << "," << T_tran.cols() << endl;
		//cout << vec_y_ordered.rows() << "," << vec_y_ordered.cols() << endl;
		VectorXd cx = cub_m_inverse*TT_inverse*T_tran*vec_x_ordered;
		VectorXd cy = cub_m_inverse*TT_inverse*T_tran*vec_y_ordered;
		//cout << cx << endl;cout << cy << endl;
		cubBezier->p0 = mypoint2f(cx[0], cy[0]);
		cubBezier->p1 = mypoint2f(cx[1], cy[1]);
		cubBezier->p2 = mypoint2f(cx[2], cy[2]);
		cubBezier->p3 = mypoint2f(cx[3], cy[3]);
		cubicBezier_vec.push_back(*cubBezier);
		getPListfromCubic(0, *cubBezier, 0, cbcurve);
	}
	else
	{
		mypoint2f cb_p0 = allplist->at(proMatrix_order[0].first);
		mypoint2f cb_p3 = allplist->at(proMatrix_order.back().first);
		//需要用到参数t(para_t) ,p0p3(cb_p0,cb_p3), 有序的数据点(vec_x/y_ordered),计算A_1, A_2, A_12,C_1,C_2,再计算p1p2
		double A_1 = 0; double A_2 = 0;
		double A_12 = 0;
		//double C_1 = 0; double c_2 = 0;
		Vector2d tmp_c1(0, 0);
		Vector2d tmp_c2(0, 0);
		Vector2d cb_p0_vec(cb_p0.x, cb_p0.y);
		Vector2d cb_p3_vec(cb_p3.x, cb_p3.y);
		for (int i = 1; i < plist_size - 1; ++i)
		{
			double b_0 = pow((1 - para_t[i]), 3);
			double b_1 = 3 * pow((1 - para_t[i]), 2) * para_t[i];
			double b_2 = 3 * (1 - para_t[i]) * pow(para_t[i], 2);
			double b_3 = pow(para_t[i], 3);
			A_1 = A_1 + pow(b_1, 2);
			A_2 = A_2 + pow(b_2, 2);
			A_12 = A_12 + b_1*b_2;
			//C_1 = C_1 + b_1*(vec_x_ordered[i] - b_0*cb_p0.x - b_3*cb_p3.x);
			Vector2d tmp_pi(vec_x_ordered[i], vec_y_ordered[i]);
			Vector2d mid_vec = tmp_pi - b_0*cb_p0_vec - b_3*cb_p3_vec;
			tmp_c1 = tmp_c1 + b_1*mid_vec;
			tmp_c2 = tmp_c2 + b_2*mid_vec;
		}
		double fenmu = A_1*A_2 - pow(A_12, 2);
		Vector2d fenzi = A_2*tmp_c1 - A_12*tmp_c2;
		Vector2d final_c1 = fenzi / fenmu;
		fenzi = A_1*tmp_c2 - A_12*tmp_c1;
		Vector2d final_c2 = fenzi / fenmu;
		cubBezier->p0 = cb_p0;
		cubBezier->p1 = mypoint2f(final_c1[0], final_c1[1]);
		cubBezier->p2 = mypoint2f(final_c2[0], final_c2[1]);
		cubBezier->p3 = cb_p3;
		cubicBezier_vec.push_back(*cubBezier);
		getPListfromCubic(0, *cubBezier, 0, cbcurve);
	}
}
void MyGraph::removeCurve(vector<int>* rm_curve, int cindex_plr)
{
	//在sgraph中删除arc node
	for (unsigned int i = 1; i < rm_curve->size(); ++i)
	{
		ArcNode* tmparc = sgraph.nodeList[i-1].firstEdge;
		while (tmparc != NULL)
		{
			if (tmparc->adjVertex == rm_curve->at(i))
			{
				removeArcNode(rm_curve->at(i - 1), rm_curve->at(i), 0);
				removeArcNode(rm_curve->at(i), rm_curve->at(i - 1), 0);
				break;
			}
			tmparc = tmparc->next;
		}
	}
	//在plrgraph.attachlines（outline）中删除curve
	plrGraph.attachlines[cindex_plr].type = -1;
	//vector<ICurve>::iterator iter_l = plrGraph.outlines.begin();
	//iter_l = plrGraph.outlines.erase(iter_l + cindex_plr);  ///////cheeeeck
}
bool MyGraph::combineTwoCurves_getplist(int lindex1, int lindex2, vector<mypoint2f>* newcurve, double aver_w)
{
	bool sucflag = false; 
	if (plrGraph.otherlines[lindex2].type < 0)
	{
		return sucflag;
	}
	vector<int> curve1 = plrGraph.otherlines[lindex1].pt_vec;
	vector<mypoint2f> l1_plist;
	jointMethod(&curve1, &l1_plist, 0);
	vector<int> curve2 = plrGraph.otherlines[lindex2].pt_vec;
	vector<mypoint2f> l2_plist;
	jointMethod(&curve2, &l2_plist, 0);

	int n_pindex1 = -1; int n_pindex2 = -1;
	double proxi = getProximity(&l1_plist, &l2_plist, n_pindex1, n_pindex2);//找两条curve距离最近的点，并返回他们的位置pindex
	double conti = getContinuity(&l1_plist, &l2_plist, n_pindex1, n_pindex2);//最近距离处的夹角
	if (proxi < aver_w && conti < 20 && conti >= 0)  // 夹角小于20 ？？！！！！
	{
		sucflag = true;
		vector<mypoint2f> total_plist;
		total_plist.insert(total_plist.end(), l1_plist.begin(), l1_plist.end());
		total_plist.insert(total_plist.end(), l2_plist.begin(), l2_plist.end());
		cubicBezier bz;
		fitCubBezier(&total_plist, newcurve, &bz, 1);//拟合新的Beziercurve
		if (_isnan(bz.p1.x) || !_finite(bz.p1.x))
			sucflag = false;
		if (_isnan(bz.p2.x) || !_finite(bz.p2.x))
			sucflag = false;
	}
	return sucflag;	
}
void MyGraph::checkIntersect(int lindex, vector<int>* newlvec, vector<int>* neighbor_patch)
{
	int frontpt = newlvec->front();
	int backpt = newlvec->back();
	//bool cflag1 = true;
	//ArcNode* checkarc = sgraph.nodeList[frontpt].firstEdge;
	//while (checkarc != NULL)
	//{
	//	if (!checkarc->patchesIndex.empty())
	//	{
	//		cflag1 = false;
	//		break;
	//	}
	//	checkarc = checkarc->next;
	//}
	//bool cflag2 = true;
	//checkarc = sgraph.nodeList[backpt].firstEdge;
	//while (checkarc != NULL)
	//{
	//	if (!checkarc->patchesIndex.empty())
	//	{
	//		cflag2 = false;
	//		break;
	//	}
	//	checkarc = checkarc->next;
	//}

	//if (cflag1 == false || cflag2 = false)
	//	return;


	mypoint2f f_posi = sgraph.nodeList[frontpt].position;
	mypoint2f b_posi = sgraph.nodeList[backpt].position;
	double len = 2*LineSegmentLength(f_posi, b_posi);
	vector<int> nei_mixnode;//存 lindex 的neighbor mixednode，通过MixedTopol
	if (!neighbor_patch->empty())
		nei_mixnode = *neighbor_patch;

	////广度遍历  ,距离lindex的两个端点,len 范围内的patch
	//int mixedsize = mixedTopology.allNode.size();
	//int* visited = new int[mixedsize];
	//for (int i = 0; i < mixedsize; ++i)
	//{
	//	visited[i] = 0;
	//}
	//queue<int> st;
	//st.push(lindex);
	//visited[lindex] = 1;
	//while (!st.empty())
	//{
	//	int tmpnode = st.front();
	//	st.pop();
	//	visited[tmpnode] = 2;
	//	MixedArc* marc = mixedTopology.allNode[tmpnode].firstAdjNode;
	//	while (marc)
	//	{
	//		if (visited[marc->adj_index] == 0)
	//		{
	//			if (mixedTopology.allNode[marc->adj_index].status >= 0)
	//			{
	//				if (mixedTopology.allNode[marc->adj_index].pl_type == 'p')
	//				{
	//					double clen = LineSegmentLength(plrGraph.patches[marc->adj_index].centre_position, f_posi);
	//					double cle2 = LineSegmentLength(plrGraph.patches[marc->adj_index].centre_position, b_posi);
	//					if (clen <= len)
	//					{
	//						nei_mixnode.push_back(marc->adj_index);
	//						st.push(marc->adj_index);
	//						visited[marc->adj_index] = 1;
	//					}
	//					else if (cle2 <= len)
	//					{
	//						nei_mixnode.push_back(marc->adj_index);
	//						st.push(marc->adj_index);
	//						visited[marc->adj_index] = 1;
	//					}
	//				}
	//			}
	//		}
	//		marc = marc->next_arc;
	//	}
	//}
	//delete[] visited;
	//visited = NULL;
	////删除重复的patch
	//deleRepeatPoint(&nei_mixnode);
	//寻找patch与line最近的两个点
	vector<mypoint2f> l_plist;
	jointMethod(newlvec, &l_plist, 0);
	vector<vector<int>> patch_crossindex; //第一个值是patch 的index，其后的值是patch 的pointlist中相交的位置
	vector<int> line_crossindex;
	for (unsigned int i = 0; i < nei_mixnode.size(); ++i)
	{
		int pch = nei_mixnode[i];
		vector<mypoint2f> pch_plist;
		jointMethod(&plrGraph.patches[pch].pointIndexOfPatch, &pch_plist, 0);
		int size1 = pch_plist.size() - 2;
		for (int j = 1; j < size1; ++j)//pch_plist[j] 第j个点 ，与l_plist上的所有点进行比较，找到最近的那个点k
		{
			double mini_dis = 6553500;
			int mini_index = 0;
			int k;
			int size2 = l_plist.size() - 2;
			for (k = 1; k < size2; ++k)
			{
				double tmp = LineSegmentLength(pch_plist[j], l_plist[k]);
				if (tmp < mini_dis)
				{
					mini_dis = tmp;
					mini_index = k;
				}
			}
			bool crossflag = isSegmentIntersected(pch_plist[j - 1], pch_plist[j], l_plist[k - 1], l_plist[k]);
			bool crossflag2 = isSegmentIntersected(pch_plist[j - 1], pch_plist[j], l_plist[k], l_plist[k + 1]);
			bool crossflag3 = isSegmentIntersected(pch_plist[j + 1], pch_plist[j], l_plist[k - 1], l_plist[k]);
			bool crossflag4 = isSegmentIntersected(pch_plist[j + 1], pch_plist[j], l_plist[k], l_plist[k + 1]);
			vector<int> across_vec;
			across_vec.push_back(pch);
			if (crossflag == true)
			{
				across_vec.push_back(j);
				line_crossindex.push_back(k);
			}
			else if (crossflag2 == true)
			{
				across_vec.push_back(j);
				line_crossindex.push_back(k + 1);
			}
			else if (crossflag3 == true)
			{
				across_vec.push_back(j + 1);
				line_crossindex.push_back(k);
			}
			else if (crossflag4 == true)
			{
				across_vec.push_back(j + 1);
				line_crossindex.push_back(k + 1);
			}
			if (across_vec.size() > 1)
			{
				patch_crossindex.push_back(across_vec);
			}
		}
	}
	if (!patch_crossindex.empty())
	{
		vector<vector<int>>::iterator iter_out = patch_crossindex.begin();
		while (iter_out != patch_crossindex.end())
		{
			vector<vector<int>>::iterator iter_inner = iter_out + 1;
			while (iter_inner != patch_crossindex.end())
			{
				if (iter_out->at(0) == iter_inner->at(0) && iter_out->at(1) == iter_inner->at(1))
				{
					iter_inner = patch_crossindex.erase(iter_inner);
				}
				else 
					iter_inner++;
			}
			iter_out++;
		}
	}
	
	//只取line_crossindex两端的点，中间不对patch进行分割，
	//sgraph添加一个点，一个arc
	//patch只添加一个节点，吧一个edge分成两半，然后添加一个point_line
	if (!line_crossindex.empty())
	{
		deleRepeatPoint(&line_crossindex);
		if (line_crossindex.size() > 1)
		{
			sort(line_crossindex.begin(), line_crossindex.end());
		}

		vector<int> fbpt_ofline;//line与patches相交的第一个和最后一个点 其在l_plist中的位置。最多两个（从小到大）
		if (line_crossindex.front() != line_crossindex.back())
		{
			fbpt_ofline.push_back(line_crossindex.front());
			fbpt_ofline.push_back(line_crossindex.back());
		}
		else
		{
			fbpt_ofline.push_back(line_crossindex.front());
		}
		vector<pair<int, int>> pch_ptindex; //与fbpt_ofline相对应的patch的index，和 patch_plist中的相交的位置。最多两个
		for (unsigned int i = 0; i < fbpt_ofline.size(); ++i)
		{
			double mini_length = 6553500;
			int mini_patchindex = -1;
			int	mini_ptchplistindex = -1;	
			mypoint2f apt = l_plist[fbpt_ofline[i]];
			for (unsigned int j = 0; j < patch_crossindex.size(); ++j)
			{
				//patch_crossindex第一个数是patch的index，之后的数是其plist上的点的index
				int tmp_pchindex = patch_crossindex[j].at(0);
				vector<mypoint2f> pch_plist;
				jointMethod(&plrGraph.patches[tmp_pchindex].pointIndexOfPatch, &pch_plist, 0);
				for (unsigned int k = 1; k < patch_crossindex[j].size(); ++k)//k必须使用1开始！！！
				{
					double smalllen = LineSegmentLength(pch_plist[patch_crossindex[j][k]], apt);
					if (smalllen < mini_length)
					{
						mini_length = smalllen;
						mini_patchindex = tmp_pchindex;
						mini_ptchplistindex = patch_crossindex[j][k];
					}
				}
			}
			if (mini_patchindex != -1)
			{
				pch_ptindex.push_back(make_pair(mini_patchindex, mini_ptchplistindex));//第mini_patchindex个patch 的plist中的第mini_ptchPlistIndex个点距离linede crossPoint最近
			}
		}
		for (unsigned int i = 0; i < fbpt_ofline.size(); ++i) //fbpt_ofline 与 pch_ptindex一一对应
		{
			int crossPatch = pch_ptindex[i].first;
			int crosspchlist_index = pch_ptindex[i].second;
			//计算中间交叉点 cross_pt
			mypoint2f p1_ofline = l_plist[fbpt_ofline[i]];//  lplist中的一小段
			mypoint2f p2_ofline = l_plist[fbpt_ofline[i] - 1];
			vector<mypoint2f> pch_plist;
			vector<int> pch_segPt_num;
			pch_segPt_num.push_back(0);
			for (unsigned int j = 1; j < plrGraph.patches[crossPatch].pointIndexOfPatch.size(); ++j)
			{
				int tmp_pindex1 = plrGraph.patches[crossPatch].pointIndexOfPatch[j - 1];
				int tmp_pindex2 = plrGraph.patches[crossPatch].pointIndexOfPatch[j];
				mypoint2f tmp_p1 = sgraph.nodeList[tmp_pindex1].position;
				mypoint2f tmp_p2 = sgraph.nodeList[tmp_pindex2].position;
				if (!pch_plist.empty())
					pch_plist.pop_back();
				getPListOfaCurve(tmp_p1, tmp_p2, tmp_pindex1, tmp_pindex2, &pch_plist, 0);
				pch_segPt_num.push_back(pch_plist.size() - 1);
			}
			mypoint2f p1_ofpatch = pch_plist[crosspchlist_index];//patch 的plist中的一小段
			mypoint2f p2_ofpatch = pch_plist[crosspchlist_index-1];
			mypoint2f cross_pt = getIntersectPoint(p1_ofline, p2_ofline, p1_ofpatch, p2_ofpatch);//交叉点
			//crosspt在patch的那哪两个点之间，
			int tmp_cr = crosspchlist_index - 1;  int ci = 1;
			while (tmp_cr < pch_segPt_num[ci-1] || tmp_cr > pch_segPt_num[ci])
			{
				ci++;
			}//位于ci-1 ， ci之间
			int left_index = plrGraph.patches[crossPatch].pointIndexOfPatch[ci-1];
			int right_index = plrGraph.patches[crossPatch].pointIndexOfPatch[ci];			
			//sgraph 里添加新点
			ArcNode origiArc = *getArcNode(left_index, right_index);
			EdgeOfPatch origiEdge = *getBorderofPatch(crossPatch, left_index, right_index);
			VertexNode anewnode;
			anewnode.position = cross_pt;
			anewnode.verIndex = sgraph.nodeList.size();
			anewnode.firstEdge = NULL;
			sgraph.nodeList.push_back(anewnode);
			//将该条edge分成俩段, polyline_vec里添加新的line  ,sgraph添加arc
			vector<mypoint2f> patch_seg1, patch_seg2;
				//第一段
			for (unsigned int j = pch_segPt_num[ci - 1]; j <= tmp_cr; ++j)
			{
				patch_seg1.push_back(pch_plist[j]);
			}
			patch_seg1.push_back(cross_pt);
			PolyLine add_pseg_pl1;
			add_pseg_pl1.origi_index = polyline_vec.size();
			add_pseg_pl1.origi_type = 'L';
			add_pseg_pl1.pointlist = patch_seg1;
			polyline_vec.push_back(add_pseg_pl1);
			
			ArcNode* anewarc=new ArcNode;
			anewarc->adjVertex = left_index;
			anewarc->curveIndex = polyline_vec.size() - 1;
			anewarc->curveType = 'L';
			anewarc->lineFlag = 0;
			anewarc->patchesIndex = origiArc.patchesIndex;
			anewarc->next = sgraph.nodeList.back().firstEdge;
			sgraph.nodeList.back().firstEdge = anewarc;
			ArcNode* anewarc_reve = new ArcNode;
			anewarc_reve->adjVertex = sgraph.nodeList.size()-1;
			anewarc_reve->curveIndex = polyline_vec.size() - 1;
			anewarc_reve->curveType = 'L';
			anewarc_reve->lineFlag = 0;
			anewarc_reve->patchesIndex = origiArc.patchesIndex;
			anewarc_reve->next = sgraph.nodeList[left_index].firstEdge;
			sgraph.nodeList[left_index].firstEdge = anewarc_reve;
				//第二段
			patch_seg2.push_back(cross_pt);
			for (unsigned int j = tmp_cr + 1; j < pch_segPt_num[ci]; ++j)
			{
				patch_seg2.push_back(pch_plist[i]);
			}
			PolyLine add_pseg_pl2;
			add_pseg_pl2.origi_index = polyline_vec.size();
			add_pseg_pl2.origi_type = 'L';
			add_pseg_pl2.pointlist = patch_seg2;
			polyline_vec.push_back(add_pseg_pl2);
			ArcNode* anewarc_r = new ArcNode;
			anewarc_r->adjVertex = right_index;
			anewarc_r->curveIndex = polyline_vec.size() - 1;
			anewarc_r->curveType = 'L';
			anewarc_r->lineFlag = 0;
			anewarc_r->patchesIndex = origiArc.patchesIndex;
			anewarc_r->next = sgraph.nodeList.back().firstEdge;
			sgraph.nodeList.back().firstEdge = anewarc_r;
			ArcNode* anewarc_r_reve = new ArcNode;
			anewarc_r_reve->adjVertex = sgraph.nodeList.size() - 1;
			anewarc_r_reve->curveIndex = polyline_vec.size() - 1;
			anewarc_r_reve->curveType = 'L';
			anewarc_r_reve->lineFlag = 0;
			anewarc_r_reve->patchesIndex = origiArc.patchesIndex;
			anewarc_r_reve->next = sgraph.nodeList[right_index].firstEdge;
			sgraph.nodeList[right_index].firstEdge = anewarc_r_reve;
			//patch添加新节点 和 edge
			vector<int>::iterator iter_pvec = plrGraph.patches[crossPatch].pointIndexOfPatch.begin() + ci;//在ci之前添加新节点
			plrGraph.patches[crossPatch].pointIndexOfPatch.insert(iter_pvec, sgraph.nodeList.size() - 1);
				//edge1
			EdgeOfPatch add_edge1;
			add_edge1.startp_index = left_index;
			add_edge1.endp_index = sgraph.nodeList.size() - 1;
			add_edge1.edgeflag = 0;
			add_edge1.adjPatchIndex = origiEdge.adjPatchIndex;
				//edge2
			EdgeOfPatch add_edge2;
			add_edge2.startp_index = sgraph.nodeList.size() - 1;
			add_edge2.endp_index = right_index;
			add_edge2.edgeflag = 0;
			add_edge2.adjPatchIndex = origiEdge.adjPatchIndex;
			vector<EdgeOfPatch>::iterator iter_edges=plrGraph.patches[crossPatch].edges_vec.begin();
			while (iter_edges != plrGraph.patches[crossPatch].edges_vec.end())
			{
				if (iter_edges->startp_index == right_index)
				{
					iter_edges = plrGraph.patches[crossPatch].edges_vec.insert(iter_edges, add_edge2);
					iter_edges = plrGraph.patches[crossPatch].edges_vec.insert(iter_edges, add_edge1);
					break;
				}
				iter_edges++;
			}
			//把l_plist分开，加入polyline_vec。。且patch添加point_line
			vector<mypoint2f> line_seg1;
			if (i == 0)
			{
				//l_plist分成 0-- > fbpt_ofline[i]
				int seg_i = fbpt_ofline[i] - 1;
				for (unsigned int j = 0; j <= seg_i; ++j)
				{
					line_seg1.push_back(l_plist[j]);
				}
				line_seg1.push_back(cross_pt);
			}
			else
			{
				//分成 fbpt_ofline[i]-- > end
				line_seg1.push_back(cross_pt);
				for (unsigned int j = fbpt_ofline[i]; j < l_plist.size(); ++j)
				{
					line_seg1.push_back(l_plist[j]);
				}
			}
			PolyLine add_pl1;
			add_pl1.origi_index = polyline_vec.size();
			add_pl1.origi_type = 'L';
			add_pl1.pointlist = line_seg1;
			polyline_vec.push_back(add_pl1);
			ICurve icv;
			icv.type = 0;
			icv.pt_vec.push_back(sgraph.nodeList.size() - 1);
			if (sgraph.nodeList[newlvec->front()].position == line_seg1.front())
			{
				icv.pt_vec.push_back(newlvec->front());

			}
			else if (sgraph.nodeList[newlvec->back()].position == line_seg1.back())
			{
				icv.pt_vec.push_back(newlvec->back());
			}
			plrGraph.otherlines.push_back(icv);
			ArcNode* arc_l = new ArcNode;
			arc_l->adjVertex = icv.pt_vec[0];
			arc_l->curveIndex = polyline_vec.size() - 1;
			arc_l->curveType = 'L';
			arc_l->lineFlag = 0;
			//arc_l->patchesIndex = origiArc.patchesIndex;
			arc_l->next = sgraph.nodeList[icv.pt_vec[1]].firstEdge;
			sgraph.nodeList[icv.pt_vec[1]].firstEdge = arc_l;
			ArcNode* arc_l_reve = new ArcNode;
			arc_l_reve->adjVertex = icv.pt_vec[1];
			arc_l_reve->curveIndex = polyline_vec.size() - 1;
			arc_l_reve->curveType = 'L';
			arc_l_reve->lineFlag = 0;
			//arc_l_reve->patchesIndex = origiArc.patchesIndex;
			arc_l_reve->next = sgraph.nodeList[icv.pt_vec[0]].firstEdge;
			sgraph.nodeList[icv.pt_vec[0]].firstEdge = arc_l_reve;
			if (origiEdge.adjPatchIndex.size() > 1)
				cout << "?" << endl;
			//mixedTopology
			MixedNode mnode;
			mnode.pl_index = mixedTopology.allNode.size();
			mnode.pl_index = plrGraph.otherlines.size();
			mnode.pl_type = 'l';
			mnode.status = 0;
			mixedTopology.allNode.push_back(mnode);
			MixedArc* marc = new MixedArc;
			marc->adj_index = crossPatch;
			marc->weight = 0;
			marc->next_arc = mixedTopology.allNode[mixedTopology.allNode.size() - 1].firstAdjNode;
			mixedTopology.allNode[mixedTopology.allNode.size() - 1].firstAdjNode = marc;
			MixedArc* marc_reve = new MixedArc;
			marc_reve->adj_index = mixedTopology.allNode.size() - 1;
			marc_reve->weight = 0;
			marc_reve->next_arc = mixedTopology.allNode[crossPatch].firstAdjNode;
			mixedTopology.allNode[crossPatch].firstAdjNode = marc_reve;
		}

	}
}



/*调整流程，混合排序 同 method2,*/
double MyGraph::getAverageWidth(double &w, double &h, double ratio)
{
	int pchsize = plrGraph.patches.size();
	int linesize = plrGraph.isolatedlines.size() + plrGraph.attachlines.size();
	double aver_w = 0;

	/*第一种：对角线*0.05*/
	double w1 = 6553500;//左
	double w2 = -6553500;//右
	double h1 = 6553500;//下
	double h2 = -6553500;//上
	for (int i = 0; i < pchsize; ++i)
	{
		vector<mypoint2f> pplist;
		jointMethod(&plrGraph.patches[i].pointIndexOfPatch, &pplist, 2);
		for (unsigned int j = 0; j < pplist.size(); ++j)
		{
			if (pplist[j].x < w1)
				w1 = pplist[j].x;
			if (pplist[j].x > w2)
				w2 = pplist[j].x;
			if (pplist[j].y < h1)
				h1 = pplist[j].y;
			if (pplist[j].y > h2)
				h2 = pplist[j].y;
		}
	}
	for (int i = 0; i < linesize; ++i)
	{
		vector<mypoint2f> pplist;
		jointMethod(&plrGraph.otherlines[i].pt_vec, &pplist, 2);
		for (unsigned int j = 0; j < pplist.size(); ++j)
		{
			if (pplist[j].x < w1)
				w1 = pplist[j].x;
			if (pplist[j].x > w2)
				w2 = pplist[j].x;
			if (pplist[j].y < h1)
				h1 = pplist[j].y;
			if (pplist[j].y > h2)
				h2 = pplist[j].y;
		}
	}
	//cout << "w1=" << w1 << "," << "w2=" << w2 << endl;
	//cout << "h1=" << h1 << "," << "h2=" << h2 << endl;
	double dui = sqrt(pow((h2 - h1), 2) + pow((w2 - w1), 2));
	aver_w = dui*ratio;
	h = h2 - h1;
	w - w2 - w1;
	cout << "宽=" << h2 - h1 << "," << "长=" << w2 - w1 << endl;
	cout << "对角线=" << dui;
	cout << "对角线* ratio=" << aver_w << endl;
	/*第二种：patch的平均宽度*/
	//vector<double> patch_w;   //顺便获取三个feature，在.h文件中有定义vector<pair<int, vector<double>>> patchfeatures;
	//for (int i = 0; i < pchsize; ++i)
	//{
	//	vector<double> tmp;
	//	patchfeatures.push_back(make_pair(i, tmp));
	//	if (plrGraph.patches[i].type >= 0)
	//	{
	//		vector<mypoint2f> total_plist;
	//		jointMethod(&plrGraph.patches[i].pointIndexOfPatch, &total_plist, 0);
	//		//2、长度，方法二，计算出最大矩形bounding box
	//		double elongation = 0;
	//		vector<mypoint2f> plist_box;
	//		mypoint2f rec_center(0, 0);
	//		getBoundRectangle(&total_plist, rec_center, &plist_box);
	//		double box_w = fabs(plist_box[0].x - plist_box[1].x);
	//		double box_h = fabs(plist_box[1].y - plist_box[2].y);
	//		if (box_w < box_h)
	//		{
	//			elongation = box_w / box_h;
	//			patch_w.push_back(box_w);
	//		}
	//		else
	//		{
	//			elongation = box_h / box_w;
	//			patch_w.push_back(box_h);
	//		}
	//	}
	//}
	//aver_w = GetAverange(patch_w);
	//cout << "average width:" << aver_w << endl;

	return aver_w;
}
double MyGraph::ConstructMixedTopology(double w_ratio)//线的宽度设置为对角线的0.05/patch的平均宽度，角度设置为30，若宽度=0默认是30
{
	cout << "ConstructMixedTopology" << endl;
	//MixedGraph mixedTopology,构建混合拓扑结构，在mixedTopology.allNode中存放的先是patch 然后是attachline,然后是isolateline
	/*1.先建立patch的topology（同上），（patch与patch）*/
	int pchsize = plrGraph.patches.size();
	int attachsize = plrGraph.attachlines.size();
	int isolatedsize = plrGraph.isolatedlines.size();
	int linesize = attachsize + isolatedsize;   //即otherline的大小
	mixedTopology.mixnode_number = pchsize + attachsize + isolatedsize;
		
	for (int i = 0; i < pchsize; ++i)
	{
		//注意不管patch的type是不是-1 ，都存入mix topology(vec)中。若是type=-1（unavailiable），这就是一个孤立的点，可以直接用patchindex找到
		MixedNode mnode;
		mnode.pl_type = 'p';
		mnode.pl_index = plrGraph.patches[i].patch_index;////?????
		mnode.firstAdjNode = NULL;
		mnode.status = 0;
		if (plrGraph.patches[i].type >= 0)
		{
			vector<int> adjpatch;
			for (unsigned int j = 0; j < plrGraph.patches[i].edges_vec.size(); ++j)
			{
				adjpatch.insert(adjpatch.end(), plrGraph.patches[i].edges_vec[j].adjPatchIndex.begin(), plrGraph.patches[i].edges_vec[j].adjPatchIndex.end());
			}
			deleRepeatPoint(&adjpatch);
			for (unsigned int j = 0; j < adjpatch.size(); ++j)
			{
				if (adjpatch[j] != i)
				{
					MixedArc *anarc = new MixedArc;
					anarc->adj_index = adjpatch[j];////?????
					anarc->weight = -1;
					anarc->next_arc = mnode.firstAdjNode;
					mnode.firstAdjNode = anarc;
				}
			}
		}
		else
		{
			mnode.status = -1;
		}
		mixedTopology.allNode.push_back(mnode);
	}
	/*2.面线，将attachline加入到mixedTopology。（patch与attachline）*/
	//遍历每个patch，查看point_line是否为空，若不为空说明有attachline,将attachline加入到mixedTopology
	double pic_w = 0; double pic_h = 0;
	double patch_averWidth = getAverageWidth(pic_w, pic_h, w_ratio);//用patch的平均width当做（1）两条线最短距离的阈值patch_averWidth（2）attachline的长度小于这个值就删除（status=-1）？？？
	if (patch_averWidth <= 0)
		patch_averWidth = 20;   //如果没有average width=0，默认设置为30
	mycircle width_cir;
	width_cir.cx =  file_width / 2;// pic_w / 2;
	width_cir.cy = file_height / 2;// pic_h / 2;
	width_cir.radius = patch_averWidth/2;
	circle_vec.push_back(width_cir);  //（画圆）patch 的平均宽度
	for (int i = 0; i < attachsize; ++i)
	{
		MixedNode mnode;
		mnode.pl_type = 'l';
		mnode.pl_index = i + pchsize;//!!
		mnode.firstAdjNode = NULL;

		vector<mypoint2f> tmp_plist;
		jointMethod(&plrGraph.attachlines[i].pt_vec, &tmp_plist, 0);
		double perimeter_l=getPerimeterofVec(&tmp_plist);
		if (perimeter_l < 1) //attachline的长度小于这个值就删除（status=-1）!!
		{
			mnode.status = -1;
			plrGraph.otherlines[i].type = -1;//!!!!!!!!!!!!!!!
			plrGraph.attachlines[i].type = -1;
		}
		else
			mnode.status = 0;
		mixedTopology.allNode.push_back(mnode);
	}
	for (int i = 0; i < pchsize; ++i)
	{
		if (plrGraph.patches[i].type >= 0)
		{
			if (!plrGraph.patches[i].point_line.empty())
			{
				for (unsigned int j = 0; j < plrGraph.patches[i].point_line.size(); ++j)
				{
					int l_index = plrGraph.patches[i].point_line[j].second;//attachline ，在attachline中的index
					if (mixedTopology.allNode[l_index + pchsize].status >= 0)//status<0说明这条attachlin短，小于阈值patch_averWidth！！！
					{
						//line（l_index + pchsize）-->patch（i）
						MixedArc *marc1 = new MixedArc;
						marc1->adj_index = i;                   //patch 的下标
						marc1->weight = -1;
						marc1->next_arc = mixedTopology.allNode[l_index + pchsize].firstAdjNode;
						mixedTopology.allNode[l_index + pchsize].firstAdjNode = marc1;
						//patch(i)-->line(l_index + pchsize)
						MixedArc *marc2 = new MixedArc;
						marc2->adj_index = l_index + pchsize;   //line 的下标
						marc2->weight = -1;
						marc2->next_arc = mixedTopology.allNode[i].firstAdjNode;
						mixedTopology.allNode[i].firstAdjNode = marc2;
					}			
				}
			}
		}
	}
	/*3、线线相交+相近 （attach与attach相交，attach与attach相近, attach与isolated相近）*/
	if (!plrGraph.isolatedlines.empty())
	{
		for (int i = 0; i < isolatedsize; ++i)
		{
			MixedNode mnode;
			mnode.pl_type = 'i';
			mnode.pl_index = i + pchsize + attachsize;//!!
			mnode.firstAdjNode = NULL;
			mnode.status = 0;
			mixedTopology.allNode.push_back(mnode);
		}
	}
	for (int i = 0; i < linesize - 1; ++i)//n^2
	{
		if (mixedTopology.allNode[i + pchsize].status >= 0)  //status >= 0
		{
			vector<int> l1_vec = plrGraph.otherlines[i].pt_vec;
			vector<mypoint2f> l1_plist;
			jointMethod(&plrGraph.otherlines[i].pt_vec, &l1_plist, 0);
			for (int j = i + 1; j < linesize; ++j)
			{
				if (mixedTopology.allNode[j + pchsize].status >= 0)  //status >= 0
				{
					vector<int> l2_vec = plrGraph.otherlines[j].pt_vec;
					vector<mypoint2f> l2_plist;
					jointMethod(&plrGraph.otherlines[j].pt_vec, &l2_plist, 0);
					int n_pindex1 = -1; int n_pindex2 = -1;
					double proxi = getProximity(&l1_plist, &l2_plist, n_pindex1, n_pindex2);//找两条curve距离最近的点，并返回他们的位置pindex
					double conti = getContinuity(&l1_plist, &l2_plist, n_pindex1, n_pindex2);//最近距离处的夹角
					if (proxi < patch_averWidth && conti<20 && conti>=0)  // 夹角小于20 ？？！！！！
					{
						MixedArc *marc1 = new MixedArc;
						marc1->adj_index = i + pchsize;//attachline(i) 的下标
						marc1->weight = -1;
						marc1->next_arc = mixedTopology.allNode[j + pchsize].firstAdjNode;//j+i
						mixedTopology.allNode[j + pchsize].firstAdjNode = marc1;
						MixedArc *marc2 = new MixedArc;
						marc2->adj_index = j + pchsize;//attachline(j) 的下标
						marc2->weight = -1;
						marc2->next_arc = mixedTopology.allNode[i + pchsize].firstAdjNode;//i+j
						mixedTopology.allNode[i + pchsize].firstAdjNode = marc2;
					}
				}			
			}
		}
	}

	//总结：mixedTopology.allNode由patch，attachline,isolatedline按顺序组成。
	//mixnode.firstAdjNode=NULL : patch.type=-1或者patch 没有neighbor，line同理
	//mixnode.pl_type："p,l,i"
	cout << "ConstructMixedTopology  done" << endl;
	return patch_averWidth;
}
MixedArc* MyGraph::getMixedArc(int mindex, int nei_mindex)
{
	MixedArc *tmparc = NULL;
	MixedArc* arc = mixedTopology.allNode[mindex].firstAdjNode;
	while (arc != NULL)
	{
		if (arc->adj_index == nei_mindex)
		{
			tmparc = arc;
			break;
		}
		arc = arc->next_arc;
	}
	return tmparc;
}
void MyGraph::removeMixedNode(int hostnode, int delnode)
{
	MixedArc *arcnode = mixedTopology.allNode[hostnode].firstAdjNode;
	MixedArc *prior_arc = NULL;
	while (arcnode)
	{
		if (arcnode->adj_index == delnode)
		{
			if (prior_arc == NULL)
			{
				mixedTopology.allNode[hostnode].firstAdjNode = arcnode->next_arc;
			}
			else
			{
				prior_arc->next_arc = arcnode->next_arc;
			}
			delete arcnode;
			arcnode = NULL;
			break;
		}
		prior_arc = arcnode;
		arcnode = arcnode->next_arc;
	}
}
void MyGraph::getMixedNeighbor(int mindex, vector<int>* m_nei)
{
	MixedArc* arc = mixedTopology.allNode[mindex].firstAdjNode;
	while (arc != NULL)
	{
		m_nei->push_back(arc->adj_index);
		arc = arc->next_arc;
	}
}

void MyGraph::MixSimilarityRange(int simply_num, double patch_averWidth, int similarity_type,double similarity_t)//similarity_type=0/1
{
	////划分mixGraph的子图
	//int mixsize = mixedTopology.mixnode_number;
	//vector<pair<int, string>> mix_subgraphs;
	//int *visitStatus = new int[mixsize];//0:unvisited,1:push in queue,2: out of queue
	//for (int i = 0; i < mixsize; ++i)
	//{
	//	visitStatus[i] = 0;
	//}
	//for (int i = 0; i < mixsize; ++i)//breadth-first traversal  BFT
	//{
	//	if (visitStatus[i] == 0 && patchTopology.allpatches[i].patchType >= 0)
	//	{
	//		if (patchTopology.allpatches[i].firstAdjPatch != NULL)
	//		{
	//			int arcCount = 0;
	//			int verCount = 0;
	//			int tmp_ver = 0;
	//			queue<int> verQueue;
	//			verQueue.push(i);
	//			visitStatus[i] = 1;
	//			while (!verQueue.empty())
	//			{
	//				tmp_ver = verQueue.front();
	//				verQueue.pop();
	//				verCount = verCount + 1;
	//				visitStatus[tmp_ver] = 2;
	//				PatchArc *tmp_arc = patchTopology.allpatches[tmp_ver].firstAdjPatch;
	//				while (tmp_arc)
	//				{
	//					if (visitStatus[tmp_arc->adj_patchindex] == 0)
	//					{
	//						verQueue.push(tmp_arc->adj_patchindex);
	//						visitStatus[tmp_arc->adj_patchindex] = 1;
	//						arcCount = arcCount + 1;
	//					}
	//					else if (visitStatus[tmp_arc->adj_patchindex] == 1)
	//					{
	//						arcCount = arcCount + 1;
	//					}
	//					tmp_arc = tmp_arc->nextPatch;
	//				}
	//			}
	//			if (verCount - 1 == arcCount)
	//			{
	//				mix_subgraphs.push_back(make_pair(i, "tree"));
	//			}
	//			else
	//			{
	//				mix_subgraphs.push_back(make_pair(i, "loop"));
	//			}
	//		}
	//		else
	//			mix_subgraphs.push_back(make_pair(i, "node"));
	//	}
	//}
	//delete[] visitStatus;
	//visitStatus = NULL;
	
	cout << "MixSimilarityRange" << endl;  
	/*1、计算 patch和line 的FD*/
	int pchsize = plrGraph.patches.size();
	int attachsize = plrGraph.attachlines.size();
	int isolatedsize = plrGraph.isolatedlines.size();
	int linesize = attachsize + isolatedsize;   //即otherline的大小
	int mixsize = pchsize + attachsize + isolatedsize;/////////////////////////////////////
	//int mixsize = pchsize;/////////////////////////////

	vector<vector<double>> mix_feature;//顺序 与mixedtopology相同，即先patch，attachline，isolatedline
	vector<vector<double>> mix_FDs;
	for (int i = 0; i < mixsize; ++i)
	{
		vector<double> tmp;
		mix_feature.push_back(tmp);     //赋初值 三个feature
		mix_FDs.push_back(tmp);          //赋初值 FD
		if (mixedTopology.allNode[i].status >= 0)
		{
			if (mixedTopology.allNode[i].pl_type == 'p')
			{
				vector<mypoint2f> total_plist;
				jointMethod(&plrGraph.patches[i].pointIndexOfPatch, &total_plist, 0);
				total_plist.pop_back();  //!!!	
				vector<double> pchFD; //patch 的FD
				getFDFeature(&total_plist, &plrGraph.patches[i].pointIndexOfPatch,i, simply_num, similarity_type, &pchFD);//该函数内部包含了plist的化简，参数：第0种FD计算方法，或者第1种
				mix_FDs[i] = pchFD;

				//vector<double> fvec;//特征向量   /////////////////////////////
				//getThreeFeature(&fvec, &total_plist, plrGraph.patches[i].area, plrGraph.patches[i].perimeter);////////////////////
				//mix_feature[i] = fvec;////////////////////////////////
			}
			else
			{
				vector<mypoint2f> attach_plist;
				jointMethod(&plrGraph.otherlines[i-pchsize].pt_vec, &attach_plist, 0);  //attach + isolate
				//vector<double> line_tf;////////////////////////////////////
				//getThreeFeature_forline(&line_tf, &attach_plist);//////////////////////
				//mix_feature[i] = line_tf;////////////////////////////////
				
				vector<mypoint2f> patch_plist;//point list围成patch
				patch_plist.insert(patch_plist.end(), attach_plist.begin(), attach_plist.end());
				for (int ppi = (int)attach_plist.size()-2; ppi >0; --ppi)//patch_plist 的 front！=back
				{
					patch_plist.push_back(attach_plist[ppi]);
				}
				vector<int> line_pvec;//point vec围成patch
				line_pvec = plrGraph.otherlines[i - pchsize].pt_vec;
				for (int j = plrGraph.otherlines[i - pchsize].pt_vec.size()-2; j >=0; --j)//line_pvec的front = back
				{
					line_pvec.push_back(plrGraph.otherlines[i - pchsize].pt_vec[j]);
				}
				vector<double> l_FD; //line 的FD
				getFDFeature(&patch_plist, &line_pvec, i, simply_num, similarity_type, &l_FD);//该函数内部包含了plist的化简。第三个参数写i或者0都可以。参数similarity_type：第0种FD计算方法，或者第1种
				mix_FDs[i] = (l_FD);
			}
		}
	}
	/*2、计算mixtopology.allnode的weight。分开计算: patch与patch合并，patch与line合并，line和line合并*/
	//计算每个patch与neighbor合并后的特征值, dis(center_ft,new_ft)->patchTopology.allpatches[i].arc.weight中 (//patch_subgraphs.push_back(make_pair(i, "tree"));  )
	ofstream ps_file1("E:/svgfile/p+l_equal/p_similarity_dis.txt");
	ps_file1.is_open();
	vector<double> pdis;
	ofstream ls_file2("E:/svgfile/p+l_equal/l_similarity_dis.txt");
	ls_file2.is_open();
	vector<double> ldis;
	ofstream pl_file3("E:/svgfile/p+l_equal/p+l_similarity_dis.txt");
	pl_file3.is_open();
	vector<double> pl_dis;

	vector<double> tt_w;//记录权值
	vector<int> weight_pindex;//记录某个patch的最小weight对应的patchindex
	for (int i = 0; i < mixsize; ++i)  //mixedTopology 从0->pchsize
	{
		tt_w.push_back(-1);
		weight_pindex.push_back(-1);
		if (mixedTopology.allNode[i].status < 0)	//if (plrGraph.patches[i].type >= 0)
		{
			continue;
		}
		if (mixedTopology.allNode[i].pl_type == 'p')  // 1、这个点是patch
		{
			MixedArc* marc = mixedTopology.allNode[i].firstAdjNode;
			double min_dis = 6553500;
			int min_pindex = -1;
			while (marc != NULL)//每一个arc 是一个neighbor patch, 中心patch与neighborpatch合并，然后计算合并后的feature，然后计算该newfeature与patch_tmpfeature[i].second的二范数，存入arc 的weight
			{
				if (marc->weight >= 0)
				{
					if (marc->weight < min_dis)
					{
						min_dis = marc->weight;
						min_pindex = marc->adj_index;
					}
					marc = marc->next_arc;
					continue;
				}
				char nodeType = mixedTopology.allNode[marc->adj_index].pl_type;//判断neighbor 的类型，
				if (nodeType == 'p')	    /*1.1、 p p 合并*/
				{
					vector<int> combine_pvec;//合并的pointlist
					combineTwoPatch_getpvec(i, marc->adj_index, &combine_pvec);//combine_plist  front()！= back()	
					if (!combine_pvec.empty())//是NULL的话说明不能合并
					{
						vector<mypoint2f> combine_plist;
						jointMethod(&combine_pvec, &combine_plist, 0);

						//vector<double> threeft;//////////////////////////////////////
						//getThreeFeature(&threeft, &combine_plist, 0, 0);//////////////////////////////////

						//计算新patch的周长+FD
						vector<double> new_ft;
						combine_plist.pop_back();  //front()！= back()	
						getFDFeature(&combine_plist, &combine_pvec, -1, simply_num, similarity_type, &new_ft);				
						double dis1 = getEuclidDistance(&mix_FDs[i], &new_ft);
						double dis2 = getEuclidDistance(&mix_FDs[marc->adj_index], &new_ft);
						//double dis1 = getEuclidDistance(&mix_FDs[i], &new_ft)*0.4+ getEuclidDistance(&mix_feature[i], &threeft)*0.6;//////////////////////////////////
						//double dis2 = getEuclidDistance(&mix_FDs[marc->adj_index], &new_ft)*0.4 + getEuclidDistance(&mix_feature[marc->adj_index], &threeft)*0.6;//////////////////////////
						if (marc->weight < 0)//parc->weight初始化是赋值为-1，在选取最小值的注意把-1排除!
						{
							marc->weight = dis1 < dis2 ? dis1 : dis2;//取较小的那个值
						}
						else
						{
							////取parc->weight,dis1,dis2三者的最小值
								//if (dis1 < parc->weight)
								//{
								//	if (dis2 < dis1)
								//		parc->weight = dis2;
								//	else
								//		parc->weight = dis1;
								//}
								//else 
								//{
								//	if (dis2 < parc->weight)
								//		parc->weight = dis2;
								//}
							cout << "??" << endl;//不可能有这种情况
						}
						MixedArc* parc_reve = getMixedArc(marc->adj_index, i);//patchTopology.allpatches[parc->adj_patchindex].firstAdjPatch;
						if (parc_reve != NULL)
							parc_reve->weight = marc->weight;
						else
							cout << "patch arc=NULL" << endl;
						if (marc->weight < min_dis)
						{
							min_dis = marc->weight;
							min_pindex = marc->adj_index;
						}
						pdis.push_back(marc->weight);  //输出
						marc = marc->next_arc;
					}
					else
					{
						//虽然两patch相邻，但是不能merge，从neighbor 中删除
						int tmp_nindex = marc->adj_index;
						marc = marc->next_arc;
						removeMixedNode(i, tmp_nindex);
						removeMixedNode(tmp_nindex, i);
					}
				}
				else if (nodeType == 'l')		/* 1.2、p l 合并  how?????????? */
				{
					vector<mypoint2f> modi_patch;
					bool pl_flag = combinePL_getplist(i, marc->adj_index - pchsize, &modi_patch);
					if (pl_flag == true)//是NULL的话说明不能合并
					{
						modi_patch.pop_back();  //front()！= back()	
						vector<int> empty_vec;
						vector<double> new_ft;
						getFDFeature(&modi_patch, &empty_vec, i, simply_num, similarity_type, &new_ft);//从p+l得到的FD，注意第三个参数是原patch的index
						double dis1 = getEuclidDistance(&mix_FDs[i], &new_ft);
						marc->weight = dis1;//marc->type是一个line，且现在weight=0，不用判断line和newpatch的dis,直接赋值dis1
						MixedArc* parc_reve = getMixedArc(marc->adj_index, i);//patchTopology.allpatches[parc->adj_patchindex].firstAdjPatch;
						if (parc_reve != NULL)
							parc_reve->weight = marc->weight;
						if (marc->weight < min_dis)
						{
							min_dis = marc->weight;
							min_pindex = marc->adj_index;
						}
						pl_dis.push_back(marc->weight);  //输出patch与line的merge后的dis
						marc = marc->next_arc;
					}
					else
					{
						//虽然patch 和line相邻，但是不能merge，令权值最大
						marc->weight = 6553500;//marc->type是一个line，且现在weight=0，不用判断line和newpatch的dis,直接赋值dis1
						MixedArc* parc_reve = getMixedArc(marc->adj_index, i);//patchTopology.allpatches[parc->adj_patchindex].firstAdjPatch;
						if (parc_reve != NULL)
							parc_reve->weight = 6553500;
						marc = marc->next_arc;
					}
				}
				//parc = parc->next_arc;
			} //while done
			if (min_dis == 6553500)//patch 没有neighbor，或者有neighbor但是neighbor是line
			{
			}
			tt_w[i] = min_dis;
			weight_pindex[i] = min_pindex;
		}
		else if (mixedTopology.allNode[i].pl_type == 'l' || mixedTopology.allNode[i].pl_type=='i')  // 2、这个点是attachline
		{
			MixedArc* parc = mixedTopology.allNode[i].firstAdjNode;
			double min_dis = 6553500;
			int min_pindex = -1;
			while (parc != NULL)//每一个arc 是一个neighbor patch, 中心patch与neighborpatch合并，然后计算合并后的feature，然后计算该newfeature与patch_tmpfeature[i].second的二范数，存入arc 的weight
			{
				if (parc->weight >= 0)
				{
					if (parc->weight < min_dis)
					{
						min_dis = parc->weight;
						min_pindex = parc->adj_index;
					}
					parc = parc->next_arc;
					continue;
				}
				char nodeType = mixedTopology.allNode[parc->adj_index].pl_type;//判断neighbor 的类型，
				if (nodeType == 'p')		//2.1  l+p
				{
					vector<mypoint2f> modi_patch;
					bool pl_flag = combinePL_getplist(i, parc->adj_index-pchsize, &modi_patch);
					if (pl_flag == true)//是NULL的话说明不能合并
					{
						modi_patch.pop_back();  //front()！= back()	
						vector<int> empty_vec;
						vector<double> new_ft;
						getFDFeature(&modi_patch, &empty_vec, parc->adj_index, simply_num, similarity_type, &new_ft);//从p+l得到的FD，注意第三个参数是原patch的index
						double dis1 = getEuclidDistance(&mix_FDs[i], &new_ft);
						parc->weight = dis1;
						MixedArc* parc_reve = getMixedArc(parc->adj_index, i);//patchTopology.allpatches[parc->adj_patchindex].firstAdjPatch;
						if (parc_reve != NULL)
							parc_reve->weight = parc->weight;
						if (parc->weight < min_dis)
						{
							min_dis = parc->weight;
							min_pindex = parc->adj_index;
						}
						pl_dis.push_back(parc->weight);  //输出patch与line的merge后的dis
						parc = parc->next_arc;
					}
					else
					{
						//虽然patch 和line相邻，但是不能merge，令权值最大
						parc->weight = 6553500;//parc->type是一个line，且现在weight=0，不用判断line和newpatch的dis,直接赋值dis1
						MixedArc* parc_reve = getMixedArc(parc->adj_index, i);//patchTopology.allpatches[parc->adj_patchindex].firstAdjPatch;
						if (parc_reve != NULL)
							parc_reve->weight = 6553500;
						parc = parc->next_arc;
					}
				}
				else	//2.2  l+l /l+i
				{
					vector<mypoint2f> newplist;
					bool c_flag=combineTwoCurves_getplist(i - pchsize, parc->adj_index - pchsize, &newplist, patch_averWidth);
					if (c_flag == true)
					{
						vector<mypoint2f> newpatch_plist;
						newpatch_plist.insert(newpatch_plist.end(), newplist.begin(), newplist.end());
						for (int ppi = (int)newplist.size() - 2; ppi >0; --ppi) //front！=back
						{
							newpatch_plist.push_back(newplist[ppi]);
						}
						vector<int> line_pvec;//line的patch的pvec,没有是空的
						vector<double> new_FD; //line 的FD
						getFDFeature(&newpatch_plist, &line_pvec, -1, simply_num, similarity_type, &new_FD); //参数：第0种FD计算方法，或者第1种
						double dis1 = getEuclidDistance(&mix_FDs[i], &new_FD);
						double dis2 = getEuclidDistance(&mix_FDs[parc->adj_index], &new_FD);
						parc->weight = dis1 < dis2 ? dis1 : dis2;//取较小的那个值
						MixedArc* parc_reve = getMixedArc(parc->adj_index, i);//patchTopology.allpatches[parc->adj_patchindex].firstAdjPatch;
						if (parc_reve != NULL)
							parc_reve->weight = parc->weight;
						else
							cout << "patch arc=NULL" << endl;
						if (parc->weight < min_dis)
						{
							min_dis = parc->weight;
							min_pindex = parc->adj_index;
						}
						ldis.push_back(parc->weight);  //输出
						parc = parc->next_arc;
					}
					else
					{				
						//虽然line 和line相邻，但是不能merge
						int tmp_adj = parc->adj_index;
						parc = parc->next_arc;
						removeMixedNode(i, tmp_adj);
						removeMixedNode(tmp_adj, i);
					}
				}
				//parc = parc->next_arc;
			}
			if (min_dis == 6553500)//line 没有neighbor，或者有neighbor但是 是patch
			{
			}
			tt_w[i] = min_dis;
			weight_pindex[i] = min_pindex;
		}				
	}
	//输出patch 和line 的similaritydistance
	vector<double>::iterator iter_tp = pdis.begin();
	while (iter_tp != pdis.end())
	{
		vector<double>::iterator iter_tp_inn = iter_tp + 1;
		while (iter_tp_inn != pdis.end())
		{
			if (*iter_tp == *iter_tp_inn)
				iter_tp_inn = pdis.erase(iter_tp_inn);
			else
				iter_tp_inn++;
		}
		iter_tp++;
	}
	iter_tp = ldis.begin();
	while (iter_tp != ldis.end())
	{
		vector<double>::iterator iter_tp_inn = iter_tp + 1;
		while (iter_tp_inn != ldis.end())
		{
			if (*iter_tp == *iter_tp_inn)
				iter_tp_inn = ldis.erase(iter_tp_inn);
			else
				iter_tp_inn++;
		}
		iter_tp++;
	}
	iter_tp = pl_dis.begin();
	while (iter_tp != pl_dis.end())
	{
		vector<double>::iterator iter_tp_inn = iter_tp + 1;
		while (iter_tp_inn != pl_dis.end())
		{
			if (*iter_tp == *iter_tp_inn)
				iter_tp_inn = pl_dis.erase(iter_tp_inn);
			else
				iter_tp_inn++;
		}
		iter_tp++;
	}
	sort(pdis.begin(), pdis.end());
	sort(ldis.begin(), ldis.end());
	sort(pl_dis.begin(), pl_dis.end());
	cout << "patch的mini distance个数：" << pdis.size() << endl;
	cout << "line的mini distance个数：" << ldis.size() << endl;
	cout << "patch+line的mini distance个数：" << pl_dis.size() << endl;
	for (unsigned int i = 0; i < pdis.size(); ++i)
	{
		ps_file1 << pdis[i] << " ";
	}
	for (unsigned int i = 0; i < ldis.size(); ++i)
	{
		ls_file2 << ldis[i] << " ";
	}
	for (unsigned int i = 0; i < pl_dis.size(); ++i)
	{
		pl_file3 << pl_dis[i] << " ";
	}
	ps_file1.close();
	ls_file2.close();
	pl_file3.close();
	//tt_w 排序
	//建堆，tt_w保存了从第0个点到底n个点weight值，堆heap_vec保存每个patch的权值在tt_w中的下标，也是patch 的index
	vector<int> heap_vec;
	vector<int> ele_index;
	for (unsigned int i = 0; i < tt_w.size(); ++i)
	{
		ele_index.push_back(-1);
		if (mixedTopology.allNode[i].status >= 0)
		{
			heap_vec.push_back(i);
		}
		//heap_vec.push_back(i);
		//ele_index.push_back(0);
	}
	//建堆
	int h_size = heap_vec.size();
	for (int i = h_size / 2; i >= 1; i--)
	{
		heap_siftdown(i, &heap_vec, &tt_w, &ele_index);
	}
	//建立索引，知道某个权值的下标，可以直接找到在heap的位置，位置是从1开始，注意更新
	for (int i = 0; i < h_size; ++i)
	{
		ele_index[heap_vec[i]] = i;// ele_index
	}
	//////////////////
	//h_size是目前堆中的patch的数量，每次合并一个patch,h_size--，这是以数量为停止条件
	int dele_count_line = 0;
	int dele_count_patch = 0;
	int dele_count_pl = 0;
	while (tt_w[heap_vec.front()] < similarity_t && tt_w[heap_vec.front()]>=0)
	{
		int pindex = heap_vec[0];
		if (mixedTopology.allNode[pindex].status<0)
		{
			int tmp_back = heap_vec.back();
			int tmp_mini = heap_vec[0];
			tt_w[tmp_mini] = -1;
			ele_index[tmp_mini] = -1;
			ele_index[tmp_back] = 0;//!!!
			heap_vec[0] = tmp_back;
			heap_siftdown(1, &heap_vec, &tt_w, &ele_index);//下沉调整
			heap_vec.pop_back();
			h_size--;
			continue;
		}
		if (tt_w[pindex] == 6553500)
		{
			cout << "left:" << h_size <<  endl;
			break;//若final_pnum设置的过小，使weight=65535，说明剩下的patch都没有neighbor,可以直接跳出循环，此时最小weight对应的patchindex， weight_pindex[pindex]=-1
		}
		int combine_pindex = -1;
		MixedArc* marc = mixedTopology.allNode[pindex].firstAdjNode;
		while (marc != NULL)
		{
			if (marc->weight == tt_w[pindex])
			{
				combine_pindex = marc->adj_index;				
				break;
			}
			marc = marc->next_arc;
		}
		int ttt = weight_pindex[pindex];
		if (ttt != combine_pindex)
			cout << "combine_pindex error";
		combine_pindex = weight_pindex[pindex];//第pindex个patch最小权值对应的那个patch的index
		
		//1、merge 处理
		if (mixedTopology.allNode[pindex].pl_type == 'p')
		{
			if (mixedTopology.allNode[combine_pindex].pl_type == 'p')//pp
			{
				dele_count_patch++;
				bool succ = combineTwoPatch(combine_pindex, pindex);//combine_pindex的周长，面积，center也会更新
				if (succ == false)
					cout << "combine patch error" << endl;
			}
			else//pl
			{
				string fileadd = "pl_";
				fileadd.append(to_string(dele_count_pl));
				string storagePath = "E:/svgfile/p+l_equal/p+l/";
				vector<double> tmp_v;
				tmp_v.push_back(0); tmp_v.push_back(0); tmp_v.push_back(file_width); tmp_v.push_back(file_height);
				//jointCurveTPatch();
				//jointPolylines();

				combinePatchLine(pindex, combine_pindex - pchsize);//改变的patch在该函数中 被记录下来，还有line
				dele_count_pl++;
	
				//writeSVGFile_step(fileadd, storagePath, file_width, file_height, tmp_v, file_view);
			}
		}
		else//注意line的index，需要减去pchsize
		{
			if (mixedTopology.allNode[combine_pindex].pl_type == 'p')//lp
			{
				string fileadd = "pl_";
				fileadd.append(to_string(dele_count_pl));
				string storagePath = "E:/svgfile/p+l_equal/p+l/";
				vector<double> tmp_v;
				tmp_v.push_back(0); tmp_v.push_back(0); tmp_v.push_back(file_width); tmp_v.push_back(file_height);
				//jointCurveTPatch();
				//jointPolylines();

				combinePatchLine(combine_pindex, pindex - pchsize);//改变的patch在该函数中 被记录下来，还有line
				dele_count_pl++;
				
				//writeSVGFile_step(fileadd, storagePath, file_width, file_height, tmp_v, file_view);
			}
			else//ll
			{
				dele_count_line++;
				vector<mypoint2f> nplist;
				combineTwoCurves(combine_pindex - pchsize, pindex - pchsize, &nplist);
				//checkIntersect(pindex - pchsize, &plrGraph.otherlines[pindex - pchsize].pt_vec, &origi_neiPatch);
			}
		}		
		//2.先删除堆顶的pindex
		if (mixedTopology.allNode[pindex].pl_type == 'p' && mixedTopology.allNode[combine_pindex].pl_type != 'p') //pindex-->patch , combine_pindex-->l/i
		{
			//交换pindex和combine_pindex  但是不从堆中pop掉
			int exchange_index = pindex;
			pindex = combine_pindex;
			combine_pindex = exchange_index;
		}
		else
		{
			int tmp_back = heap_vec.back();
			int tmp_mini = heap_vec[0];
			tt_w[tmp_mini] = -1;
			ele_index[tmp_mini] = -1;
			ele_index[tmp_back] = 0;//!!!
			heap_vec[0] = tmp_back;
			heap_siftdown(1, &heap_vec, &tt_w, &ele_index);//下沉调整
			heap_vec.pop_back();
			h_size--;
		}
		//int tmp_back = heap_vec.back();
		//int tmp_mini = heap_vec[0];
		//tt_w[tmp_mini] = -1;
		//ele_index[tmp_mini] = -1;
		//ele_index[tmp_back] = 0;//!!!
		//heap_vec[0] = tmp_back;
		//heap_siftdown(1, &heap_vec, &tt_w, &ele_index);//下沉调整
		//heap_vec.pop_back();
		//h_size--;

		//3.更新combine_pindex的FD /combine_pindex-pchsize	
		vector<mypoint2f> total_plist;
		vector<double> tmppchFD;
		if (mixedTopology.allNode[combine_pindex].pl_type == 'p')
		{
			jointMethod(&plrGraph.patches[combine_pindex].pointIndexOfPatch, &total_plist, 0); //total_plist首 = 尾
			total_plist.pop_back();			
			getFDFeature(&total_plist, &plrGraph.patches[combine_pindex].pointIndexOfPatch, combine_pindex, simply_num, similarity_type, &tmppchFD);
			
			//vector<double> nw_ft;//////////////////////////////////
			//getThreeFeature(&nw_ft,&total_plist,0,0);///////////////////////////
			//mix_feature[combine_pindex].clear();////////////////////////////////////
			//mix_feature[combine_pindex] = nw_ft;///////////////////////////////////////
		}
		else
		{
			vector<mypoint2f> attach_plist;
			jointMethod(&plrGraph.otherlines[combine_pindex-pchsize].pt_vec, &attach_plist, 0);  //attach + isolate
			total_plist.insert(total_plist.end(), attach_plist.begin(), attach_plist.end());
			for (int ppi = (int)attach_plist.size() - 2; ppi > 0; --ppi) //front!=back
			{
				total_plist.push_back(attach_plist[ppi]);
			}
			vector<int> line_pvec;
			line_pvec = plrGraph.otherlines[combine_pindex - pchsize].pt_vec;
			for (int j = plrGraph.otherlines[combine_pindex - pchsize].pt_vec.size() - 2; j >= 0; --j)
			{
				line_pvec.push_back(plrGraph.otherlines[combine_pindex - pchsize].pt_vec[j]);
			}
			getFDFeature(&total_plist, &line_pvec, -1, simply_num, similarity_type, &tmppchFD);//该函数内部包含了plist的化简，参数：第0种FD计算方法，或者第1种
		}
		////1、圆度
		//double circularity = 4 * 3.1415926*plrGraph.patches[combine_pindex].area / pow(plrGraph.patches[combine_pindex].perimeter, 2);
		////2、长度，计算出最大矩形bounding box
		//double elongation = 0;
		//vector<mypoint2f> plist_box;
		//mypoint2f rec_center(0, 0);
		//getBoundRectangle(&total_plist, rec_center, &plist_box);
		//double box_w = fabs(plist_box[0].x - plist_box[1].x);
		//double box_h = fabs(plist_box[1].y - plist_box[2].y);
		//if (box_w < box_h)
		//	elongation = box_w / box_h;
		//else
		//	elongation = box_h / box_w;
		////3、面积比，object_area/min-bounding rectangle_area,  convex or concave?
		//double are_rect = box_w*box_h;
		//double are_ratio = plrGraph.patches[combine_pindex].area / are_rect;
		//vector<double> fvec;//特征向量
		//fvec.push_back(circularity);//圆度
		//fvec.push_back(elongation);//长度
		//fvec.push_back(are_ratio);//面积
		//patch_tmpfeature[combine_pindex].second = fvec;

		mix_FDs[combine_pindex].clear(); mix_FDs[combine_pindex].swap(vector<double>());
		mix_FDs[combine_pindex] = tmppchFD;

		//3. pindex 的 status变为-1，且把其neighbor中的pindex换成combine_pindex，
		mixedTopology.allNode[pindex].status = -1;
		MixedArc* pindex_nei = mixedTopology.allNode[pindex].firstAdjNode;
		while (pindex_nei != NULL)
		{
			MixedArc* tmp_arc = getMixedArc(pindex_nei->adj_index, combine_pindex);
			if (tmp_arc != NULL)
			{
				removeMixedNode(pindex_nei->adj_index, pindex);
			}
			else
			{
				tmp_arc = getMixedArc(pindex_nei->adj_index, pindex);//////////////////////??????????????????
				tmp_arc->adj_index = combine_pindex;
			}
			pindex_nei = pindex_nei->next_arc;
		}//注意结果是：combine_pindex中的pindex改成了combine_pindex，pindex中还连着各个neighbor		
		//修改combine_pindex的 topology ，即清空combine_pindex的neighbor，重新添加neighbor
		vector<int> nei;
		MixedArc* o_arc = mixedTopology.allNode[pindex].firstAdjNode;//pindex 的neighbor删除了pindex，但是pindex还保留着这些neighbor的连接
		while (o_arc != NULL)
		{
			if (o_arc->adj_index != combine_pindex)//!!
				nei.push_back(o_arc->adj_index);
			o_arc = o_arc->next_arc;
		}
		o_arc = mixedTopology.allNode[combine_pindex].firstAdjNode;
		while (o_arc != NULL)
		{
			if (o_arc->adj_index != combine_pindex)//!!
				nei.push_back(o_arc->adj_index);
			o_arc = o_arc->next_arc;
		}
		deleRepeatPoint(&nei);
		mixedTopology.allNode[combine_pindex].firstAdjNode = NULL;////清空combine_pindex的neighbor，重新添加
		for (unsigned int ii = 0; ii < nei.size(); ++ii)
		{
			MixedArc* neiarc = new MixedArc;
			neiarc->adj_index = nei[ii];
			neiarc->weight = -1;
			neiarc->next_arc = mixedTopology.allNode[combine_pindex].firstAdjNode;
			mixedTopology.allNode[combine_pindex].firstAdjNode = neiarc;
			MixedArc* ensure_arc = getMixedArc(nei[ii], combine_pindex);
			if (ensure_arc == NULL)
			{
				MixedArc* neiarc_rever = new MixedArc;
				neiarc_rever->adj_index = combine_pindex;
				neiarc_rever->weight = -1;
				neiarc_rever->next_arc = mixedTopology.allNode[nei[ii]].firstAdjNode;
				mixedTopology.allNode[nei[ii]].firstAdjNode = neiarc_rever;
			}
		}
		//!!!
		//char pindex_ttype = mixedTopology.allNode[pindex].pl_type;
		//char combine_Pindex_ttype = mixedTopology.allNode[combine_pindex].pl_type;
		//if ((pindex_ttype == 'p' && combine_Pindex_ttype != 'p') || (pindex_ttype != 'p' && combine_Pindex_ttype == 'p'))
		//{
		//	removeMixedNode(pindex, combine_pindex);
		//	removeMixedNode(combine_pindex, pindex);
		//	MixedArc* pindex_nei = mixedTopology.allNode[pindex].firstAdjNode;
		//	while (pindex_nei != NULL)
		//	{
		//		removeMixedNode(pindex_nei->adj_index, pindex);
		//		pindex_nei = pindex_nei->next_arc;
		//	}
		//}
		//else
		//{
		//	MixedArc* pindex_nei = mixedTopology.allNode[pindex].firstAdjNode;
		//	while (pindex_nei != NULL)
		//	{
		//		MixedArc* tmp_arc = getMixedArc(pindex_nei->adj_index, combine_pindex);
		//		if (tmp_arc != NULL)
		//		{
		//			removeMixedNode(pindex_nei->adj_index, pindex);
		//		}
		//		else
		//		{
		//			tmp_arc = getMixedArc(pindex_nei->adj_index, pindex);//////////////////////??????????????????
		//			tmp_arc->adj_index = combine_pindex;
		//		}
		//		pindex_nei = pindex_nei->next_arc;
		//	}//注意结果是：combine_pindex中的pindex改成了combine_pindex，pindex中还连着各个neighbor		
		//	//修改combine_pindex的 topology ，即清空combine_pindex的neighbor，重新添加neighbor
		//	vector<int> nei;
		//	MixedArc* o_arc = mixedTopology.allNode[pindex].firstAdjNode;//pindex 的neighbor删除了pindex，但是pindex还保留着这些neighbor的连接
		//	while (o_arc != NULL)
		//	{
		//		if (o_arc->adj_index != combine_pindex)//!!
		//			nei.push_back(o_arc->adj_index);
		//		o_arc = o_arc->next_arc;
		//	}
		//	o_arc = mixedTopology.allNode[combine_pindex].firstAdjNode;
		//	while (o_arc != NULL)
		//	{
		//		if (o_arc->adj_index != combine_pindex)//!!
		//			nei.push_back(o_arc->adj_index);
		//		o_arc = o_arc->next_arc;
		//	}
		//	deleRepeatPoint(&nei);
		//	mixedTopology.allNode[combine_pindex].firstAdjNode = NULL;////清空combine_pindex的neighbor，重新添加
		//	for (unsigned int ii = 0; ii < nei.size(); ++ii)
		//	{
		//		MixedArc* neiarc = new MixedArc;
		//		neiarc->adj_index = nei[ii];
		//		neiarc->weight = 0;
		//		neiarc->next_arc = mixedTopology.allNode[combine_pindex].firstAdjNode;
		//		mixedTopology.allNode[combine_pindex].firstAdjNode = neiarc;
		//		MixedArc* ensure_arc = getMixedArc(nei[ii], combine_pindex);
		//		if (ensure_arc == NULL)
		//		{
		//			MixedArc* neiarc_rever = new MixedArc;
		//			neiarc_rever->adj_index = combine_pindex;
		//			neiarc_rever->weight = 0;
		//			neiarc_rever->next_arc = mixedTopology.allNode[nei[ii]].firstAdjNode;
		//			mixedTopology.allNode[nei[ii]].firstAdjNode = neiarc_rever;
		//		}
		//	}
		//}
		
		//4.重新计算combine_patch 与所有neighbor的相似度weight，tt_w[],更新在堆中的位置
		if (mixedTopology.allNode[combine_pindex].pl_type == 'p')
		{
			MixedArc* cmb_adjnode = mixedTopology.allNode[combine_pindex].firstAdjNode;
			double mini_simi = 6553500;
			int simi_pindex = -1;
			while (cmb_adjnode != NULL)
			{
				if (mixedTopology.allNode[cmb_adjnode->adj_index].pl_type !='p')  //p+l
				{
					vector<mypoint2f> modi_patch;
					bool pl_flag = combinePL_getplist(combine_pindex, cmb_adjnode->adj_index - pchsize, &modi_patch);
					int record_adjindex = cmb_adjnode->adj_index;   //记录一下adj index
					if (pl_flag == true)//是NULL的话说明不能合并
					{
						modi_patch.pop_back();  //front()！= back()	
						vector<int> empty_vec;
						vector<double> new_ft;
						getFDFeature(&modi_patch, &empty_vec, combine_pindex, simply_num, similarity_type, &new_ft);//从p+l得到的FD，注意第三个参数是原patch的index
						double dis1 = getEuclidDistance(&mix_FDs[combine_pindex], &new_ft);
						cmb_adjnode->weight = dis1;
						MixedArc* rever_ptcharc = getMixedArc(cmb_adjnode->adj_index, combine_pindex);
						rever_ptcharc->weight = dis1;
						if (dis1 < mini_simi)
						{
							mini_simi = dis1;
							simi_pindex = cmb_adjnode->adj_index;//
						}
						cmb_adjnode = cmb_adjnode->next_arc;
					}
					else
					{
						bool cun_flag = false; //看看patch里有没有这条线
						for (unsigned int k = 0; k < plrGraph.patches[combine_pindex].point_line.size(); ++k)
						{
							if (plrGraph.patches[combine_pindex].point_line[k].second == cmb_adjnode->adj_index - pchsize)
							{
								cun_flag = true;
								break;
							}
						}
						if (cun_flag == true)//有这条线，但是不能p+l,赋最大值，
						{
							cmb_adjnode->weight = 6553500;
							MixedArc* rever_ptcharc = getMixedArc(cmb_adjnode->adj_index, combine_pindex);
							rever_ptcharc->weight = 6553500;
							cmb_adjnode = cmb_adjnode->next_arc;
						}
						else
						{
							int tmp_dele = cmb_adjnode->adj_index;
							cmb_adjnode = cmb_adjnode->next_arc;
							removeMixedNode(tmp_dele, combine_pindex);
							removeMixedNode(combine_pindex, tmp_dele);
						}
						/*cmb_adjnode->weight = 6553500;
						MixedArc* rever_ptcharc = getMixedArc(cmb_adjnode->adj_index, combine_pindex);
						rever_ptcharc->weight = 6553500;
						cmb_adjnode = cmb_adjnode->next_arc;*/
					}
					//重新选cmb_adjnode->adj_patchindex的最小weight
					double arc_minweight = 6553500;
					int arc_minipindex = -1;
					MixedArc* arc_neighbor = mixedTopology.allNode[record_adjindex].firstAdjNode;
					while (arc_neighbor != NULL)
					{
						if (arc_neighbor->weight < arc_minweight && arc_neighbor->weight >= 0)
						{
							arc_minweight = arc_neighbor->weight;
							arc_minipindex = arc_neighbor->adj_index;
						}
						arc_neighbor = arc_neighbor->next_arc;
					}
					//跟新
					//if (tt_w[record_adjindex] != arc_minweight || weight_pindex[record_adjindex] != arc_minipindex)
					//{
						heap_update(ele_index[record_adjindex] + 1, arc_minweight, &heap_vec, &tt_w, &ele_index);
						weight_pindex[record_adjindex] = arc_minipindex;
					//}	
				}
				else  //p+p
				{
					vector<int> combine_pvec;//合并的pointlist
					combineTwoPatch_getpvec(combine_pindex, cmb_adjnode->adj_index, &combine_pvec);//combine_plist  front()！= back()
					int record_adjindex = cmb_adjnode->adj_index;   //记录一下adj index
					if (!combine_pvec.empty())
					{
						vector<mypoint2f> combine_plist;//合并的pointlist
						jointMethod(&combine_pvec, &combine_plist, 0);

						//vector<double> threeft;//////////////////////////////////////
						//getThreeFeature(&threeft, &combine_plist, 0, 0);//////////////////////////////////

						vector<double> newpchFD;//FD
						combine_plist.pop_back();  //front()！= back()	
						getFDFeature(&combine_plist, &combine_pvec, -1, simply_num, similarity_type, &newpchFD);
						double dis1 = getEuclidDistance(&newpchFD, &mix_FDs[combine_pindex]);
						double dis2 = getEuclidDistance(&newpchFD, &mix_FDs[cmb_adjnode->adj_index]);
						//double dis1 = getEuclidDistance(&newpchFD, &mix_FDs[combine_pindex])*0.4 + getEuclidDistance(&threeft,&mix_feature[combine_pindex])*0.6;///////////////////////////////////////////////
						//double dis2 = getEuclidDistance(&newpchFD, &mix_FDs[cmb_adjnode->adj_index])*0.4 + getEuclidDistance(&threeft, &mix_feature[cmb_adjnode->adj_index])*0.6;//////////////////////////////////


						double mini_dis = dis1 < dis2 ? dis1 : dis2;
						cmb_adjnode->weight = mini_dis;
						MixedArc* rever_ptcharc = getMixedArc(cmb_adjnode->adj_index, combine_pindex);
						rever_ptcharc->weight = mini_dis;
						if (mini_dis < mini_simi)
						{
							mini_simi = mini_dis;
							simi_pindex = cmb_adjnode->adj_index;//
						}
						//vector<double> ft_vec;//3ge特征
						//getThreeFeature(&ft_vec, &combine_plist, 0, 0);
						//double dis1 = getEuclidDistance(&ft_vec, &mix_FDs[combine_pindex]);
						//parc->weight = dis1;
						//if (dis1 < mini_simi)
						//{
						//	mini_simi = dis1;
						//	simi_pindex = parc->adj_index;//
						//}
						//double dis2 = getEuclidDistance(&ft_vec, &mix_FDs[parc->adj_index]);
						//PatchArc* rever_ptcharc = getPatchArc(parc->adj_index, combine_pindex);
						//rever_ptcharc->weight = dis2;

						cmb_adjnode = cmb_adjnode->next_arc;
					}
					else
					{
						int tmp_dele = cmb_adjnode->adj_index;
						cmb_adjnode = cmb_adjnode->next_arc;
						removeMixedNode(tmp_dele, combine_pindex);
						removeMixedNode(combine_pindex, tmp_dele);
					}
					//重新选cmb_adjnode->adj_patchindex的最小weight, 
					//cmb_adjnode已经改变，用record_adjindex来代替cmb_adjnode->adj_patchindex!!!
					double arc_minweight = 6553500;
					int arc_minipindex = -1;
					MixedArc* arc_neighbor = mixedTopology.allNode[record_adjindex].firstAdjNode;
					while (arc_neighbor != NULL)
					{
						if (arc_neighbor->weight < arc_minweight && arc_neighbor->weight >= 0)
						{
							arc_minweight = arc_neighbor->weight;
							arc_minipindex = arc_neighbor->adj_index;
						}
						arc_neighbor = arc_neighbor->next_arc;
					}
					//if (tt_w[record_adjindex] != arc_minweight || weight_pindex[record_adjindex] != arc_minipindex)
					//{//跟新
						heap_update(ele_index[record_adjindex] + 1, arc_minweight, &heap_vec, &tt_w, &ele_index);
						weight_pindex[record_adjindex] = arc_minipindex;
					//}
				}
				//cmb_adjnode = cmb_adjnode->nextPatch;
			}
			heap_update(ele_index[combine_pindex] + 1, mini_simi, &heap_vec, &tt_w, &ele_index);//注意+1	tt_w[combine_pindex] = dis1;	
			weight_pindex[combine_pindex] = simi_pindex;
		}
		else  //combine_pindex 是一个line
		{
			MixedArc* cmb_adjnode = mixedTopology.allNode[combine_pindex].firstAdjNode;
			double mini_simi = 6553500;
			int simi_pindex = -1;
			while (cmb_adjnode != NULL)
			{
				int record_adjindex = cmb_adjnode->adj_index;   //记录一下adj index
				if (mixedTopology.allNode[cmb_adjnode->adj_index].pl_type == 'p')  //l+p
				{
					vector<mypoint2f> modi_patch;
					bool pl_flag = combinePL_getplist(cmb_adjnode->adj_index, combine_pindex - pchsize, &modi_patch);
					if (pl_flag == true)//是NULL的话说明不能合并
					{
						modi_patch.pop_back();  //front()！= back()	
						vector<int> empty_vec;
						vector<double> new_ft;
						getFDFeature(&modi_patch, &empty_vec, cmb_adjnode->adj_index, simply_num, similarity_type, &new_ft);//从p+l得到的FD，注意第三个参数是原patch的index
						double dis1 = getEuclidDistance(&mix_FDs[combine_pindex], &new_ft);
						cmb_adjnode->weight = dis1;
						MixedArc* rever_ptcharc = getMixedArc(cmb_adjnode->adj_index, combine_pindex);
						rever_ptcharc->weight = dis1;
						if (dis1 < mini_simi)
						{
							mini_simi = dis1;
							simi_pindex = cmb_adjnode->adj_index;//
						}
						cmb_adjnode = cmb_adjnode->next_arc;
					}
					else  
					{
						bool cun_flag = false; //看看patch里有没有这条线
						for (unsigned int k = 0; k < plrGraph.patches[cmb_adjnode->adj_index].point_line.size(); ++k)
						{
							if (plrGraph.patches[cmb_adjnode->adj_index].point_line[k].second == combine_pindex - pchsize)
							{
								cun_flag = true;
								break;
							}
						}
						if (cun_flag == true)
						{
							cmb_adjnode->weight = 6553500;
							MixedArc* rever_ptcharc = getMixedArc(cmb_adjnode->adj_index, combine_pindex);
							rever_ptcharc->weight = 6553500;
							cmb_adjnode = cmb_adjnode->next_arc;
						}
						else
						{
							int tmp_dele = cmb_adjnode->adj_index;
							cmb_adjnode = cmb_adjnode->next_arc;
							removeMixedNode(tmp_dele, combine_pindex);
							removeMixedNode(combine_pindex, tmp_dele);
						}
						/*cmb_adjnode->weight = 6553500;
						MixedArc* rever_ptcharc = getMixedArc(cmb_adjnode->adj_index, combine_pindex);
						rever_ptcharc->weight = 6553500;
						cmb_adjnode = cmb_adjnode->next_arc;*/
					}
				}
				else  //l+l
				{
					vector<mypoint2f> newplist;
					bool lsucflag = combineTwoCurves_getplist(combine_pindex - pchsize, cmb_adjnode->adj_index - pchsize, &newplist, patch_averWidth);
					if (lsucflag == true)
					{
						vector<mypoint2f> newpatch_plist;
						newpatch_plist.insert(newpatch_plist.end(), newplist.begin(), newplist.end());
						for (int ppi = (int)newplist.size() - 2; ppi > 0; --ppi) //newpatch_plist.front != back
						{
							newpatch_plist.push_back(newplist[ppi]);
						}
						vector<int> line_pvec;//line的patch的pvec,没有是空的
						vector<double> new_FD; //line 的FD
						getFDFeature(&newpatch_plist, &line_pvec,-1, simply_num, similarity_type, &new_FD);//参数：第0种FD计算方法，或者第1种
						double dis1 = getEuclidDistance(&mix_FDs[combine_pindex], &new_FD);
						double dis2 = getEuclidDistance(&mix_FDs[cmb_adjnode->adj_index], &new_FD);

						cmb_adjnode->weight = dis1 < dis2 ? dis1 : dis2;//取较小的那个值
						MixedArc* parc_reve = getMixedArc(cmb_adjnode->adj_index, combine_pindex);//patchTopology.allpatches[cmb_adjnode->adj_patchindex].firstAdjPatch;
						if (parc_reve != NULL)
							parc_reve->weight = cmb_adjnode->weight;
						else
							cout << "patch arc=NULL" << endl;
						if (cmb_adjnode->weight < mini_simi)
						{
							mini_simi = cmb_adjnode->weight;
							simi_pindex = cmb_adjnode->adj_index;
						}
						cmb_adjnode = cmb_adjnode->next_arc;
					}
					else   //两条line 不能merge
					{ 
						int tmp_dele = cmb_adjnode->adj_index;
						cmb_adjnode = cmb_adjnode->next_arc;
						removeMixedNode(tmp_dele, combine_pindex);
						removeMixedNode(combine_pindex, tmp_dele);		
					}			
				}//l+l done	
				//重新选cmb_adjnode->adj_patchindex的最小weight  ,,,
				double arc_minweight = 6553500;
				int arc_minipindex = -1;
				MixedArc* arc_neighbor = mixedTopology.allNode[record_adjindex].firstAdjNode;
				while (arc_neighbor != NULL)
				{
					if (arc_neighbor->weight < arc_minweight && arc_neighbor->weight >= 0)
					{
						arc_minweight = arc_neighbor->weight;
						arc_minipindex = arc_neighbor->adj_index;
					}
					arc_neighbor = arc_neighbor->next_arc;
				}
				//if (tt_w[record_adjindex] != arc_minweight || weight_pindex[record_adjindex] != arc_minipindex)
				//{
					//combine_pindex 的neighbor 的更新
					heap_update(ele_index[record_adjindex] + 1, arc_minweight, &heap_vec, &tt_w, &ele_index);
					weight_pindex[record_adjindex] = arc_minipindex;
				//}
			}//while循环	
			//combine_pindex的更新
			heap_update(ele_index[combine_pindex] + 1, mini_simi, &heap_vec, &tt_w, &ele_index);//注意+1	tt_w[combine_pindex] = dis1;	
			weight_pindex[combine_pindex] = simi_pindex;
		}
	}
	cout << "mix process done" << endl; 
	cout << "line merge operation:" << dele_count_line << endl;
	cout << "patch merge operation:" << dele_count_patch << endl;
	cout << "patch+line merge operation:" << dele_count_pl << endl;
	cout << endl;
}
void MyGraph::test_MixSimilarityRange(int simply_num, double patch_averWidth)
{
	cout << "MixSimilarityRange" << endl;
	/*1、计算 patch和line 的FD*/
	int pchsize = plrGraph.patches.size();
	int attachsize = plrGraph.attachlines.size();
	int isolatedsize = plrGraph.isolatedlines.size();
	int linesize = attachsize + isolatedsize;   //即otherline的大小
	int mixsize = pchsize + attachsize + isolatedsize;

	vector<vector<double>> mix_feature;//顺序 与mixedtopology相同，即先patch，attachline，isolatedline
	vector<vector<double>> mix_FDs;
	for (int i = 0; i < mixsize; ++i)
	{
		vector<double> tmp;
		mix_feature.push_back(tmp);     //赋初值 三个feature
		mix_FDs.push_back(tmp);          //赋初值 FD
		if (mixedTopology.allNode[i].status >= 0)
		{
			if (mixedTopology.allNode[i].pl_type == 'p')
			{
				vector<mypoint2f> total_plist;
				jointMethod(&plrGraph.patches[i].pointIndexOfPatch, &total_plist, 0);
				total_plist.pop_back();  //!!!	
				vector<double> pchFD; //patch 的FD
				getFDFeature(&total_plist, &plrGraph.patches[i].pointIndexOfPatch,i, simply_num,1, &pchFD);//该函数内部包含了plist的化简，参数：第0种FD计算方法，或者第1种
				mix_FDs[i] = pchFD;
				//vector<double> fvec;//特征向量
				//getThreeFeature(&fvec, &total_plist, plrGraph.patches[i].area, plrGraph.patches[i].perimeter);
				//mix_feature[i] = fvec;
			}
			//else
			//{
			//	vector<mypoint2f> attach_plist;
			//	jointMethod(&plrGraph.otherlines[i - pchsize].pt_vec, &attach_plist, 0);  //attach + isolate
			//	vector<mypoint2f> patch_plist;
			//	patch_plist.insert(patch_plist.end(), attach_plist.begin(), attach_plist.end());
			//	for (int ppi = (int)attach_plist.size() - 2; ppi >0; --ppi)//patch_plist . front！=back
			//	{
			//		patch_plist.push_back(attach_plist[ppi]);
			//	}
			//	vector<int> line_pvec;
			//	line_pvec = plrGraph.otherlines[i - pchsize].pt_vec;
			//	for (int j = plrGraph.otherlines[i - pchsize].pt_vec.size() - 2; j >= 0; --j)
			//	{
			//		line_pvec.push_back(plrGraph.otherlines[i - pchsize].pt_vec[j]);
			//	}
			//	vector<double> l_FD; //line 的FD
			//	getFDFeature(&patch_plist, &line_pvec, -1, simply_num, 0, &l_FD);//该函数内部包含了plist的化简，参数：第0种FD计算方法，或者第1种
			//	mix_FDs[i] = (l_FD);
			//}
		}
	}
	///*2、计算mixtopology.allnode的weight。分开计算: patch与patch合并，patch与line合并，line和line合并*/
	////计算每个patch与neighbor合并后的特征值, dis(center_ft,new_ft)->patchTopology.allpatches[i].arc.weight中 (//patch_subgraphs.push_back(make_pair(i, "tree"));  )
	//ofstream ps_file1("E:/svgfile/p+l/p_similarity_dis.txt");
	//ps_file1.is_open();
	//vector<double> pdis;
	//ofstream ls_file2("E:/svgfile/p+l/l_similarity_dis.txt");
	//ls_file2.is_open();
	//vector<double> ldis;

	//vector<double> tt_w;//记录权值
	//vector<int> weight_pindex;//记录某个patch的最小weight对应的patchindex
	//for (int i = 0; i < mixsize; ++i)  //mixedTopology 从0->pchsize
	//{
	//	tt_w.push_back(-1);
	//	weight_pindex.push_back(-1);
	//	if (mixedTopology.allNode[i].status < 0)	//if (plrGraph.patches[i].type >= 0)
	//	{
	//		continue;
	//	}
	//	if (mixedTopology.allNode[i].pl_type == 'p')  // 1、这个点是patch
	//	{
	//		MixedArc* parc = mixedTopology.allNode[i].firstAdjNode;
	//		double min_dis = 6553500;
	//		int min_pindex = -1;
	//		while (parc != NULL)//每一个arc 是一个neighbor patch, 中心patch与neighborpatch合并，然后计算合并后的feature，然后计算该newfeature与patch_tmpfeature[i].second的二范数，存入arc 的weight
	//		{
	//			if (parc->weight > 0)
	//			{
	//				if (parc->weight < min_dis)
	//				{
	//					min_dis = parc->weight;
	//					min_pindex = parc->adj_index;
	//				}
	//				parc = parc->next_arc;
	//				continue;
	//			}
	//			char nodeType = mixedTopology.allNode[parc->adj_index].pl_type;//判断neighbor 的类型，
	//			if (nodeType == 'p')	    /*1.1、 p p 合并*/
	//			{
	//				vector<int> combine_pvec;//合并的pointlist
	//				combineTwoPatch_getpvec(i, parc->adj_index, &combine_pvec);//combine_plist  front()！= back()
	//				if (!combine_pvec.empty())//是NULL的话说明不能合并
	//				{
	//					vector<mypoint2f> combine_plist;//合并的pointlist
	//					jointMethod(&combine_pvec, &combine_plist, 0);
	//					//计算新patch的周长+FD
	//					vector<double> new_ft;
	//					combine_plist.pop_back();  //front()！= back()	
	//					getFDFeature(&combine_plist, &combine_pvec,-1, simply_num, 0, &new_ft);
	//					double dis1 = getEuclidDistance(&mix_FDs[i], &new_ft);
	//					double dis2 = getEuclidDistance(&mix_FDs[parc->adj_index], &new_ft);
	//					if (parc->weight == 0)//parc->weight初始化是赋值为0，在选取最小值的注意把0排除!
	//					{
	//						parc->weight = dis1 < dis2 ? dis1 : dis2;//取较小的那个值
	//					}
	//					else
	//					{
	//						cout << "??" << endl;//不可能有这种情况
	//					}
	//					MixedArc* parc_reve = getMixedArc(parc->adj_index, i);//patchTopology.allpatches[parc->adj_patchindex].firstAdjPatch;
	//					if (parc_reve != NULL)
	//						parc_reve->weight = parc->weight;
	//					else
	//						cout << "patch arc=NULL" << endl;
	//					if (parc->weight < min_dis)
	//					{
	//						min_dis = parc->weight;
	//						min_pindex = parc->adj_index;
	//					}
	//					pdis.push_back(parc->weight);  //输出
	//					parc = parc->next_arc;
	//				}
	//				else
	//				{
	//					//虽然两patch相邻，但是不能merge，从neighbor 中删除
	//					int tmp_nindex = parc->adj_index;
	//					parc = parc->next_arc;
	//					removeMixedNode(i, tmp_nindex);
	//					removeMixedNode(tmp_nindex, i);
	//				}
	//			}
	//			else if (nodeType == 'l')		/* 1.2、p l 合并  how?????????? */
	//			{
	//				parc->weight = 6553500;   //还没有考虑好 patch和line怎么merge，先将arc的weight赋一个最大值，，，
	//				getMixedArc(parc->adj_index, i)->weight = 6553500;
	//				min_pindex = parc->adj_index;
	//				parc = parc->next_arc;
	//			}
	//			//parc = parc->next_arc;
	//		} //while done
	//		if (min_dis == 6553500)//patch 没有neighbor，或者有neighbor但是neighbor是line
	//		{
	//		}
	//		tt_w[i] = min_dis;
	//		weight_pindex[i] = min_pindex;
	//	}
	//	else if (mixedTopology.allNode[i].pl_type == 'l' || mixedTopology.allNode[i].pl_type == 'i')  // 2、这个点是attachline
	//	{
	//		MixedArc* parc = mixedTopology.allNode[i].firstAdjNode;
	//		double min_dis = 6553500;
	//		int min_pindex = -1;
	//		while (parc != NULL)//每一个arc 是一个neighbor patch, 中心patch与neighborpatch合并，然后计算合并后的feature，然后计算该newfeature与patch_tmpfeature[i].second的二范数，存入arc 的weight
	//		{
	//			if (parc->weight > 0)
	//			{
	//				if (parc->weight < min_dis)
	//				{
	//					min_dis = parc->weight;
	//					min_pindex = parc->adj_index;
	//				}
	//				parc = parc->next_arc;
	//				continue;
	//			}
	//			char nodeType = mixedTopology.allNode[parc->adj_index].pl_type;//判断neighbor 的类型，
	//			if (nodeType == 'p')		//2.1  l+p
	//			{
	//				parc->weight = 6553500;   //还没有考虑好 patch和line怎么merge，先将arc的weight赋一个最大值，，，
	//				getMixedArc(parc->adj_index, i)->weight = 6553500;
	//				min_pindex = parc->adj_index;
	//			}
	//			else	//2.2  l+l /l+i
	//			{
	//				vector<mypoint2f> newplist;
	//				bool c_flag = combineTwoCurves_getplist(i - pchsize, parc->adj_index - pchsize, &newplist, patch_averWidth);

	//				vector<mypoint2f> newpatch_plist;
	//				newpatch_plist.insert(newpatch_plist.end(), newplist.begin(), newplist.end());
	//				for (int ppi = (int)newplist.size() - 2; ppi >0; --ppi) //front！=back
	//				{
	//					newpatch_plist.push_back(newplist[ppi]);
	//				}
	//				vector<int> line_pvec;//line的patch的pvec,没有是空的
	//				vector<double> new_FD; //line 的FD
	//				getFDFeature(&newpatch_plist, &line_pvec, -1, simply_num, 0, &new_FD); //参数：第0种FD计算方法，或者第1种
	//				double dis1 = getEuclidDistance(&mix_FDs[i], &new_FD);
	//				double dis2 = getEuclidDistance(&mix_FDs[parc->adj_index], &new_FD);

	//				parc->weight = dis1 < dis2 ? dis1 : dis2;//取较小的那个值
	//				MixedArc* parc_reve = getMixedArc(parc->adj_index, i);//patchTopology.allpatches[parc->adj_patchindex].firstAdjPatch;
	//				if (parc_reve != NULL)
	//					parc_reve->weight = parc->weight;
	//				else
	//					cout << "patch arc=NULL" << endl;
	//				if (parc->weight < min_dis)
	//				{
	//					min_dis = parc->weight;
	//					min_pindex = parc->adj_index;
	//				}
	//				ldis.push_back(parc->weight);//输出
	//			}
	//			parc = parc->next_arc;
	//		}
	//		if (min_dis == 6553500)//line 没有neighbor，或者有neighbor但是 是patch
	//		{
	//		}
	//		tt_w[i] = min_dis;
	//		weight_pindex[i] = min_pindex;
	//	}
	//}
	//vector<double>::iterator iter_tp = pdis.begin();
	//while (iter_tp != pdis.end())
	//{
	//	vector<double>::iterator iter_tp_inn = iter_tp + 1;
	//	while (iter_tp_inn != pdis.end())
	//	{
	//		if (*iter_tp == *iter_tp_inn)
	//			iter_tp_inn = pdis.erase(iter_tp_inn);
	//		else
	//			iter_tp_inn++;
	//	}
	//	iter_tp++;
	//}
	//iter_tp = ldis.begin();
	//while (iter_tp != ldis.end())
	//{
	//	vector<double>::iterator iter_tp_inn = iter_tp + 1;
	//	while (iter_tp_inn != ldis.end())
	//	{
	//		if (*iter_tp == *iter_tp_inn)
	//			iter_tp_inn = ldis.erase(iter_tp_inn);
	//		else
	//			iter_tp_inn++;
	//	}
	//	iter_tp++;
	//}
	//sort(pdis.begin(), pdis.end());
	//sort(ldis.begin(), ldis.end());
	//cout << "patch的mini distance个数：" << pdis.size() << endl;
	//cout << "line的mini distance个数：" << ldis.size() << endl;
	//for (unsigned int i = 0; i < pdis.size(); ++i)
	//{
	//	ps_file1 << pdis[i] << " ";
	//}
	//for (unsigned int i = 0; i < ldis.size(); ++i)
	//{
	//	ls_file2 << ldis[i] << " ";
	//}
}
void MyGraph::origin_MixSimilarityRange(int simply_num, double patch_averWidth, int similarity_type, double similarity_t)
{
	/*1、计算 patch和line 的FD*/
	int pchsize = plrGraph.patches.size();
	int attachsize = plrGraph.attachlines.size();
	int isolatedsize = plrGraph.isolatedlines.size();
	int linesize = attachsize + isolatedsize;   //即otherline的大小
	int mixsize = pchsize + attachsize + isolatedsize;

	vector<vector<double>> mix_feature;//顺序 与mixedtopology相同，即先patch，attachline，isolatedline
	vector<vector<double>> mix_FDs;
	for (int i = 0; i < mixsize; ++i)
	{
		vector<double> tmp;
		mix_feature.push_back(tmp);     //赋初值 三个feature
		mix_FDs.push_back(tmp);          //赋初值 FD
		if (mixedTopology.allNode[i].status >= 0)
		{
			if (mixedTopology.allNode[i].pl_type == 'p')
			{
				vector<mypoint2f> total_plist;
				jointMethod(&plrGraph.patches[i].pointIndexOfPatch, &total_plist, 0);
				total_plist.pop_back();  //!!!	
				vector<double> pchFD; //patch 的FD
				getFDFeature(&total_plist, &plrGraph.patches[i].pointIndexOfPatch, i, simply_num, similarity_type, &pchFD);//该函数内部包含了plist的化简，参数：第0种FD计算方法，或者第1种
				mix_FDs[i] = pchFD;
				//vector<double> fvec;//特征向量
				//getThreeFeature(&fvec, &total_plist, plrGraph.patches[i].area, plrGraph.patches[i].perimeter);
				//mix_feature[i] = fvec;
			}
			else
			{
				vector<mypoint2f> attach_plist;
				jointMethod(&plrGraph.otherlines[i - pchsize].pt_vec, &attach_plist, 0);  //attach + isolate
				vector<mypoint2f> patch_plist;//point list围成patch
				patch_plist.insert(patch_plist.end(), attach_plist.begin(), attach_plist.end());
				for (int ppi = (int)attach_plist.size() - 2; ppi >0; --ppi)//patch_plist 的 front！=back
				{
					patch_plist.push_back(attach_plist[ppi]);
				}
				vector<int> line_pvec;//point vec围成patch
				line_pvec = plrGraph.otherlines[i - pchsize].pt_vec;
				for (int j = plrGraph.otherlines[i - pchsize].pt_vec.size() - 2; j >= 0; --j)//line_pvec的front = back
				{
					line_pvec.push_back(plrGraph.otherlines[i - pchsize].pt_vec[j]);
				}
				vector<double> l_FD; //line 的FD
				getFDFeature(&patch_plist, &line_pvec, -1, simply_num, similarity_type, &l_FD);//该函数内部包含了plist的化简，参数：第0种FD计算方法，或者第1种
				mix_FDs[i] = (l_FD);
			}
		}
	}

	/*2、计算mixtopology.allnode的weight。分开计算: patch与patch合并，patch与line合并，line和line合并*/
	//计算每个patch与neighbor合并后的特征值, dis(center_ft,new_ft)->patchTopology.allpatches[i].arc.weight中 (//patch_subgraphs.push_back(make_pair(i, "tree"));  )
	ofstream ps_file1("E:/svgfile/p+l/p_similarity_dis.txt");
	ps_file1.is_open();
	vector<double> pdis;
	ofstream ls_file2("E:/svgfile/p+l/l_similarity_dis.txt");
	ls_file2.is_open();
	vector<double> ldis;
	ofstream pl_file3("E:/svgfile/p+l/p+l_similarity_dis.txt");
	pl_file3.is_open();
	vector<double> pl_dis;

	vector<double> tt_w;//记录权值
	vector<int> weight_pindex;//记录某个patch的最小weight对应的patchindex
	for (int i = 0; i < mixsize; ++i)  //mixedTopology 从0->pchsize
	{
		tt_w.push_back(-1);
		weight_pindex.push_back(-1);
		if (mixedTopology.allNode[i].status < 0)	//if (plrGraph.patches[i].type >= 0)
		{
			continue;
		}
		if (mixedTopology.allNode[i].pl_type == 'p')  // 1、这个点是patch
		{
			MixedArc* parc = mixedTopology.allNode[i].firstAdjNode;
			double min_dis = 6553500;
			int min_pindex = -1;
			while (parc != NULL)//每一个arc 是一个neighbor patch, 中心patch与neighborpatch合并，然后计算合并后的feature，然后计算该newfeature与patch_tmpfeature[i].second的二范数，存入arc 的weight
			{
				if (parc->weight > 0)
				{
					if (parc->weight < min_dis)
					{
						min_dis = parc->weight;
						min_pindex = parc->adj_index;
					}
					parc = parc->next_arc;
					continue;
				}
				char nodeType = mixedTopology.allNode[parc->adj_index].pl_type;//判断neighbor 的类型，
				if (nodeType == 'p')	    /*1.1、 p p 合并*/
				{
					vector<int> combine_pvec;//合并的pointlist
					combineTwoPatch_getpvec(i, parc->adj_index, &combine_pvec);//combine_plist  front()！= back()	
					if (!combine_pvec.empty())//是NULL的话说明不能合并
					{
						vector<mypoint2f> combine_plist;
						jointMethod(&combine_pvec, &combine_plist, 0);
						//计算新patch的周长+FD
						vector<double> new_ft;
						combine_plist.pop_back();  //front()！= back()	
						getFDFeature(&combine_plist, &combine_pvec, -1, simply_num, similarity_type, &new_ft);
						double dis1 = getEuclidDistance(&mix_FDs[i], &new_ft);
						double dis2 = getEuclidDistance(&mix_FDs[parc->adj_index], &new_ft);
						if (parc->weight == 0)//parc->weight初始化是赋值为0，在选取最小值的注意把0排除!
						{
							parc->weight = dis1 < dis2 ? dis1 : dis2;//取较小的那个值
						}
						else
						{
							////取parc->weight,dis1,dis2三者的最小值
							//if (dis1 < parc->weight)
							//{
							//	if (dis2 < dis1)
							//		parc->weight = dis2;
							//	else
							//		parc->weight = dis1;
							//}
							//else 
							//{
							//	if (dis2 < parc->weight)
							//		parc->weight = dis2;
							//}
							cout << "??" << endl;//不可能有这种情况
						}
						MixedArc* parc_reve = getMixedArc(parc->adj_index, i);//patchTopology.allpatches[parc->adj_patchindex].firstAdjPatch;
						if (parc_reve != NULL)
							parc_reve->weight = parc->weight;
						else
							cout << "patch arc=NULL" << endl;
						if (parc->weight < min_dis)
						{
							min_dis = parc->weight;
							min_pindex = parc->adj_index;
						}
						pdis.push_back(parc->weight);  //输出
						parc = parc->next_arc;
					}
					else
					{
						//虽然两patch相邻，但是不能merge，从neighbor 中删除
						int tmp_nindex = parc->adj_index;
						parc = parc->next_arc;
						removeMixedNode(i, tmp_nindex);
						removeMixedNode(tmp_nindex, i);
					}
				}
				else if (nodeType == 'l')		/* 1.2、p l 合并  how?????????? */
				{
					vector<mypoint2f> modi_patch;
					bool pl_flag=combinePL_getplist(i, parc->adj_index, &modi_patch);
					if (pl_flag==true)//是NULL的话说明不能合并
					{
						modi_patch.pop_back();  //front()！= back()	
						vector<int> empty_vec;
						vector<double> new_ft;
						getFDFeature(&modi_patch, &empty_vec,-1, simply_num, similarity_type, &new_ft);
						double dis1 = getEuclidDistance(&mix_FDs[i], &new_ft);
						if (parc->weight<dis1)//parc->weight初始化是赋值为0，在选取最小值的注意把0排除!
						{
							parc->weight = dis1 ;//取较小的那个值
						}
						MixedArc* parc_reve = getMixedArc(parc->adj_index, i);//patchTopology.allpatches[parc->adj_patchindex].firstAdjPatch;
						if (parc_reve != NULL)
							parc_reve->weight = parc->weight;
						if (parc->weight < min_dis)
						{
							min_dis = parc->weight;
							min_pindex = parc->adj_index;
						}
						pl_dis.push_back(parc->weight);  //输出patch与line的merge后的dis
						parc = parc->next_arc;
					}
					else
					{
						//虽然patch 和line相邻，但是不能merge，从neighbor 中删除
						int tmp_nindex = parc->adj_index;
						parc = parc->next_arc;
						removeMixedNode(i, tmp_nindex);
						removeMixedNode(tmp_nindex, i);
					}
				}
				//parc = parc->next_arc;
			} //while done
			if (min_dis == 6553500)//patch 没有neighbor，或者有neighbor但是neighbor是line
			{
			}
			tt_w[i] = min_dis;
			weight_pindex[i] = min_pindex;
		}
		else if (mixedTopology.allNode[i].pl_type == 'l' || mixedTopology.allNode[i].pl_type == 'i')  // 2、这个点是attachline
		{
			MixedArc* parc = mixedTopology.allNode[i].firstAdjNode;
			double min_dis = 6553500;
			int min_pindex = -1;
			while (parc != NULL)//每一个arc 是一个neighbor patch, 中心patch与neighborpatch合并，然后计算合并后的feature，然后计算该newfeature与patch_tmpfeature[i].second的二范数，存入arc 的weight
			{
				if (parc->weight > 0)
				{
					if (parc->weight < min_dis)
					{
						min_dis = parc->weight;
						min_pindex = parc->adj_index;
					}
					parc = parc->next_arc;
					continue;
				}
				/*从这以下还没有改*/
				char nodeType = mixedTopology.allNode[parc->adj_index].pl_type;//判断neighbor 的类型，
				if (nodeType == 'p')		//2.1  l+p
				{
					parc->weight = 6553500;   //还没有考虑好 patch和line怎么merge，先将arc的weight赋一个最大值，，，
					getMixedArc(parc->adj_index, i)->weight = 6553500;
					min_pindex = parc->adj_index;
				}
				else	//2.2  l+l /l+i
				{
					vector<mypoint2f> newplist;
					bool c_flag = combineTwoCurves_getplist(i - pchsize, parc->adj_index - pchsize, &newplist, patch_averWidth);

					vector<mypoint2f> newpatch_plist;
					newpatch_plist.insert(newpatch_plist.end(), newplist.begin(), newplist.end());
					for (int ppi = (int)newplist.size() - 2; ppi >0; --ppi) //front！=back
					{
						newpatch_plist.push_back(newplist[ppi]);
					}
					vector<int> line_pvec;//line的patch的pvec,没有是空的
					vector<double> new_FD; //line 的FD
					getFDFeature(&newpatch_plist, &line_pvec, -1, simply_num, similarity_type, &new_FD); //参数：第0种FD计算方法，或者第1种
					double dis1 = getEuclidDistance(&mix_FDs[i], &new_FD);
					double dis2 = getEuclidDistance(&mix_FDs[parc->adj_index], &new_FD);

					parc->weight = dis1 < dis2 ? dis1 : dis2;//取较小的那个值
					MixedArc* parc_reve = getMixedArc(parc->adj_index, i);//patchTopology.allpatches[parc->adj_patchindex].firstAdjPatch;
					if (parc_reve != NULL)
						parc_reve->weight = parc->weight;
					else
						cout << "patch arc=NULL" << endl;
					if (parc->weight < min_dis)
					{
						min_dis = parc->weight;
						min_pindex = parc->adj_index;
					}
					ldis.push_back(parc->weight);  //输出
				}
				parc = parc->next_arc;
			}
			if (min_dis == 6553500)//line 没有neighbor，或者有neighbor但是 是patch
			{
			}
			tt_w[i] = min_dis;
			weight_pindex[i] = min_pindex;
		}
	}
	//输出patch 和line 的similaritydistance
	vector<double>::iterator iter_tp = pdis.begin();
	while (iter_tp != pdis.end())
	{
		vector<double>::iterator iter_tp_inn = iter_tp + 1;
		while (iter_tp_inn != pdis.end())
		{
			if (*iter_tp == *iter_tp_inn)
				iter_tp_inn = pdis.erase(iter_tp_inn);
			else
				iter_tp_inn++;
		}
		iter_tp++;
	}
	iter_tp = ldis.begin();
	while (iter_tp != ldis.end())
	{
		vector<double>::iterator iter_tp_inn = iter_tp + 1;
		while (iter_tp_inn != ldis.end())
		{
			if (*iter_tp == *iter_tp_inn)
				iter_tp_inn = ldis.erase(iter_tp_inn);
			else
				iter_tp_inn++;
		}
		iter_tp++;
	}
	sort(pdis.begin(), pdis.end());
	sort(ldis.begin(), ldis.end());
	cout << "patch的mini distance个数：" << pdis.size() << endl;
	cout << "line的mini distance个数：" << ldis.size() << endl;
	for (unsigned int i = 0; i < pdis.size(); ++i)
	{
		ps_file1 << pdis[i] << " ";
	}
	for (unsigned int i = 0; i < ldis.size(); ++i)
	{
		ls_file2 << ldis[i] << " ";
	}
	ps_file1.close();
	ls_file2.close();
	//tt_w 排序
	//建堆，tt_w保存了从第0个点到底n个点weight值，堆heap_vec保存每个patch的权值在tt_w中的下标，也是patch 的index
	vector<int> heap_vec;
	vector<int> ele_index;
	for (unsigned int i = 0; i < tt_w.size(); ++i)
	{
		ele_index.push_back(-1);
		if (mixedTopology.allNode[i].status >= 0)
		{
			heap_vec.push_back(i);
		}
		//heap_vec.push_back(i);
		//ele_index.push_back(0);
	}
	//建堆
	int h_size = heap_vec.size();
	for (int i = h_size / 2; i >= 1; i--)
	{
		heap_siftdown(i, &heap_vec, &tt_w, &ele_index);
	}
	//建立索引，知道某个权值的下标，可以直接找到在heap的位置，位置是从1开始，注意更新
	for (int i = 0; i < h_size; ++i)
	{
		ele_index[heap_vec[i]] = i;// ele_index
	}
	//////////////////
	//h_size是目前堆中的patch的数量，每次合并一个patch,h_size--，这是以数量为停止条件
	int dele_count_line = 0;
	int dele_count_patch = 0;
	while (tt_w[heap_vec.front()] < similarity_t && tt_w[heap_vec.front()] >= 0)
	{
		int pindex = heap_vec[0];
		if (tt_w[pindex] == 6553500)
		{
			cout << "left:" << h_size << endl;
			break;//若final_pnum设置的过小，使weight=65535，说明剩下的patch都没有neighbor,可以直接跳出循环，此时最小weight对应的patchindex， weight_pindex[pindex]=-1
		}
		int combine_pindex = -1;
		MixedArc* parc = mixedTopology.allNode[pindex].firstAdjNode;
		while (parc != NULL)
		{
			if (parc->weight == tt_w[pindex])
			{
				combine_pindex = parc->adj_index;
				break;
			}
			parc = parc->next_arc;
		}
		int ttt = weight_pindex[pindex];
		if (ttt != combine_pindex)
			cout << "error";
		combine_pindex = weight_pindex[pindex];//第pindex个patch最小权值对应的那个patch的index
											   //1、merge 处理
		if (mixedTopology.allNode[pindex].pl_type == 'p')
		{
			if (mixedTopology.allNode[combine_pindex].pl_type == 'p')//pp
			{
				dele_count_patch++;
				bool succ = combineTwoPatch(combine_pindex, pindex);//combine_pindex的周长，面积，center也会更新
				if (succ == false)
					cout << "combine patch error" << endl;
			}
			else//pl
			{
				cout << "?" << endl;
			}
		}
		else//注意line的index，需要减去pchsize
		{
			if (mixedTopology.allNode[combine_pindex].pl_type == 'p')//lp
			{
				cout << "?" << endl;
			}
			else//ll
			{
				dele_count_line++;
				vector<mypoint2f> nplist;
				combineTwoCurves(combine_pindex - pchsize, pindex - pchsize, &nplist);
			}
		}
		//2.先删除堆顶的pindex
		int tmp_back = heap_vec.back();
		int tmp_mini = heap_vec[0];
		tt_w[tmp_mini] = -1;
		ele_index[tmp_mini] = -1;
		ele_index[tmp_back] = 0;//!!!
		heap_vec[0] = tmp_back;
		heap_siftdown(1, &heap_vec, &tt_w, &ele_index);//下沉调整
		heap_vec.pop_back();
		h_size--;
		//3.更新combine_pindex的FD /combine_pindex-pchsize
		vector<mypoint2f> total_plist;
		vector<double> tmppchFD;
		if (mixedTopology.allNode[combine_pindex].pl_type == 'p')
		{
			jointMethod(&plrGraph.patches[combine_pindex].pointIndexOfPatch, &total_plist, 0); //total_plist首 = 尾
			total_plist.pop_back();
			getFDFeature(&total_plist, &plrGraph.patches[combine_pindex].pointIndexOfPatch, combine_pindex, simply_num, similarity_type, &tmppchFD);
		}
		else
		{
			vector<mypoint2f> attach_plist;
			jointMethod(&plrGraph.otherlines[combine_pindex - pchsize].pt_vec, &attach_plist, 0);  //attach + isolate
			total_plist.insert(total_plist.end(), attach_plist.begin(), attach_plist.end());
			for (int ppi = (int)attach_plist.size() - 2; ppi > 0; --ppi) //front!=back
			{
				total_plist.push_back(attach_plist[ppi]);
			}
			vector<int> line_pvec;
			line_pvec = plrGraph.otherlines[combine_pindex - pchsize].pt_vec;
			for (int j = plrGraph.otherlines[combine_pindex - pchsize].pt_vec.size() - 2; j >= 0; --j)
			{
				line_pvec.push_back(plrGraph.otherlines[combine_pindex - pchsize].pt_vec[j]);
			}
			getFDFeature(&total_plist, &line_pvec, -1, simply_num, similarity_type, &tmppchFD);//该函数内部包含了plist的化简，参数：第0种FD计算方法，或者第1种
		}

		////1、圆度
		//double circularity = 4 * 3.1415926*plrGraph.patches[combine_pindex].area / pow(plrGraph.patches[combine_pindex].perimeter, 2);
		////2、长度，计算出最大矩形bounding box
		//double elongation = 0;
		//vector<mypoint2f> plist_box;
		//mypoint2f rec_center(0, 0);
		//getBoundRectangle(&total_plist, rec_center, &plist_box);
		//double box_w = fabs(plist_box[0].x - plist_box[1].x);
		//double box_h = fabs(plist_box[1].y - plist_box[2].y);
		//if (box_w < box_h)
		//	elongation = box_w / box_h;
		//else
		//	elongation = box_h / box_w;
		////3、面积比，object_area/min-bounding rectangle_area,  convex or concave?
		//double are_rect = box_w*box_h;
		//double are_ratio = plrGraph.patches[combine_pindex].area / are_rect;
		//vector<double> fvec;//特征向量
		//fvec.push_back(circularity);//圆度
		//fvec.push_back(elongation);//长度
		//fvec.push_back(are_ratio);//面积
		//patch_tmpfeature[combine_pindex].second = fvec;

		mix_FDs[combine_pindex].clear(); mix_FDs[combine_pindex].swap(vector<double>());
		mix_FDs[combine_pindex] = tmppchFD;

		//3. pindex 的 status变为-1，且把其neighbor中的pindex换成combine_pindex，
		mixedTopology.allNode[pindex].status = -1;//!!!
		MixedArc* pindex_nei = mixedTopology.allNode[pindex].firstAdjNode;
		while (pindex_nei != NULL)
		{
			MixedArc* tmp_arc = getMixedArc(pindex_nei->adj_index, combine_pindex);
			if (tmp_arc != NULL)
			{
				removeMixedNode(pindex_nei->adj_index, pindex);
			}
			else
			{
				tmp_arc = getMixedArc(pindex_nei->adj_index, pindex);
				tmp_arc->adj_index = combine_pindex;
			}
			pindex_nei = pindex_nei->next_arc;
		}//注意结果是：combine_pindex中的pindex改成了combine_pindex，pindex中还连着各个neighbor
		 //修改combine_pindex的 topology ，即清空combine_pindex的neighbor，重新添加
		vector<int> nei;
		MixedArc* o_arc = mixedTopology.allNode[pindex].firstAdjNode;//pindex 的neighbor删除了pindex，但是pindex还保留着这些neighbor的连接
		while (o_arc != NULL)
		{
			if (o_arc->adj_index != combine_pindex)//!!
				nei.push_back(o_arc->adj_index);
			o_arc = o_arc->next_arc;
		}
		o_arc = mixedTopology.allNode[combine_pindex].firstAdjNode;
		while (o_arc != NULL)
		{
			if (o_arc->adj_index != combine_pindex)//!!
				nei.push_back(o_arc->adj_index);
			o_arc = o_arc->next_arc;
		}
		deleRepeatPoint(&nei);
		mixedTopology.allNode[combine_pindex].firstAdjNode = NULL;////清空combine_pindex的neighbor，重新添加
		for (unsigned int ii = 0; ii < nei.size(); ++ii)
		{
			MixedArc* neiarc = new MixedArc;
			neiarc->adj_index = nei[ii];
			neiarc->weight = 0;
			neiarc->next_arc = mixedTopology.allNode[combine_pindex].firstAdjNode;
			mixedTopology.allNode[combine_pindex].firstAdjNode = neiarc;
		}
		//原来的方法，patch 和line 分开修改
		//if (mixedTopology.allNode[combine_pindex].pl_type == 'p')
		//{
		//	mixedTopology.allNode[combine_pindex].firstAdjNode = NULL;////清空neighbor，重新添加
		//	vector<int> compch_nei = getNeighborPatches(combine_pindex, 0);
		//	for (unsigned int j = 0; j < compch_nei.size(); ++j)
		//	{
		//		MixedArc* neiarc = new MixedArc;
		//		neiarc->adj_index = compch_nei[j];
		//		neiarc->weight = 0;
		//		neiarc->next_arc = mixedTopology.allNode[combine_pindex].firstAdjNode;
		//		mixedTopology.allNode[combine_pindex].firstAdjNode = neiarc;
		//	}
		//}
		//else
		//{
		//	vector<int> nei;
		//	MixedArc* o_arc = mixedTopology.allNode[pindex].firstAdjNode;//pindex 的neighbor删除了pindex，但是pindex还保留着这些neighbor的连接
		//	while (o_arc != NULL)
		//	{
		//		if (o_arc->adj_index != combine_pindex)//!!
		//			nei.push_back(o_arc->adj_index);
		//		o_arc = o_arc->next_arc;
		//	}
		//	o_arc = mixedTopology.allNode[combine_pindex].firstAdjNode;
		//	while (o_arc != NULL)
		//	{
		//		if (o_arc->adj_index != combine_pindex)//!!
		//			nei.push_back(o_arc->adj_index);
		//		o_arc = o_arc->next_arc;
		//	}
		//	deleRepeatPoint(&nei);
		//	mixedTopology.allNode[combine_pindex].firstAdjNode = NULL;////清空combine_pindex的neighbor，重新添加
		//	for (unsigned int ii = 0; ii < nei.size(); ++ii)
		//	{				
		//		MixedArc* neiarc = new MixedArc;
		//		neiarc->adj_index = nei[ii];
		//		neiarc->weight = 0;
		//		neiarc->next_arc = mixedTopology.allNode[combine_pindex].firstAdjNode;
		//		mixedTopology.allNode[combine_pindex].firstAdjNode = neiarc;		
		//	}	
		//}

		//4.计算combine_patch 与所有neighbor的相似度weight，tt_w[],更新在堆中的位置
		if (mixedTopology.allNode[combine_pindex].pl_type == 'p')
		{
			parc = mixedTopology.allNode[combine_pindex].firstAdjNode;
			double mini_simi = 6553500;
			int simi_pindex = -1;
			while (parc != NULL)
			{
				if (mixedTopology.allNode[parc->adj_index].pl_type == 'l' || mixedTopology.allNode[parc->adj_index].pl_type == 'i')
				{
					parc->weight = 6553500;
					MixedArc* rever_ptcharc = getMixedArc(parc->adj_index, combine_pindex);
					rever_ptcharc->weight = 6553500;
					parc = parc->next_arc;
				}
				else
				{
					vector<int> combine_pvec;//合并的pointlist
					combineTwoPatch_getpvec(combine_pindex, parc->adj_index, &combine_pvec);//combine_plist  front()！= back()
					if (!combine_pvec.empty())
					{
						vector<mypoint2f> combine_plist;//合并的pointlist
						jointMethod(&combine_pvec, &combine_plist, 0);
						vector<double> newpchFD;//FD
						combine_plist.pop_back();  //front()！= back()	
						getFDFeature(&combine_plist, &combine_pvec, -1, simply_num, similarity_type, &newpchFD);
						double dis1 = getEuclidDistance(&newpchFD, &mix_FDs[combine_pindex]);
						double dis2 = getEuclidDistance(&newpchFD, &mix_FDs[parc->adj_index]);
						double mini_dis = dis1 < dis2 ? dis1 : dis2;
						parc->weight = mini_dis;
						MixedArc* rever_ptcharc = getMixedArc(parc->adj_index, combine_pindex);
						rever_ptcharc->weight = mini_dis;
						if (mini_dis < mini_simi)
						{
							mini_simi = mini_dis;
							simi_pindex = parc->adj_index;//
						}
						//vector<double> ft_vec;//3ge特征
						//getThreeFeature(&ft_vec, &combine_plist, 0, 0);
						//double dis1 = getEuclidDistance(&ft_vec, &mix_FDs[combine_pindex]);
						//parc->weight = dis1;
						//if (dis1 < mini_simi)
						//{
						//	mini_simi = dis1;
						//	simi_pindex = parc->adj_index;//
						//}
						//double dis2 = getEuclidDistance(&ft_vec, &mix_FDs[parc->adj_index]);
						//PatchArc* rever_ptcharc = getPatchArc(parc->adj_index, combine_pindex);
						//rever_ptcharc->weight = dis2;

						//重新选parc->adj_patchindex的最小weight
						double arc_minweight = 6553500;
						int arc_minipindex = -1;
						MixedArc* arc_neighbor = mixedTopology.allNode[parc->adj_index].firstAdjNode;
						while (arc_neighbor != NULL)
						{
							if (arc_neighbor->weight < arc_minweight && arc_neighbor->weight != 0)
							{
								arc_minweight = arc_neighbor->weight;
								arc_minipindex = arc_neighbor->adj_index;
							}
							arc_neighbor = arc_neighbor->next_arc;
						}
						//跟新
						heap_update(ele_index[parc->adj_index] + 1, arc_minweight, &heap_vec, &tt_w, &ele_index);
						weight_pindex[parc->adj_index] = arc_minipindex;
						parc = parc->next_arc;
					}
					else
					{
						int tmp_dele = parc->adj_index;
						parc = parc->next_arc;
						removeMixedNode(tmp_dele, combine_pindex);
						removeMixedNode(combine_pindex, tmp_dele);
						double arc_minweight = 6553500;
						int arc_minipindex = -1;
						MixedArc* arc_neighbor = mixedTopology.allNode[tmp_dele].firstAdjNode;
						while (arc_neighbor != NULL)
						{
							if (arc_neighbor->weight < arc_minweight && arc_neighbor->weight != 0)
							{
								arc_minweight = arc_neighbor->weight;
								arc_minipindex = arc_neighbor->adj_index;
							}
							arc_neighbor = arc_neighbor->next_arc;
						}
						//跟新
						heap_update(ele_index[tmp_dele] + 1, arc_minweight, &heap_vec, &tt_w, &ele_index);
						weight_pindex[tmp_dele] = arc_minipindex;
					}
				}
				//parc = parc->nextPatch;
			}
			heap_update(ele_index[combine_pindex] + 1, mini_simi, &heap_vec, &tt_w, &ele_index);//注意+1	tt_w[combine_pindex] = dis1;	
			weight_pindex[combine_pindex] = simi_pindex;
		}
		else
		{
			parc = mixedTopology.allNode[combine_pindex].firstAdjNode;
			double mini_simi = 6553500;
			int simi_pindex = -1;
			while (parc != NULL)
			{
				if (mixedTopology.allNode[parc->adj_index].pl_type == 'p')
				{
					parc->weight = 6553500;
					MixedArc* rever_ptcharc = getMixedArc(parc->adj_index, combine_pindex);
					rever_ptcharc->weight = 6553500;
					parc = parc->next_arc;
				}
				else
				{
					vector<mypoint2f> newplist;
					bool lsucflag = combineTwoCurves_getplist(combine_pindex - pchsize, parc->adj_index - pchsize, &newplist, patch_averWidth);
					if (lsucflag == true)
					{
						vector<mypoint2f> newpatch_plist;
						newpatch_plist.insert(newpatch_plist.end(), newplist.begin(), newplist.end());
						for (int ppi = (int)newplist.size() - 2; ppi > 0; --ppi) //newpatch_plist.front != back
						{
							newpatch_plist.push_back(newplist[ppi]);
						}
						vector<int> line_pvec;//line的patch的pvec,没有是空的
						vector<double> new_FD; //line 的FD
						getFDFeature(&newpatch_plist, &line_pvec, -1, simply_num, similarity_type, &new_FD);//参数：第0种FD计算方法，或者第1种
						double dis1 = getEuclidDistance(&mix_FDs[combine_pindex], &new_FD);
						double dis2 = getEuclidDistance(&mix_FDs[parc->adj_index], &new_FD);

						parc->weight = dis1 < dis2 ? dis1 : dis2;//取较小的那个值
						MixedArc* parc_reve = getMixedArc(parc->adj_index, combine_pindex);//patchTopology.allpatches[parc->adj_patchindex].firstAdjPatch;
						if (parc_reve != NULL)
							parc_reve->weight = parc->weight;
						else
							cout << "patch arc=NULL" << endl;
						if (parc->weight < mini_simi)
						{
							mini_simi = parc->weight;
							simi_pindex = parc->adj_index;
						}
						//重新选parc->adj_patchindex的最小weight
						double arc_minweight = 6553500;
						int arc_minipindex = -1;
						MixedArc* arc_neighbor = mixedTopology.allNode[parc->adj_index].firstAdjNode;
						while (arc_neighbor != NULL)
						{
							if (arc_neighbor->weight < arc_minweight && arc_neighbor->weight != 0)
							{
								arc_minweight = arc_neighbor->weight;
								arc_minipindex = arc_neighbor->adj_index;
							}
							arc_neighbor = arc_neighbor->next_arc;
						}
						//更新
						heap_update(ele_index[parc->adj_index] + 1, arc_minweight, &heap_vec, &tt_w, &ele_index);
						weight_pindex[parc->adj_index] = arc_minipindex;
						parc = parc->next_arc;
					}
					else   //两条line 不能merge
					{
						int tmp_dele = parc->adj_index;
						parc = parc->next_arc;
						removeMixedNode(tmp_dele, combine_pindex);
						removeMixedNode(combine_pindex, tmp_dele);
						double arc_minweight = 6553500;
						int arc_minipindex = -1;
						MixedArc* arc_neighbor = mixedTopology.allNode[tmp_dele].firstAdjNode;
						while (arc_neighbor != NULL)
						{
							if (arc_neighbor->weight < arc_minweight && arc_neighbor->weight != 0)
							{
								arc_minweight = arc_neighbor->weight;
								arc_minipindex = arc_neighbor->adj_index;
							}
							arc_neighbor = arc_neighbor->next_arc;
						}
						//跟新
						heap_update(ele_index[tmp_dele] + 1, arc_minweight, &heap_vec, &tt_w, &ele_index);
						weight_pindex[tmp_dele] = arc_minipindex;
					}
				}
				//while一次循环结束
			}
			heap_update(ele_index[combine_pindex] + 1, mini_simi, &heap_vec, &tt_w, &ele_index);//注意+1	tt_w[combine_pindex] = dis1;	
			weight_pindex[combine_pindex] = simi_pindex;
		}
	}
}


bool MyGraph::combinePL_getplist(int pindex, int lindex, vector<mypoint2f>* transform_pch)
{
	//pindex和lindex分别是plrGraph.patch/otherlines中的。 lindex=mixindex-pachsize
	//满足条件 才合并：接触点不是拐点 + line的方向和patch某条边的方向一致
	if (plrGraph.otherlines[lindex].type < 0)
	{
		return false;
	}
	//1、判断line和patch的位置
	int touchpoint = -1;
	int pch_size = plrGraph.patches[pindex].pointIndexOfPatch.size();
	for (unsigned int i = 0; i < plrGraph.patches[pindex].point_line.size(); ++i)
	{
		if (plrGraph.patches[pindex].point_line[i].second == lindex)//lindex既可以是attachlines的index ，也是otherline 的index
		{
			touchpoint = plrGraph.patches[pindex].point_line[i].first;
			break;
		}
	}
	bool flag = false;
	if (touchpoint < 0)
		return flag;
	vector<int> seg_1;  //若连接处不是拐点，则记录下该点两边的线
	vector<int> seg_2;
	double theta_threshold = 140;
	double theta = 0;  //touchpoint点处的拐角
	vector<int> pindex_patch = plrGraph.patches[pindex].pointIndexOfPatch;
	for (int i = 0; i < pch_size; ++i)
	{
		if (touchpoint == pindex_patch[i])
		{
			if (i > 0 && i < (pch_size - 1))//若touch 点在两端
			{
				//若touch 点不在两端,调整到两端
				if (i >(pch_size / 2))
				{
					//拆后面的，补到前面
					int stop_index = pindex_patch[i];
					pindex_patch.pop_back();
					vector<int>::iterator iter_pch = pindex_patch.begin();
					while (pindex_patch.back() != stop_index)
					{
						int pop_index = pindex_patch.back();
						pindex_patch.pop_back();
						iter_pch = pindex_patch.insert(iter_pch, pop_index);
					}
					iter_pch = pindex_patch.insert(iter_pch, stop_index);
				}
				else
				{
					//删除前面的，push到后面
					int stop_index = pindex_patch[i];
					vector<int>::iterator iter_pch = pindex_patch.begin();
					iter_pch = pindex_patch.erase(iter_pch);///????
					while (pindex_patch.front() != stop_index)
					{
						int erase_index = pindex_patch.front();
						iter_pch = pindex_patch.erase(iter_pch);
						pindex_patch.push_back(erase_index);
					}
					pindex_patch.push_back(stop_index);
				}
				i = 0;
			}
			int p1 = pindex_patch[i];
			mypoint2f p1_posi = sgraph.nodeList[p1].position;
			int p2 = pindex_patch[i + 1];
			int p3 = pindex_patch[pch_size - 2];
			vector<int> ttpvec;
			vector<mypoint2f> ttplist;
			ttpvec.push_back(p1); ttpvec.push_back(p2);
			jointMethod(&ttpvec, &ttplist, 0);
			mypoint2f p2_posi = ttplist[1];
			ttpvec[1] = p3;
			ttplist.swap(vector<mypoint2f>());
			jointMethod(&ttpvec, &ttplist, 0);
			mypoint2f p3_posi = ttplist[1];
			double costheta = getCosineAngle(p1_posi, p2_posi, p1_posi, p3_posi);
			theta = acos(costheta) * 180 / 3.1415926;
			int adj_pch_seg1 = -1;
			int adj_pch_seg2 = -1;
			EdgeOfPatch * anedge = getBorderofPatch(pindex, p1, p2);
			if (anedge->adjPatchIndex.size() > 1)
				adj_pch_seg1 = anedge->adjPatchIndex.back();
			anedge = getBorderofPatch(pindex, p1, p3);
			if (anedge->adjPatchIndex.size() > 1)
				adj_pch_seg2 = anedge->adjPatchIndex.back();
			if (theta > theta_threshold)//|| theta <45)//设置若角度大于140说明不是拐角，可以merge /////150
			{
				flag = true;
				seg_1.push_back(p1);
				double stop_theta = 180;
				int i_tmp = i + 1;
				mypoint2f last_posi = p1_posi;
				while (stop_theta > 160 && i_tmp<pch_size - 1)///////////165
				{
					int p_tmp = pindex_patch[i_tmp];
					seg_1.push_back(p_tmp);
					mypoint2f p_tmp_posi = sgraph.nodeList[p_tmp].position;
					int p_prior = pindex_patch[i_tmp + 1];
					mypoint2f prior_posi = sgraph.nodeList[p_prior].position;
					int next_adjpch = -1;
					EdgeOfPatch * nnedge = getBorderofPatch(pindex, p_tmp, p_prior);
					if (nnedge->adjPatchIndex.size() > 1)
						next_adjpch = nnedge->adjPatchIndex.back();
					if (next_adjpch != adj_pch_seg1)
						break;
					vector<int> sv1, sv2;
					sv1.push_back(pindex_patch[i_tmp]);
					sv1.push_back(pindex_patch[i_tmp - 1]);
					sv2.push_back(pindex_patch[i_tmp]);
					sv2.push_back(pindex_patch[i_tmp + 1]);
					vector<mypoint2f> sp1, sp2;
					jointMethod(&sv1, &sp1, 0);
					jointMethod(&sv2, &sp2, 0);
					double stp_cos = getCosineAngle(sp1[0], sp1[1], sp2[0], sp2[1]);//double stp_cos = getCosineAngle(p_tmp_posi, prior_posi, p_tmp_posi, last_posi);
					stop_theta = acos(stp_cos) * 180 / 3.1415926;
					last_posi = p_tmp_posi;
					i_tmp++;
				}
				seg_2.push_back(p1);//逆序
				stop_theta = 180;
				i_tmp = pch_size - 2;
				last_posi = p1_posi;
				while (stop_theta > 160 && i_tmp >0) ///////////165
				{
					int p_tmp = pindex_patch[i_tmp];
					seg_2.push_back(p_tmp);
					mypoint2f ptmp_posi = sgraph.nodeList[p_tmp].position;
					int p_tmp_prior = pindex_patch[i_tmp - 1];
					mypoint2f prior_posi = sgraph.nodeList[p_tmp_prior].position;
					int next_adjpch = -1;
					EdgeOfPatch * nnedge = getBorderofPatch(pindex, p_tmp, p_tmp_prior);
					if (nnedge->adjPatchIndex.size() > 1)
						next_adjpch = nnedge->adjPatchIndex.back();
					if (next_adjpch != adj_pch_seg2)
						break;
					vector<int> sv1, sv2;
					sv1.push_back(pindex_patch[i_tmp]);
					sv1.push_back(pindex_patch[i_tmp - 1]);
					sv2.push_back(pindex_patch[i_tmp]);
					sv2.push_back(pindex_patch[i_tmp + 1]);
					vector<mypoint2f> sp1, sp2;
					jointMethod(&sv1, &sp1, 0);
					jointMethod(&sv2, &sp2, 0);
					double stp_cos = getCosineAngle(sp1[0], sp1[1], sp2[0], sp2[1]);//double stp_cos = getCosineAngle(p_tmp_posi, prior_posi, p_tmp_posi, last_posi);
					stop_theta = acos(stp_cos) * 180 / 3.1415926;
					last_posi = ptmp_posi;
					i_tmp--;
				}
			}
			else
			{
				flag = false;//非拐角，非平坦区域。接近直角？
			}
			break;
		}
	}
	if (flag == false)
	{
		return flag;
	}
	//根据theta分情况处理
	if (theta > theta_threshold)//第一种处理
	{
		//2.判断seg_1或seg_2的方向，seg_1 seg_2都是从touchpoint出发的线段
		vector<mypoint2f> l_plist;
		vector<int> l_pvec = plrGraph.otherlines[lindex].pt_vec;//plrGraph.attachlines[lindex].pt_vec;
		if (l_pvec.front() == touchpoint)
			jointMethod(&l_pvec, &l_plist, 0);
		else
		{
			reverse(l_pvec.begin(), l_pvec.end());
			jointMethod(&l_pvec, &l_plist, 0);
		}
		double l_peri = getPerimeterofVec(&l_plist);
		vector<mypoint2f> seg_plist;
		vector<int> seg_pvec;
		seg_pvec = seg_1;
		jointMethod(&seg_1, &seg_plist, 0);
		double cos_theta1 = getCosineAngle(seg_plist.front(), seg_plist[1], l_plist[0], l_plist[1]);
		double theta1 = acos(cos_theta1) * 180 / 3.1415926;
		if (theta1 > (180 - theta_threshold))/////
		{
			seg_plist.swap(vector<mypoint2f>());
			seg_pvec.swap(vector<int>());
			jointMethod(&seg_2, &seg_plist, 0);
			seg_pvec = seg_2;
			double cos_theta2 = getCosineAngle(seg_plist.front(), seg_plist[1], l_plist[0], l_plist[1]);
			double theta2 = acos(cos_theta2) * 180 / 3.1415926;
			//double theta2 = getCosineAngle(seg_plist.front(), seg_plist[2], l_plist[0], l_plist[1]) * 180 / 3.1415926;
			if (theta2 > (180 - theta_threshold))////////
			{
				flag = false;//!!
				return flag;
			}
		}
		if (seg_plist.size() < 5 || l_plist.size()<5)//被影响的patch的边太小，不能合并patch+line////////////////从这往下不要？？？？
		{
			flag = false;//!!
			return flag;
		}
		double seg_peri = getPerimeterofVec(&seg_plist);

		//3.判断line的长度,curve fittng ,
		vector<mypoint2f> new_plist;
		if (l_peri > (seg_peri / 4 * 3))
		{
			//把l_plist变短为 seg_plist的3/4
			double accu_peri = 0;
			double shortperi = seg_peri / 4 * 3;
			int ppi = 1;
			while (accu_peri < shortperi && ppi < l_plist.size())
			{
				double tmplen = LineSegmentLength(l_plist[ppi - 1], l_plist[ppi]);
				accu_peri = accu_peri + tmplen;
				ppi++;
			}
			l_plist.erase(l_plist.begin() + ppi, l_plist.end());////?????		
		}
		vector<mypoint2f> tt_plist;
		tt_plist.insert(tt_plist.end(), l_plist.begin(), l_plist.end());
		tt_plist.insert(tt_plist.end(), seg_plist.begin(), seg_plist.end());
		cubicBezier bb;
		fitCubBezier(&tt_plist, &new_plist, &bb, 1);
		if (_isnan(bb.p1.x) || !_finite(bb.p1.x))
		{
			flag = false;
			return flag;
		}
		//4.判断两端点是否正确,不正确的话，返回false。否则 连接成patch
		if (!isPointRoughEqual(bb.p0, seg_plist[0]) && !isPointRoughEqual(bb.p0, seg_plist.back()))
		{
			flag = false;
		}
		else if (!isPointRoughEqual(bb.p3, seg_plist[0]) && !isPointRoughEqual(bb.p3, seg_plist.back()))
		{
			flag = false;
		}
		else
		{
			if (seg_pvec[1] == pindex_patch[pch_size - 2])
			{
				reverse(pindex_patch.begin(), pindex_patch.end());
			}
			transform_pch->insert(transform_pch->end(), new_plist.begin(), new_plist.end());
			int j = 0;
			while (j<seg_pvec.size())
			{
				if (pindex_patch[j] == seg_pvec[j])
					j++;
				else
					break;
			}
			vector<int> left_pvec;
			for (int i = j - 1; i < pch_size; ++i)
			{
				left_pvec.push_back(pindex_patch[i]);
			}
			jointMethod(&left_pvec, transform_pch, 0);
		}
	}
	//else if (theta<45)//第二种处理，x y方向放大一定倍数在平移，与line连接起来形成新patch
	//{
	//	//2.判断seg_1或seg_2的方向，seg_1 seg_2都是从touchpoint出发的线段
	//	vector<mypoint2f> l_plist;
	//	vector<int> l_pvec = plrGraph.otherlines[lindex].pt_vec;//plrGraph.attachlines[lindex].pt_vec;	
	//	if (l_pvec.front() == touchpoint)
	//		jointMethod(&l_pvec, &l_plist, 0);
	//	else
	//	{
	//		reverse(l_pvec.begin(), l_pvec.end());
	//		jointMethod(&l_pvec, &l_plist, 0);
	//	}
	//	//取一半的l_plist
	//	vector<mypoint2f> l_plist_half;
	//	int halfsize = l_plist.size() / 2;
	//	if (halfsize < 2)
	//		halfsize = 2;
	//	for (int i = 0; i < halfsize; ++i)
	//	{
	//		l_plist_half.push_back(l_plist[i]);
	//	}
	//	double l_peri = getPerimeterofVec(&l_plist_half);
	//	
	//	//计算夹角
	//	vector<mypoint2f> seg1_plist;
	//	jointMethod(&seg_1, &seg1_plist, 0);
	//	double cos_theta1 = getCosineAngle(seg1_plist.front(), seg1_plist[1], l_plist_half[0], l_plist_half[1]);
	//	double theta1 = acos(cos_theta1) * 180 / 3.1415926;
	//	vector<mypoint2f> seg2_plist;
	//	jointMethod(&seg_2, &seg2_plist, 0);
	//	double cos_theta2 = getCosineAngle(seg2_plist.front(), seg2_plist[1], l_plist_half[0], l_plist_half[1]);
	//	double theta2 = acos(cos_theta2) * 180 / 3.1415926;
	//	//3.选择夹角较小的那个seg
	//	vector<mypoint2f> seg_plist;
	//	vector<int> seg_pvec;
	//	if (theta1 < theta2)/////
	//	{
	//		if (theta1 < 150)  //if ((theta1+theta2) < 150)//	
	//		{
	//			flag = false;
	//			return flag;
	//		}
	//		else
	//		{
	//			seg_plist = seg1_plist;
	//			seg_pvec = seg_1;
	//		}
	//	}
	//	else
	//	{
	//		if (theta2 < 150) //if((theta1+theta2) < 150)//
	//		{
	//			flag = false;
	//			return flag;
	//		}
	//		else
	//		{
	//			seg_plist = seg2_plist;
	//			seg_pvec = seg_2;
	//		}
	//	}
	//	//4.若seg_pvec[0]->[1]有邻接的patch，不能merge ，令flag=false
	//	EdgeOfPatch* bder = getBorderofPatch(pindex, seg_pvec[0], seg_pvec[1]);
	//	ArcNode* aarc = getArcNode(seg_pvec[0], seg_pvec[1]);
	//	if (bder->adjPatchIndex.size()>1 || !aarc->patchesIndex.size()>1)
	//	{
	//		flag = false;
	//	}
	//	else
	//	{
	//		//seg_plist.back->seg_plist.front构成的矩形， 变成 seg_plist.back->l_plist_half.back构成的矩形
	//		double rect1_w = abs(seg_plist.back().x - seg_plist.front().x);
	//		double rect1_h = abs(seg_plist.back().y - seg_plist.front().y);
	//		double rect2_w = abs(seg_plist.back().x - l_plist_half.back().x);
	//		double rect2_h = abs(seg_plist.back().y - l_plist_half.back().y);
	//		double w_radio = 0;// 1 + (rect2_w - rect1_w) / rect1_w;
	//		double h_radio = 0;// 1 + (rect2_h - rect1_h) / rect1_h;
	//		if (rect1_w == 0)
	//			w_radio = rect2_w - rect1_w;
	//		else 
	//			w_radio = 1 + (rect2_w - rect1_w) / rect1_w;
	//		if(rect1_h==0)
	//			h_radio = rect2_h - rect1_h;
	//		else
	//			h_radio = 1 + (rect2_h - rect1_h) / rect1_h;
	//		vector<mypoint2f> new_plist;
	//		for (unsigned int i = 0; i < seg_plist.size(); ++i)
	//		{
	//			mypoint2f tmp_pt(seg_plist[i].x*w_radio, seg_plist[i].y*h_radio);
	//			new_plist.push_back(tmp_pt);
	//		}
	//		mypoint2f move_pt = new_plist.front() - l_plist_half.back();
	//		for (unsigned int i = 0; i < new_plist.size(); ++i)
	//		{
	//			new_plist[i] = new_plist[i] - move_pt;
	//		}
	//		//组成新的patch
	//		if (seg_pvec[1] == pindex_patch[pch_size - 2])
	//		{
	//			reverse(pindex_patch.begin(), pindex_patch.end());
	//		}
	//		for (int i = 0; i< l_plist_half.size(); ++i)
	//		{
	//			transform_pch->push_back(l_plist_half[i]);
	//		}
	//		transform_pch->pop_back();
	//		transform_pch->insert(transform_pch->end(), new_plist.begin(), new_plist.end());
	//		int j = 0;
	//		while (j<seg_pvec.size())
	//		{
	//			if (pindex_patch[j] == seg_pvec[j])
	//				j++;
	//			else
	//				break;
	//		}
	//		vector<int> left_pvec;
	//		for (int i = j - 1; i < pch_size; ++i)
	//		{
	//			left_pvec.push_back(pindex_patch[i]);
	//		}
	//		jointMethod(&left_pvec, transform_pch, 0);
	//	}
	//}
	return flag;
}
void MyGraph::combinePatchLine(int pindex, int lindex)
{
	//pindex和lindex分别是plrGraph.patch/otherlines中的。
	//满足条件 才合并：接触点不是拐点 + line的方向和patch某条边的方向一致
	//1、判断line和patch的位置
	vector<mypoint2f> origi_plist_line;
	jointMethod(&plrGraph.otherlines[lindex].pt_vec, &origi_plist_line, 0);
	pandl_combine_line.push_back(origi_plist_line);

	int touchpoint = -1;
	int pch_size = plrGraph.patches[pindex].pointIndexOfPatch.size();
	for (unsigned int i = 0; i < plrGraph.patches[pindex].point_line.size(); ++i)
	{
		if (plrGraph.patches[pindex].point_line[i].second == lindex)//lindex既可以是attachlines的index ，也是otherline 的index
		{
			touchpoint = plrGraph.patches[pindex].point_line[i].first;
			break;
		}
	}
	bool flag = false;
	if (touchpoint < 0)
		return;
	vector<int> seg_1;  //若连接处不是拐点，则记录下该点两边的线
	vector<int> seg_2;
	double theta_threshold = 140;
	double theta = 0;
	vector<int> pindex_patch = plrGraph.patches[pindex].pointIndexOfPatch;
	for (int i = 0; i < pch_size; ++i)
	{
		if (touchpoint == pindex_patch[i])
		{
			if (i > 0 && i < (pch_size - 1))//若touch 点不在两端,调整到两端
			{
				if (i >(pch_size / 2))
				{
					//拆后面的，补到前面
					int stop_index = pindex_patch[i];
					pindex_patch.pop_back();
					vector<int>::iterator iter_pch = pindex_patch.begin();
					while (pindex_patch.back() != stop_index)
					{
						int pop_index = pindex_patch.back();
						pindex_patch.pop_back();
						iter_pch = pindex_patch.insert(iter_pch, pop_index);
					}
					iter_pch = pindex_patch.insert(iter_pch, stop_index);
				}
				else
				{
					//删除前面的，push到后面
					int stop_index = pindex_patch[i];
					vector<int>::iterator iter_pch = pindex_patch.begin();
					iter_pch = pindex_patch.erase(iter_pch);///????
					while (pindex_patch.front() != stop_index)
					{
						int erase_index = pindex_patch.front();
						iter_pch = pindex_patch.erase(iter_pch);
						pindex_patch.push_back(erase_index);
					}
					pindex_patch.push_back(stop_index);
				}
				i = 0;
			}
			int p1 = pindex_patch[i];
			mypoint2f p1_posi = sgraph.nodeList[p1].position;
			int p2 = pindex_patch[i + 1];
			int p3 = pindex_patch[pch_size - 2];
			vector<int> ttpvec;
			vector<mypoint2f> ttplist;
			ttpvec.push_back(p1); ttpvec.push_back(p2);
			jointMethod(&ttpvec, &ttplist, 0);
			mypoint2f p2_posi = ttplist[1];
			ttpvec[1] = p3;
			ttplist.swap(vector<mypoint2f>());
			jointMethod(&ttpvec, &ttplist, 0);
			mypoint2f p3_posi = ttplist[1];
			double costheta = getCosineAngle(p1_posi, p2_posi, p1_posi, p3_posi);
			theta = acos(costheta) * 180 / 3.1415926;
			int adj_pch_seg1 = -1;
			int adj_pch_seg2 = -1;
			EdgeOfPatch * anedge = getBorderofPatch(pindex, p1, p2);
			if (anedge->adjPatchIndex.size() > 1)
				adj_pch_seg1 = anedge->adjPatchIndex.back();
			anedge = getBorderofPatch(pindex, p1, p3);
			if (anedge->adjPatchIndex.size() > 1)
				adj_pch_seg2 = anedge->adjPatchIndex.back();
			if (theta > theta_threshold)// || theta <45)//设置若角度大于150说明不是拐角，可以merge
			{
				flag = true;
				seg_1.push_back(p1);
				double stop_theta = 180;
				int i_tmp = i + 1;
				mypoint2f last_posi = p1_posi;
				while (stop_theta > 160 && i_tmp<pch_size - 1)
				{
					int p_tmp = pindex_patch[i_tmp];
					seg_1.push_back(p_tmp);
					mypoint2f ptmp_posi = sgraph.nodeList[p_tmp].position;
					int p_tmp_prior = pindex_patch[i_tmp + 1];
					mypoint2f prior_posi = sgraph.nodeList[p_tmp_prior].position;
					int next_adjpch = -1;
					EdgeOfPatch * nnedge = getBorderofPatch(pindex, p_tmp, p_tmp_prior);
					if (nnedge->adjPatchIndex.size() > 1)
						next_adjpch = nnedge->adjPatchIndex.back();
					if (next_adjpch != adj_pch_seg1)
						break;
					vector<int> sv1, sv2;
					sv1.push_back(pindex_patch[i_tmp]);
					sv1.push_back(pindex_patch[i_tmp - 1]);
					sv2.push_back(pindex_patch[i_tmp]);
					sv2.push_back(pindex_patch[i_tmp + 1]);
					vector<mypoint2f> sp1, sp2;
					jointMethod(&sv1, &sp1, 0);
					jointMethod(&sv2, &sp2, 0);
					double cos_st_theta = getCosineAngle(sp1[0], sp1[1], sp2[0], sp2[1]);//double cos_st_theta = getCosineAngle(ptmp_posi, prior_posi, ptmp_posi, last_posi);
					stop_theta = acos(cos_st_theta) * 180 / 3.1415926;
					last_posi = ptmp_posi;
					i_tmp++;
				}
				seg_2.push_back(p1);//逆序
				stop_theta = 180;
				i_tmp = pch_size - 2;
				last_posi = p1_posi;
				while (stop_theta > 160 && i_tmp >0)
				{
					int p_tmp = pindex_patch[i_tmp];
					seg_2.push_back(p_tmp);
					mypoint2f ptmp_posi = sgraph.nodeList[p_tmp].position;
					int p_tmp_prior = pindex_patch[i_tmp - 1];
					mypoint2f prior_posi = sgraph.nodeList[p_tmp_prior].position;
					int next_adjpch = -1;
					EdgeOfPatch * nnedge = getBorderofPatch(pindex, p_tmp, p_tmp_prior);
					if (nnedge->adjPatchIndex.size() > 1)
						next_adjpch = nnedge->adjPatchIndex.back();
					if (next_adjpch != adj_pch_seg2)
						break;
					vector<int> sv1, sv2;
					sv1.push_back(pindex_patch[i_tmp]);
					sv1.push_back(pindex_patch[i_tmp - 1]);
					sv2.push_back(pindex_patch[i_tmp]);
					sv2.push_back(pindex_patch[i_tmp + 1]);
					vector<mypoint2f> sp1, sp2;
					jointMethod(&sv1, &sp1, 0);
					jointMethod(&sv2, &sp2, 0);
					double cos_st_theta = getCosineAngle(sp1[0], sp1[1], sp2[0], sp2[1]);//double cos_st_theta = getCosineAngle(ptmp_posi, prior_posi, ptmp_posi, last_posi);
					stop_theta = acos(cos_st_theta) * 180 / 3.1415926;
					last_posi = ptmp_posi;
					i_tmp--;
				}
			}
			else
			{
				flag = false;//非拐角，非平坦区域。接近直角？
			}
			break;
		}
	}
	if (flag == false)
	{
		return;
	}
	//根据theta分情况处理
	vector<mypoint2f> transform_pch;  //变形后的newpatch  的 pointlist
	vector<mypoint2f> new_plist;     //变形的那条线的plist
	vector<mypoint2f> l_plist;
	vector<int> l_pvec = plrGraph.otherlines[lindex].pt_vec;//plrGraph.attachlines[lindex].pt_vec;
	if (l_pvec.front() == touchpoint)
		jointMethod(&l_pvec, &l_plist, 0);
	else
	{
		reverse(l_pvec.begin(), l_pvec.end());
		jointMethod(&l_pvec, &l_plist, 0);
	}
	double l_peri = getPerimeterofVec(&l_plist);

	vector<mypoint2f> seg_plist;
	vector<int> seg_pvec;
	double seg_peri = 0;
	if (theta > theta_threshold)//第一种处理
	{
		//2.判断seg_1或seg_2的方向		
		seg_pvec = seg_1;
		jointMethod(&seg_1, &seg_plist, 0);
		double cos_theta1 = getCosineAngle(seg_plist.front(), seg_plist[1], l_plist[0], l_plist[1]);
		double theta1 = acos(cos_theta1) * 180 / 3.1415926;
		if (theta1 > (180 - theta_threshold))
		{
			seg_plist.swap(vector<mypoint2f>());
			seg_pvec.swap(vector<int>());
			jointMethod(&seg_2, &seg_plist, 0);
			seg_pvec = seg_2;
			double cos_theta2 = getCosineAngle(seg_plist.front(), seg_plist[1], l_plist[0], l_plist[1]);
			double theta2 = acos(cos_theta2) * 180 / 3.1415926;
			if (theta2 > (180 - theta_threshold))
			{
				flag = false;
				return;
			}
		}
		seg_peri = getPerimeterofVec(&seg_plist);
		//3.判断line的长度,curve fittng ,
		//vector<mypoint2f> new_plist;
		if (l_peri > (seg_peri / 4 * 3))
		{
			//把l_plist变短为 seg_plist的3/4
			double accu_peri = 0;
			double shortperi = seg_peri / 4 * 3;
			int ppi = 1;
			while (accu_peri < shortperi && ppi < l_plist.size())
			{
				double tmplen = LineSegmentLength(l_plist[ppi - 1], l_plist[ppi]);
				accu_peri = accu_peri + tmplen;
				ppi++;
			}
			l_plist.erase(l_plist.begin() + ppi, l_plist.end());////!!!!!		
		}
		vector<mypoint2f> tt_plist;
		tt_plist.insert(tt_plist.end(), l_plist.begin(), l_plist.end());
		tt_plist.insert(tt_plist.end(), seg_plist.begin(), seg_plist.end());
		cubicBezier bb;
		fitCubBezier(&tt_plist, &new_plist, &bb, 1);
		//4.判断两端点是否正确,不正确的话，返回false。否则 修改各个变量 sgraph.nodeList ， plrGraph.patches ,plrGraph.attachlines/otherline
		if (!isPointRoughEqual(bb.p0, seg_plist[0]) && !isPointRoughEqual(bb.p0, seg_plist.back()))
		{
			flag = false;
		}
		else if (!isPointRoughEqual(bb.p3, seg_plist[0]) && !isPointRoughEqual(bb.p3, seg_plist.back()))
		{
			flag = false;
		}
		else
		{
			//删除sgraph中的 arc
			vector<int> rm_line = plrGraph.otherlines[lindex].pt_vec;
			for (size_t ti = 1; ti < rm_line.size(); ++ti)
			{
				removeArcNode(rm_line[ti - 1], rm_line[ti], 0);
				removeArcNode(rm_line[ti], rm_line[ti - 1], 0);
			}
			//重新连成patch，transform_pch.  计算面积 ,周长 ,中心
			if (seg_pvec[1] == pindex_patch[pch_size - 2])
			{
				reverse(pindex_patch.begin(), pindex_patch.end());
			}
			transform_pch.insert(transform_pch.end(), new_plist.begin(), new_plist.end());
			int pj = 0;
			while (pj<seg_pvec.size())
			{
				if (pindex_patch[pj] == seg_pvec[pj])
					pj++;
				else
					break;
			}
			vector<int> left_pvec;
			for (int i = pj - 1; i < pch_size; ++i)
			{
				left_pvec.push_back(pindex_patch[i]);
			}
			jointMethod(&left_pvec, &transform_pch, 0);//done
													   //修改sgraph.nodeList[。].firstEdge->curveType;   sgraph.nodeList[。].firstEdge->curveIndex;  polyline_vec[.].pointlist;
			if (seg_pvec.size() > 2)
			{
				vector<double> peri_ratio;
				vector<mypoint2f> sub_plist;
				for (unsigned int i = 1; i < seg_pvec.size(); ++i)
				{
					int sp1 = seg_pvec[i - 1];
					int sp2 = seg_pvec[i];
					mypoint2f posi1 = sgraph.nodeList[sp1].position;
					mypoint2f posi2 = sgraph.nodeList[sp2].position;
					if (!sub_plist.empty())
						sub_plist.pop_back();
					getPListOfaCurve(posi1, posi2, sp1, sp2, &sub_plist, 0);
					double sub_peri = getPerimeterofVec(&sub_plist);
					peri_ratio.push_back(sub_peri / seg_peri);
				}
				int last_subi = 0;
				for (unsigned int i = 0; i < peri_ratio.size(); ++i)//seg_pvec有3个，peri_ratio有2两个
				{
					int sub_i = new_plist.size()*peri_ratio[i] - 1;  //sub_i把new_plist按比例大体分成几段
					if (sub_i <= 0)
						sub_i = 1;
					vector<mypoint2f> modi_plist;
					for (int j = last_subi; j <= sub_i; ++j)
					{
						modi_plist.push_back(new_plist[j]);
					}
					ArcNode* anarc = getArcNode(seg_pvec[i], seg_pvec[i + 1]);
					exchangeline(&modi_plist, anarc);//polyline_vec中替换新的line：modi_plist
					last_subi = sub_i;
					////原来，直接修改点的位置
					//sgraph.nodeList[seg_pvec[i + 1]].position = new_plist[sub_i];/////!!!!!!!
					////现在，patch 中point_line中有seg_pvec[i+1]点所连接的线，整体平移
					mypoint2f move_dic = new_plist[sub_i] - sgraph.nodeList[seg_pvec[i + 1]].position;
					sgraph.nodeList[seg_pvec[i + 1]].position = new_plist[sub_i];
					if (!plrGraph.patches[pindex].point_line.empty())
					{
						for (unsigned int lli = 0; lli < plrGraph.patches[pindex].point_line.size(); ++lli)
						{
							if (plrGraph.patches[pindex].point_line[lli].first == seg_pvec[i + 1])
							{
								int move_lindex = plrGraph.patches[pindex].point_line[lli].second;
								vector<int> move_lvec = plrGraph.otherlines[move_lindex].pt_vec;
								for (unsigned int k = 1; k < move_lvec.size(); ++k)
								{
									ArcNode* llarc = getArcNode(move_lvec[k - 1], move_lvec[k]);
									vector<int> tmplvec;
									tmplvec.push_back(move_lvec[k - 1]);
									tmplvec.push_back(move_lvec[k]);
									vector<mypoint2f> move_plist;
									jointMethod(&tmplvec, &move_plist, 0);
									for (unsigned int g = 0; g < move_plist.size(); ++g)
										move_plist[g] = move_plist[g] + move_dic;
									exchangeline(&move_plist, llarc);
								}
							}
						}
					}

				}
			}
			else
			{
				ArcNode* anarc = getArcNode(seg_pvec[0], seg_pvec[1]);
				exchangeline(&new_plist, anarc);//polyline_vec中替换新的line：modi_plist，代替下面的代码
												//int modi_lindex = anarc->curveIndex;
												//char modi_ltype = anarc->curveType;
												//if (modi_ltype == 'L')//?????
												//{
												//	polyline_vec[modi_lindex].pointlist.swap(vector<mypoint2f>());
												//	polyline_vec[modi_lindex].pointlist = new_plist;/////!!!!!!!!
												//}
												//else if (modi_ltype == 'C' || modi_ltype == 'Q')
												//{
												//	anarc->curveType = 'L';
												//	PolyLine pll;
												//	pll.origi_index = polyline_vec.size();
												//	pll.origi_type = 'L';
												//	pll.pointlist = new_plist;
												//	polyline_vec.push_back(pll);
												//	anarc->curveIndex = polyline_vec.size() - 1;
												//}
			}
		}
	}
	//else if (theta < 45)
	//{
	//	//取一半的l_plist
	//	vector<mypoint2f> l_halfplist;
	//	int halfsize = l_plist.size() / 2;
	//	if (halfsize < 2)
	//		halfsize = 2;
	//	for (int i = 0; i < halfsize; ++i)
	//	{
	//		l_halfplist.push_back(l_plist[i]);
	//	}
	//	if (l_pvec.size() > 2)
	//	{
	//		vector<double> peri_ratio;
	//		vector<mypoint2f> sub_plist;
	//		for (unsigned int i = 1; i < l_pvec.size(); ++i)
	//		{
	//			int sp1 = l_pvec[i - 1];
	//			int sp2 = l_pvec[i];
	//			mypoint2f posi1 = sgraph.nodeList[sp1].position;
	//			mypoint2f posi2 = sgraph.nodeList[sp2].position;
	//			if (!sub_plist.empty())
	//				sub_plist.pop_back();
	//			getPListOfaCurve(posi1, posi2, sp1, sp2, &sub_plist, 0);
	//			double sub_peri = getPerimeterofVec(&sub_plist);
	//			peri_ratio.push_back(sub_peri / l_peri);
	//		}
	//		int last_subi = 0;
	//		for (unsigned int i = 0; i < peri_ratio.size(); ++i)//l_pvec有3个，peri_ratio有2两个
	//		{
	//			int sub_i = l_halfplist.size()*peri_ratio[i] - 1;  //sub_i把l_halfplist按比例大体分成几段
	//			if (sub_i <= 0)
	//				sub_i = 1;
	//			sgraph.nodeList[l_pvec[i + 1]].position = l_halfplist[sub_i];/////!!!!!!!
	//			vector<mypoint2f> modi_plist;
	//			for (int j = last_subi; j <= sub_i; ++j)
	//			{
	//				modi_plist.push_back(l_halfplist[j]);
	//			}
	//			ArcNode* anarc = getArcNode(l_pvec[i], l_pvec[i + 1]);
	//			exchangeline(&modi_plist, anarc);//polyline_vec中替换新的line：modi_plist，代替下面的代码
	//			last_subi = sub_i;
	//		}
	//	}
	//	else
	//	{
	//		sgraph.nodeList[l_pvec.back()].position = l_halfplist.back();
	//		ArcNode* anarc = getArcNode(l_pvec[0], l_pvec[1]);
	//		exchangeline(&l_halfplist, anarc);
	//	}
	//	//计算夹角
	//	//2.判断seg_1或seg_2的方向，seg_1 seg_2都是从touchpoint出发的线段
	//	vector<mypoint2f> seg1_plist;
	//	jointMethod(&seg_1, &seg1_plist, 0);
	//	double cos_theta1 = getCosineAngle(seg1_plist.front(), seg1_plist[1], l_halfplist[0], l_halfplist[1]);
	//	double theta1 = acos(cos_theta1) * 180 / 3.1415926;
	//	vector<mypoint2f> seg2_plist;
	//	jointMethod(&seg_2, &seg2_plist, 0);
	//	double cos_theta2 = getCosineAngle(seg2_plist.front(), seg2_plist[1], l_halfplist[0], l_halfplist[1]);
	//	double theta2 = acos(cos_theta2) * 180 / 3.1415926;
	//	//选择夹角较小的那个seg
	//	if (theta1 < theta2)/////
	//	{
	//		seg_plist = seg1_plist;
	//		seg_pvec = seg_1;
	//	}
	//	else
	//	{
	//		seg_plist = seg2_plist;
	//		seg_pvec = seg_2;
	//	}
	//	seg_peri = getPerimeterofVec(&seg_plist);
	//	//seg_plist.back->seg_plist.front构成的矩形， 变成 seg_plist.back->l_halfplist.back构成的矩形
	//	double rect1_w = abs(seg_plist.back().x - seg_plist.front().x);
	//	double rect1_h = abs(seg_plist.back().y - seg_plist.front().y);
	//	double rect2_w = abs(seg_plist.back().x - l_halfplist.back().x);
	//	double rect2_h = abs(seg_plist.back().y - l_halfplist.back().y);
	//	double w_radio = 1 + (rect2_w - rect1_w) / rect1_w;
	//	double h_radio = 1 + (rect2_h - rect1_h) / rect1_h;
	//	//vector<mypoint2f> new_plist;
	//	for (unsigned int i = 0; i < seg_plist.size(); ++i)
	//	{
	//		mypoint2f tmp_pt(seg_plist[i].x*w_radio, seg_plist[i].y*h_radio);
	//		new_plist.push_back(tmp_pt);
	//	}
	//	mypoint2f move_pt = new_plist.front() - l_halfplist.back();
	//	for (unsigned int i = 0; i < new_plist.size(); ++i)
	//	{
	//		new_plist[i] = new_plist[i] - move_pt;
	//	}
	//	//组成新的patch,(transform_pch)
	//	if (seg_pvec[1] == pindex_patch[pch_size - 2])
	//	{
	//		reverse(pindex_patch.begin(), pindex_patch.end());
	//	}
	//	for (int i = 0; i< l_halfplist.size(); ++i)
	//	{
	//		transform_pch.push_back(l_halfplist[i]);
	//	}
	//	transform_pch.pop_back();
	//	transform_pch.insert(transform_pch.end(), new_plist.begin(), new_plist.end());
	//	int k = 0;
	//	while (k<seg_pvec.size())
	//	{
	//		if (pindex_patch[k] == seg_pvec[k])
	//			k++;
	//		else
	//			break;
	//	}
	//	vector<int> left_pvec;
	//	for (int i = k - 1; i < pch_size; ++i)
	//	{
	//		left_pvec.push_back(pindex_patch[i]);
	//	}
	//	jointMethod(&left_pvec, &transform_pch, 0);
	//	//向newpatch 中添加point 和edge
	//	for (int i = 0; i < pch_size; ++i)
	//	{
	//		if (touchpoint == plrGraph.patches[pindex].pointIndexOfPatch[i])
	//		{
	//			if (i == 0)
	//			{
	//				if (plrGraph.patches[pindex].pointIndexOfPatch[1] == seg_pvec[1])
	//				{
	//					vector<int>::iterator iter_tmp = plrGraph.patches[pindex].pointIndexOfPatch.begin() + 1;
	//					for (unsigned int j = 1; j < l_pvec.size(); ++j)
	//					{
	//						iter_tmp = plrGraph.patches[pindex].pointIndexOfPatch.insert(iter_tmp, l_pvec[j]);
	//						iter_tmp++;
	//					}
	//				}
	//				else
	//				{
	//					vector<int>::iterator iter_tmp = plrGraph.patches[pindex].pointIndexOfPatch.begin() + pch_size - 1;//back
	//					for (unsigned int j = 1; j < l_pvec.size(); ++j)
	//					{
	//						iter_tmp = plrGraph.patches[pindex].pointIndexOfPatch.insert(iter_tmp, l_pvec[j]);
	//					}
	//				}
	//			}
	//			else
	//			{
	//				if (plrGraph.patches[pindex].pointIndexOfPatch[i - 1] == seg_pvec[1])
	//				{
	//					vector<int>::iterator iter_tmp = plrGraph.patches[pindex].pointIndexOfPatch.begin() + i;
	//					for (unsigned int j = 1; j < l_pvec.size(); ++j)
	//					{
	//						iter_tmp = plrGraph.patches[pindex].pointIndexOfPatch.insert(iter_tmp, l_pvec[j]);
	//					}
	//				}
	//				else
	//				{
	//					vector<int>::iterator iter_tmp = plrGraph.patches[pindex].pointIndexOfPatch.begin() + i + 1;
	//					for (unsigned int j = 1; j < l_pvec.size(); ++j)
	//					{
	//						iter_tmp = plrGraph.patches[pindex].pointIndexOfPatch.insert(iter_tmp, l_pvec[j]);
	//						iter_tmp++;
	//					}
	//				}
	//			}
	//			break;
	//		}
	//	}
	//	vector<int> adjpatch_vec = getBorderofPatch(pindex, seg_pvec[0], seg_pvec[1])->adjPatchIndex;
	//	vector<EdgeOfPatch> add_edge;
	//	for (unsigned int i = 1; i < plrGraph.patches[pindex].pointIndexOfPatch.size(); ++i)
	//	{
	//		int ss = plrGraph.patches[pindex].pointIndexOfPatch[i - 1];
	//		int ee = plrGraph.patches[pindex].pointIndexOfPatch[i];
	//		EdgeOfPatch* border=getBorderofPatch(pindex, ss, ee);
	//		if (border == NULL)
	//		{
	//			EdgeOfPatch newborder;
	//			newborder.startp_index = ss;
	//			newborder.endp_index = ee;
	//			newborder.edgeflag = 0;
	//			newborder.adjPatchIndex = adjpatch_vec;
	//			add_edge.push_back(newborder);
	//		}
	//	}
	//	vector<EdgeOfPatch>::iterator iter_border = plrGraph.patches[pindex].edges_vec.begin();
	//	while (iter_border != plrGraph.patches[pindex].edges_vec.end())
	//	{
	//		if (iter_border->startp_index == seg_pvec[0] && iter_border->endp_index == seg_pvec[1])//先删除1个edge 再插入
	//		{
	//			iter_border = plrGraph.patches[pindex].edges_vec.erase(iter_border);
	//			plrGraph.patches[pindex].edges_vec.insert(iter_border, add_edge.begin(), add_edge.end());
	//			break;
	//		}
	//		else if (iter_border->startp_index == seg_pvec[1] && iter_border->endp_index == seg_pvec[0])//先删除 在插入
	//		{
	//			iter_border = plrGraph.patches[pindex].edges_vec.erase(iter_border);
	//			plrGraph.patches[pindex].edges_vec.insert(iter_border, add_edge.begin(), add_edge.end());
	//			break;
	//		}
	//		iter_border++;
	//	}
	//	//原来的seg_vec[0,1,2] ,变成了[l_pvec.back,1,2]!! 
	//	//修改sgraph.nodeList[。].firstEdge->curveType;   sgraph.nodeList[。].firstEdge->curveIndex;  polyline_vec[.].pointlist;
	//	if (seg_pvec.size() > 2)
	//	{
	//		vector<double> peri_ratio;
	//		vector<mypoint2f> sub_plist;
	//		for (unsigned int i = 1; i < seg_pvec.size(); ++i)
	//		{
	//			int sp1 = seg_pvec[i - 1];
	//			int sp2 = seg_pvec[i];
	//			mypoint2f posi1 = sgraph.nodeList[sp1].position;
	//			mypoint2f posi2 = sgraph.nodeList[sp2].position;
	//			if (!sub_plist.empty())
	//				sub_plist.pop_back();
	//			getPListOfaCurve(posi1, posi2, sp1, sp2, &sub_plist, 0);
	//			double sub_peri = getPerimeterofVec(&sub_plist);
	//			peri_ratio.push_back(sub_peri / seg_peri);
	//		}
	//		int last_subi = 0;
	//		ArcNode store_arc = *getArcNode(seg_pvec[0], seg_pvec[1]);//!!!!!!!!!!!!!!
	//		removeArcNode(seg_pvec[0], seg_pvec[1], 0);
	//		removeArcNode(seg_pvec[1], seg_pvec[0], 0);
	//		//原来的seg_vec[0,1] ,变成了[l_pvec.back,1]!! ,
	//		ArcNode* anarc = new ArcNode;
	//		anarc->adjVertex = seg_pvec[1];
	//		anarc->curveType = 'L';
	//		anarc->curveIndex = polyline_vec.size();
	//		anarc->lineFlag = 0;
	//		anarc->patchesIndex = store_arc.patchesIndex;
	//		anarc->next = sgraph.nodeList[l_pvec.back()].firstEdge;
	//		sgraph.nodeList[l_pvec.back()].firstEdge = anarc;
	//		ArcNode* arc_rever = new ArcNode;
	//		arc_rever->adjVertex = l_pvec.back();
	//		arc_rever->curveType = 'L';
	//		arc_rever->curveIndex = polyline_vec.size();
	//		arc_rever->lineFlag = 0;
	//		arc_rever->patchesIndex = store_arc.patchesIndex;
	//		arc_rever->next = sgraph.nodeList[seg_pvec[1]].firstEdge;
	//		sgraph.nodeList[seg_pvec[1]].firstEdge = arc_rever;
	//		PolyLine apl;
	//		apl.origi_index = polyline_vec.size();
	//		apl.origi_type = 'L';
	//		apl.pointlist = new_plist;
	//		polyline_vec.push_back(apl);
	//		seg_pvec[0] = l_pvec.back();//原来的seg_vec[0,1] ,变成了[l_pvec.back,1]!!!!!!,
	//		for (unsigned int i = 0; i < peri_ratio.size(); ++i)//seg_pvec有3个，peri_ratio有2两个
	//		{
	//			int sub_i = new_plist.size()*peri_ratio[i] - 1;  //sub_i把new_plist按比例大体分成几段
	//			if (sub_i <= 0)
	//			{
	//				dataIncreasing(&new_plist, new_plist.size() * 2);
	//				sub_i = new_plist.size()*peri_ratio[i] - 1;
	//			}
	//			sgraph.nodeList[seg_pvec[i + 1]].position = new_plist[sub_i];/////!!!!!!!
	//			vector<mypoint2f> modi_plist;
	//			for (int j = last_subi; j <= sub_i; ++j)
	//			{
	//				modi_plist.push_back(new_plist[j]);
	//			}
	//			ArcNode* anarc = getArcNode(seg_pvec[i], seg_pvec[i + 1]);
	//			exchangeline(&modi_plist, anarc);
	//			
	//			//int modi_lindex = anarc->curveIndex;
	//			//char modi_ltype = anarc->curveType;
	//			//if (modi_ltype == 'L')//?????
	//			//{
	//			//	polyline_vec[modi_lindex].pointlist.swap(vector<mypoint2f>());
	//			//	polyline_vec[modi_lindex].pointlist = modi_plist;/////!!!!!!!!
	//			//}
	//			//else if (modi_ltype == 'C' || modi_ltype == 'Q')
	//			//{
	//			//	anarc->curveType = 'L';
	//			//	PolyLine pll;
	//			//	pll.origi_index = polyline_vec.size();
	//			//	pll.origi_type = 'L';
	//			//	pll.pointlist = modi_plist;
	//			//	polyline_vec.push_back(pll);
	//			//	anarc->curveIndex = polyline_vec.size() - 1;
	//			//}
	//			last_subi = sub_i;
	//		}
	//	}
	//	else
	//	{
	//		ArcNode store_arc = *getArcNode(seg_pvec[0], seg_pvec[1]);//!!!!!!!!!!!!!!
	//		removeArcNode(seg_pvec[0], seg_pvec[1], 0);
	//		removeArcNode(seg_pvec[1], seg_pvec[0], 0);
	//		//原来的seg_vec[0,1] ,变成了[l_pvec.back,1]!! ,
	//		ArcNode* anarc = new ArcNode;
	//		anarc->adjVertex = seg_pvec[1];
	//		anarc->curveType = 'L';
	//		anarc->curveIndex = polyline_vec.size();
	//		anarc->lineFlag = 0;
	//		anarc->patchesIndex = store_arc.patchesIndex;
	//		anarc->next = sgraph.nodeList[l_pvec.back()].firstEdge;
	//		sgraph.nodeList[l_pvec.back()].firstEdge = anarc;
	//		ArcNode* arc_rever = new ArcNode;
	//		arc_rever->adjVertex = l_pvec.back();
	//		arc_rever->curveType = 'L';
	//		arc_rever->curveIndex = polyline_vec.size();
	//		arc_rever->lineFlag = 0;
	//		arc_rever->patchesIndex = store_arc.patchesIndex;
	//		arc_rever->next = sgraph.nodeList[seg_pvec[1]].firstEdge;
	//		sgraph.nodeList[seg_pvec[1]].firstEdge = arc_rever;
	//		PolyLine apl;
	//		apl.origi_index = polyline_vec.size();
	//		apl.origi_type = 'L';
	//		apl.pointlist = new_plist;
	//		polyline_vec.push_back(apl);
	//	}
	//}
	if (flag == true)
	{
		//patch 中删除point_line
		vector<pair<int, int>>::iterator iter_pl = plrGraph.patches[pindex].point_line.begin();
		while (iter_pl != plrGraph.patches[pindex].point_line.end())
		{
			if (iter_pl->first == touchpoint && iter_pl->second == lindex)
				iter_pl = plrGraph.patches[pindex].point_line.erase(iter_pl);
			else
				iter_pl++;
		}
		//对attachline/otherline的处理
		plrGraph.otherlines[lindex].type = -1;  //	plrGraph.attachlines[lindex].type = -1;

												//计算面积 周长 中心点
		double new_area = GetPolygonArea(&transform_pch);
		double new_peri = getPerimeterofVec(&transform_pch);
		mypoint2f sum_pt(0, 0);
		for (unsigned int i = 0; i < transform_pch.size(); ++i)
		{
			sum_pt = sum_pt + transform_pch[i];
		}
		sum_pt.x = sum_pt.x / transform_pch.size();
		sum_pt.y = sum_pt.y / transform_pch.size();
		plrGraph.patches[pindex].area = new_area;
		plrGraph.patches[pindex].perimeter = new_peri;
		plrGraph.patches[pindex].centre_position = sum_pt;
		//测试
		pandl_combine.push_back(transform_pch);//////////////////////////////
											   /*
											   sgraph.nodeList[0].firstEdge->curveType; sgraph.nodeList[0].firstEdge->curveIndex; polyline_vec[0].pointlist;
											   plrGraph.patches[0].area;
											   plrGraph.patches[0].centre_position;
											   plrGraph.patches[0].perimeter;
											   plrGraph.patches[0].point_line.erase();
											   plrGraph.attachlines[0].type = -1;
											   plrGraph.otherlines[0].type = -1;*/
	}
}
bool MyGraph::mergePL_getplist(int pindex, int lindex, vector<mypoint2f>* transform_pch)
{
	//pindex和lindex分别是plrGraph.patch/otherlines中的index
	//满足条件 才合并：接触点不是拐点 
	if (plrGraph.otherlines[lindex].type < 0)
	{
		return false;
	}
	//1、判断line和patch的位置
	int touchpoint = -1;
	int pch_size = plrGraph.patches[pindex].pointIndexOfPatch.size();
	for (unsigned int i = 0; i < plrGraph.patches[pindex].point_line.size(); ++i)
	{
		if (plrGraph.patches[pindex].point_line[i].second == lindex)//lindex既可以是attachlines的index ，也是otherline 的index
		{
			touchpoint = plrGraph.patches[pindex].point_line[i].first;
			break;
		}
	}
	bool flag = false;
	if (touchpoint < 0)
		return flag;
	vector<int> seg_1;  //若连接处不是拐点，则记录下该点 的左右两边
	vector<int> seg_2;
	double theta_thread = 140;//广角>140 ， line与patch 的edge的角度<40度
	double theta = 0;  //touchpoint点处的拐角
	vector<int> pindex_patch = plrGraph.patches[pindex].pointIndexOfPatch;
	for (int i = 0; i < pch_size; ++i)
	{
		if (touchpoint == pindex_patch[i])
		{
			if (i > 0 && i < (pch_size - 1)) //若touch 点不在两端,调整到两端
			{
				if (i >(pch_size / 2))
				{
					//拆后面的，补到前面
					int stop_index = pindex_patch[i];
					pindex_patch.pop_back();
					vector<int>::iterator iter_pch = pindex_patch.begin();
					while (pindex_patch.back() != stop_index)
					{
						int pop_index = pindex_patch.back();
						pindex_patch.pop_back();
						iter_pch = pindex_patch.insert(iter_pch, pop_index);
					}
					iter_pch = pindex_patch.insert(iter_pch, stop_index);
				}
				else
				{
					//删除前面的，push到后面
					int stop_index = pindex_patch[i];
					vector<int>::iterator iter_pch = pindex_patch.begin();
					iter_pch = pindex_patch.erase(iter_pch);///????
					while (pindex_patch.front() != stop_index)
					{
						int erase_index = pindex_patch.front();
						iter_pch = pindex_patch.erase(iter_pch);
						pindex_patch.push_back(erase_index);
					}
					pindex_patch.push_back(stop_index);
				}
				i = 0;
			}
			int p1 = pindex_patch[i];
			mypoint2f p1_posi = sgraph.nodeList[p1].position;
			int p2 = pindex_patch[i + 1];
			int p3 = pindex_patch[pch_size - 2];
			vector<int> ttpvec;
			vector<mypoint2f> ttplist;
			ttpvec.push_back(p1); ttpvec.push_back(p2);
			jointMethod(&ttpvec, &ttplist, 0);
			mypoint2f p2_posi = ttplist[1];
			ttpvec[1] = p3;
			ttplist.swap(vector<mypoint2f>());
			jointMethod(&ttpvec, &ttplist, 0);
			mypoint2f p3_posi = ttplist[1];
			double costheta = getCosineAngle(p1_posi, p2_posi, p1_posi, p3_posi);
			theta = acos(costheta) * 180 / 3.1415926;
			int adj_pch_seg1 = -1;
			int adj_pch_seg2 = -1;
			EdgeOfPatch * anedge = getBorderofPatch(pindex, p1, p2);
			if (anedge->adjPatchIndex.size() > 1)
				adj_pch_seg1 = anedge->adjPatchIndex.back();
			anedge = getBorderofPatch(pindex, p1, p3);
			if (anedge->adjPatchIndex.size() > 1)
				adj_pch_seg2 = anedge->adjPatchIndex.back();
			if (theta > theta_thread)//|| theta <45)//设置若角度大于140说明不是拐角，可以merge /////150
			{
				flag = true;
				seg_1.push_back(p1);
				double stop_theta = 180;
				int i_tmp = i + 1;
				mypoint2f last_posi = p1_posi;
				while (stop_theta > 160 && i_tmp<pch_size - 1)///////////165
				{
					int p_tmp = pindex_patch[i_tmp];
					seg_1.push_back(p_tmp);
					mypoint2f ptmp_posi = sgraph.nodeList[p_tmp].position;
					int p_tmp_prior = pindex_patch[i_tmp + 1];
					mypoint2f prior_posi = sgraph.nodeList[p_tmp_prior].position;
					int next_adjpch = -1;
					EdgeOfPatch * nnedge = getBorderofPatch(pindex, p_tmp, p_tmp_prior);
					if (nnedge->adjPatchIndex.size() > 1)
						next_adjpch = nnedge->adjPatchIndex.back();
					if (next_adjpch != adj_pch_seg1)
						break;
					vector<int> sv1, sv2;
					sv1.push_back(pindex_patch[i_tmp]);
					sv1.push_back(pindex_patch[i_tmp - 1]);
					sv2.push_back(pindex_patch[i_tmp]);
					sv2.push_back(pindex_patch[i_tmp + 1]);
					vector<mypoint2f> sp1, sp2;
					jointMethod(&sv1, &sp1, 0);
					jointMethod(&sv2, &sp2, 0);
					double cos_st_theta = getCosineAngle(sp1[0], sp1[1], sp2[0], sp2[1]);//double cos_st_theta = getCosineAngle(ptmp_posi, prior_posi, ptmp_posi, last_posi);
					stop_theta = acos(cos_st_theta) * 180 / 3.1415926;
					last_posi = ptmp_posi;
					i_tmp++;
				}
				seg_2.push_back(p1);//逆序
				stop_theta = 180;
				i_tmp = pch_size - 2;
				last_posi = p1_posi;
				while (stop_theta > 160 && i_tmp >0) ///////////165
				{
					int p_tmp = pindex_patch[i_tmp];
					seg_2.push_back(p_tmp);
					mypoint2f ptmp_posi = sgraph.nodeList[p_tmp].position;
					int p_tmp_prior = pindex_patch[i_tmp - 1];
					mypoint2f prior_posi = sgraph.nodeList[p_tmp_prior].position;
					int next_adjpch = -1;
					EdgeOfPatch * nnedge = getBorderofPatch(pindex, p_tmp, p_tmp_prior);
					if (nnedge->adjPatchIndex.size() > 1)
						next_adjpch = nnedge->adjPatchIndex.back();
					if (next_adjpch != adj_pch_seg2)
						break;
					vector<int> sv1, sv2;
					sv1.push_back(pindex_patch[i_tmp]);
					sv1.push_back(pindex_patch[i_tmp - 1]);
					sv2.push_back(pindex_patch[i_tmp]);
					sv2.push_back(pindex_patch[i_tmp + 1]);
					vector<mypoint2f> sp1, sp2;
					jointMethod(&sv1, &sp1, 0);
					jointMethod(&sv2, &sp2, 0);
					double cos_st_theta = getCosineAngle(sp1[0], sp1[1], sp2[0], sp2[1]);//double cos_st_theta = getCosineAngle(ptmp_posi, prior_posi, ptmp_posi, last_posi);
					stop_theta = acos(cos_st_theta) * 180 / 3.1415926;
					last_posi = ptmp_posi;
					i_tmp--;
				}
			}
			else
			{
				flag = false;//非拐角，非平坦区域。接近直角？
			}
			break;
		}
	}
	if (flag == false)
	{
		return flag;
	}
	//根据theta分情况处理
	if (theta > theta_thread)//第一种处理
	{
		//2.判断seg_1或seg_2的方向，seg_1 seg_2 line都是从touchpoint出发的线段
		vector<mypoint2f> l_plist;
		vector<int> l_pvec = plrGraph.otherlines[lindex].pt_vec;//plrGraph.attachlines[lindex].pt_vec;
		if (l_pvec.front() == touchpoint)
			jointMethod(&l_pvec, &l_plist, 0);
		else
		{
			reverse(l_pvec.begin(), l_pvec.end());
			jointMethod(&l_pvec, &l_plist, 0);
		}
		double l_peri = getPerimeterofVec(&l_plist);
		vector<mypoint2f> seg_plist;  //被影响的edge
		vector<int> seg_pvec;  //被影响的edge
		seg_pvec = seg_1;
		jointMethod(&seg_1, &seg_plist, 0);
		double cos_theta1 = getCosineAngle(seg_plist.front(), seg_plist[1], l_plist[0], l_plist[1]);
		double theta1 = acos(cos_theta1) * 180 / 3.1415926;
		if (theta1 > (180 - theta_thread))/////180-140   edge1与line的角度>40
		{
			seg_plist.swap(vector<mypoint2f>());
			seg_pvec.swap(vector<int>());
			jointMethod(&seg_2, &seg_plist, 0);
			seg_pvec = seg_2;
			double cos_theta2 = getCosineAngle(seg_plist.front(), seg_plist[1], l_plist[0], l_plist[1]);
			double theta2 = acos(cos_theta2) * 180 / 3.1415926;
			//double theta2 = getCosineAngle(seg_plist.front(), seg_plist[2], l_plist[0], l_plist[1]) * 180 / 3.1415926;
			if (theta2 > (180 - theta_thread))//////edge2与line的角度>40
			{
				flag = false;//!!
				return flag;
			}
		}
		if (seg_plist.size() < 5 || l_plist.size()<5)//被影响的patch的边太小，不能合并patch+line////////////////从这往下不要？？？？
		{
			flag = false;//!!
			return flag;
		}
		//vector<mypoint2f> seg_plist;  vector<int> seg_pvec;  //受影响的edge
		//vector<mypoint2f> l_plist; 	vector<int> l_pvec;   //施加影响的line  
		vector<mypoint2f> new_plist;//被影响的edge ，变形之后的新线段
		double seg_peri = getPerimeterofVec(&seg_plist);

		new_plist.push_back(seg_plist.front());
		int record_i = 1;
		for (unsigned int j = 1; j<seg_plist.size() - 1; ++j)//遍历edge上的每一点（除了两端），寻找距离line的最近的3个点
		{
			double mini_dis = 65535;
			int mini_i = 0;
			int listend = l_plist.size() - 1;
			for (unsigned int i = record_i; i < listend; ++i)
			{
				double tmp_len = LineSegmentLength(seg_plist[j], l_plist[i]);
				if (mini_dis > tmp_len)
				{
					mini_dis = tmp_len;
					mini_i = i;
					(i == 0) ? record_i = i : record_i = i - 1;
					record_i = i - 1;
				}
				else
				{
					break;
				}
			}
			vector<mypoint2f> normal_dirc;
			normal_dirc.push_back(normalization(seg_plist[j], l_plist[mini_i]));
			normal_dirc.push_back(normalization(seg_plist[j], l_plist[mini_i + 1]));
			normal_dirc.push_back(normalization(seg_plist[j], l_plist[mini_i - 1]));
			normal_dirc.push_back(normalization(seg_plist[j], seg_plist[j + 1]));
			normal_dirc.push_back(normalization(seg_plist[j], seg_plist[j - 1]));
			vector<double> each_dis;
			each_dis.push_back(LineSegmentLength(seg_plist[j], l_plist[mini_i]));
			each_dis.push_back(LineSegmentLength(seg_plist[j], l_plist[mini_i + 1]));
			each_dis.push_back(LineSegmentLength(seg_plist[j], l_plist[mini_i - 1]));
			each_dis.push_back(LineSegmentLength(seg_plist[j], seg_plist[j + 1]));
			each_dis.push_back(LineSegmentLength(seg_plist[j], seg_plist[j - 1]));
			double disSum = 0;
			for (int k = 0; k < 5; ++k)//距离的倒数
			{
				disSum += 1 / each_dis[k];
			}
			for (int k = 0; k < 5; ++k)//归一化权值,并与单位向量相乘
			{
				each_dis[k] = (1 / each_dis[k]) / disSum;
				normal_dirc[k].x = normal_dirc[k].x*each_dis[k];
				normal_dirc[k].y = normal_dirc[k].y*each_dis[k];
			}
			mypoint2f addpt(0, 0);
			for (int k = 0; k < 5; ++k)
			{
				addpt = addpt + normal_dirc[k];
			}
			addpt = addpt + seg_plist[j];
			new_plist.push_back(addpt);
			//mypoint2f addpt = l_plist[mini_i] + l_plist[mini_i + 1] + l_plist[mini_i - 1] + seg_plist[j - 1] + seg_plist[j + 1];
			//addpt.x = addpt.x / 5;
			//addpt.y = addpt.y / 5;
			//new_plist.push_back(addpt);
		}
		new_plist.push_back(seg_plist.back());
		//new_plist与patch的其他edge 组成transform_pch（point list）
		if (seg_pvec[1] == pindex_patch[pch_size - 2])
		{
			reverse(pindex_patch.begin(), pindex_patch.end());
		}
		transform_pch->insert(transform_pch->end(), new_plist.begin(), new_plist.end());
		int j = 0;
		while (j<seg_pvec.size())
		{
			if (pindex_patch[j] == seg_pvec[j])
				j++;
			else
				break;
		}
		vector<int> left_pvec;
		for (int i = j - 1; i < pch_size; ++i)
		{
			left_pvec.push_back(pindex_patch[i]);
		}
		jointMethod(&left_pvec, transform_pch, 0);
	}
	return  flag;
}
void MyGraph::mergePatchLine(int pindex, int lindex)
{
	//pindex和lindex分别是plrGraph.patch/otherlines中的index
	//满足条件 才合并：接触点不是拐点
	//1、判断line和patch的位置
	vector<mypoint2f> origi_plist_line;
	jointMethod(&plrGraph.otherlines[lindex].pt_vec, &origi_plist_line, 0);
	pandl_combine_line.push_back(origi_plist_line);

	int touchpoint = -1;
	int pch_size = plrGraph.patches[pindex].pointIndexOfPatch.size();
	for (unsigned int i = 0; i < plrGraph.patches[pindex].point_line.size(); ++i)
	{
		if (plrGraph.patches[pindex].point_line[i].second == lindex)//lindex既可以是attachlines的index ，也是otherline 的index
		{
			touchpoint = plrGraph.patches[pindex].point_line[i].first;
			break;
		}
	}
	bool flag = false;
	if (touchpoint < 0)
		return;
	vector<int> seg_1;  //若连接处不是拐点，则记录下该点两边的线
	vector<int> seg_2;
	double theta_thread = 140;//广角>140 ， line与patch 的edge的角度<40度
	double theta = 0;
	vector<int> pindex_patch = plrGraph.patches[pindex].pointIndexOfPatch;
	for (int i = 0; i < pch_size; ++i)
	{
		if (touchpoint == pindex_patch[i])
		{
			if (i > 0 && i < (pch_size - 1))//若touch 点不在两端,调整到两端
			{
				if (i >(pch_size / 2))
				{
					//拆后面的，补到前面
					int stop_index = pindex_patch[i];
					pindex_patch.pop_back();
					vector<int>::iterator iter_pch = pindex_patch.begin();
					while (pindex_patch.back() != stop_index)
					{
						int pop_index = pindex_patch.back();
						pindex_patch.pop_back();
						iter_pch = pindex_patch.insert(iter_pch, pop_index);
					}
					iter_pch = pindex_patch.insert(iter_pch, stop_index);
				}
				else
				{
					//删除前面的，push到后面
					int stop_index = pindex_patch[i];
					vector<int>::iterator iter_pch = pindex_patch.begin();
					iter_pch = pindex_patch.erase(iter_pch);///????
					while (pindex_patch.front() != stop_index)
					{
						int erase_index = pindex_patch.front();
						iter_pch = pindex_patch.erase(iter_pch);
						pindex_patch.push_back(erase_index);
					}
					pindex_patch.push_back(stop_index);
				}
				i = 0;
			}
			int p1 = pindex_patch[i];
			mypoint2f p1_posi = sgraph.nodeList[p1].position;
			int p2 = pindex_patch[i + 1];
			int p3 = pindex_patch[pch_size - 2];
			vector<int> ttpvec;
			vector<mypoint2f> ttplist;
			ttpvec.push_back(p1); ttpvec.push_back(p2);
			jointMethod(&ttpvec, &ttplist, 0);
			mypoint2f p2_posi = ttplist[1];
			ttpvec[1] = p3;
			ttplist.swap(vector<mypoint2f>());
			jointMethod(&ttpvec, &ttplist, 0);
			mypoint2f p3_posi = ttplist[1];
			double costheta = getCosineAngle(p1_posi, p2_posi, p1_posi, p3_posi);
			theta = acos(costheta) * 180 / 3.1415926;
			int adj_pch_seg1 = -1;
			int adj_pch_seg2 = -1;
			EdgeOfPatch * anedge = getBorderofPatch(pindex, p1, p2);
			if (anedge->adjPatchIndex.size() > 1)
				adj_pch_seg1 = anedge->adjPatchIndex.back();
			anedge = getBorderofPatch(pindex, p1, p3);
			if (anedge->adjPatchIndex.size() > 1)
				adj_pch_seg2 = anedge->adjPatchIndex.back();
			if (theta > theta_thread)// || theta <45)//设置若角度大于150说明不是拐角，可以merge
			{
				flag = true;
				seg_1.push_back(p1);
				double stop_theta = 180;
				int i_tmp = i + 1;
				mypoint2f last_posi = p1_posi;
				while (stop_theta > 160 && i_tmp<pch_size - 1)
				{
					int p_tmp = pindex_patch[i_tmp];
					seg_1.push_back(p_tmp);
					mypoint2f ptmp_posi = sgraph.nodeList[p_tmp].position;
					int p_tmp_prior = pindex_patch[i_tmp + 1];
					mypoint2f prior_posi = sgraph.nodeList[p_tmp_prior].position;
					int next_adjpch = -1;
					EdgeOfPatch * nnedge = getBorderofPatch(pindex, p_tmp, p_tmp_prior);
					if (nnedge->adjPatchIndex.size() > 1)
						next_adjpch = nnedge->adjPatchIndex.back();
					if (next_adjpch != adj_pch_seg1)
						break;
					vector<int> sv1, sv2;
					sv1.push_back(pindex_patch[i_tmp]);
					sv1.push_back(pindex_patch[i_tmp - 1]);
					sv2.push_back(pindex_patch[i_tmp]);
					sv2.push_back(pindex_patch[i_tmp + 1]);
					vector<mypoint2f> sp1, sp2;
					jointMethod(&sv1, &sp1, 0);
					jointMethod(&sv2, &sp2, 0);
					double cos_st_theta = getCosineAngle(sp1[0], sp1[1], sp2[0], sp2[1]);//double cos_st_theta = getCosineAngle(ptmp_posi, prior_posi, ptmp_posi, last_posi);
					stop_theta = acos(cos_st_theta) * 180 / 3.1415926;
					last_posi = ptmp_posi;
					i_tmp++;
				}
				seg_2.push_back(p1);//逆序
				stop_theta = 180;
				i_tmp = pch_size - 2;
				last_posi = p1_posi;
				while (stop_theta > 160 && i_tmp >0)
				{
					int p_tmp = pindex_patch[i_tmp];
					seg_2.push_back(p_tmp);
					mypoint2f ptmp_posi = sgraph.nodeList[p_tmp].position;
					int p_tmp_prior = pindex_patch[i_tmp - 1];
					mypoint2f prior_posi = sgraph.nodeList[p_tmp_prior].position;
					int next_adjpch = -1;
					EdgeOfPatch * nnedge = getBorderofPatch(pindex, p_tmp, p_tmp_prior);
					if (nnedge->adjPatchIndex.size() > 1)
						next_adjpch = nnedge->adjPatchIndex.back();
					if (next_adjpch != adj_pch_seg2)
						break;
					vector<int> sv1, sv2;
					sv1.push_back(pindex_patch[i_tmp]);
					sv1.push_back(pindex_patch[i_tmp - 1]);
					sv2.push_back(pindex_patch[i_tmp]);
					sv2.push_back(pindex_patch[i_tmp + 1]);
					vector<mypoint2f> sp1, sp2;
					jointMethod(&sv1, &sp1, 0);
					jointMethod(&sv2, &sp2, 0);
					double cos_st_theta = getCosineAngle(sp1[0], sp1[1], sp2[0], sp2[1]);//double cos_st_theta = getCosineAngle(ptmp_posi, prior_posi, ptmp_posi, last_posi);
					stop_theta = acos(cos_st_theta) * 180 / 3.1415926;
					//stop_theta = getCosineAngle(ptmp_posi, prior_posi, ptmp_posi, last_posi) * 180 / 3.1415926;
					last_posi = ptmp_posi;
					i_tmp--;
				}
			}
			else
			{
				flag = false;//非拐角，非平坦区域。接近直角？
			}
			break;
		}
	}
	if (flag == false)
	{
		return;
	}
	//根据theta分情况处理
	vector<mypoint2f> transform_pch;  //变形后的newpatch  的 pointlist
	vector<mypoint2f> new_plist;     //变形的那条线的plist
	vector<mypoint2f> l_plist;
	vector<int> l_pvec = plrGraph.otherlines[lindex].pt_vec;//plrGraph.attachlines[lindex].pt_vec;
	if (l_pvec.front() == touchpoint)
		jointMethod(&l_pvec, &l_plist, 0);
	else
	{
		reverse(l_pvec.begin(), l_pvec.end());
		jointMethod(&l_pvec, &l_plist, 0);
	}
	double l_peri = getPerimeterofVec(&l_plist);

	vector<mypoint2f> seg_plist;
	vector<int> seg_pvec;
	double seg_peri = 0;
	if (theta > theta_thread)//第一种处理
	{
		//2.判断seg_1或seg_2的方向		
		seg_pvec = seg_1;
		jointMethod(&seg_1, &seg_plist, 0);
		double cos_theta1 = getCosineAngle(seg_plist.front(), seg_plist[1], l_plist[0], l_plist[1]);
		double theta1 = acos(cos_theta1) * 180 / 3.1415926;
		if (theta1 > (180 - theta_thread))
		{
			seg_plist.swap(vector<mypoint2f>());
			seg_pvec.swap(vector<int>());
			jointMethod(&seg_2, &seg_plist, 0);
			seg_pvec = seg_2;
			double cos_theta2 = getCosineAngle(seg_plist.front(), seg_plist[1], l_plist[0], l_plist[1]);
			double theta2 = acos(cos_theta2) * 180 / 3.1415926;
			if (theta2 > (180 - theta_thread))
			{
				flag = false;
				return;
			}
		}
		//vector<mypoint2f> seg_plist;  vector<int> seg_pvec;  //受影响的edge
		//vector<mypoint2f> l_plist; 	vector<int> l_pvec;   //施加影响的line  
		//vector<mypoint2f> new_plist;//被影响的edge ，变形之后的新线段
		seg_peri = getPerimeterofVec(&seg_plist);
		new_plist.push_back(seg_plist.front());
		int record_i = 1;
		for (unsigned int j = 1; j<seg_plist.size() - 1; ++j)//遍历edge上的每一点（除了两端），寻找距离line的最近的3个点
		{
			double mini_dis = 65535;
			int mini_i = 0;
			int listend = l_plist.size() - 1;
			for (unsigned int i = record_i; i < listend; ++i)
			{
				double tmp_len = LineSegmentLength(seg_plist[j], l_plist[i]);
				if (mini_dis > tmp_len)
				{
					mini_dis = tmp_len;
					mini_i = i;
					record_i = i - 1;
				}
				else
				{
					break;
				}
			}
			vector<mypoint2f> normal_dirc;
			normal_dirc.push_back(normalization(seg_plist[j], l_plist[mini_i]));
			normal_dirc.push_back(normalization(seg_plist[j], l_plist[mini_i + 1]));
			normal_dirc.push_back(normalization(seg_plist[j], l_plist[mini_i - 1]));
			normal_dirc.push_back(normalization(seg_plist[j], seg_plist[j + 1]));
			normal_dirc.push_back(normalization(seg_plist[j], seg_plist[j - 1]));
			vector<double> each_dis;
			each_dis.push_back(LineSegmentLength(seg_plist[j], l_plist[mini_i]));
			each_dis.push_back(LineSegmentLength(seg_plist[j], l_plist[mini_i + 1]));
			each_dis.push_back(LineSegmentLength(seg_plist[j], l_plist[mini_i - 1]));
			each_dis.push_back(LineSegmentLength(seg_plist[j], seg_plist[j + 1]));
			each_dis.push_back(LineSegmentLength(seg_plist[j], seg_plist[j - 1]));
			double disSum = 0;
			for (int k = 0; k < 5; ++k)//距离的倒数
			{
				disSum += 1 / each_dis[k];
			}
			for (int k = 0; k < 5; ++k)//归一化权值,并与单位向量相乘
			{
				each_dis[k] = (1 / each_dis[k]) / disSum;
				normal_dirc[k].x = normal_dirc[k].x*each_dis[k];
				normal_dirc[k].y = normal_dirc[k].y*each_dis[k];
			}
			mypoint2f addpt(0, 0);
			for (int k = 0; k < 5; ++k)
			{
				addpt = addpt + normal_dirc[k];
			}
			addpt = addpt + seg_plist[j];
			new_plist.push_back(addpt);
			//mypoint2f addpt = l_plist[mini_i] + l_plist[mini_i + 1] + l_plist[mini_i - 1] + seg_plist[j - 1] + seg_plist[j + 1];
			//addpt.x = addpt.x / 5;
			//addpt.y = addpt.y / 5;
			//new_plist.push_back(addpt);
		}
		new_plist.push_back(seg_plist.back());
		//new_plist与patch的其他edge 组成transform_pch（point list）
		if (seg_pvec[1] == pindex_patch[pch_size - 2])
		{
			reverse(pindex_patch.begin(), pindex_patch.end());
		}
		transform_pch.insert(transform_pch.end(), new_plist.begin(), new_plist.end());
		int j = 0;
		while (j<seg_pvec.size())
		{
			if (pindex_patch[j] == seg_pvec[j])
				j++;
			else
				break;
		}
		vector<int> left_pvec;
		for (int i = j - 1; i < pch_size; ++i)
		{
			left_pvec.push_back(pindex_patch[i]);
		}
		jointMethod(&left_pvec, &transform_pch, 0);// transform_pch done
												   //删除sgraph中的 arc//????
		vector<int> rm_line = plrGraph.otherlines[lindex].pt_vec;
		for (size_t ti = 1; ti < rm_line.size(); ++ti)
		{
			removeArcNode(rm_line[ti - 1], rm_line[ti], 0);
			removeArcNode(rm_line[ti], rm_line[ti - 1], 0);
		}
		//修改sgraph.nodeList[。].firstEdge->curveType/curveIndex;  polyline_vec[.].pointlist;
		if (seg_pvec.size() > 2)
		{
			vector<double> peri_ratio;
			vector<mypoint2f> sub_plist;
			for (unsigned int i = 1; i < seg_pvec.size(); ++i)
			{
				int sp1 = seg_pvec[i - 1];
				int sp2 = seg_pvec[i];
				mypoint2f posi1 = sgraph.nodeList[sp1].position;
				mypoint2f posi2 = sgraph.nodeList[sp2].position;
				if (!sub_plist.empty())
					sub_plist.pop_back();
				getPListOfaCurve(posi1, posi2, sp1, sp2, &sub_plist, 0);
				double sub_peri = getPerimeterofVec(&sub_plist);
				peri_ratio.push_back(sub_peri / seg_peri);
			}
			int last_subi = 0;
			for (unsigned int i = 0; i < peri_ratio.size(); ++i)//seg_pvec有3个，peri_ratio有2两个
			{
				int sub_i = new_plist.size()*peri_ratio[i] - 1;  //sub_i把new_plist按比例大体分成几段
				if (sub_i <= 0)
					sub_i = 1;
				vector<mypoint2f> modi_plist;
				for (int j = last_subi; j <= sub_i; ++j)
				{
					modi_plist.push_back(new_plist[j]);
				}
				ArcNode* anarc = getArcNode(seg_pvec[i], seg_pvec[i + 1]);
				exchangeline(&modi_plist, anarc);//polyline_vec中替换新的line：modi_plist
				last_subi = sub_i;
				////原来，直接修改点的位置
				//sgraph.nodeList[seg_pvec[i + 1]].position = new_plist[sub_i];/////!!!!!!!
				////现在，patch 中point_line中有seg_pvec[i+1]点所连接的线，整体平移
				mypoint2f move_dic = new_plist[sub_i] - sgraph.nodeList[seg_pvec[i + 1]].position;
				sgraph.nodeList[seg_pvec[i + 1]].position = new_plist[sub_i];
				if (!plrGraph.patches[pindex].point_line.empty())
				{
					for (unsigned int lli = 0; lli < plrGraph.patches[pindex].point_line.size(); ++lli)
					{
						if (plrGraph.patches[pindex].point_line[lli].first == seg_pvec[i + 1])
						{
							int move_lindex = plrGraph.patches[pindex].point_line[lli].second;
							vector<int> move_lvec = plrGraph.otherlines[move_lindex].pt_vec;
							for (unsigned int k = 1; k < move_lvec.size(); ++k)
							{
								ArcNode* llarc = getArcNode(move_lvec[k - 1], move_lvec[k]);
								vector<int> tmplvec;
								tmplvec.push_back(move_lvec[k - 1]);
								tmplvec.push_back(move_lvec[k]);
								vector<mypoint2f> move_plist;
								jointMethod(&tmplvec, &move_plist, 0);
								for (unsigned int g = 0; g < move_plist.size(); ++g)
									move_plist[g] = move_plist[g] + move_dic;
								exchangeline(&move_plist, llarc);
							}
						}
					}
				}

			}
		}
		else
		{
			ArcNode* anarc = getArcNode(seg_pvec[0], seg_pvec[1]);
			exchangeline(&new_plist, anarc);
		}
	}
	if (flag == true)
	{
		//patch 中删除point_line
		vector<pair<int, int>>::iterator iter_pl = plrGraph.patches[pindex].point_line.begin();
		while (iter_pl != plrGraph.patches[pindex].point_line.end())
		{
			if (iter_pl->first == touchpoint && iter_pl->second == lindex)
				iter_pl = plrGraph.patches[pindex].point_line.erase(iter_pl);
			else
				iter_pl++;
		}
		//对attachline/otherline的处理
		plrGraph.otherlines[lindex].type = -1;  //	plrGraph.attachlines[lindex].type = -1;

												//计算面积 周长 中心点
		double new_area = GetPolygonArea(&transform_pch);
		double new_peri = getPerimeterofVec(&transform_pch);
		mypoint2f sum_pt(0, 0);
		for (unsigned int i = 0; i < transform_pch.size(); ++i)
		{
			sum_pt = sum_pt + transform_pch[i];
		}
		sum_pt.x = sum_pt.x / transform_pch.size();
		sum_pt.y = sum_pt.y / transform_pch.size();
		plrGraph.patches[pindex].area = new_area;
		plrGraph.patches[pindex].perimeter = new_peri;
		plrGraph.patches[pindex].centre_position = sum_pt;
		//测试
		pandl_combine.push_back(transform_pch);//////////////////////////////
											   /*
											   sgraph.nodeList[0].firstEdge->curveType; sgraph.nodeList[0].firstEdge->curveIndex; polyline_vec[0].pointlist;
											   plrGraph.patches[0].area;
											   plrGraph.patches[0].centre_position;
											   plrGraph.patches[0].perimeter;
											   plrGraph.patches[0].point_line.erase();
											   plrGraph.attachlines[0].type = -1;
											   plrGraph.otherlines[0].type = -1;*/
	}
}
void MyGraph::exchangeline(vector<mypoint2f>* nplist, ArcNode* exc_node)
{
	int modi_lindex = exc_node->curveIndex;
	char modi_ltype = exc_node->curveType;
	if (modi_ltype == 'L')//?????
	{
		polyline_vec[modi_lindex].pointlist.swap(vector<mypoint2f>());
		polyline_vec[modi_lindex].pointlist = *nplist;/////!!!!!!!!
	}
	else if (modi_ltype == 'C' || modi_ltype == 'Q')
	{
		exc_node->curveType = 'L';
		PolyLine pll;
		pll.origi_index = polyline_vec.size();
		pll.origi_type = 'L';
		pll.pointlist = *nplist;
		polyline_vec.push_back(pll);
		exc_node->curveIndex = polyline_vec.size() - 1;
	}
}


void MyGraph::globalSmoothing()
{
	//check mixedTopology中A是不是有B B是不是有A
	int mixsize = mixedTopology.allNode.size();
	int pchsize = plrGraph.patches.size();
	int lsize = plrGraph.otherlines.size();
	for (int i = 0; i < mixsize; ++i)
	{
		if (mixedTopology.allNode[i].status >= 0)
		{
			if (mixedTopology.allNode[i].pl_type != 'p')  //line
			{
				int ml_index = mixedTopology.allNode[i].pl_index - pchsize;
				int ss = plrGraph.otherlines[ml_index].type;
				if (ss < 0)//mix node 的status=-1，单对应的patch/line 的type！=-1
				{
					cout << "mixedTopology.allNode[i]" << mixedTopology.allNode[i].pl_index << endl;
					mixedTopology.allNode[i].status = -1;
					MixedArc* mmarc = mixedTopology.allNode[i].firstAdjNode;
					while (mmarc != NULL)
					{
						int neiarc = mmarc->adj_index;
						mmarc = mmarc->next_arc;
						removeMixedNode(i, neiarc);
						removeMixedNode(neiarc, i);
					}
					vector<int> rm_line = plrGraph.otherlines[ml_index].pt_vec;
					for (size_t ti = 1; ti < rm_line.size(); ++ti)
					{
						removeArcNode(rm_line[ti - 1], rm_line[ti], 0);
						removeArcNode(rm_line[ti], rm_line[ti - 1], 0);
					}
				}
			}
			else  //patch
			{
				int ss = plrGraph.patches[i].type;
				if (ss < 0)//mix node 的status=-1，单对应的patch/line 的type！=-1
				{
					cout << "error" << endl;
				}
			}
		}
	}	
	for (int i = 0; i < pchsize; ++i)
	{
		if (plrGraph.patches[i].type >= 0)
		{
			if (mixedTopology.allNode[i].status < 0)
				cout << "error" << endl;
			vector<int> pch = plrGraph.patches[i].pointIndexOfPatch;
			for (size_t j = 1; j < pch.size(); ++j)
			{
				ArcNode* narc = getArcNode(pch[j - 1], pch[j]);
				if (narc == NULL)
				{
					cout << "error" << endl;
				}
				narc = getArcNode(pch[j], pch[j - 1]);
				if (narc == NULL)
				{
					cout << "error" << endl;
				}
			}
		}
	}
	for (int i = 0; i < lsize; ++i)
	{
		if (plrGraph.otherlines[i].type >= 0)
		{
			if(mixedTopology.allNode[i + pchsize].status < 0)
				cout << "error" << endl;
			vector<int> pch = plrGraph.otherlines[i].pt_vec;
			for (size_t j = 1; j < pch.size(); ++j)
			{
				ArcNode* narc = getArcNode(pch[j - 1], pch[j]);
				if (narc == NULL)
				{
					cout << "error" << endl;
				}
				narc = getArcNode(pch[j], pch[j - 1]);
				if (narc == NULL)
				{
					cout << "error" << endl;
				}
			}
		}
	}

	////遍历每个patch， otherline， innerline
//	vector<cubicBezier> total_pw_cbcurve;
	double fit_error = 5;
	vector<int> visited((int)sgraph.nodeList.size(),0);
	for(int i = 0; i < pchsize;++i)
	{
		if (plrGraph.patches[i].type >= 0)
		{
			vector<int> pch_vec_adjust;	
			bool degree_flag = false;  //若所有的点度=2，degree_flag=false;
			if (countEdgeNum(plrGraph.patches[i].pointIndexOfPatch[0], 0) == 2)//第一个点度=2，调整一下 另度为3的点做第一个点
			{
				vector<int> pch_vec = plrGraph.patches[i].pointIndexOfPatch;
				int store_j = -1;
				vector<int> store_vec;
				for (size_t j = 1; j < pch_vec.size(); ++j)
				{
					int edgenum = countEdgeNum(pch_vec[j], 0);
					if (edgenum == 2)
					{
						store_vec.push_back(pch_vec[j]);
					}
					else
					{
						degree_flag = true;
						store_j = j;
						break;
					}
				}
				if (store_j >= 0)
				{
					for (int k = store_j; k < pch_vec.size(); ++k)
					{
						pch_vec_adjust.push_back(pch_vec[k]);
					}
					if (!store_vec.empty())
						pch_vec_adjust.insert(pch_vec_adjust.end(), store_vec.begin(), store_vec.end());
					pch_vec_adjust.push_back(pch_vec[store_j]);
				}		
			}
			else
			{
				pch_vec_adjust = plrGraph.patches[i].pointIndexOfPatch;
			}


			if (degree_flag == false)//所有的点度 = 2，
			{
				pch_vec_adjust = plrGraph.patches[i].pointIndexOfPatch;
				vector<cubicBezier> tmp_cb;
				piecewiseCurveFittting(&pch_vec_adjust,fit_error,&tmp_cb);//////////////////
				total_pw_cbcurve.push_back(tmp_cb);
			}
			else
			{
				vector<int> seg_vec;
				seg_vec.push_back(pch_vec_adjust[0]);
				for (size_t j = 1; j < pch_vec_adjust.size(); ++j)
				{
					if (countEdgeNum(pch_vec_adjust[j], 0) !=2)
					{
						seg_vec.push_back(pch_vec_adjust[j]);
					//	visited[pch_vec_adjust[j]] == 1;
						if (seg_vec.size() > 1)
						{
							vector<cubicBezier> tmp_cb;
							piecewiseCurveFittting(&seg_vec, fit_error, &tmp_cb);//////////////////
							total_pw_cbcurve.push_back(tmp_cb);
						}
						seg_vec.swap(vector<int>());
						seg_vec.push_back(pch_vec_adjust[j]);
					}
					else
					{
						if (visited[pch_vec_adjust[j]] == 0)
						{
							if(seg_vec.empty())
								seg_vec.push_back(pch_vec_adjust[j-1]);
							seg_vec.push_back(pch_vec_adjust[j]);
							visited[pch_vec_adjust[j]] == 1;
						}
						else
						{
							seg_vec.push_back(pch_vec_adjust[j]);
							visited[pch_vec_adjust[j]] == 1;
							if (seg_vec.size() > 1)
							{
								vector<cubicBezier> tmp_cb;
								piecewiseCurveFittting(&seg_vec, fit_error, &tmp_cb);//////////////////
								total_pw_cbcurve.push_back(tmp_cb);
							}
							seg_vec.swap(vector<int>());
						}
					}
					
				}
			}
		}
	}
	for (int i = 0; i < lsize; ++i)
	{
		if (plrGraph.otherlines[i].type >= 0)
		{
			vector<int> lvec = plrGraph.otherlines[i].pt_vec;
			//visited[lvec[j]] == 1;
			vector<cubicBezier> tmp_cb;
			piecewiseCurveFittting(&lvec, fit_error, &tmp_cb);//////////////////
			total_pw_cbcurve.push_back(tmp_cb);
		}
	}
	int il_size = plrGraph.innerlines.size();
	if (il_size > 0)
	{
		for (size_t i = 0; i < il_size; ++i)
		{
			vector<int> lvec = plrGraph.innerlines[i].pt_vec;
			//visited[lvec[j]] == 1;
			vector<cubicBezier> tmp_cb;
			piecewiseCurveFittting(&lvec, fit_error, &tmp_cb);//////////////////
			total_pw_cbcurve.push_back(tmp_cb);
		}
	}
}
void MyGraph::piecewiseCurveFittting(vector<int>* pt_vec, double error, vector<cubicBezier>* pw_curve)
{
	if (pt_vec->size() <= 2)
	{
		vector<mypoint2f> ptdata;
		jointMethod(pt_vec, &ptdata, 0);
		FitCurves(&ptdata, ptdata.size(), error, pw_curve);
	}
	else
	{
		int starti = 1;
		vector<vector<int>> all_segments;
		int vecsize = pt_vec->size();
		while (starti < vecsize - 1)
		{
			vector<int> seg_vec;
			for (int i = starti; i < vecsize - 1; ++i)
			{
				vector<int> tmp1;
				vector<mypoint2f> ptdata1;
				tmp1.push_back(pt_vec->at(i - 1));
				tmp1.push_back(pt_vec->at(i));
				jointMethod(&tmp1, &ptdata1, 0);
				vector<int> tmp2;
				vector<mypoint2f> ptdata2;
				tmp2.push_back(pt_vec->at(i));
				tmp2.push_back(pt_vec->at(i + 1));
				jointMethod(&tmp2, &ptdata2, 0);

				double bend_cos = getCosineAngle(ptdata2[0], ptdata2[1], ptdata1.back(), ptdata1[ptdata1.size() - 2]);
				double bend_angle = acos(bend_cos) * 180 / 3.1415926;
				if (bend_angle < 135)
				{
					for (int j = starti - 1; j <= i; ++j)
					{
						seg_vec.push_back(pt_vec->at(j));
					}
					all_segments.push_back(seg_vec);
					starti = i + 1;//!!!!!!!!!
					break;
				}
			}
			if (seg_vec.empty())
			{
				for (int j = starti - 1; j <vecsize; ++j)
				{
					seg_vec.push_back(pt_vec->at(j));
				}
				all_segments.push_back(seg_vec);
				break;
			}
		}
		for (size_t i = 0; i < all_segments.size(); ++i)
		{
			vector<int> tmpseg = all_segments[i];
			vector<mypoint2f> ptdata;
			jointMethod(&tmpseg, &ptdata, 0);
			FitCurves(&ptdata, ptdata.size(), error, pw_curve);
		}
	}

	////原来的，没有在拐点分段
	//vector<mypoint2f> ptdata;
	//jointMethod(pt_vec, &ptdata, 0);
	//FitCurves(&ptdata, ptdata.size(), error, pw_curve);
}
void MyGraph::clearStructure()
{
	//清空原来的数据结构
	cubicBezier_vec_origin.swap(vector<cubicBezier>());
	quadraticBezier_vec_origin.swap(vector<quaBezier>());
	polyline_vec_origin.swap(vector<PolyLine>());

	cubicBezier_vec.swap(vector<cubicBezier>());
	quadraticBezier_vec.swap(vector<quaBezier>());
	polyline_vec.swap(vector<PolyLine>());
	added_polyline.swap(vector<PolyLine>());

	pointTopl.swap(vector<pair<mypoint2f, vector<TerminalPoint>>>());
	myEdgeGraph.edgeNumber = 0;
	myEdgeGraph.edgesOfGraph.swap(vector<EdgeStruct>());
	sgraph.arcNum = 0;
	sgraph.vertexNum = 0;
	sgraph.nodeList.swap(vector<VertexNode>());
	subgraph_map.swap(map<int, string>());
	plrGraph.attachlines.swap(vector<ICurve>());
	plrGraph.innerlines.swap(vector<ICurve>());
	plrGraph.isolatedlines.swap(vector<ICurve>());
	plrGraph.otherlines.swap(vector<ICurve>());
	plrGraph.patches.swap(vector<Patch>());

	pointListOfPatches.swap(vector<pair<int, vector<mypoint2f>>>());
	pointListOfLines.swap(vector<PolyLine>());

	patchfeatures.swap(vector<pair<int, vector<double>>>());
	mixedTopology.mixnode_number = 0;
	mixedTopology.allNode.swap(vector<MixedNode>());
	//赋初值
	for (size_t i = 0; i < total_pw_cbcurve.size(); ++i)
	{
		cubicBezier_vec_origin.insert(cubicBezier_vec_origin.end(), total_pw_cbcurve[i].begin(), total_pw_cbcurve[i].end());
	}
	total_pw_cbcurve.swap(vector<vector<cubicBezier>>());

	/*
	vector<pair<mypoint2f, vector<TerminalPoint>>> pointTopl;
	EdgeGraph myEdgeGraph;
	SvgGraph sgraph;
	static double errorDeviation;//没有用到，，删
	map<int, string> subgraph_map;//record tree's or loop's start index
	PLRGraph plrGraph;
	vector<pair<int, vector<mypoint2f>>> pointListOfPatches;
	vector<PolyLine> pointListOfLines;

	vector<pair<int, string>> patch_subgraphs;
	PatchGraph patchTopology;
	vector<pair<int, string>> patch_subgraphs_group;
	vector<pair<int, vector<double>>> patchfeatures;//int对应patchindex，  vector<double>对应特征向量

	vector<PolyLine> neighbor_plist_group;//测试 获取一条线的neighborline，借助topology，在getNeiborLines()中用到
	vector<pair<int, vector<mypoint2f>>> patch_simplified_plist;//测试，所有patch具有固定数量的点的plist
	vector<vector<mypoint2f>> pandl_combine;
	vector<vector<mypoint2f>> pandl_combine_line;

	MixedGraph mixedTopology;
	vector<vector<cubicBezier>> total_pw_cbcurve;
	*/
}


//输出svgfile 着色
void MyGraph::writeSVGFile(string fname, string outpath, double fwidth, double fheight, vector<double> fview, string view_str)
{
	/*vector<pair<int,vector<mypoint2f>>> pointListOfPatches;
	vector<PolyLine> pointListOfLines;*/
	TiXmlDocument *writeDoc = new TiXmlDocument(); //xml文档指针  
	//XML 文件声明。 standalone=no 属性:该属性规定此 SVG 文件是否是“独立的”，或含有对外部文件的引用。 
	//standalone = "no" 意味着 SVG 文档会引用一个外部文件-->在这里，是 DTD 文件。
	TiXmlDeclaration *decl = new TiXmlDeclaration("1.0", "UTF-8", "no");
	writeDoc->LinkEndChild(decl); //写入文档  
	//有错
	//TiXmlElement *refer_file = new TiXmlElement("!DOCTYPE");
	//refer_file->SetAttribute("svg PUBLIC", "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd");
	//writeDoc->LinkEndChild(refer_file);
	//根元素svg
	TiXmlElement *rootele = new TiXmlElement("svg");
	rootele->SetAttribute("version", "1.1");
	rootele->SetAttribute("xmlns", "http://www.w3.org/2000/svg");
	rootele->SetAttribute("xmlns:xlink", "http://www.w3.org/1999/xlink");
	rootele->SetAttribute("width",fwidth);
	rootele->SetAttribute("height", fheight);
	rootele->SetAttribute("viewBox", view_str.c_str());
	rootele->SetAttribute("xml:space", "preserve");
	writeDoc->LinkEndChild(rootele);
	//一个path是一个多边形面
	//例子：
	//注意 fill-rule属性：nonzero, evenodd。 style = "fill:white;stroke:red;stroke-width:2"
	/* <path d="M 100 100 L 300 100 L 200 300 z" fill="red" stroke="blue" stroke-width="3" />*/
	/*<polygon points="220,100 300,210 170,250" style="fill:#cccccc;stroke:#000000;stroke-width:1"/>*/
	vector<char*> color_vec;
	char *c1 = "ff0000"; color_vec.push_back(c1);
	char *c2 = "00ff00"; color_vec.push_back(c2);
	char *c3 = "0000ff"; color_vec.push_back(c3);
	char *c4 = "ffff00"; color_vec.push_back(c4);
	char *c5 = "00ffff"; color_vec.push_back(c5);
	char *c6 = "ff00ff"; color_vec.push_back(c6);
	char *c7 = "b400a0"; color_vec.push_back(c7);
	char *c8 = "ff7878"; color_vec.push_back(c8);
	char *c9 = "c8ff00"; color_vec.push_back(c9);	
	char *c10 = "ff96ff"; color_vec.push_back(c10);
	char *c11 = "3cc8ff"; color_vec.push_back(c11);
	char *c12 = "ffa000"; color_vec.push_back(c12);
	char *c13 = "9655e6"; color_vec.push_back(c13);
	char *c14 = "509600"; color_vec.push_back(c14);
	char *c15 = "fff0c8"; color_vec.push_back(c15);
	char *c16 = "9acd32"; color_vec.push_back(c16);
	char *c17 = "ff6eb4"; color_vec.push_back(c17);
	char *c18 = "4169ff"; color_vec.push_back(c18);
	int patch_count = 0;
	for (int i = 0; i < pointListOfPatches.size(); ++i)     //patch_simplified_plist   pointListOfPatches
	{
		if (pointListOfPatches[i].first >= 0)//不能=-1
		{
			TiXmlElement* apath = new TiXmlElement("polygon");
			//apath->SetAttribute("style", "fill:#aaaaaa;stroke:#000000;stroke-width:2;stroke:red");
			//apath->SetAttribute("points", "110,100 200,210 70,250");
			vector<mypoint2f> p_vec = pointListOfPatches[i].second;
			int psize = p_vec.size();
			string path_str;       //points="220,100 300,210 170,250"
			for (unsigned int j = 0; j < psize; ++j)
			{
				path_str.append(to_string(p_vec[j].x));
				path_str.append(",");
				path_str.append(to_string(ginstance->file_height-p_vec[j].y)); // 转置！！！！！
				path_str.append(" ");
			}
			path_str.pop_back();
			apath->SetAttribute("points", path_str.c_str());
			//string style_str = "stroke-opacity:1;stroke-width:1;stroke:#606060;fill:#a0a0a0";//灰黑边，内灰
			string style_str = "stroke-opacity:1;stroke-width:1;stroke:#";//stroke的颜色和patch保持一致
			int cnum = patch_count % 18;
			style_str.append(color_vec[cnum]);
			style_str.append(";fill:#");
			style_str.append(color_vec[cnum]);
			apath->SetAttribute("style", style_str.c_str());
			rootele->LinkEndChild(apath);
			patch_count++;
		}
	}

	//for (unsigned int i = 0; i < pandl_combine.size(); ++i)
	//{
	//	TiXmlElement* apath = new TiXmlElement("polygon");
	//	//apath->SetAttribute("style", "fill:#aaaaaa;stroke:#000000;stroke-width:2;stroke:red");
	//	//apath->SetAttribute("points", "110,100 200,210 70,250");
	//	vector<mypoint2f> p_vec = pandl_combine[i];
	//	int psize = p_vec.size();
	//	string path_str;       //points="220,100 300,210 170,250"
	//	for (unsigned int j = 0; j < psize; ++j)
	//	{
	//		path_str.append(to_string(p_vec[j].x));
	//		path_str.append(",");
	//		path_str.append(to_string(ginstance->file_height - p_vec[j].y)); // 转置！！！！！
	//		path_str.append(" ");
	//	}
	//	path_str.pop_back();
	//	apath->SetAttribute("points", path_str.c_str());
	//	string style_str = "stroke-opacity:1;stroke-width:1;stroke:#494949;fill:#ff0000";//灰黑边，内红
	//	apath->SetAttribute("style", style_str.c_str());
	//	rootele->LinkEndChild(apath);
	//}

	for (unsigned int i = 0; i < pointListOfLines.size(); ++i)
	{
		TiXmlElement* apath = new TiXmlElement("path");
		apath->SetAttribute("style", "fill:none;stroke-width:1;stroke:#606060");//线条是灰色
		//apath->SetAttribute("d", "220,100 300,210 170,250");
		//polyline 中点的转换
		string path_str;       //d="M 1062.5,587.134 C 1025.09,592.19 978.75,614.75"
		path_str.append("M");
		path_str.append(" ");
		path_str.append(to_string(pointListOfLines[i].pointlist.front().x));//不用TransfloatToStr()，直接用to_string()也行
		path_str.append(",");
		path_str.append(to_string(ginstance->file_height - pointListOfLines[i].pointlist.front().y));
		path_str.append("L");
		for (unsigned int j = 0; j < pointListOfLines[i].pointlist.size(); ++j)
		{
			path_str.append(" ");
			path_str.append(to_string(pointListOfLines[i].pointlist[j].x));
			path_str.append(",");
			path_str.append(to_string(ginstance->file_height - pointListOfLines[i].pointlist[j].y));
		}
		apath->SetAttribute("d", path_str.c_str());
		rootele->LinkEndChild(apath);
	}
	
	//int cbsize = total_pw_cbcurve.size();  	//	vector<vector<cubicBezier>> total_pw_cbcurve;
	//for (int i = 0; i < cbsize; ++i)
	//{
	//	TiXmlElement* apath = new TiXmlElement("path");
	//	apath->SetAttribute("style", "fill:none;stroke-width:1;stroke:#606060");//线条是灰色
	//	//apath->SetAttribute("d", "220,100 300,210 170,250");
	//	//polyline 中点的转换
	//	vector<cubicBezier> path_cb= total_pw_cbcurve[i];
	//	cubicBezier headcb = path_cb[0];
	//	string path_str;       //d="M 1062.5,587.134 C 1025.09,592.19 978.75,614.75"
	//	path_str.append("M");
	//	path_str.append(" ");
	//	path_str.append(to_string(headcb.p0.x));//不用TransfloatToStr()，直接用to_string()也行
	//	path_str.append(",");
	//	path_str.append(to_string(ginstance->file_height - headcb.p0.y));
	//	path_str.append("C");
	//	for (size_t j = 0; j < path_cb.size(); ++j)
	//	{
	//		cubicBezier tmpcb1 = path_cb[j];
	//		path_str.append(" ");
	//		path_str.append(to_string(tmpcb1.p1.x));
	//		path_str.append(",");
	//		path_str.append(to_string(ginstance->file_height - tmpcb1.p1.y));
	//		path_str.append(" ");
	//		path_str.append(to_string(tmpcb1.p2.x));
	//		path_str.append(",");
	//		path_str.append(to_string(ginstance->file_height - tmpcb1.p2.y));
	//		path_str.append(" ");
	//		path_str.append(to_string(tmpcb1.p3.x));
	//		path_str.append(",");
	//		path_str.append(to_string(ginstance->file_height - tmpcb1.p3.y));
	//		path_str.append("C");
	//	}
	//	path_str.pop_back();
	//	apath->SetAttribute("d", path_str.c_str());
	//	rootele->LinkEndChild(apath);
	//}

	string savefile;
	savefile.append(outpath);
	savefile.append(fname);
	savefile.append("_output.svg");
	writeDoc->SaveFile(savefile.c_str());
	delete writeDoc;
}
void MyGraph::writeSVGFile_step(string fname, string outpath, double fwidth, double fheight, vector<double> fview, string view_str)
{
	/*vector<pair<int,vector<mypoint2f>>> pointListOfPatches;
	vector<PolyLine> pointListOfLines;*/
	TiXmlDocument *writeDoc = new TiXmlDocument(); //xml文档指针  
	//XML 文件声明。 standalone=no 属性:该属性规定此 SVG 文件是否是“独立的”，或含有对外部文件的引用。 
	//standalone = "no" 意味着 SVG 文档会引用一个外部文件-->在这里，是 DTD 文件。
	TiXmlDeclaration *decl = new TiXmlDeclaration("1.0", "UTF-8", "no");
	writeDoc->LinkEndChild(decl); //写入文档  
	//有错
	//TiXmlElement *refer_file = new TiXmlElement("!DOCTYPE");
	//refer_file->SetAttribute("svg PUBLIC", "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd");
	//writeDoc->LinkEndChild(refer_file);
	//根元素svg
	TiXmlElement *rootele = new TiXmlElement("svg");
	rootele->SetAttribute("version", "1.1");
	rootele->SetAttribute("xmlns", "http://www.w3.org/2000/svg");
	rootele->SetAttribute("xmlns:xlink", "http://www.w3.org/1999/xlink");
	rootele->SetAttribute("width", fwidth);
	rootele->SetAttribute("height", fheight);
	rootele->SetAttribute("viewBox", view_str.c_str());
	rootele->SetAttribute("xml:space", "preserve");
	writeDoc->LinkEndChild(rootele);
	//一个path是一个多边形面
	//例子：
	//注意 fill-rule属性：nonzero, evenodd。 style = "fill:white;stroke:red;stroke-width:2"
	/* <path d="M 100 100 L 300 100 L 200 300 z" fill="red" stroke="blue" stroke-width="3" />*/
	/*<polygon points="220,100 300,210 170,250" style="fill:#cccccc;stroke:#000000;stroke-width:1"/>*/
	vector<char*> color_vec;
	char *c1 = "ff0000"; color_vec.push_back(c1);
	char *c2 = "00ff00"; color_vec.push_back(c2);
	char *c3 = "0000ff"; color_vec.push_back(c3);
	char *c4 = "ffff00"; color_vec.push_back(c4);
	char *c5 = "00ffff"; color_vec.push_back(c5);
	char *c6 = "ff00ff"; color_vec.push_back(c6);
	char *c7 = "b400a0"; color_vec.push_back(c7);
	char *c8 = "ff7878"; color_vec.push_back(c8);
	char *c9 = "c8ff00"; color_vec.push_back(c9);
	char *c10 = "ff96ff"; color_vec.push_back(c10);
	char *c11 = "3cc8ff"; color_vec.push_back(c11);
	char *c12 = "ffa000"; color_vec.push_back(c12);
	char *c13 = "9655e6"; color_vec.push_back(c13);
	char *c14 = "509600"; color_vec.push_back(c14);
	char *c15 = "fff0c8"; color_vec.push_back(c15);
	char *c16 = "9acd32"; color_vec.push_back(c16);
	char *c17 = "ff6eb4"; color_vec.push_back(c17);
	char *c18 = "4169ff"; color_vec.push_back(c18);
	int patch_count = 0;
	for (int i = 0; i < pointListOfPatches.size(); ++i)     //patch_simplified_plist   pointListOfPatches
	{
		if (pointListOfPatches[i].first >= 0)//不能=-1
		{
			TiXmlElement* apath = new TiXmlElement("polygon");
			//apath->SetAttribute("style", "fill:#aaaaaa;stroke:#000000;stroke-width:2;stroke:red");
			//apath->SetAttribute("points", "110,100 200,210 70,250");
			vector<mypoint2f> p_vec = pointListOfPatches[i].second;
			int psize = p_vec.size();
			string path_str;       //points="220,100 300,210 170,250"
			for (unsigned int j = 0; j < psize; ++j)
			{
				path_str.append(to_string(p_vec[j].x));
				path_str.append(",");
				path_str.append(to_string(ginstance->file_height - p_vec[j].y)); // 转置！！！！！
				path_str.append(" ");
			}
			path_str.pop_back();
			apath->SetAttribute("points", path_str.c_str());
			string style_str = "stroke-opacity:1;stroke-width:1;stroke:#606060;fill:#a0a0a0";//灰黑边，内灰
			//string style_str = "stroke-opacity:1;stroke-width:1;stroke:#";//stroke的颜色和patch保持一致
			//int cnum = patch_count % 18;
			//style_str.append(color_vec[cnum]);
			//style_str.append(";fill:#");
			//style_str.append(color_vec[cnum]);
			apath->SetAttribute("style", style_str.c_str());
			rootele->LinkEndChild(apath);
			patch_count++;
		}
	}
	for (unsigned int i = 0; i < pointListOfLines.size(); ++i)
	{
		TiXmlElement* apath = new TiXmlElement("path");
		apath->SetAttribute("style", "fill:none;stroke-width:1;stroke:#606060");//线条是灰色
		//apath->SetAttribute("d", "220,100 300,210 170,250");
		//polyline 中点的转换
		string path_str;       //d="M 1062.5,587.134 C 1025.09,592.19 978.75,614.75"
		path_str.append("M");
		path_str.append(" ");
		path_str.append(to_string(pointListOfLines[i].pointlist.front().x));//不用TransfloatToStr()，直接用to_string()也行
		path_str.append(",");
		path_str.append(to_string(ginstance->file_height - pointListOfLines[i].pointlist.front().y));
		path_str.append("L");
		for (unsigned int j = 0; j < pointListOfLines[i].pointlist.size(); ++j)
		{
			path_str.append(" ");
			path_str.append(to_string(pointListOfLines[i].pointlist[j].x));
			path_str.append(",");
			path_str.append(to_string(ginstance->file_height - pointListOfLines[i].pointlist[j].y));
		}
		apath->SetAttribute("d", path_str.c_str());
		rootele->LinkEndChild(apath);
	}

	for (unsigned int i = 0; i < pandl_combine.size(); ++i)
	{
		TiXmlElement* apath = new TiXmlElement("polygon");
		//apath->SetAttribute("style", "fill:#aaaaaa;stroke:#000000;stroke-width:2;stroke:red");
		//apath->SetAttribute("points", "110,100 200,210 70,250");
		vector<mypoint2f> p_vec = pandl_combine[i];
		int psize = p_vec.size();
		string path_str;       //points="220,100 300,210 170,250"
		for (unsigned int j = 0; j < psize; ++j)
		{
			path_str.append(to_string(p_vec[j].x));
			path_str.append(",");
			path_str.append(to_string(ginstance->file_height - p_vec[j].y)); // 转置！！！！！
			path_str.append(" ");
		}
		path_str.pop_back();
		apath->SetAttribute("points", path_str.c_str());
		string style_str = "stroke-opacity:.6;stroke-width:1;stroke:#ff0000;fill-opacity:.6;fill:#ff0000";//红边，内红  透明度设置
		apath->SetAttribute("style", style_str.c_str());
		rootele->LinkEndChild(apath);
	}
	for (unsigned int i = 0; i < pandl_combine_line.size(); ++i)
	{
		TiXmlElement* apath = new TiXmlElement("path");
		apath->SetAttribute("style", "fill:none;stroke-width:1;stroke:#00ff00");//线条是绿色
		//apath->SetAttribute("d", "220,100 300,210 170,250");
		//polyline 中点的转换
		string path_str;       //d="M 1062.5,587.134 C 1025.09,592.19 978.75,614.75"
		path_str.append("M");
		path_str.append(" ");
		path_str.append(to_string(pandl_combine_line[i].front().x));//不用TransfloatToStr()，直接用to_string()也行
		path_str.append(",");
		path_str.append(to_string(ginstance->file_height - pandl_combine_line[i].front().y));
		path_str.append("L");
		for (unsigned int j = 0; j < pandl_combine_line[i].size(); ++j)
		{
			path_str.append(" ");
			path_str.append(to_string(pandl_combine_line[i][j].x));
			path_str.append(",");
			path_str.append(to_string(ginstance->file_height - pandl_combine_line[i][j].y));
		}
		apath->SetAttribute("d", path_str.c_str());
		rootele->LinkEndChild(apath);
	}
	pointListOfPatches.swap(vector<pair<int, vector<mypoint2f>>>());
	pointListOfLines.swap(vector<PolyLine>());
	pandl_combine.swap(vector<vector<mypoint2f>>());
	pandl_combine_line.swap(vector<vector<mypoint2f>>());

	string savefile;
	savefile.append(outpath);
	savefile.append(fname);
	savefile.append("_output.svg");
	writeDoc->SaveFile(savefile.c_str());
	delete writeDoc;
}
void MyGraph::writesmoothFile(string fname, string outpath, double fwidth, double fheight, vector<double> fview, string view_str)
{
	/*vector<pair<int,vector<mypoint2f>>> pointListOfPatches;
	vector<PolyLine> pointListOfLines;*/
	TiXmlDocument *writeDoc = new TiXmlDocument(); //xml文档指针  
												   //XML 文件声明。 standalone=no 属性:该属性规定此 SVG 文件是否是“独立的”，或含有对外部文件的引用。 
												   //standalone = "no" 意味着 SVG 文档会引用一个外部文件-->在这里，是 DTD 文件。
	TiXmlDeclaration *decl = new TiXmlDeclaration("1.0", "UTF-8", "no");
	writeDoc->LinkEndChild(decl); //写入文档  
								  //有错
								  //TiXmlElement *refer_file = new TiXmlElement("!DOCTYPE");
								  //refer_file->SetAttribute("svg PUBLIC", "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd");
								  //writeDoc->LinkEndChild(refer_file);
								  //根元素svg
	TiXmlElement *rootele = new TiXmlElement("svg");
	rootele->SetAttribute("version", "1.1");
	rootele->SetAttribute("xmlns", "http://www.w3.org/2000/svg");
	rootele->SetAttribute("xmlns:xlink", "http://www.w3.org/1999/xlink");
	rootele->SetAttribute("width", fwidth);
	rootele->SetAttribute("height", fheight);
	rootele->SetAttribute("viewBox", view_str.c_str());
	rootele->SetAttribute("xml:space", "preserve");
	writeDoc->LinkEndChild(rootele);
	//一个path是一个多边形面
	//例子：
	//注意 fill-rule属性：nonzero, evenodd。 style = "fill:white;stroke:red;stroke-width:2"
	/* <path d="M 100 100 L 300 100 L 200 300 z" fill="red" stroke="blue" stroke-width="3" />*/
	/*<polygon points="220,100 300,210 170,250" style="fill:#cccccc;stroke:#000000;stroke-width:1"/>*/
	vector<char*> color_vec;
	char *c1 = "ff0000"; color_vec.push_back(c1);
	char *c2 = "00ff00"; color_vec.push_back(c2);
	char *c3 = "0000ff"; color_vec.push_back(c3);
	char *c4 = "ffff00"; color_vec.push_back(c4);
	char *c5 = "00ffff"; color_vec.push_back(c5);
	char *c6 = "ff00ff"; color_vec.push_back(c6);
	char *c7 = "b400a0"; color_vec.push_back(c7);
	char *c8 = "ff7878"; color_vec.push_back(c8);
	char *c9 = "c8ff00"; color_vec.push_back(c9);
	char *c10 = "ff96ff"; color_vec.push_back(c10);
	char *c11 = "3cc8ff"; color_vec.push_back(c11);
	char *c12 = "ffa000"; color_vec.push_back(c12);
	char *c13 = "9655e6"; color_vec.push_back(c13);
	char *c14 = "509600"; color_vec.push_back(c14);
	char *c15 = "fff0c8"; color_vec.push_back(c15);
	char *c16 = "9acd32"; color_vec.push_back(c16);
	char *c17 = "ff6eb4"; color_vec.push_back(c17);
	char *c18 = "4169ff"; color_vec.push_back(c18);
	int patch_count = 0;

	int cbsize = total_pw_cbcurve.size();  	//	vector<vector<cubicBezier>> total_pw_cbcurve;
	for (int i = 0; i < cbsize; ++i)
	{
		TiXmlElement* apath = new TiXmlElement("path");
		apath->SetAttribute("style", "fill:none;stroke-width:1;stroke:#606060");//线条是灰色
		//apath->SetAttribute("d", "220,100 300,210 170,250");
		//polyline 中点的转换
		vector<cubicBezier> path_cb= total_pw_cbcurve[i];
		cubicBezier headcb = path_cb[0];
		string path_str;       //d="M 1062.5,587.134 C 1025.09,592.19 978.75,614.75"
		path_str.append("M");
		path_str.append(" ");
		path_str.append(to_string(headcb.p0.x));//不用TransfloatToStr()，直接用to_string()也行
		path_str.append(",");
		path_str.append(to_string(ginstance->file_height - headcb.p0.y));
		path_str.append("C");
		for (size_t j = 0; j < path_cb.size(); ++j)
		{
			cubicBezier tmpcb1 = path_cb[j];
			path_str.append(" ");
			path_str.append(to_string(tmpcb1.p1.x));
			path_str.append(",");
			path_str.append(to_string(ginstance->file_height - tmpcb1.p1.y));
			path_str.append(" ");
			path_str.append(to_string(tmpcb1.p2.x));
			path_str.append(",");
			path_str.append(to_string(ginstance->file_height - tmpcb1.p2.y));
			path_str.append(" ");
			path_str.append(to_string(tmpcb1.p3.x));
			path_str.append(",");
			path_str.append(to_string(ginstance->file_height - tmpcb1.p3.y));
			path_str.append("C");
		}
		path_str.pop_back();
		apath->SetAttribute("d", path_str.c_str());
		rootele->LinkEndChild(apath);
	}

	string savefile;
	savefile.append(outpath);
	savefile.append(fname);
	savefile.append("_Smooth.svg");
	writeDoc->SaveFile(savefile.c_str());
	delete writeDoc;
}


//about opengl
void MyGraph::testShowGraph()
{
	setGraphInstance();
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	glutInitWindowPosition(50, 100);
	glutInitWindowSize(winWidth, winHeight);
	glutCreateWindow("graph");
	init();
	glutDisplayFunc(graphDisplay);
	glutReshapeFunc(winreshape);
	glutSpecialFunc(myspecialkey);
	glutMainLoop();
}
void MyGraph::setGraphInstance()
{
	//cout << "MyGraph::instance = this;" << endl;
	MyGraph::ginstance = this;
}
void MyGraph::graphDisplay()
{
	// GLfloat view_y = 1052;
	glClearColor(1.0, 1.0, 1.0, 0.0);
	glClear(GL_COLOR_BUFFER_BIT);
	glEnable(GL_BLEND);  //开启混合模式  
	//源因子（将要画上去）目标因子（原来的颜色）
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); //glBlendFunc(GL_ONE, GL_ZERO);  //
	glDisable(GL_DEPTH_TEST);  //关闭深度测试 

	glLineWidth(2);//设置线段宽度
	glBegin(GL_LINES);//x axis ,y axis
	//glColor4f(1.0, 0.0, 0.0, 0.5);
	glVertex2f(-10.0f, 0.0f);
	glVertex2f(100.0f, 0.0f);
	glVertex2f(0.0f, -10.0f);
	glVertex2f(0.0f, 100.0f);
	glEnd();

	//showEdgeGraph();
	showLines();//paintPolyline(ginstance->pointListOfLines);
	showpatches();
	//showBoundingBox(svginstance->polygon_vec);
	//showPatchTopology(ginstance->patch_subgraphs);
	//showClassifiedPatch("cluster");
	//showPatchTopology(ginstance->patch_subgraphs_group);
	//showClassifiedPatch("group");

	paintCtlPoint(ginstance->circle_vec);//circle black
	//paintCubicBezierCurve(ginstance->cubicBezier_vec);//red
	//paintQuaBezierCurve(ginstance->quadraticBezier_vec);//red
	//glColor4f(0.0, 0.0, 0.0, 0.5);
	//paintPolyline(ginstance->polyline_vec);//polyline green

	glFlush();
}
void MyGraph::showEdgeGraph()
{
	//vector<PolyLine> poly_vec;	
	//for (unsigned int i = 0; i < ginstance->myEdgeGraph.edgesOfGraph.size(); i++)
	//{
	//	PolyLine p1;
	//	vector<mypoint2f> plist;
	//	ginstance->getPListOfaCurve_two(i, ginstance->myEdgeGraph.edgesOfGraph[i].T1.p_position, ginstance->myEdgeGraph.edgesOfGraph[i].T2.p_position, &plist);
	//	p1.pointlist.insert(p1.pointlist.end(), plist.begin(), plist.end());
	//	poly_vec.push_back(p1);
	//}
	//glColor4f(0.0, 0.0, 0.0, 0.5);
	//paintPolyline(poly_vec);

	//show original pic  (sgraph.nodeList，findline()时删掉了一些节点，不全，用pointTopl)
	vector<PolyLine> poly_vec;
	int count = 0;
	int *visitStatus = new int[ginstance->sgraph.vertexNum];
	//0:unvisited,1:push in queue,2: out of queue
	for (int i = 0; i < ginstance->sgraph.vertexNum; ++i)
	{
		visitStatus[i] = 0;
	}
	for (int i = 0; i < ginstance->sgraph.vertexNum; ++i)//breadth-first traversal  BFT
	{
		if (visitStatus[i] == 0 && ginstance->countEdgeNum(i,0)>0)
		{
			int tmp_ver = 0;
			queue<int> verQueue;
			verQueue.push(i);
			visitStatus[i] = 1;
			while (!verQueue.empty())
			{
				tmp_ver = verQueue.front();
				verQueue.pop();
				visitStatus[tmp_ver] = 2;
				ArcNode *tmp_arc = ginstance->sgraph.nodeList[tmp_ver].firstEdge;
				while (tmp_arc)
				{
					if (visitStatus[tmp_arc->adjVertex] != 2)
					{
						if (visitStatus[tmp_arc->adjVertex] == 0)
						{
							verQueue.push(tmp_arc->adjVertex);
							visitStatus[tmp_arc->adjVertex] = 1;
						}
						////show inner line
						//ArcNode *tt = ginstance->getEdgeNode(tmp_arc->adjVertex, tmp_ver);
						//if (tmp_arc->patchesIndex.empty() && tt->patchesIndex.empty())
						//	//if (tmp_arc->patchesIndex.empty() && tt->patchesIndex.empty() && tmp_arc->innerLine == false && tt->innerLine == false)
						//{
						//	PolyLine p1;
						//	ginstance->getPListOfaCurve(ginstance->sgraph.nodeList[tmp_ver].position, ginstance->sgraph.nodeList[tmp_arc->adjVertex].position, tmp_ver, 0, &p1.pointlist);
						//	poly_vec.push_back(p1);
						//	count++;
						//}
						PolyLine p1;
						ginstance->getPListOfaCurve(ginstance->sgraph.nodeList[tmp_ver].position, ginstance->sgraph.nodeList[tmp_arc->adjVertex].position, tmp_ver, tmp_arc->adjVertex, &p1.pointlist,0);
						poly_vec.push_back(p1);
						count++;
					}				
					tmp_arc = tmp_arc->next;
				}
			}
		}
	}
	delete[]visitStatus;
	visitStatus = NULL;
	glColor4f(0.0, 0.0, 0.0, 0.5);
	paintPolyline(poly_vec);
}
void MyGraph::showLines()
{
	glColor4f(0.0, 0.0, 0.0, 0.5);
	paintPolyline(ginstance->pointListOfLines);
}
void MyGraph::showpatches()
{	
	//GLfloat view_y = 1052.0;
	GLfloat colors[][4] = {
		{ 1.0, 0.0, 0.0, 0.7 },
		{ 0.0, 0.0, 1.0, 0.7 },
		{ 0.0, 1.0, 0.0, 0.7},
		{ 1.0, 1.0, 0.0, 0.7 },//黄
		{ 1.0, 0.0, 1.0, 0.7 },//洋红
		{ 0.0, 1.0, 1.0, 0.7 },//青

		{ 1.0, 0.5, 0.0, 0.7 },//橙
		{ 1.0, 0.0, 0.5, 0.7 },//桃红
		{ 0.5, 1.0, 0.0, 0.7 },//草绿
		{ 0.0, 1.0, 0.5, 0.7 },//湖绿
		{ 0.5, 0.0, 1.0, 0.7 },//紫
		{ 0.0, 0.5, 1.0, 0.7 },//孔蓝

		{ 1.0, 0.5, 0.5, 0.7 },//肤粉
		{ 0.0, 0.7, 0.7, 0.7},//湖蓝
		{ 0.6, 1.0, 0.6, 0.7 },//浅绿
		{ 0.6, 1.0, 0.6, 0.7 },//紫红
		{ 0.5, 0.5, 1.0, 0.7 },//淡紫
		{ 1.0, 0.5, 0.0, 0.7 },//橘红

		/*{ 1.0, 0.1, 0.1, 0.5 },
		{ 0.1, 0.1, 1.0, 0.5 },
		{ 0.1, 1.0, 0.1, 0.5 },
		{ 1.0, 1.0, 0.1, 0.5 },
		{ 1.0, 0.1, 1.0, 0.5 },
		{ 0.1, 1.0, 1.0, 0.5 },
		{ 0.0, 0.0, 0.0, 0.5 },
		*/
	};
	//draw patches boundary-->GL_LINE_STRIP
	for (unsigned int i = 0; i < ginstance->pointListOfPatches.size(); ++i)  //patch_simplified_plist  pointListOfPatches
	{
		vector<mypoint2f> tmpvec = ginstance->pointListOfPatches[i].second;
		glColor4fv(colors[(i%12)]);
		glBegin(GL_LINE_STRIP);//GL_POLYGON GL_LINE_STRIP
		for (unsigned int j = 0; j < tmpvec.size(); j++)
		{
			glVertex2f((GLfloat)tmpvec[j].x, (GLfloat)tmpvec[j].y);
		}
		glEnd();
	}	

	//paint patches with colors
	GLUtesselator *tobj = gluNewTess();
	if (!tobj)
	{
		return;
	}
	gluTessCallback(tobj, GLU_TESS_VERTEX, (void (CALLBACK *)())vertexCallback);
	gluTessCallback(tobj, GLU_TESS_BEGIN, (void (CALLBACK *)())beginCallback);
	gluTessCallback(tobj, GLU_TESS_END, (void (CALLBACK *)())endCallback);
	for (unsigned int i = 0; i < ginstance->pointListOfPatches.size(); ++i)  ////////
	{
		vector<mypoint2f> tmpvec = ginstance->pointListOfPatches[i].second;//////////
		GLdouble **pointarre = new GLdouble*[tmpvec.size()];
		for (int j = 0; j < tmpvec.size(); ++j)
		{
			pointarre[j] = new GLdouble[3];
		}
		for (int j = 0; j < tmpvec.size(); ++j)
		{
			//pointarre[i] = new GLdouble[3];
			//pointarre[i] = { (GLdouble)tmpvec[j].x, (GLdouble)tmpvec[j].y, 0.0 };
			pointarre[j][0] = (GLdouble)tmpvec[j].x;
			pointarre[j][1] = (GLdouble)tmpvec[j].y;
			pointarre[j][2] = (GLdouble)0.0;
		}
		gluTessBeginPolygon(tobj, NULL);
		gluTessBeginContour(tobj);
		//if (ginstance->pointListOfPatches[i].first==0)
		//	glColor4fv(colors[(i % 12)]);
		//else if (ginstance->pointListOfPatches[i].first == 1)
		//	glColor4fv(colors[13]);
		//else if (ginstance->pointListOfPatches[i].first == 2)
		//	glColor4fv(colors[14]);
		//else if (ginstance->pointListOfPatches[i].first == 3)
		//	glColor4fv(colors[15]);
		//else if (ginstance->pointListOfPatches[i].first == 4)
		//	glColor4fv(colors[16]);
		//else if (ginstance->pointListOfPatches[i].first == 5)
		//	glColor4fv(colors[17]);

	    glColor4fv(colors[(i % 18)]);   //18 color
		for (int k = 0; k < tmpvec.size(); k++)
		{
			gluTessVertex(tobj, pointarre[k], pointarre[k]);
		}
		gluTessEndContour(tobj);
		gluTessEndPolygon(tobj);
		delete[]pointarre;
		pointarre = NULL;
	}
	gluDeleteTess(tobj);

}
void MyGraph::showPatchTopology(vector<pair<int, string>> patch_subgraph)
{
	int *visited = new int[ginstance->patchTopology.patch_number];
	for (int i = 0; i < ginstance->patchTopology.patch_number; ++i)
	{
		visited[i] = 0;	//0:unvisited,1:push in queue,2: out of queue
	}
	vector<mypoint2f> pchsite;
	vector<vector<mypoint2f>> lsegment;
	for (unsigned int i = 0; i < patch_subgraph.size(); ++i)
	{
		if (patch_subgraph[i].second != "node")
		{
			int startp = patch_subgraph[i].first;
			queue<int> subg_q;
			subg_q.push(startp);
			visited[startp] = 1;
			while (!subg_q.empty())
			{
				int tmp_ver = 0;
				tmp_ver = subg_q.front();
				subg_q.pop();
				visited[tmp_ver] = 2;
				PatchArc *tmp_arc = ginstance->patchTopology.allpatches[tmp_ver].firstAdjPatch;
				while (tmp_arc != NULL)
				{
					if (visited[tmp_arc->adj_patchindex] != 2)
					{
						if (visited[tmp_arc->adj_patchindex] == 0)
						{
							subg_q.push(tmp_arc->adj_patchindex);
							visited[tmp_arc->adj_patchindex] = 1;
						}
						vector<mypoint2f> pvec;//用质心画出patch的拓扑图，patchtopology的个数和plrgraph.patch是一样的
						pvec.push_back(ginstance->plrGraph.patches[tmp_ver].centre_position);
						pvec.push_back(ginstance->plrGraph.patches[tmp_arc->adj_patchindex].centre_position);
						pchsite.push_back(pvec[0]);
						pchsite.push_back(pvec[1]);
						lsegment.push_back(pvec);
					}
					tmp_arc = tmp_arc->nextPatch;
				}
			}
		}
		else
			pchsite.push_back(ginstance->plrGraph.patches[patch_subgraph[i].first].centre_position);
	}
	delete[]visited;
	visited = NULL;

	//GLfloat view_y = 1052.0;
	glPointSize(5.0);
	glColor4f(0.0, 0.0, 0.0, 0.5);
	glBegin(GL_POINTS);
	for (unsigned int i = 0; i < pchsite.size(); ++i)
	{
		glVertex2f(pchsite[i].x, pchsite[i].y);//第二种计算中心方法是用 ：glVertex2f(pchsite[i].x, pchsite[i].y);//
	}
	glEnd();

	glColor4f(0.0, 0.0, 0.0, 0.5);
	glBegin(GL_LINES);
	for (unsigned int i = 0; i < lsegment.size(); ++i)
	{
		glVertex2f(lsegment[i].at(0).x, lsegment[i].at(0).y);
		glVertex2f(lsegment[i].at(1).x, lsegment[i].at(1).y);
		/*glVertex2f(lsegment[i].at(0).x,lsegment[i].at(0).y);
		glVertex2f(lsegment[i].at(1).x, lsegment[i].at(1).y);*/
	}
	glEnd();
}
void MyGraph::showClassifiedPatch(string str)
{
	//GLfloat view_y = 1052.0;
	GLfloat colors[][4] = {
		{ 1.0, 0.0, 0.0, 0.6 },
		{ 0.0, 0.0, 1.0, 0.6 },
		{ 0.0, 1.0, 0.0, 0.6 },
		{ 1.0, 1.0, 0.0, 0.6 },
		{ 1.0, 0.0, 1.0, 0.6 },
		{ 0.0, 1.0, 1.0, 0.6 },

		{ 1.0, 0.5, 0.0, 0.6 },
		{ 1.0, 0.0, 0.5, 0.6 },
		{ 0.5, 1.0, 0.0, 0.6 },
		{ 0.0, 1.0, 0.5, 0.6 },
		{ 0.5, 0.0, 1.0, 0.6 },
		{ 0.0, 0.5, 1.0, 0.6 },

		//{ 0.0, 0.2, 0.8, 0.9 },//深蓝
		//{ 0.2, 0.16, 0.68, 0.9 },
		//{ 0.40, 0.12, 0.56, 0.9 },
		//{ 0.6, 0.08, 0.44, 0.9 },
		//{ 0.8, 0.04, 0.32, 0.9 },
		//{ 1.0, 0.0, 0.2, 0.9 },//红

	};
	//draw patches boundary-->GL_LINE_STRIP
	for (unsigned int i = 0; i < ginstance->pointListOfPatches.size(); ++i)
	{
		vector<mypoint2f> tmpvec = ginstance->pointListOfPatches[i].second;
		glColor4f(0.0,0.0,0.0,0.5);//glColor4fv(colors[(i % 12)]);
		glBegin(GL_LINE_STRIP);//GL_POLYGON GL_LINE_STRIP
		for (unsigned int j = 0; j < tmpvec.size(); j++)
		{
			glVertex2f((GLfloat)tmpvec[j].x, tmpvec[j].y);
		}
		glEnd();
	}

	//paint patches with colors
	GLUtesselator *tobj = gluNewTess();
	if (!tobj)
	{
		return;
	}
	gluTessCallback(tobj, GLU_TESS_VERTEX, (void (CALLBACK *)())vertexCallback);
	gluTessCallback(tobj, GLU_TESS_BEGIN, (void (CALLBACK *)())beginCallback);
	gluTessCallback(tobj, GLU_TESS_END, (void (CALLBACK *)())endCallback);
	for (unsigned int i = 0; i < ginstance->pointListOfPatches.size(); ++i)
	{
		vector<mypoint2f> tmpvec = ginstance->pointListOfPatches[i].second;
		GLdouble **pointarre = new GLdouble*[tmpvec.size()];
		for (int j = 0; j < tmpvec.size(); ++j)
		{
			pointarre[j] = new GLdouble[3];
		}
		for (int j = 0; j < tmpvec.size(); ++j)
		{
			//pointarre[i] = new GLdouble[3];
			//pointarre[i] = { (GLdouble)tmpvec[j].x, (GLdouble)tmpvec[j].y, 0.0 };
			pointarre[j][0] = (GLdouble)tmpvec[j].x;
			pointarre[j][1] = (GLdouble)tmpvec[j].y;
			pointarre[j][2] = (GLdouble)0.0;
		}
		gluTessBeginPolygon(tobj, NULL);
		gluTessBeginContour(tobj);
		if (str == "group") //group  cluster
		{
			int ttp = ginstance->pointListOfPatches[i].first % 12;
			if (ttp==0)
				glColor4f(0.0, 0.0, 0.0, 0.5);
			else
				glColor4fv(colors[ttp-1]);
		}
		else if (str == "cluster")
		{
			if (ginstance->pointListOfPatches[i].first == 0)//灰色 
				glColor4f(0.0, 0.0, 0.0, 0.5);
			if (ginstance->pointListOfPatches[i].first == 1)//红
				glColor4f(1.0, 0.0, 0.0, 0.8);
			if (ginstance->pointListOfPatches[i].first == 2)//蓝色
				glColor4f(0.0, 0.0, 1.0, 0.8);
			if (ginstance->pointListOfPatches[i].first == 3)//绿色(噪声点)
				glColor4f(0.0, 1.0, 0.0, 0.8);
			if (ginstance->pointListOfPatches[i].first == 4)//?色
				glColor4f(1.0, 1.0, 0.0, 0.8);
		}	
		for (int k = 0; k < tmpvec.size(); k++)
		{
			gluTessVertex(tobj, pointarre[k], pointarre[k]);
		}
		gluTessEndContour(tobj);
		gluTessEndPolygon(tobj);
		delete[]pointarre;
		pointarre = NULL;
	}
	gluDeleteTess(tobj);
}
void MyGraph::showBoundingBox(vector<vector<mypoint2f>> plist)
{
	GLfloat px = 0; GLfloat py = 0;	//GLfloat view_y = 1052.0;
	glColor4f(0.0, 0.0, 0.0, 0.5);  //黑色
	for (unsigned int i = 0; i < plist.size(); i++)
	{	
		glBegin(GL_LINE_STRIP);
		for (unsigned int j = 0; j < plist[i].size(); j++)
		{
			px = (GLfloat)plist[i].at(j).x;
			py = (GLfloat) plist[i].at(j).y;
			glVertex2f(px, py);
		}
		glEnd();
	}
}


void CALLBACK MyGraph::beginCallback(GLenum type)
{
	glBegin(type);
}
void CALLBACK MyGraph::vertexCallback(GLdouble * vertex)
{
	const GLdouble *pointer = (GLdouble *)vertex;
	//glColor4f(1.0,0.1,0.1,0.5);
	glVertex3dv(pointer);
}
void CALLBACK MyGraph::endCallback()
{
	glEnd();
}

///////////////////////////////////////////////////////////////////////////////////////// edge Graph


/*边图edgeGraph做的尝试，只是找线，还没有找面和innerline*/
void MyGraph::addAssitantP_two(int edge_i, int edge_j, bool &delflag)
{
	vector<mypoint2f> plist_i, plist_j;
	getPListOfaCurve_two(edge_i, myEdgeGraph.edgesOfGraph[edge_i].T1.p_position, myEdgeGraph.edgesOfGraph[edge_i].T2.p_position, &plist_i);
	getPListOfaCurve_two(edge_j, myEdgeGraph.edgesOfGraph[edge_i].T1.p_position, myEdgeGraph.edgesOfGraph[edge_i].T2.p_position, &plist_j);
	if (plist_i.size() == plist_j.size() && plist_i.size() == 2)
	{
		////delete plist_j
		//delflag = true;
		int rmEdgeIndex = myEdgeGraph.edgesOfGraph[edge_j].index;
		EdgeStruct* eg = &myEdgeGraph.edgesOfGraph[edge_i];
		removeEdgeOfTP(eg, myEdgeGraph.edgesOfGraph[edge_i].T1.p_position, rmEdgeIndex);
		/*vector<int>::iterator iter = myEdgeGraph.edgesOfGraph[edge_i].T1.adjEdgeIndex_vec.begin();
		while (iter != myEdgeGraph.edgesOfGraph[edge_i].T1.adjEdgeIndex_vec.end())
		{
			if (*iter == myEdgeGraph.edgesOfGraph[edge_j].index)
			{
				myEdgeGraph.edgesOfGraph[edge_i].T1.adjEdgeIndex_vec.erase(iter);
				break;
			}
			iter++;
		}*/
		removeEdgeOfTP(eg, myEdgeGraph.edgesOfGraph[edge_i].T2.p_position, rmEdgeIndex);
		/*iter = myEdgeGraph.edgesOfGraph[edge_i].T2.adjEdgeIndex_vec.begin();
		while (iter != myEdgeGraph.edgesOfGraph[edge_i].T2.adjEdgeIndex_vec.end())
		{
			if (*iter == myEdgeGraph.edgesOfGraph[edge_j].index)
			{
				myEdgeGraph.edgesOfGraph[edge_i].T2.adjEdgeIndex_vec.erase(iter);
				break;
			}
			iter++;
		}*/

		/*for (unsigned int i = 0; i < myEdgeGraph.edgesOfGraph[edge_i].T1.adjEdgeIndex_vec.size(); ++i)
		{
			int e_index = myEdgeGraph.edgesOfGraph[edge_i].T1.adjEdgeIndex_vec[i];
			EdgeStruct *e_p = &myEdgeGraph.edgesOfGraph[e_index];
			removeEdgeOfTP(e_p, myEdgeGraph.edgesOfGraph[edge_i].T1.p_position, rmEdgeIndex);
		}
		for (unsigned int i = 0; i < myEdgeGraph.edgesOfGraph[edge_i].T2.adjEdgeIndex_vec.size(); ++i)
		{
			int e_index = myEdgeGraph.edgesOfGraph[edge_i].T2.adjEdgeIndex_vec[i];
			EdgeStruct *e_p = &myEdgeGraph.edgesOfGraph[e_index];
			removeEdgeOfTP(e_p, myEdgeGraph.edgesOfGraph[edge_i].T2.p_position, rmEdgeIndex);
		}*/
	}
	else if (plist_j.size()>2)
	{
		//split plist_j
		if (isPointRoughEqual(myEdgeGraph.edgesOfGraph[edge_j].T2.p_position, myEdgeGraph.edgesOfGraph[edge_i].T1.p_position))
		{
			plist_j.clear();
			plist_j.swap(vector<mypoint2f>());
			getPListOfaCurve_two(edge_j, myEdgeGraph.edgesOfGraph[edge_i].T2.p_position, myEdgeGraph.edgesOfGraph[edge_i].T1.p_position, &plist_j);
		}
		int mid_index = plist_j.size() / 2;

		PolyLine apl1;
		apl1.pointlist.insert(apl1.pointlist.end(), plist_j.begin(), plist_j.begin()+mid_index);
		polyline_vec.push_back(apl1);
		myEdgeGraph.edgesOfGraph[edge_j].curveType = 'l';
		myEdgeGraph.edgesOfGraph[edge_j].curveIndex = polyline_vec.size() - 1;
		myEdgeGraph.edgesOfGraph[edge_j].T2.p_position = apl1.pointlist.back();
		mypoint2f midp = apl1.pointlist.back();///////////////////////////添加点-》circle_vec
		mycircle acir;
		acir.cx = midp.x;
		acir.cy = midp.y;
		acir.radius = 1;
		circle_vec.push_back(acir);
		PolyLine apl2;
		apl2.pointlist.insert(apl2.pointlist.end(), plist_j.begin() + mid_index-1, plist_j.end());
		polyline_vec.push_back(apl2);//polylineVec.push_back(seg_curve1);
		EdgeStruct anedge;
		anedge.index = myEdgeGraph.edgesOfGraph.back().index + 1;
		anedge.curveIndex = polyline_vec.size() - 1;
		anedge.curveType = 'l';
		anedge.T1.p_position = apl2.pointlist.front();
		anedge.T2.p_position = apl2.pointlist.back();
		myEdgeGraph.edgesOfGraph.push_back(anedge);
		myEdgeGraph.edgeNumber++;
		if (isPointRoughEqual(myEdgeGraph.edgesOfGraph[edge_j].T1.p_position, myEdgeGraph.edgesOfGraph[edge_i].T1.p_position))
		{
			if (!myEdgeGraph.edgesOfGraph[edge_i].T2.adjEdgeIndex_vec.empty())
			{
				EdgeStruct *e_p = &myEdgeGraph.edgesOfGraph[edge_i];
				removeEdgeOfTP(e_p, myEdgeGraph.edgesOfGraph[edge_i].T2.p_position, myEdgeGraph.edgesOfGraph[edge_j].index);
				/*vector<int>::iterator iter_edge=myEdgeGraph.edgesOfGraph[edge_i].T2.adjEdgeIndex_vec.begin();
				while (iter_edge != myEdgeGraph.edgesOfGraph[edge_i].T2.adjEdgeIndex_vec.end())
				{
					if (*iter_edge == myEdgeGraph.edgesOfGraph[edge_j].index)
					{
						myEdgeGraph.edgesOfGraph[edge_i].T2.adjEdgeIndex_vec.erase(iter_edge);
						break;
					}
					iter_edge++;
				}*/
			}
			myEdgeGraph.edgesOfGraph[edge_i].T2.adjEdgeIndex_vec.push_back(myEdgeGraph.edgesOfGraph.back().index);
		}
		else
		{
			if (!myEdgeGraph.edgesOfGraph[edge_i].T1.adjEdgeIndex_vec.empty())
			{
				EdgeStruct *e_p = &myEdgeGraph.edgesOfGraph[edge_i];
				removeEdgeOfTP(e_p, myEdgeGraph.edgesOfGraph[edge_i].T1.p_position, myEdgeGraph.edgesOfGraph[edge_j].index);
				/*vector<int>::iterator iter_edge = myEdgeGraph.edgesOfGraph[edge_i].T1.adjEdgeIndex_vec.begin();
				while (iter_edge != myEdgeGraph.edgesOfGraph[edge_i].T1.adjEdgeIndex_vec.end())
				{
					if (*iter_edge == myEdgeGraph.edgesOfGraph[edge_j].index)
					{
						myEdgeGraph.edgesOfGraph[edge_i].T1.adjEdgeIndex_vec.erase(iter_edge);
						break;
					}
					iter_edge++;
				}*/
			}
			myEdgeGraph.edgesOfGraph[edge_i].T1.adjEdgeIndex_vec.push_back(myEdgeGraph.edgesOfGraph.back().index);
		}
		if (!myEdgeGraph.edgesOfGraph[edge_j].T2.adjEdgeIndex_vec.empty())
		{
			EdgeStruct *e_p = &myEdgeGraph.edgesOfGraph[edge_j];
			removeEdgeOfTP(e_p, myEdgeGraph.edgesOfGraph[edge_j].T2.p_position, myEdgeGraph.edgesOfGraph[edge_i].index);
			/*vector<int>::iterator iter_edge_j = myEdgeGraph.edgesOfGraph[edge_j].T2.adjEdgeIndex_vec.begin();
			while (iter_edge_j != myEdgeGraph.edgesOfGraph[edge_j].T2.adjEdgeIndex_vec.end())
			{
				if (*iter_edge_j == myEdgeGraph.edgesOfGraph[edge_i].index)
				{
					myEdgeGraph.edgesOfGraph[edge_j].T2.adjEdgeIndex_vec.erase(iter_edge_j);
					break;
				}
				iter_edge_j++;
			}*/
			myEdgeGraph.edgesOfGraph[edge_j].T2.adjEdgeIndex_vec.push_back(myEdgeGraph.edgesOfGraph.back().index);
		}
	}
	else
	{
		//split plist_i
		int mid_index = plist_i.size() / 2;

		PolyLine apl1;
		apl1.pointlist.insert(apl1.pointlist.end(), plist_i.begin(), plist_i.begin() + mid_index);
		polyline_vec.push_back(apl1);
		myEdgeGraph.edgesOfGraph[edge_i].curveType = 'l';
		myEdgeGraph.edgesOfGraph[edge_i].curveIndex = polyline_vec.size() - 1;
		myEdgeGraph.edgesOfGraph[edge_i].T2.p_position = apl1.pointlist.back();
		mypoint2f midp = apl1.pointlist.back();///////////////////////////添加点-》circle_vec
		mycircle acir;
		acir.cx = midp.x;
		acir.cy = midp.y;
		acir.radius = 1;
		circle_vec.push_back(acir);
		PolyLine apl2;
		apl2.pointlist.insert(apl2.pointlist.end(), plist_i.begin() + mid_index-1, plist_i.end());
		polyline_vec.push_back(apl2);//polylineVec.push_back(seg_curve1);
		EdgeStruct anedge;
		anedge.index = myEdgeGraph.edgesOfGraph.back().index + 1;
		anedge.curveIndex = polyline_vec.size() - 1;
		anedge.curveType = 'l';
		anedge.T1.p_position = apl2.pointlist.front();
		anedge.T2.p_position = apl2.pointlist.back();
		myEdgeGraph.edgesOfGraph.push_back(anedge);
		myEdgeGraph.edgeNumber++;
		if (!myEdgeGraph.edgesOfGraph[edge_i].T2.adjEdgeIndex_vec.empty())
		{
			EdgeStruct *e_p = &myEdgeGraph.edgesOfGraph[edge_i];
			removeEdgeOfTP(e_p, myEdgeGraph.edgesOfGraph[edge_i].T2.p_position, myEdgeGraph.edgesOfGraph[edge_j].index);
		/*	vector<int>::iterator iter_edge = myEdgeGraph.edgesOfGraph[edge_i].T2.adjEdgeIndex_vec.begin();
			while (iter_edge != myEdgeGraph.edgesOfGraph[edge_i].T2.adjEdgeIndex_vec.end())
			{
				if (*iter_edge == myEdgeGraph.edgesOfGraph[edge_j].index)
				{
					myEdgeGraph.edgesOfGraph[edge_i].T2.adjEdgeIndex_vec.erase(iter_edge);
					break;
				}
				iter_edge++;
			}*/
		}
		myEdgeGraph.edgesOfGraph[edge_i].T2.adjEdgeIndex_vec.push_back(myEdgeGraph.edgesOfGraph.back().index);

		if (isPointRoughEqual(myEdgeGraph.edgesOfGraph[edge_j].T1.p_position, myEdgeGraph.edgesOfGraph[edge_i].T1.p_position))
		{
			if (!myEdgeGraph.edgesOfGraph[edge_j].T2.adjEdgeIndex_vec.empty())
			{
				EdgeStruct *e_p = &myEdgeGraph.edgesOfGraph[edge_j];
				removeEdgeOfTP(e_p, myEdgeGraph.edgesOfGraph[edge_j].T2.p_position, myEdgeGraph.edgesOfGraph[edge_i].index);
				/*vector<int>::iterator iter_edge_j = myEdgeGraph.edgesOfGraph[edge_j].T2.adjEdgeIndex_vec.begin();
				while (iter_edge_j != myEdgeGraph.edgesOfGraph[edge_j].T2.adjEdgeIndex_vec.end())
				{
					if (*iter_edge_j == myEdgeGraph.edgesOfGraph[edge_i].index)
					{
						myEdgeGraph.edgesOfGraph[edge_j].T2.adjEdgeIndex_vec.erase(iter_edge_j);
						break;
					}
					iter_edge_j++;
				}*/
				myEdgeGraph.edgesOfGraph[edge_j].T2.adjEdgeIndex_vec.push_back(myEdgeGraph.edgesOfGraph.back().index);
			}
		}
		else
		{
			if (!myEdgeGraph.edgesOfGraph[edge_j].T1.adjEdgeIndex_vec.empty())
			{
				EdgeStruct *e_p = &myEdgeGraph.edgesOfGraph[edge_j];
				removeEdgeOfTP(e_p, myEdgeGraph.edgesOfGraph[edge_j].T1.p_position, myEdgeGraph.edgesOfGraph[edge_i].index);
				/*vector<int>::iterator iter_edge_j = myEdgeGraph.edgesOfGraph[edge_j].T1.adjEdgeIndex_vec.begin();
				while (iter_edge_j != myEdgeGraph.edgesOfGraph[edge_j].T1.adjEdgeIndex_vec.end())
				{
					if (*iter_edge_j == myEdgeGraph.edgesOfGraph[edge_i].index)
					{
						myEdgeGraph.edgesOfGraph[edge_j].T1.adjEdgeIndex_vec.erase(iter_edge_j);
						break;
					}
					iter_edge_j++;
				}*/
				myEdgeGraph.edgesOfGraph[edge_j].T1.adjEdgeIndex_vec.push_back(myEdgeGraph.edgesOfGraph.back().index);
			}
		}
	}
}
void MyGraph::getPListOfaCurve_two(int edge_index, mypoint2f from_p, mypoint2f to_p, vector<mypoint2f> *vec)
{
	int order = -1;//positive order=0  , negative order=1
	char ctype = myEdgeGraph.edgesOfGraph[edge_index].curveType;
	int cindex = myEdgeGraph.edgesOfGraph[edge_index].curveIndex;
	if (isPointRoughEqual(from_p, myEdgeGraph.edgesOfGraph[edge_index].T1.p_position))//positive order
	{
		order = 0;
	}
	else if (isPointRoughEqual(from_p, myEdgeGraph.edgesOfGraph[edge_index].T2.p_position))
	{
		order = 1;
	}
	else
		return;
	if (ctype == 'c')
	{
		getPListfromCubic(order, cubicBezier_vec[cindex], 0,vec);
	}
	else if (ctype == 'q')
	{
		getPListfromQuadratic(order, quadraticBezier_vec[cindex],0, vec);
	}
	else if (ctype == 'l')
	{
		getPListfromPolyline(order, polyline_vec[cindex],0, vec);
	}
}
void MyGraph::findLines_two()
{
	for (unsigned int i = 0; i < myEdgeGraph.edgesOfGraph.size(); ++i)
	{
		if (myEdgeGraph.edgesOfGraph[i].T1.adjEdgeIndex_vec.size() == 0)
		{			
			vector<int> linepath;//store edge's index

			EdgeStruct *anedge = &myEdgeGraph.edgesOfGraph[i];
			TerminalVertex_E *termiP = getTPofEdge(*anedge, myEdgeGraph.edgesOfGraph[i].T1.p_position);
			TerminalVertex_E *termiP_other=NULL;
			int adjEdgeCount = termiP->adjEdgeIndex_vec.size();
			while (adjEdgeCount == 0)
			{
				linepath.push_back(anedge->index);
				termiP_other=getOtherTPofEdge(*anedge, termiP->p_position);
				if (!termiP_other->adjEdgeIndex_vec.empty())
				{
					EdgeStruct *next_anedge = &myEdgeGraph.edgesOfGraph[termiP_other->adjEdgeIndex_vec[0]];
					termiP = getTPofEdge(*next_anedge, termiP_other->p_position);

					removeWholeEdge(anedge);
					anedge = next_anedge;
					adjEdgeCount = termiP->adjEdgeIndex_vec.size();
				}
				else
					break;
			}
			if (!linepath.empty())
			{
				myplrGraph.lines.push_back(linepath);
			}
		}
		else if (myEdgeGraph.edgesOfGraph[i].T2.adjEdgeIndex_vec.size() == 0)
		{
			vector<int> linepath;
			EdgeStruct anedge = myEdgeGraph.edgesOfGraph[i];
			TerminalVertex_E *termiP = getTPofEdge(anedge, myEdgeGraph.edgesOfGraph[i].T2.p_position);
			TerminalVertex_E *termiP_other = NULL;
			int adjEdgeCount = termiP->adjEdgeIndex_vec.size();
			while (adjEdgeCount == 0)
			{
				linepath.push_back(anedge.index);
				termiP_other = getOtherTPofEdge(anedge, termiP->p_position);
				if (!termiP_other->adjEdgeIndex_vec.empty())
				{
					EdgeStruct next_anedge = myEdgeGraph.edgesOfGraph[termiP_other->adjEdgeIndex_vec[0]];
					termiP = getTPofEdge(next_anedge, termiP_other->p_position);

					removeWholeEdge(&anedge);
					anedge = next_anedge;
					adjEdgeCount = termiP->adjEdgeIndex_vec.size();
				}
				else
					break;
			}
			if (!linepath.empty())
			{
				myplrGraph.lines.push_back(linepath);
			}
		}
	}
}
void MyGraph::removeWholeEdge(EdgeStruct *e1)
{
	int rmEdge_index = e1->index;
	mypoint2f p1 = e1->T1.p_position;
	for (unsigned int i = 0; i < e1->T1.adjEdgeIndex_vec.size(); ++i)
	{
		EdgeStruct *e_t1_i=&myEdgeGraph.edgesOfGraph[e1->T1.adjEdgeIndex_vec[i]];
		removeEdgeOfTP(e_t1_i, p1, rmEdge_index);
		/*TerminalVertex_E *t1 = getTPofEdge(e_t1_i, p1);
		vector<int>::iterator iter_e = t1->adjEdgeIndex_vec.begin();
		while (iter_e != t1->adjEdgeIndex_vec.end())
		{
			if (*iter_e == rmEdge_index)
			{
				t1->adjEdgeIndex_vec.erase(iter_e);
				break;
			}
			iter_e++;
		}*/
	}
	mypoint2f p2 = e1->T2.p_position;
	for (unsigned int i = 0; i < e1->T2.adjEdgeIndex_vec.size(); ++i)
	{
		EdgeStruct *e_t2_i = &myEdgeGraph.edgesOfGraph[e1->T2.adjEdgeIndex_vec[i]];
		removeEdgeOfTP(e_t2_i, p2, rmEdge_index);
		/*TerminalVertex_E *t2 = getTPofEdge(e_t2_i, p2);
		vector<int>::iterator iter_e = t2->adjEdgeIndex_vec.begin();
		while (iter_e != t2->adjEdgeIndex_vec.end())
		{
			if (*iter_e == rmEdge_index)
			{
				t2->adjEdgeIndex_vec.erase(iter_e);
				break;
			}
			iter_e++;
		}*/
	}
	e1->T1.adjEdgeIndex_vec.clear();
	e1->T2.adjEdgeIndex_vec.clear();
}
void MyGraph::removeEdgeOfTP(EdgeStruct *eg, mypoint2f tp, int rmEdgeIndex)
{
	if (isPointRoughEqual(eg->T1.p_position, tp))
	{
		vector<int>::iterator iter = eg->T1.adjEdgeIndex_vec.begin();
		while (iter != eg->T1.adjEdgeIndex_vec.end())
		{
			if (*iter == rmEdgeIndex)
			{
				eg->T1.adjEdgeIndex_vec.erase(iter);
				break;
			}
			iter++;
		}
	}
	else if (isPointRoughEqual(eg->T2.p_position, tp))
	{
		vector<int>::iterator iter = eg->T2.adjEdgeIndex_vec.begin();
		while (iter != eg->T2.adjEdgeIndex_vec.end())
		{
			if (*iter == rmEdgeIndex)
			{
				eg->T2.adjEdgeIndex_vec.erase(iter);
				break;
			}
			iter++;
		}
	}
}
TerminalVertex_E* MyGraph::getTPofEdge(EdgeStruct e1, mypoint2f tp)
{
	if (isPointRoughEqual(myEdgeGraph.edgesOfGraph[e1.index].T1.p_position, tp))
		return &myEdgeGraph.edgesOfGraph[e1.index].T1;
	else if (isPointRoughEqual(myEdgeGraph.edgesOfGraph[e1.index].T2.p_position, tp))
		return &myEdgeGraph.edgesOfGraph[e1.index].T2;
}
TerminalVertex_E* MyGraph::getOtherTPofEdge(EdgeStruct e1, mypoint2f tp)
{
	if (isPointRoughEqual(e1.T1.p_position, tp))
		return &myEdgeGraph.edgesOfGraph[e1.index].T2;
	else if (isPointRoughEqual(e1.T2.p_position, tp))
		return &myEdgeGraph.edgesOfGraph[e1.index].T1;
}
void MyGraph::jointPolylines_two()
{
	for (unsigned int i = 0; i < myplrGraph.lines.size(); ++i)
	{
		vector<int> aline = myplrGraph.lines[i];
		if (aline.size() == 1)
		{
			PolyLine pl;
			getPListOfaCurve_two(aline[0], myEdgeGraph.edgesOfGraph[aline[0]].T1.p_position, myEdgeGraph.edgesOfGraph[aline[0]].T2.p_position, &pl.pointlist);
			pointListOfLines.push_back(pl);
		}
		else
		{
			PolyLine pl;
			mypoint2f comm_p;
			if (isPointRoughEqual(myEdgeGraph.edgesOfGraph[aline[0]].T1.p_position, myEdgeGraph.edgesOfGraph[aline[1]].T1.p_position))
			{
				comm_p = myEdgeGraph.edgesOfGraph[aline[0]].T1.p_position;
				getPListOfaCurve_two(aline[0], myEdgeGraph.edgesOfGraph[aline[0]].T2.p_position, comm_p, &pl.pointlist);
			}
			else if (isPointRoughEqual(myEdgeGraph.edgesOfGraph[aline[0]].T2.p_position, myEdgeGraph.edgesOfGraph[aline[1]].T1.p_position))
			{
				comm_p = myEdgeGraph.edgesOfGraph[aline[0]].T2.p_position;
				getPListOfaCurve_two(aline[0], myEdgeGraph.edgesOfGraph[aline[0]].T1.p_position, comm_p, &pl.pointlist);
			}
			else if (isPointRoughEqual(myEdgeGraph.edgesOfGraph[aline[0]].T1.p_position, myEdgeGraph.edgesOfGraph[aline[1]].T2.p_position))
			{
				comm_p = myEdgeGraph.edgesOfGraph[aline[0]].T1.p_position;
				getPListOfaCurve_two(aline[0], myEdgeGraph.edgesOfGraph[aline[0]].T2.p_position, comm_p, &pl.pointlist);
			}
			else if (isPointRoughEqual(myEdgeGraph.edgesOfGraph[aline[0]].T2.p_position, myEdgeGraph.edgesOfGraph[aline[1]].T2.p_position))
			{
				comm_p = myEdgeGraph.edgesOfGraph[aline[0]].T2.p_position;
				getPListOfaCurve_two(aline[0], myEdgeGraph.edgesOfGraph[aline[0]].T1.p_position, comm_p, &pl.pointlist);
			}
			for (unsigned int j = 1; j < aline.size(); ++j)
			{
				if (isPointRoughEqual(myEdgeGraph.edgesOfGraph[aline[j]].T1.p_position, comm_p))
				{
					getPListOfaCurve_two(aline[j], comm_p, myEdgeGraph.edgesOfGraph[aline[j]].T2.p_position, &pl.pointlist);
					comm_p = myEdgeGraph.edgesOfGraph[aline[j]].T2.p_position;
				}
				else
				{
					getPListOfaCurve_two(aline[j], comm_p, myEdgeGraph.edgesOfGraph[aline[j]].T1.p_position, &pl.pointlist);
					comm_p = myEdgeGraph.edgesOfGraph[aline[j]].T1.p_position;
				}
			}
			pointListOfLines.push_back(pl);
		}
	}
}


