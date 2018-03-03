#pragma once
#include "MySvgFile.h"

#define MAX 32768.0
//using namespace std;

struct TerminalPoint{
	mypoint2f endpoint;
	int curveindex;
	char curvetype;
};
struct ArcNode{
	int adjVertex;
	ArcNode *next;
	int curveIndex;
	char curveType;  //'c' 'q' 'l' represent'cubic Bezier', 'quadratic Bezier' and 'polyline'	
	vector<int> patchesIndex;  //empty：一开始是空的，若这条边attach在一个patch上 则不为空
	int lineFlag;  //0:初始值 1：innerline  2:attachline  3:isolatedline 
};
struct VertexNode{
	int verIndex;
	mypoint2f position;
	ArcNode *firstEdge;
	vector<ArcNode> orderedArcs;
};
struct SvgGraph{
	int vertexNum,arcNum;
	//VertexNode *nodeList;
	vector<VertexNode> nodeList;
};
/////////////////////////////////////////////////////////用了前三个
struct TerminalVertex_E{
	mypoint2f p_position;
	//vector<mypoint2f> adjPoints_vec;
	vector<int> adjEdgeIndex_vec;//ordered
};
struct EdgeStruct{
	int index;
	char curveType;
	int curveIndex;
	TerminalVertex_E T1;
	TerminalVertex_E T2;
	vector<int> patchesIndex;
};
struct EdgeGraph{
	vector<EdgeStruct> edgesOfGraph;
	int edgeNumber;
};

struct Patch_two{
	int patch_index;
	vector<int> edgesOfPatch_vec;//edge list
};
struct PLRGraph_two{
	//Point-Line-Region graph
	vector<Patch_two> patches;//vector<vector<int>> patches;
	vector<vector<int>> lines;//vector<vector<mypoint2f>> lines;
};
/*---------------以patch为节点，为构建patch的拓扑结构，声明PatchGraph patchTopology;...-----------------*/
struct PatchArc{
	int adj_patchindex;
	double weight;
	PatchArc *nextPatch;
};
struct PatchNode{
	int patch_index;
	int patchType;//-1-> unavaliable 
	PatchArc* firstAdjPatch;
};
struct PatchGraph{
	vector<PatchNode> allpatches;
	int patch_number;
};
/*-----------（类似于上面的patchgraph）这个是混合排序，以patch/line为节点，为构建混合的拓扑结构，声明mixedGarph...------------*/
struct MixedArc{
	int adj_index;  
	double weight;
	MixedArc *next_arc;
};
struct MixedNode{
	int pl_index;  //1.patch在plrgraph中的index。2.line在otherline中的index + patch的size（记住attachline和isolatedline的数量，就可以判断出这条线是attach还是isolate）
	char pl_type;//n-> unavaliable  p:patch，l:attach  i:isolated
	MixedArc* firstAdjNode;
	int status;//初始值=0，-1被删除 不可用
};
struct MixedGraph{
	vector<MixedNode> allNode;
	int mixnode_number;
};

//-----------------------------------------//



class MyGraph:public MySvgFile
{

private:
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


public:
	MyGraph();
	~MyGraph();
	MyGraph(string filename);

	void constuctEdgeGraph(double w,double h,string file_v);
	void constructPointGraph();
	void combineTwoCurves(int i_ptopl,int tp_k,int tp_j);
	void addAssistantPoint(int i_ptopl, int tp_k, int tp_j);
	bool isPointRoughEqual(mypoint2f p1,mypoint2f p2);
	void createGraph();
	int getAdjVerIndex(mypoint2f point);

	void findLines();//找到outlines和isolatedlines 并从sgraph.nodelist中删除相应的arcnode,  即设置arcnode.innerline标志位为2
	int countEdgeNum(int startVerIndex, int flag);//第二个参数counttype=0表示计算所有的arcnode ，=1表示除了innerline outline isolateline其他的arc个数
	ArcNode *getArcNode(int startVerIndex,int endVerIndex);
	void removeArcNode(int startVerIndex, int endVerIndex,int flag);//从startVerIndex中删除endVerIndex（单向的）, 第三个参数 0：真正的删除,  1/2/3/4：不删除，只是修改标志位innerline为1，2，3
	//void addArcNode(int startver, int endver, int flag);

	void divideSubGraph();//breadth-first traversal  BFT
	void MakeOrderedArcs();

	void findPatches_new();//+1
	void method(int pindex1, int pindex2, int &pcount,bool &cflag,vector<int> &inn_pt);//+1
	void deleOverlayPatches(int dele_pch,int store_pch);//第二个参数没用，可以删了//+1
	void deleRepeatPoint(vector<int> *plist);//+1+1
	int findInnerLine(int graph_p);//+1

	int deleRepeatPatch();//+1 =deleOverlayPatches()
	bool isDifferentPath(vector<int> p1, vector<int> p2);//+1
	void printpathOfPatches();//+1
	void divideAttachLines();//在找完patch之后，吧isolatedline和attachline分开(在findline中统一标记为2)。并且找到attachline对应的那个patch，修改（添加）bool attach_flag;int adjCurve;等信息


	void jointPolylines();
	void jointCurveTPatch();
	void jointMethod(vector<int>* pindex_vec,vector<mypoint2f>* plist,int density);// density=1隔一个去一个点，=2隔2个取一个点，=0无间隔
	void joinCurveToSpecialPatch(vector<int> pch);//为了测试观察某几个patch用
	void getPListOfaCurve(mypoint2f point1, mypoint2f point2, int verindex1, int verindex2, vector<mypoint2f> *vec,int density);//point1->point2
	void getPListfromCubic(int order, cubicBezier aCBcurve , int density,vector<mypoint2f> *vec);
	void getPListfromQuadratic(int order, quaBezier aQBcurve, int density, vector<mypoint2f> *vec);
	void getPListfromPolyline(int order, PolyLine aPLline, int density, vector<mypoint2f> *vec);

	void detectShorterEdge();
	Patch *getAPatch(int patchindex);
	EdgeOfPatch *getBorderofPatch(int patchindex,int startp,int endp);
	void mergeEdge();
	bool getAddpFromAdjpatch(vector<int> pvec_origi, vector<int> pvec_adj, vector<int> &addp,int &startp,int &top);
	vector<int> getNeighborPatches(int patchindex,int flag);//0有相邻边的neighbor patch，1包括对角patch
	double getPerimeterofVec(vector<mypoint2f>* pvec);


	void ConstructPatchTopology();
	void removePatchNode(int hostpch, int delpch);
	PatchArc* getPatchArc(int pindex, int nei_pindex);
	int countPatchEdgeNum(int pchindex);
	void deleUndesirablePatch();//在下面的KmeansClustering()函数之后执行
	void jointCurveTPatch_subgraph();

	void getPatchFeatures();
	void KmeansClustering(int classnum);
	void KmedoidClustering(int classnum);
	void DBSCANClusering();
	void FCMClustering(int classnum,double para_m);//para_m->1隶属度区别越大

	//找到与line平行角度偏差小的neighbor line
	void getOrigiCurveList(char currt_type, int currt_index, char& origin_type, int& origin_index, vector<mypoint2f> &plist,int flag);//flag=0,只获得origin_typeh和 origin_index，flag=1,有plist
	vector<EdgeStruct> getNeiborLines(int nodeindex1, int nodeindex2, int scope, double width_threshold, double angle_threshold);
	//腐蚀成一条线，
	void Erosion_test();
	void getMidPlist(vector<mypoint2f> plist1, vector<mypoint2f> plist2, vector<mypoint2f> *midplist);
	bool isInVector(vector<int> nodevec,int node);
	void truncateAnEdge(vector<mypoint2f> origi_plist,vector<mypoint2f> *cutout, mypoint2f fromp, mypoint2f endp);
	void translateOutline(int startp,ArcNode* anodep,mypoint2f trans_dis);
	void translateEdge(char ctype, int cindex, mypoint2f trans_dis);
	void substituteEdge(char ctype, int cindex, vector<mypoint2f> newplist);
	void transformEdge(vector<mypoint2f> *new_plist, mypoint2f addp, vector<mypoint2f> *origi_plist);
	void modifyPatch_plist(Patch *pch,vector<int> origip,vector<int> targetp);
	void specialTransformEdge(int startp,ArcNode *anode,int origi_pch,vector<int> *left_p,vector<mypoint2f> *leftp_project,vector<int> *used_node);

	//局部操作 尝试 (从这到线是一个整体，先面后线处理，stopcondition：patchnumber和similarity)
	void localProcess();

	//计算相似度
	void patchSimilarity();
	//patch contour的简化，使其有固定的fixed point number 
	void dataReduction(vector<mypoint2f> *plist, vector<int>* p_vec, int pchindex, int fixednum, vector<mypoint2f> *adjust_plist);
	//patch contour的增加，数据点少，则增加 
	void dataIncreasing(vector<mypoint2f> *plist, int fixednum);
	void dataReduction_byheap(vector<mypoint2f>* plist, int fixednum, double perimeter, vector<mypoint2f> *adjust_plist);
	//当用radial function计算FD时，用等距采样来简化得到固定的pointnumber
	void equalSampling(vector<mypoint2f> *plist, vector<int>* p_vec, int pchindex, int fixednum, vector<mypoint2f> *adjust_plist);

	//方法一：计算bend angle Function,先等距重采样
	void getbendAngleFunction(vector<mypoint2f> *adjusted_plist, vector<vector<double>>* bendAf, int sample_num);
	//方法二：计算傅里叶描述子，输入contour point，输出一个一维的向量FD
	void getFourierDescriptor(vector<mypoint2f>* contour, vector<complex<double>>* FD, int pch);
	

	void localProcess_two(int simply_num);	
	//M1
	double localProcess_heap(int final_pnum, double similarity_t, int simply_num, int iteration);//用heap， 每次删除一个patch，例如把A与B合并，B的weight需要重新计算，并且更新他在堆中的位置，
	void combineTwoPatch_getpvec(int pch1, int pch2, vector<int>* totalpvec);//计算两个patch合并后的plist,density=0保持原来的采样密度
	void getFourierDescriptor(vector<mypoint2f>* contour, vector<double>* FD, int fd_type);// fd_type=0,用bendangle function然后等采样计算FD，=1用radial function计算FD
	double similarityMeasure(vector<double>* fd1, vector<double>* fd2);
	void patchToCurve(int pchindex, vector<mypoint2f>* newcurve);

	//M2
	double similarityRange(int final_pnum, double similarity_t, int simply_num,int file_num);//按合并后相似度最小的当做weight 进行排序,返回patch的平均宽度
	void getThreeFeature(vector<double>* t_ft, vector<mypoint2f> *plist, double p_area, double p_perimeter);// p_area不知道的写0 +1
	void getThreeFeature_forline(vector<double>* t_ft, vector<mypoint2f> *plist);
	void getFDFeature(vector<mypoint2f>* contour,vector<int>* p_vec, int pchindex, int fixednum, int fd_type, vector<double>* fd_feature);// contour的front!=back +1

	//M3
	void shapeSimilarityRange(int final_pnum, double similarity_t, int simply_num, int file_num);//按照，自己的形状与neighbor的形状的相似度进行排序
	bool checkMergeable(int pch1, int pch2);
	
	/*patch单一操作*/
	bool combineTwoPatch(int pch1,int pch2);//通过一条边将两个patch合并，这条边的两个端点必须是度！=2//保留pch1并修改edge_vec等等，删除pch2
	void erosionApatch(int pch);//腐蚀成一条线(删)
	bool checkErosion(int pch,int p1,int p2,int& anotherpch);//（在erosionApatch()内，删）

	/*线*/
	void attachLineProcess(double width_threshold, double angle_threshold);
	void combineTwoCurves(int lindex1, int lindex2, vector<mypoint2f>* newcurve);
	void fitCubBezier(vector<mypoint2f>* allplist, vector<mypoint2f>* cbcurve, cubicBezier* cubBezier, int flag);//flag=0表示不固定两端点，=1固定两个端点
	void removeCurve(vector<int>* rm_curve,int cindex_plr);//还没用（删）,涉及到 在sgraph中删除arc node 和plrGraph.out/inner/isolines[cindex_plr].type
	bool combineTwoCurves_getplist(int lindex1, int lindex2, vector<mypoint2f>* newcurve,double aver_w);
	void checkIntersect(int lindex,vector<int>* newlvec,vector<int>* neighbor_patch);//新的线可能会和其他patch相交


	/*调整pipeline，面和线统一用 FD 来度量相似度，面和线混合排序*/
	double getAverageWidth(double &w,double &h,double ratio);//用所有patch的width平均值
	double ConstructMixedTopology(double w_ratio);
	MixedArc* getMixedArc(int mindex, int nei_mindex);
	void removeMixedNode(int hostnode, int delnode);
	void getMixedNeighbor(int mindex, vector<int>* m_nei);

	void MixSimilarityRange(int simply_num, double patch_averWidth, int similarity_type, double similarity_t);
	void test_MixSimilarityRange(int simply_num, double patch_averWidth);
	void origin_MixSimilarityRange(int simply_num, double patch_averWidth, int similarity_type, double similarity_t);

	//line和patch的merge--->变形的patch
	bool combinePL_getplist(int pindex, int lindex, vector<mypoint2f>* transform_pch);//pindex和lindex分别是plrGraph.patch/otherlines中的,若line在patch的里面，则不能变形
	void combinePatchLine(int pindex, int lindex);
	bool mergePL_getplist(int pindex, int lindex, vector<mypoint2f>* transform_pch);
	void mergePatchLine(int pindex, int lindex);
	void exchangeline(vector<mypoint2f>* nplist,ArcNode* exc_node);

	//平滑
	void globalSmoothing();
	void piecewiseCurveFittting(vector<int>* pt_vec,double error,vector<cubicBezier>* pw_curve);
	void clearStructure();

	void writeSVGFile(string fname,string outpath,double fwidth,double fheight,vector<double> fview,string view_str);
	void writeSVGFile_step(string fname, string outpath, double fwidth, double fheight, vector<double> fview, string view_str);
	void writesmoothFile(string fname, string outpath, double fwidth, double fheight, vector<double> fview, string view_str);

//about opengl
public:
	double file_width;
	double file_height;
	string file_view;

	void testShowGraph();
	void setGraphInstance();
	static void graphDisplay();

	static void showEdgeGraph();
	static void showLines();
	static void showpatches();
	static void showPatchTopology(vector<pair<int, string>> patch_subgraph);
	static void showClassifiedPatch(string str);
	static void showBoundingBox(vector<vector<mypoint2f>> plist);

	static void CALLBACK beginCallback(GLenum type);
	static void CALLBACK vertexCallback(GLdouble * vertex);
	static void CALLBACK endCallback();

protected:
	static MyGraph *ginstance;

////////////////////////////////////////////////////////////////////   Edge Graph method:
public:
	//void constuctEdgeGraph();//放上面了
	//void constructPointGraph();//放上面了
	void addAssitantP_two(int edge_i, int edge_j, bool &delflag);
	void getPListOfaCurve_two(int edge_index, mypoint2f from_p, mypoint2f to_p, vector<mypoint2f> *vec);
	/*void getPListfromCubic(int order, cubicBezier aCBcurve, vector<mypoint2f> *vec);
	void getPListfromQuadratic(int order, quaBezier aQBcurve, vector<mypoint2f> *vec);
	void getPListfromPolyline(int order, PolyLine aPLline, vector<mypoint2f> *vec);*/
	void findLines_two();
	TerminalVertex_E* getTPofEdge(EdgeStruct e1,mypoint2f tp);
	TerminalVertex_E* getOtherTPofEdge(EdgeStruct e1, mypoint2f tp);
	void removeWholeEdge(EdgeStruct *e1);
	void removeEdgeOfTP(EdgeStruct *eg,mypoint2f tp,int rmEdgeIndex);
	void jointPolylines_two();
	
private:
	//EdgeGraph myEdgeGraph;//放上面了
	PLRGraph_two myplrGraph;
public:
	
};

