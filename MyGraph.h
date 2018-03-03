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
	vector<int> patchesIndex;  //empty��һ��ʼ�ǿյģ���������attach��һ��patch�� ��Ϊ��
	int lineFlag;  //0:��ʼֵ 1��innerline  2:attachline  3:isolatedline 
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
/////////////////////////////////////////////////////////����ǰ����
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
/*---------------��patchΪ�ڵ㣬Ϊ����patch�����˽ṹ������PatchGraph patchTopology;...-----------------*/
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
/*-----------�������������patchgraph������ǻ��������patch/lineΪ�ڵ㣬Ϊ������ϵ����˽ṹ������mixedGarph...------------*/
struct MixedArc{
	int adj_index;  
	double weight;
	MixedArc *next_arc;
};
struct MixedNode{
	int pl_index;  //1.patch��plrgraph�е�index��2.line��otherline�е�index + patch��size����סattachline��isolatedline���������Ϳ����жϳ���������attach����isolate��
	char pl_type;//n-> unavaliable  p:patch��l:attach  i:isolated
	MixedArc* firstAdjNode;
	int status;//��ʼֵ=0��-1��ɾ�� ������
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
	static double errorDeviation;//û���õ�����ɾ
	map<int, string> subgraph_map;//record tree's or loop's start index	
	PLRGraph plrGraph;
	vector<pair<int, vector<mypoint2f>>> pointListOfPatches;
	vector<PolyLine> pointListOfLines;

	vector<pair<int, string>> patch_subgraphs;
	PatchGraph patchTopology;
	vector<pair<int, string>> patch_subgraphs_group;
	vector<pair<int, vector<double>>> patchfeatures;//int��Ӧpatchindex��  vector<double>��Ӧ��������

	vector<PolyLine> neighbor_plist_group;//���� ��ȡһ���ߵ�neighborline������topology����getNeiborLines()���õ�
	vector<pair<int, vector<mypoint2f>>> patch_simplified_plist;//���ԣ�����patch���й̶������ĵ��plist
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

	void findLines();//�ҵ�outlines��isolatedlines ����sgraph.nodelist��ɾ����Ӧ��arcnode,  ������arcnode.innerline��־λΪ2
	int countEdgeNum(int startVerIndex, int flag);//�ڶ�������counttype=0��ʾ�������е�arcnode ��=1��ʾ����innerline outline isolateline������arc����
	ArcNode *getArcNode(int startVerIndex,int endVerIndex);
	void removeArcNode(int startVerIndex, int endVerIndex,int flag);//��startVerIndex��ɾ��endVerIndex������ģ�, ���������� 0��������ɾ��,  1/2/3/4����ɾ����ֻ���޸ı�־λinnerlineΪ1��2��3
	//void addArcNode(int startver, int endver, int flag);

	void divideSubGraph();//breadth-first traversal  BFT
	void MakeOrderedArcs();

	void findPatches_new();//+1
	void method(int pindex1, int pindex2, int &pcount,bool &cflag,vector<int> &inn_pt);//+1
	void deleOverlayPatches(int dele_pch,int store_pch);//�ڶ�������û�ã�����ɾ��//+1
	void deleRepeatPoint(vector<int> *plist);//+1+1
	int findInnerLine(int graph_p);//+1

	int deleRepeatPatch();//+1 =deleOverlayPatches()
	bool isDifferentPath(vector<int> p1, vector<int> p2);//+1
	void printpathOfPatches();//+1
	void divideAttachLines();//������patch֮�󣬰�isolatedline��attachline�ֿ�(��findline��ͳһ���Ϊ2)�������ҵ�attachline��Ӧ���Ǹ�patch���޸ģ���ӣ�bool attach_flag;int adjCurve;����Ϣ


	void jointPolylines();
	void jointCurveTPatch();
	void jointMethod(vector<int>* pindex_vec,vector<mypoint2f>* plist,int density);// density=1��һ��ȥһ���㣬=2��2��ȡһ���㣬=0�޼��
	void joinCurveToSpecialPatch(vector<int> pch);//Ϊ�˲��Թ۲�ĳ����patch��
	void getPListOfaCurve(mypoint2f point1, mypoint2f point2, int verindex1, int verindex2, vector<mypoint2f> *vec,int density);//point1->point2
	void getPListfromCubic(int order, cubicBezier aCBcurve , int density,vector<mypoint2f> *vec);
	void getPListfromQuadratic(int order, quaBezier aQBcurve, int density, vector<mypoint2f> *vec);
	void getPListfromPolyline(int order, PolyLine aPLline, int density, vector<mypoint2f> *vec);

	void detectShorterEdge();
	Patch *getAPatch(int patchindex);
	EdgeOfPatch *getBorderofPatch(int patchindex,int startp,int endp);
	void mergeEdge();
	bool getAddpFromAdjpatch(vector<int> pvec_origi, vector<int> pvec_adj, vector<int> &addp,int &startp,int &top);
	vector<int> getNeighborPatches(int patchindex,int flag);//0�����ڱߵ�neighbor patch��1�����Խ�patch
	double getPerimeterofVec(vector<mypoint2f>* pvec);


	void ConstructPatchTopology();
	void removePatchNode(int hostpch, int delpch);
	PatchArc* getPatchArc(int pindex, int nei_pindex);
	int countPatchEdgeNum(int pchindex);
	void deleUndesirablePatch();//�������KmeansClustering()����֮��ִ��
	void jointCurveTPatch_subgraph();

	void getPatchFeatures();
	void KmeansClustering(int classnum);
	void KmedoidClustering(int classnum);
	void DBSCANClusering();
	void FCMClustering(int classnum,double para_m);//para_m->1����������Խ��

	//�ҵ���lineƽ�нǶ�ƫ��С��neighbor line
	void getOrigiCurveList(char currt_type, int currt_index, char& origin_type, int& origin_index, vector<mypoint2f> &plist,int flag);//flag=0,ֻ���origin_typeh�� origin_index��flag=1,��plist
	vector<EdgeStruct> getNeiborLines(int nodeindex1, int nodeindex2, int scope, double width_threshold, double angle_threshold);
	//��ʴ��һ���ߣ�
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

	//�ֲ����� ���� (���⵽����һ�����壬������ߴ���stopcondition��patchnumber��similarity)
	void localProcess();

	//�������ƶ�
	void patchSimilarity();
	//patch contour�ļ򻯣�ʹ���й̶���fixed point number 
	void dataReduction(vector<mypoint2f> *plist, vector<int>* p_vec, int pchindex, int fixednum, vector<mypoint2f> *adjust_plist);
	//patch contour�����ӣ����ݵ��٣������� 
	void dataIncreasing(vector<mypoint2f> *plist, int fixednum);
	void dataReduction_byheap(vector<mypoint2f>* plist, int fixednum, double perimeter, vector<mypoint2f> *adjust_plist);
	//����radial function����FDʱ���õȾ�������򻯵õ��̶���pointnumber
	void equalSampling(vector<mypoint2f> *plist, vector<int>* p_vec, int pchindex, int fixednum, vector<mypoint2f> *adjust_plist);

	//����һ������bend angle Function,�ȵȾ��ز���
	void getbendAngleFunction(vector<mypoint2f> *adjusted_plist, vector<vector<double>>* bendAf, int sample_num);
	//�����������㸵��Ҷ�����ӣ�����contour point�����һ��һά������FD
	void getFourierDescriptor(vector<mypoint2f>* contour, vector<complex<double>>* FD, int pch);
	

	void localProcess_two(int simply_num);	
	//M1
	double localProcess_heap(int final_pnum, double similarity_t, int simply_num, int iteration);//��heap�� ÿ��ɾ��һ��patch�������A��B�ϲ���B��weight��Ҫ���¼��㣬���Ҹ������ڶ��е�λ�ã�
	void combineTwoPatch_getpvec(int pch1, int pch2, vector<int>* totalpvec);//��������patch�ϲ����plist,density=0����ԭ���Ĳ����ܶ�
	void getFourierDescriptor(vector<mypoint2f>* contour, vector<double>* FD, int fd_type);// fd_type=0,��bendangle functionȻ��Ȳ�������FD��=1��radial function����FD
	double similarityMeasure(vector<double>* fd1, vector<double>* fd2);
	void patchToCurve(int pchindex, vector<mypoint2f>* newcurve);

	//M2
	double similarityRange(int final_pnum, double similarity_t, int simply_num,int file_num);//���ϲ������ƶ���С�ĵ���weight ��������,����patch��ƽ�����
	void getThreeFeature(vector<double>* t_ft, vector<mypoint2f> *plist, double p_area, double p_perimeter);// p_area��֪����д0 +1
	void getThreeFeature_forline(vector<double>* t_ft, vector<mypoint2f> *plist);
	void getFDFeature(vector<mypoint2f>* contour,vector<int>* p_vec, int pchindex, int fixednum, int fd_type, vector<double>* fd_feature);// contour��front!=back +1

	//M3
	void shapeSimilarityRange(int final_pnum, double similarity_t, int simply_num, int file_num);//���գ��Լ�����״��neighbor����״�����ƶȽ�������
	bool checkMergeable(int pch1, int pch2);
	
	/*patch��һ����*/
	bool combineTwoPatch(int pch1,int pch2);//ͨ��һ���߽�����patch�ϲ��������ߵ������˵�����Ƕȣ�=2//����pch1���޸�edge_vec�ȵȣ�ɾ��pch2
	void erosionApatch(int pch);//��ʴ��һ����(ɾ)
	bool checkErosion(int pch,int p1,int p2,int& anotherpch);//����erosionApatch()�ڣ�ɾ��

	/*��*/
	void attachLineProcess(double width_threshold, double angle_threshold);
	void combineTwoCurves(int lindex1, int lindex2, vector<mypoint2f>* newcurve);
	void fitCubBezier(vector<mypoint2f>* allplist, vector<mypoint2f>* cbcurve, cubicBezier* cubBezier, int flag);//flag=0��ʾ���̶����˵㣬=1�̶������˵�
	void removeCurve(vector<int>* rm_curve,int cindex_plr);//��û�ã�ɾ��,�漰�� ��sgraph��ɾ��arc node ��plrGraph.out/inner/isolines[cindex_plr].type
	bool combineTwoCurves_getplist(int lindex1, int lindex2, vector<mypoint2f>* newcurve,double aver_w);
	void checkIntersect(int lindex,vector<int>* newlvec,vector<int>* neighbor_patch);//�µ��߿��ܻ������patch�ཻ


	/*����pipeline�������ͳһ�� FD ���������ƶȣ�����߻������*/
	double getAverageWidth(double &w,double &h,double ratio);//������patch��widthƽ��ֵ
	double ConstructMixedTopology(double w_ratio);
	MixedArc* getMixedArc(int mindex, int nei_mindex);
	void removeMixedNode(int hostnode, int delnode);
	void getMixedNeighbor(int mindex, vector<int>* m_nei);

	void MixSimilarityRange(int simply_num, double patch_averWidth, int similarity_type, double similarity_t);
	void test_MixSimilarityRange(int simply_num, double patch_averWidth);
	void origin_MixSimilarityRange(int simply_num, double patch_averWidth, int similarity_type, double similarity_t);

	//line��patch��merge--->���ε�patch
	bool combinePL_getplist(int pindex, int lindex, vector<mypoint2f>* transform_pch);//pindex��lindex�ֱ���plrGraph.patch/otherlines�е�,��line��patch�����棬���ܱ���
	void combinePatchLine(int pindex, int lindex);
	bool mergePL_getplist(int pindex, int lindex, vector<mypoint2f>* transform_pch);
	void mergePatchLine(int pindex, int lindex);
	void exchangeline(vector<mypoint2f>* nplist,ArcNode* exc_node);

	//ƽ��
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
	//void constuctEdgeGraph();//��������
	//void constructPointGraph();//��������
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
	//EdgeGraph myEdgeGraph;//��������
	PLRGraph_two myplrGraph;
public:
	
};

