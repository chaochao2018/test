#pragma once

#include<vector>
#include<queue>
#include<list>
#include<algorithm>
#include<math.h>
#include<stack>
#include<fstream>
#include<complex>
#include"Eigen\Dense"
#include<tbb/tbb.h>

using namespace Eigen;
using namespace std;

extern double eps;
extern double eps_isOnsegment;

struct mypoint2f{
	double x;
	double y;
	mypoint2f() :x(0), y(0){}
	mypoint2f(double xvalue, double yvalue) :x(xvalue), y(yvalue){}
	mypoint2f(const mypoint2f& ap)
	{
		x = ap.x;
		y = ap.y;
	}
	bool operator ==(const mypoint2f p1)
	{
		if (fabs(p1.x - this->x)<eps && fabs(p1.y - this->y)<eps)//if (p1.x == this->x&&p1.y == this->y)
			return true;
		else
			return false;
	}
	bool operator !=(const mypoint2f p1)
	{
		if ((p1.x != this->x)&& (p1.y != this->y))//if (p1.x == this->x&&p1.y == this->y)
			return true;
		else
			return false;
	}
	mypoint2f operator + (const mypoint2f &a) const
	{
		return mypoint2f(x + a.x, y + a.y);
	}
	mypoint2f operator - (const mypoint2f &a) const
	{
		return mypoint2f(x - a.x, y - a.y);
	}
};
struct mycircle{
	double cx;
	double cy;
	double radius;
};
struct myrect{
	double x;    double y;
	double width; double height;
};
struct cubicBezier{
	mypoint2f p0;
	mypoint2f p1;
	mypoint2f p2;
	mypoint2f p3;
	//int origi_index; 
	//char origi_type;
};
struct quaBezier{
	mypoint2f p0;
	mypoint2f p1;
	mypoint2f p2;
	//int origi_index;
	//char origi_type;
};
struct PolyLine{
	vector<mypoint2f> pointlist;
	int origi_index;//若某条曲线与其他线相交而被分割成不同的折线段，则该变量保存未分割之前的曲线的类型，对于贝塞尔曲线来说，没有这个变量
	char origi_type; //这两个变量，在mysvgfile中的d的getGraphPrimitive中初始化,在splitcurve中加入新的。例如cubicbezier_vec[5]曲线被分割，新的polyline加入polyine_vec，则origi_index=5,origi_type=c!!!
};


//struct PointofPatch{
//	int pindex;
//	vector<int> attachCurve;//empty： 没有邻接的line，1：innerline，2：attachline，3：isolate。。和ArcNode中的lineflag一样0:初始值 1：innerline  2:attachline  3:isolatedline 
//	PointofPatch(){
//		this->pindex = 0;
//	}
//	PointofPatch(int pch){
//		this->pindex = pch;
//	}
//};
struct EdgeOfPatch{
	int startp_index;//PLRGraph plrGraph.patch[i] 和 SvgGraph sgraph之间的桥梁吧  （这个startp_index对应sgraph.nodeList的下标) 同理pointIndexOfPatch
	int endp_index;
	vector<int> adjPatchIndex;//adjacent patch's index;
	int edgeflag;//初始赋值0，赋值为1表示是最短边
};
struct Patch{
	int patch_index;//patch index
	int type;//初始赋值0，若等于-1说明不可用（重复或patch已经被合并）
	vector<int> pointIndexOfPatch;//point list，组成patch的点,且每个点都有一个是否有attachline的标志位，和attachline的下标
	vector<EdgeOfPatch> edges_vec;//edge list，组成patch的edge
	vector<pair<int,int>> point_line;//point对应点的下标，line是这个点连接有一条突出的线，line的index赋值时用的是attachline的下标，也可以用otherlines来查找,因为otherline=先存attachline，再存isolatedline

	double perimeter;
	double area;
	mypoint2f centre_position;
};
struct ICurve{//跟patch是平等的 ，在findline()，findinnerline()和中找，在findpatch()中修改
	vector<int> pt_vec;
	int type;//初始值=0   删除了=-1
};
struct PLRGraph{
	//Point-Line-Region graph
	vector<Patch> patches;//vector<vector<int>> patches;
	vector<mypoint2f> points;

	vector<ICurve> attachlines;//vector<vector<mypoint2f>> lines;
	vector<ICurve> innerlines;
	vector<ICurve> isolatedlines;
	vector<ICurve> otherlines;//先存attachline，再存isolatedline，因为patch的point_line中line的 index是 attachline的下标
};
//用于dbscan聚类的数据结构
struct DataPoint{
	int dataID;          //数据点ID
	vector<double> pch_ft;        //数据
	int clusterId;       //所属聚类ID
	int datatype;        //1表示噪声点，2表示边界，3表示核心对象
	bool visited;         //是否已访问
	vector<int> neigbor_ct;   //在半径radius 之内的neighbor 的个数
};
int dbscanclustering(vector<DataPoint> &dataset, double radius, int minpts);//返回聚类的个数
//void kmeans(vector<pair<int, double>> &type_dis, int classnum, vector<pair<int, vector<double>>> pchFeathers);



int dcmp(double x);
double dot(mypoint2f a, mypoint2f b);
double cross(mypoint2f a, mypoint2f b);
double getEuclidDistance(vector<double>* vect1,vector<double>* vect2);
mypoint2f normalization(mypoint2f from_a, mypoint2f to_b);
double getCosineAngle(mypoint2f from_a1, mypoint2f to_a2, mypoint2f from_b1, mypoint2f to_b2);
bool isSegmentIntersected(mypoint2f a1, mypoint2f a2, mypoint2f b1, mypoint2f b2);
bool isOnSegment(mypoint2f p, mypoint2f a1, mypoint2f a2);
double pointToLineDistance(mypoint2f p, mypoint2f lp1, mypoint2f lp2);
double getAtan(mypoint2f p);
double getAtan2(mypoint2f p);
double getProximity(vector<mypoint2f>* plist1, vector<mypoint2f>* plist2,int &posi1,int &posi2);
double getContinuity(vector<mypoint2f>* plist1, vector<mypoint2f>* plist2,int posi1,int posi2);



double Area(const mypoint2f& p1, const mypoint2f& p2, const mypoint2f& p3);
double fArea(const mypoint2f& p1, const mypoint2f& p2, const mypoint2f& p3);
mypoint2f getIntersectPoint(const mypoint2f& p1, const mypoint2f& p2, const mypoint2f& p3, const mypoint2f& p4);


vector<double> quadricRealPolynomial(double coeff_a, double coeff_b, double coeff_c);
vector<double> cubicRealPolynomial(double coeff_a, double coeff_b, double coeff_c, double coeff_d);

double LineSegmentLength(mypoint2f p1,mypoint2f p2);
void CubicBezrCoefficient(double t,double &a1,double &a2,double &a3,double &a4);
void CubicBezrDerivative(double t, double &b1, double &b2, double &b3, double &b4);
void QuaBezrCoefficient(double t, double &a1, double &a2, double &a3);
void QuaBezrDerivative(double t, double &b1, double &b2, double &b3);
//mypoint2f GetPointFromCBezr(double a1,double a2,double a3,double a4,double p0,double p1,double p2,double p3);
//mypoint2f GetPointFromQBezr(double t, mypoint2f &point);


double GetPolygonArea(vector<mypoint2f>* plist);

double GetAverange(vector<double> vec);
double GetVariance(vector<double> vec1, double averange1, vector<double> vec2, double averange2);
double GetVariance_xx(vector<double> vec1, double averange1);

bool mycompare(const pair<int, double>& v1, const pair<int, double>& v2);


bool checkAcuteTriangle(mypoint2f p1,mypoint2f p2,mypoint2f p3);
void makeAcuteTrangle(vector<mypoint2f>* plist,mypoint2f *p1, mypoint2f *p2, mypoint2f *p3);
void maxInscribedCircle(vector<mypoint2f> *polygon_vec,mypoint2f &center,double &radius);

//自定义小顶堆 操作
//上浮，下沉，updata
//a表示需要上调的data的位置，a取值范围1...N, 更新索引表
int heap_siftup(int a, vector<int>* b,vector<double>* ele,vector<int>* h_index);
//a表示需要下调的data的位置，a取值范围1...N, 返回新位置
int heap_siftdown(int a, vector<int>* b, vector<double>* ele, vector<int>* h_index);
//更新位置index处的data为newdata，index取值范围1...N , 返回新位置
int heap_update(int index, double newdata, vector<int>* b, vector<double>* ele, vector<int>* h_index);


void getBoundRectangle(vector<mypoint2f>* total_plist, mypoint2f &centerpt, vector<mypoint2f>* four_corner);

//piecewise curve fitting
mypoint2f V2Scale(mypoint2f *p, double newlen);
mypoint2f V2ScaleIII(mypoint2f v, double s);
mypoint2f V2AddII(mypoint2f a, mypoint2f b);
mypoint2f V2SubII(mypoint2f a, mypoint2f b);
mypoint2f *V2Negate(mypoint2f *v);

mypoint2f ComputeLeftTangent(vector<mypoint2f>* d, int end);//计算d->at(end + 1) - d->at(end)的单位切向量;
mypoint2f ComputeRightTangent(vector<mypoint2f>* d, int end);
mypoint2f ComputeCenterTangent(vector<mypoint2f>* d, int center);

void ChordLengthParameterize(vector<mypoint2f>* d, int first, int last, vector<double> *u);//返回参数u
double B0(double u);
double B1(double u);
double B2(double u);
double B3(double u);
cubicBezier GenerateBezier(vector<mypoint2f>* d, int first, int last, vector<double>* uPrime, mypoint2f tHat1, mypoint2f tHat2);
mypoint2f BezierII(int degree, vector<mypoint2f>* V, double t);
double ComputeMaxError(vector<mypoint2f>* d, int first, int last, cubicBezier bezCurve, vector<double>* u, int* splitPoint);
double NewtonRaphsonRootFind(cubicBezier Q, mypoint2f P, double u);
void Reparameterize(vector<mypoint2f>* d, int first, int last, vector<double>* u, cubicBezier bezCurve, vector<double> *uPrime);

void FitCubic(vector<mypoint2f>* d, int first, int last, mypoint2f tHat1, mypoint2f tHat2, double error, vector<cubicBezier>* piecewise_cb);
void FitCurves(vector<mypoint2f>* d, int nPts, double error, vector<cubicBezier>* piecewise_cb);
