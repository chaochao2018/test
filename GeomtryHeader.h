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
	int origi_index;//��ĳ���������������ཻ�����ָ�ɲ�ͬ�����߶Σ���ñ�������δ�ָ�֮ǰ�����ߵ����ͣ����ڱ�����������˵��û���������
	char origi_type; //��������������mysvgfile�е�d��getGraphPrimitive�г�ʼ��,��splitcurve�м����µġ�����cubicbezier_vec[5]���߱��ָ�µ�polyline����polyine_vec����origi_index=5,origi_type=c!!!
};


//struct PointofPatch{
//	int pindex;
//	vector<int> attachCurve;//empty�� û���ڽӵ�line��1��innerline��2��attachline��3��isolate������ArcNode�е�lineflagһ��0:��ʼֵ 1��innerline  2:attachline  3:isolatedline 
//	PointofPatch(){
//		this->pindex = 0;
//	}
//	PointofPatch(int pch){
//		this->pindex = pch;
//	}
//};
struct EdgeOfPatch{
	int startp_index;//PLRGraph plrGraph.patch[i] �� SvgGraph sgraph֮���������  �����startp_index��Ӧsgraph.nodeList���±�) ͬ��pointIndexOfPatch
	int endp_index;
	vector<int> adjPatchIndex;//adjacent patch's index;
	int edgeflag;//��ʼ��ֵ0����ֵΪ1��ʾ����̱�
};
struct Patch{
	int patch_index;//patch index
	int type;//��ʼ��ֵ0��������-1˵�������ã��ظ���patch�Ѿ����ϲ���
	vector<int> pointIndexOfPatch;//point list�����patch�ĵ�,��ÿ���㶼��һ���Ƿ���attachline�ı�־λ����attachline���±�
	vector<EdgeOfPatch> edges_vec;//edge list�����patch��edge
	vector<pair<int,int>> point_line;//point��Ӧ����±꣬line�������������һ��ͻ�����ߣ�line��index��ֵʱ�õ���attachline���±꣬Ҳ������otherlines������,��Ϊotherline=�ȴ�attachline���ٴ�isolatedline

	double perimeter;
	double area;
	mypoint2f centre_position;
};
struct ICurve{//��patch��ƽ�ȵ� ����findline()��findinnerline()�����ң���findpatch()���޸�
	vector<int> pt_vec;
	int type;//��ʼֵ=0   ɾ����=-1
};
struct PLRGraph{
	//Point-Line-Region graph
	vector<Patch> patches;//vector<vector<int>> patches;
	vector<mypoint2f> points;

	vector<ICurve> attachlines;//vector<vector<mypoint2f>> lines;
	vector<ICurve> innerlines;
	vector<ICurve> isolatedlines;
	vector<ICurve> otherlines;//�ȴ�attachline���ٴ�isolatedline����Ϊpatch��point_line��line�� index�� attachline���±�
};
//����dbscan��������ݽṹ
struct DataPoint{
	int dataID;          //���ݵ�ID
	vector<double> pch_ft;        //����
	int clusterId;       //��������ID
	int datatype;        //1��ʾ�����㣬2��ʾ�߽磬3��ʾ���Ķ���
	bool visited;         //�Ƿ��ѷ���
	vector<int> neigbor_ct;   //�ڰ뾶radius ֮�ڵ�neighbor �ĸ���
};
int dbscanclustering(vector<DataPoint> &dataset, double radius, int minpts);//���ؾ���ĸ���
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

//�Զ���С���� ����
//�ϸ����³���updata
//a��ʾ��Ҫ�ϵ���data��λ�ã�aȡֵ��Χ1...N, ����������
int heap_siftup(int a, vector<int>* b,vector<double>* ele,vector<int>* h_index);
//a��ʾ��Ҫ�µ���data��λ�ã�aȡֵ��Χ1...N, ������λ��
int heap_siftdown(int a, vector<int>* b, vector<double>* ele, vector<int>* h_index);
//����λ��index����dataΪnewdata��indexȡֵ��Χ1...N , ������λ��
int heap_update(int index, double newdata, vector<int>* b, vector<double>* ele, vector<int>* h_index);


void getBoundRectangle(vector<mypoint2f>* total_plist, mypoint2f &centerpt, vector<mypoint2f>* four_corner);

//piecewise curve fitting
mypoint2f V2Scale(mypoint2f *p, double newlen);
mypoint2f V2ScaleIII(mypoint2f v, double s);
mypoint2f V2AddII(mypoint2f a, mypoint2f b);
mypoint2f V2SubII(mypoint2f a, mypoint2f b);
mypoint2f *V2Negate(mypoint2f *v);

mypoint2f ComputeLeftTangent(vector<mypoint2f>* d, int end);//����d->at(end + 1) - d->at(end)�ĵ�λ������;
mypoint2f ComputeRightTangent(vector<mypoint2f>* d, int end);
mypoint2f ComputeCenterTangent(vector<mypoint2f>* d, int center);

void ChordLengthParameterize(vector<mypoint2f>* d, int first, int last, vector<double> *u);//���ز���u
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
