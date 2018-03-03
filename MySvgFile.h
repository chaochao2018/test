#pragma once

#include"tinystr.h"
#include"tinyxml.h"
#include<iostream>
#include <string.h>
#include<sstream>
#include<map>
//#include<stack>
#include<ctime>
#include<algorithm>
#include<GL\glut.h>
#include"GeomtryHeader.h"
//#include <GL\freeglut.h>

//using namespace std;

//
//struct mypoint2f{
//	double x;
//	double y;
//	mypoint2f():x(0), y(0){}
//	mypoint2f(double xvalue, double yvalue) :x(xvalue), y(yvalue){}
//	mypoint2f(const mypoint2f& ap)
//	{
//		x = ap.x;
//		y = ap.y;
//	}
//	bool operator ==(const mypoint2f p1) 
//	{
//		if (p1.x == this->x&&p1.y == this->y)
//		{
//			return true;
//		}
//		else
//			return false;
//	}
//};
//struct mycircle{
//	double cx;
//	double cy;
//	double radius;
//};
//struct myrect{
//	double x;    double y;
//	double width; double height;
//};

struct DiffCurve{
	int index;
	char curveType;
	int curveIndex;
	vector<mypoint2f> split_point;
	vector<double> split_t;
	vector<int> itsct_cur;
};


class MySvgFile
{
public:
	MySvgFile();
	~MySvgFile();
	MySvgFile(string name);
	void getGraphPrimitive(double *fwidth,double *fheight,vector<double> *view,string &view_str);
	void getDifferentGraphics(TiXmlElement *graphelement,char *g_type);//get different graphics from path,circle,rect,line....And print their type(for checking)
	void getCrossCurve();
	//show the graph 
	void testShowDiffvec();

	//测试PCA
	void TestPCA();   //PCA + 3次 Bezier曲线
	//vector<mypoint2f> ordered_plist;
	void LaplacianEigenmap();
	//average curve
	void TestAverageCurve();
	mypoint2f algorithm_one(vector<vector<mypoint2f>> normal_list, vector<vector<mypoint2f>> line_list, mypoint2f apoint,vector<int> &start_GRP_index,mypoint2f& N_newline);
	vector<int> computeGRP(mypoint2f point, vector<int> startindex, vector<vector<mypoint2f>> line_list);
	mypoint2f snap(vector<vector<mypoint2f>> normal_list, vector<vector<mypoint2f>> line_list, mypoint2f p, vector<int>& grp_vec, mypoint2f& normal_p);
	//测试分段Beziercurvefitting
	void piecewiseBezieFitting();
	void writefile();

private:
	double getFileSize(string str);
	int shapeStrTNum(char *shapename);
	void getTransform(TiXmlElement *graphelement);
	vector<mypoint2f> transOrigiCooridate(vector<mypoint2f> tp_vec);
	mypoint2f transOrigiCooridate(mypoint2f tp);
	vector<double> strToNum(string str,int startp,int endp);
	string trimSpace(string str);
	void divideShapeFromPath(TiXmlElement *graphelement);
	//<path>
	mypoint2f path_CubicBezier(mypoint2f &startp, vector<double> number_vec, char type);
	void path_SmoothCubBezier(mypoint2f &startp, mypoint2f lastp, vector<double> number_vec, char type);
	mypoint2f path_QuadraticBezier(mypoint2f &startp, vector<double> number_vec, char type);
	void path_SmoothQuaBezier(mypoint2f &startp, mypoint2f lastp, vector<double> number_vec, char type);
	void path_LineTo(mypoint2f &startp, vector<double> number_vec, char type);
	void path_HorizontalTo(mypoint2f &startp, vector<double> number_vec, char type);
	void path_VerticalTo(mypoint2f &startp, vector<double> number_vec, char type);
	//void path_MoveTo();
	//<line>
	void shape_Line(TiXmlElement *graphelement);//only two points
	//<circle>
	void shape_Circle(TiXmlElement *graphelement);
	//<rect>
	void shape_Rect(TiXmlElement *graphelement);
	//<polyline>
	void shape_Polyline(TiXmlElement *graphelement);
	//<polygon>
	void shape_Polygon(TiXmlElement *graphelement);


	//compute crossing point between different curves
	void coarseJudgeByConvex(int c1,int c2, vector<DiffCurve> *diffcur);
	vector<mypoint2f> getConvexPoints(DiffCurve curve);
	bool isBoudEncircled(vector<mypoint2f> convex1, vector<mypoint2f> convex2);
	vector<pair<mypoint2f, mypoint2f>> getConvexSegment(DiffCurve curve);

	//bool isInCrossCurvePair(DiffCurve curve1, DiffCurve curve2, vector<pair<DiffCurve, DiffCurve>> crosscur_vec);
	//bool isOnPolylineSegmt(DiffCurve polycur, DiffCurve cmpcurve);

	void getCrossCurvePosi(vector<DiffCurve> *diffcurlist, vector<bool> *isSplited);
	bool isTerpointOnCurve(DiffCurve curve1, DiffCurve curve2, vector<double> &split_t, vector<mypoint2f> &split_p, double &mini_dis);
	vector<mypoint2f> getTerpoint(DiffCurve curve);
	void getExtendedLine(DiffCurve cross_curve, DiffCurve ext_curve, mypoint2f ext_p, PolyLine &ext_line,int scalor);
	void addSplitPosi(vector<DiffCurve> *diffcurlist, int cur1, int cur2, vector<double> split_t1, vector<double> split_t2, vector<mypoint2f> split_p1, vector<mypoint2f> split_p2);

	bool twoDiffCurveIntersct(DiffCurve cur1, DiffCurve cur2, vector<double> &split_t1, vector<double> &split_t2, vector<mypoint2f> &split_p1, vector<mypoint2f> &split_p2);
	bool checkTValidity(vector<double> &split_t);
	bool checkPValidity(vector<mypoint2f> &split_p, mypoint2f sp, mypoint2f ep);
	bool twoCubicBzrIntersect(cubicBezier cubBcur1, cubicBezier cubBcur2,vector<double> &split_t1,vector<double> &split_t2);
	bool twoQuaBzrIntersect(quaBezier quaBcur1, quaBezier quaBcur2, vector<double> &split_t1, vector<double> &split_t2);
	bool CBandQBIntersect(cubicBezier Cbcur, quaBezier Qbcur, vector<double> &split_t1, vector<double> &split_t2);
	bool CBandPLIntersect(cubicBezier Cbcur, vector<mypoint2f> polyline, vector<double> &split_t1, vector<mypoint2f> &split_p);
	bool QBandPLIntersect(quaBezier qbcur, vector<mypoint2f> polyline, vector<double> &split_t1, vector<mypoint2f> &split_p);
	bool twoPolyIntersect(vector<mypoint2f> p1, vector<mypoint2f> p2, vector<mypoint2f> &split_p1, vector<mypoint2f> &split_p2);
	//并行的计算相交
	bool twoCubicBzrIntersect_para(int para_i, cubicBezier cubBcur1, cubicBezier cubBcur2, vector<double> &split_t1, vector<double> &split_t2);

	//void splitACurve(DiffCurve curve, double split_t, mypoint2f split_p);
	void splitACurve(DiffCurve curve,vector<double> &split_t,vector<mypoint2f> &split_p,int density);
	void splitCubBzrs(vector<double> split_t, cubicBezier cbcurve,int origi_index, vector<PolyLine> &add_polyline, int density);
	void splitQuaBzrs(vector<double> split_t, quaBezier qbcurve , int origi_index, vector<PolyLine> &add_polyline, int density);
	void splitPolyline(vector<mypoint2f> splitps, vector<mypoint2f> origiline, int origi_index, vector<PolyLine> &add_polyline);

	void deleOrigiLongCurve(vector<DiffCurve> *diffcurlist, vector<bool> *marksplit);
	void makeTerpointToCircle();

	
private:
	string filename;
	TiXmlDocument *filedoc;
	
	int tree_depth;
	//static double errorDeviation;	
	vector<pair<string, vector<double>>> transform_vec;//map<string, vector<double>> transform_vec;

public:
	vector<cubicBezier> cubicBezier_vec_origin;//这三个是从svg文件中解析出来的 最初的曲线。可能会在另一个文件中用到，public类型
	vector<quaBezier> quadraticBezier_vec_origin;
	vector<PolyLine> polyline_vec_origin;

	vector<PolyLine> added_polyline;//相交曲线分段，分段后的折线存储在这个变量中

	vector<cubicBezier> cubicBezier_vec;//这三个是将相交的曲线分段后的存储不同曲线的变量，即删除了已经分段的贝塞尔曲线，增添了新的折线
	vector<quaBezier> quadraticBezier_vec;//以后基本上都是使用这三个变量进行操作。若想找到被分割以前的原始线条 用上面的xxx_origin
	vector<PolyLine> polyline_vec;

	vector<vector< mypoint2f>> polygon_vec;//这些变量可能用不到
	vector<mycircle> circle_vec;
	vector<myrect> rect_vec;

	vector<PolyLine> polyline_extend;// 测试用，
	vector<cubicBezier> piecewise_Bezier;  //测试用

//about opengl 
public:
	static float zoomup;
	static float zoomdown;
	static int winWidth;
	static float winHeight;
	static float w_x1, w_x2, w_y1, w_y2;
	/*const static float zoomup;
	const static float zoomdown;
	const static int winWidth;
	const static float winHeight;
	const static float w_x1, w_x2, w_y1, w_y2;*/
protected:
	static MySvgFile *svginstance;
public:
	//opengl show graph
	double svg_width;
	double svg_height;
	string svg_view;

	void setInstance();
	void init();
	static void winreshape(GLint newWidth, GLint newHeight);
	static void mydisplay();
	static void myspecialkey(int key, int x, int y);
	static void paintCtlConvex();
	static void paintCtlPoint(vector<mycircle> circles);
	static void paintCubicBezierCurve(vector<cubicBezier> cbcurve);
	static void paintQuaBezierCurve(vector<quaBezier> qbcurve);
	static void paintPolyline(vector<PolyLine> pline);

	static void test_paintPCABezierCurve(vector<cubicBezier> cbcurve);

public:
		//void getPointList_BezierCurve();
};

