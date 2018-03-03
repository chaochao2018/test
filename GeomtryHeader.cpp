#include"GeomtryHeader.h"

double eps = 1e-5;
//double eps_isOnsegment = eps*100;//0.001  
double eps_isOnsegment = 0.00141421214815;  //����xy�����С������ȷ��0.001 ��������֮��ľ��������eps_isOnsegment->0.00141421214815�� 0.01��Ӧ0.0141421214815  

int dcmp(double x)
{
	if (fabs(x)<eps)
		return 0;
	else
		return x>0 ? 1 : -1;
}
double dot(mypoint2f a, mypoint2f b)
{
	return a.x*b.x + a.y*b.y;
}
double cross(mypoint2f a, mypoint2f b)
{
	return a.x*b.y - a.y*b.x;
}
double getEuclidDistance(vector<double>* vect1, vector<double>* vect2)
{
	double sum = 0;
	int size = vect1->size();
	for (int i = 0; i<size; i++)
	{
		sum += (vect1->at(i) - vect2->at(i))*(vect1->at(i) - vect2->at(i));
	}
	sum = sqrt(sum);
	return sum;
}
mypoint2f normalization(mypoint2f from_a, mypoint2f to_b)
{
	mypoint2f p = to_b - from_a;
	double vector_length = LineSegmentLength(to_b , from_a);
	p.x = p.x / vector_length;
	p.y = p.y / vector_length;
	return p;
}
double getCosineAngle(mypoint2f from_a1, mypoint2f to_a2, mypoint2f from_b1, mypoint2f to_b2)
{
	double aa = dot(to_a2 - from_a1, to_b2 - from_b1);
	double bb = (LineSegmentLength(from_a1, to_a2) * LineSegmentLength(from_b1, to_b2));
	double result = aa / bb;
	if (result > 1)
		result = 1;
	if (result<-1)
		result = -1;
	return result;
}
bool isSegmentIntersected(mypoint2f a1, mypoint2f a2, mypoint2f b1, mypoint2f b2)
{
	if (max(a1.x, a2.x) > min(b1.x, b2.x) && max(a1.y, a2.y) > min(b1.y, b2.y) && max(b1.x, b2.x) > min(a1.x, a2.x) && max(b1.y, b2.y) > min(a1.y, a2.y))
	{
		double c1 = cross(a2 - a1, b1 - a1);
		double c2 = cross(a2 - a1, b2 - a1);
		double c3 = cross(b2 - b1, a1 - b1);
		double c4 = cross(b2 - b1, a2 - b1);
		return dcmp(c1)*dcmp(c2) <= 0 && dcmp(c3)*dcmp(c4) <= 0;
	}
	else
		return false;
}
bool isOnSegment(mypoint2f p, mypoint2f a1, mypoint2f a2)
{
	double cross_r = cross(a1 - p, a2 - p);
	double dot_r = dot(a1 - p, a2 - p);
	if (fabs(cross_r) < (eps*10) && dcmp(dot_r) < 0)//���Ӧ�õ����㣬���ǿ��ǵ����������ཻ����ʵ��û�н���������
		return true;
	else
		return false;

	//return dcmp(cross(a1 - p, a2 - p)) == 0 && dcmp(dot(a1 - p, a2 - p))<0;

	/*if (p.x >= min(a1.x, a2.x) && p.x <= max(a1.x, a2.x) && p.y <= min(a1.y, a2.y) && p.y <= max(a1.y, a2.y))
	{
	return dcmp(cross(a1 - p, a2 - p)) == 0 && dcmp(dot(a1 - p, a2 - p))<0;
	}
	else
	return false;*/
}
double pointToLineDistance(mypoint2f p, mypoint2f lp1, mypoint2f lp2)
{
	double k, b;
	k = (lp1.y - lp2.y) / (lp1.x - lp2.x);
	b = lp1.y - k*lp1.x;
	double A=0, B=0, C=0;
	A = k;
	B = -1;
	C = b;
	double dis = abs(A*p.x + B*p.y + C) / sqrt((A*A) + (B*B));
	return dis;
}
double getAtan(mypoint2f p)
{
	double a = 0;
	if (p.x < eps)
		p.x = eps * 10;
	a = atan(p.y / p.x);
	return a;
}
double getAtan2(mypoint2f p)
{
	double a = 0;
	a = atan2(p.y, p.x);
	return a;
}
double getProximity(vector<mypoint2f>* plist1, vector<mypoint2f>* plist2, int &posi1, int &posi2)
{
	double minimum = 655350000;
	posi1 = 0; posi2 = 0;
	for (unsigned int i = 0; i < plist1->size(); ++i)
	{
		mypoint2f p1 = plist1->at(i);
		double last_len = 655350000;
		unsigned int j = 0;
		for ( j = 0; j < plist2->size(); ++j)
		{
			mypoint2f p2 = plist2->at(j);
			double tmp_len = LineSegmentLength(p1, p2);
			if (tmp_len < last_len)
			{
				last_len = tmp_len;
			}
			else
				break;
		}
		if (last_len < minimum)
		{
			minimum = last_len;
			posi1 = (int)i;
			posi2 = (int)j-1;
		}
	}
	return minimum;
}
double getContinuity(vector<mypoint2f>* plist1, vector<mypoint2f>* plist2, int posi1, int posi2)
{
	//���������������Ȼ��1-|dot(,)|  ֵԽСԽ���,[0��1]��
	/*mypoint2f tmp_p = plist1.front();
	plist1.insert(plist1.begin(), tmp_p);
	tmp_p = plist1.back();
	plist1.push_back(tmp_p);
	tmp_p = plist2.front();
	plist2.insert(plist2.begin(), tmp_p);
	tmp_p = plist2.back();
	plist2.push_back(tmp_p);
	posi1++; posi2++;
	mypoint2f vec1 = plist1[posi1] - plist1[posi1 - 1];
	mypoint2f vec2 = plist1[posi1+1] - plist1[posi1];
	mypoint2f vec3((vec1.x + vec2.x) / 2, (vec1.y + vec2.y) / 2);
	mypoint2f vec4 = plist2[posi2] - plist2[posi2 - 1];
	mypoint2f vec5 = plist2[posi2 + 1] - plist2[posi2];
	mypoint2f vec6((vec4.x + vec5.x) / 2, (vec4.y + vec5.y) / 2);
	double contn = 1 - abs(dot(vec3, vec6));
	return contn;*/

	mypoint2f vec1, vec2;
	if (posi1 != (int)plist1->size() - 1)
		vec1 = normalization(plist1->at(posi1), plist1->at(posi1 + 1));
	else
		vec1 = normalization(plist1->at(posi1-1), plist1->at(posi1));
	if (posi2 != (int)plist2->size() - 1)
		vec2 = normalization(plist2->at(posi2), plist2->at(posi2 + 1));
	else
		vec2 = normalization(plist2->at(posi2 - 1), plist2->at(posi2));
	//double contn = 1 - abs(dot(vec1, vec2)); //method 1
	double costheta = abs(getCosineAngle(mypoint2f(0, 0), vec1, mypoint2f(0, 0), vec2));   //method 2 ��|atan2|ֵ
	double contn = acos(costheta) * 180 / 3.1415926;   //ֵ��0-pi

	return contn;
}





double Area(const mypoint2f& p1, const mypoint2f& p2, const mypoint2f& p3)
{
	return cross(p1 - p2, p1 - p3);
}
double fArea(const mypoint2f& p1, const mypoint2f& p2, const mypoint2f& p3)
{
	return fabs(Area(p1, p2, p3));
}
mypoint2f getIntersectPoint(const mypoint2f& p1, const mypoint2f& p2, const mypoint2f& p3, const mypoint2f& p4)//�ж�����ֱ�ߵĽ��㣬�������߶Σ���
{
	//double k = fArea(p1, p2, p3) / fArea(p1, p2, p4);
	//return mypoint2f((p3.x + k*p4.x) / (1 + k), (p3.y + k*p4.y) / (1 + k));
	double s1 = fArea(p1, p2, p3);
	double s2 = fArea(p1, p2, p4);
	return mypoint2f((p4.x*s1 + p3.x*s2) / (s1 + s2), (p4.y*s1 + p3.y*s2) / (s1 + s2));
}

vector<double> quadricRealPolynomial(double coeff_a, double coeff_b, double coeff_c)
{
	vector<double> realroot;
	double a, b, c;
	double X1, X2;
	double delt = 0;
	a = coeff_a; b = coeff_b; c = coeff_c;
	if (a == 0)
	{
		X1 = -c / b;
		realroot.push_back(X1);	
	}
	else
	{
		delt = b*b - 4 * a*c;
		if (delt == 0)
		{
			X1 = -b / (2 * a);
			realroot.push_back(X1);
		}
		else if (delt > 0)
		{
			if (c == 0)
			{
				X1 = 0;
				X2 = -b / a;
				realroot.push_back(X1);
				realroot.push_back(X2);
			}
			else
			{
				X1 = (-b + sqrt(delt)) / (2 * a);
				X2 = (-b - sqrt(delt)) / (2 * a);
				realroot.push_back(X1);
				realroot.push_back(X2);
			}
		}
	}
	return realroot;
}
vector<double> cubicRealPolynomial(double coeff_a, double coeff_b, double coeff_c, double coeff_d)
{
	vector<double> realroot;
	double a, b, c, d;//����ϵ��  
	double A, B, C;
	double delt;//�б�ʽ  
	double X1, X2, X3;
	a = coeff_a; b = coeff_b; c = coeff_c; d = coeff_d;
	if (a == 0)
	{
		if (b == 0)
		{
			if (c != 0)
			{
				X1 = -d / c;
				realroot.push_back(X1);
			}
		}
		else
		{
			quadricRealPolynomial(coeff_b, coeff_c, coeff_d);
		}
	}
	else
	{
		A = pow(b, 2) - 3 * a*c;
		B = b*c - 9 * a*d;
		C = pow(c, 2) - 3 * b*d;
		delt = pow(B, 2) - 4 * A*C;
		//��A=B=0ʱ����ʽ1.+����6  
		if ((A == 0 && B == 0) || (delt == 0 && A == 0))
		{
			X1 = -b / (3 * a);//-c/b  or  -3d/c
			X2 = X1;
			X3 = X1;
			realroot.push_back(X1);
		}
		if (delt>0)//��=B2��4AC>0ʱ����ʽ2  
		{
			double Y1 = A*b + 3 * a*((-B + pow(delt, 0.5)) / 2);
			double Y2 = A*b + 3 * a*((-B - pow(delt, 0.5)) / 2);
			//�������������������  
			double Y1_three, Y2_three;
			if (Y1<0)
				Y1_three = -pow(-Y1, 1.0 / 3.0);
			else
				Y1_three = pow(Y1, 1.0 / 3.0);
			if (Y2<0)
				Y2_three = -pow(-Y2, 1.0 / 3.0);
			else
				Y2_three = pow(Y2, 1.0 / 3.0);
			X1 = (-b - (Y1_three + Y2_three)) / (3 * a);
			//X3,X2Ϊ���������ﲻ��Ҫ��δʵ�֡�  
			realroot.push_back(X1);
		}
		else if (delt == 0 && (A != 0))//����=B2��4AC=0ʱ����ʽ3 
		{
			double K = B / A;
			X1 = -b / a + K;
			X2 = -K / 2.0;
			X3 = X2;
			realroot.push_back(X1);
			realroot.push_back(X2);
		}
		else if (delt<0)//����=B2��4AC<0ʱ����ʽ4  
		{
			double T = (2 * A*b - 3 * a*B) / (2 * A*pow(A, 0.5));//��A>0����1<T<1��  
			double theita = acos(T);
			X1 = (-b - 2 * pow(A, 0.5)*cos(theita / 3.0)) / (3 * a);
			X2 = (-b + pow(A, 0.5)*(cos(theita / 3.0) + pow(3, 0.5)*sin(theita / 3.0))) / (3 * a);
			X3 = (-b + pow(A, 0.5)*(cos(theita / 3.0) - pow(3, 0.5)*sin(theita / 3.0))) / (3 * a);
			realroot.push_back(X1);
			realroot.push_back(X2);
			realroot.push_back(X3);
		}
	}
	return realroot;
}

double LineSegmentLength(mypoint2f p1, mypoint2f p2)
{
	double len = sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2));
	return len;
}
void CubicBezrCoefficient(double t, double &a1, double &a2, double &a3, double &a4)
{
	a1 = pow((1 - t), 3);
	a2 = pow((1 - t), 2) * 3 * t;
	a3 = 3 * t*t*(1 - t);
	a4 = t*t*t;
}
void CubicBezrDerivative(double t, double &b1, double &b2, double &b3, double &b4)
{
	b1 = -3 * pow((1 - t), 2);
	b2 = 3 * (3 * pow(t, 2) - 4 * t + 1);
	b3 = 3 * (2 * t - 3 * pow(t, 2));
	b4 = 3 * pow(t, 2);
}
void QuaBezrCoefficient(double t, double &a1, double &a2, double &a3)
{
	a1 = pow((1 - t), 2);
	a2 = 2 * t*(1 - t);
	a3 = pow(t, 2);
}
void QuaBezrDerivative(double t, double &b1, double &b2, double &b3)
{
	b1 = -2 * (1 - t);
	b2 = 2 - 4 * t;
	b3 = 2 * t;
}


double GetPolygonArea(vector<mypoint2f>* plist)//plist������β��ͬ
{
	////����һ�����ĵ��ǣ�0��0����patch����svg����ϵͳ������ʱ��ģ�ת����opengl��˳ʱ���,�����svg����ϵ�£�����Ǹ��ģ���opengl����ϵ�ֶ�����Ǹ��ģ�
	//float view_y = 1052;
	//for (unsigned int i = 0; i < plist.size(); ++i)//svg---->opengl������ת��
	//{
	//	plist[i].y = view_y - plist[i].y;
	//}

	double area = 0;
	double temp = 0;
	for (unsigned int i = 1; i < plist->size(); ++i)
	{
		temp = cross(plist->at(i-1), plist->at(i));
		area = area + temp;
	}
	//��svg����ϵ�£�area����Ǹ���ֵ,����-1��Ϊ��ֵ��ע��������area����������������ģ�˵������overlap patch���� repeatpatch
	area = area / 2*(-1);
	return area;

	//mypoint2f startp = plist.front();
	//plist.push_back(startp);
	//double area = 0;
	//double temp = 0;
	//for (int i = plist.size()-1; i >=1; --i)
	//{
	//	temp = cross(plist[i], plist[i-1]);
	//	area = area + temp;
	//}
	//area = area / 2;
	//return area;

	//������ patch����ʱ��ģ��õ���patch���<0,���������˳��Ϊ+
	//mypoint2f startp = plist.front();
	//double area = 0;
	//double temp = 0;
	//for (unsigned int i = 1; i < plist.size()-1; ++i)
	//{
	//	temp = cross(plist[i]-startp, plist[i+1]-startp);
	//	area = area + temp;
	//}
	//area = area / 2;
	//return area;

	//mypoint2f startp = plist.front();
	//double area = 0;
	//double temp = 0;
	//for (int i = plist.size() - 1; i >= 1; --i)
	//{
	//	temp = cross(plist[i]-startp, plist[i-1]-startp);
	//	area = area + temp;
	//}
	//area = area / 2;
	//return area;

}

int dbscanclustering(vector<DataPoint> &dataset, double radius, int minpts)
{
	int datasize = dataset.size();
	for (int i = 0; i < datasize; ++i)
	{
		vector<double> ft1 = dataset[i].pch_ft;
		for (int j = i+1; j < datasize; ++j)
		{
			vector<double> ft2 = dataset[j].pch_ft;
			double dis = getEuclidDistance(&ft1, &ft2);
			if (dis <= radius)
			{
				dataset[i].neigbor_ct.push_back(j);//dataset[j].dataID=j
				dataset[j].neigbor_ct.push_back(i);
			}
		}
	}
	for (int i = 0; i < datasize; ++i)
	{
		if (dataset[i].neigbor_ct.size() >= minpts)
		{
			dataset[i].datatype = 3;
		}
	}

	int clusterindex = 1;
	for (int i = 0; i < datasize; ++i)
	{
		if (!dataset[i].visited && dataset[i].datatype == 3)
		{
			stack<int> core_stack;
			core_stack.push(i);
			dataset[i].clusterId = clusterindex;
			while (!core_stack.empty())
			{
				int tmp = core_stack.top();
				core_stack.pop();		
				dataset[i].visited = true;
				for (unsigned int j = 0; j < dataset[tmp].neigbor_ct.size(); ++j)
				{
					int pindex = dataset[tmp].neigbor_ct[j];
					if (dataset[pindex].datatype == 3)
					{
						if (!dataset[pindex].visited)
						{
							core_stack.push(pindex);
							dataset[pindex].visited = true;//!!
							dataset[pindex].clusterId = clusterindex;
						}
					}
					else
					{
						if (!dataset[pindex].visited)
						{
							dataset[pindex].visited = true;//!!
							dataset[pindex].datatype = 2;//�߽��
							dataset[pindex].clusterId = clusterindex;
						}
					}
					//if (dataset[pindex].datatype == 3 && !dataset[pindex].visited)
					//{
					//	core_stack.push(pindex);
					//	dataset[pindex].visited = true;//!!
					//	dataset[pindex].clusterId = clusterindex;
					//}
				}
			}
			clusterindex++;
		}
	}
	return clusterindex;
}
//void kmeans(vector<pair<int, double>> &type_dis, int classnum, vector<pair<int, vector<double>>> pchFeathers)
//{
//	int ft_dimention = pchFeathers[0].second.size();//ft_dimension ������������ά��
//	double maxlong = -1; double minlong = 1; double maxArea = 0;
//
//	//����patch index������������������������
//	vector<vector<double>> centroids;
//	srand((unsigned)time(NULL));
//	for (int i = 0; i < classnum; i++)
//	{
//		int rnum = rand() % (patchfeatures.size() - 1);
//		cout << rnum << ' ';
//		centroids.push_back(patchfeatures[rnum].second);
//	}
//}


double GetAverange(vector<double> vec)
{
	int size = vec.size();
	if (size == 0)
		return 0;
	else
	{
		double sum = 0;
		for (int i = 0; i < size; ++i)
		{
			sum = vec[i] + sum;
		}
		return sum / size;
	}
}
double GetVariance(vector<double> vec1, double averange1, vector<double> vec2, double averange2)
{
	double aver1 = 0; double aver2 = 0;
	if (averange1 >= 0)
		aver1 = averange1;
	else
		aver1 = GetAverange(vec1);
	if (averange2 >= 0)
		aver2 = averange2;
	else
		aver2 = GetAverange(vec2);

	int size = vec1.size();
	double sum = 0;
	for (int i = 0; i < size; ++i)
	{
		sum = (vec1[i] - aver1)*(vec2[i] - aver2) + sum;
	}
	return sum / (size-1);
}
double GetVariance_xx(vector<double> vec1, double averange1)
{
	double aver1 = 0; double aver2 = 0;
	if (averange1 >= 0)
		aver1 = averange1;
	else
		aver1 = GetAverange(vec1);
	int size = vec1.size();
	double sum = 0;
	for (int i = 0; i < size; ++i)
	{
		sum = pow((vec1[i] - aver1),2) + sum;
	}
	return sum / (size - 1);
}

bool mycompare(const pair<int, double>& v1, const pair<int, double>& v2)
{
	return v1.second<v2.second; //��С��������
}



bool checkAcuteTriangle(mypoint2f p1, mypoint2f p2, mypoint2f p3)
{
	bool flag = false;
	//���ж��ǲ���������
	double a = LineSegmentLength(p1, p2);//�������߶εĳ���
	double b = LineSegmentLength(p2, p3);
	double c = LineSegmentLength(p1, p3);
	if ((a + b)>c && (a + c)>b && (b + c) > a)
	{
		//���ܹ��������Σ� ���ж��ǲ������������
		double an1 = getCosineAngle(p1, p2, p1, p3);
		double an2 = getCosineAngle(p2, p3, p2, p1);
		double an3 = getCosineAngle(p3, p1, p3, p2);
		if (an1 < 0 || an2<0 || an3<0)
		{
			flag = false;
		}
		else
			flag = true;
	}
	return flag;
}
void makeAcuteTrangle(vector<mypoint2f>* plist, mypoint2f *p1, mypoint2f *p2, mypoint2f *p3)
{
	////����һ����pointlistƽ���ֳ�3�ݣ�
	//int divide = plist->size() / 3;
	//*p1 = plist->at(0);
	//*p2 = plist->at(divide);
	//*p3 = plist->at(divide * 2);
	////�ж��ǲ������������ ,������
	//double an1 = getCosineAngle(*p1, *p2, *p1, *p3);
	//double an2 = getCosineAngle(*p2, *p3, *p2, *p1);
	//double an3 = getCosineAngle(*p3, *p1, *p3, *p2);
	//int begin_i = 0;
	//int end_i = 0;
	//if (an1 >= -1 && an1 <= 0)
	//{
	//	//p1���Ƕ۽ǣ��̶�p1p2,��p3�ƶ�
	//	double tmp_an = getCosineAngle(*p1, *p2, *p1, plist->at(divide * 2 + 1));
	//	if (tmp_an > an1)
	//	{
	//		bool change = false;
	//		begin_i = divide * 2 + 1;
	//		end_i = plist->size() - 1;
	//		for (begin_i; begin_i <= end_i; ++begin_i)
	//		{
	//			tmp_an = getCosineAngle(*p1, *p2, *p1, plist->at(begin_i));
	//			if (tmp_an > 0 && tmp_an < 1)
	//			{
	//				*p3 = plist->at(begin_i);
	//				change = true;
	//				break;
	//			}
	//		}
	//		if (change == false)//�̶�p1p3,��p2�ƶ�
	//		{
	//			tmp_an = getCosineAngle(*p1, *p3, *p1, plist->at(divide + 1));
	//			if (tmp_an > an2)
	//			{
	//				begin_i = divide * 2 + 1;
	//				end_i = plist->size() - 1;
	//				for (begin_i; begin_i <= end_i; ++begin_i)
	//				{
	//					tmp_an = getCosineAngle(*p1, *p2, *p1, plist->at(begin_i));
	//					if (tmp_an > 0 && tmp_an < 1)
	//					{
	//						*p3 = plist->at(begin_i);
	//						change = true;
	//						break;
	//					}
	//				}
	//			}
	//			
	//		}
	//	}
	//	else
	//	{
	//		begin_i = divide * 2 - 1;
	//		end_i = divide + 1;
	//		for (begin_i; begin_i >= end_i; --begin_i)
	//		{
	//			tmp_an = getCosineAngle(*p1, *p2, *p1, plist->at(begin_i));
	//			if (tmp_an > 0 && tmp_an < 1)
	//			{
	//				*p3 = plist->at(begin_i);
	//				break;
	//			}
	//		}
	//	}
	//}
	//else
	//{//p2���Ƕ۽ǣ��̶�p2p1,��p3�ƶ�
	//	if (an2 >= -1 && an2 <= 0)
	//	{
	//		double tmp_an = getCosineAngle(*p2, *p1, *p2, plist->at(divide * 2 + 1));
	//		if (tmp_an > an2)
	//		{
	//			begin_i = divide * 2 + 1;
	//			end_i = plist->size() - 1;
	//			for (begin_i; begin_i <= end_i; ++begin_i)
	//			{
	//				tmp_an = getCosineAngle(*p2, *p1, *p2, plist->at(begin_i));
	//				if (tmp_an > 0 && tmp_an < 1)
	//				{
	//					*p3 = plist->at(begin_i);
	//					break;
	//				}
	//			}
	//		}
	//		else
	//		{
	//			begin_i = divide * 2 - 1;
	//			end_i = divide + 1;
	//			for (begin_i; begin_i >= end_i; --begin_i)
	//			{
	//				tmp_an = getCosineAngle(*p2, *p1, *p2, plist->at(begin_i));
	//				if (tmp_an > 0 && tmp_an < 1)
	//				{
	//					*p3 = plist->at(begin_i);
	//					break;
	//				}
	//			}
	//		}
	//	}
	//	else
	//	{
	//		if (an3 >= -1 && an3 <= 0)
	//		{//p3���Ƕ۽ǣ��̶�p3p1,��p2�ƶ�
	//			double tmp_an = getCosineAngle(*p3, *p1, *p3, plist->at(divide + 1));
	//			if (tmp_an > an3)
	//			{
	//				begin_i = divide + 1;
	//				end_i = divide * 2 - 1;
	//				for (begin_i; begin_i <= end_i; ++begin_i)
	//				{
	//					tmp_an = getCosineAngle(*p3, *p1, *p3, plist->at(begin_i));
	//					if (tmp_an > 0 && tmp_an < 1)
	//					{
	//						*p2 = plist->at(begin_i);
	//						break;
	//					}
	//				}
	//			}
	//			else
	//			{
	//				begin_i = divide - 1;
	//				end_i = 1;
	//				for (begin_i; begin_i >= end_i; --begin_i)
	//				{
	//					tmp_an = getCosineAngle(*p3, *p1, *p3, plist->at(begin_i));
	//					if (tmp_an > 0 && tmp_an < 1)
	//					{
	//						*p2 = plist->at(begin_i);
	//						break;
	//					}
	//				}
	//			}
	//		}
	//	}
	//}

	//ƽ���ֳ�4��	
//int divide = plist->size() / 4;
	//*p1 = plist->at(divide);
	//*p2 = plist->at(divide*2);
	//*p3 = plist->at(divide * 3);
	////�ж��ǲ������������ ,������
	//double an1 = getCosineAngle(*p1, *p2, *p1, *p3);
	//double an2 = getCosineAngle(*p2, *p3, *p2, *p1);
	//double an3 = getCosineAngle(*p3, *p1, *p3, *p2);
	//int begin_i = 0;
	//int end_i = 0;
	//if (an1 >= -1 && an1 <= 0)
	//{
	//	//p1���Ƕ۽ǣ��̶�p1p2,��p3�ƶ�
	//	double tmp_an = getCosineAngle(*p1, *p2, *p1, plist->at(divide * 3 + 1));
	//	if (tmp_an > an1)
	//	{
	//		begin_i = divide * 3 + 1;
	//		end_i = plist->size() - 1;
	//		for (begin_i; begin_i <= end_i; ++begin_i)
	//		{
	//			tmp_an = getCosineAngle(*p1, *p2, *p1, plist->at(begin_i));
	//			if (tmp_an > 0 && tmp_an < 1)
	//			{
	//				*p3 = plist->at(begin_i);
	//				break;
	//			}
	//		}
	//	}
	//	else
	//	{
	//		begin_i = divide * 3 - 1;
	//		end_i = divide*2 + 1;
	//		for (begin_i; begin_i >= end_i; --begin_i)
	//		{
	//			tmp_an = getCosineAngle(*p1, *p2, *p1, plist->at(begin_i));
	//			if (tmp_an > 0 && tmp_an < 1)
	//			{
	//				*p3 = plist->at(begin_i);
	//				break;
	//			}
	//		}
	//	}
	//}
	//else
	//{//p2���Ƕ۽ǣ��̶�p2p1,��p3�ƶ�
	//	if (an2 >= -1 && an2 <= 0)
	//	{
	//		double tmp_an = getCosineAngle(*p2, *p1, *p2, plist->at(divide * 3 + 1));
	//		if (tmp_an > an2)
	//		{
	//			begin_i = divide * 3 + 1;
	//			end_i = plist->size() - 1;
	//			for (begin_i; begin_i <= end_i; ++begin_i)
	//			{
	//				tmp_an = getCosineAngle(*p2, *p1, *p2, plist->at(begin_i));
	//				if (tmp_an > 0 && tmp_an < 1)
	//				{
	//					*p3 = plist->at(begin_i);
	//					break;
	//				}
	//			}
	//		}
	//		else
	//		{
	//			begin_i = divide * 3 - 1;
	//			end_i = divide*2 + 1;
	//			for (begin_i; begin_i >= end_i; --begin_i)
	//			{
	//				tmp_an = getCosineAngle(*p2, *p1, *p2, plist->at(begin_i));
	//				if (tmp_an > 0 && tmp_an < 1)
	//				{
	//					*p3 = plist->at(begin_i);
	//					break;
	//				}
	//			}
	//		}
	//	}
	//	else
	//	{
	//		if (an3 >= -1 && an3 <= 0)
	//		{//p3���Ƕ۽ǣ��̶�p3p1,��p2�ƶ�
	//			double tmp_an = getCosineAngle(*p3, *p1, *p3, plist->at(divide*2 + 1));
	//			if (tmp_an > an3)
	//			{
	//				begin_i = divide*2 + 1;
	//				end_i = divide * 3 - 1;
	//				for (begin_i; begin_i <= end_i; ++begin_i)
	//				{
	//					tmp_an = getCosineAngle(*p3, *p1, *p3, plist->at(begin_i));
	//					if (tmp_an > 0 && tmp_an < 1)
	//					{
	//						*p2 = plist->at(begin_i);
	//						break;
	//					}
	//				}
	//			}
	//			else
	//			{
	//				begin_i = divide*2 - 1;
	//				end_i = divide+1;
	//				for (begin_i; begin_i >= end_i; --begin_i)
	//				{
	//					tmp_an = getCosineAngle(*p3, *p1, *p3, plist->at(begin_i));
	//					if (tmp_an > 0 && tmp_an < 1)
	//					{
	//						*p2 = plist->at(begin_i);
	//						break;
	//					}
	//				}
	//			}
	//		}
	//	}
	//}

	//����2 ���ҵ��յ㣬��ȡ���յ���е㣬������һ��patch�����������յ�
	vector<mypoint2f> tmpplist =*plist;//plist�Ǳ�ԭ����һ���ĵ㣨���ܼ���
	mypoint2f pfront = plist->front();
	mypoint2f pback = plist->back();
	tmpplist.insert(tmpplist.begin(), pback);
	tmpplist.push_back(pfront);
	int turn_count = 0;
	vector<int> record_tp;
	for (unsigned int j = 1; j < tmpplist.size() - 1; ++j)
	{
		mypoint2f vec1 = tmpplist.at(j) - tmpplist.at(j - 1);
		mypoint2f vec2 = tmpplist.at(j + 1) - tmpplist.at(j);
		double cos_theta = getCosineAngle(tmpplist.at(j - 1), tmpplist.at(j), tmpplist.at(j), tmpplist.at(j + 1));
		if (cos_theta >= -1 && cos_theta <= 0.707106)//ת����45��180��Ϊ��һ���յ�   135-->-0.7071  150--> -0.866
		{
			turn_count++; //cout << cos_theta << endl;
			record_tp.push_back(j-1);
		}
	}
	tmpplist.pop_back();
	tmpplist.erase(tmpplist.begin());
	if (turn_count > 1)//patch ����Բ��
	{
		if (record_tp[0] != 0)
		{
			int posi = record_tp[0];
			for (unsigned int i = 0; i < record_tp.size(); ++i)
			{
				record_tp[i] = record_tp[i] - posi;
			}
			vector<mypoint2f>::iterator iter = tmpplist.begin();
			while (posi > 0)
			{
				mypoint2f pp = tmpplist.front();
				tmpplist.push_back(pp);
				iter = tmpplist.erase(iter);
				posi--;
			}
		}
		if (turn_count == 2)
		{
			int mid_posi = (record_tp[1] - record_tp[0]) / 2;
			*p1 = tmpplist[mid_posi];
			mid_posi = record_tp[1] + (tmpplist.size() - record_tp[1]) / 2;
			*p2 = tmpplist[mid_posi];
			bool findflag = false;
			for (unsigned int i = 0; i < tmpplist.size(); ++i)
			{
				if (tmpplist[i] != *p1 && tmpplist[i] != *p2)
				{
					mypoint2f tmp_p3 = tmpplist[i];
					bool acu_tri = checkAcuteTriangle(*p1, *p2, tmp_p3);
					if (acu_tri == true)
					{
						*p3 = tmp_p3;
						findflag = true;
						break;
					}
				}
			}
			if (findflag == false)
			{
				//û���ҵ���Ӧ��p3,
				*p3 = mypoint2f(-1, -1);
			}
		}
		else if (turn_count == 3)
		{
			//ѡ�ϳ���������
			int mid_posi = (record_tp[1] - record_tp[0]) / 2;
			mypoint2f tmp_p1 = tmpplist[mid_posi];
			mid_posi = record_tp[1] + (record_tp[2] - record_tp[1]) / 2;
			mypoint2f tmp_p2 = tmpplist[mid_posi];
			mid_posi = record_tp[2] + (tmpplist.size() - record_tp[2]) / 2;
			mypoint2f tmp_p3 = tmpplist[mid_posi];
			double len1 = LineSegmentLength(tmp_p1, tmp_p2);
			double len2 = LineSegmentLength(tmp_p2, tmp_p3);
			double len3 = LineSegmentLength(tmp_p1, tmp_p3);
			if (len1 < len2 && len1 < len3)
			{
				*p1 = tmp_p1;
				*p2 = tmp_p2;
			}
			else if (len2 < len1 && len2 < len3)
			{
				*p1 = tmp_p2;
				*p2 = tmp_p3;
			}
			else if (len3 < len1 &&len3 < len2)
			{
				*p1 = tmp_p1;
				*p2 = tmp_p3;
			}
			bool findflag = false;
			for (unsigned int i = 0; i < tmpplist.size(); ++i)
			{
				if (tmpplist[i] != *p1 && tmpplist[i] != *p2)
				{
					tmp_p3 = tmpplist[i];
					bool acu_tri = checkAcuteTriangle(*p1, *p2, tmp_p3);
					if (acu_tri == true)
					{
						*p3 = tmp_p3;
						findflag = true;
						break;
					}
				}
			}
			if (findflag == false)
			{
				//û���ҵ���Ӧ��p3,
				*p3 = mypoint2f(-1, -1);
			}
		}
		else
		{
			////����һ��ƽ���ֳ�3��
			//int divide = plist->size() / 3;
			//*p1 = plist->at(0);
			//*p2 = plist->at(divide);
			//*p3 = plist->at(divide * 2);
			////�ж��ǲ������������ ,������
			//double an1 = getCosineAngle(*p1, *p2, *p1, *p3);
			//double an2 = getCosineAngle(*p2, *p3, *p2, *p1);
			//double an3 = getCosineAngle(*p3, *p1, *p3, *p2);
			//int begin_i = 0;
			//int end_i = 0;
			//if ( an1 <= 0)
			//{
			//	//p1���Ƕ۽ǣ��̶�p1p2,��p3�ƶ�
			//	double tmp_an = getCosineAngle(*p1, *p2, *p1, plist->at(divide * 2 + 1));
			//	bool change = false;
			//	if (tmp_an > an1)
			//	{
			//		bool change = false;
			//		begin_i = divide * 2 + 1;
			//		end_i = plist->size() - 1;
			//		for (begin_i; begin_i <= end_i; ++begin_i)
			//		{
			//			tmp_an = getCosineAngle(*p1, *p2, *p1, plist->at(begin_i));
			//			if (tmp_an > 0 && tmp_an < 1)
			//			{
			//				*p3 = plist->at(begin_i);
			//				change = true;
			//				break;
			//			}
			//		}
			//		if (change == false)//�̶�p1p3,��p2�ƶ�
			//		{
			//			tmp_an = getCosineAngle(*p1, *p3, *p1, plist->at(divide + 1));
			//			if (tmp_an > an2)
			//			{
			//				begin_i = divide * 2 + 1;
			//				end_i = plist->size() - 1;
			//				for (begin_i; begin_i <= end_i; ++begin_i)
			//				{
			//					tmp_an = getCosineAngle(*p1, *p2, *p1, plist->at(begin_i));
			//					if (tmp_an > 0 && tmp_an < 1)
			//					{
			//						*p3 = plist->at(begin_i);
			//						change = true;
			//						break;
			//					}
			//				}
			//			}
			//			if (change == false)
			//			{
			//				*p3 = mypoint2f(-1, -1);
			//			}
			//		}
			//	}
			//	else
			//	{
			//		bool cflag = false;
			//		begin_i = divide * 2 - 1;
			//		end_i = divide + 1;
			//		for (begin_i; begin_i >= end_i; --begin_i)
			//		{
			//			tmp_an = getCosineAngle(*p1, *p2, *p1, plist->at(begin_i));
			//			if (tmp_an > 0 && tmp_an < 1)
			//			{
			//				*p3 = plist->at(begin_i);
			//				cflag = true;
			//				break;
			//			}
			//		}
			//		if (cflag==false)
			//			*p3 = mypoint2f(-1, -1);
			//	}
			//}
			//else
			//{//��p2���Ƕ۽ǣ��̶�p2p1,��p3�ƶ�
			//	if (an2 >= -1 && an2 <= 0)
			//	{
			//		double tmp_an = getCosineAngle(*p2, *p1, *p2, plist->at(divide * 2 + 1));
			//		bool cflag = false;
			//		if (tmp_an > an2)
			//		{
			//			begin_i = divide * 2 + 1;
			//			end_i = plist->size() - 1;
			//			for (begin_i; begin_i <= end_i; ++begin_i)
			//			{
			//				tmp_an = getCosineAngle(*p2, *p1, *p2, plist->at(begin_i));
			//				if (tmp_an > 0 && tmp_an < 1)
			//				{
			//					*p3 = plist->at(begin_i);
			//					cflag = true;
			//					break;
			//				}
			//			}
			//		}
			//		else
			//		{
			//			begin_i = divide * 2 - 1;
			//			end_i = divide + 1;
			//			for (begin_i; begin_i >= end_i; --begin_i)
			//			{
			//				tmp_an = getCosineAngle(*p2, *p1, *p2, plist->at(begin_i));
			//				if (tmp_an > 0 && tmp_an < 1)
			//				{
			//					*p3 = plist->at(begin_i);
			//					cflag = true;
			//					break;
			//				}
			//			}
			//		}
			//		if (cflag ==false)
			//			*p3 = mypoint2f(-1, -1);
			//	}
			//	else
			//	{
			//		if (an3 >= -1 && an3 <= 0)
			//		{//p3���Ƕ۽ǣ��̶�p3p1,��p2�ƶ�
			//			double tmp_an = getCosineAngle(*p3, *p1, *p3, plist->at(divide + 1));
			//			bool cflag = false;
			//			if (tmp_an > an3)
			//			{
			//				begin_i = divide + 1;
			//				end_i = divide * 2 - 1;
			//				for (begin_i; begin_i <= end_i; ++begin_i)
			//				{
			//					tmp_an = getCosineAngle(*p3, *p1, *p3, plist->at(begin_i));
			//					if (tmp_an > 0 && tmp_an < 1)
			//					{
			//						*p2 = plist->at(begin_i);
			//						cflag = true;
			//						break;
			//					}
			//				}
			//			}
			//			else
			//			{
			//				begin_i = divide - 1;
			//				end_i = 1;
			//				for (begin_i; begin_i >= end_i; --begin_i)
			//				{
			//					tmp_an = getCosineAngle(*p3, *p1, *p3, plist->at(begin_i));
			//					if (tmp_an > 0 && tmp_an < 1)
			//					{
			//						*p2 = plist->at(begin_i);
			//						cflag = true;
			//						break;
			//					}
			//				}
			//			}
			//			if (cflag == false)
			//				*p3 = mypoint2f(-1, -1);
			//		}
			//	}
			//}
			//������
			int mid_posi1 = (record_tp[1] - record_tp[0]) / 2;
			*p1 = tmpplist[mid_posi1];
			int mid_posi2 = record_tp[1] + (record_tp[2] - record_tp[1]) / 2;
			*p2 = tmpplist[mid_posi2];
			int mid_posi3 = record_tp[2] + (tmpplist.size() - record_tp[2]) / 2;
			*p3 = tmpplist[mid_posi3];
			double an1 = getCosineAngle(*p1, *p2, *p1, *p3);
			double an2 = getCosineAngle(*p2, *p3, *p2, *p1);
			double an3 = getCosineAngle(*p3, *p1, *p3, *p2);
			if (an1 < 0 )//p1���Ƕ۽�
			{
				bool change = false;
				for (unsigned int i = 0; i<plist->size(); ++i)//�̶�p1p2 �ƶ�p3
				{
					if (plist->at(i) != *p1 &&plist->at(i) != *p2)
					{
						double tmp_an = getCosineAngle(*p1, *p2, *p1, plist->at(i));
						if (tmp_an > 0 && tmp_an < 1)
						{
							*p3 = plist->at(i);
							change = true;
							break;
						}
					}
				}
				if (change == false)//�̶�p1p3 �ƶ�p2
				{
					for (unsigned int i = 0; i<plist->size(); ++i)
					{
						if (plist->at(i) != *p1 &&plist->at(i) != *p3)
						{
							double tmp_an = getCosineAngle(*p1, *p3, *p1, plist->at(i));
							if (tmp_an > 0 && tmp_an < 1)
							{
								*p2 = plist->at(i);
								change = true;
								break;
							}
						}
					}
				}
				if (change == false)
					*p3 = mypoint2f(-1, -1);
			}
			else
			{
				//��p2���Ƕ۽�
				if (an2 <0)
				{
					bool change = false;
					for (unsigned int i = 0; i<plist->size(); ++i)//�̶�p2p1,��p3�ƶ�
					{
						if (plist->at(i) != *p1 &&plist->at(i) != *p2)
						{
							double tmp_an = getCosineAngle(*p2, *p1, *p2, plist->at(i));
							if (tmp_an > 0 && tmp_an < 1)
							{
								*p3 = plist->at(i);
								change = true;
								break;
							}
						}
					}
					if (change == false)//�̶�p2p3,��p1�ƶ�
					{
						for (unsigned int i = 0; i<plist->size(); ++i)
						{
							if (plist->at(i) != *p2 &&plist->at(i) != *p3)
							{
								double tmp_an = getCosineAngle(*p2, *p3, *p2, plist->at(i));
								if (tmp_an > 0 && tmp_an < 1)
								{
									*p1 = plist->at(i);
									change = true;
									break;
								}
							}
						}
					}
					if (change == false)
						*p3 = mypoint2f(-1, -1);
				}
				else
				{
					//��p3���Ƕ۽�
					if (an3 < 0)
					{
						bool change = false;
						for (unsigned int i = 0; i<plist->size(); ++i)//�̶�p3p1,��p2�ƶ�
						{
							if (plist->at(i) != *p1 &&plist->at(i) != *p3)
							{
								double tmp_an = getCosineAngle(*p3, *p1, *p3, plist->at(i));
								if (tmp_an > 0 && tmp_an < 1)
								{
									*p2 = plist->at(i);
									change = true;
									break;
								}
							}
						}
						if (change == false)//�̶�p3p2,��p1�ƶ�
						{
							for (unsigned int i = 0; i<plist->size(); ++i)
							{
								if (plist->at(i) != *p3 &&plist->at(i) != *p2)
								{
									double tmp_an = getCosineAngle(*p3, *p2, *p3, plist->at(i));
									if (tmp_an > 0 && tmp_an < 1)
									{
										*p1 = plist->at(i);
										change = true;
										break;
									}
								}
							}
						}
						if (change == false)
							*p3 = mypoint2f(-1, -1);
					}
				}
			}
		}
	}
	else//Բ �����
	{
		int divide = plist->size() / 3;
		*p1 = plist->at(0);
		*p2 = plist->at(divide);
		*p3 = plist->at(divide * 2);
		bool acute = checkAcuteTriangle(*p1, *p2,*p3);
		if (acute==false)
			*p3 = mypoint2f(-1, -1);
	}
}
void maxInscribedCircle(vector<mypoint2f> *polygon, mypoint2f &center, double &radius)
{
	mypoint2f p1, p2, p3;
	if (polygon->size()<=10)//polygon��������٣����Ƽ������ڽ�Բ�뾶��Բ��
	{
		mypoint2f tmp(0, 0);
		for (unsigned int i = 0; i < polygon->size(); ++i)
		{
			tmp = tmp + polygon->at(i);
		}
		tmp.x = tmp.x / polygon->size();
		tmp.y = tmp.y / polygon->size();
		double mini_r = 65535;
		for (unsigned int i = 0; i < polygon->size(); ++i)
		{
			double dis = LineSegmentLength(polygon->at(i), tmp);
			if (dis < mini_r)
			{
				mini_r = dis;
				center = polygon->at(i);
			}
		}
		radius = mini_r;
	}
	else
	{
		//����һЩ�㣬����
		vector<mypoint2f>::iterator iter_pl = polygon->begin()+1;
		while (iter_pl != polygon->end())
		{
			mypoint2f spt = *iter_pl;
			mypoint2f ept = *(iter_pl - 1);
			mypoint2f mpt((spt.x + ept.x)/2, (spt.y + ept.y) / 2);
			iter_pl = polygon->insert(iter_pl,mpt);
			iter_pl = iter_pl + 2;//+2?
		}
		makeAcuteTrangle(polygon, &p1, &p2, &p3);//�õ���ʼp1p2p3�� ����������Σ�
		if (p3 == mypoint2f(-1, -1))//û���ҵ����ʵ���������Σ���if�еķ���
		{
			mypoint2f tmp(0, 0);
			for (unsigned int i = 0; i < polygon->size(); ++i)
			{
				tmp = tmp + polygon->at(i);
			}
			tmp.x = tmp.x / polygon->size();
			tmp.y = tmp.y / polygon->size();
			double mini_r = 65535;
			for (unsigned int i = 0; i < polygon->size(); ++i)
			{
				double dis = LineSegmentLength(polygon->at(i), tmp);
				if (dis < mini_r)
				{
					mini_r = dis;
					center = polygon->at(i);
				}
			}
			radius = mini_r;
		}
		else
		{
			double mini_dis = 65535;
			double tmp_r = 1; //�뾶
			mypoint2f tmp_center(0, 0);//Բ��
			double cha = mini_dis - tmp_r;
			while (cha > 0.01)//ֱ���뾶Լ������С����
			{
				//��֪p1p2p3,����Բ���� 
				//tmp_center = mypoint2f(0, 0);
				//mypoint2f pm((p1.x + p3.x) / 2, (p1.y + p3.y) / 2);
				//mypoint2f pn((p2.x + p3.x) / 2, (p2.y + p3.y) / 2);
				//double fenmu = (p3.x - p1.x)*(p3.y - p2.y) - (p3.x - p2.x)*(p3.y - p1.y);
				//tmp_center.x = (pm.y - pn.y)*(p3.y - p1.y)*(p3.y - p2.y) + pm.x*(p3.x - p1.x)*(p3.y - p2.y) - pn.x*(p3.x - p2.x)*(p3.y - p1.y);
				//tmp_center.x = tmp_center.x / fenmu;
				//tmp_center.y = (pn.x - pm.x)*(p3.x - p1.x)*(p3.x - p2.x) - pm.y*(p3.y - p1.y)*(p3.x - p2.x) + pn.y*(p3.y - p2.y)*(p3.x - p1.x);
				//tmp_center.y = tmp_center.y / fenmu;	
				double u1 = (pow(p2.x, 2) - pow(p1.x, 2) + pow(p2.y, 2) - pow(p1.y, 2)) / 2;
				double u2 = (pow(p3.x, 2) - pow(p1.x, 2) + pow(p3.y, 2) - pow(p1.y, 2)) / 2;
				double d11 = p2.x - p1.x; double d12 = p2.y - p1.y; double d21 = p3.x - p1.x; double d22 = p3.y - p1.y;
				tmp_center.x = (u1*d22 - u2*d12) / (d11*d22 - d21*d12);
				tmp_center.y = (u2*d11 - u1*d21) / (d11*d22 - d21*d12);
				//�ж��ǲ���������
				double last_rmp = 0;
				for (unsigned int i = 1; i < polygon->size(); ++i)
				{
					mypoint2f vec1 = polygon->at(i) - polygon->at(i - 1);
					mypoint2f vec2 = tmp_center - polygon->at(i - 1);
					double bbb = cross(vec1, vec2);///�� P �� Q > 0 , ��P��Q��˳ʱ�뷽�� =0����
					if (last_rmp*bbb < 0)
					{
						//cout << "center point��patch ����" << endl;
						break;
					}
					else
						last_rmp = bbb;
				}

				//����뾶
				tmp_r = LineSegmentLength(p1, tmp_center);
				//���㵱ǰԲ�ĵ�polygon���������С����
				mypoint2f newpoint(0, 0);
				mini_dis = 65535;
				for (unsigned int i = 0; i < polygon->size(); ++i)
				{
					double dis = LineSegmentLength(polygon->at(i), tmp_center);
					if (dis < mini_dis)
					{
						mini_dis = dis;
						newpoint = polygon->at(i);
					}
				}
				//������С�����ǲ��ǵ��ڰ뾶
				cha = fabs(tmp_r - mini_dis);
				//����С�����Ӧ�Ķ������p1p2p3����һ�������·���p1p2p3
				//1:newpoin1+p2 p3
				double angle1 = getCosineAngle(newpoint, p2, newpoint, p3);
				double angle2 = getCosineAngle(p2, newpoint, p2, p3);
				double angle3 = getCosineAngle(p3, newpoint, p3, p1);
				if ((angle1>0 && angle1 < 1) && (angle2>0 && angle2 < 1) && (angle3>0 && angle3 < 1))
				{
					p1 = newpoint;
				}
				else
				{
					//2:newpoin1+p1 p3
					angle1 = getCosineAngle(newpoint, p1, newpoint, p3);
					angle2 = getCosineAngle(p1, newpoint, p1, p3);
					angle3 = getCosineAngle(p3, newpoint, p3, p1);
					if ((angle1 > 0 && angle1 < 1) && (angle2>0 && angle2 < 1) && (angle3>0 && angle3 < 1))
					{
						p2 = newpoint;
					}
					else
					{
						//3:newpoin1+p1 p2
						angle1 = getCosineAngle(newpoint, p1, newpoint, p2);
						angle2 = getCosineAngle(p1, newpoint, p1, p2);
						angle3 = getCosineAngle(p2, newpoint, p2, p1);
						if ((angle1 > 0 && angle1 < 1) && (angle2>0 && angle2 < 1) && (angle3>0 && angle3 < 1))
						{
							p3 = newpoint;
						}
						else
							cha = 0.0001;
					}
				}
			}
			//while done
			radius = tmp_r;
			center = tmp_center;
		}
		
	}
}


//a��ʾ��Ҫ�ϵ���data��λ�ã�aȡֵ��Χ1...N ����������h_index
int heap_siftup(int a, vector<int>* b, vector<double>* ele, vector<int>* h_index)
{
	int size = a;
	int flag = 0; //��������Ƿ���Ҫ�������ϵ���
	if (size == 1)
		return a; //����ǶѶ����ͷ��أ�����Ҫ������    
	//���ڶѶ� ���� ��ǰ���i��ֵ�ȸ����С��ʱ��������ϵ��� 
	while (size > 1 && flag == 0)
	{
		//�ж��Ƿ�ȸ�����С 
		int ff = size / 2;
		if (ele->at(b->at(size - 1)) < ele->at(b->at(ff - 1)))// <С����
		{
			//�ȸ���������h_index
			int tmp = h_index->at(b->at(size - 1));
			h_index->at(b->at(size - 1)) = h_index->at(b->at(ff - 1));
			h_index->at(b->at(ff - 1)) = tmp;
			//swap(size, size / 2);//����heap�е�λ��(�������ְֵ�λ�� )
			tmp = b->at(size - 1);
			b->at(size - 1) = b->at(ff - 1);
			b->at(ff - 1) = tmp;
			
			size = size / 2; //���±��iΪ�������ı�ţ��Ӷ�������һ�μ������ϵ��� 
		}
		else
			flag = 1;//��ʾ�Ѿ�����Ҫ�����ˣ���ǰ����ֵ�ȸ�����ֵҪ�� 
		//size = size / 2; //���±��iΪ�������ı�ţ��Ӷ�������һ�μ������ϵ��� 
	}
	return size;//���ص�ǰ���µ�λ��
}
//a��ʾ��Ҫ�µ���data��λ�ã�aȡֵ��Χ1...N ����������h_index
int heap_siftdown(int a, vector<int>* b, vector<double>* ele, vector<int>* h_index)
{
	int size = b->size();
	int i = a;//i����Ҫ�����Ľڵ�λ��  1...N
	if (size >1)
	{
		int t, flag = 0;//flag��������Ƿ���Ҫ�������µ��� 
		//��i����ж��ӵ�ʱ����ʵ������������ӵ�����£���������Ҫ����������ʱ��ѭ����ִ��
		while (i * 2 <= size && flag == 0)
		{
			//�����ж�����������ӵĹ�ϵ������t��¼ֵ��С�Ľ���� 
			if (ele->at(b->at(i - 1)) > ele->at(b->at(i * 2 - 1)))// >С����
				t = i * 2;
			else
				t = i;
			//��������Ҷ��ӵ�����£��ٶ��Ҷ��ӽ������� 
			if (i * 2 + 1 <= size)
			{
				//����Ҷ��ӵ�ֵ��С�����½�С�Ľ����  
				if (ele->at(b->at(t - 1)) >ele->at(b->at(i * 2)))//  >С����
					t = i * 2 + 1;
			}
			//���������С�Ľ���Ų����Լ���˵���ӽ�����бȸ�����С��  
			if (t != i)
			{
				//�Ƚ���������h_index
				int tmp = h_index->at(b->at(t - 1));
				//h_index->at(b->at(t - 1)) = h_index->at(b->at(i - 1));
				//h_index->at(b->at(i - 1)) = tmp;
				h_index->at(b->at(t - 1)) = i - 1;
				h_index->at(b->at(i - 1)) = t - 1;
				//swap(t, i);//��������
				tmp = b->at(t - 1);
				b->at(t - 1) = b->at(i - 1);
				b->at(i - 1) = tmp;

				i = t;//����iΪ�ղ����������Ķ��ӽ��ı�ţ����ڽ������������µ��� 
			}
			else
				flag = 1;//���˵����ǰ�ĸ�����Ѿ��������ӽ�㶼ҪС�ˣ�����Ҫ�ڽ��е����� 
		}
	}
	return i;
}
//����λ��index����dataΪnewdata��indexȡֵ��Χ1...N , ������λ��
int heap_update(int index, double newdata, vector<int>* b, vector<double>* ele, vector<int>* h_index)//indexȡֵ��Χ1...N
{
	int newposi = 0;
	ele->at(b->at(index - 1)) = newdata;
	//int flag = 0;//0��ʾ�����ƶ���1��ʾ�ϸ����ȸ��ڵ�С����2��ʾ�³���������һ���ӽڵ�� С����
	if (index > 1)
	{
		double d1 = ele->at(b->at(index - 1));
		double d2 = ele->at(b->at(index / 2 - 1));
		if (d1 < d2)//if (ele->at(b->at(index - 1)) < ele->at(b->at(index / 2 - 1)))�ȸ��ڵ� С���ϸ�
		{
			newposi = heap_siftup(index, b, ele, h_index);
		}
		else
		{
			//�������������ӽڵ�Ƚ�
			newposi = heap_siftdown(index, b, ele, h_index);
		}
	}
	else
	{
		newposi = heap_siftdown(index, b, ele, h_index);//index�Ǹ��ڵ㣬ֱ���³�
	}
	return newposi;
}


void getBoundRectangle(vector<mypoint2f>* total_plist, mypoint2f &centerpt, vector<mypoint2f>* four_corner)
{
	//PCA ���ɷַ����� ֱ������ϵ�� ���ȼ���Э�������A ,��ʼ�� ���ݵ����MatrixXd B��2��n�� ���ݰ������У���
	int plist_size = total_plist->size();
	vector<double> vec_x;	//vec_x��vec_y��B ����ֱ������ϵ�µģ�
	vector<double> vec_y;
	MatrixXd B = MatrixXd::Zero(2, plist_size);
	MatrixXd B_revolve = MatrixXd::Zero(2, plist_size);
	for (int j = 0; j < plist_size; ++j)
	{
		vec_x.push_back(total_plist->at(j).x);
		vec_y.push_back(total_plist->at(j).y);
		B(0, j) = total_plist->at(j).x;
		B(1, j) = total_plist->at(j).y;
	}
	//cout << B[0] << endl;
	double aver_x = GetAverange(vec_x);
	double aver_y = GetAverange(vec_y);
	double var_cxx = GetVariance(vec_x, aver_x, vec_x, aver_x);
	double var_cyy = GetVariance(vec_y, aver_y, vec_y, aver_y);
	double var_cxy = GetVariance(vec_x, aver_x, vec_y, aver_y);
	Matrix2d A;
	A << var_cxx, var_cxy, var_cxy, var_cyy;  //cout << "Here is a covariance 2*2 matrix, A:" << endl << A << endl << endl;
	//����Э������������ֵ ��������
	EigenSolver<MatrixXd> es(A);
	//cout << "The eigenvalues of A are:" << endl << es.eigenvalues() << endl;//ע����������������ֵ���� �������ͣ���complex
	//cout << "The eigenvectors is:" << endl << es.eigenvectors() << endl;//���������ǰ������е�,���Ѿ��ǵ�λ����
	//�������������double���͵�matrix��
	VectorXd an_vec;
	Matrix2d eigenvec_double = MatrixXd::Zero(2, 2);
	for (int m_i = 0; m_i<2; m_i++)
	{
		an_vec = es.eigenvectors().col(m_i).real();
		eigenvec_double.col(m_i) = an_vec;
	}
	// B�е����ݵ�����������������ת�ã���������任
	B_revolve = eigenvec_double.transpose()*B;  //(2*2) * (2*n)�� eigenvec_double.transpose()����ת�ú�ľ���eigenvec_double����û�б�
	//����������ϵ�µ� �������� �߽� 
	Vector2d apt_vec;//�Ƚ�ĳһ���������任����
	apt_vec << total_plist->at(0).x, total_plist->at(0).y;
	Vector2d apt_revolve = eigenvec_double.transpose()*apt_vec;  //����������������ת�ã�
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
	//����������ϵ�µ� ԭ��,
	mypoint2f revol_center_pt(0, 0);
	revol_center_pt.x = (posi_left + posi_right) / 2;
	revol_center_pt.y = (posi_bottom + posi_top) / 2;
	//���� ���� �� ����
	double box_w = fabs(posi_right - posi_left);
	double box_h = fabs(posi_top - posi_bottom);
	double elongation = 0;
	if (box_w < box_h)
		elongation = box_w / box_h;
	else
		elongation = box_h / box_w;
	/*showBoundingBox(),����bounding box�ĸ���*/
	vector<mypoint2f> plist_box;
	mypoint2f bp(revol_center_pt.x - box_w / 2, revol_center_pt.y + box_h / 2);//���Ͻ�
	plist_box.push_back(bp);
	bp.x = revol_center_pt.x + box_w / 2;  bp.y = revol_center_pt.y + box_h / 2;//���Ͻ�
	plist_box.push_back(bp);
	bp.x = revol_center_pt.x + box_w / 2;  bp.y = revol_center_pt.y - box_h / 2;//���½�
	plist_box.push_back(bp);
	bp.x = revol_center_pt.x - box_w / 2;  bp.y = revol_center_pt.y - box_h / 2;//���½�
	plist_box.push_back(bp);
	//4���������任�� ԭ������ϵ��
	MatrixXd box_matrix = MatrixXd::Zero(2, 5);//�ȴ浽������
	for (int j = 0; j < 4; j++)
	{
		box_matrix(0, j) = plist_box[j].x;
		box_matrix(1, j) = plist_box[j].y;
	}
	MatrixXd box_matrix_o = eigenvec_double*box_matrix;//����� ������������
	for (int j = 0; j < 4; j++)
	{
		mypoint2f pt(box_matrix_o(0, j), box_matrix_o(1, j));
		four_corner->push_back(pt);
		//plist_box[j].x = box_matrix_o(0, j);
		//plist_box[j].y = box_matrix_o(1, j);
	}
	//ԭ��ת���� ԭ��������ϵ
	Vector2d revol_center_vec;
	revol_center_vec << revol_center_pt.x, revol_center_pt.y;
	Vector2d center_vec1 = eigenvec_double*revol_center_vec;//��� ������������
	centerpt.x = center_vec1[0];
	centerpt.y = center_vec1[1];
}



//piecewise curve fitting
mypoint2f V2Scale(mypoint2f *p, double newlen)
{
	double len = sqrt(p->x*p->x + p->y*p->y);  //len=1 ��Ϊp�ǵ�λ����
	double ss = (newlen / len);
	if (len != 0.0) { p->x = p->x *ss;   p->y = p->y*ss; }
	
	mypoint2f v=*p;
	return(v);
}
mypoint2f V2ScaleIII(mypoint2f v, double s)
{
	mypoint2f result;
	result.x = v.x * s; result.y = v.y * s;
	return (result);
}
mypoint2f V2AddII(mypoint2f a, mypoint2f b)
{
	mypoint2f	c;
	c.x = a.x + b.x;  c.y = a.y + b.y;
	return (c);
}
mypoint2f V2SubII(mypoint2f a, mypoint2f b)
{
	mypoint2f	c;
	c.x = a.x - b.x; c.y = a.y - b.y;
	return (c);
}
mypoint2f *V2Negate(mypoint2f *v)
{
	v->x = -v->x;  v->y = -v->y;
	return(v);
}

mypoint2f ComputeLeftTangent(vector<mypoint2f>* d, int end)
/*	Point2	*d;			  Digitized points*/
/*int		end;		  Index to "left" end of region */
{
	mypoint2f tHat1(0, 0);
	//tHat1 = d->at(end + 1) - d->at(end); //V2SubII(d[end + 1], d[end]);
	tHat1 = normalization(d->at(end), d->at(end + 1));
	return tHat1;
}
mypoint2f ComputeRightTangent(vector<mypoint2f>* d, int end)
/*	Point2	*d;			  Digitized points		*/
/*int		end;		 Index to "right" end of region */
{
	mypoint2f tHat2;
	//tHat2 = V2SubII(d[end - 1], d[end]);
	tHat2 = normalization(d->at(end), d->at(end - 1));
	return tHat2;
}
mypoint2f ComputeCenterTangent(vector<mypoint2f>* d, int center)
/*	Point2	*d;			  Digitized points			*/
/*int		center;		  Index to point inside region	*/
{
	mypoint2f	V1(0, 0), V2(0, 0), tHatCenter(0, 0);
	V1 = normalization(d->at(center), d->at(center - 1));  //V2SubII(d[center - 1], d[center]);
	V2 = normalization(d->at(center + 1), d->at(center)); //V2SubII(d[center], d[center + 1]);
	tHatCenter.x = (V1.x + V2.x) / 2.0;
	tHatCenter.y = (V1.y + V2.y) / 2.0;
	double len = LineSegmentLength(mypoint2f(0, 0), tHatCenter);

	tHatCenter.x = tHatCenter.x / len;
	tHatCenter.y = tHatCenter.y / len;
	return tHatCenter;
}

void ChordLengthParameterize(vector<mypoint2f>* d, int first, int last, vector<double> *u)
/*	Point2	*d;			 Array of digitized points */
/*int		first, last;		 Indices defining region	*/
{
	//double	*u;			// Parameterization
	//u = (double *)malloc((unsigned)(last - first + 1) * sizeof(double));
	//u[0] = 0.0;
	int vecsize = last - first + 1;
	for (int i = 0; i < vecsize; i++)
	{
		u->push_back(0.0);
	}
	for (int i = first + 1; i <= last; i++) 
	{
		u->at(i - first) = u->at(i - first - 1) + LineSegmentLength(d->at(i), d->at(i - 1));
	}
	for (int i = first + 1; i <= last; i++)
	{
		u->at(i - first) = u->at(i - first) / u->at(last - first);
	}
	//int vecsize = last - first + 1;
	//u->push_back(0.0);
	//double accu_len = 0;
	//for (int i = 1; i < vecsize; ++i)
	//{
	//	accu_len = accu_len + LineSegmentLength(d->at(i+first-1), d->at(i + first));
	//	u->push_back(accu_len);
	//}
	//for (int i =0; i < vecsize; i++)
	//{
	//	u->at(i) = u->at(i) / accu_len;
	//}
}

double B0(double u)
{
	double tmp = 1.0 - u;
	return (tmp * tmp * tmp);
}
double B1(double u)
{
	double tmp = 1.0 - u;
	return (3 * u * (tmp * tmp));
}
double B2(double u)
{
	double tmp = 1.0 - u;
	return (3 * u * u * tmp);
}
double B3(double u)
{
	return (u * u * u);
}

/*
*  GenerateBezier :
*  Use least-squares method to find Bezier control points for region.
*/
cubicBezier GenerateBezier(vector<mypoint2f>* d, int first, int last, vector<double>* uPrime, mypoint2f tHat1, mypoint2f tHat2)
/*	Point2	*d;			  Array of digitized points	*/
/*	int		first, last;		  Indices defining region	*/
/*	double	*uPrime;		  Parameter values for region */
/*	Vector2	tHat1, tHat2;	  Unit tangents at endpoints	*/
{
	int 	i;
	vector<vector<mypoint2f>> A; //Vector_pt A[MAXPOINTS][2]; /* Precomputed rhs for eqn�� A�� ��nPts = last - first + 1���У� 2�� �ľ���*/
	int 	nPts;			/* Number of pts in sub-curve */
	double 	C[2][2];			/* Matrix C	 Э�������	*/
	double 	X[2];			/* Matrix X			*/
	double 	det_C0_C1,		/* Determinants of matrices	*/
		det_C0_X,
		det_X_C1;
	double 	alpha_l,		/* Alpha values, left and right	*/
		alpha_r;
	cubicBezier	bezCurve;	/* RETURN bezier curve ctl pts	*/
	double  segLength;
	double  epsilon;

	nPts = last - first + 1;

	/* Compute the A's	*/
	for (i = 0; i < nPts; i++) {
		mypoint2f v1, v2;
		v1 = tHat1;
		v2 = tHat2;
		V2Scale(&v1, B1(uPrime->at(i)));
		V2Scale(&v2, B2(uPrime->at(i)));
		//A[i][0] = v1;
		//A[i][1] = v2;
		vector<mypoint2f> tmp_vec;
		tmp_vec.push_back(v1); tmp_vec.push_back(v2);
		A.push_back(tmp_vec);
	}

	/* Create the C and X matrices	*/
	C[0][0] = 0.0;
	C[0][1] = 0.0;
	C[1][0] = 0.0;
	C[1][1] = 0.0;
	X[0] = 0.0;
	X[1] = 0.0;

	for (i = 0; i < nPts; i++) {
		C[0][0] += dot(A[i][0], A[i][0]); //V2Dot(&A[i][0], &A[i][0]);
		C[0][1] += dot(A[i][0], A[i][1]); //V2Dot(&A[i][0], &A[i][1]);
		C[1][0] = C[0][1]; 		/*	C[1][0] += V2Dot(&A[i][0], &A[i][1]);*/
		C[1][1] += dot(A[i][1], A[i][1]);//V2Dot(&A[i][1], &A[i][1]);

		mypoint2f tmp = V2SubII(d->at(first + i),
			V2AddII(
				V2ScaleIII(d->at(first), B0(uPrime->at(i))),
				V2AddII(
					V2ScaleIII(d->at(first), B1(uPrime->at(i))),
					V2AddII(
						V2ScaleIII(d->at(last), B2(uPrime->at(i))),
						V2ScaleIII(d->at(last), B3(uPrime->at(i)))))));

		//mypoint2f tmp1 = V2ScaleIII(d->at(last), B2(uPrime->at(i))) + V2ScaleIII(d->at(last), B3(uPrime->at(i)));
		//mypoint2f tmp2 = V2ScaleIII(d->at(first), B1(uPrime->at(i))) + tmp1;
		//mypoint2f tmp3 = V2ScaleIII(d->at(first), B0(uPrime->at(i))) + tmp2;
		//mypoint2f tmp4 = d->at(first + i) - tmp3;  // tmp4=tmp

		X[0] += dot(A[i][0], tmp);  //V2Dot(&A[i][0], &tmp);
		X[1] += dot(A[i][1], tmp);  //V2Dot(&A[i][1], &tmp);
	}

	/* Compute the determinants of C and X	*/
	det_C0_C1 = C[0][0] * C[1][1] - C[1][0] * C[0][1];
	det_C0_X = C[0][0] * X[1] - C[1][0] * X[0];
	det_X_C1 = X[0] * C[1][1] - X[1] * C[0][1];

	/* Finally, derive alpha values	*/
	alpha_l = (det_C0_C1 == 0) ? 0.0 : det_X_C1 / det_C0_C1;
	alpha_r = (det_C0_C1 == 0) ? 0.0 : det_C0_X / det_C0_C1;

	/* If alpha negative, use the Wu/Barsky heuristic (see text) */
	/* (if alpha is 0, you get coincident control points that lead to
	* divide by zero in any subsequent NewtonRaphsonRootFind() call. */
	segLength = LineSegmentLength(d->at(last), d->at(first)); //2DistanceBetween2Points(&d[last], &d[first]);
	epsilon = 1.0e-6 * segLength;
	if (alpha_l < epsilon || alpha_r < epsilon)
	{
		/* fall back on standard (probably inaccurate) formula, and subdivide further if needed. */
		double dist = segLength / 3.0;
		bezCurve.p0 = d->at(first);  //bezCurve[0] = d[first];
		bezCurve.p3 = d->at(last);   //bezCurve[3] = d[last];
		bezCurve.p1 = bezCurve.p0 + V2Scale(&tHat1, dist);  //V2Add(&bezCurve[0], V2Scale(&tHat1, dist), &bezCurve[1]);
		bezCurve.p2 = bezCurve.p3 + V2Scale(&tHat2, dist);  //V2Add(&bezCurve[3], V2Scale(&tHat2, dist), &bezCurve[2]);
		return (bezCurve);
	}

	/*  First and last control points of the Bezier curve are */
	/*  positioned exactly at the first and last data points */
	/*  Control points 1 and 2 are positioned an alpha distance out */
	/*  on the tangent vectors, left and right, respectively */
	bezCurve.p0 = d->at(first);  //bezCurve[0] = d[first];
	bezCurve.p3 = d->at(last);   //bezCurve[3] = d[last];
	bezCurve.p1 = bezCurve.p0 + V2Scale(&tHat1, alpha_l);  //	V2Add(&bezCurve[0], V2Scale(&tHat1, alpha_l), &bezCurve[1]);
	bezCurve.p2 = bezCurve.p3 + V2Scale(&tHat2, alpha_r);  //	V2Add(&bezCurve[3], V2Scale(&tHat2, alpha_r), &bezCurve[2]);
	return (bezCurve);
}
/*
*  Bezier :
*  	Evaluate a Bezier curve at a particular parameter value
*/
mypoint2f BezierII(int degree, vector<mypoint2f>* V, double t)
/*	int		degree;		 The degree of the bezier curve	*/
/*Point2 	*V;		 Array of control points		*/
/*double 	t;		 Parametric value to find point for	*/
{
	mypoint2f Q;	        /* Point on curve at parameter t	*/
	vector<mypoint2f> Vtemp;	/* Local copy of control points		*/
	Vtemp = *V;                 /* Copy array	*/
	/* Triangle computation	*/
	for (int i = 1; i <= degree; i++) {
		for (int j = 0; j <= degree - i; j++) {
			Vtemp[j].x = (1.0 - t) * Vtemp[j].x + t * Vtemp[j + 1].x;
			Vtemp[j].y = (1.0 - t) * Vtemp[j].y + t * Vtemp[j + 1].y;
		}
	}
	Q = Vtemp[0];
	return Q;
}
/*
*  ComputeMaxError :
*	Find the maximum squared distance of digitized points
*	to fitted curve.
*/
double ComputeMaxError(vector<mypoint2f>* d, int first, int last, cubicBezier bezCurve, vector<double>* u, int* splitPoint)
/*	Point2	*d;			  Array of digitized points	*/
/*int		first, last;		  Indices defining region	*/
/*BezierCurve	bezCurve;		 Fitted Bezier curve		*/
/*double	*u;			  Parameterization of points	*/
/*int		*splitPoint;		  Point of maximum error	*/
{
	int		i;
	double	maxDist;		/*  Maximum error		*/
	double	dist;		/*  Current error		*/
	mypoint2f P;			/*  Point on curve		*/
	mypoint2f v;			/*  Vector from point to curve	*/

	*splitPoint = (last - first + 1) / 2;
	maxDist = 0.0;
	for (i = first + 1; i < last; i++)
	{
		vector<mypoint2f> bc;
		bc.push_back(bezCurve.p0);
		bc.push_back(bezCurve.p1);
		bc.push_back(bezCurve.p2);
		bc.push_back(bezCurve.p3);

		P = BezierII(3, &bc, u->at(i - first));
		v = V2SubII(P, d->at(i));
		dist = v.x*v.x + v.y*v.y;  //V2SquaredLength(&v);
		if (dist >= maxDist) 
		{
			maxDist = dist;
			*splitPoint = i;
		}
	}
	return (maxDist);
}

/*
*  NewtonRaphsonRootFind :
*	Use Newton-Raphson iteration to find better root.
*/
double NewtonRaphsonRootFind(cubicBezier Q, mypoint2f P, double u)
/*	BezierCurve	Q;			  Current fitted curve	*/
/*Point2 		P;		  Digitized point		*/
/*double 		u;		  Parameter value for "P"	*/
{
	double 		numerator, denominator;
	vector<mypoint2f> Q1, Q2;      	//Point2 Q1[3], Q2[2];	/*  Q' and Q''	*/
	mypoint2f Q_u, Q1_u, Q2_u;      //	Point2	Q_u, Q1_u, Q2_u; /*u evaluated at Q, Q', & Q''	*/
	double 		uPrime;		/*  Improved u			*/

	/* Compute Q(u)	*/
	vector<mypoint2f> tmp_Q;
	tmp_Q.push_back(Q.p0);
	tmp_Q.push_back(Q.p1);
	tmp_Q.push_back(Q.p2);
	tmp_Q.push_back(Q.p3);
	Q_u = BezierII(3, &tmp_Q, u);

	/* Generate control vertices for Q'	���浽Q1��*/
	for (int i = 0; i <= 2; i++) 
	{
		mypoint2f pt((tmp_Q.at(i + 1).x - tmp_Q.at(i).x) * 3.0, (tmp_Q.at(i + 1).y - tmp_Q.at(i).y) * 3.0);
		Q1.push_back(pt);
		//Q1[i].x = (Q[i + 1].x - Q[i].x) * 3.0;
		//Q1[i].y = (Q[i + 1].y - Q[i].y) * 3.0;
	}

	/* Generate control vertices for Q'' ���浽Q2��*/
	for (int i = 0; i <= 1; i++)
	{
		mypoint2f pt((Q1.at(i + 1).x - Q1.at(i).x) * 2.0, (Q1.at(i + 1).y - Q1.at(i).y) * 2.0);
		Q2.push_back(pt);
		//Q2[i].x = (Q1[i + 1].x - Q1[i].x) * 2.0;
		//Q2[i].y = (Q1[i + 1].y - Q1[i].y) * 2.0;
	}

	/* Compute Q'(u) and Q''(u)	*/
	Q1_u = BezierII(2, &Q1, u);
	Q2_u = BezierII(1, &Q2, u);

	/* Compute f(u)/f'(u) */
	numerator = (Q_u.x - P.x) * (Q1_u.x) + (Q_u.y - P.y) * (Q1_u.y);
	denominator = (Q1_u.x) * (Q1_u.x) + (Q1_u.y) * (Q1_u.y) +
		(Q_u.x - P.x) * (Q2_u.x) + (Q_u.y - P.y) * (Q2_u.y);
	if (denominator == 0.0f)
		return u;
	/* u = u - f(u)/f'(u) */
	uPrime = u - (numerator / denominator);
	return (uPrime);
}
/*
*  Reparameterize:
*	Given set of points and their parameterization, try to find
*   a better parameterization.
*
*/
void Reparameterize(vector<mypoint2f>* d, int first, int last, vector<double>* u, cubicBezier bezCurve, vector<double> *uPrime)
/*	Point2	*d;			  Array of digitized points	*/
/*int		first, last;		  Indices defining region	*/
/*double	*u;			  Current parameter values	*/
/*BezierCurve	bezCurve;	  Current fitted curve	*/
{
	int 	nPts = last - first + 1;
	int 	i;
	//vector<double> uPrime;		/*  New parameter values	*/
	for (i = first; i <= last; i++) 
	{
		//uPrime[i - first] = NewtonRaphsonRootFind(bezCurve, d[i], u[i - first]);
		double tmp= NewtonRaphsonRootFind(bezCurve, d->at(i), u->at(i - first));
		uPrime->push_back(tmp);
	}
}



void FitCubic(vector<mypoint2f>* d, int first, int last, mypoint2f tHat1, mypoint2f tHat2, double error, vector<cubicBezier>* piecewise_cb)
/*	Point2	*d;			  Array of digitized points */
/*	int		first, last;	 Indices of first and last pts in region */
/*	Vector2	tHat1, tHat2;	 Unit tangent vectors at endpoints */
/*	double	error;		 User-defined error squared	   */
{
	cubicBezier	bezCurve; /*Control points of fitted Bezier curve*/
	vector<double> u; //double	*u;		/*  Parameter values for point  */
	//vector<double> uPrime;	/*  Improved parameter values */
	double	maxError;	/*  Maximum fitting error	 */
	int	splitPoint;	/*  Point to split point set at	 */
	int nPts;		/*  Number of points in subset  */
	double iterationError; /*Error below which you try iterating  */
	int maxIterations = 4; /*  Max times to try iterating  */
	mypoint2f tHatCenter;   	/* Unit tangent vector at splitPoint */

	iterationError = error * error;
	nPts = last - first + 1;

	/*  Use heuristic if region only has two points in it */
	if (nPts == 2) {
		double dist = LineSegmentLength(d->at(last), d->at(first)) / 3.0;

		bezCurve.p0 = d->at(first);
		bezCurve.p3 = d->at(last);
		bezCurve.p1 = bezCurve.p0 + V2Scale(&tHat1, dist);  //V2Add(&bezCurve[0], V2Scale(&tHat1, dist), &bezCurve[1]);
		bezCurve.p2 = bezCurve.p3 + V2Scale(&tHat2, dist); //V2Add(&bezCurve[3], V2Scale(&tHat2, dist), &bezCurve[2]);
		piecewise_cb->push_back(bezCurve);      //DrawBezierCurve(3, bezCurve);
		return;
	}
	/*  Parameterize points, and attempt to fit curve */
	ChordLengthParameterize(d, first, last, &u); //u = ChordLengthParameterize(d, first, last);
	bezCurve = GenerateBezier(d, first, last, &u, tHat1, tHat2);
	/*  Find max deviation of points to fitted curve */
	maxError = ComputeMaxError(d, first, last, bezCurve, &u, &splitPoint);
	if (maxError < error) {
		piecewise_cb->push_back(bezCurve); //DrawBezierCurve(3, bezCurve);
		//free((void *)u);
		//free((void *)bezCurve);/////////////////////ԭ���У���ע�͵���
		return;
	}

	/*  If error not too large, try some reparameterization  */
	/*  and iteration */
	if (maxError < iterationError)
	{
		vector<double> uPrime;	/*  Improved parameter values */
		for (int i = 0; i < maxIterations; i++)
		{
			Reparameterize(d, first, last, &u, bezCurve, &uPrime);
			//free((void *)bezCurve);
			bezCurve = GenerateBezier(d, first, last, &uPrime, tHat1, tHat2);
			maxError = ComputeMaxError(d, first, last,bezCurve, &uPrime, &splitPoint);
			if (maxError < error)
			{
				piecewise_cb->push_back(bezCurve);//DrawBezierCurve(3, bezCurve);
				//free((void *)u);
				//free((void *)bezCurve);
				//free((void *)uPrime);
				return;
			}
			u.swap(vector<double>());  //free((void *)u);
			u = uPrime;  //!!!??
		}
	}

	/* Fitting failed -- split at max error point and fit recursively */
	//free((void *)u);
	//free((void *)bezCurve);
	tHatCenter = ComputeCenterTangent(d, splitPoint);
	FitCubic(d, first, splitPoint, tHat1, tHatCenter, error, piecewise_cb);
	V2Negate(&tHatCenter);
	FitCubic(d, splitPoint, last, tHatCenter, tHat2, error, piecewise_cb);
}

void FitCurves(vector<mypoint2f>* d, int nPts, double error, vector<cubicBezier>* piecewise_cb)
/*	vector<mypoint2f>	*d;			  Array of digitized points	*/
/*  int		nPts;		  Number of digitized points	*/
/*  double	error;		  User-defined error squared	*/
{
	mypoint2f tHat1, tHat2;	/*  Unit tangent vectors at endpoints */
	tHat1 = ComputeLeftTangent(d, 0);
	tHat2 = ComputeRightTangent(d, nPts - 1);

	FitCubic(d, 0, nPts - 1, tHat1, tHat2, error,piecewise_cb);
}



