#include "MySvgFile.h"

float MySvgFile::zoomup=0.9;
float MySvgFile::zoomdown=1.1;
int MySvgFile::winWidth = 1400;
float MySvgFile::winHeight=800;
float MySvgFile::w_x1 = 0;
float MySvgFile::w_x2 = winWidth;
float MySvgFile::w_y1 = 0;
float MySvgFile::w_y2 = winHeight;
//double MySvgFile::errorDeviation = 0.5;
MySvgFile *MySvgFile::svginstance = NULL;

MySvgFile::MySvgFile() 
{
	filename = " ";
	filedoc=NULL;
	svg_width=0;
	svg_height=0;
	 svg_view = " ";
	//originCoordinate.x = 0;
	//originCoordinate.y = 0;
	tree_depth = 0;
	//zoomup = 0.9;  zoomdown = 1.1;
	//winWidth = 750;  winHeight = 1000;
	//w_x1 = 0;
	//w_x2 = winWidth;//= winWidth;
	//w_y1 = 0;
	//w_y2 = winHeight;// = winHeight;
}
MySvgFile::~MySvgFile()
{
}
MySvgFile::MySvgFile(string name) 
{
	filename = name;
	filedoc = new TiXmlDocument(name.c_str());
	svg_width = 0;
	svg_height = 0;
	 svg_view=" ";

	setInstance();
	//originCoordinate.x = 0;
	//originCoordinate.y = 0;
	bool isload = filedoc->LoadFile();
	if (!isload)
	{
		cout << "Could not load test file:" << filename<<",Error="<< filedoc->ErrorDesc() << endl;
	}
	tree_depth = 0;
	//zoomup = 0.9;  zoomdown = 1.1;
	//winWidth = 750;  winHeight = 1000;
	//w_x1 = 0;
	//w_x2 = winWidth;//= winWidth;
	//w_y1 = 0;
	//w_y2 = winHeight;// = winHeight;
}

void MySvgFile::getGraphPrimitive(double *fwidth, double *fheight, vector<double> *view, string &view_str)
{
	TiXmlElement *rootElement = filedoc->RootElement();
	cout << "rootelement.name=" << rootElement->Value() << endl;
	//TiXmlElement *ele_g = rootElement->FirstChildElement("g");
	//cout << "ele_g.name=" << ele_g->Value() << ",attribute=" << ele_g->FirstAttribute()->Name() << "," << ele_g->FirstAttribute()->Value()<< endl;
	//找到<svg  />属性里的width height viewbox:
	const char *const_view = rootElement->Attribute("viewBox");//name="transform",value=graphelement->Attribute("transform")
	if (const_view)
	{
		string str_view = const_cast<char *>(const_view);
		view_str = str_view;
		str_view = trimSpace(str_view);
		*view = strToNum(str_view, 0, str_view.size() - 1);
	}
	else
	{
		view->push_back(0); view->push_back(0); view->push_back(1500); view->push_back(1500);
		view_str = "0 0 1500 1500";
	}
	const char *const_w = rootElement->Attribute("width");//name="transform",value=graphelement->Attribute("transform")
	if (const_w)
	{
		string str_w = const_cast<char *>(const_w);
		str_w = trimSpace(str_w);
		*fwidth = getFileSize(str_w);
	}
	else
		*fwidth = 1500;
	const char *const_h = rootElement->Attribute("height");//name="transform",value=graphelement->Attribute("transform")
	if (const_h)
	{
		string str_h = const_cast<char *>(const_h);
		str_h = trimSpace(str_h);
		*fheight = getFileSize(str_h);
	}
	else
		*fheight == 1500;
	svg_width = *fwidth;
	svg_height = *fheight;
	svg_view = view_str;

	char *lastGraphName = " ";
	stack<TiXmlElement*> ele_stack;
	ele_stack.push(rootElement);//push root into the stack
	tree_depth = 0;
	while (!ele_stack.empty())
	{
		TiXmlElement *tmpelement = ele_stack.top();
		ele_stack.pop();
		TiXmlElement *childnode = tmpelement->FirstChildElement();
		tree_depth = tree_depth + 1;
		while (childnode)
		{		
			if (childnode->ValueTStr() == "g")//if it is <g>,put into stack
			{
				ele_stack.push(childnode);
				getTransform(childnode);
				//transOrigiCooridate(childnode);//Generally "transform" or "matric"attribute exits in <g>,if it has,then change the origin coordinate!
			}
			else//path or line or circle...
			{
				//count_path = count_path + 1;
				//lastGraphName = childnode->Value();
				//cout << childnode->Value() << endl;
				try
				{
					getDifferentGraphics(childnode, lastGraphName);//get Different Graphics according to their tag
				}
				catch (int)
				{
					cout << "Do not find conrresponding shape in 'shape_ELement[6]' " << endl;
				}
				catch (double)
				{
					cout << "Not find a double attribute in tag<line> or other tag" << endl;
				}
			}
			childnode = childnode->NextSiblingElement();
		}
	}
	cout << endl;
}
double MySvgFile::getFileSize(string str)
{
	string numstr;
	for (int i = 0;i< str.length(); i++)
	{
		int chr = str.at(i);
		if (chr >= 48 && chr <= 57)
			numstr.push_back(str[i]);
		else
			break;
	}
	double tmp = atof(numstr.c_str());
	return tmp;
}
void MySvgFile::getTransform(TiXmlElement *graphelement)
{
	const char *const_trans = graphelement->Attribute("transform");//name="transform",value=graphelement->Attribute("transform")
	string tran_name;
	vector<double> tran_number;
	if (const_trans)
	{		
		string str_trans = const_cast<char *>(const_trans);
		//2 situations:
		// <g transform="translate(5,10)">
		//<g transform="matrix(0.71808638,0,0,0.71808638,-2301.4416,-2305.7758)">
		cout << "str_transform=" << str_trans.c_str() << endl;
		if (str_trans.find("translate") != string::npos)
		{
			tran_name = "translate";
			int leftp = str_trans.find("(");
			int midComma = str_trans.find(",");
			int rightp = str_trans.find(")");
			tran_number.push_back(atof(str_trans.substr(leftp + 1, midComma - leftp - 1).c_str()));
			tran_number.push_back(atof(str_trans.substr(midComma + 1, rightp - midComma - 1).c_str()));
			transform_vec.push_back(make_pair(tran_name, tran_number));
			//originCoordinate.x = originCoordinate.x + atof(str_trans.substr(leftp + 1, midComma - leftp - 1).c_str());//substr(start,length);start included
			//originCoordinate.y = originCoordinate.y + atof(str_trans.substr(midComma + 1, rightp - midComma - 1).c_str());
		}
		else if (str_trans.find("matrix") != string::npos)
		{
			tran_name = "matrix";
			int leftp = str_trans.find("(");
			int rightp = str_trans.find(")");
			string middlestr = str_trans.substr(leftp + 1, rightp - leftp - 1);		
			tran_number = strToNum(middlestr, 0, middlestr.size() - 1);//store parameters"a,b,c,d,e,f" from matrix
			transform_vec.push_back(make_pair(tran_name, tran_number));
			/*double tt_x = Matrix_vec[0] * originCoordinate.x + Matrix_vec[2] * originCoordinate.y + Matrix_vec[4];
			double tt_y = Matrix_vec[1] * originCoordinate.x + Matrix_vec[3] * originCoordinate.y + Matrix_vec[5];
			originCoordinate.x = tt_x;
			originCoordinate.y = tt_y;*/
		}
	}
	else
	{
		tran_name=" ";
		transform_vec.push_back(make_pair(tran_name, tran_number));//tran_name=""tran_numbe is empty
	}
}
vector<mypoint2f> MySvgFile::transOrigiCooridate(vector<mypoint2f> tp_vec)
{
	//int tree_depth;
	//map<string, vector<double>> transform_vec;
	if (!transform_vec.empty())
	{
		double temp_x = 0;
		double temp_y = 0;
		vector<pair<string, vector<double>>>::iterator iter = transform_vec.end();
		for (unsigned int i = 0; i < tree_depth - 1; i++)
		{
			--iter;
			if (iter->first == "translate")
			{
				for (int j = 0; j < tp_vec.size(); j++)
				{
					tp_vec[j].x = tp_vec[j].x + iter->second.at(0);
					tp_vec[j].y = tp_vec[j].y + iter->second.at(1);
				}
			}
			else if (iter->first == "matrix")
			{
				for (int j = 0; j < tp_vec.size(); j++)
				{
					temp_x = tp_vec[j].x * iter->second.at(0) + tp_vec[j].y*iter->second.at(2) + iter->second.at(4);
					temp_y = tp_vec[j].x * iter->second.at(1) + tp_vec[j].y*iter->second.at(3) + iter->second.at(5);
					tp_vec[j].x = temp_x;
					tp_vec[j].y = temp_y;
				}
			}
		}
		/*map<string, vector<double>>::iterator iter = transform_vec.begin();
		for (unsigned int i = 0; i < tree_depth-1; i++)
		{
			if (iter->first == "translate")
			{
				for (int j = 0; j < tp_vec.size(); j++)
				{
					tp_vec[j].x = tp_vec[j].x + iter->second.at(0);
					tp_vec[j].y = tp_vec[j].y + iter->second.at(1);
				}
			}
			else if (iter->first == "matrix")
			{
				for (int j = 0; j < tp_vec.size(); j++)
				{
					temp_x = tp_vec[j].x * iter->second.at(0) + tp_vec[j].y*iter->second.at(2) + iter->second.at(4);
					temp_y = tp_vec[j].x * iter->second.at(1) + tp_vec[j].y*iter->second.at(3) + iter->second.at(5);
					tp_vec[j].x = temp_x;
					tp_vec[j].y = temp_y;
				}
			}
			if (i >= transform_vec.size()-1)
				break;
			else
				++iter;
		}*/
	}
	return tp_vec;
}
mypoint2f MySvgFile::transOrigiCooridate(mypoint2f tp)
{
	if (!transform_vec.empty())
	{
		double temp_x = 0;
		double temp_y = 0;
		vector<pair<string, vector<double>>>::iterator iter = transform_vec.end();
		for (unsigned int i = 0; i < tree_depth - 1; i++)
		{
			--iter;
			if (iter->first == "translate")
			{
				tp.x = tp.x + iter->second.at(0);
				tp.y = tp.y + iter->second.at(1);
			}
			else if (iter->first == "matrix")
			{
				temp_x = tp.x * iter->second.at(0) + tp.y*iter->second.at(2) + iter->second.at(4);
				temp_y = tp.x * iter->second.at(1) + tp.y*iter->second.at(3) + iter->second.at(5);
				tp.x = temp_x;
				tp.y = temp_y;
			}
		}
	}
	return tp;
}
//void MySvgFile::transOrigiCooridate(TiXmlElement *graphelement)
//{
//	//原 ，计算原点坐标的变换，originCoordinate.y
//	const char *const_trans = graphelement->Attribute("transform");//name="transform",value=graphelement->Attribute("transform")
//	if (const_trans)
//	{
//		string str_trans = const_cast<char *>(const_trans);
//		//2 situations:
//		// <g transform="translate(5,10)">
//		//<g transform="matrix(0.71808638,0,0,0.71808638,-2301.4416,-2305.7758)">
//		cout << "str_transform=" << str_trans.c_str() << endl;
//		if (str_trans.find("translate") != string::npos)
//		{
//			int leftp = str_trans.find("(");
//			int midComma = str_trans.find(",");
//			int rightp = str_trans.find(")");
//			originCoordinate.x = originCoordinate.x + atof(str_trans.substr(leftp + 1, midComma - leftp - 1).c_str());//substr(start,length);start included
//			originCoordinate.y = originCoordinate.y + atof(str_trans.substr(midComma + 1, rightp - midComma - 1).c_str());
//		}
//		else if (str_trans.find("matrix") != string::npos)
//		{
//			int leftp = str_trans.find("(");
//			int rightp = str_trans.find(")");
//			string middlestr = str_trans.substr(leftp + 1, rightp - leftp - 1);
//			vector<double> transMatrix_vec;//store parameters"a,b,c,d,e,f" from matrix
//			transMatrix_vec = strToNum(middlestr, 0, middlestr.size() - 1);
//			double tt_x = transMatrix_vec[0] * originCoordinate.x + transMatrix_vec[2] * originCoordinate.y + transMatrix_vec[4];
//			double tt_y = transMatrix_vec[1] * originCoordinate.x + transMatrix_vec[3] * originCoordinate.y + transMatrix_vec[5];
//			originCoordinate.x = tt_x;
//			originCoordinate.y = tt_y;
//			cout << transMatrix_vec[0] << "," << transMatrix_vec[2] << "," << transMatrix_vec[4] << endl;
//		}
//		cout << "original coordinate is changed:" << originCoordinate.x << "," << originCoordinate.y << endl;
//	}
//}

void MySvgFile::getDifferentGraphics(TiXmlElement *graphelement, char *g_type)
{
	const char * const_graph_type=graphelement->Value();
	char *graph_type = const_cast<char *>(const_graph_type);
	//if (strcmp(graph_type, g_type) != 0)//test
	//{
	//	g_type = const_cast<char *>(const_graph_type);
	//	//cout << "graph type:" << g_type << endl;
	//}
	//
	int typenum = shapeStrTNum(graph_type);
	if (typenum < 0)
	{
		throw typenum;
	}
	else
	{
		switch (typenum)  //const char *shapeType[6] = {"path","line","circle","rect","polyline","polygon"}
		{
			case 0://path
			{
				divideShapeFromPath(graphelement);//subdivision..<path d="M...L...Q...C..." /> continue split graphics from path
				break;
			}
			case 1://line
			{
				shape_Line(graphelement);
				break;
			}
			case 2://circle
			{
				shape_Circle(graphelement);
				break;
			}
			case 3://rect
			{
				shape_Rect(graphelement);
				break;
			}
			case 4://polyline
			{
				shape_Polyline(graphelement);
				break;
			}
			case 5:
			{
				shape_Polygon(graphelement);
				break;
			}
		}
	}
	
}
int MySvgFile::shapeStrTNum(char *shapename)
{
	char *shape_ELement[6] = { "path", "line", "circle", "rect", "polyline", "polygon" };
	int shapenum = -1;
	for (int i = 0; i < 6; ++i)
	{
		if (strcmp(shape_ELement[i], shapename) == 0)
		{
			shapenum = i;
		}
	}
	return shapenum;
}

void MySvgFile::divideShapeFromPath(TiXmlElement *graphelement)
{
	//<path d = "M...L...Q...C..." / > for example: d="m 1171.1176,585.56556 q 140.4006,229.92032 175.9721,532.56874" 

	const char *d_const = graphelement->Attribute("d");
	string d_str = d_const;
	char *typestr = "mMlLqQtTcCsShHvVzZ";
	int found_pos = d_str.find_first_of(typestr);//<path d="m/M..."> must start with m/M. Get character 'm'/'M' firstly.
	int next_index = d_str.find_first_of(typestr, found_pos + 1);

	int temp_index = found_pos;
	vector<pair<char,int>> subElement;
	while (temp_index != string::npos)
	{
		subElement.push_back(make_pair(d_str[temp_index], temp_index));
		temp_index = d_str.find_first_of(typestr, temp_index + 1);
	}
	if (subElement.back().first == 'z'|| subElement.back().first == 'Z')
		subElement.pop_back();
	if (subElement.size() == 1)
	{
		//get the first start point
		string leftstr = d_str.substr(found_pos + 1);//substract the middle str which is the start point
		leftstr = trimSpace(leftstr);//remove the first and last space
		vector<double> xy_vec = strToNum(leftstr, 0, leftstr.size() - 1);
		mypoint2f startpoint;
		startpoint.x = xy_vec[0];//startpoint.x = xy_vec[0] + originCoordinate.x;
		startpoint.y = xy_vec[1];//startpoint.y = xy_vec[1] + originCoordinate.y;

		vector<double>::iterator iter = xy_vec.begin();
		xy_vec.erase(iter); 
		iter = xy_vec.begin();
		xy_vec.erase(iter);
		cout << "only 'm/M'" << endl;
		path_LineTo(startpoint, xy_vec, d_str[found_pos]);//the left point consists of polyline
	}
	else
	{
		//get the first start point
		string midstr = d_str.substr(found_pos + 1, next_index-1);//substract the middle str which is the start point
		midstr = trimSpace(midstr);//remove the first and last space
		vector<double> xy_vec = strToNum(midstr, 0, midstr.size() - 1);
		mypoint2f startpoint;
		startpoint.x = xy_vec[0];//startpoint.x = xy_vec[0] + originCoordinate.x;
		startpoint.y = xy_vec[1];//startpoint.y = xy_vec[1] + originCoordinate.y;
		mypoint2f lastpoint_cb(0,0);
		mypoint2f lastpoint_qb(0,0);
		// find the next graph character after 'm'
		for (unsigned int i = 1; i < subElement.size(); i++)
		{
			char gtype = subElement[i].first;
			int gindex = subElement[i].second;
			string subPointstr;
			if ((i + 1) >= subElement.size())
			{
				subPointstr = d_str.substr(gindex + 1);				
			}
			else
			{
				subPointstr = d_str.substr(gindex + 1, subElement[i + 1].second - gindex - 1);
			}
			subPointstr = trimSpace(subPointstr);
			vector<double> subpvec = strToNum(subPointstr, 0, subPointstr.size() - 1);
			if (gtype == 'c' || gtype == 'C')
				lastpoint_cb = path_CubicBezier(startpoint, subpvec, gtype);
			else if (gtype == 's' || gtype == 'S')
				path_SmoothCubBezier(startpoint, lastpoint_cb, subpvec, gtype);
			else if (gtype == 'q' || gtype == 'Q')
				lastpoint_qb = path_QuadraticBezier(startpoint, subpvec, gtype);
			else if (gtype == 't' || gtype == 'T')
				path_SmoothQuaBezier(startpoint, lastpoint_qb, subpvec, gtype);
			else if (gtype == 'l' || gtype == 'L')
				path_LineTo(startpoint, subpvec, gtype);
			else if (gtype == 'h' || gtype == 'H')
				path_HorizontalTo(startpoint, subpvec, gtype);
			else if (gtype == 'v' || gtype == 'V')
				path_VerticalTo(startpoint, subpvec, gtype);			
		}
	}
}
string MySvgFile::trimSpace(string str)
{
	if (!str.empty())
	{
		str.erase(0, str.find_first_not_of(" "));
		str.erase(str.find_last_not_of(" ") + 1);
		str.erase(0, str.find_first_not_of('\t'));
		str.erase(str.find_last_not_of('\t') + 1);
		str.erase(0, str.find_first_not_of('\n'));
		str.erase(str.find_last_not_of('\n') + 1);
    }
	return str;
}
vector<double> MySvgFile::strToNum(string str, int startp, int endp)//example:str="12,33 44,22 5,2 34,11" strToNum(str,0,str.size()-1)
{
	vector<double> number_vec;
	double temp_num = 0;
	int current_index = startp-1;
	int next_index = 0;
	int i=startp;
	for (i; i <= endp; i++)
	{
		if (str.at(i) == ' ' || str.at(i) == ',')
		{
			next_index = i;
			string midstr = str.substr(current_index + 1, next_index - current_index-1);
			int midstr_size = next_index - current_index - 1;
			int curr_nega_pos = 0;
			int next_nega_pos = 0;
			for (int j = 0; j <midstr_size; ++j)
			{		
				if (midstr.at(j) == '-' && j!=0)
				{
					next_nega_pos = j;
					temp_num = atof(midstr.substr(curr_nega_pos, next_nega_pos - curr_nega_pos).c_str());
					number_vec.push_back(temp_num);
					curr_nega_pos = next_nega_pos;
				}
			}
			if (curr_nega_pos < midstr_size)
			{
				temp_num = atof(midstr.substr(curr_nega_pos, midstr_size - curr_nega_pos).c_str());//put the last number into vector;
				number_vec.push_back(temp_num);
			}
			current_index = next_index;
		}
	}
	if (next_index < endp)
	{
		string midstr = str.substr(current_index + 1, endp - current_index);
		int midstr_size = endp - current_index;
		int curr_nega_pos = 0;
		int next_nega_pos = 0;
		for (int j = 0; j < midstr_size; ++j)
		{
			if (midstr.at(j) == '-' && j != 0)
			{
				next_nega_pos = j;
				temp_num = atof(midstr.substr(curr_nega_pos, next_nega_pos - curr_nega_pos).c_str());
				number_vec.push_back(temp_num);
				curr_nega_pos = next_nega_pos;
			}
		}
		if (curr_nega_pos < midstr_size)
		{
			temp_num = atof(midstr.substr(curr_nega_pos, midstr_size - curr_nega_pos).c_str());//put the last number into vector;
			number_vec.push_back(temp_num);
		}
	}
	//if (str.at(startp) == '-')
	//{
	//	i = startp + 1;
	//	current_index = startp;
	//}
	//else
	//{
	//	i = startp;
	//	current_index = startp - 1;
	//}
	//for (i ; i <= endp; i++)
	//{
	//	if (str.at(i) == ' ' || str.at(i) == ',' || str.at(i) == '-')
	//	{
	//		next_index = i;
	//		if (current_index>= 0 && str.at(current_index) == '-')
	//		{			
	//			temp_num = atof(str.substr(current_index , next_index - current_index).c_str());
	//			number_vec.push_back(temp_num);
	//		}
	//		else
	//		{
	//			temp_num = atof(str.substr(current_index + 1, next_index - current_index - 1).c_str());
	//			number_vec.push_back(temp_num);
	//		}	
	//		current_index = next_index;
	//	}
	//}
	//if (next_index < endp)//check if it is over the border
	//{
	//	if (str.at(next_index) == '-')
	//	{//temp_num = atof(str.substr(next_index).c_str());
	//		temp_num = atof(str.substr(next_index, endp - next_index+1).c_str());//put the last number into vector;
	//		number_vec.push_back(temp_num);
	//	}
	//	else
	//	{//temp_num = atof(str.substr(next_index+1).c_str());
	//		temp_num = atof(str.substr(next_index + 1, endp - next_index).c_str());//put the last number into vector;
	//		number_vec.push_back(temp_num);
	//	}
	//}
	//else
	//{
	//	cout << "the last seperator is at the end of the string " << endl;
	//}	
	return number_vec;
}

mypoint2f MySvgFile::path_CubicBezier(mypoint2f &startp, vector<double> number_vec, char type)
{
	vector<mypoint2f> cubicBezierList;
	cubicBezierList.push_back(startp);
	mypoint2f apoint;
	mypoint2f temp_startp = startp;
	if (type == 'C')//absolute distance
	{
		for (unsigned int i = 0; i < number_vec.size(); i++)
		{
			apoint.x = number_vec[i] ;//apoint.x = number_vec[i];
			apoint.y = number_vec[i + 1];//apoint.y = number_vec[i + 1];
			i = i + 1;
			cubicBezierList.push_back(apoint);
		}
	}
	else//relative distance  
	{
		unsigned int i;
		for ( i = 0; i < number_vec.size(); i++)
		{
			if (i % 6 == 0 && i != 0)
			{
				temp_startp.y = temp_startp.y +number_vec[i - 1];
				temp_startp.x = temp_startp.x + number_vec[i - 2];
			}
			apoint.x = number_vec[i] +temp_startp.x;
			apoint.y = number_vec[i + 1] +temp_startp.y;
			i = i + 1;
			cubicBezierList.push_back(apoint);
		}
	}
	startp = cubicBezierList.back();//修改起始坐标
	cubicBezierList = transOrigiCooridate(cubicBezierList);//坐标transform
	for (unsigned int j = 0; j < cubicBezierList.size(); j += 3)//分段
	{
		cubicBezier acubicBezier;
		acubicBezier.p0 = cubicBezierList[j];
		acubicBezier.p1 = cubicBezierList[j + 1];
		acubicBezier.p2 = cubicBezierList[j + 2];
		acubicBezier.p3 = cubicBezierList[j + 3];
		acubicBezier.p0.y = svg_height - acubicBezier.p0.y;
		acubicBezier.p1.y = svg_height - acubicBezier.p1.y;
		acubicBezier.p2.y = svg_height - acubicBezier.p2.y;
		acubicBezier.p3.y = svg_height - acubicBezier.p3.y;

		cubicBezier_vec_origin.push_back(acubicBezier);
		if (j + 4 >= cubicBezierList.size())
		{
			break;
		}
	}
	return cubicBezierList.at(cubicBezierList.size() - 2);
}
mypoint2f MySvgFile::path_QuadraticBezier(mypoint2f &startp, vector<double> number_vec, char type)
{
	vector<mypoint2f> quadBezierList;
	quadBezierList.push_back(startp);
	mypoint2f apoint;
	mypoint2f temp_startp = startp;
	if (type == 'Q')//absolute distance
	{
		for (unsigned int i = 0; i < number_vec.size(); i++)
		{
			apoint.x = number_vec[i];// +originCoordinate.x;
			apoint.y = number_vec[i + 1];// +originCoordinate.y;
			i = i + 1;
			quadBezierList.push_back(apoint);
		}
	}
	else//relative distance  
	{
		for (unsigned int i = 0; i < number_vec.size(); i++)
		{
			if (i % 4 == 0 && i != 0)
			{
				temp_startp.y = temp_startp.y + number_vec[i - 1];
				temp_startp.x = temp_startp.x + number_vec[i - 2];
			}
			apoint.x = number_vec[i] + temp_startp.x;
			apoint.y = number_vec[i + 1] + temp_startp.y;
			//apoint.x = number_vec[i] + startp.x;
			//apoint.y = number_vec[i + 1] + startp.y;
			i = i + 1;
			quadBezierList.push_back(apoint);
		}
	}
	startp = quadBezierList.back();//修改起始坐标
	quadBezierList = transOrigiCooridate(quadBezierList);
	for (unsigned int j = 0; j < quadBezierList.size(); j += 2)//分段
	{
		quaBezier aquaBezier;
		aquaBezier.p0 = quadBezierList[j];
		aquaBezier.p1 = quadBezierList[j + 1];
		aquaBezier.p2 = quadBezierList[j + 2];
		aquaBezier.p0.y = svg_height - aquaBezier.p0.y;
		aquaBezier.p1.y = svg_height - aquaBezier.p1.y;
		aquaBezier.p2.y = svg_height - aquaBezier.p2.y;

		quadraticBezier_vec_origin.push_back(aquaBezier);
		if (j + 3 >= quadBezierList.size())
		{
			break;
		}
	}
	return quadBezierList.at(quadBezierList.size() - 2);
}
void MySvgFile::path_SmoothCubBezier(mypoint2f &startp, mypoint2f lastp, vector<double> number_vec, char type)
{
	vector<mypoint2f> cubicBezierList;
	cubicBezierList.push_back(startp);
	mypoint2f apoint;
	mypoint2f temp_startp = startp;
	mypoint2f tmp_lastp = startp + (startp - lastp);
	if (type == 'S')//absolute distance
	{
		for (unsigned int i = 1; i < number_vec.size(); i++)
		{
			apoint.x = number_vec[i-1];//apoint.x = number_vec[i];
			apoint.y = number_vec[i];//apoint.y = number_vec[i + 1];
			cubicBezierList.push_back(apoint);
		}
	}
	else//relative distance  
	{
		unsigned int i;
		for (i = 0; i < number_vec.size(); i++)
		{
			if (i % 4 == 0 && i != 0)
			{
				temp_startp.y = temp_startp.y + number_vec[i - 1];
				temp_startp.x = temp_startp.x + number_vec[i - 2];
			}
			apoint.x = number_vec[i] + temp_startp.x;
			apoint.y = number_vec[i + 1] + temp_startp.y;
			i = i + 1;
			cubicBezierList.push_back(apoint);
		}
	}
	vector<mypoint2f>::iterator iter = cubicBezierList.begin()+1;
	for (iter; iter != cubicBezierList.end(); iter+=3)
	{
		iter=cubicBezierList.insert(iter, tmp_lastp);
		tmp_lastp = *(iter + 2) + (*(iter + 2) - *(iter + 1));
	}
	startp = cubicBezierList.back();//修改起始坐标
	cubicBezierList = transOrigiCooridate(cubicBezierList);//坐标transform
	for (unsigned int j = 0; j < cubicBezierList.size(); j += 3)//分段
	{
		cubicBezier acubicBezier;
		acubicBezier.p0 = cubicBezierList[j];
		acubicBezier.p1 = cubicBezierList[j + 1];
		acubicBezier.p2 = cubicBezierList[j + 2];
		acubicBezier.p3 = cubicBezierList[j + 3];
		acubicBezier.p0.y = svg_height - acubicBezier.p0.y;
		acubicBezier.p1.y = svg_height - acubicBezier.p1.y;
		acubicBezier.p2.y = svg_height - acubicBezier.p2.y;
		acubicBezier.p3.y = svg_height - acubicBezier.p3.y;

		cubicBezier_vec_origin.push_back(acubicBezier);
		if (j + 4 >= cubicBezierList.size())
		{
			break;
		}
	}
}
void MySvgFile::path_SmoothQuaBezier(mypoint2f &startp, mypoint2f lastp, vector<double> number_vec, char type)
{
	vector<mypoint2f> quadBezierList;
	quadBezierList.push_back(startp);
	mypoint2f apoint;
	mypoint2f temp_startp = startp;
	mypoint2f tmp_lastp = startp + (startp - lastp);
	if (type == 'Q')//absolute distance
	{
		for (unsigned int i = 0; i < number_vec.size(); i++)
		{
			apoint.x = number_vec[i];// +originCoordinate.x;
			apoint.y = number_vec[i + 1];// +originCoordinate.y;
			i = i + 1;
			quadBezierList.push_back(apoint);
		}
	}
	else//relative distance  
	{
		for (unsigned int i = 0; i < number_vec.size(); i++)
		{
			if (i % 2 == 0 && i != 0)
			{
				temp_startp.y = temp_startp.y + number_vec[i - 1];
				temp_startp.x = temp_startp.x + number_vec[i - 2];
			}
			apoint.x = number_vec[i] + temp_startp.x;
			apoint.y = number_vec[i + 1] + temp_startp.y;
			i = i + 1;
			quadBezierList.push_back(apoint);
		}
	}
	vector<mypoint2f>::iterator iter = quadBezierList.begin()+1;
	for (iter; iter != quadBezierList.end(); iter += 3)
	{
		iter = quadBezierList.insert(iter, tmp_lastp);
		tmp_lastp = *(iter + 2) + (*(iter + 2) - *(iter + 1));
	}
	startp = quadBezierList.back();//修改起始坐标
	quadBezierList = transOrigiCooridate(quadBezierList);
	for (unsigned int j = 0; j < quadBezierList.size(); j += 2)//分段
	{
		quaBezier aquaBezier;
		aquaBezier.p0 = quadBezierList[j];
		aquaBezier.p1 = quadBezierList[j + 1];
		aquaBezier.p2 = quadBezierList[j + 2];
		aquaBezier.p0.y = svg_height - aquaBezier.p0.y;
		aquaBezier.p1.y = svg_height - aquaBezier.p1.y;
		aquaBezier.p2.y = svg_height - aquaBezier.p2.y;

		quadraticBezier_vec_origin.push_back(aquaBezier);
		if (j + 3 >= quadBezierList.size())
		{
			break;
		}
	}
}
void MySvgFile::path_LineTo(mypoint2f &startp, vector<double> number_vec, char type)/////只适用于从noris文章中的输出文件 
{
	vector<mypoint2f> aPolyline;
	aPolyline.push_back(startp);
	mypoint2f apoint;
	if (type == 'L' || type == 'M')//absolute distance
	{
		for (unsigned int i = 0; i < number_vec.size(); i++)
		{
			apoint.x = number_vec[i];// +originCoordinate.x;
			apoint.y = number_vec[i + 1];// +originCoordinate.y;
			i = i + 1;
			aPolyline.push_back(apoint);
		}
	}
	else//relative distance  
	{
		mypoint2f lastp = startp;
		for (unsigned int i = 0; i < number_vec.size(); i++)
		{
			apoint.x = number_vec[i] + lastp.x;
			apoint.y = number_vec[i + 1] + lastp.y;
			i = i + 1;
			aPolyline.push_back(apoint);
			lastp = apoint;
		}
	}
	startp = aPolyline.back();

	aPolyline = transOrigiCooridate(aPolyline);	
	for (unsigned int i = 1; i < aPolyline.size(); ++i)//divide aPolyline into segments
	{
		PolyLine pl;
		pl.pointlist.push_back(aPolyline[i - 1]);
		pl.pointlist.back().y = svg_height - pl.pointlist.back().y;
		pl.pointlist.push_back(aPolyline[i]);
		pl.pointlist.back().y = svg_height - pl.pointlist.back().y;
		//vector<mypoint2f> tmpvec;
		//tmpvec.push_back(aPolyline[i - 1]);
		//tmpvec.push_back(aPolyline[i]);
		pl.origi_index = (int)polyline_vec_origin.size();
		pl.origi_type = 'l';
		polyline_vec_origin.push_back(pl);
	}
}
void MySvgFile::path_HorizontalTo(mypoint2f &startp, vector<double> number_vec, char type)
{
	vector<mypoint2f> aPolyline;
	aPolyline.push_back(startp);
	mypoint2f apoint;
	if (type == 'H')//absolute distance
	{
		for (unsigned int i = 0; i < number_vec.size(); i++)
		{
			apoint.x = number_vec[i] + startp.x;
			apoint.y = startp.y;
			aPolyline.push_back(apoint);
		}
	}
	else//relative distance  
	{
		for (unsigned int i = 0; i < number_vec.size(); i++)
		{
			apoint.x = number_vec[i] + aPolyline.back().x;
			apoint.y = aPolyline.back().y;
			aPolyline.push_back(apoint);
		}
	}
	startp = aPolyline.back();
	aPolyline = transOrigiCooridate(aPolyline);
	for (unsigned int i = 1; i < aPolyline.size(); ++i)
	{
		//vector<mypoint2f> tmpvec;
		//tmpvec.push_back(aPolyline[i - 1]);
		//tmpvec.push_back(aPolyline[i]);
		//polyline_vec.push_back(tmpvec);
		PolyLine pl;
		pl.pointlist.push_back(aPolyline[i - 1]);
		pl.pointlist.back().y = svg_height - pl.pointlist.back().y;
		pl.pointlist.push_back(aPolyline[i]);
		pl.pointlist.back().y = svg_height - pl.pointlist.back().y;
		pl.origi_index = (int)polyline_vec_origin.size();
		pl.origi_type = 'l';
		polyline_vec_origin.push_back(pl);
	}
}
void MySvgFile::path_VerticalTo(mypoint2f &startp, vector<double> number_vec, char type)
{
	vector<mypoint2f> aPolyline;
	aPolyline.push_back(startp);
	mypoint2f apoint;
	if (type == 'V')//absolute distance
	{
		for (unsigned int i = 0; i < number_vec.size(); i++)
		{
			apoint.y = number_vec[i] + startp.y;
			apoint.x = startp.x;
			aPolyline.push_back(apoint);
		}
	}
	else//relative distance  
	{
		for (unsigned int i = 0; i < number_vec.size(); i++)
		{
			apoint.y = number_vec[i] + aPolyline.back().y;
			apoint.x = aPolyline.back().x;
			aPolyline.push_back(apoint);
		}
	}
	startp = aPolyline.back();
	aPolyline = transOrigiCooridate(aPolyline);
	for (unsigned int i = 1; i < aPolyline.size(); ++i)
	{
		//vector<mypoint2f> tmpvec;
		//tmpvec.push_back(aPolyline[i - 1]);
		//tmpvec.push_back(aPolyline[i]);
		//polyline_vec.push_back(tmpvec);
		PolyLine pl;
		pl.pointlist.push_back(aPolyline[i - 1]);
		pl.pointlist.back().y = svg_height - pl.pointlist.back().y;
		pl.pointlist.push_back(aPolyline[i]);
		pl.pointlist.back().y = svg_height - pl.pointlist.back().y;

		pl.origi_index = (int)polyline_vec_origin.size();
		pl.origi_type = 'l';
		polyline_vec_origin.push_back(pl);
	}
}
void MySvgFile::shape_Line(TiXmlElement *graphelement)
{
	mypoint2f p1,p2;
	PolyLine pl;
	//vector<mypoint2f> aline;
	double temp_num=0;
	if (graphelement->QueryDoubleAttribute("x1", &temp_num) == TIXML_SUCCESS)
		p1.x = temp_num;
	else
		throw temp_num;
	if (graphelement->QueryDoubleAttribute("y1", &temp_num) == TIXML_SUCCESS)
		p1.y = svg_height - temp_num;
	else
		throw temp_num;
	if (graphelement->QueryDoubleAttribute("x2", &temp_num) == TIXML_SUCCESS)
		p2.x = temp_num;
	else
		throw temp_num;
	if (graphelement->QueryDoubleAttribute("y2", &temp_num) == TIXML_SUCCESS)
		p2.y = svg_height - temp_num;
	else
		throw temp_num;
	p1 = transOrigiCooridate(p1);
	p2 = transOrigiCooridate(p2);
	pl.pointlist.push_back(p1);
	pl.pointlist.push_back(p2);
	pl.origi_index = (int)polyline_vec_origin.size();
	pl.origi_type = 'l';
	polyline_vec_origin.push_back(pl);
	//aline.push_back(p1);
	//aline.push_back(p2);
	//polyline_vec.push_back(aline);
}
void MySvgFile::shape_Circle(TiXmlElement *graphelement)
{
	mycircle acircle;
	mypoint2f center_circle;
	double temp_num=0;
	if (graphelement->QueryDoubleAttribute("cx", &temp_num) == TIXML_SUCCESS)
		center_circle.x = temp_num;// +originCoordinate.x;
	else
		throw temp_num;
	if (graphelement->QueryDoubleAttribute("cy", &temp_num) == TIXML_SUCCESS)
		center_circle.y = temp_num;
	else
		throw temp_num;
	if (graphelement->QueryDoubleAttribute("r", &temp_num) == TIXML_SUCCESS)
		acircle.radius = temp_num;
	else
		throw temp_num;
	center_circle = transOrigiCooridate(center_circle);
	acircle.cx = center_circle.x;
	acircle.cy = svg_height - center_circle.y;
	circle_vec.push_back(acircle);
}
void MySvgFile::shape_Rect(TiXmlElement *graphelement)
{
	myrect arect;
	mypoint2f rec_xy;
	double temp_num=0;
	if (graphelement->QueryDoubleAttribute("x", &temp_num) == TIXML_SUCCESS)
		rec_xy.x = temp_num;// +originCoordinate.x;
	else
		throw temp_num;
	if (graphelement->QueryDoubleAttribute("y", &temp_num) == TIXML_SUCCESS)
		rec_xy.y = temp_num;// +originCoordinate.y;
	else
		throw temp_num;
	if (graphelement->QueryDoubleAttribute("width", &temp_num) == TIXML_SUCCESS)
		arect.width = temp_num;
	else
		throw temp_num;
	if (graphelement->QueryDoubleAttribute("height", &temp_num) == TIXML_SUCCESS)
		arect.height = temp_num;
	else
		throw temp_num;
	rec_xy = transOrigiCooridate(rec_xy);
	arect.x = rec_xy.x;
	arect.y = svg_height - rec_xy.y;
	rect_vec.push_back(arect);
}
void MySvgFile::shape_Polyline(TiXmlElement *graphelement)
{
	const char *c_pointset=graphelement->Attribute("points");
	string str_pointset = c_pointset;
	vector<double> xy_vec = strToNum(str_pointset, 0, str_pointset.size() - 1);

	vector<mypoint2f> a_polyline_vec;
	mypoint2f apoint;
	for (unsigned int i = 1; i < xy_vec.size(); i++)
	{
		apoint.x = xy_vec.at(i - 1);// +originCoordinate.x;
		apoint.y = xy_vec.at(i);// +originCoordinate.y;
		i = i + 1;
		a_polyline_vec.push_back(apoint);
	}
	a_polyline_vec = transOrigiCooridate(a_polyline_vec);
	for (unsigned int j = 1; j < a_polyline_vec.size(); ++j)//cut the long polyline into small line segments
	{
		//vector<mypoint2f> tmpvec;
		//tmpvec.push_back(a_polyline_vec[j - 1]);
		//tmpvec.push_back(a_polyline_vec[j]);
		//polyline_vec.push_back(tmpvec);
		PolyLine pl;
		pl.pointlist.push_back(a_polyline_vec[j - 1]);
		pl.pointlist.back().y = svg_height - pl.pointlist.back().y;
		pl.pointlist.push_back(a_polyline_vec[j]);
		pl.pointlist.back().y = svg_height - pl.pointlist.back().y;
		pl.origi_index = (int)polyline_vec_origin.size();
		pl.origi_type = 'l';
		polyline_vec_origin.push_back(pl);
	}
}
void MySvgFile::shape_Polygon(TiXmlElement *graphelement)//similar with shape_Polyline()
{
	const char *c_pointset = graphelement->Attribute("points");
	string str_pointset = c_pointset;
	vector<double> xy_vec = strToNum(str_pointset, 0, str_pointset.size() - 1);
	double temp_num = 0;
	//put the first two number into vector again
	temp_num = xy_vec.at(0);
	xy_vec.push_back(temp_num);
	temp_num = xy_vec.at(1);
	xy_vec.push_back(temp_num);

	vector<mypoint2f> a_polygon_vec;
	mypoint2f apoint;
	for (unsigned int i = 0; i < xy_vec.size(); i++)
	{
		apoint.x = xy_vec.at(i);// +originCoordinate.x;
		apoint.y = xy_vec.at(i + 1);
		i = i + 1;
		a_polygon_vec.push_back(apoint);
	}
	a_polygon_vec = transOrigiCooridate(a_polygon_vec);
	for (unsigned int i = 0; i < a_polygon_vec.size(); ++i)
	{
		a_polygon_vec[i].y = svg_height - a_polygon_vec[i].y;
	}
	polygon_vec.push_back(a_polygon_vec);
}

//about opengl
void MySvgFile::testShowDiffvec()
{
	cout << "quadraticBezier:" << quadraticBezier_vec.size() << endl;
	//for (unsigned int i = 0; i < quadraticBezier_vec.size(); i++)
	//{
	//	cout << quadraticBezier_vec[i].p0.x << "," << quadraticBezier_vec[i].p0.y << " ";
	//	cout << quadraticBezier_vec[i].p1.x << "," << quadraticBezier_vec[i].p1.y << " ";
	//	cout << quadraticBezier_vec[i].p2.x << "," << quadraticBezier_vec[i].p2.y << " ";
	//	cout << endl;
	//}
	cout << "cubicBezier_vec:" << cubicBezier_vec.size() << endl;
	//for (unsigned int i = 0; i < cubicBezier_vec.size(); i++)
	//{
	//	cout << cubicBezier_vec[i].p0.x << "," << cubicBezier_vec[i].p0.y << " ";
	//	cout << cubicBezier_vec[i].p1.x << "," << cubicBezier_vec[i].p1.y << " ";
	//	cout << cubicBezier_vec[i].p2.x << "," << cubicBezier_vec[i].p2.y << " ";
	//	cout << cubicBezier_vec[i].p3.x << "," << cubicBezier_vec[i].p3.y << " ";
	//	cout << endl;
	//}
	cout << "polyline_vec:" << polyline_vec.size() << endl;
	/*for (unsigned int i = 0; i < polyline_vec.size(); i++)
	{
		for (unsigned int j = 0; j < polyline_vec[i].pointlist.size(); j++)
		{
			cout << polyline_vec[i].pointlist.at(j).x << "," << polyline_vec[i].pointlist.at(j).y << " ";
		}
		cout << endl;
	}*/
	cout << "circle_vec:" << circle_vec.size() << endl;
	//for (unsigned int i = 0; i < circle_vec.size(); i++)
	//{
	//	cout << circle_vec[i].cx << "," << circle_vec[i].cy << " ,r="<<circle_vec[i].radius << endl;;
	//	cout << endl;
	//}
	setInstance();
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	glutInitWindowPosition(50, 100);
	glutInitWindowSize(winWidth, winHeight);
	glutCreateWindow("svg_controlPoint");
	init();
	glutDisplayFunc(mydisplay);
	glutReshapeFunc(winreshape);
	glutSpecialFunc(myspecialkey);
	glutMainLoop();

}
void MySvgFile::setInstance()
{
	cout << "MySvgFile::instance = this;" << endl;
	MySvgFile::svginstance = this;
}
void MySvgFile::init()
{
	glViewport(0, 0, winWidth, winHeight);
	glClearColor(1.0, 1.0, 1.0, 0.0);
	glClear(GL_COLOR_BUFFER_BIT);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(w_x1, w_x2, w_y1, w_y2);
	//gluPerspective(100, (GLfloat)winWidth / (GLfloat)winHeight, 1, 1000.0);
	//gluLookAt(winWidth / 2, winHeight / 2, 100.0, winWidth / 2, winHeight / 2, 0.0, 0.0, 1.0, 0.0);
	//gluLookAt(90.0, 90.0, 100.0, 90, 90.0, 0.0, 0.0, 1.0, 0.0);
}
void MySvgFile::winreshape(GLint newWidth, GLint newHeight)
{
	winWidth = newWidth;
	winHeight = newHeight;
	w_x1 = 0;
	w_x2 = winWidth;
	w_y1 = 0;
	w_y2 = winHeight;
	glViewport(0, 0, newWidth, newHeight);//glViewport(0, 0, (GLsizei) w, (GLsizei) h)注释掉此句窗口中的图形将不会随着窗口的拉伸而发生形变,因为视口的width和height没有改变
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	//gluPerspective(100, (GLdouble)newWidth / (GLdouble)newHeight, 10.0, 1000.0);
	//gluLookAt(winWidth / 2, winHeight / 2, 100.0, winWidth / 2, winHeight / 2, 0.0, 0.0, 1.0, 0.0);
	//gluLookAt(newWidth / 2, newHeight / 2, 100.0, newWidth / 2, newHeight / 2, 0.0, 0.0, 1.0, 0.0);
	//gluLookAt(90.0, 90.0, 100.0, 90, 90.0, 0.0, 0.0, 1.0, 0.0);
	gluOrtho2D(w_x1, w_x2, w_y1, w_y2);
}
void MySvgFile::myspecialkey(int key, int x, int y)
{
	if (key == GLUT_KEY_UP)
	{
		float keyx = (w_x2 - w_x1) / winWidth;
		float keyy = (w_y2 - w_y1) / winHeight;
		float newcenterx = x*keyx + w_x1;//实际的 在世界坐标系中的坐标
		float newcentery = (winHeight - y)*keyy + w_y1;
		//cout << "real x=" << newcenterx << "real y=" << newcentery << endl;
		float win_wid = w_x2 - w_x1;
		float win_hei = w_y2 - w_y1;
		w_x1 = newcenterx - win_wid / 2 * zoomup;//以鼠标所在点为中心，构建一个新的正交投影窗口
		w_x2 = newcenterx + win_wid / 2 * zoomup;
		w_y1 = newcentery - win_hei / 2 * zoomup;
		w_y2 = newcentery + win_hei / 2 * zoomup;
		//cout << w_x1 << "," << w_x2 << "," << w_y1 << "," << w_y2 << endl << endl;
		//glViewport(0, 0, winWidth, winHeight);
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		gluOrtho2D(w_x1, w_x2, w_y1, w_y2);
		glutPostRedisplay();
	}
	if (key == GLUT_KEY_DOWN)
	{
		float keyx = (w_x2 - w_x1) / winWidth;
		float keyy = (w_y2 - w_y1) / winHeight;
		float newcenterx = x*keyx + w_x1;
		float newcentery = (winHeight - y)*keyy + w_y1;
		//cout << "real x=" << newcenterx << "real y=" << newcentery << endl;
		float win_wid = w_x2 - w_x1;
		float win_hei = w_y2 - w_y1;
		w_x1 = newcenterx - win_wid / 2 * zoomdown;
		w_x2 = newcenterx + win_wid / 2 * zoomdown;
		w_y1 = newcentery - win_hei / 2 * zoomdown;
		w_y2 = newcentery + win_hei / 2 * zoomdown;
		//cout << w_x1 << "," << w_x2 << "," << w_y1 << "," << w_y2 << endl << endl;
		//glViewport(0, 0, winWidth, winHeight);
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		gluOrtho2D(w_x1, w_x2, w_y1, w_y2);
		glutPostRedisplay();
	}
}

void MySvgFile::mydisplay(void)
{
	glClearColor(1.0, 1.0, 1.0, 0.0);
	glClear(GL_COLOR_BUFFER_BIT); 

	glLineWidth(2);//设置线段宽度
	glBegin(GL_LINES);
	glVertex2f(-10.0f, 0.0f);
	glVertex2f(100.0f, 0.0f);
	glVertex2f(0.0f, -10.0f);
	glVertex2f(0.0f, 100.0f);
	glEnd();

	//paintCtlConvex();//blue or green
	glColor4f(0.0, 0.0, 0.0, 0.5);
	paintCtlPoint(svginstance->circle_vec); 
	paintCubicBezierCurve(svginstance->cubicBezier_vec_origin);//black
	paintQuaBezierCurve(svginstance->quadraticBezier_vec);//black
	glColor4f(0.0, 0.0, 0.0, 0.5);
	paintPolyline(svginstance->polyline_vec_origin);//black

	glColor4f(1.0, 0.0, 0.0, 0.5);
	paintCubicBezierCurve(svginstance->piecewise_Bezier);
	//paintPolyline(svginstance->polyline_vec_origin);//black
	//test_paintPCABezierCurve(svginstance->cubicBezier_vec_origin);  //red+blue
	//glColor4f(1.0, 0.0, 0.0, 0.5);
	//paintPolyline(svginstance->polyline_extend);

	////画出这一组散乱点的有序连接
	//GLfloat view_y = 1052.0;
	//glBegin(GL_LINE_STRIP);
	//for (unsigned int i = 0; i < svginstance->ordered_plist.size(); ++i)
	//{
	//	GLfloat px = (GLfloat)svginstance->ordered_plist.at(i).x;
	//	GLfloat py = (GLfloat)(view_y - svginstance->ordered_plist.at(i).y);
	//	glVertex2f(px, py);
	//}
	//glEnd();

	glFlush();
}
void MySvgFile::paintCtlConvex()
{
	GLfloat px = 0; GLfloat py = 0;
	//GLfloat view_y = 1052.0;
	//quadraticBezier and cubicBezier control polygon blue
	for (unsigned int i = 0; i < svginstance->quadraticBezier_vec.size(); i++)
	{
		glBegin(GL_LINE_STRIP);
		//glColor3f(0.0, 0.0, 1.0);
		glColor4f(0.5, 0.5, 1.0, 1.0);//light blue
		px = (GLfloat)svginstance->quadraticBezier_vec[i].p0.x;
		py = (GLfloat)svginstance->quadraticBezier_vec[i].p0.y;
		glVertex2f(px, py);
		px = (GLfloat)svginstance->quadraticBezier_vec[i].p1.x;
		py = (GLfloat)svginstance->quadraticBezier_vec[i].p1.y;
		glVertex2f(px, py);
		px = (GLfloat)svginstance->quadraticBezier_vec[i].p2.x;
		py = (GLfloat)svginstance->quadraticBezier_vec[i].p2.y;
		glVertex2f(px, py);
		glEnd();
	}
	for (unsigned int i = 0; i < svginstance->cubicBezier_vec.size(); i++)
	{
		glBegin(GL_LINE_STRIP);
		//glColor3f(0.0, 0.0, 1.0);
		glColor4f(0.5, 0.5, 1.0, 1.0);//light blue
		px = (GLfloat)svginstance->cubicBezier_vec[i].p0.x;
		py = (GLfloat)svginstance->cubicBezier_vec[i].p0.y;
		glVertex2f(px, py);
		px = (GLfloat)svginstance->cubicBezier_vec[i].p1.x;
		py = (GLfloat)svginstance->cubicBezier_vec[i].p1.y;
		glVertex2f(px, py);
		px = (GLfloat)svginstance->cubicBezier_vec[i].p2.x;
		py = (GLfloat)svginstance->cubicBezier_vec[i].p2.y;
		glVertex2f(px, py);
		px = (GLfloat)svginstance->cubicBezier_vec[i].p3.x;
		py = (GLfloat)svginstance->cubicBezier_vec[i].p3.y;
		glVertex2f(px, py);
		glEnd();
	}
}
void MySvgFile::paintCtlPoint(vector<mycircle> circles)
{
	//circle  black
	//GLfloat view_y = 1052.0;
	GLdouble pi = 3.1415926;
	for (unsigned int i = 0; i < circles.size(); ++i)
	{
		glPointSize(3.0);
		glColor3f(0.0, 0.0, 0.0);
		glBegin(GL_POINTS);
		glVertex2f(circles[i].cx, circles[i].cy);
		glEnd();

		glColor3f(0.0, 0.0, 0.0);
		glBegin(GL_LINE_LOOP);// GL_POLYGON	
		GLfloat c_x = circles[i].cx;
		GLfloat c_y = circles[i].cy;
		for (unsigned int j = 0; j < 10; ++j)
		{
			glVertex2f(c_x + circles[i].radius*cos(2 * pi * j / 10), c_y + circles[i].radius*sin(2 * pi  * j / 10));
		}
		glEnd();
	}
}
void MySvgFile::paintCubicBezierCurve(vector<cubicBezier> cbcurve)//红色
{
	for (unsigned int j = 0; j < cbcurve.size(); ++j)  //red red red
	{
		GLfloat p1[2] = { cbcurve[j].p0.x, cbcurve[j].p0.y };
		GLfloat p2[2] = { cbcurve[j].p1.x, cbcurve[j].p1.y };
		GLfloat p3[2] = { cbcurve[j].p2.x, cbcurve[j].p2.y };
		GLfloat p4[2] = { cbcurve[j].p3.x, cbcurve[j].p3.y };
		GLfloat pcurve[51][2];
		GLint i = 0;
		GLfloat step = 1 / 50;//t[i]均匀的
		//三次公式 B(t)=p0*(1-t)^3 + p1*3t(1-t)^2 + p2*3t^2(1-t) + p3*t^3
		for (GLfloat t = 0.0; t <= 1.0; t += 0.02)
		{
			GLfloat a1 = pow((1 - t), 3);
			GLfloat a2 = pow((1 - t), 2) * 3 * t;
			GLfloat a3 = 3 * t*t*(1 - t);
			GLfloat a4 = t*t*t;
			pcurve[i][0] = a1*p1[0] + a2*p2[0] + a3*p3[0] + a4*p4[0];
			pcurve[i][1] = a1*p1[1] + a2*p2[1] + a3*p3[1] + a4*p4[1];
			i = i + 1;
		}
		//glColor4f(0.0, 0.0, 0.0, 1.0);
		glBegin(GL_LINE_STRIP);
		for (int i = 0; i < 51; i++)
		{
			glVertex2fv(pcurve[i]);
		}
		glEnd();
	}
}
void MySvgFile::paintQuaBezierCurve(vector<quaBezier> qbcurve)//red色
{
	for (unsigned int j = 0; j < qbcurve.size(); ++j)
	{
		GLfloat p1[2] = { qbcurve[j].p0.x, qbcurve[j].p0.y };
		GLfloat p2[2] = { qbcurve[j].p1.x, qbcurve[j].p1.y };
		GLfloat p3[2] = { qbcurve[j].p2.x, qbcurve[j].p2.y };
		GLfloat pcurve[51][2];
		GLint i = 0;
		//GLdouble step = 1 / 50;
		//二次贝塞尔 B(t)=p0*(1-t)^2 + p1*2t(1-t) + p2*t^2
		for (GLfloat t = 0.0; t <= 1.0; t += 0.02)
		{
			GLfloat a1 = pow((1 - t), 2);
			GLfloat a2 = 2 * t*(1 - t);
			GLfloat a3 = t*t;
			pcurve[i][0] = a1*p1[0] + a2*p2[0] + a3*p3[0];
			pcurve[i][1] = a1*p1[1] + a2*p2[1] + a3*p3[1];
			i = i + 1;
		}
		glColor4f(0.0, 0.0, 0.0, 1.0);
		glBegin(GL_LINE_STRIP);
		for (int i = 0; i < 51; i++)
		{
			glVertex2fv(pcurve[i]);
		}
		glEnd();
	}
}
void MySvgFile::paintPolyline(vector<PolyLine> pline)
{
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
	};
	GLfloat px = 0; GLfloat py = 0;
	//GLfloat view_y = 1052.0;
	for (unsigned int i = 0; i < pline.size(); i++)
	{
		//glColor4fv(colors[(i % 12)]);//glColor4f(0.0, 0.0, 0.0, 0.5);
		glBegin(GL_LINE_STRIP);
		for (unsigned int j = 0; j < pline[i].pointlist.size(); j++)
		{
			px = (GLfloat)pline[i].pointlist.at(j).x;
			py = (GLfloat)pline[i].pointlist.at(j).y;
			glVertex2f(px, py);
		}
		glEnd();
	}

	////画出前n-1个线是黑色的，最后一条innerline是红色的
	//for (unsigned int i = 0; i < pline.size()-1; i++)
	//{
	//	glColor4f(0.0, 0.0, 0.0, 0.5);
	//	glBegin(GL_LINE_STRIP);	
	//	for (unsigned int j = 0; j < pline[i].pointlist.size(); j++)
	//	{
	//		px = (GLfloat)pline[i].pointlist.at(j).x;
	//		py = (GLfloat)(view_y - pline[i].pointlist.at(j).y);
	//		glVertex2f(px, py);
	//	}
	//	glEnd();
	//}
	//glColor4f(0.0, 1.0, 0.0, 0.5);
	//glBegin(GL_LINE_STRIP);
	//for (unsigned int j = 0; j < pline.back().pointlist.size(); j++)
	//{
	//	px = (GLfloat)pline.back().pointlist.at(j).x;
	//	py = (GLfloat)(view_y - pline.back().pointlist.at(j).y);
	//	glVertex2f(px, py);
	//}
	//glEnd();

	////旋转公式
	//vector<mypoint2f> aline;
	//aline.push_back(mypoint2f(1.0, 2.0));
	//aline.push_back(mypoint2f(2.0, 2.0));
	//aline.push_back(mypoint2f(3.0, 2.0));
	//////x'=(x+rtx)costheta + (y+rty)sintheta-rtx;
	//////y'=(x+rtx)sintheta + (y+rty)costheta-rty;
	//////rotate polyline red
	//GLfloat pi = 3.1415;
	//GLfloat px2 = 0.0; GLfloat py2 = 0.0;
	//GLfloat px_rotate = 0.0; GLfloat py_rotate = 0.0;
	//GLfloat rtx = (-1)*(GLfloat)aline.at(0).x;
	//GLfloat rty = (-1)*(GLfloat)aline.at(0).y;
	//glColor4f(0.0, 0.0, 1.0, 0.5);
	//glBegin(GL_LINE_STRIP);
	////glColor3f(0.0, 1.0, 0.0);		
	//for (unsigned int j = 0; j < aline.size(); j++)
	//{
	//	px2 = (GLfloat)aline.at(j).x;
	//	py2 = (GLfloat)(aline.at(j).y);
	//	px_rotate = (px2+rtx)*cos(-pi / 6) - (py2+rty)*sin(-pi / 6)-rtx;
	//	py_rotate = (py2 + rty)*cos(-pi / 6) + (px2 + rtx)*sin(-pi / 6)-rty;
	//	glVertex2f(px_rotate, py_rotate);
	//}
	//glEnd();
	
}

void MySvgFile::test_paintPCABezierCurve(vector<cubicBezier> cbcurve)//原来的线都是红色，最后一条新的beziercurve 是蓝色的
{
	for (unsigned int j = 0; j < cbcurve.size(); ++j)  //red red red
	{
		GLfloat p1[2] = { cbcurve[j].p0.x, cbcurve[j].p0.y };
		GLfloat p2[2] = { cbcurve[j].p1.x, cbcurve[j].p1.y };
		GLfloat p3[2] = { cbcurve[j].p2.x, cbcurve[j].p2.y };
		GLfloat p4[2] = { cbcurve[j].p3.x, cbcurve[j].p3.y };
		GLfloat pcurve[51][2];
		GLint i = 0;
		GLfloat step = 1 / 50;//t[i]均匀的
		//三次公式 B(t)=p0*(1-t)^3 + p1*3t(1-t)^2 + p2*3t^2(1-t) + p3*t^3
		for (GLfloat t = 0.0; t <= 1.0; t += 0.02)
		{
			GLfloat a1 = pow((1 - t), 3);
			GLfloat a2 = pow((1 - t), 2) * 3 * t;
			GLfloat a3 = 3 * t*t*(1 - t);
			GLfloat a4 = t*t*t;
			pcurve[i][0] = a1*p1[0] + a2*p2[0] + a3*p3[0] + a4*p4[0];
			pcurve[i][1] = a1*p1[1] + a2*p2[1] + a3*p3[1] + a4*p4[1];
			i = i + 1;
		}
		glColor4f(1.0, 0.0, 0.0, 1.0);
		glBegin(GL_LINE_STRIP);
		for (int i = 0; i < 51; i++)
		{
			glVertex2fv(pcurve[i]);
		}
		glEnd();
	}

	//	画出最后一条拟合的曲线  blue blue blue
	GLfloat p1[2] = { cbcurve.back().p0.x, cbcurve.back().p0.y };
	GLfloat p2[2] = { cbcurve.back().p1.x, cbcurve.back().p1.y };
	GLfloat p3[2] = { cbcurve.back().p2.x, cbcurve.back().p2.y };
	GLfloat p4[2] = { cbcurve.back().p3.x, cbcurve.back().p3.y };
	GLfloat pcurve[51][2];
	GLint i = 0;
	GLfloat step = 1 / 50;//t[i]均匀的
	//三次公式 B(t)=p0*(1-t)^3 + p1*3t(1-t)^2 + p2*3t^2(1-t) + p3*t^3
	for (GLfloat t = 0.0; t <= 1.0; t += 0.02)
	{
		GLfloat a1 = pow((1 - t), 3);
		GLfloat a2 = pow((1 - t), 2) * 3 * t;
		GLfloat a3 = 3 * t*t*(1 - t);
		GLfloat a4 = t*t*t;
		pcurve[i][0] = a1*p1[0] + a2*p2[0] + a3*p3[0] + a4*p4[0];
		pcurve[i][1] = a1*p1[1] + a2*p2[1] + a3*p3[1] + a4*p4[1];
		i = i + 1;
	}
	glColor4f(0.0, 0.0, 1.0, 1.0);
	glBegin(GL_LINE_STRIP);
	for (int i = 0; i < 51; i++)
	{
		glVertex2fv(pcurve[i]);
	}
	glEnd();
}



void MySvgFile::getCrossCurve()
{
	cout << "getCrossCurve" << endl;
	int cbsize = (int)cubicBezier_vec_origin.size();
	int qbsize = (int)quadraticBezier_vec_origin.size();
	int plsize = (int)polyline_vec_origin.size();
	int diffcursize = cbsize + qbsize + plsize;
	vector<DiffCurve> diffcurList;// = new DiffCurve[diffcursize];
	for (int i = 0; i < cbsize; ++i)
	{
		DiffCurve dtcur;
		dtcur.curveIndex = i;
		dtcur.curveType = 'c';
		dtcur.index = i;
		diffcurList.push_back(dtcur);
	}
	for (int i = 0; i < qbsize; ++i)
	{
		DiffCurve dtcur;
		dtcur.curveIndex = i;
		dtcur.curveType = 'q';
		dtcur.index = i + cbsize;
		diffcurList.push_back(dtcur);
	}
	int tmp_addNum = diffcurList.size(); //=cbsize + qbsize
	for (int i = 0; i < plsize; ++i)
	{
		DiffCurve dtcur;
		dtcur.curveIndex = i;
		dtcur.curveType = 'l';
		dtcur.index = i + tmp_addNum; //= i + cbsize + qbsize;
		diffcurList.push_back(dtcur);
	}

	diffcursize = diffcurList.size();
	vector<bool> isSplitedCurve(diffcursize, false);//因为isSplitedCurve是动态变化的 所以用vector
	/*bool *isSplitedCurve = new bool[diffcursize];
	for (int i = 0; i < diffcursize; ++i)
	{
		isSplitedCurve[i] = false;
	}*/
	for (int i = 0; i < diffcursize; ++i)
	{
		for (int j = 0; j < diffcursize; ++j)
		{
			if (i != j)//if (i != j && !isAdjOnTerminate(diffcurList[i], diffcurList[j]))//exclude the condition where two curves adjacent by a common point
			{
				coarseJudgeByConvex(i, j, &diffcurList);
			}
		}
	}
	getCrossCurvePosi(&diffcurList, &isSplitedCurve);//normal condition
	diffcursize = diffcurList.size();//！！diffcurList添加了新的extended line ，
	for (int i = 0; i < diffcursize; ++i)
	{
		if (!diffcurList[i].split_t.empty() || !diffcurList[i].split_point.empty())
		{
			splitACurve(diffcurList[i], diffcurList[i].split_t, diffcurList[i].split_point,1);//最后一个参数0表示每段曲线参数t=1/50=0.02. 1表示俺曲线大概长度去不同的t值
		}
	}

	deleOrigiLongCurve(&diffcurList,&isSplitedCurve);
	//In order to observe,make ternimal point into little circle
	//makeTerpointToCircle();

	//delete[] isSplitedCurve;
	//isSplitedCurve = NULL;
	cout << "getCrossCurve done" << endl;
}

void MySvgFile::coarseJudgeByConvex(int c1, int c2, vector<DiffCurve> *diffcur)
{  
	DiffCurve curve1 = diffcur->at(c1);// [c1];
	DiffCurve curve2 = diffcur->at(c2);// [c2];
	bool cflag = false;
	vector<pair<mypoint2f, mypoint2f>> origiConvexSt = getConvexSegment(curve1);
	vector<pair<mypoint2f, mypoint2f>> compareConvexSt = getConvexSegment(curve2);
	for (unsigned int i = 0; i < origiConvexSt.size(); ++i)//先判断convex相交
	{
		for (unsigned int j = 0; j < compareConvexSt.size(); ++j)
		{
			if (isSegmentIntersected(origiConvexSt[i].first, origiConvexSt[i].second, compareConvexSt[j].first, compareConvexSt[j].second))//convex相交
			{
				cflag = true;
				break;
			}
		}
		if (cflag == true)
			break;
	}
	if (cflag == true)
	{
		diffcur->at(c1).itsct_cur.push_back(c2);
	}
	else//convex不相交,判断是不是在内部
	{
		vector<mypoint2f> origiConvexPt = getConvexPoints(curve1);
		vector<mypoint2f> compareConvexPt = getConvexPoints(curve2);
		cflag = isBoudEncircled(origiConvexPt, compareConvexPt);
		if (cflag == true)//convex相交,记录下两个curve
		{
			diffcur->at(c1).itsct_cur.push_back(c2);
		}
	}

	//vector<mypoint2f> origiConvexPt = getConvexPoints(curve1);
	//vector<mypoint2f> compareConvexPt = getConvexPoints(curve2);
	//cflag = isConvexIntersect(origiConvexPt, compareConvexPt);
	//if (cflag == true)//convex相交,记录下两个curve
	//{
	//	diffcur[c1].itsct_cur.push_back(c2);
	//}
	//else
	//{
	//	vector<mypoint2f> divide_p;
	//	vector<double> divide_t;
	//	if (isTerpointOnCurve(curve1, curve2, divide_t, divide_p))//T型相交.(水平垂直)
	//	{
	//		if (curve1.curveType != 'l')
	//			diffcur[c1].split_t.insert(diffcur[c1].split_t.end(), divide_t.begin(), divide_t.end());
	//		else
	//			diffcur[c1].split_point.insert(diffcur[c1].split_point.end(), divide_p.begin(), divide_p.end());
	//		marksplit[curve1.index] = true;
	//	}
	//	else if (isTerpointOnCurve(curve2, curve1, divide_t, divide_p))
	//	{
	//		if (curve2.curveType != 'l')
	//			diffcur[c2].split_t.insert(diffcur[c2].split_t.end(), divide_t.begin(), divide_t.end());
	//		else
	//			diffcur[c2].split_point.insert(diffcur[c2].split_point.end(), divide_p.begin(), divide_p.end());
	//		marksplit[curve2.index] = true;
	//	}
	//}
}
vector<pair<mypoint2f, mypoint2f>> MySvgFile::getConvexSegment(DiffCurve curve)
{
	vector<pair<mypoint2f, mypoint2f>> pver;
	if (curve.curveType == 'c')
	{
		cubicBezier cbcur = cubicBezier_vec_origin[curve.curveIndex];
		pver.push_back(make_pair(cbcur.p0, cbcur.p1));
		pver.push_back(make_pair(cbcur.p1, cbcur.p2));
		pver.push_back(make_pair(cbcur.p2, cbcur.p3));
		pver.push_back(make_pair(cbcur.p3, cbcur.p0));
	}
	else if (curve.curveType == 'q')
	{
		quaBezier qbcur = quadraticBezier_vec_origin[curve.curveIndex];
		pver.push_back(make_pair(qbcur.p0, qbcur.p1));
		pver.push_back(make_pair(qbcur.p1, qbcur.p2));
		pver.push_back(make_pair(qbcur.p2, qbcur.p0));
	}
	else if (curve.curveType == 'l')
	{
		vector<mypoint2f> pline = polyline_vec_origin[curve.curveIndex].pointlist;
		for (unsigned int i = 1; i < pline.size(); ++i)
		{
			pver.push_back(make_pair(pline[i - 1], pline[i]));
		}
	}
	return pver;
}
vector<mypoint2f> MySvgFile::getConvexPoints(DiffCurve curve)
{
	vector<mypoint2f> pver;
	if (curve.curveType == 'c')
	{
		cubicBezier cbcur = cubicBezier_vec_origin[curve.curveIndex];
		pver.push_back(cbcur.p0);
		pver.push_back(cbcur.p1);
		pver.push_back(cbcur.p2);
		pver.push_back(cbcur.p3);
	}
	else if (curve.curveType == 'q')
	{
		quaBezier qbcur = quadraticBezier_vec_origin[curve.curveIndex];
		pver.push_back(qbcur.p0);
		pver.push_back(qbcur.p1);
		pver.push_back(qbcur.p2);
	}
	else if (curve.curveType == 'l')
	{
		vector<mypoint2f> pline = polyline_vec_origin[curve.curveIndex].pointlist;
		for (unsigned int i = 0; i < pline.size(); ++i)
		{
			pver.push_back(pline[i]);
		}
	}
	return pver;
}
bool MySvgFile::isBoudEncircled(vector<mypoint2f> convex1, vector<mypoint2f> convex2)
{
	double max_x1 = -65535; double min_x1 = 65535; double max_y1 = -65535; double min_y1 = 65535;
	double max_x2 = -65535; double min_x2 = 65535; double max_y2 = -65535; double min_y2 = 65535;
	for (unsigned int i = 0; i < convex1.size(); ++i)
	{
		if (max_x1 < convex1[i].x)
			max_x1 = convex1[i].x;
		if (min_x1 > convex1[i].x)
			min_x1 = convex1[i].x;
		if (max_y1 < convex1[i].y)
			max_y1 = convex1[i].y;
		if (min_y1 > convex1[i].y)
			min_y1 = convex1[i].y;
	}
	for (unsigned int j = 0; j < convex2.size(); ++j)
	{
		if (max_x2 < convex2[j].x)
			max_x2 = convex2[j].x;
		if (min_x2 > convex2[j].x)
			min_x2 = convex2[j].x;
		if (max_y2 < convex2[j].y)
			max_y2 = convex2[j].y;
		if (min_y2 > convex2[j].y)
			min_y2 = convex2[j].y;
	}
	bool iflag = false;
	for (unsigned int i = 0; i < convex1.size(); ++i)
	{
		if ((convex1[i].x > min_x2&&convex1[i].x < max_x2) && (convex1[i].y>min_y2&&convex1[i].y<max_y2))
		{
			iflag = true;
			return iflag;
		}
	}
	for (unsigned int j = 0; j < convex2.size(); ++j)
	{
		if ((convex2[j].x > min_x1&&convex2[j].x<max_x1) && (convex2[j].y>min_y1&&convex2[j].y < max_y1))
		{
			iflag = true;
			return iflag;
		}
	}
	return iflag;
}
void MySvgFile::deleOrigiLongCurve(vector<DiffCurve> *diffcurlist, vector<bool> *marksplit)
{
	int cbsize = (int)cubicBezier_vec_origin.size();
	int qbsize = (int)quadraticBezier_vec_origin.size();
	int plsize = (int)polyline_vec_origin.size();
	int diffcursize = cbsize + qbsize + plsize;
	vector<cubicBezier> ncb_vec;
	vector<quaBezier> nqb_vec;
	vector<PolyLine> npl_vec;
	int i;
	for (i = 0; i < cbsize; ++i)
	{
		if (marksplit->at(i) == false)
		{
			ncb_vec.push_back(cubicBezier_vec_origin[i]);
		}
	}
	for (i = cbsize; i < (cbsize + qbsize); ++i)
	{
		if (marksplit->at(i) == false)
		{
			nqb_vec.push_back(quadraticBezier_vec_origin[i - cbsize]);
		}
	}
	for (i = (cbsize + qbsize); i < diffcursize; ++i)
	{
		if (marksplit->at(i) == false)
		{
			npl_vec.push_back(polyline_vec_origin[i - (cbsize + qbsize)]);
		}
	}
	//cubicBezier_vec.clear(); cubicBezier_vec.swap(vector<cubicBezier>());
	//quadraticBezier_vec.clear(); quadraticBezier_vec.swap(vector<quaBezier>());
	//polyline_vec.clear(); polyline_vec.swap(vector<PolyLine>());

	cubicBezier_vec.insert(cubicBezier_vec.end(), ncb_vec.begin(), ncb_vec.end());
	quadraticBezier_vec.insert(quadraticBezier_vec.end(), nqb_vec.begin(), nqb_vec.end());
	polyline_vec.insert(polyline_vec.end(), npl_vec.begin(), npl_vec.end());
	polyline_vec.insert(polyline_vec.end(),added_polyline.begin(),added_polyline.end());
}
void MySvgFile::makeTerpointToCircle()
{
	for (unsigned int i = 0; i < added_polyline.size(); ++i)
	{
		mycircle acir;
		acir.cx = added_polyline[i].pointlist.front().x;
		acir.cy = added_polyline[i].pointlist.front().y;
		acir.radius = 0.5;
		circle_vec.push_back(acir);
		acir.cx = added_polyline[i].pointlist.back().x;
		acir.cy = added_polyline[i].pointlist.back().y;
		acir.radius = 0.5;
		circle_vec.push_back(acir);
	}

	/*for (unsigned int i = 0; i < cubicBezier_vec.size(); ++i)
	{
		mycircle acir;
		acir.cx = cubicBezier_vec[i].p0.x;
		acir.cy = cubicBezier_vec[i].p0.y;
		acir.radius = 1;
		circle_vec.push_back(acir);
		acir.cx = cubicBezier_vec[i].p3.x;
		acir.cy = cubicBezier_vec[i].p3.y;
		acir.radius = 1;
		circle_vec.push_back(acir);
	}*/

	//for (unsigned int i = 0; i < quadraticBezier_vec.size(); ++i)
	//{
	//	mycircle acir;
	//	acir.cx = quadraticBezier_vec[i].p0.x;
	//	acir.cy = quadraticBezier_vec[i].p0.y;
	//	acir.radius = 2;
	//	circle_vec.push_back(acir);
	//	acir.cx = quadraticBezier_vec[i].p2.x;
	//	acir.cy = quadraticBezier_vec[i].p2.y;
	//	acir.radius = 2;
	//	circle_vec.push_back(acir);
	//}
	//for (unsigned int i = 0; i < polyline_vec.size(); ++i)
	//{
	//	mycircle acir;
	//	acir.cx = polyline_vec[i].front().x;
	//	acir.cy = polyline_vec[i].front().y;
	//	acir.radius = 2;
	//	circle_vec.push_back(acir);
	//	acir.cx = polyline_vec[i].back().x;
	//	acir.cy = polyline_vec[i].back().y;
	//	acir.radius = 2;
	//	circle_vec.push_back(acir);
	//}
}

void MySvgFile::getCrossCurvePosi(vector<DiffCurve> *diffcurlist, vector<bool> *isSplited)
{
	/*normal condition*/
	int dcsize = diffcurlist->size();
	for (int i = 0; i < dcsize; ++i)  //size可能在isTerpointOnCurve那个地方 动态增加，需要固定遍历范围
	{
		DiffCurve curve1 = diffcurlist->at(i);
		vector<int> cur_vec = diffcurlist->at(i).itsct_cur;
		for (unsigned int j = 0; j < cur_vec.size(); ++j)
		{
			if (cur_vec[j] > i)////!!
			{
				int curve_j = cur_vec[j];
				DiffCurve curve2 = diffcurlist->at(curve_j);
				vector<double> st1, st2; vector<mypoint2f> sp1, sp2;
				
				double mini_Terminate_dis = 0;
				bool itsflag = false;
				if (isTerpointOnCurve(curve1, curve2, st1, sp1,mini_Terminate_dis))//判断curve2的两端点之一是否在curve1上，st1, sp1记录的是curve1的split
				{
					itsflag = true;
					if (!st1.empty())
						isSplited->at(curve1.index) = true;
					if (!sp1.empty())
						isSplited->at(curve1.index) = true;
				}
				else if (isTerpointOnCurve(curve2, curve1, st2, sp2, mini_Terminate_dis))//curve1的两端点之一在curve2上...
				{
					itsflag = true;
					if (!st2.empty())
						isSplited->at(curve2.index) = true;
					if (!sp2.empty())
						isSplited->at(curve2.index) = true;
				}
				if (itsflag == true)
					addSplitPosi(diffcurlist, i, cur_vec[j], st1, st2, sp1, sp2);
				st1.clear(); st2.clear(); sp1.clear(); sp2.clear();
				itsflag = twoDiffCurveIntersct(curve1, curve2, st1, st2, sp1, sp2);
				if (itsflag == true)
				{
					if (curve1.curveType != 'l'&&curve2.curveType != 'l')
					{
						if (!st1.empty())
							isSplited->at(curve1.index) = true;
						if (!st2.empty())
							isSplited->at(curve2.index) = true;
					}
					else if (curve1.curveType == 'l'&&curve2.curveType == 'l')
					{
						if (!sp1.empty())
							isSplited->at(curve1.index) = true;
						if (!sp2.empty())
							isSplited->at(curve2.index) = true;
					}
					else
					{
						if (curve1.curveType == 'l')
						{
							if (!sp1.empty())
								isSplited->at(curve1.index) = true;
							if (!st1.empty())
								isSplited->at(curve2.index) = true;
						}
						else//curve2='l'
						{
							if (!sp1.empty())
								isSplited->at(curve2.index) = true;
							if (!st1.empty())
								isSplited->at(curve1.index) = true;
						}
					}
					addSplitPosi(diffcurlist, i, cur_vec[j], st1, st2, sp1, sp2);
				}

				////先判断 点是不是在线上
				//double tm_dis1 = 0, tm_dis2 = 0; //curve2两端点距离curve1 的距离
				//double mini_dis = 0;
				//bool onflag = false;
				//if (isTerpointOnCurve(curve1, curve2, st1, sp1, tm_dis1))//判断curve2的两端点之一是否在curve1上，st1, sp1记录的是curve1的split
				//{
				//	onflag = true;
				//	if (!st1.empty() || !sp1.empty())
				//		isSplited->at(curve1.index) = true;//curve1被分开，置为true
				//}
				//else if (isTerpointOnCurve(curve2, curve1, st2, sp2, tm_dis2))//curve1的两端点之一在curve2上...
				//{
				//	onflag = true;
				//	if (!st2.empty() || !sp2.empty())
				//		isSplited->at(curve2.index) = true;
				//}
				//if (onflag == false)//若端点都不在线上，则判断 距离是不是<1,
				//{
				//	DiffCurve extend_cur, cross_cur;
				//	mypoint2f ext_point(0, 0);
				//	if (tm_dis1 < tm_dis2)  //curve1被分割，求curve2的延长线
				//	{
				//		mini_dis = tm_dis1;
				//		extend_cur = curve2;
				//		cross_cur = curve1;
				//		ext_point = sp1.front();
				//	}
				//	else  //求curve1的延长线
				//	{
				//		mini_dis = tm_dis2;
				//		extend_cur = curve1;
				//		cross_cur = curve2;
				//		ext_point = sp2.front();
				//	}
				//	if (mini_dis < 2)//距离>1, 可以考虑延长curve一个单位长度，延长的过程中计算交叉点,添加新延长的short line
				//	{
				//		PolyLine ext_l;
				//		getExtendedLine(cross_cur, extend_cur, ext_point, ext_l,3);//计算延长线ext_l,从ext_point到向外延长的点			
				//		polyline_vec_origin.push_back(ext_l);
				//		DiffCurve dtcur;
				//		dtcur.curveIndex = ext_l.origi_index;
				//		dtcur.curveType = 'l';
				//		dtcur.index = (int)diffcurlist->size();
				//		diffcurlist->push_back(dtcur);
				//		//计算cross_curve 与延长线 的交点。若相交，记得设置：onflag=true!  isSplited[curve2.index] = true!!
				//		if (!st1.empty())
				//			st1.swap(vector<double>());
				//		if (!st2.empty())
				//			st2.swap(vector<double>());
				//		if (!sp1.empty())
				//			sp1.swap(vector<mypoint2f>());
				//		if (!sp2.empty())
				//			sp2.swap(vector<mypoint2f>());
				//		bool extflag = twoDiffCurveIntersct(cross_cur, dtcur, st1, st2, sp1, sp2);//cross_cur对应st1 sp1
				//		if (extflag == true)
				//		{
				//			//cross_cur 可能是bezier/line ， dtcur一定是 line
				//			if (cross_cur.curveType == 'l')
				//			{
				//				if (!sp1.empty())
				//					isSplited->at(cross_cur.index) = true;
				//				if (!sp2.empty())
				//				{
				//					isSplited->push_back(false);//isSplited[dtcur.index] = true;
				//					polyline_vec_origin.back().pointlist.pop_back(); 		//把延长线，截取有用的一段，剩下的不要,所以设置标志位=false
				//					polyline_vec_origin.back().pointlist.push_back(sp2.front());
				//				}
				//			}
				//			else
				//			{
				//				if (!st1.empty())
				//					isSplited->at(cross_cur.index) = true;
				//				if (!sp1.empty())
				//				{
				//					isSplited->push_back(false);//isSplited[dtcur.index] = true;
				//					polyline_vec_origin.back().pointlist.pop_back(); 		//把延长线，截取有用的一段，剩下的不要,所以设置标志位=false
				//					polyline_vec_origin.back().pointlist.push_back(sp1.front());
				//				}
				//			}
				//			addSplitPosi(diffcurlist, cross_cur.index, dtcur.index, st1, st2, sp1, sp2);// 延长相交时,添加到相应的vector中
				//			polyline_extend.push_back(polyline_vec_origin.back());//测试用
				//		}
				//		else
				//		{
				//			polyline_vec_origin.pop_back();
				//			diffcurlist->pop_back();
				//		}
				//	}
				//}
				//else//正常oncurve时，添加到相应的vector中
				//{
				//	if (curve1.curveType != 'l'&&curve2.curveType != 'l')
				//	{
				//		if (!st1.empty())
				//			isSplited->at(curve1.index) = true;
				//		if (!st2.empty())
				//			isSplited->at(curve2.index) = true;
				//	}
				//	else if (curve1.curveType == 'l'&&curve2.curveType == 'l')
				//	{
				//		if (!sp1.empty())
				//			isSplited->at(curve1.index) = true;
				//		if (!sp2.empty())
				//			isSplited->at(curve2.index) = true;
				//	}
				//	else
				//	{
				//		if (curve1.curveType == 'l')
				//		{
				//			if (!sp1.empty())
				//				isSplited->at(curve1.index) = true;
				//			if (!st1.empty())
				//				isSplited->at(curve2.index) = true;
				//		}
				//		else//curve2='l'
				//		{
				//			if (!sp1.empty())
				//				isSplited->at(curve2.index) = true;
				//			if (!st1.empty())
				//				isSplited->at(curve1.index) = true;
				//		}
				//	}
				//	addSplitPosi(diffcurlist, i, curve_j, st1, st2, sp1, sp2);
				//}
				////再判断是不是相交
				//st1.swap(vector<double>()); st2.swap(vector<double>()); sp1.swap(vector<mypoint2f>()); sp2.swap(vector<mypoint2f>());
				//bool itsflag = twoDiffCurveIntersct(curve1, curve2, st1, st2, sp1, sp2);
				//if (itsflag == true)
				//{
				//	if (curve1.curveType != 'l'&&curve2.curveType != 'l')
				//	{
				//		if (!st1.empty())
				//			isSplited->at(curve1.index) = true;
				//		if (!st2.empty())
				//			isSplited->at(curve2.index) = true;
				//	}
				//	else if (curve1.curveType == 'l'&&curve2.curveType == 'l')
				//	{
				//		if (!sp1.empty())
				//			isSplited->at(curve1.index) = true;
				//		if (!sp2.empty())
				//			isSplited->at(curve2.index) = true;
				//	}
				//	else
				//	{
				//		if (curve1.curveType == 'l')
				//		{
				//			if (!sp1.empty())
				//				isSplited->at(curve1.index) = true;
				//			if (!st1.empty())
				//				isSplited->at(curve2.index) = true;
				//		}
				//		else//curve2='l'
				//		{
				//			if (!sp1.empty())
				//				isSplited->at(curve2.index) = true;
				//			if (!st1.empty())
				//				isSplited->at(curve1.index) = true;
				//		}
				//	}
				//	addSplitPosi(diffcurlist, i, curve_j, st1, st2, sp1, sp2);
				//}
				////if (itsflag == true)
				////addSplitPosi(diffcurlist, i, curve_j, st1, st2, sp1, sp2);
			}
		}
	}
	cout << "add short extend line: " << diffcurlist->size() - dcsize << endl;
}
bool MySvgFile::isTerpointOnCurve(DiffCurve curve1, DiffCurve curve2, vector<double> &split_t, vector<mypoint2f> &split_p, double &mini_dis)
{
	//split_t,split_p在相交的情况下返回相交的点。在不相交的时候split_p返回curve2 curve2 curve2的其中一个端点！！！
	//check whether curve2's terminal points are on curve1
	bool flag = false;
	vector<mypoint2f> terpoint = getTerpoint(curve2);//判断curve2的两个端点是否在curve1上
	double diss1 = 6553500; //curve2的其中一个端点到curve1 的距离
	double diss2 = 6553500;//curve2的另一个端点到curve1 的距离
	if (curve1.curveType == 'l')//check whether curve2's terminate are on curve1. curve1 is a polyline
	{
		vector<mypoint2f> plline = polyline_vec_origin[curve1.curveIndex].pointlist;
		mypoint2f p1 = plline.front();
		mypoint2f p2 = plline.back();
		diss1 = pointToLineDistance(terpoint[0], p1, p2);
		diss2 = pointToLineDistance(terpoint[1], p1, p2);
		if (isOnSegment(terpoint[0], p1, p2))// && LineSegmentLength(terpoint[0], p1)>eps_isOnsegment && LineSegmentLength(terpoint[0], p2)>eps_isOnsegment)
		{
			split_p.push_back(terpoint[0]);
			flag = true;
		}
		else if (isOnSegment(terpoint[1], p1, p2))// && LineSegmentLength(terpoint[1], p1)>eps_isOnsegment && LineSegmentLength(terpoint[1], p2)>eps_isOnsegment)
		{
			split_p.push_back(terpoint[1]);
			flag = true;
		}
	}
	else if (curve1.curveType == 'q')
	{
		quaBezier qbcur = quadraticBezier_vec_origin[curve1.curveIndex];
		double a = qbcur.p0.x - 2 * qbcur.p1.x + qbcur.p2.x;
		double b = -2 * qbcur.p0.x + 2 * qbcur.p1.x;
		double c = qbcur.p0.x - terpoint[0].x;
		//double a = qbcur.p0.x - 2 * qbcur.p1.x + qbcur.p2.x + qbcur.p0.y - 2 * qbcur.p1.y + qbcur.p2.y;//double a = qbcur.p0.x - 2 * qbcur.p1.x + qbcur.p2.x;
		//double b = -2 * qbcur.p0.x + 2 * qbcur.p1.x -2 * qbcur.p0.y + 2 * qbcur.p1.y;//double b = -2 * qbcur.p0.x + 2 * qbcur.p1.x;
		//double c = qbcur.p0.x - terpoint[0].x + qbcur.p0.y - terpoint[0].y;//double c = qbcur.p0.x - terpoint[0].x;
		vector<double> t_root = quadricRealPolynomial(a, b, c);
		for (unsigned int i = 0; i < t_root.size(); ++i)
		{
			if (t_root[i] > 0 && t_root[i] < 1)//if (t_root[i] > eps && (1 - t_root[i])>eps)
			{
				double y1 = pow((1 - t_root[i]), 2)*qbcur.p0.y + 2 * t_root[i] * (1 - t_root[i])*qbcur.p1.y + pow(t_root[i], 2)*qbcur.p2.y;
				double x1 = pow((1 - t_root[i]), 2)*qbcur.p0.x + 2 * t_root[i] * (1 - t_root[i])*qbcur.p1.x + pow(t_root[i], 2)*qbcur.p2.x;//terpoint[0].x;
				mypoint2f interp(x1, y1);
				diss1 = LineSegmentLength(interp, terpoint[0]);
				if (diss1< eps_isOnsegment)//&& LineSegmentLength(interp, qbcur.p0)>eps_isOnsegment && LineSegmentLength(interp, qbcur.p2)>eps_isOnsegment)//在线段上&&不在曲线两个端点上
				{
					split_t.push_back(t_root[i]);
					flag = true;
				}
			}
		}
		//if (flag == false)
		//{
		t_root.clear(); t_root.swap(vector<double>());
		c = qbcur.p0.x - terpoint[1].x;//another terminate point //c = qbcur.p0.x - terpoint[1].x + qbcur.p0.y - terpoint[1].y;//
		t_root = quadricRealPolynomial(a, b, c);
		for (unsigned int i = 0; i < t_root.size(); ++i)
		{
			if (t_root[i] > 0 && t_root[i] < 1)//if (t_root[i] > eps && (1 - t_root[i])>eps)
			{
				double y1 = pow(1 - t_root[i], 2)*qbcur.p0.y + 2 * t_root[i] * (1 - t_root[i])*qbcur.p1.y + pow(t_root[i], 2)*qbcur.p2.y;
				double x1 = pow(1 - t_root[i], 2)*qbcur.p0.x + 2 * t_root[i] * (1 - t_root[i])*qbcur.p1.x + pow(t_root[i], 2)*qbcur.p2.x;//terpoint[1].x;
				mypoint2f interp(x1, y1);
				diss2 = LineSegmentLength(interp, terpoint[1]);
				if (diss2< eps_isOnsegment)// && LineSegmentLength(interp, qbcur.p0)>eps_isOnsegment && LineSegmentLength(interp, qbcur.p2)>eps_isOnsegment)///
				{
					split_t.push_back(t_root[i]);
					flag = true;
				}
			}
		}
		//}
		//diss1 < diss2 ? mini_dis = diss1 : mini_dis = diss2;
	}
	else if (curve1.curveType == 'c')
	{
		cubicBezier cbcurve = cubicBezier_vec_origin[curve1.curveIndex];
		//double A = -cbcur.p0.x + 3 * cbcur.p1.x - 3 * cbcur.p2.x + cbcur.p3.x -cbcur.p0.y + 3 * cbcur.p1.y - 3 * cbcur.p2.y + cbcur.p3.y;
		//double B = 3 * cbcur.p0.x - 6 * cbcur.p1.x + 3 * cbcur.p2.x + 3 * cbcur.p0.y - 6 * cbcur.p1.y + 3 * cbcur.p2.y;
		//double C = -3 * cbcur.p0.x + 3 * cbcur.p1.x -3 * cbcur.p0.y + 3 * cbcur.p1.y;
		//double D = cbcur.p0.x - terpoint[0].x + cbcur.p0.y - terpoint[0].y;
		double A = -cbcurve.p0.x + 3 * cbcurve.p1.x - 3 * cbcurve.p2.x + cbcurve.p3.x;
		double B = 3 * cbcurve.p0.x - 6 * cbcurve.p1.x + 3 * cbcurve.p2.x;
		double C = -3 * cbcurve.p0.x + 3 * cbcurve.p1.x;
		double D = cbcurve.p0.x - terpoint[0].x;
		vector<double> t_root = cubicRealPolynomial(A, B, C, D);
		for (unsigned int i = 0; i < t_root.size(); ++i)
		{
			if (t_root[i] > 0 && t_root[i] < 1)//if (t_root[i] > eps && (1 - t_root[i])>eps)
			{
				double y1 = pow(1 - t_root[i], 3)*cbcurve.p0.y + 3 * t_root[i] * pow((1 - t_root[i]), 2)*cbcurve.p1.y + 3 * pow(t_root[i], 2)*(1 - t_root[i])*cbcurve.p2.y + pow(t_root[i], 3)*cbcurve.p3.y;
				double x1 = pow(1 - t_root[i], 3)*cbcurve.p0.x + 3 * t_root[i] * pow((1 - t_root[i]), 2)*cbcurve.p1.x + 3 * pow(t_root[i], 2)*(1 - t_root[i])*cbcurve.p2.x + pow(t_root[i], 3)*cbcurve.p3.x;
				mypoint2f interp(x1, y1);
				diss1 = LineSegmentLength(interp, terpoint[0]);
				if (diss1 < eps_isOnsegment)// && LineSegmentLength(interp, cbcur.p0)>eps_isOnsegment && LineSegmentLength(interp, cbcur.p3)>eps_isOnsegment)///
				{
					split_t.push_back(t_root[i]);
					flag = true;
				}
			}
		}
		//if (flag == false)
		//{
		D = cbcurve.p0.x - terpoint[1].x;//D = cbcur.p0.x - terpoint[1].x + cbcur.p0.y - terpoint[1].y;//
		t_root.clear(); t_root.swap(vector<double>());
		t_root = cubicRealPolynomial(A, B, C, D);
		for (unsigned int i = 0; i < t_root.size(); ++i)
		{
			if (t_root[i] > 0 && t_root[i] < 1)//if (t_root[i] > eps && (1 - t_root[i])>eps)
			{
				double y1 = pow(1 - t_root[i], 3)*cbcurve.p0.y + 3 * t_root[i] * pow((1 - t_root[i]), 2)*cbcurve.p1.y + 3 * pow(t_root[i], 2)*(1 - t_root[i])*cbcurve.p2.y + pow(t_root[i], 3)*cbcurve.p3.y;
				double x1 = pow(1 - t_root[i], 3)*cbcurve.p0.x + 3 * t_root[i] * pow((1 - t_root[i]), 2)*cbcurve.p1.x + 3 * pow(t_root[i], 2)*(1 - t_root[i])*cbcurve.p2.x + pow(t_root[i], 3)*cbcurve.p3.x;
				mypoint2f interp(x1, y1);
				diss2 = LineSegmentLength(interp, terpoint[1]);
				if (diss2 < eps_isOnsegment)// && LineSegmentLength(interp, cbcur.p0)>eps_isOnsegment && LineSegmentLength(interp, cbcur.p3)>eps_isOnsegment)///
				{
					split_t.push_back(t_root[i]);
					flag = true;
				}
			}
		}
		//}
		//diss1 < diss2 ? mini_dis = diss1 : mini_dis = diss2;
	}
	if (diss1 < diss2)
	{
		mini_dis = diss1;
		if (flag == false)//在不相交的时候split_p返回curve2 curve2 curve2的其中一个端点！！！
		{
			split_p.push_back(terpoint[0]);
		}
	}
	else
	{
		mini_dis = diss2;
		if (flag == false)
		{
			split_p.push_back(terpoint[1]);
		}
	}
	return flag;
}
vector<mypoint2f> MySvgFile::getTerpoint(DiffCurve curve)
{
	vector<mypoint2f> point_vec;
	if (curve.curveType == 'c')
	{
		cubicBezier cbcurve = cubicBezier_vec_origin[curve.curveIndex];
		point_vec.push_back(cbcurve.p0);
		point_vec.push_back(cbcurve.p3);
	}
	else if (curve.curveType == 'q')
	{
		quaBezier qbcurve = quadraticBezier_vec_origin[curve.curveIndex];
		point_vec.push_back(qbcurve.p0);
		point_vec.push_back(qbcurve.p2);
	}
	else if (curve.curveType == 'l')
	{
		point_vec.push_back(polyline_vec_origin[curve.curveIndex].pointlist.front());
		point_vec.push_back(polyline_vec_origin[curve.curveIndex].pointlist.back());
	}
	return point_vec;
}
void MySvgFile::getExtendedLine(DiffCurve cross_curve, DiffCurve ext_curve, mypoint2f ext_pt, PolyLine &ext_line, int scalor)
{
	//先计算ext_curve切线方向的单位向量 （作延长线）
	mypoint2f tan_direct(0, 0);
	mypoint2f ext_terminate = ext_pt;
	vector<mypoint2f> ter_pt = getTerpoint(ext_curve);//获得p0,p3 | p0,p2 | front back
	if (ext_curve.curveType == 'c')
	{
		cubicBezier acub = cubicBezier_vec_origin[ext_curve.curveIndex];
		if (ter_pt[0] == ext_pt)  //可以用==判断
		{
			tan_direct = normalization(acub.p1, acub.p0);
		}
		else
		{
			tan_direct = normalization(acub.p2, acub.p3);
		}
	}
	else if (ext_curve.curveType == 'q')
	{
		quaBezier qbcurve = quadraticBezier_vec_origin[ext_curve.curveIndex];
		if (ter_pt[0] == ext_pt)
		{
			tan_direct = normalization(qbcurve.p1, qbcurve.p0);
		}
		else
		{
			tan_direct = normalization(qbcurve.p1, qbcurve.p2);
		}
	}
	else if (ext_curve.curveType == 'l')
	{
		PolyLine apl = polyline_vec_origin[ext_curve.curveIndex];
		if (ter_pt[0] == ext_pt)  //可以用==判断
		{
			tan_direct = normalization(apl.pointlist.back(), apl.pointlist.front());
		}
		else
		{
			tan_direct = normalization(apl.pointlist.front(), apl.pointlist.back());
		}
	}
	//线
	tan_direct.x = tan_direct.x*scalor;
	tan_direct.y = tan_direct.y*scalor;
	ext_line.origi_index = (int)polyline_vec_origin.size();//不要push进去，再该函数之外push
	ext_line.origi_type = 'l';
	ext_line.pointlist.push_back(ext_terminate);  //有顺序的，要先push起始点
	ext_line.pointlist.push_back(tan_direct + ext_terminate);
}

void MySvgFile::addSplitPosi(vector<DiffCurve> *diffcurlist, int cur1, int cur2, vector<double> split_t1, vector<double> split_t2, vector<mypoint2f> split_p1, vector<mypoint2f> split_p2)
{
	if (diffcurlist->at(cur1).curveType != 'l'&&diffcurlist->at(cur2).curveType != 'l')
	{
		diffcurlist->at(cur1).split_t.insert(diffcurlist->at(cur1).split_t.end(), split_t1.begin(), split_t1.end()); //push_back(split_t1);
		diffcurlist->at(cur2).split_t.insert(diffcurlist->at(cur2).split_t.end(), split_t2.begin(), split_t2.end());
	}
	else if (diffcurlist->at(cur1).curveType != 'l'&&diffcurlist->at(cur2).curveType == 'l')
	{
		diffcurlist->at(cur1).split_t.insert(diffcurlist->at(cur1).split_t.end(), split_t1.begin(), split_t1.end());
		diffcurlist->at(cur2).split_point.insert(diffcurlist->at(cur2).split_point.end(), split_p1.begin(), split_p1.end());
	}
	else if (diffcurlist->at(cur1).curveType == 'l'&&diffcurlist->at(cur2).curveType != 'l')
	{
		diffcurlist->at(cur1).split_point.insert(diffcurlist->at(cur1).split_point.end(), split_p1.begin(), split_p1.end());
		diffcurlist->at(cur2).split_t.insert(diffcurlist->at(cur2).split_t.end(), split_t1.begin(), split_t1.end());
	}
	else if (diffcurlist->at(cur1).curveType == 'l'&&diffcurlist->at(cur2).curveType == 'l')
	{
		diffcurlist->at(cur1).split_point.insert(diffcurlist->at(cur1).split_point.end(), split_p1.begin(), split_p1.end());
		diffcurlist->at(cur2).split_point.insert(diffcurlist->at(cur2).split_point.end(), split_p2.begin(), split_p2.end());
	}
}
bool MySvgFile::twoDiffCurveIntersct(DiffCurve curve1, DiffCurve curve2, vector<double> &split_t1, vector<double> &split_t2, vector<mypoint2f> &split_p1, vector<mypoint2f> &split_p2)
{
	int ret_i = 0;
	bool iflag = false;
	/*if (curve1.curveType != 'l'&&curve2.curveType != 'l')
	{
		if (curve1.curveType == 'c'&&curve2.curveType == 'c')
		{
			cubicBezier cubBcur1 = cubicBezier_vec[curve1.curveIndex];
			cubicBezier cubBcur2 = cubicBezier_vec[curve2.curveIndex];
			iflag = twoCubicBzrIntersect(cubBcur1, cubBcur2, split_t1, split_t2);
		}
		else if (curve1.curveType == 'q'&&curve2.curveType == 'q')
		{
			quaBezier QBcur1 = quadraticBezier_vec[curve1.curveIndex];
			quaBezier QBcur2 = quadraticBezier_vec[curve2.curveIndex];
			iflag = twoQuaBzrIntersect(QBcur1, QBcur2, split_t1, split_t2);
		}
		else if (curve1.curveType == 'c'&&curve2.curveType == 'q')
		{
			cubicBezier CBcur = cubicBezier_vec[curve1.curveIndex];
			quaBezier QBcur = quadraticBezier_vec[curve2.curveIndex];
			iflag = CBandQBIntersect(CBcur, QBcur, split_t1, split_t2);
		}
		else if (curve1.curveType == 'q'&&curve2.curveType == 'c')
		{
			cubicBezier CBcur = cubicBezier_vec[curve2.curveIndex];
			quaBezier QBcur = quadraticBezier_vec[curve1.curveIndex];
			iflag = CBandQBIntersect(CBcur, QBcur, split_t1, split_t2);
		}
		if (iflag == true)
		{
			if (checkTValidity(split_t1))
			{
				if (checkTValidity(split_t2))
					ret_i = 3;
				else
					ret_i = 1;
			}
			else
			{
				if (checkTValidity(split_t2))
					ret_i = 2;
			}			
		}
	}
	else if (curve1.curveType == 'l'&&curve2.curveType == 'l')
	{
		vector<mypoint2f> plLine1 = polyline_vec[curve1.curveIndex].pointlist;
		vector<mypoint2f> plLine2 = polyline_vec[curve2.curveIndex].pointlist;
		iflag = twoPolyIntersect(plLine1, plLine2, split_p);
	}
	else 
	{
		mypoint2f startp, endp;
		if (curve1.curveType == 'c'&&curve2.curveType == 'l')
		{
			cubicBezier CBcur = cubicBezier_vec[curve1.curveIndex];
			vector<mypoint2f> plLine = polyline_vec[curve2.curveIndex].pointlist;
			startp = plLine.front(); endp = plLine.back();
			iflag = CBandPLIntersect(CBcur, plLine, split_t1, split_p);
		}
		else if (curve1.curveType == 'l'&&curve2.curveType == 'c')
		{
			cubicBezier CBcur = cubicBezier_vec[curve2.curveIndex];
			vector<mypoint2f> plLine = polyline_vec[curve1.curveIndex].pointlist;
			startp = plLine.front(); endp = plLine.back();
			iflag = CBandPLIntersect(CBcur, plLine, split_t1, split_p);
		}
		else if (curve1.curveType == 'q'&&curve2.curveType == 'l')
		{
			quaBezier QBcur = quadraticBezier_vec[curve1.curveIndex];
			vector<mypoint2f> plLine = polyline_vec[curve2.curveIndex].pointlist;
			startp = plLine.front(); endp = plLine.back();
			iflag = QBandPLIntersect(QBcur, plLine, split_t1, split_p);
		}
		else if (curve1.curveType == 'l'&&curve2.curveType == 'q')
		{
			quaBezier QBcur = quadraticBezier_vec[curve2.curveIndex];
			vector<mypoint2f> plLine = polyline_vec[curve1.curveIndex].pointlist;
			startp = plLine.front(); endp = plLine.back();
			iflag = QBandPLIntersect(QBcur, plLine, split_t1, split_p);
		}
	}*/		

	mypoint2f startp, endp;
	if (curve1.curveType == 'l'&&curve2.curveType == 'l')
	{
		vector<mypoint2f> plLine1 = polyline_vec_origin[curve1.curveIndex].pointlist;
		vector<mypoint2f> plLine2 = polyline_vec_origin[curve2.curveIndex].pointlist;
		iflag = twoPolyIntersect(plLine1, plLine2, split_p1, split_p2);
		if (iflag == true)
		{
			startp = plLine1.front(); endp = plLine1.back();
			checkPValidity(split_p1, startp, endp);
			startp = plLine2.front(); endp = plLine2.back();
			checkPValidity(split_p2, startp, endp);
		}		
	}
	else
	{
		if (curve1.curveType == 'c'&&curve2.curveType == 'c')
		{
			cubicBezier cubBcur1 = cubicBezier_vec_origin[curve1.curveIndex];
			cubicBezier cubBcur2 = cubicBezier_vec_origin[curve2.curveIndex];
			iflag = twoCubicBzrIntersect(cubBcur1, cubBcur2, split_t1, split_t2);
		}
		else if (curve1.curveType == 'q'&&curve2.curveType == 'q')
		{
			quaBezier QBcur1 = quadraticBezier_vec_origin[curve1.curveIndex];
			quaBezier QBcur2 = quadraticBezier_vec_origin[curve2.curveIndex];
			iflag = twoQuaBzrIntersect(QBcur1, QBcur2, split_t1, split_t2);
		}
		else if (curve1.curveType == 'c'&&curve2.curveType == 'q')
		{
			cubicBezier CBcur = cubicBezier_vec_origin[curve1.curveIndex];
			quaBezier QBcur = quadraticBezier_vec_origin[curve2.curveIndex];
			iflag = CBandQBIntersect(CBcur, QBcur, split_t1, split_t2);
		}
		else if (curve1.curveType == 'q'&&curve2.curveType == 'c')
		{
			cubicBezier CBcur = cubicBezier_vec_origin[curve2.curveIndex];
			quaBezier QBcur = quadraticBezier_vec_origin[curve1.curveIndex];
			iflag = CBandQBIntersect(CBcur, QBcur, split_t1, split_t2);
		}
		else if (curve1.curveType == 'c'&&curve2.curveType == 'l')
		{
			cubicBezier CBcur = cubicBezier_vec_origin[curve1.curveIndex];
			vector<mypoint2f> plLine = polyline_vec_origin[curve2.curveIndex].pointlist;
			startp = plLine.front(); endp = plLine.back();
			iflag = CBandPLIntersect(CBcur, plLine, split_t1, split_p1);
		}
		else if (curve1.curveType == 'l'&&curve2.curveType == 'c')
		{
			cubicBezier CBcur = cubicBezier_vec_origin[curve2.curveIndex];
			vector<mypoint2f> plLine = polyline_vec_origin[curve1.curveIndex].pointlist;
			startp = plLine.front(); endp = plLine.back();
			iflag = CBandPLIntersect(CBcur, plLine, split_t1, split_p1);
		}
		else if (curve1.curveType == 'q'&&curve2.curveType == 'l')
		{
			quaBezier QBcur = quadraticBezier_vec_origin[curve1.curveIndex];
			vector<mypoint2f> plLine = polyline_vec_origin[curve2.curveIndex].pointlist;
			startp = plLine.front(); endp = plLine.back();
			iflag = QBandPLIntersect(QBcur, plLine, split_t1, split_p1);
		}
		else if (curve1.curveType == 'l'&&curve2.curveType == 'q')
		{
			quaBezier QBcur = quadraticBezier_vec_origin[curve2.curveIndex];
			vector<mypoint2f> plLine = polyline_vec_origin[curve1.curveIndex].pointlist;
			startp = plLine.front(); endp = plLine.back();
			iflag = QBandPLIntersect(QBcur, plLine, split_t1, split_p1);
		}
		if (iflag == true)
		{
			checkTValidity(split_t1);
			checkTValidity(split_t2);
			checkPValidity(split_p1, startp, endp);
		}		
	}	
	return iflag;
}
bool MySvgFile::checkTValidity(vector<double> &split_t)
{
	bool vali_flag = false;
	//if (split_t.size() == 1)
	//{
	//	if ((split_t[0] - 0 > eps) && (1-split_t[0]> eps))
	//		vali_flag = true;
	//}
	//else if(split_t.size()>1)
	//{
	//	sort(split_t.begin(), split_t.end());
	//	if (split_t.front() - 0 < eps)
	//		split_t.erase(split_t.begin());
	//	if (1 - split_t.back() < eps)
	//		split_t.erase(split_t.end()-1);
	//	if (!split_t.empty())
	//	{
	//		vali_flag = true;
	//		if (split_t.size() > 1)
	//		{
	//			vector<double>::iterator iter = split_t.begin();
	//			while (iter != split_t.end())
	//			{
	//				vector<double>::iterator iter_inn = iter + 1;
	//				while (iter_inn != split_t.end())
	//				{
	//					if ((*iter_inn - *iter)<eps)//if (*iter == *iter_inn),(约等于)
	//						iter_inn = split_t.erase(iter_inn);
	//					else
	//						iter_inn++;
	//				}
	//				++iter;
	//			}
	//		}						
	//	}
	//}
	if (!split_t.empty())
	{
		vector<double>::iterator iter = split_t.begin();
		while (iter != split_t.end())
		{
			if ((*iter - 0 > eps) && (1 - *iter> eps))
			{
				vali_flag = true;
				++iter;
			}
			else
			{
				if (split_t.size()>1)//保留一个极小值
					iter = split_t.erase(iter);
				else
					++iter;
			}
			//else
				//iter = split_t.erase(iter);
		}
	}	
	return vali_flag;
}
bool MySvgFile::checkPValidity(vector<mypoint2f> &split_p, mypoint2f sp, mypoint2f ep)
{
	bool vali_flag = false;
	if (!split_p.empty())
	{
		vector<mypoint2f>::iterator iter = split_p.begin();
		while (iter != split_p.end())
		{
			if (LineSegmentLength(*iter, sp) > eps && LineSegmentLength(*iter, ep) > eps)
			{
				vali_flag = true;
				++iter;
			}
			else
			{
				if (split_p.size()>1)//保留一个极小值
					iter = split_p.erase(iter);
				else
					++iter;
				//iter = split_p.erase(iter);
			}
		}
	}	
	return vali_flag;
}

//bool MySvgFile::twoCubicBzrIntersect(cubicBezier cubBcur1, cubicBezier cubBcur2, vector<double> &split_t1, vector<double> &split_t2)
//{
//	bool isIntersect = false;
//	double start_t = 0;
//	while (start_t <= 1)
//	{
//		double t = start_t;
//		double t_prior = 0;
//		double step = 0.02;
//		cubicBezier tmp_cur1 = cubBcur1;
//		cubicBezier tmp_cur2 = cubBcur2;
//		double deviation = 1000;
//		int iteration = 0;
//		double store_t; bool findflag = false; bool turnflag = false;
//		while (t <= 1 && t >= 0 && deviation > eps && iteration < 30)
//		{
//			double a1_cur1, a2_cur1, a3_cur1, a4_cur1;
//			double b1_cur1, b2_cur1, b3_cur1, b4_cur1;
//			CubicBezrCoefficient(t, a1_cur1, a2_cur1, a3_cur1, a4_cur1);
//			CubicBezrDerivative(t, b1_cur1, b2_cur1, b3_cur1, b4_cur1);
//
//			double x_cur1 = tmp_cur1.p0.x*a1_cur1 + tmp_cur1.p1.x*a2_cur1 + tmp_cur1.p2.x*a3_cur1 + tmp_cur1.p3.x*a4_cur1;
//			double y_cur1 = tmp_cur1.p0.y*a1_cur1 + tmp_cur1.p1.y*a2_cur1 + tmp_cur1.p2.y*a3_cur1 + tmp_cur1.p3.y*a4_cur1;
//			double dx = tmp_cur1.p0.x*b1_cur1 + tmp_cur1.p1.x*b2_cur1 + tmp_cur1.p2.x*b3_cur1 + tmp_cur1.p3.x*b4_cur1;
//			double dy = tmp_cur1.p0.y*b1_cur1 + tmp_cur1.p1.y*b2_cur1 + tmp_cur1.p2.y*b3_cur1 + tmp_cur1.p3.y*b4_cur1;
//			double a = dy / dx;
//			double b = -1;
//			double c = y_cur1 - a*x_cur1;
//
//			double A = a*(-tmp_cur2.p0.x + 3 * tmp_cur2.p1.x - 3 * tmp_cur2.p2.x + tmp_cur2.p3.x) +
//				b*(-tmp_cur2.p0.y + 3 * tmp_cur2.p1.y - 3 * tmp_cur2.p2.y + tmp_cur2.p3.y);
//			double B = a*(3 * tmp_cur2.p0.x - 6 * tmp_cur2.p1.x + 3 * tmp_cur2.p2.x) +
//				b*(3 * tmp_cur2.p0.y - 6 * tmp_cur2.p1.y + 3 * tmp_cur2.p2.y);
//			double C = a*(-3 * tmp_cur2.p0.x + 3 * tmp_cur2.p1.x) +
//				b*(-3 * tmp_cur2.p0.y + 3 * tmp_cur2.p1.y);
//			double D = a*tmp_cur2.p0.x + b*tmp_cur2.p0.y + c;
//			vector<double> t_root = cubicRealPolynomial(A, B, C, D);
//
//				vector<double>::iterator iter = t_root.begin();
//				while (iter != t_root.end())
//				{
//					if (*iter > 0 && *iter < 1)
//						iter++;
//					else
//						iter = t_root.erase(iter);
//				}
//				if (!t_root.empty())
//				{
//					bool firstflag = false;
//					for (unsigned int i = 0; i < t_root.size(); ++i)
//					{
//						if (t_root[i]>0 && t_root[i]<1)//if (t_root[i]>eps && (1 - t_root[i])>eps)//
//						{
//							double a1_cur2, a2_cur2, a3_cur2, a4_cur2;
//							CubicBezrCoefficient(t_root[i], a1_cur2, a2_cur2, a3_cur2, a4_cur2);
//							double x_cur2 = tmp_cur2.p0.x*a1_cur2 + tmp_cur2.p1.x*a2_cur2 + tmp_cur2.p2.x*a3_cur2 + tmp_cur2.p3.x*a4_cur2;
//							double y_cur2 = tmp_cur2.p0.y*a1_cur2 + tmp_cur2.p1.y*a2_cur2 + tmp_cur2.p2.y*a3_cur2 + tmp_cur2.p3.y*a4_cur2;
//							double  tmpd = sqrt(pow(x_cur1 - x_cur2, 2) + pow(y_cur1 - y_cur2, 2));
//							if (t_root.size() > 1)
//							{
//								if (firstflag == false)
//								{
//									firstflag = true;
//									deviation = tmpd;
//									t_prior = t;
//									t = t_root[i];
//								}
//								else
//								{
//									if (tmpd < deviation)//considering there are two more intersecting points,and select the nearest one
//									{
//										deviation = tmpd;
//										t_prior = t;
//										t = t_root[i];
//									}
//								}
//							}
//							else
//							{
//								deviation = tmpd;
//								t_prior = t;
//								t = t_root[i];
//							}
//						}
//					}
//					cubicBezier tmpcb = tmp_cur1;
//					tmp_cur1 = tmp_cur2;
//					tmp_cur2 = tmpcb;
//					findflag = true;////////////////			
//				}
//				else
//				{
//					t = t + step;
//				}
//
//			iteration = iteration + 1;
//			if (t > 1 || t < 0)//step=-0.02或0.02  只翻转一次/////////////////////
//			{
//				if (findflag == false && turnflag == false)
//				{
//					step = step*(-1);
//					t = store_t;
//					turnflag = true;
//				}
//			}/////////////////////
//		}
//
//		if (turnflag == true && split_t1.empty())/////////////////////
//			break;/////////////////////
//		if (deviation <= eps)
//		{
//			if (tmp_cur1.p0 == cubBcur1.p0 && tmp_cur1.p1 == cubBcur1.p1)//判断tmp_cur1是不是cubBcur1
//			{
//				split_t1.push_back(t); //split_t1 = t;//t-> tmp_cur1;
//				split_t2.push_back(t_prior);//split_t2 = t_prior;//t_prior-> tmp_cur2;
//			}
//			else
//			{
//				split_t1.push_back(t_prior);//split_t1 = t_prior;
//				split_t2.push_back(t);//split_t2 = t;
//			}
//			isIntersect = true;
//		}
//		start_t = start_t + 0.2;
//		if (start_t > 0.5&&split_t1.empty())
//			break;
//	}
//	return isIntersect;
//}

bool MySvgFile::twoCubicBzrIntersect(cubicBezier cubBcur1, cubicBezier cubBcur2, vector<double> &split_t1, vector<double> &split_t2)
{
	//double cflag;
	//tbb::parallel_for(tbb::blocked_range<size_t>(0, 10), [&](tbb::blocked_range<size_t>& r){for (size_t ii = r.begin(); ii != r.end(); ++ii){twoCubicBzrIntersect_para(ii, cubBcur1, cubBcur2, split_t1, split_t2); } });
	//if (split_t1.empty() && split_t2.empty())
	//	cflag = false;
	//else
	//	cflag = true;
	//return cflag;

	bool isIntersect = false;
	double start_t = 0;
	while (start_t <= 1)
	{
		bool first_t = false;////////////////////////////////////
		double t = start_t;
		double t_prior = 0;
		double step = 0.02;
		cubicBezier tmp_cur1 = cubBcur1;
		cubicBezier tmp_cur2 = cubBcur2;
		double deviation = 1000;
		int iteration = 0;
		while (t <= 1 && t >=0 && deviation > eps && iteration < 30)
		{
			double a1_cur1, a2_cur1, a3_cur1, a4_cur1;
			double b1_cur1, b2_cur1, b3_cur1, b4_cur1;
			CubicBezrCoefficient(t, a1_cur1, a2_cur1, a3_cur1, a4_cur1);
			CubicBezrDerivative(t, b1_cur1, b2_cur1, b3_cur1, b4_cur1);

			double x_cur1 = tmp_cur1.p0.x*a1_cur1 + tmp_cur1.p1.x*a2_cur1 + tmp_cur1.p2.x*a3_cur1 + tmp_cur1.p3.x*a4_cur1;
			double y_cur1 = tmp_cur1.p0.y*a1_cur1 + tmp_cur1.p1.y*a2_cur1 + tmp_cur1.p2.y*a3_cur1 + tmp_cur1.p3.y*a4_cur1;
			double dx = tmp_cur1.p0.x*b1_cur1 + tmp_cur1.p1.x*b2_cur1 + tmp_cur1.p2.x*b3_cur1 + tmp_cur1.p3.x*b4_cur1;
			double dy = tmp_cur1.p0.y*b1_cur1 + tmp_cur1.p1.y*b2_cur1 + tmp_cur1.p2.y*b3_cur1 + tmp_cur1.p3.y*b4_cur1;
			double a = dy / dx;
			double b = -1;
			double c = y_cur1 - a*x_cur1;

			double A = a*(-tmp_cur2.p0.x + 3 * tmp_cur2.p1.x - 3 * tmp_cur2.p2.x + tmp_cur2.p3.x) +
				b*(-tmp_cur2.p0.y + 3 * tmp_cur2.p1.y - 3 * tmp_cur2.p2.y + tmp_cur2.p3.y);
			double B = a*(3 * tmp_cur2.p0.x - 6 * tmp_cur2.p1.x + 3 * tmp_cur2.p2.x) +
				b*(3 * tmp_cur2.p0.y - 6 * tmp_cur2.p1.y + 3 * tmp_cur2.p2.y);
			double C = a*(-3 * tmp_cur2.p0.x + 3 * tmp_cur2.p1.x) +
				b*(-3 * tmp_cur2.p0.y + 3 * tmp_cur2.p1.y);
			double D = a*tmp_cur2.p0.x + b*tmp_cur2.p0.y + c;
			vector<double> t_root = cubicRealPolynomial(A, B, C, D);
			if (!t_root.empty())//t_root求出来的根是tmp_cur2上的t值
			{
				int suit_t = 0;
				for (unsigned int i = 0; i < t_root.size(); ++i)
				{
					if (t_root[i]>0 && t_root[i]<1)//if (t_root[i]>eps && (1 - t_root[i])>eps)//
						suit_t = suit_t + 1;
				}
				bool firstflag = false;
				double tmp_prior_t = 0;
				for (unsigned int i = 0; i < t_root.size(); ++i)
				{
					if (t_root[i]>0 && t_root[i]<1)//if (t_root[i]>eps && (1 - t_root[i])>eps)//
					{
						if (first_t == false && t_root[i]>0.95)
						{
							step = -0.02;
							first_t = true;
						}
						double a1_cur2, a2_cur2, a3_cur2, a4_cur2;
						CubicBezrCoefficient(t_root[i], a1_cur2, a2_cur2, a3_cur2, a4_cur2);
						double x_cur2 = tmp_cur2.p0.x*a1_cur2 + tmp_cur2.p1.x*a2_cur2 + tmp_cur2.p2.x*a3_cur2 + tmp_cur2.p3.x*a4_cur2;
						double y_cur2 = tmp_cur2.p0.y*a1_cur2 + tmp_cur2.p1.y*a2_cur2 + tmp_cur2.p2.y*a3_cur2 + tmp_cur2.p3.y*a4_cur2;
						double  tmpd = sqrt(pow(x_cur1 - x_cur2, 2) + pow(y_cur1 - y_cur2, 2));
						if (suit_t > 1)
						{
							if (firstflag == false)
							{
								firstflag = true;
								deviation = tmpd;
								tmp_prior_t = t;
								t_prior = t;
								t = t_root[i];
							}
							else
							{
								if (tmpd < deviation)//considering there are two more intersecting points,and select the nearest one
								{
									deviation = tmpd;
									t_prior = tmp_prior_t;
									//t_prior = t;
									t = t_root[i];
								}
							}
						}
						else
						{
							deviation = tmpd;
							t_prior = t;
							t = t_root[i];
						}
					}
				}
				if (suit_t > 0)
				{
					cubicBezier tmpcb = tmp_cur1;
					tmp_cur1 = tmp_cur2;
					tmp_cur2 = tmpcb;
				}
				else
				{
					t = t + step;
				}
			}
			else
			{
				t = t + step;
			}
			iteration = iteration + 1;
		}
		if (deviation <= eps)
		{
			if (tmp_cur1.p0 == cubBcur1.p0 && tmp_cur1.p1 == cubBcur1.p1)//判断tmp_cur1是不是cubBcur1
			{
				split_t1.push_back(t); //split_t1 = t;//t-> tmp_cur1;
				split_t2.push_back(t_prior);//split_t2 = t_prior;//t_prior-> tmp_cur2;
			}
			else
			{
				split_t1.push_back(t_prior);//split_t1 = t_prior;
				split_t2.push_back(t);//split_t2 = t;
			}
			isIntersect = true;
		}
		start_t = start_t + 0.1;//+0.2的话，当两条线挨得很近话，可能检测不出来交叉点
		//if (start_t > 0.5&&split_t1.empty())
		//	break;
	}
	return isIntersect;

	
}
bool MySvgFile::twoQuaBzrIntersect(quaBezier quaBcur1, quaBezier quaBcur2, vector<double> &split_t1, vector<double> &split_t2)
{
	bool isIntersect = false;
	double start_t = 0;
	while (start_t <= 1)
	{
		double t = start_t;
		double t_prior = 0;
		double step = 0.02;
		quaBezier tmp_cur1 = quaBcur1;
		quaBezier tmp_cur2 = quaBcur2;
		double deviation = 1000;
		int iteration = 0;
		while (t <= 1 && deviation > eps && iteration < 30)
		{
			double a1_cur1, a2_cur1, a3_cur1;
			double b1_cur1, b2_cur1, b3_cur1;
			QuaBezrCoefficient(t, a1_cur1, a2_cur1, a3_cur1);
			QuaBezrDerivative(t, b1_cur1, b2_cur1, b3_cur1);

			double x_cur1 = tmp_cur1.p0.x*a1_cur1 + tmp_cur1.p1.x*a2_cur1 + tmp_cur1.p2.x*a3_cur1;
			double y_cur1 = tmp_cur1.p0.y*a1_cur1 + tmp_cur1.p1.y*a2_cur1 + tmp_cur1.p2.y*a3_cur1;
			double dx = tmp_cur1.p0.x*b1_cur1 + tmp_cur1.p1.x*b2_cur1 + tmp_cur1.p2.x*b3_cur1;
			double dy = tmp_cur1.p0.y*b1_cur1 + tmp_cur1.p1.y*b2_cur1 + tmp_cur1.p2.y*b3_cur1;
			double a = dy / dx;
			double b = -1;
			double c = y_cur1 - a*x_cur1;

			double A = a*(tmp_cur2.p0.x - 2 * tmp_cur2.p1.x + tmp_cur2.p2.x) +
				b*(tmp_cur2.p0.y - 2 * tmp_cur2.p1.y + tmp_cur2.p2.y);
			double B = a*(-2 * tmp_cur2.p0.x + 2 * tmp_cur2.p1.x) +
				b*(-2 * tmp_cur2.p0.y + 2 * tmp_cur2.p1.y);
			double C = a* tmp_cur2.p0.x + b*tmp_cur2.p0.y + c;
			vector<double> t_root = quadricRealPolynomial(A, B, C);
			if (!t_root.empty())
			{
				int suit_t = 0;
				for (unsigned int i = 0; i < t_root.size(); ++i)
				{
					if (t_root[i]>0 && t_root[i]<1)//if (t_root[i]>eps && (1-t_root[i])>eps)//
						suit_t = suit_t + 1;
				}
				bool firstflag = false;
				double tmp_prior_t = 0;
				for (unsigned int i = 0; i < t_root.size(); ++i)
				{
					if (t_root[i]>0 && t_root[i]<1)//if (t_root[i]>eps && (1 - t_root[i])>eps)//
					{
						double a1_cur2, a2_cur2, a3_cur2;
						QuaBezrCoefficient(t_root[i], a1_cur2, a2_cur2, a3_cur2);
						double x_cur2 = tmp_cur2.p0.x*a1_cur2 + tmp_cur2.p1.x*a2_cur2 + tmp_cur2.p2.x*a3_cur2;
						double y_cur2 = tmp_cur2.p0.y*a1_cur2 + tmp_cur2.p1.y*a2_cur2 + tmp_cur2.p2.y*a3_cur2;
						double  tmpd = sqrt(pow(x_cur1 - x_cur2, 2) + pow(y_cur1 - y_cur2, 2));
						if (suit_t > 1)
						{
							if (firstflag == false)
							{
								firstflag = true;
								deviation = tmpd;
								tmp_prior_t = t;
								t_prior = t;
								t = t_root[i];
							}
							else
							{
								if (tmpd < deviation)//considering there are two intersecting points,and select the nearest one
								{
									deviation = tmpd;
									t_prior = tmp_prior_t;
									//t_prior = t;
									t = t_root[i];
								}
							}
						}
						else
						{
							deviation = tmpd;
							t_prior = t;
							t = t_root[i];
						}
					}
				}
				if (suit_t > 0)
				{
					quaBezier tmpcb = tmp_cur1;
					tmp_cur1 = tmp_cur2;
					tmp_cur2 = tmpcb;
				}
				else
					t = t + step;
			}
			else
			{
				t = t + step;
			}
			iteration = iteration + 1;
		}
		if (deviation <= eps)
		{
			if (tmp_cur1.p0 == quaBcur1.p0)
			{
				split_t1.push_back(t);//split_t1 = t;//t-> tmp_cur1;
				split_t2.push_back(t_prior);//split_t2 = t_prior;	//t_prior-> tmp_cur2;
			}
			else
			{
				split_t1.push_back(t_prior);//split_t1 = t_prior;
				split_t2.push_back(t);//split_t2 = t;
			}
			isIntersect = true;
		}
		start_t = start_t + 0.1;
		//if (start_t > 0.5&&split_t1.empty())
		//	break;
	}
	//bool isIntersect=false;
	//double t = 0;
	//double t_prior = 0;
	//double step = 0.1;
	//quaBezier tmp_cur1 = quaBcur1;
	//quaBezier tmp_cur2 = quaBcur2;
	//double deviation = 1000;
	//int iteration = 0;
	//while (t <= 1 && deviation > eps && iteration <= 15)
	//{
	//	double a1_cur1, a2_cur1, a3_cur1;
	//	double b1_cur1, b2_cur1, b3_cur1;
	//	QuaBezrCoefficient(t, a1_cur1, a2_cur1, a3_cur1);
	//	QuaBezrDerivative(t, b1_cur1, b2_cur1, b3_cur1);
	//	double x_cur1 = tmp_cur1.p0.x*a1_cur1 + tmp_cur1.p1.x*a2_cur1 + tmp_cur1.p2.x*a3_cur1;
	//	double y_cur1 = tmp_cur1.p0.y*a1_cur1 + tmp_cur1.p1.y*a2_cur1 + tmp_cur1.p2.y*a3_cur1;
	//	double dx = tmp_cur1.p0.x*b1_cur1 + tmp_cur1.p1.x*b2_cur1 + tmp_cur1.p2.x*b3_cur1;
	//	double dy = tmp_cur1.p0.y*b1_cur1 + tmp_cur1.p1.y*b2_cur1 + tmp_cur1.p2.y*b3_cur1;
	//	double a = dy / dx;
	//	double b = -1;
	//	double c = y_cur1 - a*x_cur1;
	//	double A = a*(tmp_cur2.p0.x - 2 * tmp_cur2.p1.x + tmp_cur2.p2.x) +
	//		b*(tmp_cur2.p0.y - 2 * tmp_cur2.p1.y + tmp_cur2.p2.y);
	//	double B = a*(-2 * tmp_cur2.p0.x + 2 * tmp_cur2.p1.x) +
	//		b*(-2 * tmp_cur2.p0.y + 2 * tmp_cur2.p1.y);
	//	double C = a* tmp_cur2.p0.x + b*tmp_cur2.p0.y + c;
	//	vector<double> t_root = quadricRealPolynomial(A, B, C);
	//	if (!t_root.empty())
	//	{
	//		int suit_t = 0;
	//		for (unsigned int i = 0; i < t_root.size(); ++i)
	//		{
	//			if (t_root[i] > 0 && t_root[i] < 1)
	//				suit_t = suit_t + 1;
	//		}
	//		bool firstflag = false;
	//		for (unsigned int i = 0; i < t_root.size(); ++i)
	//		{
	//			if (t_root[i] > 0 && t_root[i] < 1)
	//			{
	//				double a1_cur2, a2_cur2, a3_cur2;
	//				QuaBezrCoefficient(t_root[i], a1_cur2, a2_cur2, a3_cur2);
	//				double x_cur2 = tmp_cur2.p0.x*a1_cur2 + tmp_cur2.p1.x*a2_cur2 + tmp_cur2.p2.x*a3_cur2;
	//				double y_cur2 = tmp_cur2.p0.y*a1_cur2 + tmp_cur2.p1.y*a2_cur2 + tmp_cur2.p2.y*a3_cur2;
	//				double  tmpd = sqrt(pow(x_cur1 - x_cur2, 2) + pow(y_cur1 - y_cur2, 2));
	//				if (suit_t > 1)
	//				{
	//					if (firstflag == false)
	//					{
	//						firstflag = true;
	//						deviation = tmpd;
	//						t_prior = t;
	//						t = t_root[i];
	//					}
	//					else
	//					{
	//						if (tmpd < deviation)//considering there are two intersecting points,and select the nearest one
	//						{
	//							deviation = tmpd;
	//							t_prior = t;
	//							t = t_root[i];
	//						}
	//					}
	//				}
	//				else
	//				{
	//					deviation = tmpd;
	//					t_prior = t;
	//					t = t_root[i];
	//				}
	//			}
	//		}
	//		if (suit_t > 0)
	//		{
	//			quaBezier tmpcb = tmp_cur1;
	//			tmp_cur1 = tmp_cur2;
	//			tmp_cur2 = tmpcb;
	//		}
	//		else
	//			t = t + step;
	//	}
	//	else
	//	{
	//		t = t + step;
	//	}
	//	iteration = iteration + 1;
	//}
	//if (deviation <= eps)
	//{
	//	if (tmp_cur1.p0 == quaBcur1.p0)
	//	{
	//		split_t1.push_back(t);//split_t1 = t;//t-> tmp_cur1;
	//		split_t2.push_back(t_prior);//split_t2 = t_prior;	//t_prior-> tmp_cur2;
	//	}
	//	else
	//	{
	//		split_t1.push_back(t_prior);//split_t1 = t_prior;
	//		split_t2.push_back(t);//split_t2 = t;
	//	}
	//	isIntersect = true;
	//}
	//return isIntersect;
	return isIntersect;
}
bool MySvgFile::CBandQBIntersect(cubicBezier Cbcur, quaBezier Qbcur, vector<double> &split_t1, vector<double> &split_t2)
{
	bool isIntersect = false;
	double start_t = 0;
	while (start_t <= 1)
	{
		double t = start_t;
		double t_prior = 0;
		double step = 0.02;
		int changflag = 0;
		cubicBezier tmp_cur1 = Cbcur;
		quaBezier tmp_cur2 = Qbcur;
		double deviation = 1000;
		int iteration = 0;
		while (t <= 1 && deviation > eps && iteration < 30)
		{
			if (changflag == 0)
			{
				double a1_cur1, a2_cur1, a3_cur1, a4_cur1;
				double b1_cur1, b2_cur1, b3_cur1, b4_cur1;
				CubicBezrCoefficient(t, a1_cur1, a2_cur1, a3_cur1, a4_cur1);
				CubicBezrDerivative(t, b1_cur1, b2_cur1, b3_cur1, b4_cur1);

				double x_cur1 = tmp_cur1.p0.x*a1_cur1 + tmp_cur1.p1.x*a2_cur1 + tmp_cur1.p2.x*a3_cur1 + tmp_cur1.p3.x*a4_cur1;
				double y_cur1 = tmp_cur1.p0.y*a1_cur1 + tmp_cur1.p1.y*a2_cur1 + tmp_cur1.p2.y*a3_cur1 + tmp_cur1.p3.y*a4_cur1;
				double dx = tmp_cur1.p0.x*b1_cur1 + tmp_cur1.p1.x*b2_cur1 + tmp_cur1.p2.x*b3_cur1 + tmp_cur1.p3.x*b4_cur1;
				double dy = tmp_cur1.p0.y*b1_cur1 + tmp_cur1.p1.y*b2_cur1 + tmp_cur1.p2.y*b3_cur1 + tmp_cur1.p3.y*b4_cur1;
				double a = dy / dx;
				double b = -1;
				double c = y_cur1 - a*x_cur1;

				double A = a*(tmp_cur2.p0.x - 2 * tmp_cur2.p1.x + tmp_cur2.p2.x) +
					b*(tmp_cur2.p0.y - 2 * tmp_cur2.p1.y + tmp_cur2.p2.y);
				double B = a*(-2 * tmp_cur2.p0.x + 2 * tmp_cur2.p1.x) +
					b*(-2 * tmp_cur2.p0.y + 2 * tmp_cur2.p1.y);
				double C = a* tmp_cur2.p0.x + b*tmp_cur2.p0.y + c;
				vector<double> t_root = quadricRealPolynomial(A, B, C);
				if (!t_root.empty())
				{
					int suit_t = 0;
					for (unsigned int i = 0; i < t_root.size(); ++i)
					{
						if (t_root[i]>0 && t_root[i]<1)//if (t_root[i]>eps && (1 - t_root[i])>eps)//
						{
							suit_t = suit_t + 1;
						}
					}
					bool first_flag = false;
					double tmp_prior_t = 0;
					for (unsigned int i = 0; i < t_root.size(); ++i)
					{
						if (t_root[i]>0 && t_root[i]<1)//if (t_root[i]>eps && (1 - t_root[i])>eps)//
						{
							double a1_cur2, a2_cur2, a3_cur2;
							QuaBezrCoefficient(t_root[i], a1_cur2, a2_cur2, a3_cur2);
							double x_cur2 = tmp_cur2.p0.x*a1_cur2 + tmp_cur2.p1.x*a2_cur2 + tmp_cur2.p2.x*a3_cur2;
							double y_cur2 = tmp_cur2.p0.y*a1_cur2 + tmp_cur2.p1.y*a2_cur2 + tmp_cur2.p2.y*a3_cur2;
							double  tmpd = sqrt(pow(x_cur1 - x_cur2, 2) + pow(y_cur1 - y_cur2, 2));
							if (suit_t > 1)
							{
								if (first_flag == false)
								{
									first_flag = true;
									deviation = tmpd;
									tmp_prior_t = t;
									t_prior = t;
									t = t_root[i];
								}
								else
								{
									if (tmpd < deviation)//considering there are two intersecting points,and select the nearest one
									{
										deviation = tmpd;
										t_prior = tmp_prior_t;
										//t_prior = t;
										t = t_root[i];
									}
								}
							}
							else
							{
								deviation = tmpd;
								t_prior = t;
								t = t_root[i];
							}
						}
					}
					if (suit_t > 0)
					{
						changflag = 1;
					}
					else
						t = t + step;
				}
				else
				{
					t = t + step;
				}
			}
			else
			{
				double a1_cur2, a2_cur2, a3_cur2;
				double b1_cur2, b2_cur2, b3_cur2;
				QuaBezrCoefficient(t, a1_cur2, a2_cur2, a3_cur2);
				QuaBezrDerivative(t, b1_cur2, b2_cur2, b3_cur2);

				double x_cur2 = tmp_cur2.p0.x*a1_cur2 + tmp_cur2.p1.x*a2_cur2 + tmp_cur2.p2.x*a3_cur2;
				double y_cur2 = tmp_cur2.p0.y*a1_cur2 + tmp_cur2.p1.y*a2_cur2 + tmp_cur2.p2.y*a3_cur2;
				double dx = tmp_cur2.p0.x*b1_cur2 + tmp_cur2.p1.x*b2_cur2 + tmp_cur2.p2.x*b3_cur2;
				double dy = tmp_cur2.p0.y*b1_cur2 + tmp_cur2.p1.y*b2_cur2 + tmp_cur2.p2.y*b3_cur2;
				double a = dy / dx;
				double b = -1;
				double c = y_cur2 - a*x_cur2;

				double A = a*(-tmp_cur1.p0.x + 3 * tmp_cur1.p1.x - 3 * tmp_cur1.p2.x + tmp_cur1.p3.x) +
					b*(-tmp_cur1.p0.y + 3 * tmp_cur1.p1.y - 3 * tmp_cur1.p2.y + tmp_cur1.p3.y);
				double B = a*(3 * tmp_cur1.p0.x - 6 * tmp_cur1.p1.x + 3 * tmp_cur1.p2.x) +
					b*(3 * tmp_cur1.p0.y - 6 * tmp_cur1.p1.y + 3 * tmp_cur1.p2.y);
				double C = a*(-3 * tmp_cur1.p0.x + 3 * tmp_cur1.p1.x) +
					b*(-3 * tmp_cur1.p0.y + 3 * tmp_cur1.p1.y);
				double D = a*tmp_cur1.p0.x + b*tmp_cur1.p0.y + c;
				vector<double> t_root = cubicRealPolynomial(A, B, C, D);
				if (!t_root.empty())
				{
					int suit_t = 0;
					for (unsigned int i = 0; i < t_root.size(); ++i)
					{
						if (t_root[i]>0 && t_root[i]<1)//if (t_root[i]>eps && (1 - t_root[i])>eps)//
						{
							suit_t = suit_t + 1;
						}
					}
					bool first_flag = false;
					double tmp_prior_t = 0;
					for (unsigned int i = 0; i < t_root.size(); ++i)
					{
						if (t_root[i]>0 && t_root[i]<1)//if (t_root[i]>eps && (1 - t_root[i])>eps)//
						{
							double a1_cur1, a2_cur1, a3_cur1, a4_cur1;
							CubicBezrCoefficient(t_root[i], a1_cur1, a2_cur1, a3_cur1, a4_cur1);
							double x_cur1 = tmp_cur1.p0.x*a1_cur1 + tmp_cur1.p1.x*a2_cur1 + tmp_cur1.p2.x*a3_cur1 + tmp_cur1.p3.x*a4_cur1;
							double y_cur1 = tmp_cur1.p0.y*a1_cur1 + tmp_cur1.p1.y*a2_cur1 + tmp_cur1.p2.y*a3_cur1 + tmp_cur1.p3.y*a4_cur1;
							double  tmpd = sqrt(pow(x_cur1 - x_cur2, 2) + pow(y_cur1 - y_cur2, 2));
							if (suit_t > 1)
							{
								if (first_flag == false)
								{
									first_flag = true;
									deviation = tmpd;
									tmp_prior_t = t;
									t_prior = t;
									t = t_root[i];
								}
								else
								{
									if (tmpd < deviation)//considering there are two intersecting points,and select the nearest one
									{
										deviation = tmpd;
										t_prior = tmp_prior_t;
										//t_prior = t;
										t = t_root[i];
									}
								}
							}
							else
							{
								deviation = tmpd;
								t_prior = t;
								t = t_root[i];
							}
						}
					}
					if (suit_t > 0)
					{
						changflag = 0;
					}
					else
						t = t + step;
				}
				else
				{
					t = t + step;
				}
			}
			iteration = iteration + 1;
		}
		if (deviation <= eps)
		{
			if (changflag == 0)
			{
				split_t1.push_back(t);//split_t1 = t; //splitCubBzrs(t, tmp_cur1); //t-> tmp_cur1;
				split_t2.push_back(t_prior);//split_t2 = t_prior; //splitQuaBzrs(t_prior, tmp_cur2);//t_prior-> tmp_cur2;		
			}
			else
			{
				split_t1.push_back(t_prior);//split_t1 = t_prior; //splitCubBzrs(t_prior, tmp_cur1);//t_prior-> tmp_cur1;
				split_t2.push_back(t);//split_t2 = t; //splitQuaBzrs(t, tmp_cur2);//t-> tmp_cur2;			
			}
			isIntersect = true;
		}
		start_t = start_t + 0.1;
		//if (start_t > 0.5&&split_t1.empty())
		//	break;
	}	
	return isIntersect;
}
bool MySvgFile::CBandPLIntersect(cubicBezier Cbcur, vector<mypoint2f> polyline, vector<double> &split_t1, vector<mypoint2f> &split_p)
{
	bool flag = false;
	mypoint2f p1 = polyline.front();
	mypoint2f p2 = polyline.back();
	double a = p2.y - p1.y;
	double b = p1.x - p2.x;
	double c = p2.x*p1.y - p1.x*p2.y;

	double A = a*(-Cbcur.p0.x + 3 * Cbcur.p1.x - 3 * Cbcur.p2.x + Cbcur.p3.x) +
		b*(-Cbcur.p0.y + 3 * Cbcur.p1.y - 3 * Cbcur.p2.y + Cbcur.p3.y);
	double B = a*(3 * Cbcur.p0.x - 6 * Cbcur.p1.x + 3 * Cbcur.p2.x) +
		b*(3 * Cbcur.p0.y - 6 * Cbcur.p1.y + 3 * Cbcur.p2.y);
	double C = a*(-3 * Cbcur.p0.x + 3 * Cbcur.p1.x) +
		b*(-3 * Cbcur.p0.y + 3 * Cbcur.p1.y);
	double D = a*Cbcur.p0.x + b*Cbcur.p0.y + c;
	vector<double> t_root = cubicRealPolynomial(A, B, C, D);
	if (!t_root.empty())
	{
		for (unsigned int i = 0; i < t_root.size(); ++i)
		{
			if (t_root[i]>0 && t_root[i]<1)//if (t_root[i]>eps && (1 - t_root[i])>eps)
			{
				double a1_cur, a2_cur, a3_cur, a4_cur;
				CubicBezrCoefficient(t_root[i], a1_cur, a2_cur, a3_cur, a4_cur);
				double x_cur = Cbcur.p0.x*a1_cur + Cbcur.p1.x*a2_cur + Cbcur.p2.x*a3_cur + Cbcur.p3.x*a4_cur;
				double y_cur = Cbcur.p0.y*a1_cur + Cbcur.p1.y*a2_cur + Cbcur.p2.y*a3_cur + Cbcur.p3.y*a4_cur;
				mypoint2f iterp(x_cur, y_cur);				
				if (isOnSegment(iterp, p1, p2))//&& LineSegmentLength(iterp, Cbcur.p0)>eps &&LineSegmentLength(iterp, Cbcur.p3)>eps)
				{
					split_t1.push_back(t_root[i]);
					split_p.push_back(iterp);
					flag = true;
				}
			}
		}
		return flag;
	}
	else
		return flag;
}
bool MySvgFile::QBandPLIntersect(quaBezier qbcur, vector<mypoint2f> polyline, vector<double> &split_t1, vector<mypoint2f> &split_p)
{
	bool flag = false;
	mypoint2f p1 = polyline.front();
	mypoint2f p2 = polyline.back();
	double a = p2.y - p1.y;
	double b = p1.x - p2.x;
	double c = p2.x*p1.y - p1.x*p2.y;

	double A = a*(qbcur.p0.x - 2 * qbcur.p1.x + qbcur.p2.x) +
		b*(qbcur.p0.y - 2 * qbcur.p1.y + qbcur.p2.y);
	double B = a*(-2 * qbcur.p0.x + 2 * qbcur.p1.x) +
		b*(-2 * qbcur.p0.y + 2 * qbcur.p1.y);
	double C = a* qbcur.p0.x + b*qbcur.p0.y + c;
	vector<double> t_root = quadricRealPolynomial(A, B, C);
	if (!t_root.empty())
	{
		for (unsigned int i = 0; i < t_root.size(); ++i)
		{		
			if (t_root[i]>0 && t_root[i]<1)//if (t_root[i]>eps && (1 - t_root[i])>eps)
			{
				double a1_cur2, a2_cur2, a3_cur2;
				QuaBezrCoefficient(t_root[i], a1_cur2, a2_cur2, a3_cur2);
				double x_cur2 = qbcur.p0.x*a1_cur2 + qbcur.p1.x*a2_cur2 + qbcur.p2.x*a3_cur2;
				double y_cur2 = qbcur.p0.y*a1_cur2 + qbcur.p1.y*a2_cur2 + qbcur.p2.y*a3_cur2;
				mypoint2f iterp(x_cur2, y_cur2);
				if (isOnSegment(iterp, p1, p2) )//&& LineSegmentLength(iterp, qbcur.p0)>eps &&LineSegmentLength(iterp, qbcur.p2)>eps)
				{
					split_t1.push_back(t_root[i]);
					split_p.push_back(iterp);
					flag = true;
				}
			}
		}
		return flag;
	}
	else
		return flag;
}
bool MySvgFile::twoPolyIntersect(vector<mypoint2f> pline1, vector<mypoint2f> pline2, vector<mypoint2f> &split_p1, vector<mypoint2f> &split_p2)
{
	mypoint2f p1 = pline1.front();
	mypoint2f p2 = pline1.back();
	mypoint2f p3 = pline2.front();
	mypoint2f p4 = pline2.back();
	mypoint2f inpt = getIntersectPoint(p1, p2, p3, p4);
	if (isOnSegment(inpt, p1, p2) && isOnSegment(inpt, p3, p4))
	{
		split_p1.push_back(inpt);//split_p = inpt;
		split_p2.push_back(inpt);
		return true;
	}
	else
	{
		return false;
	}	
}
//并行的处理计算相交
bool MySvgFile::twoCubicBzrIntersect_para(int para_i, cubicBezier cubBcur1, cubicBezier cubBcur2, vector<double> &split_t1, vector<double> &split_t2)
{
	bool isIntersect = false;
	double start_t = 0 + 0.1*(double)para_i;  //以start_t为cub1的起始点 寻找cross_t
	//while (start_t <= 1)
	//{
		bool first_t = false;////////////////////////////////////
		double t = start_t;
		double t_prior = 0;
		double step = 0.02;
		cubicBezier tmp_cur1 = cubBcur1;
		cubicBezier tmp_cur2 = cubBcur2;
		double deviation = 1000;
		int iteration = 0;
		while (t <= 1 && t >= 0 && deviation > eps && iteration < 30)
		{
			double a1_cur1, a2_cur1, a3_cur1, a4_cur1;
			double b1_cur1, b2_cur1, b3_cur1, b4_cur1;
			CubicBezrCoefficient(t, a1_cur1, a2_cur1, a3_cur1, a4_cur1);
			CubicBezrDerivative(t, b1_cur1, b2_cur1, b3_cur1, b4_cur1);

			double x_cur1 = tmp_cur1.p0.x*a1_cur1 + tmp_cur1.p1.x*a2_cur1 + tmp_cur1.p2.x*a3_cur1 + tmp_cur1.p3.x*a4_cur1;
			double y_cur1 = tmp_cur1.p0.y*a1_cur1 + tmp_cur1.p1.y*a2_cur1 + tmp_cur1.p2.y*a3_cur1 + tmp_cur1.p3.y*a4_cur1;
			double dx = tmp_cur1.p0.x*b1_cur1 + tmp_cur1.p1.x*b2_cur1 + tmp_cur1.p2.x*b3_cur1 + tmp_cur1.p3.x*b4_cur1;
			double dy = tmp_cur1.p0.y*b1_cur1 + tmp_cur1.p1.y*b2_cur1 + tmp_cur1.p2.y*b3_cur1 + tmp_cur1.p3.y*b4_cur1;
			double a = dy / dx;
			double b = -1;
			double c = y_cur1 - a*x_cur1;

			double A = a*(-tmp_cur2.p0.x + 3 * tmp_cur2.p1.x - 3 * tmp_cur2.p2.x + tmp_cur2.p3.x) +
				b*(-tmp_cur2.p0.y + 3 * tmp_cur2.p1.y - 3 * tmp_cur2.p2.y + tmp_cur2.p3.y);
			double B = a*(3 * tmp_cur2.p0.x - 6 * tmp_cur2.p1.x + 3 * tmp_cur2.p2.x) +
				b*(3 * tmp_cur2.p0.y - 6 * tmp_cur2.p1.y + 3 * tmp_cur2.p2.y);
			double C = a*(-3 * tmp_cur2.p0.x + 3 * tmp_cur2.p1.x) +
				b*(-3 * tmp_cur2.p0.y + 3 * tmp_cur2.p1.y);
			double D = a*tmp_cur2.p0.x + b*tmp_cur2.p0.y + c;
			vector<double> t_root = cubicRealPolynomial(A, B, C, D);
			if (!t_root.empty())//t_root求出来的根是tmp_cur2上的t值
			{
				int suit_t = 0;
				for (unsigned int i = 0; i < t_root.size(); ++i)
				{
					if (t_root[i]>0 && t_root[i]<1)//if (t_root[i]>eps && (1 - t_root[i])>eps)//
						suit_t = suit_t + 1;
				}
				bool firstflag = false;
				double tmp_prior_t = 0;
				for (unsigned int i = 0; i < t_root.size(); ++i)
				{
					if (t_root[i]>0 && t_root[i]<1)//if (t_root[i]>eps && (1 - t_root[i])>eps)//
					{
						if (first_t == false && t_root[i]>0.95)
						{
							step = -0.02;
							first_t = true;
						}
						double a1_cur2, a2_cur2, a3_cur2, a4_cur2;
						CubicBezrCoefficient(t_root[i], a1_cur2, a2_cur2, a3_cur2, a4_cur2);
						double x_cur2 = tmp_cur2.p0.x*a1_cur2 + tmp_cur2.p1.x*a2_cur2 + tmp_cur2.p2.x*a3_cur2 + tmp_cur2.p3.x*a4_cur2;
						double y_cur2 = tmp_cur2.p0.y*a1_cur2 + tmp_cur2.p1.y*a2_cur2 + tmp_cur2.p2.y*a3_cur2 + tmp_cur2.p3.y*a4_cur2;
						double  tmpd = sqrt(pow(x_cur1 - x_cur2, 2) + pow(y_cur1 - y_cur2, 2));
						if (suit_t > 1)
						{
							if (firstflag == false)
							{
								firstflag = true;
								deviation = tmpd;
								tmp_prior_t = t;
								t_prior = t;
								t = t_root[i];
							}
							else
							{
								if (tmpd < deviation)//considering there are two more intersecting points,and select the nearest one
								{
									deviation = tmpd;
									t_prior = tmp_prior_t;
									//t_prior = t;
									t = t_root[i];
								}
							}
						}
						else
						{
							deviation = tmpd;
							t_prior = t;
							t = t_root[i];
						}
					}
				}
				if (suit_t > 0)
				{
					cubicBezier tmpcb = tmp_cur1;
					tmp_cur1 = tmp_cur2;
					tmp_cur2 = tmpcb;
				}
				else
				{
					t = t + step;
				}
			}
			else
			{
				t = t + step;
			}
			iteration = iteration + 1;
		}
		if (deviation <= eps)
		{
			if (tmp_cur1.p0 == cubBcur1.p0 && tmp_cur1.p1 == cubBcur1.p1)//判断tmp_cur1是不是cubBcur1
			{
				tbb::mutex mu;
				mu.lock();
				split_t1.push_back(t); //split_t1 = t;//t-> tmp_cur1;
				split_t2.push_back(t_prior);//split_t2 = t_prior;//t_prior-> tmp_cur2;
				mu.unlock();
			}
			else
			{
				tbb::mutex mu;
				mu.lock();
				split_t1.push_back(t_prior);//split_t1 = t_prior;
				split_t2.push_back(t);//split_t2 = t;
				mu.unlock();
			}
			isIntersect = true;
		}
		//start_t = start_t + 0.1;//+0.2的话，当两条线挨得很近话，可能检测不出来交叉点
	//}
	return isIntersect;
}


void MySvgFile::splitACurve(DiffCurve curve, vector<double> &split_t, vector<mypoint2f> &split_p, int density)
{
	//delete repeated split points
	if (split_t.size()>1)
	{
		vector<double>::iterator iter = split_t.begin();
		while (iter != split_t.end())
		{
			vector<double>::iterator iter_inn = iter + 1;
			while (iter_inn != split_t.end())
			{
				if (fabs(*iter - *iter_inn)<eps)//if (*iter == *iter_inn),(约等于)
					iter_inn=split_t.erase(iter_inn);
				else
					iter_inn++;
			}
			++iter;
		}
	}
	if (split_p.size()>1)
	{
		vector<mypoint2f>::iterator iter = split_p.begin();
		while (iter != split_p.end())
		{
			vector<mypoint2f>::iterator iter_inn = iter + 1;
			while (iter_inn != split_p.end())
			{
				if (LineSegmentLength(*iter, *iter_inn)<eps)//if (isPointRoughEqual(*iter,*iter_inn))//if (*iter == *iter_inn)
					iter_inn = split_p.erase(iter_inn);
				else
					iter_inn++;
			}
			++iter;
		}	
	}
	if (curve.curveType == 'c'&& !split_t.empty())
	{
		splitCubBzrs(split_t, cubicBezier_vec_origin[curve.curveIndex], curve.curveIndex, this->added_polyline, density);
	}
	else if (curve.curveType == 'q' && !split_t.empty())
	{
		splitQuaBzrs(split_t, quadraticBezier_vec_origin[curve.curveIndex], curve.curveIndex, this->added_polyline, density);
	}
	else if (curve.curveType == 'l' && !split_p.empty())
	{
		splitPolyline(split_p, polyline_vec_origin[curve.curveIndex].pointlist, curve.curveIndex, this->added_polyline);
	}
}
void MySvgFile::splitCubBzrs(vector<double> split_t, cubicBezier cbcurve, int origi_index, vector<PolyLine> &add_polyline, int density)
{
 	double p1[2] = { cbcurve.p0.x, cbcurve.p0.y };
	double p2[2] = { cbcurve.p1.x, cbcurve.p1.y };
	double p3[2] = { cbcurve.p2.x, cbcurve.p2.y };
	double p4[2] = { cbcurve.p3.x, cbcurve.p3.y };	
	double totallong = 0;
	totallong = LineSegmentLength(cbcurve.p0, cbcurve.p1) + LineSegmentLength(cbcurve.p1, cbcurve.p2) + LineSegmentLength(cbcurve.p2, cbcurve.p3);
	
	sort(split_t.begin(), split_t.end());
	if (split_t.back()!=1)
		split_t.push_back(1);
	if (split_t.front() == 0)
		split_t.erase(split_t.begin());
	double start_t = 0;
	PolyLine apl;//vector<mypoint2f> seg_curve1;
	double a1, a2, a3, a4;
	mypoint2f apoint;
	if (density == 0)
	{
		double step = 0.02;
		//三次公式 B(t)=p0*(1-t)^3 + p1*3t(1-t)^2 + p2*3t^2(1-t) + p3*t^3
		for (unsigned int i = 0; i < split_t.size(); ++i)
		{
			step = 0.02;
			if (split_t[i] - start_t < step)
				step = abs(split_t[i] - start_t);
			else
				step = 0.02;
			for (double t = start_t; t < split_t[i]; t += step)
			{
				a1 = pow((1 - t), 3);
				a2 = pow((1 - t), 2) * 3 * t;
				a3 = 3 * t*t*(1 - t);
				a4 = t*t*t;
				apoint.x = a1*p1[0] + a2*p2[0] + a3*p3[0] + a4*p4[0];
				apoint.y = a1*p1[1] + a2*p2[1] + a3*p3[1] + a4*p4[1];
				apl.pointlist.push_back(apoint);
				//mycircle acir;
				//acir.cx = apoint.x;
				//acir.cy = apoint.y;
				//acir.radius = 0.2;
				//circle_vec.push_back(acir);
			}
			a1 = pow((1 - split_t[i]), 3);
			a2 = pow((1 - split_t[i]), 2) * 3 * split_t[i];
			a3 = 3 * pow(split_t[i], 2)*(1 - split_t[i]);
			a4 = pow(split_t[i], 3);
			apoint.x = a1*p1[0] + a2*p2[0] + a3*p3[0] + a4*p4[0];
			apoint.y = a1*p1[1] + a2*p2[1] + a3*p3[1] + a4*p4[1];
			if (fabs(apoint.x - apl.pointlist.back().x)>eps * 10 || fabs(apoint.y - apl.pointlist.back().y)>eps * 10)
				apl.pointlist.push_back(apoint);
			else
			{
				apl.pointlist.pop_back();
				apl.pointlist.push_back(apoint);
			}
			apl.origi_index = origi_index;
			apl.origi_type = 'c';
			add_polyline.push_back(apl);
			apl.pointlist.clear();
			apl.pointlist.swap(vector<mypoint2f>());
			start_t = split_t[i];
		}
	}
	else if (density == 1)
	{
		double step = 0;
		//三次公式 B(t)=p0*(1-t)^3 + p1*3t(1-t)^2 + p2*3t^2(1-t) + p3*t^3
		for (unsigned int i = 0; i < split_t.size(); ++i)
		{
			int diviNum = (split_t[i] - start_t)*totallong /0.8;//0.8
			if (diviNum < 1)
				step = split_t[i] - start_t;
			else
				step = (split_t[i] - start_t) / diviNum;

			for (double t = start_t; t < split_t[i]; t += step)
			{
				a1 = pow((1 - t), 3);
				a2 = pow((1 - t), 2) * 3 * t;
				a3 = 3 * t*t*(1 - t);
				a4 = t*t*t;
				apoint.x = a1*p1[0] + a2*p2[0] + a3*p3[0] + a4*p4[0];
				apoint.y = a1*p1[1] + a2*p2[1] + a3*p3[1] + a4*p4[1];
				apl.pointlist.push_back(apoint);
				//mycircle acir;
				//acir.cx = apoint.x;
				//acir.cy = apoint.y;
				//acir.radius = 0.2;
				//circle_vec.push_back(acir);
			}
			a1 = pow((1 - split_t[i]), 3);
			a2 = pow((1 - split_t[i]), 2) * 3 * split_t[i];
			a3 = 3 * pow(split_t[i], 2)*(1 - split_t[i]);
			a4 = pow(split_t[i], 3);
			apoint.x = a1*p1[0] + a2*p2[0] + a3*p3[0] + a4*p4[0];
			apoint.y = a1*p1[1] + a2*p2[1] + a3*p3[1] + a4*p4[1];
			if (fabs(apoint.x - apl.pointlist.back().x)>eps * 10 || fabs(apoint.y - apl.pointlist.back().y)>eps * 10)
				apl.pointlist.push_back(apoint);
			else
			{
				apl.pointlist.pop_back();
				apl.pointlist.push_back(apoint);
			}
			apl.origi_index = origi_index;
			apl.origi_type = 'c';
			add_polyline.push_back(apl);
			apl.pointlist.clear();
			apl.pointlist.swap(vector<mypoint2f>());
			start_t = split_t[i];
		}
	}
}
void MySvgFile::splitQuaBzrs(vector<double> split_t, quaBezier qbcurve, int origi_index, vector<PolyLine> &add_polyline, int density)
{
	double p1[2] = { qbcurve.p0.x, qbcurve.p0.y };
	double p2[2] = { qbcurve.p1.x, qbcurve.p1.y };
	double p3[2] = { qbcurve.p2.x, qbcurve.p2.y };
	double totallong = 0;
	totallong = LineSegmentLength(qbcurve.p0, qbcurve.p1) + LineSegmentLength(qbcurve.p1, qbcurve.p2);
	sort(split_t.begin(), split_t.end());
	if (split_t.back() != 1)
		split_t.push_back(1);
	if (split_t.front() == 0)
		split_t.erase(split_t.begin());
	double start_t = 0;
	PolyLine apl;
	double a1, a2, a3;
	mypoint2f apoint;
	if (density == 0)
	{
		double step = 0.02;
		for (unsigned int i = 0; i < split_t.size(); ++i)
		{
			step = 0.02;
			if (split_t[i] - start_t < step)
				step = abs(split_t[i] - start_t);
			else
				step = 0.02;
			for (double t = start_t; t < split_t[i]; t += step)
			{
				a1 = pow((1 - t), 2);
				a2 = 2 * t*(1 - t);
				a3 = t*t;
				apoint.x = a1*p1[0] + a2*p2[0] + a3*p3[0];
				apoint.y = a1*p1[1] + a2*p2[1] + a3*p3[1];
				apl.pointlist.push_back(apoint);
			}
			a1 = pow((1 - split_t[i]), 2);
			a2 = 2 * split_t[i] * (1 - split_t[i]);
			a3 = pow(split_t[i], 2);
			apoint.x = a1*p1[0] + a2*p2[0] + a3*p3[0];
			apoint.y = a1*p1[1] + a2*p2[1] + a3*p3[1];
			if (fabs(apoint.x - apl.pointlist.back().x)>eps || fabs(apoint.y - apl.pointlist.back().y)>eps)
				apl.pointlist.push_back(apoint);
			else
			{
				apl.pointlist.pop_back();
				apl.pointlist.push_back(apoint);
			}
			apl.origi_index = origi_index;
			apl.origi_type = 'q';
			add_polyline.push_back(apl);
			apl.pointlist.swap(vector<mypoint2f>());
			start_t = split_t[i];
		}
	}
	else if (density == 1)
	{
		double step = 0;
		for (unsigned int i = 0; i < split_t.size(); ++i)
		{
			int diviNum = (split_t[i] - start_t)*totallong / 0.8;//0.8
			if (diviNum < 1)
				step = split_t[i] - start_t;
			else
				step = (split_t[i] - start_t) / diviNum;
			for (double t = start_t; t < split_t[i]; t += step)
			{
				a1 = pow((1 - t), 2);
				a2 = 2 * t*(1 - t);
				a3 = t*t;
				apoint.x = a1*p1[0] + a2*p2[0] + a3*p3[0];
				apoint.y = a1*p1[1] + a2*p2[1] + a3*p3[1];
				apl.pointlist.push_back(apoint);
			}
			a1 = pow((1 - split_t[i]), 2);
			a2 = 2 * split_t[i] * (1 - split_t[i]);
			a3 = pow(split_t[i], 2);
			apoint.x = a1*p1[0] + a2*p2[0] + a3*p3[0];
			apoint.y = a1*p1[1] + a2*p2[1] + a3*p3[1];
			if (fabs(apoint.x - apl.pointlist.back().x)>eps || fabs(apoint.y - apl.pointlist.back().y)>eps)
				apl.pointlist.push_back(apoint);
			else
			{
				apl.pointlist.pop_back();
				apl.pointlist.push_back(apoint);
			}
			apl.origi_index = origi_index;
			apl.origi_type = 'q';
			add_polyline.push_back(apl);
			apl.pointlist.swap(vector<mypoint2f>());
			start_t = split_t[i];
		}
	}
}
void MySvgFile::splitPolyline(vector<mypoint2f> splitps, vector<mypoint2f> origiline, int origi_index, vector<PolyLine> &add_polyline)
{
	mypoint2f startp, endp;
	if (fabs(origiline.front().x - origiline.back().x) > fabs(origiline.front().y - origiline.back().y))
	{
		if (origiline.front().x < origiline.back().x)
		{
			startp = origiline.front();
			endp = origiline.back();
		}
		else
		{
			startp = origiline.back();
			endp = origiline.front();
		}	
		splitps.push_back(endp);
		for (unsigned int i = 0; i < splitps.size(); ++i)
		{
			for (unsigned int j = i; j < splitps.size(); ++j)
			{
				if (splitps[i].x > splitps[j].x)
				{
					mypoint2f tmp = splitps[i];
					splitps[i] = splitps[j];
					splitps[j] = tmp;
				}
			}
		}
	}
	else
	{
		if (origiline.front().y < origiline.back().y)
		{
			startp = origiline.front();
			endp = origiline.back();
		}
		else
		{
			startp = origiline.back();
			endp = origiline.front();
		}
		splitps.push_back(endp);
		for (unsigned int i = 0; i < splitps.size(); ++i)
		{
			for (unsigned int j = i; j < splitps.size(); ++j)
			{
				if (splitps[i].y > splitps[j].y)
				{
					mypoint2f tmp = splitps[i];
					splitps[i] = splitps[j];
					splitps[j] = tmp;
				}
			}
		}
	}
	//splitps.push_back(endp);
	//for (unsigned int i = 0; i < splitps.size(); ++i)
	//{
	//	for (unsigned int j = i; j < splitps.size(); ++j)
	//	{
	//		if (splitps[i].x > splitps[j].x)
	//		{
	//			mypoint2f tmp = splitps[i];
	//			splitps[i] = splitps[j];
	//			splitps[j] = tmp;
	//		}
	//	}
	//}
	
	for (unsigned int i = 0; i < splitps.size(); ++i)
	{
		//vector<mypoint2f> tp;
		PolyLine apl;
		double sub_len = LineSegmentLength(startp, splitps[i]);
		int sub_num = sub_len / 0.8;
		if (sub_num < 2)
			sub_num = 2;
		mypoint2f direct_vec = splitps[i] - startp;
		direct_vec.x = direct_vec.x / sub_num;
		direct_vec.y = direct_vec.y / sub_num;
		apl.pointlist.push_back(startp);
		mypoint2f tmp_pt = startp;
		for (int j = 0; j < sub_num; ++j)
		{
			tmp_pt = tmp_pt + direct_vec;
			apl.pointlist.push_back(tmp_pt);
		}
		apl.pointlist.pop_back();
		apl.pointlist.push_back(splitps[i]);
		//apl.pointlist.push_back(startp);
		//apl.pointlist.push_back(splitps[i]);
		apl.origi_index = origi_index;
		apl.origi_type = 'l';
		add_polyline.push_back(apl);
		startp = splitps[i];
	}
}

void MySvgFile::writefile()
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
	rootele->SetAttribute("width", svg_width);
	rootele->SetAttribute("height", svg_height);
	rootele->SetAttribute("viewBox", svg_view.c_str());
	rootele->SetAttribute("xml:space", "preserve");
	writeDoc->LinkEndChild(rootele);
	//一个path是一个多边形面
	//例子：
	//注意 fill-rule属性：nonzero, evenodd。 style = "fill:white;stroke:red;stroke-width:2"
	/* <path d="M 100 100 L 300 100 L 200 300 z" fill="red" stroke="blue" stroke-width="3" />*/
	/*<polygon points="220,100 300,210 170,250" style="fill:#cccccc;stroke:#000000;stroke-width:1"/>*/
	vector<char*> color_vec;
	char *c0 = "787878"; color_vec.push_back(c0);
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

	vector<vector<mypoint2f>> plist_diffcurve;
	for (unsigned int i = 0; i < cubicBezier_vec_origin.size(); ++i)
	{
		vector<mypoint2f> cubplist;
		double step = 0.02;
		cubicBezier tmpcub = cubicBezier_vec_origin[i];
		double t = 0;
		for (t = 0; t < 1; t += step)
		{
			mypoint2f apoint(0, 0);
			double a1 = pow((1 - t), 3);
			double a2 = pow((1 - t), 2) * 3 * t;
			double a3 = 3 * t*t*(1 - t);
			double a4 = t*t*t;
			apoint.x = a1*tmpcub.p0.x + a2*tmpcub.p1.x + a3*tmpcub.p2.x + a4*tmpcub.p3.x;
			apoint.y = a1*tmpcub.p0.y + a2*tmpcub.p1.y + a3*tmpcub.p2.y + a4*tmpcub.p3.y;
			cubplist.push_back(apoint);
		}
		if (cubplist.back() != tmpcub.p3)
			cubplist.push_back(tmpcub.p3);
		plist_diffcurve.push_back(cubplist);
	}
	for (unsigned int i = 0; i < quadraticBezier_vec_origin.size(); ++i)
	{
		quaBezier tmpqua = quadraticBezier_vec_origin[i];
		vector < mypoint2f> quaplist;
		double step = 0.02;
		for (double t = 0.0; t <= 1.0; t += step)
		{
			mypoint2f pt;
			GLfloat a1 = pow((1 - t), 2);
			GLfloat a2 = 2 * t*(1 - t);
			GLfloat a3 = t*t;
			pt.x = a1*tmpqua.p0.x + a2*tmpqua.p1.x + a3*tmpqua.p2.x;
			pt.y= a1*tmpqua.p0.y + a2*tmpqua.p1.y + a3*tmpqua.p2.y;
			quaplist.push_back(pt);
		}
		plist_diffcurve.push_back(quaplist);
	}
	for (unsigned int i = 0; i < polyline_vec_origin.size(); ++i)
	{
		PolyLine tmppl = polyline_vec_origin[i];
		plist_diffcurve.push_back(tmppl.pointlist);
	}
	vector<vector<mypoint2f>> new_cubplist;
	for (unsigned int i = 0; i < piecewise_Bezier.size(); ++i)
	{
		vector<mypoint2f> cubplist;
		double step = 0.02;
		cubicBezier tmpcub = piecewise_Bezier[i];
		for (double t = 0; t < 1; t += step)
		{
			mypoint2f apoint(0, 0);
			double a1 = pow((1 - t), 3);
			double a2 = pow((1 - t), 2) * 3 * t;
			double a3 = 3 * t*t*(1 - t);
			double a4 = t*t*t;
			apoint.x = a1*tmpcub.p0.x + a2*tmpcub.p1.x + a3*tmpcub.p2.x + a4*tmpcub.p3.x;
			apoint.y = a1*tmpcub.p0.y + a2*tmpcub.p1.y + a3*tmpcub.p2.y + a4*tmpcub.p3.y;
			cubplist.push_back(apoint);
		}
		cubplist.push_back(tmpcub.p3);
		new_cubplist.push_back(cubplist);
	}
	//先写灰色的plist_diffcurve
	for (unsigned int i = 0; i < plist_diffcurve.size(); ++i)
	{
		TiXmlElement* apath = new TiXmlElement("path");
		string style_str = "fill:none;stroke-width:1;stroke:#";
		int cnum = 0;
		style_str.append(color_vec[cnum]);  // piecewise cubic curve： 底线是灰色的787878， newcubiccurve是red
		apath->SetAttribute("style", style_str.c_str());
		//apath->SetAttribute("d", "220,100 300,210 170,250");
		//polyline 中点的转换
		string path_str;       //d="M 1062.5,587.134 C 1025.09,592.19 978.75,614.75"
		path_str.append("M");
		path_str.append(" ");
		path_str.append(to_string(plist_diffcurve[i].front().x));//不用TransfloatToStr()，直接用to_string()也行
		path_str.append(",");
		path_str.append(to_string(svg_height - plist_diffcurve[i].front().y));
		path_str.append("L");
		for (unsigned int j = 0; j < plist_diffcurve[i].size(); ++j)
		{
			path_str.append(" ");
			path_str.append(to_string(plist_diffcurve[i].at(j).x));
			path_str.append(",");
			path_str.append(to_string(svg_height-plist_diffcurve[i].at(j).y));
		}
		apath->SetAttribute("d", path_str.c_str());
		rootele->LinkEndChild(apath);
	}
	//在写 红色的新new_cubplist
	for (unsigned int i = 0; i < new_cubplist.size(); ++i)
	{
		TiXmlElement* apath = new TiXmlElement("path");
		string style_str = "fill:none;stroke-width:1;stroke:#";
		int cnum = 1;
		style_str.append(color_vec[cnum]);  // piecewise cubic curve： 底线是灰色的787878， newcubiccurve是red
		apath->SetAttribute("style", style_str.c_str());
		//apath->SetAttribute("d", "220,100 300,210 170,250");
		//polyline 中点的转换
		string path_str;       //d="M 1062.5,587.134 C 1025.09,592.19 978.75,614.75"
		path_str.append("M");
		path_str.append(" ");
		path_str.append(to_string(new_cubplist[i].front().x));//不用TransfloatToStr()，直接用to_string()也行
		path_str.append(",");
		path_str.append(to_string(svg_height-new_cubplist[i].front().y));
		path_str.append("L");
		for (unsigned int j = 0; j < new_cubplist[i].size(); ++j)
		{
			path_str.append(" ");
			path_str.append(to_string(new_cubplist[i].at(j).x));
			path_str.append(",");
			path_str.append(to_string(svg_height-new_cubplist[i].at(j).y));
		}
		apath->SetAttribute("d", path_str.c_str());
		rootele->LinkEndChild(apath);
	}
	string savefile;
	savefile.append("E:/svgfile/piecewise/");
	savefile.append("house");
	savefile.append("_output.svg");
	writeDoc->SaveFile(savefile.c_str());
	delete writeDoc;
}


//PCA降维算法
void MySvgFile::TestPCA()
{
	int density = 1;
	vector<mypoint2f> total_plist;//记录L条线中所有的N个点,一定会有重复的点
	for (unsigned int i = 0; i < cubicBezier_vec_origin.size(); ++i)
	{
		cubicBezier cbcurve = cubicBezier_vec_origin[i];
		double p1[2] = { cbcurve.p0.x, cbcurve.p0.y };
		double p2[2] = { cbcurve.p1.x, cbcurve.p1.y };
		double p3[2] = { cbcurve.p2.x, cbcurve.p2.y };
		double p4[2] = { cbcurve.p3.x, cbcurve.p3.y };
		double totallong = 0;
		totallong = LineSegmentLength(cbcurve.p0, cbcurve.p1) + LineSegmentLength(cbcurve.p1, cbcurve.p2) + LineSegmentLength(cbcurve.p2, cbcurve.p3);
		
		PolyLine apl;//vector<mypoint2f> seg_curve1;
		double a1, a2, a3, a4;
		mypoint2f apoint;
		double step = 0;
		int diviNum = 0;
		if (density == 0)
		{
			step = 0.05;		
		}
		else if (density == 1)
		{
			diviNum =totallong / 0.8;
			step = 1 / (double)diviNum;
			////三次公式 B(t)=p0*(1-t)^3 + p1*3t(1-t)^2 + p2*3t^2(1-t) + p3*t^3
			//double t=0;
			//for (t = 0; t < 1; t += step)
			//{
			//	a1 = pow((1 - t), 3);
			//	a2 = pow((1 - t), 2) * 3 * t;
			//	a3 = 3 * t*t*(1 - t);
			//	a4 = t*t*t;
			//	apoint.x = a1*p1[0] + a2*p2[0] + a3*p3[0] + a4*p4[0];
			//	apoint.y = a1*p1[1] + a2*p2[1] + a3*p3[1] + a4*p4[1];
			//	apl.pointlist.push_back(apoint);
			//		mycircle acir;
			//		acir.cx = apoint.x;
			//		acir.cy = apoint.y;
			//		acir.radius = 0.2;
			//		circle_vec.push_back(acir);
			//}
			//t = 1;
			//a1 = pow((1 - t), 3);
			//a2 = pow((1 - t), 2) * 3 * t;
			//a3 = 3 * pow(t, 2)*(1 - t);
			//a4 = pow(t, 3);
			//apoint.x = a1*p1[0] + a2*p2[0] + a3*p3[0] + a4*p4[0];
			//apoint.y = a1*p1[1] + a2*p2[1] + a3*p3[1] + a4*p4[1];
			//if (fabs(apoint.x - apl.pointlist.back().x)>eps * 10 || fabs(apoint.y - apl.pointlist.back().y)>eps * 10)
			//	apl.pointlist.push_back(apoint);
			//else
			//{
			//	apl.pointlist.pop_back();
			//	apl.pointlist.push_back(apoint);
			//}
			//total_plist.insert(total_plist.end(), apl.pointlist.begin(), apl.pointlist.end());
			//apl.pointlist.clear();
			//apl.pointlist.swap(vector<mypoint2f>());
		}
		//三次公式 B(t)=p0*(1-t)^3 + p1*3t(1-t)^2 + p2*3t^2(1-t) + p3*t^3
		double t = 0;
		for (t = 0; t < 1; t += step)
		{
			a1 = pow((1 - t), 3);
			a2 = pow((1 - t), 2) * 3 * t;
			a3 = 3 * t*t*(1 - t);
			a4 = t*t*t;
			apoint.x = a1*p1[0] + a2*p2[0] + a3*p3[0] + a4*p4[0];
			apoint.y = a1*p1[1] + a2*p2[1] + a3*p3[1] + a4*p4[1];
			apl.pointlist.push_back(apoint);
			mycircle acir;
			acir.cx = apoint.x;
			acir.cy = apoint.y;
			acir.radius = 0.2;
			circle_vec.push_back(acir);
		}
		if (fabs(p4[0] - apl.pointlist.back().x)>eps * 10 || fabs(p4[1] - apl.pointlist.back().y)>eps * 10)
			apl.pointlist.push_back(mypoint2f(p4[0],p4[1]));
		else
		{
			apl.pointlist.pop_back();
			apl.pointlist.push_back(mypoint2f(p4[0], p4[1]));
		}
		total_plist.insert(total_plist.end(), apl.pointlist.begin(), apl.pointlist.end());
		apl.pointlist.clear();
		apl.pointlist.swap(vector<mypoint2f>());
	}
	for (unsigned int i = 0; i < quadraticBezier_vec_origin.size(); ++i)
	{
		cout << "not complish in quadraticBezier_vec" << endl;
	}
	for (unsigned int i = 0; i < polyline_vec_origin.size(); ++i)
	{
		PolyLine pl = polyline_vec_origin[i];
		double line_long=LineSegmentLength(pl.pointlist.front(),pl.pointlist.back());
		double diviNum = line_long / 0.08;
		if (diviNum <= 1)
			diviNum = 2;
		double step = 1 / diviNum;
		mypoint2f direction_vec = pl.pointlist.back() - pl.pointlist.front();
		vector<mypoint2f> new_plist;
		for (int i = 0; i <= diviNum; ++i)
		{
			mypoint2f add_p(step*i*direction_vec.x, step*i*direction_vec.y);
			mypoint2f accum_p = pl.pointlist.front() + add_p;
			new_plist.push_back(accum_p);		

			mycircle acir;
			acir.cx = accum_p.x;
			acir.cy = accum_p.y;
			acir.radius = 0.2;
			circle_vec.push_back(acir);
		}
		total_plist.insert(total_plist.end(), new_plist.begin(), new_plist.end());
	}
	//删除重复的数据点
	vector<mypoint2f>::iterator iter_p = total_plist.begin();
	while (iter_p != total_plist.end())
	{
		vector<mypoint2f>::iterator iter_p_inner = iter_p+1;
		while (iter_p_inner != total_plist.end())
		{
			if (LineSegmentLength(*iter_p, *iter_p_inner) < eps)
			{
				iter_p_inner = total_plist.erase(iter_p_inner);
			}
			else
				iter_p_inner++;
		}
		iter_p++;
	}
	//PCA 降维  直角坐标系
	//先计算协方差矩阵A ,初始化 数据点矩阵MatrixXd B（2，n） 数据按列排列！！
	int plist_size = total_plist.size();
	//vec_x，vec_y，B 都是直角坐标系下的，
	vector<double> vec_x;
	vector<double> vec_y;
	MatrixXd B = MatrixXd::Zero(2, plist_size);	
	for (int i = 0; i < plist_size; ++i)
	{
		vec_x.push_back(total_plist[i].x);
		vec_y.push_back(total_plist[i].y);
		B(0, i) = total_plist[i].x;
		B(1, i) = total_plist[i].y;
	}
	double aver_x = GetAverange(vec_x);
	double aver_y = GetAverange(vec_y);
	double var_cxx = GetVariance(vec_x, aver_x, vec_x, aver_x);
	double var_cyy = GetVariance(vec_y, aver_y, vec_y, aver_y);
	double var_cxy = GetVariance(vec_x, aver_x, vec_y, aver_y);
	Matrix2d A;
	A << var_cxx, var_cxy, var_cxy, var_cyy;
	cout << "Here is a covariance 2*2 matrix, A:" << endl << A << endl << endl;
	//计算协方差矩阵的特征值 特征向量
	EigenSolver<MatrixXd> es(A);
	cout << "The eigenvalues of A are:" << endl << es.eigenvalues() << endl;//注意特征向量和特征值都是 复数类型！！complex
	cout << "The eigenvectors is:" << endl << es.eigenvectors() << endl;
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
	}
	cout << "max_evalue="<<max_evalue << endl;
	cout << "max_v=" << max_vec << endl;
	//projection_m是所有数据点 投影到主特征向量上的结果（投影点），二维降到一维
	VectorXd projection_m = max_vec.transpose()*B;
	//把vectorXd 类型转换成 vector<double>, 并sort()从小到大排序，得到有序点
	vector<pair<int,double>> proMatrix_order; //int对应元数据点在total_plist中的下标，double 是projection_m中的值
	for (int i = 0; i < plist_size; ++i)
	{
		proMatrix_order.push_back(make_pair(i, projection_m[i]));
	}
	sort(proMatrix_order.begin(), proMatrix_order.end(), mycompare);

	/*方法二：固定了前后端点，最小二乘法，参考论文 vectorization of hand-drawn image using piecewise cubic Bezier curves fitting!!！！*/
	mypoint2f cb_p0 = total_plist[proMatrix_order[0].first];
	mypoint2f cb_p3 = total_plist[proMatrix_order.back().first];
	//初始化有序坐标点vec_x_ordered , vec_y_ordered
	VectorXd vec_x_ordered = VectorXd::Zero(plist_size, 1);//即x坐标值 列向量
	VectorXd vec_y_ordered = VectorXd::Zero(plist_size, 1);
	for (int i = 0; i < plist_size; ++i)
	{
		vec_x_ordered[i] = total_plist[proMatrix_order[i].first].x;
		vec_y_ordered[i] = total_plist[proMatrix_order[i].first].y;
		//ordered_plist.push_back(mypoint2f(vec_x_ordered[i], vec_y_ordered[i]));
	}
	//初始化参数矩阵T，弦长参数化？proMatrix_order
	double startp = proMatrix_order.front().second;
	double endp = proMatrix_order.back().second;
	double arclen =abs(endp - startp);
	vector<double> para_t;
	for (int i = 0; i < plist_size; ++i)
	{
		para_t.push_back((proMatrix_order[i].second - startp) / arclen);
	}
	//需要用到参数t(para_t) ,p0p3(cb_p0,cb_p3), 有序的数据点(vec_x/y_ordered),计算A_1, A_2, A_12,C_1,C_2,再计算p1p2
	double A_1 = 0; double A_2 = 0;
	double A_12 = 0;
	//double C_1 = 0; double c_2 = 0;
	Vector2d tmp_c1(0,0);
	Vector2d tmp_c2(0,0);
	Vector2d cb_p0_vec(cb_p0.x, cb_p0.y);
	Vector2d cb_p3_vec(cb_p3.x, cb_p3.y);
	for (int i = 1; i < plist_size-1; ++i)
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
	 cubicBezier newcub;
	 newcub.p0 = cb_p0; //cb_p0; //
	 newcub.p1 = mypoint2f(final_c1[0], final_c1[1]);
	 newcub.p2 = mypoint2f(final_c2[0], final_c2[1]);
	 newcub.p3 = cb_p3;  //cb_p3;  //
	 cubicBezier_vec_origin.push_back(newcub);
	/*方法一：（没有固定前后端点）最小二乘法 拟合cubic Bezier curve ，变量的命名参考自http://jimherold.com/2012/04/20/least-squares-bezier-fit/*/
	////初始化有序坐标点vec_x_ordered , vec_y_ordered
	//// ordered_plist;  //在.h文件里
	//VectorXd vec_x_ordered = VectorXd::Zero(plist_size, 1);//即x坐标值 列向量
	//VectorXd vec_y_ordered = VectorXd::Zero(plist_size, 1);
	//for (int i = 0; i < plist_size; ++i)
	//{
	//	vec_x_ordered[i] = total_plist[proMatrix_order[i].first].x;
	//	vec_y_ordered[i] = total_plist[proMatrix_order[i].first].y;
	//	ordered_plist.push_back(mypoint2f(vec_x_ordered[i], vec_y_ordered[i]));
	//}
	////初始化矩阵M, 赋值符号<<是一行一行的赋值
	//Matrix4d cub_m = Matrix4d::Zero(4, 4);
	//cub_m << -1, 3, -3, 1, 3, -6, 3, 0, -3, 3, 0, 0, 1, 0, 0, 0;
	//Matrix4d cub_m_inverse = cub_m.inverse();
	////初始化参数矩阵T，弦长参数化？proMatrix_order
	//double startp = proMatrix_order.front().second;
	//double endp = proMatrix_order.back().second;
	//double arclen = endp-startp;
	//vector<double> para_t;
	//for (int i = 0; i < plist_size; ++i)
	//{
	//	para_t.push_back((proMatrix_order[i].second - startp) / arclen);
	//}
	//MatrixXd T = MatrixXd::Zero(plist_size, 4);
	//for (int i = 0; i < plist_size; ++i)
	//{
	//	double tmp_t = para_t[i];
	//	T(i, 0) = pow(tmp_t, 3);
	//	T(i, 1) = pow(tmp_t, 2);
	//	T(i, 2) = tmp_t;
	//	T(i, 3) = 1;
	//}
	//MatrixXd T_tran = T.transpose();
	//MatrixXd TT_inverse = (T_tran*T).inverse();
	//cout << cub_m_inverse.rows() << "," << cub_m_inverse.cols() << endl;
	//cout << TT_inverse.rows() << "," << TT_inverse.cols() << endl;
	//cout << T_tran.rows() << "," << T_tran.cols() << endl;
	//cout << vec_y_ordered.rows() << "," << vec_y_ordered.cols() << endl;
	////VectorXd cx
	//VectorXd cx = cub_m_inverse*TT_inverse*T_tran*vec_x_ordered;
	//VectorXd cy = cub_m_inverse*TT_inverse*T_tran*vec_y_ordered;
	// cout << cx<<endl;
	// cout << cy << endl;
	// cubicBezier newcub;
	// mypoint2f cb_p0 = total_plist[proMatrix_order[0].first];
	// mypoint2f cb_p3 = total_plist[proMatrix_order.back().first];
	// newcub.p0 = mypoint2f(cx[0], cy[0]);// cb_p0; // 
	// newcub.p1 = mypoint2f(cx[1], cy[1]);
	// newcub.p2 = mypoint2f(cx[2], cy[2]);
	// newcub.p3 =  mypoint2f(cx[3], cy[3]);   //cb_p3; //
	// cubicBezier_vec_origin.push_back(newcub);


	////PCA 降维  极坐标系 p theta
	////先计算协方差矩阵A ,初始化 数据点矩阵MatrixXd B_polar（2，n） 数据按列排列！！
	//int plist_size = total_plist.size();
	////以下是极坐标系下的变量
	//vector<double> vec_p_polar;
	//vector<double> vec_theta_polar;
	//MatrixXd B_polar = MatrixXd::Zero(2, plist_size);
	////输出极坐标，在matlab里看一下
	//ofstream outfile1("E:/polar_p.txt");
	//ofstream outfile2("E:/polar_theta.txt");
	//outfile1.is_open();
	//outfile2.is_open();
	//for (int i = 0; i < plist_size; ++i)
	//{
	//	double len=LineSegmentLength(total_plist[i], mypoint2f(0, 0));
	//	vec_p_polar.push_back(len);
	//	double tan_theta = atan2(total_plist[i].y, total_plist[i].x);//tan(total_plist[i].y / total_plist[i].x); 
	//	if (tan_theta < 0)
	//		tan_theta = 3.1415926 * 2 - abs(tan_theta);
	//	vec_theta_polar.push_back(tan_theta);
	//	B_polar(0, i) = len;
	//	B_polar(1, i) = tan_theta;
	//	outfile1 << len;
	//	outfile1 << "\n";
	//	outfile2 << tan_theta;
	//	outfile2 << "\n";
	//}
	//outfile1.close();
	//outfile2.close();
	//double aver_p = GetAverange(vec_p_polar);
	//double aver_theta = GetAverange(vec_theta_polar);
	//double var_cxx = GetVariance(vec_p_polar, aver_p, vec_p_polar, aver_p);
	//double var_cyy = GetVariance(vec_theta_polar, aver_theta, vec_theta_polar, aver_theta);
	//double var_cxy = GetVariance(vec_p_polar, aver_p, vec_theta_polar, aver_theta);
	//Matrix2d A;
	//A << var_cxx, var_cxy, var_cxy, var_cyy;
	//cout << "Here is a covariance 2*2 matrix, A:" << endl << A << endl << endl;
	////计算协方差矩阵的特征值 特征向量
	//EigenSolver<MatrixXd> es(A);
	//cout << "The eigenvalues of A are:" << endl << es.eigenvalues() << endl;//注意特征向量和特征值都是 复数类型！！complex
	//cout << "The eigenvectors is:" << endl << es.eigenvectors() << endl;
	////选择较大的那个特征值和其对应的特征向量，作为PCA的主方向
	//double max_evalue = 0;
	//VectorXd max_vec;
	//if (es.eigenvalues()[0].real()> es.eigenvalues()[1].real())
	//{
	//	max_evalue = es.eigenvalues()[0].real();
	//	max_vec = es.eigenvectors().col(0).real();
	//}
	//else
	//{
	//	max_evalue = es.eigenvalues()[1].real();
	//	max_vec = es.eigenvectors().col(1).real();
	//}
	//cout << "max_evalue="<<max_evalue << endl;
	//cout << "max_v=" << max_vec << endl;
	////projection_m是所有数据点 投影到主特征向量上的结果（投影点），二维降到一维
	//VectorXd projection_m = max_vec.transpose()*B_polar;
	////把vectorXd 类型转换成 vector<double>, 并sort()从小到大排序，得到有序点
	//vector<pair<int,double>> proMatrix_order; //int对应元数据点在total_plist中的下标，double 是projection_m中的值
	//for (int i = 0; i < plist_size; ++i)
	//{
	//	proMatrix_order.push_back(make_pair(i, projection_m[i]));
	//}
	//sort(proMatrix_order.begin(), proMatrix_order.end(), mycompare);
	////最小二乘法 拟合cubic Bezier curve ，变量的命名参考自http://jimherold.com/2012/04/20/least-squares-bezier-fit/
	////初始化有序坐标点vec_x_ordered , vec_y_ordered
	//ordered_plist;  //在.h文件里
	//VectorXd vec_x_ordered = VectorXd::Zero(plist_size, 1);//即 x坐标值 列向量
	//VectorXd vec_y_ordered = VectorXd::Zero(plist_size, 1);//y值，列向量
	//for (int i = 0; i < plist_size; ++i)
	//{
	//	vec_x_ordered[i] = total_plist[proMatrix_order[i].first].x;
	//	vec_y_ordered[i] = total_plist[proMatrix_order[i].first].y;
	//	ordered_plist.push_back(mypoint2f(vec_x_ordered[i], vec_y_ordered[i]));
	//}
	////初始化矩阵M, 赋值符号<<是一行一行的赋值
	//Matrix4d cub_m = Matrix4d::Zero(4, 4);
	//cub_m << -1, 3, -3, 1, 3, -6, 3, 0, -3, 3, 0, 0, 1, 0, 0, 0;
	//Matrix4d cub_m_inverse = cub_m.inverse();
	////初始化参数矩阵T，弦长参数化？proMatrix_order
	//double startp = proMatrix_order.front().second;
	//double endp = proMatrix_order.back().second;
	//double arclen = endp-startp;
	//vector<double> para_t;
	//for (int i = 0; i < plist_size; ++i)
	//{
	//	para_t.push_back((proMatrix_order[i].second - startp) / arclen);
	//}
	//MatrixXd T = MatrixXd::Zero(plist_size, 4);
	//for (int i = 0; i < plist_size; ++i)
	//{
	//	double tmp_t = para_t[i];
	//	T(i, 0) = pow(tmp_t, 3);
	//	T(i, 1) = pow(tmp_t, 2);
	//	T(i, 2) = tmp_t;
	//	T(i, 3) = 1;
	//}
	//MatrixXd T_tran = T.transpose();
	//MatrixXd TT_inverse = (T_tran*T).inverse();
	//cout << cub_m_inverse.rows() << "," << cub_m_inverse.cols() << endl;
	//cout << TT_inverse.rows() << "," << TT_inverse.cols() << endl;
	//cout << T_tran.rows() << "," << T_tran.cols() << endl;
	//cout << vec_y_ordered.rows() << "," << vec_y_ordered.cols() << endl;
	////VectorXd cx
	//VectorXd cx = cub_m_inverse*TT_inverse*T_tran*vec_x_ordered;
	//VectorXd cy = cub_m_inverse*TT_inverse*T_tran*vec_y_ordered;
	// cubicBezier newcub;
	// newcub.p0 = mypoint2f(cx[0], cy[0]);//为了跟之前的cubic bezier保持一致，y坐标在翻转一次
	// newcub.p1 = mypoint2f(cx[1], cy[1]);
	// newcub.p2 = mypoint2f(cx[2], cy[2]);
	// newcub.p3 = mypoint2f(cx[3], cy[3]);
	// cubicBezier_vec_origin.push_back(newcub);

	/*例子
	MatrixXd A = MatrixXd::Random(6,6);
	cout << "Here is a random 6x6 matrix, A:" << endl << A << endl << endl;
	EigenSolver<MatrixXd> es(A);
	cout << "The eigenvalues of A are:" << endl << es.eigenvalues() << endl;
	cout << "The matrix of eigenvectors, V, is:" << endl << es.eigenvectors() << endl << endl;
	complex<double> lambda = es.eigenvalues()[0];
	cout << "Consider the first eigenvalue, lambda = " << lambda << endl;
	VectorXcd v = es.eigenvectors().col(0);
	cout << "If v is the corresponding eigenvector, then lambda * v = " << endl << lambda * v << endl;
	cout << "... and A * v = " << endl << A.cast<complex<double> >() * v << endl << endl;
	MatrixXcd D = es.eigenvalues().asDiagonal();
	MatrixXcd V = es.eigenvectors();
	cout << "Finally, V * D * V^(-1) = " << endl << V * D * V.inverse() << endl;
	*/
	/*
	Matrix2d am;
	am << 2, 1, 1, 3;
	EigenSolver<Matrix2d> es(am);
	Matrix2d D = es.pseudoEigenvalueMatrix();
	Matrix2d V = es.pseudoEigenvectors();
	cout << D << endl;
	cout << V << endl;
	*/
	/*
	PCA的算法步骤： 
	设有m条n维数据。 
	1）将原始数据 按列 组成n行m列矩阵X 
	2）将X的每一行（代表一个属性字段）进行零均值化，即减去这一行的均值 
	3）求出协方差矩阵C=XXTC 
	4）求出协方差矩阵的特征值及对应的特征向量 
	5）将特征向量按对应特征值大小从上到下 按行 排列成矩阵，取前k行组成矩阵P 
	6）Y=PX  Y=PX即为降维到k维后的数据
	*/
}

//未完，有错误，怎么排序？
void MySvgFile::LaplacianEigenmap()
{
	int density = 0;
	vector<mypoint2f> total_plist;//记录L条线中所有的的N个点,一定会有重复的点
	for (unsigned int i = 0; i < cubicBezier_vec_origin.size(); ++i)
	{
		cubicBezier cbcurve = cubicBezier_vec_origin[i];
		double p1[2] = { cbcurve.p0.x, cbcurve.p0.y };
		double p2[2] = { cbcurve.p1.x, cbcurve.p1.y };
		double p3[2] = { cbcurve.p2.x, cbcurve.p2.y };
		double p4[2] = { cbcurve.p3.x, cbcurve.p3.y };
		double totallong = 0;
		totallong = LineSegmentLength(cbcurve.p0, cbcurve.p1) + LineSegmentLength(cbcurve.p1, cbcurve.p2) + LineSegmentLength(cbcurve.p2, cbcurve.p3);

		PolyLine apl;//vector<mypoint2f> seg_curve1;
		double a1, a2, a3, a4;
		mypoint2f apoint;
		if (density == 0)
		{
			double step = 0.02;
			//三次公式 B(t)=p0*(1-t)^3 + p1*3t(1-t)^2 + p2*3t^2(1-t) + p3*t^3
			double t = 0;
			for (t = 0; t < 1; t += step)
			{
				a1 = pow((1 - t), 3);
				a2 = pow((1 - t), 2) * 3 * t;
				a3 = 3 * t*t*(1 - t);
				a4 = t*t*t;
				apoint.x = a1*p1[0] + a2*p2[0] + a3*p3[0] + a4*p4[0];
				apoint.y = a1*p1[1] + a2*p2[1] + a3*p3[1] + a4*p4[1];
				apl.pointlist.push_back(apoint);
				mycircle acir;
				acir.cx = apoint.x;
				acir.cy = apoint.y;
				acir.radius = 0.2;
				circle_vec.push_back(acir);
			}
			t = 1;
			a1 = pow((1 - t), 3);
			a2 = pow((1 - t), 2) * 3 * t;
			a3 = 3 * pow(t, 2)*(1 - t);
			a4 = pow(t, 3);
			apoint.x = a1*p1[0] + a2*p2[0] + a3*p3[0] + a4*p4[0];
			apoint.y = a1*p1[1] + a2*p2[1] + a3*p3[1] + a4*p4[1];
			if (fabs(apoint.x - apl.pointlist.back().x)>eps * 10 || fabs(apoint.y - apl.pointlist.back().y)>eps * 10)
				apl.pointlist.push_back(apoint);
			else
			{
				apl.pointlist.pop_back();
				apl.pointlist.push_back(apoint);
			}
			total_plist.insert(total_plist.end(), apl.pointlist.begin(), apl.pointlist.end());
			apl.pointlist.clear();
			apl.pointlist.swap(vector<mypoint2f>());
		}
		else if (density == 1)
		{
			int diviNum = totallong / 0.8;
			double step = 1 / (double)diviNum;
			//三次公式 B(t)=p0*(1-t)^3 + p1*3t(1-t)^2 + p2*3t^2(1-t) + p3*t^3
			double t = 0;
			for (t = 0; t < 1; t += step)
			{
				a1 = pow((1 - t), 3);
				a2 = pow((1 - t), 2) * 3 * t;
				a3 = 3 * t*t*(1 - t);
				a4 = t*t*t;
				apoint.x = a1*p1[0] + a2*p2[0] + a3*p3[0] + a4*p4[0];
				apoint.y = a1*p1[1] + a2*p2[1] + a3*p3[1] + a4*p4[1];
				apl.pointlist.push_back(apoint);
				mycircle acir;
				acir.cx = apoint.x;
				acir.cy = apoint.y;
				acir.radius = 0.2;
				circle_vec.push_back(acir);
			}
			t = 1;
			a1 = pow((1 - t), 3);
			a2 = pow((1 - t), 2) * 3 * t;
			a3 = 3 * pow(t, 2)*(1 - t);
			a4 = pow(t, 3);
			apoint.x = a1*p1[0] + a2*p2[0] + a3*p3[0] + a4*p4[0];
			apoint.y = a1*p1[1] + a2*p2[1] + a3*p3[1] + a4*p4[1];
			if (fabs(apoint.x - apl.pointlist.back().x)>eps * 10 || fabs(apoint.y - apl.pointlist.back().y)>eps * 10)
				apl.pointlist.push_back(apoint);
			else
			{
				apl.pointlist.pop_back();
				apl.pointlist.push_back(apoint);
			}
			total_plist.insert(total_plist.end(), apl.pointlist.begin(), apl.pointlist.end());
			apl.pointlist.clear();
			apl.pointlist.swap(vector<mypoint2f>());
		}
	}

	for (unsigned int i = 0; i < quadraticBezier_vec_origin.size(); ++i)
	{
		cout << "not complish in quadraticBezier_vec" << endl;
	}
	for (unsigned int i = 0; i < polyline_vec_origin.size(); ++i)
	{
		PolyLine pl = polyline_vec_origin[i];
		total_plist.insert(total_plist.end(), pl.pointlist.begin(), pl.pointlist.end());
	}
	//删除重复的数据点
	vector<mypoint2f>::iterator iter_p = total_plist.begin();
	while (iter_p != total_plist.end())
	{
		vector<mypoint2f>::iterator iter_p_inner = iter_p + 1;
		while (iter_p_inner != total_plist.end())
		{
			if (LineSegmentLength(*iter_p, *iter_p_inner) < eps)
			{
				iter_p_inner = total_plist.erase(iter_p_inner);
			}
			else
				iter_p_inner++;
		}
		iter_p++;
	}

	//laplacian Eigenmap降维  
	////先计算协方差矩阵A ,初始化 数据点矩阵MatrixXd B_polar（2，n） 数据按列排列！！
	//for (unsigned int i = 0; i < total_plist.size(); ++i)  //先把y翻转过来
	//{
	//	total_plist[i].y = 1052.0 - total_plist[i].y;
	//}
	int plist_size = total_plist.size();
	//计算关联矩阵W + 高斯核函数。  即 W(ij)=exp(-dij^2/(alpha*sigma)), 根据论文，W(ij)>0.75则保留，否则为0
	MatrixXd W = MatrixXd::Zero(plist_size, plist_size);
	vector<double> allpairs_len;//记录每对数据点对的距离，为了下一步计算标准差用
	for (int i = 0; i < plist_size; ++i)
	{
		for (int j = i; j < plist_size; ++j)
		{
			double len_square = pow(total_plist[i].x - total_plist[j].x, 2) + pow(total_plist[i].y - total_plist[j].y, 2);
			W(i, j) = len_square;         //记录每对数据点对的距离的平方dij^2 存入矩阵W
			W(j, i) = len_square;
			allpairs_len.push_back(sqrt(len_square));  
		}
	}
	double aver = GetAverange(allpairs_len);//均值
	double std_var = sqrt(GetVariance_xx(allpairs_len, aver));//标准差
	double deno = 0.5*std_var;  //分母=系数0.5*标准差
	for (int i = 0; i < plist_size; ++i)
	{
		for (int j = i; j < plist_size; ++j)
		{
			double tmp = exp(-W(i, j) / deno);
			if (tmp > 0.75)
			{
				W(i, j) = tmp;
				W(j, i) = tmp;
			}
			else
			{
				W(i, j) = 0;
				W(j, i) = 0;
			}
		}
	}
	//计算对角矩阵D，D(ii)=W矩阵第i行元素之和。
	MatrixXd D = MatrixXd::Zero(plist_size, plist_size);
	for (int i = 0; i < plist_size; ++i)
	{
		double tmp_sum = 0;
		for (int j = 0; j < plist_size; ++j)
		{
			tmp_sum = tmp_sum + W(i, j);
		}
		D(i, i) = tmp_sum;
	}
	//矩阵L=D-W，计算广义矩阵特征值分解
	MatrixXd L = MatrixXd::Zero(plist_size, plist_size);
	L = D - W;
	//cout << L << endl;
	MatrixXd dl_matrix = D.inverse()*L;
	EigenSolver<MatrixXd> es(dl_matrix);
	cout << "The eigenvalues of A are:" << endl;// es.eigenvalues()是从小到大排序的，es.eigenvectors().col(0)一列是一个特征向量
	cout<< es.eigenvalues()[0] << endl;
	cout << es.eigenvalues()[1] << endl;
	cout << es.eigenvalues()[2] << endl;
	cout << es.eigenvalues()[plist_size - 1] << endl;
	cout << es.eigenvalues()[plist_size - 2] << endl;
	cout << es.eigenvalues()[plist_size - 3] << endl;
	cout << "The matrix of eigenvectors, V, is:" << endl;
	//cout<< es.eigenvectors()[0].real() << endl;
	VectorXcd v1 = es.eigenvectors().col(2);//用前两列作为 新spectral domain中的坐标  plist_size - 1  plist_size - 2
	VectorXcd v2 = es.eigenvectors().col(3);

	//把vectorXd 类型转换成 vector<double>, 并sort()从小到大排序，得到有序点
	vector<pair<int, double>> horizontal_order; //int对应元数据点在total_plist中的下标，double 是projection_m中的值
	for (int i = 0; i < plist_size; ++i)
	{
		horizontal_order.push_back(make_pair(i, v1[i].real()));
	}
	sort(horizontal_order.begin(), horizontal_order.end(), mycompare);
	cout << endl;
	//方法一：（没有固定前后端点）最小二乘法 拟合cubic Bezier curve ，变量的命名参考自http://jimherold.com/2012/04/20/least-squares-bezier-fit/
	//初始化有序坐标点vec_x_ordered , vec_y_ordered
	// ordered_plist;  //在.h文件里
	VectorXd vec_x_ordered = VectorXd::Zero(plist_size, 1);//即x坐标值 列向量
	VectorXd vec_y_ordered = VectorXd::Zero(plist_size, 1);
	for (int i = 0; i < plist_size; ++i)
	{
		vec_x_ordered[i] = total_plist[horizontal_order[i].first].x;
		vec_y_ordered[i] = total_plist[horizontal_order[i].first].y;
		//ordered_plist.push_back(mypoint2f(vec_x_ordered[i], vec_y_ordered[i]));
	}
	//初始化矩阵M, 赋值符号<<是一行一行的赋值
	Matrix4d cub_m = Matrix4d::Zero(4, 4);
	cub_m << -1, 3, -3, 1, 3, -6, 3, 0, -3, 3, 0, 0, 1, 0, 0, 0;
	Matrix4d cub_m_inverse = cub_m.inverse();
	//初始化参数矩阵T，弦长参数化？proMatrix_order
	double arclen = 1.0 / (plist_size-1);
	vector<double> para_t;
	for (int i = 0; i < plist_size; ++i)
	{
		para_t.push_back(arclen*i);
	}
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
	cout << cub_m_inverse.rows() << "," << cub_m_inverse.cols() << endl;
	cout << TT_inverse.rows() << "," << TT_inverse.cols() << endl;
	cout << T_tran.rows() << "," << T_tran.cols() << endl;
	cout << vec_y_ordered.rows() << "," << vec_y_ordered.cols() << endl;
	//VectorXd cx
	VectorXd cx = cub_m_inverse*TT_inverse*T_tran*vec_x_ordered;
	VectorXd cy = cub_m_inverse*TT_inverse*T_tran*vec_y_ordered;
	cout << cx << endl;
	cout << cy << endl;
	cubicBezier newcub;
	mypoint2f cb_p0 = total_plist[horizontal_order[0].first];
	mypoint2f cb_p3 = total_plist[horizontal_order.back().first];
	newcub.p0 = mypoint2f(cx[0], cy[0]);// cb_p0; //
	newcub.p1 = mypoint2f(cx[1], cy[1]);
	newcub.p2 = mypoint2f(cx[2], cy[2]);
	newcub.p3 = mypoint2f(cx[3], cy[3]);   // cb_p3; // 
	cubicBezier_vec_origin.push_back(newcub);

	


	//EigenSolver<MatrixXd> es(L);    //稀疏矩阵 求特征值和特征向量 速度非常非常慢！！！怎么计算？
	//MatrixXd ei_val = es.pseudoEigenvalueMatrix();
	//MatrixXd ei_vec = es.pseudoEigenvectors();
	//cout << ei_val << endl;
	//cout << ei_vec << endl;
	/*例子
	MatrixXd A = MatrixXd::Random(6,6);
	cout << "Here is a random 6x6 matrix, A:" << endl << A << endl << endl;
	EigenSolver<MatrixXd> es(A);
	cout << "The eigenvalues of A are:" << endl << es.eigenvalues() << endl;
	cout << "The matrix of eigenvectors, V, is:" << endl << es.eigenvectors() << endl << endl;
	complex<double> lambda = es.eigenvalues()[0];
	cout << "Consider the first eigenvalue, lambda = " << lambda << endl;
	VectorXcd v = es.eigenvectors().col(0);
	cout << "If v is the corresponding eigenvector, then lambda * v = " << endl << lambda * v << endl;
	cout << "... and A * v = " << endl << A.cast<complex<double> >() * v << endl << endl;
	MatrixXcd D = es.eigenvalues().asDiagonal();
	MatrixXcd V = es.eigenvectors();
	cout << "Finally, V * D * V^(-1) = " << endl << V * D * V.inverse() << endl;
	*/
	/*
	Matrix2d am;
	am << 2, 1, 1, 3;
	EigenSolver<Matrix2d> es(am);
	Matrix2d D = es.pseudoEigenvalueMatrix();
	Matrix2d V = es.pseudoEigenvectors();
	cout << D << endl;
	cout << V << endl;
	*/
	
}
//平均曲线
void MySvgFile::TestAverageCurve()
{
	int density = 1;

	vector<vector<mypoint2f>> total_plist;//记录L条线中每条线的采样点
	vector<vector<mypoint2f>> total_normal;//normal
	vector<vector<mypoint2f>> total_tangent;
	for (unsigned int i = 0; i < cubicBezier_vec_origin.size(); ++i)
	{
		vector<mypoint2f> normallist;//(y'(t),-x'(t))
		vector<mypoint2f> tangelist;//(x'(t),y'(t))
		PolyLine apl;//vector<mypoint2f> seg_curve1;

		cubicBezier cbcurve = cubicBezier_vec_origin[i];
		double p1[2] = { cbcurve.p0.x, cbcurve.p0.y };
		double p2[2] = { cbcurve.p1.x, cbcurve.p1.y };
		double p3[2] = { cbcurve.p2.x, cbcurve.p2.y };
		double p4[2] = { cbcurve.p3.x, cbcurve.p3.y };
		double totallong = 0;
		totallong = LineSegmentLength(cbcurve.p0, cbcurve.p1) + LineSegmentLength(cbcurve.p1, cbcurve.p2) + LineSegmentLength(cbcurve.p2, cbcurve.p3);
		double a1, a2, a3, a4;
		
		int diviNum = 0;
		double step = 0;
		mypoint2f apoint;
		if (density == 0)
		{
			step = 0.02;
		}
		else if (density == 1)
		{
			diviNum = totallong / 0.8;
			step = 1 / (double)diviNum;
			////三次公式 B(t)=p0*(1-t)^3 + p1*3t(1-t)^2 + p2*3t^2(1-t) + p3*t^3
			//double t = 0;
			//for (t = 0; t < 1; t += step)
			//{
			//	a1 = pow((1 - t), 3);
			//	a2 = pow((1 - t), 2) * 3 * t;
			//	a3 = 3 * t*t*(1 - t);
			//	a4 = t*t*t;
			//	apoint.x = a1*p1[0] + a2*p2[0] + a3*p3[0] + a4*p4[0];
			//	apoint.y = a1*p1[1] + a2*p2[1] + a3*p3[1] + a4*p4[1];
			//	apl.pointlist.push_back(apoint);
			//	mycircle acir;
			//	acir.cx = apoint.x;
			//	acir.cy = apoint.y;
			//	acir.radius = 0.2;
			//	circle_vec.push_back(acir);
			//}
			//t = 1;
			//a1 = pow((1 - t), 3);
			//a2 = pow((1 - t), 2) * 3 * t;
			//a3 = 3 * pow(t, 2)*(1 - t);
			//a4 = pow(t, 3);
			//apoint.x = a1*p1[0] + a2*p2[0] + a3*p3[0] + a4*p4[0];
			//apoint.y = a1*p1[1] + a2*p2[1] + a3*p3[1] + a4*p4[1];
			//if (fabs(apoint.x - apl.pointlist.back().x)>eps * 10 || fabs(apoint.y - apl.pointlist.back().y)>eps * 10)
			//	apl.pointlist.push_back(apoint);
			//else
			//{
			//	apl.pointlist.pop_back();
			//	apl.pointlist.push_back(apoint);
			//}
			//total_plist.push_back(apl.pointlist);
			//apl.pointlist.clear();
			//apl.pointlist.swap(vector<mypoint2f>());
		}
		//三次公式 B(t)=p0*(1-t)^3 + p1*3t(1-t)^2 + p2*3t^2(1-t) + p3*t^3
		double t = 0;
		for (t = 0; t < 1; t += step)
		{
			a1 = pow((1 - t), 3);
			a2 = pow((1 - t), 2) * 3 * t;
			a3 = 3 * t*t*(1 - t);
			a4 = t*t*t;
			apoint.x = a1*p1[0] + a2*p2[0] + a3*p3[0] + a4*p4[0];
			apoint.y = a1*p1[1] + a2*p2[1] + a3*p3[1] + a4*p4[1];
			apl.pointlist.push_back(apoint);//采样点
			double A = pow(t, 3)*(-p1[0] + 3 * p2[0] - 3 * p3[0] + p4[0]) + pow(t, 2)*(3 * p1[0] - 6 * p2[0] + 3 * p3[0]) + t*(-3 * p1[0] + 3 * p2[0]) + p1[0];
			double B = pow(t, 3)*(-p1[1] + 3 * p2[1] - 3 * p3[1] + p4[1]) + pow(t, 2)*(3 * p1[1] - 6 * p2[1] + 3 * p3[1]) + t*(-3 * p1[1] + 3 * p2[1]) + p1[1];
			mypoint2f tmp_p = mypoint2f(A, B);
			double vector_length = LineSegmentLength(tmp_p, mypoint2f(0, 0));
			tmp_p.x = tmp_p.x / vector_length;
			tmp_p.y = tmp_p.y / vector_length;
			tangelist.push_back(tmp_p);//切向量（归一化的）
			mypoint2f tmp_p2 =  mypoint2f(tmp_p.y, -tmp_p.x);
			normallist.push_back(tmp_p2);//法向量
		}
		t = 1;
		a1 = pow((1 - t), 3);
		a2 = pow((1 - t), 2) * 3 * t;
		a3 = 3 * pow(t, 2)*(1 - t);
		a4 = pow(t, 3);
		apoint.x = a1*p1[0] + a2*p2[0] + a3*p3[0] + a4*p4[0];
		apoint.y = a1*p1[1] + a2*p2[1] + a3*p3[1] + a4*p4[1];
		double A = pow(t, 3)*(-p1[0] + 3 * p2[0] - 3 * p3[0] + p4[0]) + pow(t, 2)*(3 * p1[0] - 6 * p2[0] + 3 * p3[0]) + t*(-3 * p1[0] + 3 * p2[0]) + p1[0];
		double B = pow(t, 3)*(-p1[1] + 3 * p2[1] - 3 * p3[1] + p4[1]) + pow(t, 2)*(3 * p1[1] - 6 * p2[1] + 3 * p3[1]) + t*(-3 * p1[1] + 3 * p2[1]) + p1[1];
		mypoint2f tmp_p = mypoint2f(A, B);
		double vector_length = LineSegmentLength(tmp_p, mypoint2f(0, 0));
		tmp_p.x = tmp_p.x / vector_length;
		tmp_p.y = tmp_p.y / vector_length;
		mypoint2f tmp_p2 =mypoint2f(tmp_p.y, -tmp_p.x);
		if (fabs(apoint.x - apl.pointlist.back().x) > eps * 10 || fabs(apoint.y - apl.pointlist.back().y) > eps * 10)
		{
			apl.pointlist.push_back(apoint);
			tangelist.push_back(tmp_p);
			normallist.push_back(tmp_p2);
		}
		else
		{
			apl.pointlist.pop_back();
			apl.pointlist.push_back(apoint);
			tangelist.pop_back();
			tangelist.push_back(tmp_p);
			normallist.pop_back();
			normallist.push_back(tmp_p2);
		}
		total_plist.push_back(apl.pointlist);
		total_tangent.push_back(tangelist);
		total_normal.push_back(normallist);
	}

	for (unsigned int i = 0; i < quadraticBezier_vec_origin.size(); ++i)
	{
		cout << "not complish in quadraticBezier_vec" << endl;
	}
	for (unsigned int i = 0; i < polyline_vec_origin.size(); ++i)
	{
		PolyLine pl = polyline_vec_origin[i];
		total_plist.push_back(pl.pointlist);
	}
	//所有笔画起点的重心
	mypoint2f centroid_p_start(0, 0);
	mypoint2f centroid_p_end(0, 0);
	for (unsigned int i = 0; i < total_plist.size(); ++i)
	{
		centroid_p_start = centroid_p_start + total_plist[i].front();
		centroid_p_end = centroid_p_end + total_plist[i].back();
	}
	centroid_p_end.x = centroid_p_end.x / (double)total_plist.size();
	centroid_p_end.y = centroid_p_end.y / (double)total_plist.size();
	centroid_p_start.x = centroid_p_start.x / (double)total_plist.size();
	centroid_p_start.y = centroid_p_start.y / (double)total_plist.size();
	vector<mypoint2f> final_plist;
	vector<int> start_grp(total_plist.size(), 0);
	vector<int> grp_origi=computeGRP(centroid_p_start, start_grp, total_plist);
	mypoint2f normal_origi(0, 0);
	mypoint2f tmp_p = snap(total_normal, total_plist, centroid_p_start, grp_origi, normal_origi);
	final_plist.push_back(tmp_p);
	mypoint2f prior_point(0, 0);
	double stop_flag = 65535;
	double step_e = 0.01;//step_e越大 越平滑！
	while (stop_flag > 10)
	{
		prior_point = final_plist.back();
		/*vector<int> Ui = computeGRP(prior_point, grp_origi, total_plist);
		mypoint2f tangent_p = snap(total_tangent, total_plist, prior_point, Ui);
		grp_origi = Ui;
		tangent_p.x = tangent_p.x*step_e;
		tangent_p.y = tangent_p.y*step_e;*/
		mypoint2f tangent_p(-normal_origi.y, normal_origi.x);///!!!!!
		tangent_p.x = tangent_p.x*step_e;
		tangent_p.y = tangent_p.y*step_e;
		mypoint2f next_p = prior_point + tangent_p;
		tmp_p = snap(total_normal, total_plist, next_p, grp_origi, normal_origi);
		final_plist.push_back(tmp_p);
		stop_flag = LineSegmentLength(tmp_p, centroid_p_end);
	}
	PolyLine apl;
	apl.pointlist = final_plist;
	polyline_vec.push_back(apl);
}
vector<int> MySvgFile::computeGRP(mypoint2f point, vector<int> startindex, vector<vector<mypoint2f>> line_list)
{
	vector<int> GRP_index;
	for (unsigned int i = 0; i < line_list.size(); ++i)
	{
		mypoint2f begin_p = line_list[i].at(startindex[i]);
		double mini = 65565;
		unsigned int j = 0;
		for (j = startindex[i]; j < line_list[i].size(); ++j)
		{
			begin_p = line_list[i].at(j);
			double tmp = LineSegmentLength(point, begin_p);
			if (tmp < mini)
			{
				mini = tmp;
			}
			else
				break;
		}
		GRP_index.push_back(j - 1);
	}
	return GRP_index;
}
mypoint2f MySvgFile::algorithm_one(vector<vector<mypoint2f>> normal_list, vector<vector<mypoint2f>> line_list, mypoint2f apoint, vector<int>& start_GRP_index, mypoint2f& N_newline)
{
	vector<int> GRP_index = computeGRP(apoint, start_GRP_index, line_list);
	for (unsigned int i = 0; i < start_GRP_index.size(); ++i)
	{
		start_GRP_index[i] = GRP_index[i];
	}
	double a = 0;
	double b = 0;
	double c = 0;
	double n = 0;
	double d = 0;
	for (unsigned int i = 0; i < normal_list.size(); ++i)
	{
		a = a + pow(normal_list[i].at(GRP_index[i]).x, 2);
		b = b + pow(normal_list[i].at(GRP_index[i]).y, 2);
		c = c + (normal_list[i].at(GRP_index[i]).x*normal_list[i].at(GRP_index[i]).y);
	}
	double alpha = atan2(2 * c, (a - b)) / 2;//有两个点 point(x1,y1), 和 point(x2,y2);  float angle = atan2( y2-y1, x2-x1 );
	mypoint2f N(cos(alpha), sin(alpha));
	N_newline = N;
	for (unsigned int i = 0; i < line_list.size(); ++i)
	{
		double ppi_ni = dot(line_list[i].at(GRP_index[i]) - apoint, normal_list[i].at(GRP_index[i]));
		double n_ni = dot(N, normal_list[i].at(GRP_index[i]));
		n = n + ppi_ni*n_ni;
		d = d + pow(n_ni, 2);
	}
	N.x = N.x*(n / d);
	N.y = N.y*(n / d);
	mypoint2f V = apoint + N;
	return V;
}
mypoint2f MySvgFile::snap(vector<vector<mypoint2f>> normal_list, vector<vector<mypoint2f>> line_list, mypoint2f p, vector<int>& grp_vec, mypoint2f& normal_p)
{
	mypoint2f projection_point=p;
	mypoint2f proir_point = p;
	mypoint2f normal_tmp_p(0,0);
	vector<int> prior_start_grp = grp_vec;// (line_list.size(), 0);
	for (int i = 0; i < 3; ++i)
	{
		projection_point = algorithm_one(normal_list, line_list, proir_point, prior_start_grp, normal_tmp_p);//返回投影点，且prior_start_grp被更新了
		proir_point = projection_point;
	}
	for (unsigned int i = 0; i < prior_start_grp.size(); ++i)
	{
		grp_vec[i]=prior_start_grp[i];
	}
	normal_p = normal_tmp_p;//projection_point处的法向量
	return projection_point;
}

//piecewise Besier curve fitting ，
void MySvgFile::piecewiseBezieFitting()
{
	int density = 1;
	vector<mypoint2f> total_plist;//记录L条线中所有的N个点,一定会有重复的点
	for (unsigned int i = 0; i < cubicBezier_vec_origin.size(); ++i)
	{
		cubicBezier cbcurve = cubicBezier_vec_origin[i];
		double p1[2] = { cbcurve.p0.x, cbcurve.p0.y };
		double p2[2] = { cbcurve.p1.x, cbcurve.p1.y };
		double p3[2] = { cbcurve.p2.x, cbcurve.p2.y };
		double p4[2] = { cbcurve.p3.x, cbcurve.p3.y };
		double totallong = 0;
		totallong = LineSegmentLength(cbcurve.p0, cbcurve.p1) + LineSegmentLength(cbcurve.p1, cbcurve.p2) + LineSegmentLength(cbcurve.p2, cbcurve.p3);

		PolyLine apl;//vector<mypoint2f> seg_curve1;
		double a1, a2, a3, a4;
		mypoint2f apoint;
		double step = 0;
		int diviNum = 0;
		if (density == 0)
		{
			step = 0.05;
		}
		else if (density == 1)
		{
			diviNum = totallong / 0.8;
			step = 1 / (double)diviNum;
			////三次公式 B(t)=p0*(1-t)^3 + p1*3t(1-t)^2 + p2*3t^2(1-t) + p3*t^3
			//double t=0;
			//for (t = 0; t < 1; t += step)
			//{
			//	a1 = pow((1 - t), 3);
			//	a2 = pow((1 - t), 2) * 3 * t;
			//	a3 = 3 * t*t*(1 - t);
			//	a4 = t*t*t;
			//	apoint.x = a1*p1[0] + a2*p2[0] + a3*p3[0] + a4*p4[0];
			//	apoint.y = a1*p1[1] + a2*p2[1] + a3*p3[1] + a4*p4[1];
			//	apl.pointlist.push_back(apoint);
			//		mycircle acir;
			//		acir.cx = apoint.x;
			//		acir.cy = apoint.y;
			//		acir.radius = 0.2;
			//		circle_vec.push_back(acir);
			//}
			//t = 1;
			//a1 = pow((1 - t), 3);
			//a2 = pow((1 - t), 2) * 3 * t;
			//a3 = 3 * pow(t, 2)*(1 - t);
			//a4 = pow(t, 3);
			//apoint.x = a1*p1[0] + a2*p2[0] + a3*p3[0] + a4*p4[0];
			//apoint.y = a1*p1[1] + a2*p2[1] + a3*p3[1] + a4*p4[1];
			//if (fabs(apoint.x - apl.pointlist.back().x)>eps * 10 || fabs(apoint.y - apl.pointlist.back().y)>eps * 10)
			//	apl.pointlist.push_back(apoint);
			//else
			//{
			//	apl.pointlist.pop_back();
			//	apl.pointlist.push_back(apoint);
			//}
			//total_plist.insert(total_plist.end(), apl.pointlist.begin(), apl.pointlist.end());
			//apl.pointlist.clear();
			//apl.pointlist.swap(vector<mypoint2f>());
		}
		//三次公式 B(t)=p0*(1-t)^3 + p1*3t(1-t)^2 + p2*3t^2(1-t) + p3*t^3
		double t = 0;
		for (t = 0; t < 1; t += step)
		{
			a1 = pow((1 - t), 3);
			a2 = pow((1 - t), 2) * 3 * t;
			a3 = 3 * t*t*(1 - t);
			a4 = t*t*t;
			apoint.x = a1*p1[0] + a2*p2[0] + a3*p3[0] + a4*p4[0];
			apoint.y = a1*p1[1] + a2*p2[1] + a3*p3[1] + a4*p4[1];
			apl.pointlist.push_back(apoint);
			mycircle acir;
			acir.cx = apoint.x;
			acir.cy = apoint.y;
			acir.radius = 0.2;
			circle_vec.push_back(acir);
		}
		if (fabs(p4[0] - apl.pointlist.back().x)>(eps * 10) || fabs(p4[1] - apl.pointlist.back().y)>(eps * 10))
			apl.pointlist.push_back(mypoint2f(p4[0], p4[1]));
		else
		{
			apl.pointlist.pop_back();
			apl.pointlist.push_back(mypoint2f(p4[0], p4[1]));
		}
		total_plist.insert(total_plist.end(), apl.pointlist.begin(), apl.pointlist.end());
		apl.pointlist.clear();
		apl.pointlist.swap(vector<mypoint2f>());
	}
	for (unsigned int i = 0; i < quadraticBezier_vec_origin.size(); ++i)
	{
		cout << "not complish in quadraticBezier_vec" << endl;
	}
	for (unsigned int i = 0; i < polyline_vec_origin.size(); ++i)
	{
		PolyLine pl = polyline_vec_origin[i];
		double line_long = LineSegmentLength(pl.pointlist.front(), pl.pointlist.back());
		double diviNum = line_long / 0.08;
		if (diviNum <= 1)
			diviNum = 2;
		double step = 1 / diviNum;
		mypoint2f direction_vec = pl.pointlist.back() - pl.pointlist.front();
		vector<mypoint2f> new_plist;
		for (int i = 0; i <= diviNum; ++i)
		{
			mypoint2f add_p(step*i*direction_vec.x, step*i*direction_vec.y);
			mypoint2f accum_p = pl.pointlist.front() + add_p;
			new_plist.push_back(accum_p);

			mycircle acir;
			acir.cx = accum_p.x;
			acir.cy = accum_p.y;
			acir.radius = 0.2;
			circle_vec.push_back(acir);
		}
		total_plist.insert(total_plist.end(), new_plist.begin(), new_plist.end());
	}
	//删除重复的数据点
	vector<mypoint2f>::iterator iter_p = total_plist.begin();
	while (iter_p != total_plist.end())
	{
		vector<mypoint2f>::iterator iter_p_inner = iter_p + 1;
		while (iter_p_inner != total_plist.end())
		{
			if (LineSegmentLength(*iter_p, *iter_p_inner) < eps)
			{
				iter_p_inner = total_plist.erase(iter_p_inner);
			}
			else
				iter_p_inner++;
		}
		iter_p++;
	}
	//PCA 降维  直角坐标系
	//先计算协方差矩阵A ,初始化 数据点矩阵MatrixXd B（2，n） 数据按列排列！！
	int plist_size = total_plist.size();
	//vec_x，vec_y，B 都是直角坐标系下的，
	vector<double> vec_x;
	vector<double> vec_y;
	MatrixXd B = MatrixXd::Zero(2, plist_size);
	for (int i = 0; i < plist_size; ++i)
	{
		vec_x.push_back(total_plist[i].x);
		vec_y.push_back(total_plist[i].y);
		B(0, i) = total_plist[i].x;
		B(1, i) = total_plist[i].y;
	}
	double aver_x = GetAverange(vec_x);
	double aver_y = GetAverange(vec_y);
	double var_cxx = GetVariance(vec_x, aver_x, vec_x, aver_x);
	double var_cyy = GetVariance(vec_y, aver_y, vec_y, aver_y);
	double var_cxy = GetVariance(vec_x, aver_x, vec_y, aver_y);
	Matrix2d A;
	A << var_cxx, var_cxy, var_cxy, var_cyy;
	cout << "Here is a covariance 2*2 matrix, A:" << endl << A << endl << endl;
	//计算协方差矩阵的特征值 特征向量
	EigenSolver<MatrixXd> es(A);
	cout << "The eigenvalues of A are:" << endl << es.eigenvalues() << endl;//注意特征向量和特征值都是 复数类型！！complex
	cout << "The eigenvectors is:" << endl << es.eigenvectors() << endl;
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
	}
	cout << "max_evalue=" << max_evalue << endl;
	cout << "max_v=" << max_vec << endl;
	//projection_m是所有数据点 投影到主特征向量上的结果（投影点），二维降到一维
	VectorXd projection_m = max_vec.transpose()*B;
	//把vectorXd 类型转换成 vector<double>, 并sort()从小到大排序，得到有序点
	vector<pair<int, double>> proMatrix_order; //int对应元数据点在total_plist中的下标，double 是projection_m中的值
	for (int i = 0; i < plist_size; ++i)
	{
		proMatrix_order.push_back(make_pair(i, projection_m[i]));
	}
	sort(proMatrix_order.begin(), proMatrix_order.end(), mycompare);

	/*方法三：piecewise curve fitting  */
	vector<mypoint2f> ptdata;
	for (int i = 0; i < plist_size; ++i)
	{
		ptdata.push_back(total_plist[proMatrix_order[i].first]);
	}
	double	error = 5;		/*  Squared error */
							//vector<cubicBezier> store_cb;
	FitCurves(&ptdata, ptdata.size(), error, &piecewise_Bezier);		/*  Fit the Bezier curves */
	for (size_t i = 0; i < piecewise_Bezier.size(); ++i)
	{
		mycircle cc;
		cc.cx = piecewise_Bezier[i].p0.x;
		cc.cy = piecewise_Bezier[i].p0.y;
		cc.radius = 1;
		circle_vec.push_back(cc);
		//mycircle cc;
		cc.cx = piecewise_Bezier[i].p1.x;
		cc.cy = piecewise_Bezier[i].p1.y;
		cc.radius = 1;
		circle_vec.push_back(cc);
		//mycircle cc;
		cc.cx = piecewise_Bezier[i].p2.x;
		cc.cy = piecewise_Bezier[i].p2.y;
		cc.radius = 1;
		circle_vec.push_back(cc);
		//mycircle cc;
		cc.cx = piecewise_Bezier[i].p3.x;
		cc.cy = piecewise_Bezier[i].p3.y;
		cc.radius = 1;
		circle_vec.push_back(cc);
	}
	cout <<"piecewise cubic curve number:"<<piecewise_Bezier.size()<< endl;
}