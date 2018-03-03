//#include<iostream>
//#include<stdlib.h>
//#include"MySvgFile.h"
#include"MyGraph.h"


string FILENAME = "grandpa";
string FILEPATH = "F:/sketch_input/";
string OUTPUT_PATH = "E:/svgfile/p+l_equal/";
//char *fname = "F:/sketch_input/tt.svg";//spiral  clock  ring   grandpa-hat  tt  closure  girl duck  car  house  shesh  geometry 

//string filepath="E:/true2form/true2form_supplyment/InputSVG/camera.svg";
//char *fname = "E:/true2form/true2form_supplyment/InputSVG/GT_cross.svg";
//char *fname = "E:/true2form/true2form_supplyment/InputSVG/camera.svg";
//char *fname ="E:/fidelity_simplicity_supplyment/figure_1/our.svg";
//char *fname = "E:/fidelity_simplicity_supplyment/figure_2/output.svg";
//char *fname = "E:/fidelity_simplicity_supplyment/figure_9/jug_user_interaction.svg"; 
//char *fname = "E:/fidelity_simplicity_supplyment/figure_10/cat.svg";//archi,mecha ,complex_drawing , street , mechanical_piece , cat
//char *fname = "E:/fidelity_simplicity_supplyment/figure_11/v1.svg";
//char *fname = "E:/fidelity_simplicity_supplyment/figure_14/bag/our.svg";
//char *fname = "E:/fidelity_simplicity_supplyment/figure_15/car/our.svg";
//char *fname = "E:/fidelity_simplicity_supplyment/figure_15/house/our.svg";
//char *fname = "E:/fidelity_simplicity_supplyment/figure_15/cartoon/our.svg"; //duck
//char *fname = "E:/true2form/true2form_supplyment/InputSVGResampledWithDots/InputSVGResampledWithDots/toothpaste_resampled_black.svg";//    
//box_bump_cut_half_resampled_black
//cup_half_resampled_black
//GT_cross_resampled_black
//GT_saddle_resampled_black
//GT_spiral_resampled_black
//GT_revolve_resampled_black
//card_reader_resampled_black
//cone_smooth_resampled_black
//fighter_half_resampled_black
//fork_half_resampled_black
//toothpaste_resampled_black
//gamepadnew_half_resampled_black
//blender_half_resampled_black
//car_93_half_resampled_black
//guitar_resampled_black
//iron_half_resampled_black
//modem_boundary_resampled_black
//new_camera_resampled_black
//Pencil_holder_half_resampled_black
//plane_half_resampled_black
//ring_quarter_resampled_black
//sedan_half_resampled_black
//sewing_machine_half_resampled_black
//vacuum_half_resampled_black
//char *fname = "E:/true2form/true2form_supplyment/InputSVG/coffee_half.svg";

#include <igl/opengl/glfw/Viewer.h>
void main()
{
	////function���������ָ��ͼ
	//string combine_name;
	//combine_name.append(FILEPATH);
	//combine_name.append(FILENAME);
	//combine_name.append(".svg");
	//double file_width, file_height;
	//vector<double> file_view;
	//string file_view_str;

	//MySvgFile afile(combine_name);
	//afile.getGraphPrimitive(&file_width, &file_height, &file_view, file_view_str);
	//afile.getCrossCurve();
	//afile.testShowDiffvec();
	////function���鿴PCA + ���α������������ ���
	//string combine_name;
	//combine_name.append(FILEPATH);
	//combine_name.append(FILENAME);
	//combine_name.append(".svg");
	//double file_width, file_height;
	//vector<double> file_view;
	//string file_view_str;
	//MySvgFile afile(combine_name.c_str());
	//afile.getGraphPrimitive(&file_width, &file_height, &file_view, file_view_str);
	////afile.getCrossCurve();
	//afile.piecewiseBezieFitting();  //afile.TestPCA();  //afile.LaplacianEigenmap();  //   afile.TestAverageCurve(); //
	//afile.writefile();
	//afile.testShowDiffvec();//ע��ȥ�޸� mydisplay()����..._origin 

////---------------------------------------------//


	//tbb::parallel_for(0, 10, [](int num) {cout << num << ": Hello TBB!" << endl; });
	////�����ļ����б��ļ���������������ļ�
	//fstream file_list("file_list.txt", std::ios::out);
	//file_list.close();
	////д���ļ����б�file_list.txt
	//system("dir /a /b >> file_list.txt");





	//////method3
	string combine_name;
	combine_name.append(FILEPATH);
	combine_name.append(FILENAME);
	combine_name.append(".svg");
	double file_width, file_height;
	vector<double> file_view;
	string file_view_str;

	MyGraph agraph(combine_name.c_str());
	agraph.getGraphPrimitive(&file_width, &file_height, &file_view, file_view_str);
	agraph.getCrossCurve();
	agraph.constuctEdgeGraph(file_width, file_height, file_view_str);
	agraph.constructPointGraph();

	agraph.createGraph();
	agraph.findLines();
	agraph.divideSubGraph();
	agraph.MakeOrderedArcs();
	agraph.findPatches_new();//	findPatches();  findPatches_new

	//agraph.detectShorterEdge();
	//agraph.ConstructPatchTopology();
	//agraph.getPatchFeatures();
	//	//agraph.KmeansClustering(2);
	//	//agraph.KmedoidClustering(2);
	//	//agraph.DBSCANClusering();
	//	//agraph.FCMClustering(2, 2);
	//  //agraph.deleUndesirablePatch();
	////agraph.jointCurveTPatch_subgraph();

	//agraph.getNeiborLines(17,16, 50,20,0.1);// 20cm����  25.85������
	//agraph.Erosion_test(); 
	//agraph.localProcess();
	//agraph.patchSimilarity();

	/*0.��һ��ʼû��heap�ķ���*/
	//agraph.localProcess_two(64);
	/*1.ʹ��heap ��weight����sum��3ft����С�������򣬺ϲ�����Ӱ����С���Ǹ�patch�ϲ�*/
		/*�������ͣ����ֻʣ15%��patch, ���ƶ�С��2.5�Ĳſ��Ժϲ�������FDʱfixed pointnumber=64,������1û�á� //������������patch��width��ƽ��ֵ����Ϊ�ߵĴ���Ĳ�����*/
	//double aver_w = agraph.localProcess_heap(15, 2.5, 64, 1);	
	/*2.ʹ��heap ��weight��patch��neighbor�ϲ���Ӱ����С��similarity distance��ΪȨֵ���ϲ�����Ӱ����С���Ǹ�patch�ϲ�*/
		/*�������ͣ����ֻʣ20%��patch, ���ƶ�С��0.06�Ĳſ��Ժϲ�������FDʱfixed pointnumber=128,������1û�á�//������������patch��width��ƽ��ֵ����Ϊ�ߵĴ���Ĳ�����*/
	//double aver_w = agraph.similarityRange(20, 0.06, 128, 1);
	/*(ɾ)3.ʹ��heap ��weight��patch FD��neighbor FD�Ĳѡ��С�ģ����ϲ�������С�ĵ��Ǹ�patch�ϲ�*/
		/*�������ͣ����ʣ100��patch, ���ƶ�С��11�Ĳſ��Ժϲ�������FDʱfixed pointnumber=128,������1û�á�*/
	//agraph.shapeSimilarityRange(100, 11, 128, 1);	
	/*�ߵĴ���*/
	//agraph.attachLineProcess(aver_w, 20);

	double aver_w=agraph.ConstructMixedTopology(0.01);
	agraph.MixSimilarityRange(256, aver_w, 1, 0.02);
	//agraph.test_MixSimilarityRange(256, aver_w);

	agraph.jointCurveTPatch();
	agraph.jointPolylines();
	agraph.writeSVGFile(FILENAME, OUTPUT_PATH, file_width, file_height, file_view, file_view_str);
	agraph.globalSmoothing();
	agraph.writesmoothFile(FILENAME, OUTPUT_PATH, file_width, file_height, file_view, file_view_str);
	agraph.clearStructure();
	//�ظ��ڶ���
	agraph.getCrossCurve();
	agraph.constuctEdgeGraph(file_width, file_height, file_view_str);
	agraph.constructPointGraph();

	agraph.createGraph();
	agraph.findLines();
	agraph.divideSubGraph();
	agraph.MakeOrderedArcs();
	agraph.findPatches_new();

	 aver_w = agraph.ConstructMixedTopology(0.02);
	agraph.MixSimilarityRange(256, aver_w, 1, 0.03);
	//agraph.test_MixSimilarityRange(256, aver_w);

	agraph.globalSmoothing();
	agraph.writesmoothFile(FILENAME, OUTPUT_PATH, file_width, file_height, file_view, file_view_str);
	/////////////////////////////

	agraph.jointCurveTPatch();
	agraph.jointPolylines();
	agraph.writeSVGFile(FILENAME, OUTPUT_PATH, file_width, file_height,file_view,file_view_str);
	agraph.testShowGraph();
	//---------------------------------------------//
	/*MyGraph agraph(fname); //û�зָ��
	agraph.getGraphPrimitive();
	agraph.findPointAndEdge();
	agraph.createGraph();
	agraph.divideSubGraph();
	agraph.findPatches();
	agraph.jointCurveTPatch();
	agraph.testShowGraph();*/
	////////////////////////////////////////////////////������������

	//vector<string> nnn;
	//nnn.push_back("girl");
	//nnn.push_back("car");
	//nnn.push_back("geometry");
	//nnn.push_back("house");
	//nnn.push_back("kneel");
	//nnn.push_back("mouse");
	//nnn.push_back("ring");
	//nnn.push_back("shesh");
	//nnn.push_back("waiter");
	//for (unsigned int i = 0; i < nnn.size(); ++i)
	//{
	//	FILENAME = nnn[i];
	//	string combine_name;
	//	combine_name.append(FILEPATH);
	//	combine_name.append(FILENAME);
	//	combine_name.append(".svg");
	//	double file_width, file_height;
	//	vector<double> file_view;
	//	string file_view_str;
	//	MyGraph agraph(combine_name.c_str());
	//	agraph.getGraphPrimitive(&file_width, &file_height, &file_view, file_view_str);
	//	if (file_width == 0)
	//		file_width = 1500;
	//	if (file_height == 0)
	//		file_height = 1500;
	//	if (file_view.empty())
	//	{
	//		file_view.push_back(0); file_view.push_back(0); file_view.push_back(1500); file_view.push_back(1500);
	//		file_view_str = "0 0 1500 1500";
	//	}
	//	agraph.getCrossCurve();
	//	agraph.constuctEdgeGraph(file_width, file_height);
	//	agraph.constructPointGraph();
	//	agraph.createGraph();
	//	agraph.findLines();
	//	agraph.divideSubGraph();
	//	agraph.MakeOrderedArcs();
	//	agraph.findPatches_new();//	findPatches();  findPatches_new
	//	agraph.jointCurveTPatch();
	//	agraph.jointPolylines();
	//	agraph.writeSVGFile(FILENAME, OUTPUT_PATH, file_width, file_height, file_view, file_view_str);
	//}
	
	system("pause");
}