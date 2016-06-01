//#define vtkRenderingCore_AUTOINIT 4(vtkInteractionStyle,vtkRenderingFreeType,vtkRenderingFreeTypeOpenGL,vtkRenderingOpenGL)
//#define vtkRenderingVolume_AUTOINIT 1(vtkRenderingVolumeOpenGL)
#include "ControlView_PLYfileSeries.h"
#include "Model.h"

//#define FILED(var) exportMat(var, #var)
//#pragma comment(linker, "/subsystem:\"windows\" /entry:\"mainCRTStartup\"")


//template<typename T>
//int exportMat(const T &matr, const string &file)
//{
//	string fileplus(file + ".txt");
//	ofstream outfile(fileplus);
//	if (!outfile.is_open()) {
//		cout << "Can not create file." << endl;
//		return -1;
//	}
//	// Opened
//	outfile << setprecision(18);
//	outfile << matr;
//	outfile.close();
//	return 0;
//}

int main(int argc, char *argv[])
{

	vector<string> modelFiles;
	vector<string> dispFiles;
	
	shared_ptr<Model> pModel = Model::New();
	//shared_ptr<ControlView_PLYfileSeries> pControlView = ControlView_PLYfileSeries::New();
	shared_ptr<ControlView> pControlView = ControlView::New();
	//pControlView->ConvertPLY2Unstrgrid(argc, argv);
	pControlView->inputModelfiles(modelFiles, dispFiles, argc, argv);

	pControlView->setMainActor();
	pControlView->setAxesActor();
	pControlView->setTextActor();
    pControlView->setLabelActor();
    
	pControlView->setRender();

	if (!dispFiles.empty()) {
		pControlView->setContent(pModel);
		pControlView->setAnimationMethod(DEFAULT_TIMERCALLBACK, dispFiles);
		pControlView->setSliderBar();
	}
    
    // derived class member function
    //pControlView->setfileName(&modelFiles, &dispFiles);
    //pControlView->AddText();
    //
    
	pControlView->setKeyboardMethod(DEFAULT_KEYPRESSCALLBACK);
	pControlView->setWindowMethod(DEFAULT_WINDOWCALLBACK);

	pControlView->Display();

	return EXIT_SUCCESS;
}


