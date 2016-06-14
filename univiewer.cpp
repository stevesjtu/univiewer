//#define vtkRenderingCore_AUTOINIT 4(vtkInteractionStyle,vtkRenderingFreeType,vtkRenderingFreeTypeOpenGL,vtkRenderingOpenGL)
//#define vtkRenderingVolume_AUTOINIT 1(vtkRenderingVolumeOpenGL)
#include "ControlView_Converter.h"
int main(int argc, char *argv[])
{

	shared_ptr<Model> pModel = Model::New();
	shared_ptr<ControlView_Converter> pControlView = ControlView_Converter::New();
	//pControlView->ConvertPLY2UGrid(argc, argv);
	pControlView->ConvertUGrid2PLY(argc, argv);
	//pControlView->ModifyUGrid(argc, argv);

	pControlView->setMainActor();
	pControlView->setAxesActor();
	pControlView->setTextActor();
    pControlView->setLabelActor();
    
	pControlView->setRender();

	pControlView->setKeyboardMethod(DEFAULT_KEYPRESSCALLBACK);
	pControlView->setWindowMethod(DEFAULT_WINDOWCALLBACK);

	pControlView->Display();

	return EXIT_SUCCESS;
}


