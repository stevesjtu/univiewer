//#define vtkRenderingCore_AUTOINIT 4(vtkInteractionStyle,vtkRenderingFreeType,vtkRenderingFreeTypeOpenGL,vtkRenderingOpenGL)
//#define vtkRenderingVolume_AUTOINIT 1(vtkRenderingVolumeOpenGL)
#include "ControlView_addText.h"
int main(int argc, char *argv[])
{

	vector<string> modelFiles;
	vector<string> dispFiles;
	
	shared_ptr<Model> pModel = Model::New();
	shared_ptr<ControlView_addText> pControlView = ControlView_addText::New();
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
    pControlView->setfileName(&modelFiles, &dispFiles);
    pControlView->AddText();
        
	pControlView->setKeyboardMethod(DEFAULT_KEYPRESSCALLBACK);
	pControlView->setWindowMethod(DEFAULT_WINDOWCALLBACK);

	pControlView->Display();

	return EXIT_SUCCESS;
}


