#include "ControlView.h"
int main(int argc, char *argv[])
{

	vector<string> modelFiles;
	vector<string> dispFiles;
	
	shared_ptr<Model> pModel = Model::New();
	shared_ptr<ControlView> pControlView = ControlView::New();
	pControlView->inputModelfiles(modelFiles, dispFiles, argc, argv);

	pControlView->setMainActor();

	// add all the configure
	int prt = AXESLINE_PART | AXESFRAME_PART | CURRENTTIMER_PART | SLIDEBAR_PART | LOOKUPTABLE_PART; 

	pControlView->setRender();

	if (!dispFiles.empty()) {
		
		pControlView->setContent(pModel);
		pControlView->setAnimationMethod(DEFAULT_TIMERCALLBACK, dispFiles);
		
	}
	else {
		prt ^= SLIDEBAR_PART;  // delete the SLIDERBAR_PART configure
	}
    
    // derived class member function
    //pControlView->setfileName(&modelFiles, &dispFiles);
    //pControlView->AddText();
        
	pControlView->setKeyboardMethod(DEFAULT_KEYPRESSCALLBACK);
	pControlView->setWindowMethod(DEFAULT_WINDOWCALLBACK);

	pControlView->Display(prt);
	
	return EXIT_SUCCESS;
}


