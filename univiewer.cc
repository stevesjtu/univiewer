#include "ControlView.h"
int main(int argc, char *argv[])
{

	shared_ptr<ControlView> pControlView = ControlView::New();
	pControlView->setRender();

	pControlView->inputModelfiles(argc, argv);

	pControlView->setAnimationMethod(DEFAULT_TIMERCALLBACK);
	pControlView->setKeyboardMethod(DEFAULT_KEYPRESSCALLBACK);
	pControlView->setWindowMethod(DEFAULT_WINDOWCALLBACK);
	pControlView->Display();
	
	return EXIT_SUCCESS;
}


