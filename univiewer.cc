#include "ControlView.h"
int main(int argc, char *argv[])
{

	shared_ptr<ControlView> pControlView = ControlView::New();
	pControlView->setRender();

	string arg1 = "/r";
	string arg2 = "c:/Users/CNJISHI10/WorkSpace/repository/msdtk/build/bin/testSolid.dat";
	argv[1] = (char*)arg1.c_str();
	argv[2] = (char*)arg2.c_str();
	argc = 3;

	pControlView->inputModelfiles(argc, argv);

	pControlView->setAnimationMethod(DEFAULT_TIMERCALLBACK);
	pControlView->setKeyboardMethod(DEFAULT_KEYPRESSCALLBACK);
	pControlView->setWindowMethod(DEFAULT_WINDOWCALLBACK);
	pControlView->Display();
	
	return EXIT_SUCCESS;
}


