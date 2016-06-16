#pragma once
#ifndef UNIVIEWER_H
#define UNIVIEWER_H
#include "ControlView.h"
#include "Model.h"

#ifdef UNIVIEWER_EXPORTS
#define UNIVIEWER_API __declspec(dllexport)
#else
#define UNIVIEWER_API __declspec(dllimport)
#endif

class UNIVIEWER_API Univiewer
{
protected:
	vector<string> modelFiles;
	vector<string> dispFiles;

	shared_ptr<Model> pModel;
	shared_ptr<ControlView> pControlView;

public:
	Univiewer() {};
	virtual ~Univiewer() {};
	static shared_ptr<Univiewer> New() {
		return make_shared<Univiewer>();
	}
	
	int plotModel(int argc, char *argv[]);
};

int Univiewer::plotModel(int argc, char *argv[])
{
	vector<string> modelFiles;
	vector<string> dispFiles;

	shared_ptr<Model> pModel = Model::New();
	shared_ptr<ControlView> pControlView = ControlView::New();
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

	pControlView->setKeyboardMethod(DEFAULT_KEYPRESSCALLBACK);
	pControlView->setWindowMethod(DEFAULT_WINDOWCALLBACK);

	pControlView->Display();

	return EXIT_SUCCESS;
}
#endif // UNIVIEWER_H