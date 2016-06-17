#pragma once
#ifndef UNIVIEWER_H
#define UNIVIEWER_H
#include "ControlView.h"
#include "Model.h"

#ifdef __APPLE__
    #define UNIVIEWER_API
#else
    #ifdef UNIVIEWER_EXPORTS
    #define UNIVIEWER_API __declspec(dllexport)
    #else
    #define UNIVIEWER_API __declspec(dllimport)
    #endif
#endif

struct TriangleMesh
{
	TriangleMesh() {};
	TriangleMesh(MatrixXu &elem, MatrixXd &node) {
		pelem = &elem;
		pnode = &node;
	}
	MatrixXu *pelem;
	MatrixXd *pnode;
};


#ifdef __APPLE__
class Univiewer
#else
class UNIVIEWER_API Univiewer
#endif
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
    int plotModel(char file[]);
    int plotModel(MatrixXu &elem, MatrixXd &node);
	int plotModel(TriangleMesh &mesh);
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

int Univiewer::plotModel(char file[])
{
    char xx[] = "x";
    char m[] = {'/','m'};
    char *argv[] = {xx, m, file};
    return plotModel(3, argv);
}

int Univiewer::plotModel(MatrixXu &elem, MatrixXd &node)
{

    shared_ptr<Model> pModel = Model::New();
    shared_ptr<ControlView> pControlView = ControlView::New();
    pControlView->inputModel(elem, node);

    pControlView->setMainActor();
    pControlView->setAxesActor();
    pControlView->setTextActor();
    pControlView->setLabelActor();
    
    pControlView->setRender();
    
//    if (!dispFiles.empty()) {
//        pControlView->setContent(pModel);
//        pControlView->setAnimationMethod(DEFAULT_TIMERCALLBACK, dispFiles);
//        pControlView->setSliderBar();
//    }
    
    pControlView->setKeyboardMethod(DEFAULT_KEYPRESSCALLBACK);
    pControlView->setWindowMethod(DEFAULT_WINDOWCALLBACK);
    
    pControlView->Display();
    return 0;

}

int Univiewer::plotModel(TriangleMesh &mesh)
{
	return plotModel(*mesh.pelem, *mesh.pnode);
}
#endif // UNIVIEWER_H