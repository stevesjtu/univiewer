#pragma once
#ifndef UNIVIEWER_H
#define UNIVIEWER_H
#include "ControlView.h"

#ifndef USING_DLL
    #define UNIVIEWER_API
#else
    #ifdef UNIVIEWER_EXPORTS
    #define UNIVIEWER_API __declspec(dllexport)
    #else
    #define UNIVIEWER_API __declspec(dllimport)
    #endif
#endif

#ifdef __APPLE__
class Univiewer
#else
class UNIVIEWER_API Univiewer
#endif
{
protected:
	vector<string> modelFiles;

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
	int plotModel(vector<TriangleMesh> &mesh);
	int plotModel(vector<TriangleMesh> &mesh, const vector<pair<unsigned, unsigned>> &pairs);
};

int Univiewer::plotModel(int argc, char *argv[])
{
	vector<string> modelFiles;

	shared_ptr<ControlView> pControlView = ControlView::New();
	pControlView->inputModelfiles(modelFiles, argc, argv);

	pControlView->setMainActor();
	pControlView->setAxesActor();
	pControlView->setTextActor();
	pControlView->setLabelActor();

	pControlView->setRender();

	pControlView->setKeyboardMethod(DEFAULT_KEYPRESSCALLBACK);

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
    shared_ptr<ControlView> pControlView = ControlView::New();
    pControlView->inputModel(elem, node);

    pControlView->setMainActor();
    pControlView->setAxesActor();
    pControlView->setTextActor();
    pControlView->setLabelActor();
    
    pControlView->setRender();
    
    pControlView->setKeyboardMethod(DEFAULT_KEYPRESSCALLBACK);
    
    pControlView->Display();
    return 0;

}

int Univiewer::plotModel(TriangleMesh &mesh)
{
	return plotModel(*mesh.pelem, *mesh.pnode);
}

int Univiewer::plotModel(vector<TriangleMesh> &mesh)
{
	shared_ptr<ControlView> pControlView = ControlView::New();

	pControlView->inputModel(mesh);

	pControlView->setMainActor();
	pControlView->setAxesActor();
	pControlView->setTextActor();
	pControlView->setLabelActor();

	pControlView->setRender();

	pControlView->setKeyboardMethod(DEFAULT_KEYPRESSCALLBACK);

	pControlView->Display();
	return 0;
}

int Univiewer::plotModel(vector<TriangleMesh> &mesh, const vector<pair<unsigned, unsigned>> &pairs)
{
	shared_ptr<ControlView> pControlView = ControlView::New();

	pControlView->inputModel(mesh, pairs);

	pControlView->setMainActor();
	pControlView->setAxesActor();
	pControlView->setTextActor();
	pControlView->setLabelActor();

	pControlView->setRender();

	pControlView->setKeyboardMethod(DEFAULT_KEYPRESSCALLBACK);

	pControlView->Display();
	return 0;
}

#endif // UNIVIEWER_H