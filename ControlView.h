#pragma once
#ifndef CONTROLVIEW_H
#define CONTROLVIEW_H

#include "Model.h"
#include "Axespart.h"
#include "Textpart.h"
#include "Sliderbarpart.h"
#include "lutpart.h"
#include "Plotpart.h"

#include "vtkProperty.h"
// for animation
#include "vtkProgrammableFilter.h"
#include "vtkCallbackCommand.h"
// for windows render and view
#include "vtkRenderWindow.h"
#include "vtkRenderer.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkCamera.h"
#include "vtkInteractorStyleTrackballCamera.h"

void argParser(const int& argc, char* argv[], vector<string> &modelFiles,
											  vector<string> &dispFiles,
											  vector<string> &contFiles,
											  vector<string> &nodeDatafiles)
{
	if (argc == 1) {
		cout << "Type 'Univiewer /h' for more help." << endl;
		exit(0);
	}
	vector<string> *ptr_vector_name = NULL;
	for (int i = 1; i<argc; ++i) {

#ifdef __APPLE__
		if (*argv[i] == '-') {
#else
		if ((*argv[i] == '-') || (*argv[i] == '/')) {
#endif
			switch (*(argv[i] + 1)) {
			case 'm':
				ptr_vector_name = &modelFiles;
				break;
			case 'o':
				ptr_vector_name = &dispFiles;
				break;
			//case 'n':
			//	ptr_vector_name = &nodeDatafiles;
			//	break;
			case 's':
				ptr_vector_name = &nodeDatafiles;
				break;
			case 'c':
				ptr_vector_name = &contFiles;
				break;
			case 'h':
				cout << endl;
				cout << "Usage: Univiewer /m file1.xml file2.xml ... fileN.xml /o disp.dat" << endl;
				cout << endl;
				cout << "The 'fileN.xml' is a model file that is compatible with VTK API, and other programs. It is a text file using *.vtu format, but you can also add your own data." << endl;
				cout << "The 'disp.dat' is a data file that contains the displacements of nodes for all the models(file1.xml, file2.xml and so on). It is a binary file with all the data as the type of double, the detail data sequence is as follows:" << endl << endl;
				cout << "#####################################################################################################" << endl;
				cout << "time1 \n node_1_disp_x node_1_disp_y node_1_disp_z \n node_2_disp_x node_2_disp_y node_2_disp_z \n ... \n node_n_disp_x node_n_disp_y node_n_disp_z" << endl;
				cout << "time2 \n ......" << endl;
				cout << "#####################################################################################################" << endl << endl;
				cout << "NOTE: Binary files are not same as txt files which contains symbols like \\n \\t, and all the data saved as String. A Binary file contains data saved as its own type without any symbols like \\n \\t." << endl;
				exit(0);
				break;
			default:
				break;
			}
			continue;
		}
		ptr_vector_name->push_back(argv[i]);
		}
}


// Timer callback subclass
class CommandSubclass : public vtkCommand
{
public:
	vtkTypeMacro(CommandSubclass, vtkCommand);
	static vtkSmartPointer<CommandSubclass> New()
	{
		return new CommandSubclass;
	}
	void Execute(vtkObject *caller, unsigned long vtkNotUsed(eventId),
		void *vtkNotUsed(callData))
	{
		vtkRenderWindowInteractor *iren =
			static_cast<vtkRenderWindowInteractor*>(caller);
		this->ProgrammableFilter->Modified();
		iren->Render();
	}
	vtkSmartPointer<vtkProgrammableFilter> ProgrammableFilter;
};

// real callback function
void TimerCallback(void* arguments);
void WindowModifiedCallback(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData);
void KeypressCallback(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData);

class ControlView
{
protected:
	// for the centering point
	vtkSmartPointer<vtkActor> actor;
	vtkSmartPointer<vtkDataSetMapper> mapper;

	// for windows control
	vtkSmartPointer<vtkRenderWindow> renderWindow;
	vtkSmartPointer<vtkRenderer> renderer;
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor;
	vtkSmartPointer<vtkInteractorStyleTrackballCamera> style;

	// for some part of graphics
	shared_ptr<Axesline> axesline;
	shared_ptr<Axesframe> axesframe;
	shared_ptr<Sliderbar> sliderbar;
	shared_ptr<CurrentTimer> currenttimer;
	shared_ptr<LookUpTable> lookuptable;

	shared_ptr<CommandText> cText;
	vector<shared_ptr<CommandText> > cTextBodies;
	// for animation and UI
	vtkSmartPointer<vtkProgrammableFilter> programmableFilter;

	// main part of visulization
	vector<shared_ptr<Model>> pModels;
	vector<shared_ptr<ContactData>> pContacts;
	//////////////////////////////////
	//////////////////////////////////

	bool play, stepPlay;
	bool ShowMarker;
	bool ShowMesh;
    bool ShowLabel;

	vector<string> dispFiles;
	vector<string> contFiles;
	vector<string> nodeDataFiles;
	//int nextButtonState, prevButtonState;

public:
	virtual ~ControlView() {}
	ControlView(): play(false), stepPlay(false), ShowMarker(true), ShowMesh(true), ShowLabel(false) {};
	static shared_ptr<ControlView> New()
	{
		return make_shared<ControlView>();
	}
	
	inline bool & IsPlay() { return play; }
	inline bool & IsStepPlay() { return stepPlay; }
	inline bool & IsShowMarker() { return ShowMarker; }
	inline bool & IsShowMesh() { return ShowMesh; }
    inline bool & IsShowLabel() { return ShowLabel; }

	vector<shared_ptr<Model>> & getModels() { return pModels; }
	vector<shared_ptr<ContactData>> &getContactData() { return pContacts; }
	vtkSmartPointer<vtkProgrammableFilter> getProgrammableFilter() { return programmableFilter; }

	shared_ptr<CommandText> getCommandText() { return cText; }
	vector<shared_ptr<CommandText>> &getCommandTextBodies() { return cTextBodies; }
	shared_ptr<CurrentTimer> getCurrentTimer() {return currenttimer;}
	shared_ptr<Axesline> getAxesline() { return axesline; }
	shared_ptr<Sliderbar> getSliderbar() { return sliderbar; }
	shared_ptr<LookUpTable> getLookuptable() { return lookuptable; }

    vtkSmartPointer<vtkRenderer> getRenderer(){return renderer;}
    vtkSmartPointer<vtkRenderWindowInteractor> getRenderWindowInteractor(){return renderWindowInteractor;}

	virtual int inputModelfiles(const int& argc,  char* argv[]) 
	{
		vector<string> modelFiles;
		argParser(argc, argv, modelFiles, dispFiles, contFiles, nodeDataFiles);
		
		pModels.resize(modelFiles.size());
		for (unsigned i = 0; i< (unsigned)modelFiles.size(); ++i) {
			pModels[i] = Model::New();
			pModels[i]->setOffset(Model::nodeNums* 3);

			if (modelFiles[i].substr(modelFiles[i].size() - 3).compare(".md") == 0) {
				pModels[i]->ReadTxtModel(modelFiles[i]);
			} else if (modelFiles[i].substr(modelFiles[i].size() - 3).compare("xml") == 0) {
				pModels[i]->ReadXmlModel(modelFiles[i]);
			}
			
			//pModels[i]->getActor()->GetProperty()->SetColor(0.9, 0.9, 0.9);
			//pModels[i]->getActor()->GetProperty()->SetOpacity(1.0);
			pModels[i]->getActor()->GetProperty()->SetEdgeColor(0.0, 0.0, 0.0);
			pModels[i]->getActor()->GetProperty()->EdgeVisibilityOn();

			if (pModels[i]->getFEMesh()->getUGrid()->GetCellType(0) == 4)
				pModels[i]->getActor()->GetProperty()->SetLineWidth(2.5);
			else
				pModels[i]->getActor()->GetProperty()->SetLineWidth(1.0);
			
			pModels[i]->getActor()->SetScale(1.0);

			renderer->AddActor(pModels[i]->getActor());
		}

		if (!this->dispFiles.empty()) {
			readDispfile(dispFiles, pModels);

			sliderbar = Sliderbar::New();
			sliderbar->setSliderBar(renderWindowInteractor, Model::step, Model::stepNum, this->play, this->stepPlay);

			currenttimer = CurrentTimer::New();
			currenttimer->setTextActor(renderWindow);
			renderer->AddActor2D(currenttimer->getTextActor());
		}

		if (!this->contFiles.empty()) {
			readContfile(contFiles[0], pContacts, pModels);
			for (auto& pcontact : pContacts) {
				pcontact->InitializeUGrid();

				pcontact->getPActor()->GetProperty()->SetLineWidth(2.0);
				pcontact->getNActor()->GetProperty()->SetLineWidth(2.0);
				pcontact->getPActor()->GetProperty()->SetPointSize(6.0);
				pcontact->getNActor()->GetProperty()->SetPointSize(6.0);
				pcontact->getPActor()->GetProperty()->EdgeVisibilityOff();
				pcontact->getNActor()->GetProperty()->EdgeVisibilityOff();

				renderer->AddActor(pcontact->getPActor());
				renderer->AddActor(pcontact->getNActor());
			}
		}

		if (!this->nodeDataFiles.empty()) {
			readNodeDatafile(nodeDataFiles[0], pModels);
		}



		return 0;
	}
	
	virtual void setAnimationMethod( void(*f)(void*)) {
		if (!dispFiles.empty()) {
			programmableFilter->SetExecuteMethod(f, this);
			renderWindowInteractor->CreateRepeatingTimer(10);
			vtkSmartPointer<CommandSubclass> timerCallback = vtkSmartPointer<CommandSubclass>::New();
			timerCallback->ProgrammableFilter = programmableFilter;
			renderWindowInteractor->AddObserver(vtkCommand::TimerEvent, timerCallback);
		}
	}

	virtual void setKeyboardMethod(void (*f)(vtkObject* , long unsigned int vtkNotUsed(eventId), void* clientData, void* vtkNotUsed(callData))){
		vtkSmartPointer<vtkCallbackCommand> keypressCallback = vtkSmartPointer<vtkCallbackCommand>::New();
		keypressCallback->SetCallback(f);
		keypressCallback->SetClientData(this);
		renderWindowInteractor->AddObserver(vtkCommand::KeyPressEvent, keypressCallback);
	}

	virtual void setWindowMethod(void (*f)(vtkObject* caller, long unsigned int vtkNotUsed(eventId), void* clientData, void* vtkNotUsed(callData))){
		vtkSmartPointer<vtkCallbackCommand> windowCallback = vtkSmartPointer<vtkCallbackCommand>::New();
		windowCallback->SetCallback(f);
		windowCallback->SetClientData(this);
		renderWindow->AddObserver(vtkCommand::ModifiedEvent, windowCallback);
	}

	void playThisStep(unsigned &step)
	{
		stringstream ss("");
		ss << "Current Time = " << Model::stepCollection[step];
		currenttimer->getTextActor()->SetInput(ss.str().c_str());

		sliderbar->getSliderRep()->SetValue(step);

		for (auto &pmodel : pModels) {
			pmodel->updateDisp(step);
		}

		for (auto &pcontact : pContacts) {
			pcontact->UpdateUGrid(step);
		}
	}

	virtual void Update() {
			
		if (stepPlay) {
			play = false;
			sliderbar->getPlayButton()->setState(0);
			playThisStep(Model::step);
		}
		
		if (play) {
			Model::step = (Model::step == Model::stepNum - 1) ? 0 : Model::step + 1;
			playThisStep(Model::step);
		}
	}

	virtual void setRender() {
		// create an environment actor (center point)
		vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
		const double p[3] = { 0.0, 0.0, 0.0 };
		// Create the topology of the point (a vertex)
		vtkSmartPointer<vtkCellArray> vertices = vtkSmartPointer<vtkCellArray>::New();
		vtkIdType pid[1];
		pid[0] = points->InsertNextPoint(p);
		vertices->InsertNextCell(1, pid);
		// Create a polydata object
		vtkSmartPointer<vtkPolyData> point = vtkSmartPointer<vtkPolyData>::New();
		// Set the points and vertices we created as the geometry and topology of the polydata
		point->SetPoints(points);
		point->SetVerts(vertices);
		
		programmableFilter = vtkSmartPointer<vtkProgrammableFilter>::New();
		programmableFilter->SetInputData(point);
		
		mapper = vtkSmartPointer<vtkDataSetMapper>::New();
		mapper->SetInputConnection(programmableFilter->GetOutputPort());
		actor = vtkSmartPointer<vtkActor>::New();
		actor->SetMapper(mapper);
		actor->GetProperty()->SetColor(0.0, 0.0, 0.0);
		actor->GetProperty()->SetPointSize(1.0);

		// Create a renderer, render window, and interactor
		renderer = vtkSmartPointer<vtkRenderer>::New();
		renderer->SetNearClippingPlaneTolerance(1e-4);
		renderer->GetActiveCamera()->SetClippingRange(1e-4, 1000);
		renderer->SetAmbient(1.0, 1.0, 1.0);
		renderer->SetLightFollowCamera(1);

		renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
		renderWindow->AddRenderer(renderer);
		renderWindow->SetSize(800, 640);
		renderWindow->SetWindowName("MLV 2.0");
        renderWindow->Render();
		renderWindow->PointSmoothingOn();
		renderWindow->LineSmoothingOn();
		renderWindow->PolygonSmoothingOn();

		renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
		renderWindowInteractor->SetRenderWindow(renderWindow);
		style = vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
		renderWindowInteractor->SetInteractorStyle(style);

		// Initialize must be called prior to creating timer events.
		renderWindowInteractor->Initialize();

	}

	virtual void Display() {

		// command text
		cText = CommandText::New();
		cText->setRenderWindowInteractor(renderWindowInteractor);
		int *winSize = renderWindowInteractor->GetRenderWindow()->GetSize();
		cText->setCommandTextContent("System View", 1.0, 1.0, 1.0, 0.5, 1);
		//cText->setTextSizePosition(10, winSize[1] - 30, 20);
		vtkSmartPointer<systemReleaseTextCallback> systemRelease = vtkSmartPointer<systemReleaseTextCallback>::New();
		systemRelease->renderWindowInteractor = renderWindowInteractor;
		cText->setTextCallback(systemRelease);

		std::stringstream ss("");
		for (unsigned i = 0; i < pModels.size(); ++i) {
			shared_ptr<CommandText> ct = CommandText::New();
			ct->setRenderWindowInteractor(renderWindowInteractor);
			ss.str("");
			ss.clear();
			ss << "* Body_" << i;
			ct->setCommandTextContent(ss.str(), 1.0, 1.0, 1.0, 0.5, 1);
			//ct->setTextSizePosition(18, winSize[1] - 30 - 15* (i+1), 15);
			vtkSmartPointer<releaseTextCallback> release = vtkSmartPointer<releaseTextCallback>::New();
			ct->setTextCallback(release);
			cTextBodies.push_back(ct);
			cText->addLeafNode(ct->getTextWidget());
		}


		// axesline part
		axesline = Axesline::New();
		axesline->setAxesActor();
		renderer->AddActor(axesline->getAxesActor());

		// axesframe part
		axesframe = Axesframe::New();
		axesframe->setAxesWidget(renderWindowInteractor);

		// labelnodes part
		for (auto& pmodel : pModels) {
			pmodel->setLabelnode();
			renderer->AddActor2D(pmodel->getLabelactor());
			pmodel->getLabelactor()->VisibilityOff();
		}
		
		if (!nodeDataFiles.empty()) {
			// lookuptable part
			lookuptable = LookUpTable::New();
			lookuptable->setScalars(pModels);
			renderer->AddActor2D(lookuptable->getScalarBar());
		}

		// Add the actor to the scene
		renderer->AddActor(actor);
		
		renderer->GradientBackgroundOn();
		renderer->SetBackground2(13.0 / 255.0, 71.0 / 255.0, 161.0 / 255.0);
		renderer->SetBackground(144.0 / 255.0, 202.0 / 255.0, 249.0 / 255.0);
		// Render and interact

		renderWindowInteractor->Start();

		//plotpart p;
		//p.plot(renderWindow, renderer);
	}

};


void TimerCallback(void* arguments)
{
	ControlView* pCtr = static_cast<ControlView*>(arguments);
	pCtr->Update();
}

void WindowModifiedCallback(vtkObject* caller, long unsigned int vtkNotUsed(eventId), void* clientData, void* vtkNotUsed(callData))
{
	
	vtkRenderWindow* window = static_cast<vtkRenderWindow*>(caller);
	int* windowSize = window->GetSize();
	ControlView* pCtr = static_cast<ControlView*>(clientData);
	
	if(Model::stepNum!=0){

		pCtr->getCurrentTimer()->getTextActor()->SetPosition(windowSize[0] - 180, windowSize[1] - 25);

		pCtr->getSliderbar()->getSliderRep()->GetPoint1Coordinate()->SetValue(windowSize[0] - 260 - windowSize[0]/4.0, windowSize[1] - 15);
		pCtr->getSliderbar()->getSliderRep()->GetPoint2Coordinate()->SetValue(windowSize[0] - 260, windowSize[1] - 15);

		int width = 20, height = 20;
		pCtr->getSliderbar()->getPlayButton()->setButtonSizePosition(width, height, windowSize[0] - 235 - width, windowSize[1] - 5 - height);
		pCtr->getSliderbar()->getNextButton()->setButtonSizePosition(width, height, windowSize[0] - 210 - width, windowSize[1] - 5 - height);
		pCtr->getSliderbar()->getPrevButton()->setButtonSizePosition(width, height, windowSize[0] - 185 - width, windowSize[1] - 5 - height);
	}

	pCtr->getCommandText()->setTextSizePosition(10, windowSize[1] - 30, 20);
	for (unsigned i = 0; i < pCtr->getModels().size(); ++i) {
		pCtr->getCommandTextBodies()[i]->setTextSizePosition(10, windowSize[1] - 30 - 15 * (i+1), 15);
	}

	if (pCtr->getLookuptable()) {
		double w = windowSize[0] * 70.0 / 800.0;
		double h = windowSize[1] * 300.0 / 640.0;
		pCtr->getLookuptable()->getScalarBar()->SetPosition(windowSize[0] - w - 10, 10);
		pCtr->getLookuptable()->getScalarBar()->SetPosition2(w, h);
	}

}

void KeypressCallback(vtkObject* caller, long unsigned int vtkNotUsed(eventId), void* clientData, void* vtkNotUsed(callData))
{

	vtkRenderWindowInteractor *iren =
		static_cast<vtkRenderWindowInteractor*>(caller);

	ControlView* pCtr = static_cast<ControlView*>(clientData);

    if(Model::stepNum!=0){
	    if (strcmp(iren->GetKeySym(), "space") == 0){
	        pCtr->IsPlay() = !pCtr->IsPlay();
			if (pCtr->IsPlay()) {
				pCtr->getSliderbar()->getPlayButton()->setState(1);
				pCtr->IsStepPlay() = false;
			}
			else {
				pCtr->getSliderbar()->getPlayButton()->setState(0);
			}
				
	    }

		if (strcmp(iren->GetKeySym(), "b") == 0) {
			pCtr->IsPlay() = false;
			pCtr->getSliderbar()->getPlayButton()->setState(0);
			Model::step = (Model::step == Model::stepNum - 1) ? Model::stepNum - 1 : Model::step + 1;
			pCtr->IsStepPlay() = true;
		}

		if (strcmp(iren->GetKeySym(), "v") == 0) {
			pCtr->IsPlay() = false;
			pCtr->getSliderbar()->getPlayButton()->setState(0);
			Model::step = (Model::step == 0) ? 0 : Model::step - 1;
			pCtr->IsStepPlay() = true;
		}
    }
	if (strcmp(iren->GetKeySym(), "i") == 0) {
		if (pCtr->IsShowMarker())
            pCtr->getAxesline()->getAxesActor()->VisibilityOff();
		else
			pCtr->getAxesline()->getAxesActor()->VisibilityOn();

		pCtr->IsShowMarker() = !pCtr->IsShowMarker();
	}

	if (strcmp(iren->GetKeySym(), "u") == 0) {
		if (pCtr->IsShowMesh()) {
			for(auto &pmodel: pCtr->getModels())
				pmodel->getActor()->GetProperty()->EdgeVisibilityOff();
		}
		else {
			for (auto &pmodel : pCtr->getModels())
				pmodel->getActor()->GetProperty()->EdgeVisibilityOn();
		}
			
		pCtr->IsShowMesh() = !pCtr->IsShowMesh();
	}

	if (strcmp(iren->GetKeySym(), "l") == 0) {

		if (pCtr->IsShowLabel()) {
			for(auto& pmodel : pCtr->getModels())
				pmodel->getLabelactor()->VisibilityOff();
		}
		else {
			for (auto& pmodel : pCtr->getModels())
				pmodel->getLabelactor()->VisibilityOn();
		}

		pCtr->IsShowLabel() = !pCtr->IsShowLabel();
	}

	iren->Render();
}
#endif //CONTROLVIEW_H