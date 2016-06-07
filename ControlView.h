#pragma once
#ifndef CONTROLVIEW_H
#define CONTROLVIEW_H
#include "vtkXMLUnstructuredGridReader.h"
#include "vtkUnstructuredGrid.h"
#include "vtkAppendFilter.h"
#include "vtkGeometryFilter.h"
#include "vtkUnstructuredGridGeometryFilter.h"
#include "vtkSmartPointer.h"
#include "vtkSphereSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkLabeledDataMapper.h"
#include "vtkPointData.h"
#include "vtkDataSetMapper.h"
#include "vtkActor.h"
#include "vtkProperty.h"
#include "vtkExtractEdges.h"
#include "vtkRenderWindow.h"
//#include "vtkCocoaRenderWindow.h"
#include "vtkRenderer.h"
#include "vtkRendererCollection.h"
#include "vtkActorCollection.h"
#include "vtkAxesActor.h"
#include "vtkOrientationMarkerWidget.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkCamera.h"
#include "vtkProgrammableFilter.h"
#include "vtkCallbackCommand.h"
#include "vtkCommand.h"
#include "vtkInteractorStyleTrackballCamera.h"

#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkUnsignedCharArray.h"
#include "vtkPoints.h"
#include "vtkLine.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkImageData.h"
#include "vtkCoordinate.h"
#include "vtkButtonWidget.h"
#include "vtkTexturedButtonRepresentation2D.h"

#include "vtkTextActor.h"
#include "vtkTextProperty.h"
#include "vtkTextRepresentation.h"

#include "vtkProperty2D.h"
#include "vtkWidgetEvent.h"
#include "vtkWidgetEventTranslator.h"
#include "vtkSliderRepresentation2D.h"
#include "vtkSliderWidget.h"
#include "Model.h"

#define DEFAULT_TIMERCALLBACK TimerCallback
#define DEFAULT_KEYPRESSCALLBACK KeypressCallbackFunction
#define DEFAULT_WINDOWCALLBACK WindowModifiedCallback

// Timer callback proxy
void TimerCallbackFunction(vtkObject* caller, long unsigned int vtkNotUsed(eventId), void* clientData, void* vtkNotUsed(callData));

// real callback function
void TimerCallback(void* arguments);
void WindowModifiedCallback(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData);
void KeypressCallbackFunction(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData);


class vtkSliderCallback : public vtkCommand
{
protected:

public:

    static vtkSliderCallback *New()
    {
        return new vtkSliderCallback;
    }
    virtual void Execute(vtkObject *caller, unsigned long, void*)
    {
        vtkSliderWidget *sliderWidget = reinterpret_cast<vtkSliderWidget*>(caller);
        *pst = (unsigned)static_cast<vtkSliderRepresentation *>(sliderWidget->GetRepresentation())->GetValue();
    }
    unsigned *pst;
};

void CreateImagePause(vtkSmartPointer<vtkImageData> image)
{
    // Specify the size of the image data
    image->SetDimensions(10,10,1);
    image->AllocateScalars(VTK_UNSIGNED_CHAR,4);
    int* dims = image->GetDimensions();

    // Fill the image with
    for (int y = 0; y < dims[1]; y++)
    {
        for (int x = 0; x < dims[0]; x++)
        {
            unsigned char* pixel =
            static_cast<unsigned char*>(image->GetScalarPointer(x,y,0));
            if(x < 4)
            {
                pixel[0] = 255;
                pixel[1] = 255;
                pixel[2] = 255;
                pixel[3] = 50;
            }
            else if(x > 5)
            {
                pixel[0] = 255;
                pixel[1] = 255;
                pixel[2] = 255;
                pixel[3] = 50;
            }
            else{
                pixel[3] = 0;
            }
        }
    }
}

void CreateImagePlay(vtkSmartPointer<vtkImageData> image)
{
    // Specify the size of the image data
    image->SetDimensions(200,200,1);
    image->AllocateScalars(VTK_UNSIGNED_CHAR,4);
    int* dims = image->GetDimensions();

    // Fill the image with
    for (int x = 0; x < dims[0]; x++) {
        for(int y = 0; y< x/2; y++){
            unsigned char* pixel = static_cast<unsigned char*>(image->GetScalarPointer(x,y,0));
            pixel[3] = 0;
        }
        for(int y = dims[1]-x/2; y< dims[1]; y++){
            unsigned char* pixel = static_cast<unsigned char*>(image->GetScalarPointer(x,y,0));
            pixel[3] = 0;
        }
        for (int y = x/2; y < dims[1]-x/2; y++) {
            unsigned char* pixel = static_cast<unsigned char*>(image->GetScalarPointer(x,y,0));
            pixel[0] = 255;
            pixel[1] = 255;
            pixel[2] = 255;
            pixel[3] = 50;
        }
    }
}


class ControlView
{
protected:
	vtkSmartPointer<vtkActor> actor;
	vtkSmartPointer<vtkActor> axesActor;
	vtkSmartPointer<vtkTextActor> textActor;
    vtkSmartPointer<vtkActor2D> labelActor;

    vector<vtkSmartPointer<vtkXMLUnstructuredGridReader> > ugridReaders;
    vtkSmartPointer<vtkLabeledDataMapper> labelMapper;
	vtkSmartPointer<vtkDataSetMapper> mapper;

	vtkSmartPointer<vtkRenderWindow> renderWindow;
	vtkSmartPointer<vtkRenderer> renderer;
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor;
	vtkSmartPointer<vtkInteractorStyleTrackballCamera> style;
	vtkSmartPointer<vtkAxesActor> axes;
	vtkSmartPointer<vtkOrientationMarkerWidget> widget;
	vtkSmartPointer<vtkPolyData> linesPolyData;
	vtkSmartPointer<vtkPoints> pts;
	vtkSmartPointer<vtkLine> line0, line1, line2;
	vtkSmartPointer<vtkCellArray> lines;
	vtkSmartPointer<vtkUnsignedCharArray> colors;
	vtkSmartPointer<vtkPolyDataMapper> axesMapper;

	vtkSmartPointer<vtkProgrammableFilter> programmableFilter;
	vtkSmartPointer<vtkCallbackCommand> timerCallback;
	vtkSmartPointer<vtkCallbackCommand> keypressCallback;
	vtkSmartPointer<vtkCallbackCommand> windowCallback;

    vtkSmartPointer<vtkSliderRepresentation2D> sliderRep;
    vtkSmartPointer<vtkSliderWidget> sliderWidget;
    vtkSmartPointer<vtkSliderCallback> SliderCallback;

    vtkSmartPointer<vtkImageData> imagePlay;
    vtkSmartPointer<vtkImageData> imagePause;
    vtkSmartPointer<vtkTexturedButtonRepresentation2D> buttonRepresentation;
    vtkSmartPointer<vtkButtonWidget> buttonWidget;

	shared_ptr<Model> pModel;
	//////////////////////////////////
	//////////////////////////////////
	unsigned stepNum;
	unsigned step;
	bool play;
	bool ShowMarker;
	bool ShowMesh;
    bool ShowLabel;

	void parser(int argc, char** argv, vector<string> &files1, vector<string> &files2)
	{
		vector<string> *ptr_vector_name = NULL;
		for (int i = 1; i<argc; ++i) {

			if ((*argv[i] == '-') || (*argv[i] == '/')) {
				switch (*(argv[i] + 1)) {
				case 'm':
					ptr_vector_name = &files1;
					break;
				case 'o':
					ptr_vector_name = &files2;
					break;
				default:
					break;
				}
				continue;
			}
			ptr_vector_name->push_back(argv[i]);
		}
	}

public:
	virtual ~ControlView() {}
	ControlView(): stepNum(0), step(0),  play(false), ShowMarker(true), ShowMesh(true), ShowLabel(false) {};
	static shared_ptr<ControlView> New()
	{
		shared_ptr<ControlView> nw = make_shared<ControlView>();
		return nw;
	}

	unsigned &getStepNum() { return stepNum; }
	unsigned & getStep() { return step; }
	bool & IsPlay() { return play; }
	bool & IsShowMarker() { return ShowMarker; }
	bool & IsShowMesh() { return ShowMesh; }
    bool & IsShowLabel() { return ShowLabel; }

	shared_ptr<Model> & getModel() { return pModel; }
	vtkSmartPointer<vtkProgrammableFilter> & getProgrammableFilter() { return programmableFilter; }
	vtkSmartPointer<vtkActor> &getActor() { return actor; }
	vtkSmartPointer<vtkTextActor> &getTextActor() {return textActor;}
	vtkSmartPointer<vtkActor> &getAxesActor() { return axesActor; }
    vtkSmartPointer<vtkActor2D> &getLabelActor() { return labelActor; }
    vtkSmartPointer<vtkSliderRepresentation2D> &getSliderRep(){return sliderRep;}
    vtkSmartPointer<vtkRenderer> &getRenderer(){return renderer;}
    vtkSmartPointer<vtkTexturedButtonRepresentation2D> &getButtonRepresentation(){return buttonRepresentation;}
    vtkSmartPointer<vtkRenderWindowInteractor> &getRenderWindowInteractor(){return renderWindowInteractor;}
    vtkSmartPointer<vtkCallbackCommand> &getTimerCallback(){return timerCallback;}

	virtual void inputModelfiles(vector<string> &modelFiles, vector<string>&dispFiles,
			const int& argc,  char* argv[]) 
	{

		parser(argc, argv, modelFiles, dispFiles);

		ugridReaders.resize(modelFiles.size());
		auto appendFilter = vtkSmartPointer<vtkAppendFilter>::New();
		for (unsigned i = 0; i< (unsigned)modelFiles.size(); ++i) {
			
			ugridReaders[i] = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
			ugridReaders[i]->SetFileName(modelFiles[i].c_str());
			ugridReaders[i]->Update();
			appendFilter->AddInputConnection(ugridReaders[i]->GetOutputPort());
		}
		appendFilter->Update();
		programmableFilter = vtkSmartPointer<vtkProgrammableFilter>::New();
		programmableFilter->AddInputConnection(appendFilter->GetOutputPort());
	}

	virtual void setContent(shared_ptr<Model> &model) {
		pModel = model;
		pModel->setVtkpnt0(programmableFilter->GetUnstructuredGridInput()->GetPoints());
	}

	virtual void Update(){

        if(buttonRepresentation->GetState())
            play = true;
        else
            play = false;

		stringstream ss("");
		ss << "Current Time = " << pModel->getStep(step);
		textActor->SetInput(ss.str().c_str());

		if (play) {
			step = (step == stepNum - 1) ? 0 : step + 1;
			//pModel->update(step);
			//sliderRep->SetTitleText(ss.str().c_str());
		}
		sliderRep->SetValue(step);
		programmableFilter->GetUnstructuredGridOutput()->SetPoints(pModel->getvtkPnts(step));	
		
	}

	virtual void setAnimationMethod( void(*f)(void*), vector<string> &dispFiles ) {
		programmableFilter->SetExecuteMethod(f, this);
		pModel->readDispfile(dispFiles);
		stepNum = pModel->getStepNum();

		pModel->initialize();
		renderWindowInteractor->CreateRepeatingTimer(10);
		timerCallback = vtkSmartPointer<vtkCallbackCommand>::New();
		timerCallback->SetCallback(TimerCallbackFunction);
		timerCallback->SetClientData(programmableFilter);
		renderWindowInteractor->AddObserver(vtkCommand::TimerEvent, timerCallback);

	}

	virtual void setKeyboardMethod(void (*f)(vtkObject* , long unsigned int vtkNotUsed(eventId), void* clientData, void* vtkNotUsed(callData))){
		keypressCallback = vtkSmartPointer<vtkCallbackCommand>::New();
		keypressCallback->SetCallback(f);
		keypressCallback->SetClientData(this);
		renderWindowInteractor->AddObserver(vtkCommand::KeyPressEvent, keypressCallback);
	}

	virtual void setWindowMethod(void (*f)(vtkObject* caller, long unsigned int vtkNotUsed(eventId), void* clientData, void* vtkNotUsed(callData))){
		windowCallback = vtkSmartPointer<vtkCallbackCommand>::New();
		windowCallback->SetCallback(f);
		windowCallback->SetClientData(this);
		renderWindow->AddObserver(vtkCommand::ModifiedEvent, windowCallback);
	}

	virtual void setMainActor() {
		// Create a mapper and actor
		mapper = vtkSmartPointer<vtkDataSetMapper>::New();
		mapper->SetInputConnection(programmableFilter->GetOutputPort());
		actor = vtkSmartPointer<vtkActor>::New();
		actor->SetMapper(mapper);
		actor->GetProperty()->SetColor(0.9, 0.9, 0.9);
		actor->GetProperty()->SetEdgeColor(0.0, 0.0, 0.0);
		actor->GetProperty()->EdgeVisibilityOn();
		//actor->SetScale(0.1);
	}

    virtual void setLabelActor(){
        labelMapper = vtkSmartPointer<vtkLabeledDataMapper>::New();

        for (unsigned i = 0; i< (unsigned)ugridReaders.size(); ++i) {
            labelMapper->AddInputConnection(ugridReaders[i]->GetOutputPort());
        }
        labelActor = vtkSmartPointer<vtkActor2D>::New();
        labelActor->SetMapper(labelMapper);

    }

	virtual void setRender() {
		// Create a renderer, render window, and interactor
		renderer = vtkSmartPointer<vtkRenderer>::New();
		renderer->GetActiveCamera()->SetClippingRange(0.01, 1000);
		renderer->SetAmbient(1.0, 1.0, 1.0);
		renderer->SetLightFollowCamera(1);

		renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
		renderWindow->AddRenderer(renderer);
		renderWindow->SetSize(800, 640);
		renderWindow->SetWindowName("MLV 2.0");
        renderWindow->Render();

		renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
		renderWindowInteractor->SetRenderWindow(renderWindow);
		style = vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
		renderWindowInteractor->SetInteractorStyle(style);

		// Initialize must be called prior to creating timer events.
		renderWindowInteractor->Initialize();

		axes = vtkSmartPointer<vtkAxesActor>::New();
		widget = vtkSmartPointer<vtkOrientationMarkerWidget>::New();
		widget->SetOutlineColor(0.9300, 0.5700, 0.1300);
		widget->SetOrientationMarker(axes);
		widget->SetInteractor(renderWindowInteractor);
		widget->SetViewport(0.0, 0.0, 0.2, 0.2);
		widget->SetEnabled(1);
		widget->InteractiveOn();
		widget->SetInteractive(0);

		int* windowSize = renderWindow->GetSize();
		textActor->SetPosition(windowSize[0] - 170, windowSize[1] - 25);
	}

	virtual void setAxesActor() {
		linesPolyData = vtkSmartPointer<vtkPolyData>::New();

		// Create three points
		double px0[3] = { -1e2, 0.0, 0.0 };
		double px1[3] = { 1e2, 0.0, 0.0 };

		double py0[3] = { 0.0, -1e2, 0.0 };
		double py1[3] = { 0.0, 1e2, 0.0 };

		double pz0[3] = { 0.0, 0.0, -1e2 };
		double pz1[3] = { 0.0, 0.0, 1e2 };

		// Create a vtkPoints container and store the points in it
		pts = vtkSmartPointer<vtkPoints>::New();
		pts->InsertNextPoint(px0);
		pts->InsertNextPoint(px1);

		pts->InsertNextPoint(py0);
		pts->InsertNextPoint(py1);

		pts->InsertNextPoint(pz0);
		pts->InsertNextPoint(pz1);
		// Add the points to the polydata container
		linesPolyData->SetPoints(pts);

		// Create the first line (between Origin and P0)
		line0 = vtkSmartPointer<vtkLine>::New();
		line0->GetPointIds()->SetId(0, 0); // the second 0 is the index of the Origin in linesPolyData's points
		line0->GetPointIds()->SetId(1, 1); // the second 1 is the index of P0 in linesPolyData's points

										   // Create the second line (between Origin and P1)
		line1 = vtkSmartPointer<vtkLine>::New();
		line1->GetPointIds()->SetId(0, 2); // the second 0 is the index of the Origin in linesPolyData's points
		line1->GetPointIds()->SetId(1, 3); // 2 is the index of P1 in linesPolyData's points

	    line2 = vtkSmartPointer<vtkLine>::New();
		line2->GetPointIds()->SetId(0, 4); // the second 0 is the index of the Origin in linesPolyData's points
		line2->GetPointIds()->SetId(1, 5); // 3 is the index of P1 in linesPolyData's points


										   // Create a vtkCellArray container and store the lines in it
		lines = vtkSmartPointer<vtkCellArray>::New();
		lines->InsertNextCell(line0);
		lines->InsertNextCell(line1);
		lines->InsertNextCell(line2);

		// Add the lines to the polydata container
		linesPolyData->SetLines(lines);

		// Create two colors - one for each line
		unsigned char red[3] = { 255, 0, 0 };
		unsigned char green[3] = { 0, 255, 0 };
		unsigned char blue[3] = { 0, 0, 255 };
		// Create a vtkUnsignedCharArray container and store the colors in it
		colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
		colors->SetNumberOfComponents(3);
		colors->InsertNextTupleValue(red);
		colors->InsertNextTupleValue(green);
		colors->InsertNextTupleValue(blue);
		// Color the lines.
		linesPolyData->GetCellData()->SetScalars(colors);

		// Setup the visualization pipeline
		axesMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		axesMapper->SetInputData(linesPolyData);
		axesActor = vtkSmartPointer<vtkActor>::New();
		axesActor->SetMapper(axesMapper);

	}

	virtual void setTextActor() {
		textActor = vtkSmartPointer<vtkTextActor>::New();
		//textActor->SetInput("Hello world");
		//textActor->SetPosition(80, 40);
		//textActor->GetTextProperty()->SetFontSize(24);
		//textActor->GetTextProperty()->SetColor(1.0, 0.0, 0.0);

		////////////////////

		textActor->GetTextProperty()->SetFontSize(16);
		textActor->GetTextProperty()->SetColor(1.0, 1.0, 1.0);

	}

    virtual void setSliderBar(){


        sliderRep = vtkSmartPointer<vtkSliderRepresentation2D>::New();

        sliderRep->SetMinimumValue(0);
        sliderRep->SetMaximumValue(stepNum-1);
        sliderRep->SetValue(0);

        //sliderRep->SetTitleText("");
        sliderRep->SetShowSliderLabel(0);

        // Set color properties:
        // Change the color of the knob that slides
        sliderRep->GetSliderProperty()->SetColor(1.0,1.0,1.0);//red
        sliderRep->GetSliderProperty()->SetOpacity(0.8);
        // Change the color of the bar
        sliderRep->GetTubeProperty()->SetColor(1.0,1.0,1.0);
        sliderRep->GetTubeProperty()->SetOpacity(0.2);
        // Change the color of the text indicating what the slider controls
        //sliderRep->GetTitleProperty()->SetColor(1,0,0);//red
//        sliderRep->GetTitleProperty()->SetBold(0);
//        sliderRep->GetTitleProperty()->SetShadow(0);
        // Change the color of the text displaying the value
        //sliderRep->GetLabelProperty()->SetColor(1,0,0);//red
        // Change the color of the knob when the mouse is held on it
        sliderRep->GetSelectedProperty()->SetColor(0.0,0.0,0.0);
        sliderRep->GetSelectedProperty()->SetOpacity(0.9);
        // Change the color of the ends of the bar
        //sliderRep->GetCapProperty()->SetColor(1,1,0);//yellow

        sliderRep->GetPoint1Coordinate()->SetCoordinateSystemToDisplay();
        sliderRep->GetPoint2Coordinate()->SetCoordinateSystemToDisplay();
        sliderRep->GetPoint1Coordinate()->SetValue(800 - 400, 640 - 15);
        sliderRep->GetPoint2Coordinate()->SetValue(800- 200, 640 - 15);

        sliderRep->SetEndCapLength(0);

        sliderWidget = vtkSmartPointer<vtkSliderWidget>::New();
        sliderWidget->SetInteractor(renderWindowInteractor);
        sliderWidget->SetRepresentation(sliderRep);
        sliderWidget->SetAnimationModeToAnimate();
        sliderWidget->EnabledOn();

        SliderCallback = vtkSmartPointer<vtkSliderCallback>::New();
        SliderCallback->pst = &step;
        sliderWidget->AddObserver(vtkCommand::InteractionEvent,SliderCallback);

        // button widget
        imagePlay = vtkSmartPointer<vtkImageData>::New();
        imagePause = vtkSmartPointer<vtkImageData>::New();
        CreateImagePlay(imagePlay);
        CreateImagePause(imagePause);
        buttonRepresentation = vtkSmartPointer<vtkTexturedButtonRepresentation2D>::New();
        buttonRepresentation->SetNumberOfStates(2);
        buttonRepresentation->SetButtonTexture(0,imagePlay);
        buttonRepresentation->SetButtonTexture(1,imagePause);

        buttonWidget = vtkSmartPointer<vtkButtonWidget>::New();
        buttonWidget->SetInteractor(renderWindowInteractor);
        buttonWidget->SetRepresentation(buttonRepresentation);

        vtkSmartPointer<vtkCoordinate> upperRight = vtkSmartPointer<vtkCoordinate>::New();
        upperRight->SetCoordinateSystemToNormalizedDisplay();
        //upperRight->SetCoordinateSystemToDisplay();
        upperRight->SetValue(1.0, 1.0);

        double bds[6];
        double sz = 20.0;

        bds[0] = upperRight->GetComputedDisplayValue(renderer)[0] - sz - 175;
        bds[1] = bds[0] + sz;
        bds[2] = upperRight->GetComputedDisplayValue(renderer)[1] - sz - 5;
        bds[3] = bds[2] + sz;
        bds[4] = bds[5] = 0.0;

        // Scale to 1, default is .5
        buttonRepresentation->SetPlaceFactor(1.0);
        buttonRepresentation->PlaceWidget(bds);

        buttonWidget->On();

    }


	virtual void Display() {
		// Add the actor to the scene
		renderer->AddActor(actor);
		renderer->AddActor(axesActor);
		renderer->AddActor2D(textActor);
        renderer->AddActor2D(labelActor);
        labelActor->VisibilityOff();

		renderer->GradientBackgroundOn();
		renderer->SetBackground2(13.0 / 255.0, 71.0 / 255.0, 161.0 / 255.0);
		renderer->SetBackground(144.0 / 255.0, 202.0 / 255.0, 249.0 / 255.0);
		// Render and interact

		renderWindowInteractor->Start();
	}

};


void TimerCallbackFunction(vtkObject* caller, long unsigned int vtkNotUsed(eventId), void* clientData, void* vtkNotUsed(callData))
{
	auto programmableFilter = static_cast<vtkProgrammableFilter*>(clientData);
	auto *iren = static_cast<vtkRenderWindowInteractor*>(caller);
	programmableFilter->Modified();
	iren->Render();
}


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
	pCtr->getTextActor()->SetPosition(windowSize[0] - 170, windowSize[1] - 25);

    if(pCtr->getStepNum()!=0){
        pCtr->getSliderRep()->GetPoint1Coordinate()->SetValue(windowSize[0] - 200 - windowSize[0]/4.0, windowSize[1] - 15);
        pCtr->getSliderRep()->GetPoint2Coordinate()->SetValue(windowSize[0] - 200, windowSize[1] - 15);

        vtkSmartPointer<vtkCoordinate> upperRight = vtkSmartPointer<vtkCoordinate>::New();
        upperRight->SetCoordinateSystemToNormalizedDisplay();
        upperRight->SetValue(1.0, 1.0);

        double bds[6];
        double sz = 20.0;

        bds[0] = upperRight->GetComputedDisplayValue(pCtr->getRenderer())[0] - sz - 175;
        bds[1] = bds[0] + sz;
        bds[2] = upperRight->GetComputedDisplayValue(pCtr->getRenderer())[1] - sz - 5;
        bds[3] = bds[2] + sz;
        bds[4] = bds[5] = 0.0;

        // Scale to 1, default is .5
        pCtr->getButtonRepresentation()->SetPlaceFactor(1.0);
        pCtr->getButtonRepresentation()->PlaceWidget(bds);
    }

}

void KeypressCallbackFunction(vtkObject* caller, long unsigned int vtkNotUsed(eventId), void* clientData, void* vtkNotUsed(callData))
{

	vtkRenderWindowInteractor *iren =
		static_cast<vtkRenderWindowInteractor*>(caller);

	ControlView* pCtr = static_cast<ControlView*>(clientData);

    if(pCtr->getStepNum()!=0){
	    if (strcmp(iren->GetKeySym(), "space") == 0){
	        pCtr->IsPlay() = !pCtr->IsPlay();
	        pCtr->getButtonRepresentation()->SetState(!pCtr->getButtonRepresentation()->GetState());
	//        if(pCtr->IsPlay()){
	//            pCtr->IsPlay() = !pCtr->IsPlay();
	//            pCtr->getRenderWindowInteractor()->InvokeEvent(vtkCommand::DisableEvent, pCtr->getTimerCallback());
	//        }
	//        else{
	//            pCtr->IsPlay() = !pCtr->IsPlay();
	//            pCtr->getRenderWindowInteractor()->InvokeEvent(vtkCommand::EnableEvent, pCtr->getTimerCallback());
	//        }

	    }


		if (strcmp(iren->GetKeySym(), "b") == 0) {
			pCtr->IsPlay() = false;
	        pCtr->getButtonRepresentation()->SetState(0);
			pCtr->getStep() = (pCtr->getStep() == pCtr->getStepNum() - 1) ? pCtr->getStepNum() - 1 : pCtr->getStep() + 1;
		}

		if (strcmp(iren->GetKeySym(), "v") == 0) {
			pCtr->IsPlay() = false;
	        pCtr->getButtonRepresentation()->SetState(0);
			pCtr->getStep() = (pCtr->getStep() == 0) ? 0 : pCtr->getStep() - 1;
		}
    }
	if (strcmp(iren->GetKeySym(), "i") == 0) {
		//vtkRendererCollection* rendercollction = iren->GetRenderWindow()->GetRenderers();
		if (pCtr->IsShowMarker())
			//rendercollction->GetFirstRenderer()->RemoveActor(pCtr->getAxesActor());
            pCtr->getAxesActor()->VisibilityOff();
		else
//			rendercollction->GetFirstRenderer()->AddActor(pCtr->getAxesActor());
            pCtr->getAxesActor()->VisibilityOn();
		pCtr->IsShowMarker() = !pCtr->IsShowMarker();
	}

	if (strcmp(iren->GetKeySym(), "u") == 0) {
		if (pCtr->IsShowMesh())
			pCtr->getActor()->GetProperty()->EdgeVisibilityOff();
		else
			pCtr->getActor()->GetProperty()->EdgeVisibilityOn();

		pCtr->IsShowMesh() = !pCtr->IsShowMesh();
	}

	if (strcmp(iren->GetKeySym(), "l") == 0) {
        if (pCtr->IsShowLabel())
            pCtr->getLabelActor()->VisibilityOff();
        else
            pCtr->getLabelActor()->VisibilityOn();

        pCtr->IsShowLabel() = !pCtr->IsShowLabel();
	}

	iren->Render();
}
#endif //CONTROLVIEW_H