#pragma once
#pragma once
#ifndef CONTROLVIEW_PLYFILESERIES_H
#define CONTROLVIEW_PLYFILESERIES_H

#include"ControlView.h"
#include "vtkPLYReader.h"
#include "vtkPLYWriter.h"
#include "vtkXMLUnstructuredGridWriter.h"
class ControlView_PLYfileSeries : public ControlView
{
protected:
	
public:
	virtual ~ControlView_PLYfileSeries() {};
	ControlView_PLYfileSeries() {};
	static shared_ptr<ControlView_PLYfileSeries> New() {
		shared_ptr<ControlView_PLYfileSeries> nw = make_shared<ControlView_PLYfileSeries>();
		return nw;
	}

	void ConvertPLY2Unstrgrid(const int & argc, char *argv[]) {

		const string prefix(argv[1]);
		const unsigned filenum(atoi(argv[2]));
		vector<vtkSmartPointer<vtkPLYReader>> plyReaders;
		plyReaders.resize(filenum);
		stringstream ss("");
		for (unsigned i = 0; i< filenum; ++i) {
			ss << prefix << i << ".ply";
			plyReaders[i] = vtkSmartPointer<vtkPLYReader>::New();
			plyReaders[i]->SetFileName(ss.str().c_str());
			plyReaders[i]->Update();
			
			cout << "reading: " << ss.str() << endl;
			ss.str("");
		}

		vtkSmartPointer<vtkAppendFilter> appendFilter =
			vtkSmartPointer<vtkAppendFilter>::New();
		appendFilter->AddInputData(plyReaders[0]->GetOutput());
		appendFilter->Update();
		vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid =
			vtkSmartPointer<vtkUnstructuredGrid>::New();
		unstructuredGrid->ShallowCopy(appendFilter->GetOutput());

		vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
		writer->SetFileName("mdel_clothball.xml");
		writer->SetInputData(unstructuredGrid);
		writer->Update();

		cout << "writting model file complete!" << endl;

		unsigned nodeNum = plyReaders[0]->GetOutput()->GetNumberOfPoints();
		VectorXd node0(nodeNum * 3);
		node0.setZero();
		for (unsigned i = 0; i < nodeNum; ++i) {
			double *xyz = plyReaders[0]->GetOutput()->GetPoint(i);
			node0.segment<3>(3 * i) = Vector3d(*xyz, *(xyz + 1), *(xyz + 2));
		}

		ofstream out;
		out.open("out_clothball.dat", ios::out | ios::binary);
		VectorXd node(3 * nodeNum), disp(3 * nodeNum);

		for (unsigned s = 0; s < plyReaders.size(); ++s) {
			node.setZero();
			for (unsigned i = 0; i < nodeNum; ++i) {
				double *xyz = plyReaders[s]->GetOutput()->GetPoint(i);
				node.segment<3>(3 * i) = Vector3d(*xyz, *(xyz + 1), *(xyz + 2));
			}
			disp = node - node0;
			double ds = (double)s;
			out.write((char*)&ds, sizeof(double));
			out.write((char*)disp.data(), sizeof(double) * 3 * nodeNum);

			cout << "write: " << s << " step" << endl;
		}
		out.close();
		cout << "writting disp file complete!" << endl;
		exit(0);
	}

	void ConvertUnstrgrid2PLY(const int & argc, char *argv[]) {
		
		vtkSmartPointer<vtkXMLUnstructuredGridReader> ugridReader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
		ugridReader->SetFileName(argv[2]);
		ugridReader->Update();

		vtkSmartPointer<vtkGeometryFilter> geometryFilter = vtkSmartPointer<vtkGeometryFilter>::New();
		geometryFilter->SetInputData(ugridReader->GetOutput());

		vtkSmartPointer<vtkPLYWriter> plyWriter = vtkSmartPointer<vtkPLYWriter>::New();
		plyWriter->SetFileName("mdel_longplate40.ply");
		plyWriter->SetInputConnection(geometryFilter->GetOutputPort());
		plyWriter->Update();

		cout << "writting ply file complete!" << endl;
 		exit(0);

	}
};





#endif //CONTROLVIEW_PLYFILESERIES_H