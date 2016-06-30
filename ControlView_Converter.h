#pragma once
#pragma once
#ifndef CONTROLVIEW_CONVERTER_H
#define CONTROLVIEW_CONVERTER_H

#include "ControlView.h"
#include "vtkPLYReader.h"
#include "vtkPLYWriter.h"
#include "vtkXMLUnstructuredGridWriter.h"
class ControlView_Converter : public ControlView
{
protected:
	
public:
	virtual ~ControlView_Converter() {};
	ControlView_Converter() {};
	static shared_ptr<ControlView_Converter> New() {
		shared_ptr<ControlView_Converter> nw = make_shared<ControlView_Converter>();
		return nw;
	}

	///////////////////////////////////////////////////////////////////////////////////////////////
	// syntax
	// exec ./data/somefile 20 ugridfile
	///////////////////////////////////////////////////////////////////////////////////////////////
	void ConvertPLYSeries2UGridAndDisp(const int & argc, char *argv[]) 
	{
		const string prefix(argv[1]);    // 1st argument: file prefix including path and prefix
		const unsigned filenum(atoi(argv[2]));	// 2nd argument: PLY file num, 0, 1, 2, ... , num
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
		string outfile = argv[3];  // 3rd argument: the unstructuredgrid file name without postfix
		writer->SetFileName( (outfile + ".xml").c_str() );
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
		out.open( outfile + ".dat", ios::out | ios::binary);
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
		
		// display the first one
		programmableFilter = vtkSmartPointer<vtkProgrammableFilter>::New();
		programmableFilter->SetInputConnection(plyReaders[0]->GetOutputPort());

	}

	///////////////////////////////////////////////////////////////////////////////////////
	// syntax
	// exec /m ugridfile1.xml [ugridfile2.xml] /o plyfile.ply
	///////////////////////////////////////////////////////////////////////////////////////
	void ConvertUGrid2PLY(const int & argc, char *argv[]) 
	{
		vector<string> ugridfiles;
		vector<string> plyfiles;
		parser(argc, argv, ugridfiles, plyfiles);
		if (plyfiles.size() != 1 && plyfiles.size() != ugridfiles.size()) {
			printf("please input appropriate amount of the .ply file name.");
			system("pause");
			exit(0);
		}

		vector<vtkSmartPointer<vtkXMLUnstructuredGridReader> > ugridReader(ugridfiles.size());
		for (unsigned i = 0; i < ugridfiles.size(); ++i) {
			ugridReader[i] = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
			ugridReader[i]->SetFileName(ugridfiles[i].c_str());
			ugridReader[i]->Update();
		}

		if (plyfiles.size() == 1) {
			vtkSmartPointer<vtkGeometryFilter> geometryFilter = vtkSmartPointer<vtkGeometryFilter>::New();
			for (unsigned i = 0; i < ugridfiles.size(); ++i) {
				geometryFilter->AddInputConnection(ugridReader[i]->GetOutputPort());
			}
			vtkSmartPointer<vtkPLYWriter> plyWriter = vtkSmartPointer<vtkPLYWriter>::New();
			plyWriter->SetFileName(plyfiles[0].c_str());
			plyWriter->SetInputConnection(geometryFilter->GetOutputPort());
			plyWriter->Update();
		}
		else {
			vector<vtkSmartPointer<vtkPLYWriter>> plyWriter(plyfiles.size());
			vector<vtkSmartPointer<vtkGeometryFilter>> geometryFilter(plyfiles.size());
			for (unsigned i = 0; i < plyfiles.size(); ++i) {
				plyWriter[i] = vtkSmartPointer<vtkPLYWriter>::New();
				geometryFilter[i] = vtkSmartPointer<vtkGeometryFilter>::New();
				geometryFilter[i]->SetInputData(ugridReader[i]->GetOutput());

				plyWriter[i]->SetFileName(plyfiles[i].c_str());
				plyWriter[i]->SetInputConnection(geometryFilter[i]->GetOutputPort());
				plyWriter[i]->Update();
			}

		}

		cout << "writting ply file complete!" << endl;

		// display it
		programmableFilter = vtkSmartPointer<vtkProgrammableFilter>::New();
		programmableFilter->SetInputConnection(ugridReader[0]->GetOutputPort());

	}

	///////////////////////////////////////////////////////////////////////////////////////
	// syntax
	// exec /m ugridfile1.xml [ugridfile2.xml]/o ugridfile3.xml
	///////////////////////////////////////////////////////////////////////////////////////
	void ModifyUGrid(const int& argc, char *argv[])
	{
		vector<string> ugridfiles;
		vector<string> ugridfiles1;
		parser(argc, argv, ugridfiles, ugridfiles1);
		if (ugridfiles1.size() != 1 && ugridfiles1.size() != ugridfiles.size()) {
			printf("please input appropriate amount of the .xml file name.");
			system("pause");
			exit(0);
		}
		// vtkAppendFilter
		vtkSmartPointer<vtkAppendFilter> appendFilter = vtkSmartPointer<vtkAppendFilter>::New();
		vtkSmartPointer<vtkProgrammableFilter> programfilter = vtkSmartPointer<vtkProgrammableFilter>::New();
		vector<vtkSmartPointer<vtkXMLUnstructuredGridReader> > ugridReader(ugridfiles.size());
		for (unsigned i = 0; i < ugridfiles.size(); ++i) {
			ugridReader[i] = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
			ugridReader[i]->SetFileName(ugridfiles[i].c_str());
			ugridReader[i]->Update();
			appendFilter->AddInputConnection(ugridReader[i]->GetOutputPort());
		}
		appendFilter->Update();
		//vtkSmartPointer<vtkProgrammableFilter> programfilter = vtkSmartPointer<vtkProgrammableFilter>::New();
		//programfilter->SetInputConnection(ugridReader[0]->GetOutputPort());

		// the second UnstructuredGird

		programfilter->SetInputData(appendFilter->GetOutput(0));
		vtkPoints *point = programfilter->GetUnstructuredGridInput()->GetPoints();

		MatrixXd node(3, point->GetNumberOfPoints()/2);
		for (unsigned i = 0; i < point->GetNumberOfPoints()/2; ++i) {
			double *pnode = point->GetPoint(i);
			node.col(i) << pnode[0], pnode[1], pnode[2];
		}
		Matrix3d rotationZ;
		double thetaZ = M_PI / 2.0;
		
		rotationZ = (MatrixXd(3, 3) << cos(thetaZ), -sin(thetaZ), 0.0, sin(thetaZ), cos(thetaZ), 0.0, 0.0, 0.0, 1.0).finished();
		MatrixXd node1 = rotationZ* node;

		for (unsigned i = 0; i < point->GetNumberOfPoints()/2; ++i) {
			double *pnode = point->GetPoint(i);
			pnode[0] = node1(0, i) + 0.3;
			pnode[1] = node1(1, i) - 0.2;
			pnode[2] = node1(2, i) + 0.05;
			point->SetPoint(i, pnode);
		}
		//
		programfilter->GetUnstructuredGridOutput()->SetPoints(point);
		vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
		writer->SetFileName(ugridfiles1[0].c_str());
		writer->SetInputData(programfilter->GetUnstructuredGridInput());
		writer->Update();
		cout << "writting modified ugrid file complete!" << endl;

		// display it
		programmableFilter = vtkSmartPointer<vtkProgrammableFilter>::New();
		programmableFilter = programfilter;
		//programmableFilter->SetInputConnection(ugridReader[0]->GetOutputPort());
	}

	///////////////////////////////////////////////////////////////////////////////////////
	// syntax
	// exec /m plyfile1.ply [plyfile2.ply] /o ugridfile.xml
	///////////////////////////////////////////////////////////////////////////////////////
	void ConvertPLY2UGrid(const int & argc, char *argv[]) 
	{

		vector<string> plyfiles;
		vector<string> ugridfiles;
		parser(argc, argv, plyfiles, ugridfiles);
		if (ugridfiles.size() != 1 && ugridfiles.size() != plyfiles.size()) {
			printf("please input appropriate amount of the .ply file name.");
			system("pause");
			exit(0);
		}

		vector<vtkSmartPointer<vtkPLYReader>> plyReaders(plyfiles.size());
		for (unsigned i = 0; i < plyfiles.size(); ++i) {
			plyReaders[i] = vtkSmartPointer<vtkPLYReader>::New();
			plyReaders[i]->SetFileName(plyfiles[i].c_str());
			plyReaders[i]->Update();
		}

		if (ugridfiles.size() == 1) {
			vtkSmartPointer<vtkAppendFilter> appendFilter = vtkSmartPointer<vtkAppendFilter>::New();
			for (unsigned i = 0; i < plyReaders.size(); ++i) {
				appendFilter->AddInputData(plyReaders[i]->GetOutput());
			}
			appendFilter->Update();
			vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
			unstructuredGrid->ShallowCopy(appendFilter->GetOutput());
			vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
			writer->SetFileName(ugridfiles[0].c_str());
			writer->SetInputData(unstructuredGrid);
			writer->Update();
		}
		else {
			vector<vtkSmartPointer<vtkAppendFilter>> appendFilter(ugridfiles.size());
			vector<vtkSmartPointer<vtkUnstructuredGrid>> unstructuredGrid(ugridfiles.size());
			vector<vtkSmartPointer<vtkXMLUnstructuredGridWriter>> writer(ugridfiles.size());
			for (unsigned i = 0; i < plyReaders.size(); ++i) {
				appendFilter[i] = vtkSmartPointer<vtkAppendFilter>::New();
				appendFilter[i]->SetInputData(plyReaders[i]->GetOutput());
				appendFilter[i]->Update();

				unstructuredGrid[i] = vtkSmartPointer<vtkUnstructuredGrid>::New();
				unstructuredGrid[i]->ShallowCopy(appendFilter[i]->GetOutput());
				writer[i] = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
				writer[i]->SetFileName(ugridfiles[i].c_str());
				writer[i]->SetInputData(unstructuredGrid[i]);
				writer[i]->Update();
			}

		}
		
		cout << "writting model file complete!" << endl;
		// display it
		programmableFilter = vtkSmartPointer<vtkProgrammableFilter>::New();
		programmableFilter->SetInputConnection(plyReaders[0]->GetOutputPort());
	}


};

#endif //CONTROLVIEW_PLYFILESERIES_H