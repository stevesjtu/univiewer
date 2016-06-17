#include "univiewer.h"

int main(int argc, char *argv[])
{
    shared_ptr<Univiewer> uv = Univiewer::New();

	//uv->plotModel(argc, argv);
    MatrixXu elem(3, 2);
    elem.col(0) << 0, 1, 2;
    elem.col(1) << 1, 3, 2;
    
    MatrixXd node(3, 4);
    node.col(0) = Vector3d(0.0, 0.0, 0.0);
    node.col(1) = Vector3d(1.0, 0.0, 0.0);
    node.col(2) = Vector3d(0.0, 1.0, 0.0);
    node.col(3) = Vector3d(1.0, 1.0, 0.0);
    
	vector<TriangleMesh> meshes;
	MatrixXd node1 = node;
	node1.row(0).array() += 1.0;

	meshes.push_back(TriangleMesh(elem, node));
	meshes.push_back(TriangleMesh(elem, node1));

    uv->plotModel(meshes);
    cout<< "xxxx" <<endl;
    return 0;
}