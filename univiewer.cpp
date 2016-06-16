#include "univiewer.h"

int main(int argc, char *argv[])
{
    shared_ptr<Univiewer> uv = Univiewer::New();
    
    MatrixXu elem(3, 2);
    elem.col(0) << 0, 1, 2;
    elem.col(1) << 1, 3, 2;
    
    MatrixXd node(3, 4);
    node.col(0) = Vector3d(0.0, 0.0, 0.0);
    node.col(1) = Vector3d(1.0, 0.0, 0.0);
    node.col(2) = Vector3d(0.0, 1.0, 0.0);
    node.col(3) = Vector3d(1.0, 1.0, 0.0);
    
    uv->plotModel(elem, node);
    cout<< "xxxx" <<endl;
    return 0;
}