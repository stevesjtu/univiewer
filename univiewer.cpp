#include "univiewer.h"

int main(int argc, char *argv[])
{
    shared_ptr<Univiewer> uv = Univiewer::New();
    uv->plotModel(argc, argv);
    cout<< "xxxx" <<endl;
    return 0;
}