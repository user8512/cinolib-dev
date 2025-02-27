#include <cinolib/meshes/meshes.h>
#include "subdivision_cc.h"
#include <cinolib/gl/glcanvas.h>
#include <cinolib/gl/volume_mesh_controls.h>


int main(int argc, char **argv)
{
    using namespace cinolib;
    std::string s = (argc==2) ? std::string(argv[1]) : "D:/Coding/cinolib/examples/data/rockerarm.mesh";
    Hexmesh<> m(s.c_str());
    subdivision_cc(m);
    m.save("D:/Coding/cinolib/subdiv/reckerarm_subdivded.mesh");
    DrawableHexmesh<> dm("D:/Coding/cinolib/subdiv/reckerarm_subdivded.mesh");
    GLcanvas gui;
    VolumeMeshControls<DrawableHexmesh<>> menu(&dm, &gui);
    gui.push(&dm);
    gui.push(&menu);
    return gui.launch();
}
