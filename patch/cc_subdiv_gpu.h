#pragma once
#include <vector>
#include <cstdint>
#include <cinolib/geometry/vec_mat.h>

namespace cinolib {

    void subdiv_cuda(
        const std::vector<vec3d>& vertsPos,
        int patchVerts, int patchEdges, int patchFaces, int patchPolys,

        const std::vector<vec2i>& edgeVerts,   // size = ne
        const std::vector<int>& faceEdges,    // size = 4*nf, 允许负号
        const std::vector<int>& polyFaces,    // size = 6*np, 允许负号

        const std::vector<uint8_t>& faceOnSurf,   // size = patchFaces (0/1)
        const std::vector<uint8_t>& edgeOnSurf,   // size = patchEdges (0/1)
        const std::vector<uint8_t>& vertOnSurf,   // size = patchVerts (0/1)

        const std::vector<int>& facePolysOffset, // size = nf+1
        const std::vector<int>& facePolys,       // 允许负号

        const std::vector<int>& edgeFacesOffset, // size = ne+1
        const std::vector<int>& edgeFaces,       // 允许负号

        const std::vector<int>& vertEdgesOffset, // size = nv+1
        const std::vector<int>& vertEdges,       // 允许负号

        std::vector<vec3d>& newPolyVerts,    // size = patchPolys
        std::vector<vec3d>& newFaceVerts,    // size = patchFaces
        std::vector<vec3d>& newEdgeVerts,    // size = patchEdges
        std::vector<vec3d>& newVertVerts     // size = patchVerts
    );

}
