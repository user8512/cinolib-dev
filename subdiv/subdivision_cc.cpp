#include <cinolib/subdivision_barycentric.h>
#include <unordered_map>
#include <iostream>

namespace cinolib
{

template<class M, class V, class E, class F, class P>
CINO_INLINE
void subdivision_cc(Hexmesh<M,V,E,F,P> & m)
{
    std::unordered_map<uint,uint> e_map; // edge midpoints
    std::unordered_map<uint,uint> f_map; // face centroids
    std::unordered_map<uint,uint> p_map; // poly centroids
    
    for(uint eid=0; eid<m.num_edges(); ++eid) e_map[eid] = m.vert_add(m.edge_sample_at(eid,0.5));
    for(uint fid=0; fid<m.num_faces(); ++fid) f_map[fid] = m.vert_add(m.face_centroid(fid));
    for(uint pid=0; pid<m.num_polys(); ++pid) p_map[pid] = m.vert_add(m.poly_centroid(pid));
    uint np = m.num_polys();
    for(uint pid=0; pid<np; ++pid)
    {
        // hex verts
        uint v[8] =
        {
            m.poly_vert_id(pid,0),
            m.poly_vert_id(pid,1),
            m.poly_vert_id(pid,2),
            m.poly_vert_id(pid,3),
            m.poly_vert_id(pid,4),
            m.poly_vert_id(pid,5),
            m.poly_vert_id(pid,6),
            m.poly_vert_id(pid,7),
        };

        // hex faces
        uint f[6][4] =
        {
            { v[HEXA_FACES[0][0]], v[HEXA_FACES[0][1]], v[HEXA_FACES[0][2]], v[HEXA_FACES[0][3]] },
            { v[HEXA_FACES[1][0]], v[HEXA_FACES[1][1]], v[HEXA_FACES[1][2]], v[HEXA_FACES[1][3]] },
            { v[HEXA_FACES[2][0]], v[HEXA_FACES[2][1]], v[HEXA_FACES[2][2]], v[HEXA_FACES[2][3]] },
            { v[HEXA_FACES[3][0]], v[HEXA_FACES[3][1]], v[HEXA_FACES[3][2]], v[HEXA_FACES[3][3]] },
            { v[HEXA_FACES[4][0]], v[HEXA_FACES[4][1]], v[HEXA_FACES[4][2]], v[HEXA_FACES[4][3]] },
            { v[HEXA_FACES[5][0]], v[HEXA_FACES[5][1]], v[HEXA_FACES[5][2]], v[HEXA_FACES[5][3]] },
        };

        // hex edges (sorted per face)
        uint e[6][4] =
        {
            { e_map.at(m.edge_id(f[0][0], f[0][1])),
              e_map.at(m.edge_id(f[0][1], f[0][2])),
              e_map.at(m.edge_id(f[0][2], f[0][3])),
              e_map.at(m.edge_id(f[0][3], f[0][0])) },

            { e_map.at(m.edge_id(f[1][0], f[1][1])),
              e_map.at(m.edge_id(f[1][1], f[1][2])),
              e_map.at(m.edge_id(f[1][2], f[1][3])),
              e_map.at(m.edge_id(f[1][3], f[1][0])) },

            { e_map.at(m.edge_id(f[2][0], f[2][1])),
              e_map.at(m.edge_id(f[2][1], f[2][2])),
              e_map.at(m.edge_id(f[2][2], f[2][3])),
              e_map.at(m.edge_id(f[2][3], f[2][0])) },
            
            { e_map.at(m.edge_id(f[3][0], f[3][1])),
              e_map.at(m.edge_id(f[3][1], f[3][2])),
              e_map.at(m.edge_id(f[3][2], f[3][3])),
              e_map.at(m.edge_id(f[3][3], f[3][0])) },

            { e_map.at(m.edge_id(f[4][0], f[4][1])),
              e_map.at(m.edge_id(f[4][1], f[4][2])),
              e_map.at(m.edge_id(f[4][2], f[4][3])),
              e_map.at(m.edge_id(f[4][3], f[4][0])) },
            
            { e_map.at(m.edge_id(f[5][0], f[5][1])),
              e_map.at(m.edge_id(f[5][1], f[5][2])),
              e_map.at(m.edge_id(f[5][2], f[5][3])),
              e_map.at(m.edge_id(f[5][3], f[5][0])) }
        };

        // face centroids
        uint fc[6] =
        {
            f_map.at(m.face_id({f[0][0], f[0][1], f[0][2], f[0][3]})),
            f_map.at(m.face_id({f[1][0], f[1][1], f[1][2], f[1][3]})),
            f_map.at(m.face_id({f[2][0], f[2][1], f[2][2], f[2][3]})),
            f_map.at(m.face_id({f[3][0], f[3][1], f[3][2], f[3][3]})),
            f_map.at(m.face_id({f[4][0], f[4][1], f[4][2], f[4][3]})),
            f_map.at(m.face_id({f[5][0], f[5][1], f[5][2], f[5][3]}))
        };

        // hex centroid
        uint c = p_map.at(pid);
        // split i^th face
        m.poly_add({e[0][4], f[0][3], e[0][3], fc[0], fc[4], e[4][1], fc[1], c});
        m.poly_add({f[0][0], e[0][4], fc[0], e[0][1], e[4][3], fc[4], c, fc[3]});
        m.poly_add({fc[0], e[0][3], f[0][2], e[0][2], c, fc[1], e[5][3], fc[5]});
        m.poly_add({e[0][1], fc[0], e[0][2], f[0][1], fc[3], c, fc[5], e[5][1]});
        m.poly_add({fc[4], e[4][1], fc[1], c, e[2][0], f[3][1], e[2][1], fc[2]});
        m.poly_add({e[4][3], fc[4], c, fc[3], f[3][0], e[2][0], fc[2], e[2][3]});
        m.poly_add({c, fc[1], e[5][3], fc[5], fc[2], e[2][1], f[3][2], e[2][2]});
        m.poly_add({fc[3], c, fc[5], e[5][1], e[2][3], fc[2], e[2][2], f[3][3]});
    }

    // remove the old polys
    for(int pid=np-1; pid>=0; --pid) m.poly_remove(pid);
}


}
