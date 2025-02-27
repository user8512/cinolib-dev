#include <cinolib/subdivision_barycentric.h>
#include <unordered_map>
#include <iostream>

namespace cinolib
{

template<class M, class V, class E, class F, class P>
CINO_INLINE
void subdivision_cc(Hexmesh<M,V,E,F,P> & m)
{
    std::unordered_map<uint,uint> p_map; // 存储体点
    std::unordered_map<uint,uint> f_map; // 存储面点
    std::unordered_map<uint,uint> e_map; // 存储边点
    std::unordered_map<uint,uint> v_map; // 存储点点

    std::unordered_map<uint,vec3d> p_centroid; // 存储体中心
    std::unordered_map<uint,vec3d> f_centroid; // 存储面中心
    std::unordered_map<uint,vec3d> e_centroid; // 存储边中点

    // m.adj_f2p(fid).size()!=2:边界面    m.vert_is_singular(vid):边界点    m.edge_is_singular(eid):边界边 
    // m.poly_centroid(pid):体的中心    m.face_centroid(fid):面的中心    m.edge_sample_at(eid,0.5):边的中点
    // 注意cinolib用于表示坐标的vec3d只重载了右乘标量(vec3d*scalar)

    /* 
    std::cout << "Vertexs " << m.num_verts() << std::endl;
    std::cout << "Edges " << m.num_edges() << std::endl;
    std::cout << "Faces " << m.num_faces() << std::endl;
    std::cout << "Polys " << m.num_polys() << std::endl;
    */

    // 提前计算各种中心点
    for (uint pid=0; pid<m.num_polys(); ++pid) 
      p_centroid[pid] = m.poly_centroid(pid);
    for (uint fid=0; fid<m.num_faces(); ++fid) 
      f_centroid[fid] = m.face_centroid(fid);
    for (uint eid=0; eid<m.num_edges(); ++eid) 
      e_centroid[eid] = m.edge_sample_at(eid,0.5);

    // 由于新的点直接添加到原网格中，用origin_verts_num来记录原有的顶点数量，可通过v_map查找对应的点点
    uint origin_verts_num = m.num_verts();

    // 体点计算规则:
    // Regular: 八个顶点的平均，内置方法完成
    // Irregular: 拓展到任意多面体时TODO
    for (uint pid=0; pid<m.num_polys(); ++pid) 
      p_map[pid] = m.vert_add(p_centroid.at(pid));

    // 面点计算规则:
    // Regular: 相邻的两个体点和面点 (C_0 + C_1 + 2 * F)/4
    // Irregular: 四个顶点的平均值
    for (uint fid=0; fid<m.num_faces(); ++fid)
    {
      std::vector<uint> adj_polys = m.adj_f2p(fid);
      if (adj_polys.size()==2)
        f_map[fid] = m.vert_add((p_centroid.at(adj_polys[0]) + p_centroid.at(adj_polys[1]) + f_centroid.at(fid)*2) / 4);
      else 
        f_map[fid] = m.vert_add(f_centroid.at(fid));
    }

    // 边点计算规则:
    // Regular:
    //      所有包含边的体点S_i的平均值记为S 包含边的面点F_i的平均值记为F 边的中点记为E
    //      N 是边相邻的面数
    //      (S + 2*F + (N-3)E)/N
    // Irregular:
    //      边所相邻的边界面的面点F_i的平均值记为F 边的中点记为E
    //      (F + E)/2
    for (uint eid=0; eid<m.num_edges(); ++eid) 
    {
      std::vector<uint> adj_faces = m.adj_e2f(eid);
      std::vector<uint> adj_polys = m.adj_e2p(eid);

      vec3d e_midpoint =e_centroid.at(eid);

      vec3d fp_avg = vec3d(0,0,0);
      for (uint adj_f : adj_faces) 
        fp_avg += m.vert(f_map.at(adj_f));
      fp_avg /= adj_faces.size();

      if (!m.edge_is_singular(eid))
      {
        vec3d pp_avg = vec3d(0,0,0);
        for (uint adj_p : adj_polys) 
          pp_avg += m.vert(p_map.at(adj_p));
        pp_avg /= adj_polys.size();

        int N = adj_faces.size();
        e_map[eid] = m.vert_add((pp_avg + fp_avg*2 + e_midpoint*(N-3)) / N);
      }
      else 
        e_map[eid] = m.vert_add((fp_avg + e_midpoint) / 2);
    }

    // 点点计算规则
    //  Regular:
    //      所有包含点的体点S_i的平均值记为S  包含点的面点F_i的平均值记为F  包含点的边点E_i的平均值记为E
    //      (S + 3*F + 3*E + V)/8
    //  Irregular:
    //      所有包含点的边界面面点F_i的平均值记为F  包含点的边界边点E_i的平均值记为E  点记为V
    //      N 为V相邻的边界边数
    //      (F + 2*E + (N-3)V)/N
    for (uint vid=0; vid<origin_verts_num; ++vid) 
    {
      std::vector<uint> adj_edges = m.adj_v2e(vid);
      std::vector<uint> adj_faces = m.adj_v2f(vid);
      std::vector<uint> adj_polys = m.adj_v2p(vid);

      vec3d v_pos = m.vert(vid);

      vec3d ep_avg = vec3d(0,0,0);
      for (uint adj_e : adj_edges)
        ep_avg += m.vert(e_map.at(adj_e));
      ep_avg /= adj_edges.size();

      vec3d fp_avg = vec3d(0,0,0);
      for (uint adj_f : adj_faces)
        fp_avg += m.vert(f_map.at(adj_f));
      fp_avg /= adj_faces.size();

      if (!m.vert_is_singular(vid))
      {
        vec3d pp_avg = vec3d(0,0,0);
        for (uint adj_p : adj_polys)
          pp_avg += m.vert(p_map.at(adj_p));
        pp_avg /= adj_polys.size();

        v_map[vid] = m.vert_add((pp_avg + fp_avg*3 + ep_avg*3 + v_pos) / 8);
      }
      else
      {
        int N = adj_edges.size();
        v_map[vid] = m.vert_add((fp_avg + ep_avg*2 + v_pos*(N-3)) / N);
      }
    }
    uint np = m.num_polys();

    // 依照原网格的多面体，逐个构建新的多面体

    // cinolib的六面体单元序号如下(点, 边), 面的顺序是下-右-上-左-前-后
    //        7-------6        +---6---+       
    //       /|      /|       7|      5|
    //      / |     / |      / 11    / 10
    //     4-------5  |     +----4--+  |
    //     |  3----|--2     |  +--2-|--+
    //     | /     | /      8 3     9 1
    //     |/      |/       |/      |/
    //     0-------1        +---0---+

    for (uint pid=0; pid<np; ++pid)
    {
        // hex verts id
        uint v_id[8] =
        {
            m.poly_vert_id(pid,0),
            m.poly_vert_id(pid,1),
            m.poly_vert_id(pid,2),
            m.poly_vert_id(pid,3),
            m.poly_vert_id(pid,4),
            m.poly_vert_id(pid,5),
            m.poly_vert_id(pid,6),
            m.poly_vert_id(pid,7)
        };

        // verts point
        uint v[8] =
        {
          v_map.at(v_id[0]),
          v_map.at(v_id[1]),
          v_map.at(v_id[2]),
          v_map.at(v_id[3]),
          v_map.at(v_id[4]),
          v_map.at(v_id[5]),
          v_map.at(v_id[6]),
          v_map.at(v_id[7])
        };

        // edges point
        // 对应关系存于standard_elements_tables.h的HEXA_EDGES中
        uint e[12] =
        {
          e_map.at(m.edge_id(v_id[0], v_id[1])),
          e_map.at(m.edge_id(v_id[1], v_id[2])),
          e_map.at(m.edge_id(v_id[2], v_id[3])),
          e_map.at(m.edge_id(v_id[3], v_id[0])),
          e_map.at(m.edge_id(v_id[4], v_id[5])),
          e_map.at(m.edge_id(v_id[5], v_id[6])),
          e_map.at(m.edge_id(v_id[6], v_id[7])),
          e_map.at(m.edge_id(v_id[7], v_id[4])),
          e_map.at(m.edge_id(v_id[0], v_id[4])),
          e_map.at(m.edge_id(v_id[1], v_id[5])),
          e_map.at(m.edge_id(v_id[2], v_id[6])),
          e_map.at(m.edge_id(v_id[3], v_id[7]))
        };

        // faces point
        // 对应关系存于standard_elements_tables.h的HEXA_FACES中
        uint f[6] =
        {
          f_map.at(m.face_id({v_id[0], v_id[3], v_id[2], v_id[1]})),
          f_map.at(m.face_id({v_id[1], v_id[2], v_id[6], v_id[5]})),
          f_map.at(m.face_id({v_id[4], v_id[5], v_id[6], v_id[7]})),
          f_map.at(m.face_id({v_id[3], v_id[0], v_id[4], v_id[7]})),
          f_map.at(m.face_id({v_id[0], v_id[1], v_id[5], v_id[4]})),
          f_map.at(m.face_id({v_id[2], v_id[3], v_id[7], v_id[6]}))
        };

        // poly point
        uint p = p_map.at(pid);

        // 构建细分后的六面体拓扑结构
        m.poly_add({v[0], e[0], f[0], e[3], e[8], f[4], p, f[3]});
        m.poly_add({e[0], v[1], e[1], f[0], f[4], e[9], f[1], p});
        m.poly_add({e[3], f[0], e[2], v[3], f[3], p, f[5], e[11]});
        m.poly_add({f[0], e[1], v[2], e[2], p, f[1], e[10], f[5]});
        m.poly_add({e[8], f[4], p, f[3], v[4], e[4], f[2], e[7]});
        m.poly_add({f[4], e[9], f[1], p, e[4], v[5], e[5], f[2]});
        m.poly_add({f[3], p, f[5], e[11], e[7], f[2], e[6], v[7]});
        m.poly_add({p, f[1], e[10], f[5], f[2], e[5], v[6], e[6]});
    }

    // remove the old polys
    for(int pid=np-1; pid>=0; --pid) 
      m.poly_remove(pid);
}
}
