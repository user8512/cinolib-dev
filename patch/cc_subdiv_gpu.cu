#define CINO_STATIC_LIB
#include "patch.h"
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/transform.h>
#include <thrust/for_each.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/functional.h>
#include <cuda_runtime.h>
#include <vector>
#include <algorithm>
#include <cstdint>

namespace cinolib {
    // 轻量设备/主机两用的三维向量（用于 GPU 侧计算）
    struct dvec3 {
        double x, y, z;
        __host__ __device__ dvec3() : x(0), y(0), z(0) {}
        __host__ __device__ dvec3(double X, double Y, double Z) : x(X), y(Y), z(Z) {}
        __host__ __device__ dvec3& operator+=(const dvec3& o) { x += o.x; y += o.y; z += o.z; return *this; }
        __host__ __device__ dvec3& operator-=(const dvec3& o) { x -= o.x; y -= o.y; z -= o.z; return *this; }
        __host__ __device__ dvec3& operator*=(double s) { x *= s; y *= s; z *= s; return *this; }
        __host__ __device__ dvec3& operator/=(double s) { x /= s; y /= s; z /= s; return *this; }
        __host__ __device__ friend dvec3 operator+(dvec3 a, const dvec3& b) { a += b; return a; }
        __host__ __device__ friend dvec3 operator-(dvec3 a, const dvec3& b) { a -= b; return a; }
        __host__ __device__ friend dvec3 operator*(dvec3 a, double s) { a *= s; return a; }
        __host__ __device__ friend dvec3 operator*(double s, dvec3 a) { a *= s; return a; }
        __host__ __device__ friend dvec3 operator/(dvec3 a, double s) { a /= s; return a; }
    };

    __host__ __device__ inline int idx_abs(int v) { return (v < 0) ? (-v - 1) : v; }

    static inline dvec3 to_d(const cinolib::vec3d& v);
    static inline vec3d to_h(const dvec3& v);

    struct EdgeCentroidOp {
        const uint2* EV;             // edgeVerts, 每条边两个端点的局部顶点索引
        const dvec3* Vpos;           // 旧顶点坐标（局部）
        dvec3* Ecentroids;     // 输出
        __host__ __device__
            EdgeCentroidOp(const uint2* ev, const dvec3* vp, dvec3* out)
            : EV(ev), Vpos(vp), Ecentroids(out) {}
        __host__ __device__
            void operator()(const int e) const {
            uint2 ev = EV[e];
            dvec3 a = Vpos[ev.x];
            dvec3 b = Vpos[ev.y];
            Ecentroids[e] = (a + b) * 0.5;
        }
    };

    // 2) 面质心：face 是四条边的平均（边心平均）
    struct FaceCentroidOp {
        const int* faceEdges; // 长度 nf*4
        const dvec3* Ecentroids;
        dvec3* Fcentroids;
        __host__ __device__
            FaceCentroidOp(const int* fe, const dvec3* ec, dvec3* out)
            : faceEdges(fe), Ecentroids(ec), Fcentroids(out) {}
        __host__ __device__
            void operator()(const int f) const {
            dvec3 c(0, 0, 0);
#pragma unroll
            for (int k = 0; k < 4; ++k) {
                int e = idx_abs(faceEdges[f * 4 + k]);
                c += Ecentroids[e];
            }
            Fcentroids[f] = c / 4.0;
        }
    };

    // 3) 体质心：poly 是六个面的平均（面心平均）
    struct PolyCentroidOp {
        const int* polyFaces; // 长度 np*6
        const dvec3* Fcentroids;
        dvec3* Pcentroids;
        __host__ __device__
            PolyCentroidOp(const int* pf, const dvec3* fc, dvec3* out)
            : polyFaces(pf), Fcentroids(fc), Pcentroids(out) {}
        __host__ __device__
            void operator()(const int p) const {
            dvec3 c(0, 0, 0);
#pragma unroll
            for (int k = 0; k < 6; ++k) {
                int f = idx_abs(polyFaces[p * 6 + k]);
                c += Fcentroids[f];
            }
            Pcentroids[p] = c / 6.0;
        }
    };

    // 4) 新面点
    struct NewFaceVertOp {
        const uint8_t* faceOnSurf;
        const uint* facePolysOffset;
        const int* facePolys;
        const dvec3* Fcentroids;
        const dvec3* Pcentroids;
        dvec3* outNewF;
        __host__ __device__
            NewFaceVertOp(const uint8_t* fos, const uint* fpoff, const int* fp,
                const dvec3* fc, const dvec3* pc, dvec3* out)
            : faceOnSurf(fos), facePolysOffset(fpoff), facePolys(fp),
            Fcentroids(fc), Pcentroids(pc), outNewF(out) {}
        __host__ __device__
            void operator()(const int f) const {
            if (!faceOnSurf[f]) {
                uint off = facePolysOffset[f];
                int p0 = idx_abs(facePolys[off + 0]);
                int p1 = idx_abs(facePolys[off + 1]);
                dvec3 v = Pcentroids[p0] + Pcentroids[p1] + (Fcentroids[f] * 2.0);
                outNewF[f] = v / 4.0;
            }
            else {
                outNewF[f] = Fcentroids[f];
            }
        }
    };

    // 5) 新边点（含边界/内部两种规则）
    // 说明：由于要做小规模的“邻接聚合 + 去重”，在每个线程用固定小数组。
    //      阈值给 32，若真实邻接超过则截断（极少见的异常网格）。
    struct NewEdgeVertOp {
        const uint8_t* edgeOnSurf;
        const uint8_t* faceOnSurf;
        const uint* edgeFacesOffset;
        const int* edgeFaces;
        const uint* facePolysOffset;
        const int* facePolys;

        const dvec3* Ecentroids;
        const dvec3* Fcentroids;
        const dvec3* Pcentroids;

        dvec3* outNewE;

        __host__ __device__
            NewEdgeVertOp(const uint8_t* eos, const uint8_t* fos,
                const uint* efoff, const int* ef,
                const uint* fpoff, const int* fp,
                const dvec3* ec, const dvec3* fc, const dvec3* pc,
                dvec3* out)
            : edgeOnSurf(eos), faceOnSurf(fos),
            edgeFacesOffset(efoff), edgeFaces(ef),
            facePolysOffset(fpoff), facePolys(fp),
            Ecentroids(ec), Fcentroids(fc), Pcentroids(pc),
            outNewE(out) {}

        __host__ __device__
            void operator()(const int e) const {
            uint begin = edgeFacesOffset[e];
            uint end = edgeFacesOffset[e + 1];
            int N = int(end - begin);

            if (!edgeOnSurf[e]) {
                // 收集相邻面
                int faces[32]; int nf = 0;
                for (int i = 0; i < N && nf < 32; ++i) {
                    int f = idx_abs(edgeFaces[begin + i]);
                    faces[nf++] = f;
                }
                // 从相邻面收集体并去重
                int polys[32]; int np = 0;
                for (int i = 0; i < nf; ++i) {
                    int f = faces[i];
                    uint fbegin = facePolysOffset[f];
                    uint fend = facePolysOffset[f + 1];
                    for (uint k = fbegin; k < fend; ++k) {
                        int p = idx_abs(facePolys[k]);
                        bool seen = false;
                        for (int t = 0; t < np; ++t) if (polys[t] == p) { seen = true; break; }
                        if (!seen && np < 32) polys[np++] = p;
                    }
                }
                // 平均
                dvec3 faceAvg(0, 0, 0);
                for (int i = 0; i < nf; ++i) faceAvg += Fcentroids[faces[i]];
                if (nf > 0) faceAvg /= double(nf);

                dvec3 polyAvg(0, 0, 0);
                for (int i = 0; i < np; ++i) polyAvg += Pcentroids[polys[i]];
                if (np > 0) polyAvg /= double(np);

                dvec3 v = polyAvg + (faceAvg * 2.0) + (Ecentroids[e] * double(N - 3));
                outNewE[e] = v / double(N);
            }
            else {
                if (N == 1) {
                    outNewE[e] = Ecentroids[e];
                }
                else {
                    int facesB[32]; int nf = 0;
                    for (int i = 0; i < N && nf < 32; ++i) {
                        int f = idx_abs(edgeFaces[begin + i]);
                        if (faceOnSurf[f]) facesB[nf++] = f;
                    }
                    dvec3 faceAvg(0, 0, 0);
                    for (int i = 0; i < nf; ++i) faceAvg += Fcentroids[facesB[i]];
                    if (nf > 0) faceAvg /= double(nf);
                    dvec3 v = faceAvg + Ecentroids[e];
                    outNewE[e] = v * 0.5;
                }
            }
        }
    };

    // 6) 新顶点
    struct NewVertVertOp {
        const uint8_t* vertOnSurf;
        const uint8_t* edgeOnSurf;
        const uint8_t* faceOnSurf;

        const uint* vertEdgesOffset;
        const int* vertEdges;

        const uint* edgeFacesOffset;
        const int* edgeFaces;

        const uint* facePolysOffset;
        const int* facePolys;

        const dvec3* Vpos;
        const dvec3* Ecentroids;
        const dvec3* Fcentroids;
        const dvec3* Pcentroids;

        dvec3* outNewV;

        __host__ __device__
            NewVertVertOp(const uint8_t* vos, const uint8_t* eos, const uint8_t* fos,
                const uint* veoff, const int* ve,
                const uint* efoff, const int* ef,
                const uint* fpoff, const int* fp,
                const dvec3* vp, const dvec3* ec, const dvec3* fc, const dvec3* pc,
                dvec3* out)
            : vertOnSurf(vos), edgeOnSurf(eos), faceOnSurf(fos),
            vertEdgesOffset(veoff), vertEdges(ve),
            edgeFacesOffset(efoff), edgeFaces(ef),
            facePolysOffset(fpoff), facePolys(fp),
            Vpos(vp), Ecentroids(ec), Fcentroids(fc), Pcentroids(pc),
            outNewV(out) {}

        __host__ __device__
            void operator()(const int v) const {
            uint begin = vertEdgesOffset[v];
            uint end = vertEdgesOffset[v + 1];
            int N = int(end - begin);

            if (!vertOnSurf[v]) {
                // 收集边/面/体（去重）
                int edges[32]; int ne = 0;
                int faces[32]; int nf = 0;
                int polys[32]; int np = 0;

                for (int i = 0; i < N && ne < 32; ++i) {
                    int e = idx_abs(vertEdges[begin + i]);
                    edges[ne++] = e;
                    // 经边找面
                    uint eb = edgeFacesOffset[e];
                    uint ee = edgeFacesOffset[e + 1];
                    for (uint j = eb; j < ee; ++j) {
                        int f = idx_abs(edgeFaces[j]);
                        bool seen = false;
                        for (int t = 0; t < nf; ++t) if (faces[t] == f) { seen = true; break; }
                        if (!seen && nf < 32) faces[nf++] = f;
                        // 经面找体
                        uint fb = facePolysOffset[f];
                        uint fe = facePolysOffset[f + 1];
                        for (uint k = fb; k < fe; ++k) {
                            int p = idx_abs(facePolys[k]);
                            bool seenp = false;
                            for (int t = 0; t < np; ++t) if (polys[t] == p) { seenp = true; break; }
                            if (!seenp && np < 32) polys[np++] = p;
                        }
                    }
                }
                dvec3 edgeAvg(0, 0, 0);
                for (int i = 0; i < ne; ++i) edgeAvg += Ecentroids[edges[i]];
                if (ne > 0) edgeAvg /= double(ne);

                dvec3 faceAvg(0, 0, 0);
                for (int i = 0; i < nf; ++i) faceAvg += Fcentroids[faces[i]];
                if (nf > 0) faceAvg /= double(nf);

                dvec3 polyAvg(0, 0, 0);
                for (int i = 0; i < np; ++i) polyAvg += Pcentroids[polys[i]];
                if (np > 0) polyAvg /= double(np);

                dvec3 vnew = polyAvg + (faceAvg * 3.0) + (edgeAvg * 3.0) + Vpos[v];
                outNewV[v] = vnew / 8.0;
            }
            else {
                // 边界点：只考虑边界边/面
                int edgesB[32]; int ne = 0;
                int facesB[32]; int nf = 0;

                for (int i = 0; i < N && ne < 32; ++i) {
                    int e = idx_abs(vertEdges[begin + i]);
                    if (edgeOnSurf[e]) {
                        edgesB[ne++] = e;
                        uint eb = edgeFacesOffset[e];
                        uint ee = edgeFacesOffset[e + 1];
                        for (uint j = eb; j < ee; ++j) {
                            int f = idx_abs(edgeFaces[j]);
                            if (faceOnSurf[f]) {
                                bool seen = false;
                                for (int t = 0; t < nf; ++t) if (facesB[t] == f) { seen = true; break; }
                                if (!seen && nf < 32) facesB[nf++] = f;
                            }
                        }
                    }
                }

                dvec3 edgeAvg(0, 0, 0);
                for (int i = 0; i < ne; ++i) edgeAvg += Ecentroids[edgesB[i]];
                if (ne > 0) edgeAvg /= double(ne);

                dvec3 faceAvg(0, 0, 0);
                for (int i = 0; i < nf; ++i) faceAvg += Fcentroids[facesB[i]];
                if (nf > 0) faceAvg /= double(nf);

                int n = ne; // 只考虑边界边的度
                dvec3 vnew = faceAvg + (edgeAvg * 2.0) + (Vpos[v] * double(n - 3));
                if (n > 0) vnew /= double(n);
                outNewV[v] = vnew;
            }
        }
    };

    // ============ 主函数实现 ============

    void Patch::singlePatch::subdiv_cuda(std::vector<vec3d>& pos, std::vector<uint>& polys)
    {
        // ------------ 一些规模参数 ------------
        const uint ne = static_cast<uint>(edgeVerts.size());
        const uint nf = static_cast<uint>(faceEdges.size() / 4);
        const uint np = static_cast<uint>(polyFaces.size() / 6);
        const uint nv = static_cast<uint>(vertsPos.size());

        // ------------- 主机 -> 设备 拷贝/整理 -------------
        // 顶点
        thrust::host_vector<dvec3> hV(nv);
        for (uint i = 0; i < nv; ++i) hV[i] = to_d(vertsPos[i]);

        // 边的端点（转为 uint2）
        thrust::host_vector<uint2> hEV(ne);
        for (uint e = 0; e < ne; ++e) {
            auto u = edgeVerts[e].x();
            auto v = edgeVerts[e].y();
            hEV[e] = make_uint2(u, v);
        }

        // 其它索引/偏移/标记
        thrust::device_vector<dvec3> dV = hV;

        thrust::device_vector<uint2> dEV = hEV;
        thrust::device_vector<int> d_faceEdges(faceEdges.begin(), faceEdges.end());
        thrust::device_vector<int> d_polyFaces(polyFaces.begin(), polyFaces.end());
        thrust::device_vector<int> d_vertEdges(vertEdges.begin(), vertEdges.end());
        thrust::device_vector<int> d_edgeFaces(edgeFaces.begin(), edgeFaces.end());
        thrust::device_vector<int> d_facePolys(facePolys.begin(), facePolys.end());

        thrust::device_vector<uint> d_veOff(vertEdgesOffset.begin(), vertEdgesOffset.end());
        thrust::device_vector<uint> d_efOff(edgeFacesOffset.begin(), edgeFacesOffset.end());
        thrust::device_vector<uint> d_fpOff(facePolysOffset.begin(), facePolysOffset.end());

        // 布尔标记压成 uint8_t，便于设备侧使用
        auto pack_bool = [](const std::vector<bool>& v) {
            thrust::host_vector<uint8_t> out(v.size());
            for (size_t i = 0; i < v.size(); ++i) out[i] = v[i] ? 1u : 0u;
            return out;
        };
        thrust::device_vector<uint8_t> d_vOn = pack_bool(vertOnSurf);
        thrust::device_vector<uint8_t> d_eOn = pack_bool(edgeOnSurf);
        thrust::device_vector<uint8_t> d_fOn = pack_bool(faceOnSurf);

        // ------------- 设备侧输出缓冲 -------------
        thrust::device_vector<dvec3> d_EC(ne); // edge centroids
        thrust::device_vector<dvec3> d_FC(nf); // face centroids
        thrust::device_vector<dvec3> d_PC(np); // poly centroids

        // 新点：
        thrust::device_vector<dvec3> d_newPoly(patchPolys);  // 体点（取前 patchPolys 个体心）
        thrust::device_vector<dvec3> d_newFace(patchFaces);  // 面点
        thrust::device_vector<dvec3> d_newEdge(patchEdges);  // 边点
        thrust::device_vector<dvec3> d_newVert(patchVerts);  // 顶点

        // ----------- GPU 并行计算 -----------
        auto c_begin_e = thrust::make_counting_iterator<int>(0);
        auto c_begin_f = thrust::make_counting_iterator<int>(0);
        auto c_begin_p = thrust::make_counting_iterator<int>(0);
        auto c_begin_v = thrust::make_counting_iterator<int>(0);

        // (1) 边心
        thrust::for_each(c_begin_e, c_begin_e + int(ne),
            EdgeCentroidOp(thrust::raw_pointer_cast(dEV.data()),
                thrust::raw_pointer_cast(dV.data()),
                thrust::raw_pointer_cast(d_EC.data())));

        // (2) 面心
        thrust::for_each(c_begin_f, c_begin_f + int(nf),
            FaceCentroidOp(thrust::raw_pointer_cast(d_faceEdges.data()),
                thrust::raw_pointer_cast(d_EC.data()),
                thrust::raw_pointer_cast(d_FC.data())));

        // (3) 体心
        thrust::for_each(c_begin_p, c_begin_p + int(np),
            PolyCentroidOp(thrust::raw_pointer_cast(d_polyFaces.data()),
                thrust::raw_pointer_cast(d_FC.data()),
                thrust::raw_pointer_cast(d_PC.data())));

        // (4) 新体点：直接取前 patchPolys 个体心
        thrust::copy(d_PC.begin(), d_PC.begin() + patchPolys, d_newPoly.begin());

        // (5) 新面点
        thrust::for_each(c_begin_f, c_begin_f + int(patchFaces),
            NewFaceVertOp(thrust::raw_pointer_cast(d_fOn.data()),
                thrust::raw_pointer_cast(d_fpOff.data()),
                thrust::raw_pointer_cast(d_facePolys.data()),
                thrust::raw_pointer_cast(d_FC.data()),
                thrust::raw_pointer_cast(d_PC.data()),
                thrust::raw_pointer_cast(d_newFace.data())));

        // (6) 新边点
        thrust::for_each(c_begin_e, c_begin_e + int(patchEdges),
            NewEdgeVertOp(thrust::raw_pointer_cast(d_eOn.data()),
                thrust::raw_pointer_cast(d_fOn.data()),
                thrust::raw_pointer_cast(d_efOff.data()),
                thrust::raw_pointer_cast(d_edgeFaces.data()),
                thrust::raw_pointer_cast(d_fpOff.data()),
                thrust::raw_pointer_cast(d_facePolys.data()),
                thrust::raw_pointer_cast(d_EC.data()),
                thrust::raw_pointer_cast(d_FC.data()),
                thrust::raw_pointer_cast(d_PC.data()),
                thrust::raw_pointer_cast(d_newEdge.data())));

        // (7) 新顶点
        thrust::for_each(c_begin_v, c_begin_v + int(patchVerts),
            NewVertVertOp(thrust::raw_pointer_cast(d_vOn.data()),
                thrust::raw_pointer_cast(d_eOn.data()),
                thrust::raw_pointer_cast(d_fOn.data()),
                thrust::raw_pointer_cast(d_veOff.data()),
                thrust::raw_pointer_cast(d_vertEdges.data()),
                thrust::raw_pointer_cast(d_efOff.data()),
                thrust::raw_pointer_cast(d_edgeFaces.data()),
                thrust::raw_pointer_cast(d_fpOff.data()),
                thrust::raw_pointer_cast(d_facePolys.data()),
                thrust::raw_pointer_cast(dV.data()),
                thrust::raw_pointer_cast(d_EC.data()),
                thrust::raw_pointer_cast(d_FC.data()),
                thrust::raw_pointer_cast(d_PC.data()),
                thrust::raw_pointer_cast(d_newVert.data())));

        // ----------- 把新点拷回主机并写入 pos -----------
        thrust::host_vector<dvec3> h_newPoly = d_newPoly;
        thrust::host_vector<dvec3> h_newFace = d_newFace;
        thrust::host_vector<dvec3> h_newEdge = d_newEdge;
        thrust::host_vector<dvec3> h_newVert = d_newVert;

        // 计算偏移（插入前的 pos.size()）
        uint pvOffset = static_cast<uint>(pos.size());
        uint fvOffset = pvOffset + static_cast<uint>(h_newPoly.size());
        uint evOffset = fvOffset + static_cast<uint>(h_newFace.size());
        uint vvOffset = evOffset + static_cast<uint>(h_newEdge.size());

        // 追加新点
        pos.reserve(pos.size() + h_newPoly.size() + h_newFace.size() + h_newEdge.size() + h_newVert.size());
        for (const auto& v : h_newPoly) pos.push_back(to_h(v));
        for (const auto& v : h_newFace) pos.push_back(to_h(v));
        for (const auto& v : h_newEdge) pos.push_back(to_h(v));
        for (const auto& v : h_newVert) pos.push_back(to_h(v));

        // ----------- 主机侧：复用你原先的拓扑装配逻辑 -----------
        // 注意：这里使用已有的成员数组（edgeVerts/faceEdges/polyFaces/...）
        //       以及我们刚刚计算好的偏移 pvOffset/fvOffset/evOffset/vvOffset。

        std::vector<uint> tempPolys;
        std::vector<uint> tempFaces;
        std::vector<uint> tempEdges;
        std::vector<uint> tempVerts;

        std::vector<uint> FV(6);
        std::vector<uint> EV(12);
        std::vector<uint> VV(8);
        bool f1Reverse = false;

        for (int p = 0; p < int(patchPolys); ++p) {
            int f1 = polyFaces[p * 6 + 0];
            f1Reverse = false;
            if (f1 < 0) { f1 = -f1 - 1; f1Reverse = true; }

            tempEdges.clear();
            for (int edgeOff = 0; edgeOff < 4; ++edgeOff) {
                int e = faceEdges[f1 * 4 + edgeOff];
                if (e < 0) e = -e - 1;
                tempEdges.push_back(uint(e));
            }

            int f2 = -1;
            for (int faceOff = 1; faceOff < 6; ++faceOff) {
                int f = polyFaces[p * 6 + faceOff];
                if (f < 0) f = -f - 1;
                bool flag = true;
                for (int edgeOff = 0; edgeOff < 4; ++edgeOff) {
                    int e = faceEdges[f * 4 + edgeOff];
                    if (e < 0) e = -e - 1;
                    if (std::find(tempEdges.begin(), tempEdges.end(), uint(e)) != tempEdges.end()) {
                        flag = false;
                    }
                }
                if (flag) { f2 = f; break; }
            }

            tempVerts.clear();
            tempEdges.clear();
            tempFaces.clear();
            for (int faceOff = 0; faceOff < 6; ++faceOff) {
                int f = polyFaces[p * 6 + faceOff];
                if (f < 0) f = -f - 1;
                if (std::find(tempFaces.begin(), tempFaces.end(), uint(f)) == tempFaces.end()) {
                    tempFaces.push_back(uint(f));
                }
                for (int edgeOff = 0; edgeOff < 4; ++edgeOff) {
                    int e = faceEdges[f * 4 + edgeOff];
                    if (e < 0) e = -e - 1;
                    if (std::find(tempEdges.begin(), tempEdges.end(), uint(e)) == tempEdges.end()) {
                        tempEdges.push_back(uint(e));
                    }
                    auto start = edgeVerts[e].x();
                    auto end = edgeVerts[e].y();
                    if (std::find(tempVerts.begin(), tempVerts.end(), start) == tempVerts.end()) tempVerts.push_back(start);
                    if (std::find(tempVerts.begin(), tempVerts.end(), end) == tempVerts.end()) tempVerts.push_back(end);
                }
            }

            FV.clear(); EV.clear(); VV.clear();
            bool edgeReverse = false;
            int vStart, vCurrent;
            int eCurrent = faceEdges[f1 * 4 + 0];
            if (eCurrent < 0) { eCurrent = -eCurrent - 1; edgeReverse = true; }

            if ((f1Reverse ^ edgeReverse)) {
                vStart = edgeVerts[eCurrent].y();
                vCurrent = edgeVerts[eCurrent].x();
            }
            else {
                vStart = edgeVerts[eCurrent].x();
                vCurrent = edgeVerts[eCurrent].y();
            }
            VV.push_back(uint(vStart));
            EV.push_back(uint(eCurrent));

            while (vCurrent != vStart) {
                VV.push_back(uint(vCurrent));
                for (int eOff = 1; eOff < 4; ++eOff) {
                    edgeReverse = false;
                    eCurrent = faceEdges[f1 * 4 + eOff];
                    if (eCurrent < 0) { eCurrent = -eCurrent - 1; edgeReverse = true; }
                    if ((f1Reverse ^ edgeReverse) && int(edgeVerts[eCurrent].y()) == vCurrent) {
                        vCurrent = edgeVerts[eCurrent].x();
                        EV.push_back(uint(eCurrent));
                        break;
                    }
                    else if (!(f1Reverse ^ edgeReverse) && int(edgeVerts[eCurrent].x()) == vCurrent) {
                        vCurrent = edgeVerts[eCurrent].y();
                        EV.push_back(uint(eCurrent));
                        break;
                    }
                }
            }

            for (int i = 0; i < 4; ++i) {
                vCurrent = int(VV[i]);
                for (auto& e : tempEdges) {
                    if (edgeVerts[e].x() == uint(vCurrent) &&
                        std::find(EV.begin(), EV.end(), e) == EV.end()) {
                        VV.push_back(edgeVerts[e].y());
                        EV.push_back(e);
                        break;
                    }
                    else if (edgeVerts[e].y() == uint(vCurrent) &&
                        std::find(EV.begin(), EV.end(), e) == EV.end()) {
                        VV.push_back(edgeVerts[e].x());
                        EV.push_back(e);
                        break;
                    }
                }
            }

            FV.push_back(uint(f1));
            FV.push_back(uint(f2));

            for (int i = 0; i < 4; ++i) {
                eCurrent = int(EV[i]);
                for (auto& f : tempFaces) {
                    if (std::find(FV.begin(), FV.end(), f) == FV.end()) {
                        for (int eOff = 0; eOff < 4; ++eOff) {
                            int fe = faceEdges[f * 4 + eOff];
                            if (fe < 0) fe = -fe - 1;
                            if (fe == eCurrent) {
                                FV.push_back(f);
                                for (int eOff2 = 0; eOff2 < 4; ++eOff2) {
                                    fe = faceEdges[f * 4 + eOff2];
                                    if (fe < 0) fe = -fe - 1;
                                    if (std::find(EV.begin(), EV.end(), uint(fe)) == EV.end()) {
                                        EV.push_back(uint(fe));
                                        break;
                                    }
                                }
                                break;
                            }
                        }
                    }
                }
            }

            uint PV = uint(p) + pvOffset;
            for (int i = 0; i<int(FV.size()); ++i) FV[i] += fvOffset;
            for (int i = 0; i<int(EV.size()); ++i) EV[i] += evOffset;
            for (int i = 0; i<int(VV.size()); ++i) VV[i] += vvOffset;

            polys.insert(polys.end(), { VV[0], EV[0], FV[0], EV[3], EV[4], FV[2], PV,    FV[5] });
            polys.insert(polys.end(), { EV[0], VV[1], EV[1], FV[0], FV[2], EV[5], FV[3], PV });
            polys.insert(polys.end(), { FV[0], EV[1], VV[2], EV[2], PV,    FV[3], EV[6], FV[4] });
            polys.insert(polys.end(), { EV[3], FV[0], EV[2], VV[3], FV[5], PV,    FV[4], EV[7] });

            polys.insert(polys.end(), { EV[4], FV[2], PV,    FV[5], VV[4], EV[8], FV[1], EV[11] });
            polys.insert(polys.end(), { FV[2], EV[5], FV[3], PV,    EV[8], VV[5], EV[9], FV[1] });
            polys.insert(polys.end(), { PV,    FV[3], EV[6], FV[4], FV[1], EV[9], VV[6], EV[10] });
            polys.insert(polys.end(), { FV[5], PV,    FV[4], EV[7], EV[11],FV[1], EV[10],VV[7] });
        }
    }
    
    static inline dvec3 to_d(const vec3d& v) {
        return dvec3(v.x(), v.y(), v.z());
    }
    static inline vec3d to_h(const dvec3& v) {
        return cinolib::vec3d(v.x, v.y, v.z);
    }
}