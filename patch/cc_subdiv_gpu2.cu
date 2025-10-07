#define CINO_STATIC_LIB
#include "patch.h"
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/transform.h>
#include <thrust/for_each.h>
#include <thrust/copy.h>
#include <thrust/sequence.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/execution_policy.h>
#include <thrust/functional.h>
#include <cuda_runtime.h>
#include <vector>
#include <algorithm>
#include <cstdint>
#include <chrono>

namespace cinolib {

    // ======= 32B 对齐的 3D 向量（改善合并访存；w 仅为填充）=======
    struct alignas(32) dvec3a {
        double x, y, z, w;
        __host__ __device__ dvec3a() : x(0), y(0), z(0), w(0) {}
        __host__ __device__ dvec3a(double X, double Y, double Z) : x(X), y(Y), z(Z), w(0) {}
    };

    __host__ __device__ __forceinline__ dvec3a operator+(const dvec3a& a, const dvec3a& b) {
        return dvec3a(a.x + b.x, a.y + b.y, a.z + b.z);
    }
    __host__ __device__ __forceinline__ dvec3a operator-(const dvec3a& a, const dvec3a& b) {
        return dvec3a(a.x - b.x, a.y - b.y, a.z - b.z);
    }
    __host__ __device__ __forceinline__ dvec3a operator*(const dvec3a& a, double s) {
        return dvec3a(a.x * s, a.y * s, a.z * s);
    }
    __host__ __device__ __forceinline__ dvec3a operator*(double s, const dvec3a& a) {
        return dvec3a(a.x * s, a.y * s, a.z * s);
    }
    __host__ __device__ __forceinline__ dvec3a operator/(const dvec3a& a, double s) {
        return dvec3a(a.x / s, a.y / s, a.z / s);
    }
    __host__ __device__ __forceinline__ dvec3a& operator+=(dvec3a& a, const dvec3a& b) {
        a.x += b.x; a.y += b.y; a.z += b.z; return a;
    }

    __host__ __device__ __forceinline__ int idx_abs(int v) { return (v < 0) ? (-v - 1) : v; }

    static inline dvec3a to_d(const cinolib::vec3d& v);
    static inline vec3d  to_h(const dvec3a& v);

    // ========== 辅助：读取标记的谓词（用于 copy_if 分割边界/内部） ==========
    struct IsBoundaryPred {
        const uint8_t* __restrict__ flag;
        __host__ __device__ IsBoundaryPred(const uint8_t* f) : flag(f) {}
        __host__ __device__ bool operator()(const int i) const { return flag[i] != 0; }
    };
    struct IsInteriorPred {
        const uint8_t* __restrict__ flag;
        __host__ __device__ IsInteriorPred(const uint8_t* f) : flag(f) {}
        __host__ __device__ bool operator()(const int i) const { return flag[i] == 0; }
    };

    // ========== 1) 边质心 ==========
    struct EdgeCentroidOp {
        const uint2* __restrict__ EV;     // edgeVerts
        const dvec3a* __restrict__ Vpos;  // 顶点坐标
        dvec3a* __restrict__ Ecentroids;  // 输出
        __host__ __device__
            EdgeCentroidOp(const uint2* ev, const dvec3a* vp, dvec3a* out)
            : EV(ev), Vpos(vp), Ecentroids(out) {}
        __host__ __device__ __forceinline__
            void operator()(const int e) const {
            const uint2 ev = EV[e];
            const dvec3a a = Vpos[ev.x];
            const dvec3a b = Vpos[ev.y];
            Ecentroids[e] = (a + b) * 0.5;
        }
    };

    // ========== 2) 面质心（四条边心平均） ==========
    struct FaceCentroidOp {
        const int* __restrict__ faceEdges; // nf*4
        const dvec3a* __restrict__ Ecentroids;
        dvec3a* __restrict__ Fcentroids;
        __host__ __device__
            FaceCentroidOp(const int* fe, const dvec3a* ec, dvec3a* out)
            : faceEdges(fe), Ecentroids(ec), Fcentroids(out) {}
        __host__ __device__ __forceinline__
            void operator()(const int f) const {
            dvec3a c(0, 0, 0);
#pragma unroll
            for (int k = 0; k < 4; ++k) {
                const int e = idx_abs(faceEdges[f * 4 + k]);
                c += Ecentroids[e];
            }
            Fcentroids[f] = c / 4.0;
        }
    };

    // ========== 3) 体质心（六个面心平均） ==========
    struct PolyCentroidOp {
        const int* __restrict__ polyFaces; // np*6
        const dvec3a* __restrict__ Fcentroids;
        dvec3a* __restrict__ Pcentroids;
        __host__ __device__
            PolyCentroidOp(const int* pf, const dvec3a* fc, dvec3a* out)
            : polyFaces(pf), Fcentroids(fc), Pcentroids(out) {}
        __host__ __device__ __forceinline__
            void operator()(const int p) const {
            dvec3a c(0, 0, 0);
#pragma unroll
            for (int k = 0; k < 6; ++k) {
                const int f = idx_abs(polyFaces[p * 6 + k]);
                c += Fcentroids[f];
            }
            Pcentroids[p] = c / 6.0;
        }
    };

    // ========== 4) 新面点（拆分：内部/边界两套 kernel） ==========
    struct NewFaceVert_Interior_Op {
        const uint* __restrict__ facePolysOffset;
        const int* __restrict__ facePolys;
        const dvec3a* __restrict__ Fcentroids;
        const dvec3a* __restrict__ Pcentroids;
        dvec3a* __restrict__ outNewF;
        __host__ __device__
            NewFaceVert_Interior_Op(const uint* fpoff, const int* fp,
                const dvec3a* fc, const dvec3a* pc, dvec3a* out)
            : facePolysOffset(fpoff), facePolys(fp),
            Fcentroids(fc), Pcentroids(pc), outNewF(out) {}
        __host__ __device__ __forceinline__
            void operator()(const int f) const {
            const uint off = facePolysOffset[f];
            const int p0 = idx_abs(facePolys[off + 0]);
            const int p1 = idx_abs(facePolys[off + 1]);
            const dvec3a v = Pcentroids[p0] + Pcentroids[p1] + (Fcentroids[f] * 2.0);
            outNewF[f] = v / 4.0;
        }
    };

    struct NewFaceVert_Boundary_Op {
        const dvec3a* __restrict__ Fcentroids;
        dvec3a* __restrict__ outNewF;
        __host__ __device__
            NewFaceVert_Boundary_Op(const dvec3a* fc, dvec3a* out)
            : Fcentroids(fc), outNewF(out) {}
        __host__ __device__ __forceinline__
            void operator()(const int f) const {
            outNewF[f] = Fcentroids[f];
        }
    };

    // ========== 5) 新边点（拆分：内部/边界） ==========
    struct NewEdgeVert_Interior_Op {
        const uint* __restrict__ edgeFacesOffset;
        const int* __restrict__ edgeFaces;
        const uint* __restrict__ facePolysOffset;
        const int* __restrict__ facePolys;

        const dvec3a* __restrict__ Ecentroids;
        const dvec3a* __restrict__ Fcentroids;
        const dvec3a* __restrict__ Pcentroids;

        dvec3a* __restrict__ outNewE;

        __host__ __device__
            NewEdgeVert_Interior_Op(const uint* efoff, const int* ef,
                const uint* fpoff, const int* fp,
                const dvec3a* ec, const dvec3a* fc, const dvec3a* pc,
                dvec3a* out)
            : edgeFacesOffset(efoff), edgeFaces(ef),
            facePolysOffset(fpoff), facePolys(fp),
            Ecentroids(ec), Fcentroids(fc), Pcentroids(pc), outNewE(out) {}

        __host__ __device__ __forceinline__
            bool contains_i(const int* arr, int n, int v) const {
#pragma unroll
            for (int i = 0; i < n; ++i) if (arr[i] == v) return true;
            return false;
        }

        __host__ __device__ __forceinline__
            void operator()(const int e) const {
            const uint begin = edgeFacesOffset[e];
            const uint end = edgeFacesOffset[e + 1];
            const int  N = int(end - begin);

            // 收集相邻面（最多 32）
            int faces[32]; int nf = 0;
#pragma unroll
            for (int i = 0; i < N && nf < 32; ++i) {
                const int f = idx_abs(edgeFaces[begin + i]);
                faces[nf++] = f;
            }
            // 面 -> 体（去重，最多 32）
            int polys[32]; int np = 0;
#pragma unroll
            for (int i = 0; i < nf; ++i) {
                const int f = faces[i];
                const uint fb = facePolysOffset[f];
                const uint fe = facePolysOffset[f + 1];
                for (uint k = fb; k < fe; ++k) {
                    const int p = idx_abs(facePolys[k]);
                    if (!contains_i(polys, np, p) && np < 32) polys[np++] = p;
                }
            }
            dvec3a faceAvg(0, 0, 0);
#pragma unroll
            for (int i = 0; i < nf; ++i) faceAvg += Fcentroids[faces[i]];
            if (nf > 0) faceAvg = faceAvg / double(nf);

            dvec3a polyAvg(0, 0, 0);
#pragma unroll
            for (int i = 0; i < np; ++i) polyAvg += Pcentroids[polys[i]];
            if (np > 0) polyAvg = polyAvg / double(np);

            const dvec3a v = polyAvg + (faceAvg * 2.0) + (Ecentroids[e] * double(N - 3));
            outNewE[e] = v / double(N);
        }
    };

    struct NewEdgeVert_Boundary_Op {
        const uint8_t* __restrict__ faceOnSurf;
        const uint* __restrict__ edgeFacesOffset;
        const int* __restrict__ edgeFaces;

        const dvec3a* __restrict__ Ecentroids;
        const dvec3a* __restrict__ Fcentroids;

        dvec3a* __restrict__ outNewE;

        __host__ __device__
            NewEdgeVert_Boundary_Op(const uint8_t* fos,
                const uint* efoff, const int* ef,
                const dvec3a* ec, const dvec3a* fc,
                dvec3a* out)
            : faceOnSurf(fos), edgeFacesOffset(efoff), edgeFaces(ef),
            Ecentroids(ec), Fcentroids(fc), outNewE(out) {}

        __host__ __device__ __forceinline__
            void operator()(const int e) const {
            const uint begin = edgeFacesOffset[e];
            const uint end = edgeFacesOffset[e + 1];
            const int  N = int(end - begin);

            if (N == 1) { // 开放边
                outNewE[e] = Ecentroids[e];
                return;
            }
            // 只平均边界面
            int facesB[32]; int nf = 0;
            for (int i = 0; i < N && nf < 32; ++i) {
                const int f = idx_abs(edgeFaces[begin + i]);
                if (faceOnSurf[f]) facesB[nf++] = f;
            }
            dvec3a faceAvg(0, 0, 0);
            for (int i = 0; i < nf; ++i) faceAvg += Fcentroids[facesB[i]];
            if (nf > 0) faceAvg = faceAvg / double(nf);
            const dvec3a v = faceAvg + Ecentroids[e];
            outNewE[e] = v * 0.5;
        }
    };

    // ========== 6) 新顶点（拆分：内部/边界） ==========
    struct NewVertVert_Interior_Op {
        const uint* __restrict__ vertEdgesOffset;
        const int* __restrict__ vertEdges;

        const uint* __restrict__ edgeFacesOffset;
        const int* __restrict__ edgeFaces;

        const uint* __restrict__ facePolysOffset;
        const int* __restrict__ facePolys;

        const dvec3a* __restrict__ Vpos;
        const dvec3a* __restrict__ Ecentroids;
        const dvec3a* __restrict__ Fcentroids;
        const dvec3a* __restrict__ Pcentroids;

        dvec3a* __restrict__ outNewV;

        __host__ __device__
            NewVertVert_Interior_Op(const uint* veoff, const int* ve,
                const uint* efoff, const int* ef,
                const uint* fpoff, const int* fp,
                const dvec3a* vp, const dvec3a* ec, const dvec3a* fc, const dvec3a* pc,
                dvec3a* out)
            : vertEdgesOffset(veoff), vertEdges(ve),
            edgeFacesOffset(efoff), edgeFaces(ef),
            facePolysOffset(fpoff), facePolys(fp),
            Vpos(vp), Ecentroids(ec), Fcentroids(fc), Pcentroids(pc), outNewV(out) {}

        __host__ __device__ __forceinline__
            bool contains_i(const int* arr, int n, int v) const {
#pragma unroll
            for (int i = 0; i < n; ++i) if (arr[i] == v) return true;
            return false;
        }

        __host__ __device__ __forceinline__
            void operator()(const int v) const {
            const uint begin = vertEdgesOffset[v];
            const uint end = vertEdgesOffset[v + 1];
            const int  N = int(end - begin);

            int edges[32]; int ne = 0;
            int faces[32]; int nf = 0;
            int polys[32]; int np = 0;

            // 边 / 面 / 体 收集 + 去重
#pragma unroll
            for (int i = 0; i < N && ne < 32; ++i) {
                const int e = idx_abs(vertEdges[begin + i]);
                edges[ne++] = e;
                // edge -> faces
                const uint eb = edgeFacesOffset[e];
                const uint ee = edgeFacesOffset[e + 1];
                for (uint j = eb; j < ee; ++j) {
                    const int f = idx_abs(edgeFaces[j]);
                    if (!contains_i(faces, nf, f) && nf < 32) {
                        faces[nf++] = f;
                    }
                    // face -> polys
                    const uint fb = facePolysOffset[f];
                    const uint fe = facePolysOffset[f + 1];
                    for (uint k = fb; k < fe; ++k) {
                        const int p = idx_abs(facePolys[k]);
                        if (!contains_i(polys, np, p) && np < 32) {
                            polys[np++] = p;
                        }
                    }
                }
            }

            dvec3a edgeAvg(0, 0, 0);
            for (int i = 0; i < ne; ++i) edgeAvg += Ecentroids[edges[i]];
            if (ne > 0) edgeAvg = edgeAvg / double(ne);

            dvec3a faceAvg(0, 0, 0);
            for (int i = 0; i < nf; ++i) faceAvg += Fcentroids[faces[i]];
            if (nf > 0) faceAvg = faceAvg / double(nf);

            dvec3a polyAvg(0, 0, 0);
            for (int i = 0; i < np; ++i) polyAvg += Pcentroids[polys[i]];
            if (np > 0) polyAvg = polyAvg / double(np);

            const dvec3a vnew = (polyAvg + (faceAvg * 3.0) + (edgeAvg * 3.0) + Vpos[v]) / 8.0;
            outNewV[v] = vnew;
        }
    };

    struct NewVertVert_Boundary_Op {
        const uint8_t* __restrict__ edgeOnSurf;
        const uint8_t* __restrict__ faceOnSurf;

        const uint* __restrict__ vertEdgesOffset;
        const int* __restrict__ vertEdges;

        const uint* __restrict__ edgeFacesOffset;
        const int* __restrict__ edgeFaces;

        const dvec3a* __restrict__ Vpos;
        const dvec3a* __restrict__ Ecentroids;
        const dvec3a* __restrict__ Fcentroids;

        dvec3a* __restrict__ outNewV;

        __host__ __device__
            NewVertVert_Boundary_Op(const uint8_t* eos, const uint8_t* fos,
                const uint* veoff, const int* ve,
                const uint* efoff, const int* ef,
                const dvec3a* vp, const dvec3a* ec, const dvec3a* fc,
                dvec3a* out)
            : edgeOnSurf(eos), faceOnSurf(fos),
            vertEdgesOffset(veoff), vertEdges(ve),
            edgeFacesOffset(efoff), edgeFaces(ef),
            Vpos(vp), Ecentroids(ec), Fcentroids(fc), outNewV(out) {}

        __host__ __device__ __forceinline__
            void operator()(const int v) const {
            const uint begin = vertEdgesOffset[v];
            const uint end = vertEdgesOffset[v + 1];

            int edgesB[32]; int ne = 0;
            int facesB[32]; int nf = 0;

            // 只考虑边界边/面
            for (uint i = begin; i < end && ne < 32; ++i) {
                const int e = idx_abs(vertEdges[i]);
                if (edgeOnSurf[e]) {
                    edgesB[ne++] = e;
                    const uint eb = edgeFacesOffset[e];
                    const uint ee = edgeFacesOffset[e + 1];
                    for (uint j = eb; j < ee; ++j) {
                        const int f = idx_abs(edgeFaces[j]);
                        if (faceOnSurf[f]) {
                            bool seen = false;
                            for (int t = 0; t < nf; ++t) if (facesB[t] == f) { seen = true; break; }
                            if (!seen && nf < 32) facesB[nf++] = f;
                        }
                    }
                }
            }

            dvec3a edgeAvg(0, 0, 0);
            for (int i = 0; i < ne; ++i) edgeAvg += Ecentroids[edgesB[i]];
            if (ne > 0) edgeAvg = edgeAvg / double(ne);

            dvec3a faceAvg(0, 0, 0);
            for (int i = 0; i < nf; ++i) faceAvg += Fcentroids[facesB[i]];
            if (nf > 0) faceAvg = faceAvg / double(nf);

            const int n = ne; // 边界点度数
            dvec3a vnew = faceAvg + (edgeAvg * 2.0) + (Vpos[v] * double(n - 3));
            if (n > 0) vnew = vnew / double(n);
            outNewV[v] = vnew;
        }
    };

    // ========== 7) Topo 装配 ==========
    struct TopoAssembleOp {
        const int* __restrict__ polyFaces; // np*6 (带符号)
        const int* __restrict__ faceEdges; // nf*4 (带符号)
        const uint2* __restrict__ edgeVerts; // ne

        const uint pvOffset, fvOffset, evOffset, vvOffset;
        uint* __restrict__ out;              // patchPolys*64

        __host__ __device__
            TopoAssembleOp(const int* pf, const int* fe, const uint2* ev,
                uint pv, uint fv, uint evv, uint vv, uint* o)
            : polyFaces(pf), faceEdges(fe), edgeVerts(ev),
            pvOffset(pv), fvOffset(fv), evOffset(evv), vvOffset(vv), out(o) {}

        __host__ __device__ __forceinline__ static int iabs(int v) { return v < 0 ? (-v - 1) : v; }
        __host__ __device__ __forceinline__ static bool contains_u(const uint* arr, int n, uint v) {
#pragma unroll
            for (int i = 0; i < n; ++i) if (arr[i] == v) return true;
            return false;
        }

        __host__ __device__ __forceinline__
            void operator()(const int p) const {
            uint VV[8];  int nVV = 0;
            uint EV[12]; int nEV = 0;
            uint FV[6];  int nFV = 0;

            int FACES[6];
#pragma unroll
            for (int fo = 0; fo < 6; ++fo) {
                FACES[fo] = iabs(polyFaces[p * 6 + fo]);
            }

            const int raw_f1 = polyFaces[p * 6 + 0];
            const bool f1Reverse = (raw_f1 < 0);
            const int f1 = iabs(raw_f1);

            int F1E[4];
#pragma unroll
            for (int eOff = 0; eOff < 4; ++eOff) {
                F1E[eOff] = iabs(faceEdges[f1 * 4 + eOff]);
            }

            int f2 = -1;
            for (int fo = 1; fo < 6; ++fo) {
                const int f = FACES[fo];
                bool share = false;
#pragma unroll
                for (int eOff = 0; eOff < 4; ++eOff) {
                    const int e = iabs(faceEdges[f * 4 + eOff]);
                    for (int k = 0; k < 4; ++k) { if (e == F1E[k]) { share = true; break; } }
                    if (share) break;
                }
                if (!share) { f2 = f; break; }
            }

            uint EDGES[12]; int nEDGES = 0;
            for (int fo = 0; fo < 6; ++fo) {
                const int f = FACES[fo];
#pragma unroll
                for (int eOff = 0; eOff < 4; ++eOff) {
                    const uint e = (uint)iabs(faceEdges[f * 4 + eOff]);
                    if (!contains_u(EDGES, nEDGES, e) && nEDGES < 12) {
                        EDGES[nEDGES++] = e;
                    }
                }
            }

            bool edgeReverse = false;
            int eCurrent = faceEdges[f1 * 4 + 0];
            if (eCurrent < 0) { eCurrent = -eCurrent - 1; edgeReverse = true; }

            int vStart, vCurrent;
            if ((f1Reverse ^ edgeReverse)) {
                vStart = edgeVerts[(uint)eCurrent].y;
                vCurrent = edgeVerts[(uint)eCurrent].x;
            }
            else {
                vStart = edgeVerts[(uint)eCurrent].x;
                vCurrent = edgeVerts[(uint)eCurrent].y;
            }
            VV[nVV++] = (uint)vStart;
            EV[nEV++] = (uint)eCurrent;

            while (vCurrent != vStart) {
                VV[nVV++] = (uint)vCurrent;
                for (int eOff = 1; eOff < 4; ++eOff) {
                    edgeReverse = false;
                    int eTmp = faceEdges[f1 * 4 + eOff];
                    if (eTmp < 0) { eTmp = -eTmp - 1; edgeReverse = true; }
                    const uint2 ev = edgeVerts[(uint)eTmp];
                    if ((f1Reverse ^ edgeReverse) && (int)ev.y == vCurrent) {
                        vCurrent = ev.x; EV[nEV++] = (uint)eTmp; break;
                    }
                    else if (!(f1Reverse ^ edgeReverse) && (int)ev.x == vCurrent) {
                        vCurrent = ev.y; EV[nEV++] = (uint)eTmp; break;
                    }
                }
            }

            for (int i = 0; i < 4; ++i) {
                const int vtx = (int)VV[i];
                for (int t = 0; t < nEDGES; ++t) {
                    const uint e = EDGES[t];
                    if (!contains_u(EV, nEV, e)) {
                        const uint2 ev = edgeVerts[e];
                        if ((int)ev.x == vtx) { VV[nVV++] = ev.y; EV[nEV++] = e; break; }
                        else if ((int)ev.y == vtx) { VV[nVV++] = ev.x; EV[nEV++] = e; break; }
                    }
                }
            }

            FV[nFV++] = (uint)f1;
            FV[nFV++] = (uint)f2;

            for (int i = 0; i < 4; ++i) {
                const int eNeed = (int)EV[i];
                for (int fo = 0; fo < 6; ++fo) {
                    const uint f = (uint)FACES[fo];
                    if (!contains_u(FV, nFV, f)) {
                        bool hit = false;
                        for (int eOff = 0; eOff < 4; ++eOff) {
                            const int fe = iabs(faceEdges[f * 4 + eOff]);
                            if (fe == eNeed) { hit = true; break; }
                        }
                        if (hit) {
                            FV[nFV++] = f;
                            for (int eOff = 0; eOff < 4; ++eOff) {
                                const int fe = iabs(faceEdges[f * 4 + eOff]);
                                if (!contains_u(EV, nEV, (uint)fe)) {
                                    EV[nEV++] = (uint)fe;
                                    break;
                                }
                            }
                        }
                    }
                }
            }

            const uint PV = (uint)p + pvOffset;
            for (int i = 0; i < 6; ++i)  FV[i] += fvOffset;
            for (int i = 0; i < 12; ++i)  EV[i] += evOffset;
            for (int i = 0; i < 8; ++i)  VV[i] += vvOffset;

            const uint base = (uint)p * 64u;
            // 1
            out[base + 0] = VV[0]; out[base + 1] = EV[0]; out[base + 2] = FV[0]; out[base + 3] = EV[3];
            out[base + 4] = EV[4]; out[base + 5] = FV[2]; out[base + 6] = PV;    out[base + 7] = FV[5];
            // 2
            out[base + 8] = EV[0]; out[base + 9] = VV[1]; out[base + 10] = EV[1]; out[base + 11] = FV[0];
            out[base + 12] = FV[2]; out[base + 13] = EV[5]; out[base + 14] = FV[3]; out[base + 15] = PV;
            // 3
            out[base + 16] = FV[0]; out[base + 17] = EV[1]; out[base + 18] = VV[2]; out[base + 19] = EV[2];
            out[base + 20] = PV;    out[base + 21] = FV[3]; out[base + 22] = EV[6]; out[base + 23] = FV[4];
            // 4
            out[base + 24] = EV[3]; out[base + 25] = FV[0]; out[base + 26] = EV[2]; out[base + 27] = VV[3];
            out[base + 28] = FV[5]; out[base + 29] = PV;    out[base + 30] = FV[4]; out[base + 31] = EV[7];
            // 5
            out[base + 32] = EV[4]; out[base + 33] = FV[2]; out[base + 34] = PV;    out[base + 35] = FV[5];
            out[base + 36] = VV[4]; out[base + 37] = EV[8]; out[base + 38] = FV[1]; out[base + 39] = EV[11];
            // 6
            out[base + 40] = FV[2]; out[base + 41] = EV[5]; out[base + 42] = FV[3]; out[base + 43] = PV;
            out[base + 44] = EV[8]; out[base + 45] = VV[5]; out[base + 46] = EV[9]; out[base + 47] = FV[1];
            // 7
            out[base + 48] = PV;    out[base + 49] = FV[3]; out[base + 50] = EV[6]; out[base + 51] = FV[4];
            out[base + 52] = FV[1]; out[base + 53] = EV[9]; out[base + 54] = VV[6]; out[base + 55] = EV[10];
            // 8
            out[base + 56] = FV[5]; out[base + 57] = PV;    out[base + 58] = FV[4]; out[base + 59] = EV[7];
            out[base + 60] = EV[11]; out[base + 61] = FV[1]; out[base + 62] = EV[10]; out[base + 63] = VV[7];
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

        // 使用单一 stream（后续可考虑捕获 CUDA Graph）
        cudaStream_t stream;
        cudaStreamCreate(&stream);
        auto exec = thrust::cuda::par.on(stream);

        // ------------- 主机 -> 设备 拷贝/整理 -------------
        // 顶点
        thrust::host_vector<dvec3a> hV(nv);
        for (uint i = 0; i < nv; ++i) hV[i] = to_d(vertsPos[i]);

        // 边的端点（转为 uint2）
        thrust::host_vector<uint2> hEV(ne);
        for (uint e = 0; e < ne; ++e) {
            auto u = edgeVerts[e].x();
            auto v = edgeVerts[e].y();
            hEV[e] = make_uint2(u, v);
        }

        // 各种索引/偏移/标记
        thrust::device_vector<dvec3a> dV = hV;
        thrust::device_vector<uint2>  dEV = hEV;

        thrust::device_vector<int>  d_faceEdges(faceEdges.begin(), faceEdges.end());
        thrust::device_vector<int>  d_polyFaces(polyFaces.begin(), polyFaces.end());
        thrust::device_vector<int>  d_vertEdges(vertEdges.begin(), vertEdges.end());
        thrust::device_vector<int>  d_edgeFaces(edgeFaces.begin(), edgeFaces.end());
        thrust::device_vector<int>  d_facePolys(facePolys.begin(), facePolys.end());

        thrust::device_vector<uint> d_veOff(vertEdgesOffset.begin(), vertEdgesOffset.end());
        thrust::device_vector<uint> d_efOff(edgeFacesOffset.begin(), edgeFacesOffset.end());
        thrust::device_vector<uint> d_fpOff(facePolysOffset.begin(), facePolysOffset.end());

        // 布尔标记压成 uint8_t
        auto pack_bool = [](const std::vector<bool>& v) {
            thrust::host_vector<uint8_t> out(v.size());
            for (size_t i = 0; i < v.size(); ++i) out[i] = v[i] ? 1u : 0u;
            return out;
        };
        thrust::device_vector<uint8_t> d_vOn = pack_bool(vertOnSurf);
        thrust::device_vector<uint8_t> d_eOn = pack_bool(edgeOnSurf);
        thrust::device_vector<uint8_t> d_fOn = pack_bool(faceOnSurf);

        // ------------- 设备侧输出缓冲 -------------
        thrust::device_vector<dvec3a> d_EC(ne); // edge centroids
        thrust::device_vector<dvec3a> d_FC(nf); // face centroids
        thrust::device_vector<dvec3a> d_PC(np); // poly centroids

        thrust::device_vector<dvec3a> d_newPoly(patchPolys);  // 体点
        thrust::device_vector<dvec3a> d_newFace(patchFaces);  // 面点
        thrust::device_vector<dvec3a> d_newEdge(patchEdges);  // 边点
        thrust::device_vector<dvec3a> d_newVert(patchVerts);  // 顶点

        // ----------- GPU 并行计算 -----------
        auto cE = thrust::make_counting_iterator<int>(0);
        auto cF = thrust::make_counting_iterator<int>(0);
        auto cP = thrust::make_counting_iterator<int>(0);
        auto cV = thrust::make_counting_iterator<int>(0);

        // (1) 边心
        thrust::for_each(exec, cE, cE + int(ne),
            EdgeCentroidOp(
                thrust::raw_pointer_cast(dEV.data()),
                thrust::raw_pointer_cast(dV.data()),
                thrust::raw_pointer_cast(d_EC.data())
            )
        );

        // (2) 面心
        thrust::for_each(exec, cF, cF + int(nf),
            FaceCentroidOp(
                thrust::raw_pointer_cast(d_faceEdges.data()),
                thrust::raw_pointer_cast(d_EC.data()),
                thrust::raw_pointer_cast(d_FC.data())
            )
        );

        // (3) 体心
        thrust::for_each(exec, cP, cP + int(np),
            PolyCentroidOp(
                thrust::raw_pointer_cast(d_polyFaces.data()),
                thrust::raw_pointer_cast(d_FC.data()),
                thrust::raw_pointer_cast(d_PC.data())
            )
        );

        // (4) 新体点：直接取前 patchPolys 个体心
        thrust::copy_n(exec, d_PC.begin(), patchPolys, d_newPoly.begin());

        // ====== 准备边界/内部索引，减少大核中的分支 ======
        // faces
        thrust::device_vector<int> d_allFaces(patchFaces);
        thrust::sequence(exec, d_allFaces.begin(), d_allFaces.end(), 0);
        thrust::device_vector<int> d_facesInterior(patchFaces), d_facesBoundary(patchFaces);
        auto it_f_int_end = thrust::copy_if(exec, d_allFaces.begin(), d_allFaces.end(), d_facesInterior.begin(),
            IsInteriorPred(thrust::raw_pointer_cast(d_fOn.data())));
        auto it_f_bnd_end = thrust::copy_if(exec, d_allFaces.begin(), d_allFaces.end(), d_facesBoundary.begin(),
            IsBoundaryPred(thrust::raw_pointer_cast(d_fOn.data())));
        d_facesInterior.resize(it_f_int_end - d_facesInterior.begin());
        d_facesBoundary.resize(it_f_bnd_end - d_facesBoundary.begin());

        // edges
        thrust::device_vector<int> d_allEdges(patchEdges);
        thrust::sequence(exec, d_allEdges.begin(), d_allEdges.end(), 0);
        thrust::device_vector<int> d_edgesInterior(patchEdges), d_edgesBoundary(patchEdges);
        auto it_e_int_end = thrust::copy_if(exec, d_allEdges.begin(), d_allEdges.end(), d_edgesInterior.begin(),
            IsInteriorPred(thrust::raw_pointer_cast(d_eOn.data())));
        auto it_e_bnd_end = thrust::copy_if(exec, d_allEdges.begin(), d_allEdges.end(), d_edgesBoundary.begin(),
            IsBoundaryPred(thrust::raw_pointer_cast(d_eOn.data())));
        d_edgesInterior.resize(it_e_int_end - d_edgesInterior.begin());
        d_edgesBoundary.resize(it_e_bnd_end - d_edgesBoundary.begin());

        // verts
        thrust::device_vector<int> d_allVerts(patchVerts);
        thrust::sequence(exec, d_allVerts.begin(), d_allVerts.end(), 0);
        thrust::device_vector<int> d_vertsInterior(patchVerts), d_vertsBoundary(patchVerts);
        auto it_v_int_end = thrust::copy_if(exec, d_allVerts.begin(), d_allVerts.end(), d_vertsInterior.begin(),
            IsInteriorPred(thrust::raw_pointer_cast(d_vOn.data())));
        auto it_v_bnd_end = thrust::copy_if(exec, d_allVerts.begin(), d_allVerts.end(), d_vertsBoundary.begin(),
            IsBoundaryPred(thrust::raw_pointer_cast(d_vOn.data())));
        d_vertsInterior.resize(it_v_int_end - d_vertsInterior.begin());
        d_vertsBoundary.resize(it_v_bnd_end - d_vertsBoundary.begin());

        // (5) 新面点：内部/边界分开跑
        {
            // 内部
            thrust::for_each(exec, d_facesInterior.begin(), d_facesInterior.end(),
                NewFaceVert_Interior_Op(
                    thrust::raw_pointer_cast(d_fpOff.data()),
                    thrust::raw_pointer_cast(d_facePolys.data()),
                    thrust::raw_pointer_cast(d_FC.data()),
                    thrust::raw_pointer_cast(d_PC.data()),
                    thrust::raw_pointer_cast(d_newFace.data())
                )
            );
            // 边界
            thrust::for_each(exec, d_facesBoundary.begin(), d_facesBoundary.end(),
                NewFaceVert_Boundary_Op(
                    thrust::raw_pointer_cast(d_FC.data()),
                    thrust::raw_pointer_cast(d_newFace.data())
                )
            );
        }

        // (6) 新边点：内部/边界分开跑
        {
            // 内部
            thrust::for_each(exec, d_edgesInterior.begin(), d_edgesInterior.end(),
                NewEdgeVert_Interior_Op(
                    thrust::raw_pointer_cast(d_efOff.data()),
                    thrust::raw_pointer_cast(d_edgeFaces.data()),
                    thrust::raw_pointer_cast(d_fpOff.data()),
                    thrust::raw_pointer_cast(d_facePolys.data()),
                    thrust::raw_pointer_cast(d_EC.data()),
                    thrust::raw_pointer_cast(d_FC.data()),
                    thrust::raw_pointer_cast(d_PC.data()),
                    thrust::raw_pointer_cast(d_newEdge.data())
                )
            );
            // 边界
            thrust::for_each(exec, d_edgesBoundary.begin(), d_edgesBoundary.end(),
                NewEdgeVert_Boundary_Op(
                    thrust::raw_pointer_cast(d_fOn.data()),
                    thrust::raw_pointer_cast(d_efOff.data()),
                    thrust::raw_pointer_cast(d_edgeFaces.data()),
                    thrust::raw_pointer_cast(d_EC.data()),
                    thrust::raw_pointer_cast(d_FC.data()),
                    thrust::raw_pointer_cast(d_newEdge.data())
                )
            );
        }

        // (7) 新顶点：内部/边界分开跑
        {
            // 内部
            thrust::for_each(exec, d_vertsInterior.begin(), d_vertsInterior.end(),
                NewVertVert_Interior_Op(
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
                    thrust::raw_pointer_cast(d_newVert.data())
                )
            );
            // 边界
            thrust::for_each(exec, d_vertsBoundary.begin(), d_vertsBoundary.end(),
                NewVertVert_Boundary_Op(
                    thrust::raw_pointer_cast(d_eOn.data()),
                    thrust::raw_pointer_cast(d_fOn.data()),
                    thrust::raw_pointer_cast(d_veOff.data()),
                    thrust::raw_pointer_cast(d_vertEdges.data()),
                    thrust::raw_pointer_cast(d_efOff.data()),
                    thrust::raw_pointer_cast(d_edgeFaces.data()),
                    thrust::raw_pointer_cast(dV.data()),
                    thrust::raw_pointer_cast(d_EC.data()),
                    thrust::raw_pointer_cast(d_FC.data()),
                    thrust::raw_pointer_cast(d_newVert.data())
                )
            );
        }

        // ----------- 把新点拷回主机并写入 pos -----------
        thrust::host_vector<dvec3a> h_newPoly = d_newPoly;
        thrust::host_vector<dvec3a> h_newFace = d_newFace;
        thrust::host_vector<dvec3a> h_newEdge = d_newEdge;
        thrust::host_vector<dvec3a> h_newVert = d_newVert;

        // 计算偏移（插入前的 pos.size()）
        const uint pvOffset = static_cast<uint>(pos.size());
        const uint fvOffset = pvOffset + static_cast<uint>(h_newPoly.size());
        const uint evOffset = fvOffset + static_cast<uint>(h_newFace.size());
        const uint vvOffset = evOffset + static_cast<uint>(h_newEdge.size());

        // 追加新点
        pos.reserve(pos.size() + h_newPoly.size() + h_newFace.size() + h_newEdge.size() + h_newVert.size());
        for (const auto& v : h_newPoly) pos.push_back(to_h(v));
        for (const auto& v : h_newFace) pos.push_back(to_h(v));
        for (const auto& v : h_newEdge) pos.push_back(to_h(v));
        for (const auto& v : h_newVert) pos.push_back(to_h(v));

        // 设备侧输出（每 poly 64 个 uint）
        thrust::device_vector<uint> d_topo(patchPolys * 64u);

        // 并行装配（每 poly 一个线程）——复用 dEV
        auto cP2 = thrust::make_counting_iterator<int>(0);
        thrust::for_each(exec, cP2, cP2 + int(patchPolys),
            TopoAssembleOp(
                thrust::raw_pointer_cast(d_polyFaces.data()),
                thrust::raw_pointer_cast(d_faceEdges.data()),
                thrust::raw_pointer_cast(dEV.data()),
                pvOffset, fvOffset, evOffset, vvOffset,
                thrust::raw_pointer_cast(d_topo.data())
            )
        );

        // 回拷并追加到 polys（顺序稳定：p=0..patchPolys-1）
        thrust::host_vector<uint> h_topo = d_topo;
        polys.reserve(polys.size() + h_topo.size());
        polys.insert(polys.end(), h_topo.begin(), h_topo.end());

        cudaStreamDestroy(stream);
    }

    // ===== 转换函数 =====
    static inline dvec3a to_d(const vec3d& v) {
        return dvec3a(v.x(), v.y(), v.z());
    }
    static inline vec3d to_h(const dvec3a& v) {
        return vec3d(v.x, v.y, v.z);
    }

} // namespace cinolib
