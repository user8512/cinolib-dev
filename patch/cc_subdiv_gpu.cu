#define CINO_STATIC_LIB
#include "patch.h"
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/transform.h>
#include <thrust/for_each.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/functional.h>
#include <cuda_runtime.h>

#ifndef CUDA_CHECK
#define CUDA_CHECK(call) do { \
    cudaError_t _e = (call);  \
    if (_e != cudaSuccess) {  \
        fprintf(stderr, "CUDA error %s:%d: %s\n", __FILE__, __LINE__, cudaGetErrorString(_e)); \
        abort(); \
    } \
} while(0)
#endif

namespace cinolib {

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

    __host__ __device__ __forceinline__ int idx_abs(int v) { return (v < 0) ? (-v - 1) : v; }

    static inline dvec3 to_d(const cinolib::vec3d& v);
    static inline vec3d to_h(const dvec3& v);

    // -------------------- Thrust 侧算子（保持不变） --------------------

    struct EdgeCentroidOp {
        const uint2* __restrict__ EV;
        const dvec3* __restrict__ Vpos;
        dvec3* __restrict__ Ecentroids;
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

    struct FaceCentroidOp {
        const int* __restrict__ faceEdges;
        const dvec3* __restrict__ Ecentroids;
        dvec3* __restrict__ Fcentroids;
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

    struct PolyCentroidOp {
        const int* __restrict__ polyFaces;
        const dvec3* __restrict__ Fcentroids;
        dvec3* __restrict__ Pcentroids;
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

    struct NewFaceVertOp {
        const uint8_t* __restrict__ faceOnSurf;
        const uint* __restrict__ facePolysOffset;
        const int* __restrict__ facePolys;
        const dvec3* __restrict__ Fcentroids;
        const dvec3* __restrict__ Pcentroids;
        dvec3* __restrict__ outNewF;
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

    struct TopoAssembleOp {
        const int* __restrict__ polyFaces;
        const int* __restrict__ faceEdges;
        const uint2* __restrict__ edgeVerts;
        const uint    pvOffset, fvOffset, evOffset, vvOffset;
        uint* __restrict__ out;

        __host__ __device__
            TopoAssembleOp(const int* pf, const int* fe, const uint2* ev,
                uint pv, uint fv, uint evv, uint vv, uint* o)
            : polyFaces(pf), faceEdges(fe), edgeVerts(ev),
            pvOffset(pv), fvOffset(fv), evOffset(evv), vvOffset(vv), out(o) {}

        __host__ __device__ static inline int iabs(int v) { return v < 0 ? (-v - 1) : v; }
        __host__ __device__ static inline bool contains_u(const uint* arr, int n, uint v) {
            for (int i = 0; i < n; ++i) if (arr[i] == v) return true;
            return false;
        }
        __host__ __device__ static inline bool contains_i(const int* arr, int n, int v) {
            for (int i = 0; i < n; ++i) if (arr[i] == v) return true;
            return false;
        }

        __host__ __device__
            void operator()(const int p) const {
            uint VV[8];  int nVV = 0;
            uint EV[12]; int nEV = 0;
            uint FV[6];  int nFV = 0;

            int FACES[6];
#pragma unroll
            for (int fo = 0; fo < 6; ++fo) {
                FACES[fo] = iabs(polyFaces[p * 6 + fo]);
            }

            int raw_f1 = polyFaces[p * 6 + 0];
            bool f1Reverse = (raw_f1 < 0);
            int f1 = iabs(raw_f1);

            int F1E[4];
#pragma unroll
            for (int eOff = 0; eOff < 4; ++eOff) {
                F1E[eOff] = iabs(faceEdges[f1 * 4 + eOff]);
            }

            int f2 = -1;
            for (int fo = 1; fo < 6; ++fo) {
                int f = FACES[fo];
                bool share = false;
#pragma unroll
                for (int eOff = 0; eOff < 4; ++eOff) {
                    int e = iabs(faceEdges[f * 4 + eOff]);
                    for (int k = 0; k < 4; ++k) { if (e == F1E[k]) { share = true; break; } }
                    if (share) break;
                }
                if (!share) { f2 = f; break; }
            }

            uint EDGES[12]; int nEDGES = 0;
            for (int fo = 0; fo < 6; ++fo) {
                int f = FACES[fo];
#pragma unroll
                for (int eOff = 0; eOff < 4; ++eOff) {
                    uint e = (uint)iabs(faceEdges[f * 4 + eOff]);
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
                    uint2 ev = edgeVerts[(uint)eTmp];
                    if ((f1Reverse ^ edgeReverse) && (int)ev.y == vCurrent) {
                        vCurrent = ev.x; EV[nEV++] = (uint)eTmp; break;
                    }
                    else if (!(f1Reverse ^ edgeReverse) && (int)ev.x == vCurrent) {
                        vCurrent = ev.y; EV[nEV++] = (uint)eTmp; break;
                    }
                }
            }

            for (int i = 0; i < 4; ++i) {
                int vtx = (int)VV[i];
                for (int t = 0; t < nEDGES; ++t) {
                    uint e = EDGES[t];
                    if (!contains_u(EV, nEV, e)) {
                        uint2 ev = edgeVerts[e];
                        if ((int)ev.x == vtx) { VV[nVV++] = ev.y; EV[nEV++] = e; break; }
                        else if ((int)ev.y == vtx) { VV[nVV++] = ev.x; EV[nEV++] = e; break; }
                    }
                }
            }

            FV[nFV++] = (uint)f1;
            FV[nFV++] = (uint)f2;

            for (int i = 0; i < 4; ++i) {
                int eNeed = (int)EV[i];
                for (int fo = 0; fo < 6; ++fo) {
                    uint f = (uint)FACES[fo];
                    if (!contains_u(FV, nFV, f)) {
                        bool hit = false;
                        for (int eOff = 0; eOff < 4; ++eOff) {
                            int fe = iabs(faceEdges[f * 4 + eOff]);
                            if (fe == eNeed) { hit = true; break; }
                        }
                        if (hit) {
                            FV[nFV++] = f;
                            for (int eOff = 0; eOff < 4; ++eOff) {
                                int fe = iabs(faceEdges[f * 4 + eOff]);
                                if (!contains_u(EV, nEV, (uint)fe)) {
                                    EV[nEV++] = (uint)fe;
                                    break;
                                }
                            }
                        }
                    }
                }
            }

            uint PV = (uint)p + pvOffset;
            for (int i = 0; i < 6; ++i)  FV[i] += fvOffset;
            for (int i = 0; i < 12; ++i)  EV[i] += evOffset;
            for (int i = 0; i < 8; ++i)  VV[i] += vvOffset;

            uint base = (uint)p * 64u;
            out[base + 0] = VV[0]; out[base + 1] = EV[3]; out[base + 2] = FV[0]; out[base + 3] = EV[0];
            out[base + 4] = EV[4]; out[base + 5] = FV[5]; out[base + 6] = PV;    out[base + 7] = FV[2];
            out[base + 8] = EV[0]; out[base + 9] = FV[0]; out[base + 10] = EV[1]; out[base + 11] = VV[1];
            out[base + 12] = FV[2]; out[base + 13] = PV; out[base + 14] = FV[3]; out[base + 15] = EV[5];
            out[base + 16] = FV[0]; out[base + 17] = EV[2]; out[base + 18] = VV[2]; out[base + 19] = EV[1];
            out[base + 20] = PV;    out[base + 21] = FV[4]; out[base + 22] = EV[6]; out[base + 23] = FV[3];
            out[base + 24] = EV[3]; out[base + 25] = VV[3]; out[base + 26] = EV[2]; out[base + 27] = FV[0];
            out[base + 28] = FV[5]; out[base + 29] = EV[7];    out[base + 30] = FV[4]; out[base + 31] = PV;
            out[base + 32] = EV[4]; out[base + 33] = FV[5]; out[base + 34] = PV;    out[base + 35] = FV[2];
            out[base + 36] = VV[4]; out[base + 37] = EV[11]; out[base + 38] = FV[1]; out[base + 39] = EV[8];
            out[base + 40] = FV[2]; out[base + 41] = PV; out[base + 42] = FV[3]; out[base + 43] = EV[5];
            out[base + 44] = EV[8]; out[base + 45] = FV[1]; out[base + 46] = EV[9]; out[base + 47] = VV[5];
            out[base + 48] = PV;    out[base + 49] = FV[4]; out[base + 50] = EV[6]; out[base + 51] = FV[3];
            out[base + 52] = FV[1]; out[base + 53] = EV[10]; out[base + 54] = VV[6]; out[base + 55] = EV[9];
            out[base + 56] = FV[5]; out[base + 57] = EV[7];    out[base + 58] = FV[4]; out[base + 59] = PV;
            out[base + 60] = EV[11]; out[base + 61] = VV[7]; out[base + 62] = EV[10]; out[base + 63] = FV[1];
        }
    };

    // -------------------- 新增：共享内存 Kernel (方案A) --------------------

    template<int BLOCK_SIZE, int MAX_DEG = 32>
    __global__ void NewEdgeVertKernel(
        const uint8_t* __restrict__ edgeOnSurf,
        const uint8_t* __restrict__ faceOnSurf,
        const uint* __restrict__ edgeFacesOffset,
        const int* __restrict__ edgeFaces,
        const uint* __restrict__ facePolysOffset,
        const int* __restrict__ facePolys,
        const dvec3* __restrict__ Ecentroids,
        const dvec3* __restrict__ Fcentroids,
        const dvec3* __restrict__ Pcentroids,
        dvec3* __restrict__ outNewE,
        int patchEdges)
    {
        extern __shared__ int smem[];
        // 每线程两段：faces[MAX_DEG], polys[MAX_DEG]
        int* faces = smem + threadIdx.x * MAX_DEG;
        int* polys = smem + BLOCK_SIZE * MAX_DEG + threadIdx.x * MAX_DEG;

        const int e = blockIdx.x * blockDim.x + threadIdx.x;
        if (e >= patchEdges) return;

        const uint begin = edgeFacesOffset[e];
        const uint end = edgeFacesOffset[e + 1];
        const int  N = int(end - begin);

        if (!edgeOnSurf[e]) {
            // 收集相邻面（最多 MAX_DEG）
            int nf = 0;
#pragma unroll
            for (int i = 0; i < MAX_DEG; ++i) {
                int idx = (i < N) ? idx_abs(edgeFaces[begin + i]) : -1;
                faces[i] = idx;
                if (idx >= 0) ++nf;
            }
            // 从相邻面收集体并去重（线性去重；MAX_DEG 很小）
            int np = 0;
            for (int i = 0; i < nf; ++i) {
                int f = faces[i];
                uint fb = facePolysOffset[f];
                uint fe = facePolysOffset[f + 1];
                for (uint k = fb; k < fe; ++k) {
                    int p = idx_abs(facePolys[k]);
                    bool seen = false;
#pragma unroll
                    for (int t = 0; t < np; ++t) { if (polys[t] == p) { seen = true; break; } }
                    if (!seen && np < MAX_DEG) polys[np++] = p;
                }
            }
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
                // 只统计边界面
                int nf = 0;
                for (int i = 0; i < N && nf < MAX_DEG; ++i) {
                    int f = idx_abs(edgeFaces[begin + i]);
                    if (faceOnSurf[f]) faces[nf++] = f;
                }
                dvec3 faceAvg(0, 0, 0);
                for (int i = 0; i < nf; ++i) faceAvg += Fcentroids[faces[i]];
                if (nf > 0) faceAvg /= double(nf);
                outNewE[e] = (faceAvg + Ecentroids[e]) * 0.5;
            }
        }
    }

    template<int BLOCK_SIZE, int MAX_DEG = 32>
    __global__ void NewVertVertKernel(
        const uint8_t* __restrict__ vertOnSurf,
        const uint8_t* __restrict__ edgeOnSurf,
        const uint8_t* __restrict__ faceOnSurf,
        const uint* __restrict__ vertEdgesOffset,
        const int* __restrict__ vertEdges,
        const uint* __restrict__ edgeFacesOffset,
        const int* __restrict__ edgeFaces,
        const uint* __restrict__ facePolysOffset,
        const int* __restrict__ facePolys,
        const dvec3* __restrict__ Vpos,
        const dvec3* __restrict__ Ecentroids,
        const dvec3* __restrict__ Fcentroids,
        const dvec3* __restrict__ Pcentroids,
        dvec3* __restrict__ outNewV,
        int patchVerts)
    {
        extern __shared__ int smem[];
        // 每线程三段：edges[MAX_DEG], faces[MAX_DEG], polys[MAX_DEG]
        int* edges = smem + threadIdx.x * MAX_DEG;
        int* faces = smem + BLOCK_SIZE * MAX_DEG + threadIdx.x * MAX_DEG;
        int* polys = smem + 2 * BLOCK_SIZE * MAX_DEG + threadIdx.x * MAX_DEG;

        const int v = blockIdx.x * blockDim.x + threadIdx.x;
        if (v >= patchVerts) return;

        const uint begin = vertEdgesOffset[v];
        const uint end = vertEdgesOffset[v + 1];
        const int  N = int(end - begin);

        if (!vertOnSurf[v]) {
            // 收集边
            int ne = 0;
            for (int i = 0; i < N && ne < MAX_DEG; ++i) {
                int e = idx_abs(vertEdges[begin + i]);
                edges[ne++] = e;
            }
            // 通过边→面（去重）
            int nf = 0;
            for (int i = 0; i < ne; ++i) {
                int e = edges[i];
                uint eb = edgeFacesOffset[e];
                uint ee = edgeFacesOffset[e + 1];
                for (uint j = eb; j < ee; ++j) {
                    int f = idx_abs(edgeFaces[j]);
                    bool seen = false;
#pragma unroll
                    for (int t = 0; t < nf; ++t) { if (faces[t] == f) { seen = true; break; } }
                    if (!seen && nf < MAX_DEG) faces[nf++] = f;
                }
            }
            // 通过面→体（去重）
            int np = 0;
            for (int i = 0; i < nf; ++i) {
                int f = faces[i];
                uint fb = facePolysOffset[f];
                uint fe = facePolysOffset[f + 1];
                for (uint k = fb; k < fe; ++k) {
                    int p = idx_abs(facePolys[k]);
                    bool seenp = false;
#pragma unroll
                    for (int t = 0; t < np; ++t) { if (polys[t] == p) { seenp = true; break; } }
                    if (!seenp && np < MAX_DEG) polys[np++] = p;
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
            int ne = 0;
            for (int i = 0; i < N && ne < MAX_DEG; ++i) {
                int e = idx_abs(vertEdges[begin + i]);
                if (edgeOnSurf[e]) edges[ne++] = e;
            }
            int nf = 0;
            for (int i = 0; i < ne; ++i) {
                int e = edges[i];
                uint eb = edgeFacesOffset[e];
                uint ee = edgeFacesOffset[e + 1];
                for (uint j = eb; j < ee; ++j) {
                    int f = idx_abs(edgeFaces[j]);
                    if (faceOnSurf[f]) {
                        bool seen = false;
#pragma unroll
                        for (int t = 0; t < nf; ++t) { if (faces[t] == f) { seen = true; break; } }
                        if (!seen && nf < MAX_DEG) faces[nf++] = f;
                    }
                }
            }

            dvec3 edgeAvg(0, 0, 0);
            for (int i = 0; i < ne; ++i) edgeAvg += Ecentroids[edges[i]];
            if (ne > 0) edgeAvg /= double(ne);

            dvec3 faceAvg(0, 0, 0);
            for (int i = 0; i < nf; ++i) faceAvg += Fcentroids[faces[i]];
            if (nf > 0) faceAvg /= double(nf);

            int n = ne; // 只考虑边界边的度
            dvec3 vnew = faceAvg + (edgeAvg * 2.0) + (Vpos[v] * double(n - 3));
            if (n > 0) vnew /= double(n);
            outNewV[v] = vnew;
        }
    }

    // -------------------- 主流程 --------------------

    void Patch::singlePatch::subdiv_cuda(std::vector<vec3d>& pos, std::vector<uint>& polys)
    {
        const uint ne = static_cast<uint>(edgeVerts.size());
        const uint nf = static_cast<uint>(faceEdges.size() / 4);
        const uint np = static_cast<uint>(polyFaces.size() / 6);
        const uint nv = static_cast<uint>(vertsPos.size());

        // ---- Host → Device ----
        thrust::host_vector<dvec3> hV(nv);
        for (uint i = 0; i < nv; ++i) hV[i] = to_d(vertsPos[i]);

        thrust::host_vector<uint2> hEV(ne);
        for (uint e = 0; e < ne; ++e) {
            auto u = edgeVerts[e].x();
            auto v = edgeVerts[e].y();
            hEV[e] = make_uint2(u, v);
        }

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

        auto pack_bool = [](const std::vector<bool>& v) {
            thrust::host_vector<uint8_t> out(v.size());
            for (size_t i = 0; i < v.size(); ++i) out[i] = v[i] ? 1u : 0u;
            return out;
        };
        thrust::device_vector<uint8_t> d_vOn = pack_bool(vertOnSurf);
        thrust::device_vector<uint8_t> d_eOn = pack_bool(edgeOnSurf);
        thrust::device_vector<uint8_t> d_fOn = pack_bool(faceOnSurf);

        // ---- 设备侧输出缓冲 ----
        thrust::device_vector<dvec3> d_EC(ne);
        thrust::device_vector<dvec3> d_FC(nf);
        thrust::device_vector<dvec3> d_PC(np);

        thrust::device_vector<dvec3> d_newPoly(patchPolys);
        thrust::device_vector<dvec3> d_newFace(patchFaces);
        thrust::device_vector<dvec3> d_newEdge(patchEdges);
        thrust::device_vector<dvec3> d_newVert(patchVerts);

        auto c_begin_e = thrust::make_counting_iterator<int>(0);
        auto c_begin_f = thrust::make_counting_iterator<int>(0);
        auto c_begin_p = thrust::make_counting_iterator<int>(0);
        auto c_begin_v = thrust::make_counting_iterator<int>(0);

        
        auto cuda_start = std::chrono::high_resolution_clock::now();

        // (1) 边心
        thrust::for_each(c_begin_e, c_begin_e + int(ne),
            EdgeCentroidOp(
                thrust::raw_pointer_cast(dEV.data()),
                thrust::raw_pointer_cast(dV.data()),
                thrust::raw_pointer_cast(d_EC.data())));

        // (2) 面心
        thrust::for_each(c_begin_f, c_begin_f + int(nf),
            FaceCentroidOp(
                thrust::raw_pointer_cast(d_faceEdges.data()),
                thrust::raw_pointer_cast(d_EC.data()),
                thrust::raw_pointer_cast(d_FC.data())));

        // (3) 体心
        thrust::for_each(c_begin_p, c_begin_p + int(np),
            PolyCentroidOp(
                thrust::raw_pointer_cast(d_polyFaces.data()),
                thrust::raw_pointer_cast(d_FC.data()),
                thrust::raw_pointer_cast(d_PC.data())));

        // (4) 新体点 = 前 patchPolys 个体心
        thrust::copy(d_PC.begin(), d_PC.begin() + patchPolys, d_newPoly.begin());

        // (5) 新面点
        thrust::for_each(c_begin_f, c_begin_f + int(patchFaces),
            NewFaceVertOp(
                thrust::raw_pointer_cast(d_fOn.data()),
                thrust::raw_pointer_cast(d_fpOff.data()),
                thrust::raw_pointer_cast(d_facePolys.data()),
                thrust::raw_pointer_cast(d_FC.data()),
                thrust::raw_pointer_cast(d_PC.data()),
                thrust::raw_pointer_cast(d_newFace.data())));

        // (6) 新边点 —— 共享内存 Kernel（替换原 Thrust for_each）
        {
            constexpr int BLOCK = 128;
            constexpr int MAX_DEG = 32; // 若邻接上限更大，可调高并留意共享内存占用
            int grid = (int(patchEdges) + BLOCK - 1) / BLOCK;
            // faces[MAX_DEG] + polys[MAX_DEG] per thread
            size_t shmem = BLOCK * (MAX_DEG + MAX_DEG) * sizeof(int);
            NewEdgeVertKernel<BLOCK, MAX_DEG> << <grid, BLOCK, shmem >> > (
                thrust::raw_pointer_cast(d_eOn.data()),
                thrust::raw_pointer_cast(d_fOn.data()),
                thrust::raw_pointer_cast(d_efOff.data()),
                thrust::raw_pointer_cast(d_edgeFaces.data()),
                thrust::raw_pointer_cast(d_fpOff.data()),
                thrust::raw_pointer_cast(d_facePolys.data()),
                thrust::raw_pointer_cast(d_EC.data()),
                thrust::raw_pointer_cast(d_FC.data()),
                thrust::raw_pointer_cast(d_PC.data()),
                thrust::raw_pointer_cast(d_newEdge.data()),
                int(patchEdges)
                );
            CUDA_CHECK(cudaGetLastError());
        }

        // (7) 新顶点 —— 共享内存 Kernel（替换原 Thrust for_each）
        {
            constexpr int BLOCK = 128;
            constexpr int MAX_DEG = 32;
            int grid = (int(patchVerts) + BLOCK - 1) / BLOCK;
            // edges[MAX_DEG] + faces[MAX_DEG] + polys[MAX_DEG] per thread
            size_t shmem = BLOCK * (3 * MAX_DEG) * sizeof(int);
            NewVertVertKernel<BLOCK, MAX_DEG> << <grid, BLOCK, shmem >> > (
                thrust::raw_pointer_cast(d_vOn.data()),
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
                thrust::raw_pointer_cast(d_newVert.data()),
                int(patchVerts)
                );
            CUDA_CHECK(cudaGetLastError());
            CUDA_CHECK(cudaDeviceSynchronize());
        }

        auto cuda_end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> elapsed = cuda_end - cuda_start;
        double patch_elapsed = elapsed.count();

        // ---- 回拷到主机 ----
        thrust::host_vector<dvec3> h_newPoly = d_newPoly;
        thrust::host_vector<dvec3> h_newFace = d_newFace;
        thrust::host_vector<dvec3> h_newEdge = d_newEdge;
        thrust::host_vector<dvec3> h_newVert = d_newVert;

        uint pvOffset = static_cast<uint>(pos.size());
        uint fvOffset = pvOffset + static_cast<uint>(h_newPoly.size());
        uint evOffset = fvOffset + static_cast<uint>(h_newFace.size());
        uint vvOffset = evOffset + static_cast<uint>(h_newEdge.size());

        pos.reserve(pos.size() + h_newPoly.size() + h_newFace.size() + h_newEdge.size() + h_newVert.size());
        for (const auto& v : h_newPoly) pos.push_back(to_h(v));
        for (const auto& v : h_newFace) pos.push_back(to_h(v));
        for (const auto& v : h_newEdge) pos.push_back(to_h(v));
        for (const auto& v : h_newVert) pos.push_back(to_h(v));

        thrust::device_vector<uint> d_topo(patchPolys * 64u);

        thrust::host_vector<uint2> hEV2(edgeVerts.size());
        for (size_t e = 0; e < edgeVerts.size(); ++e) {
            hEV2[e] = make_uint2(edgeVerts[e].x(), edgeVerts[e].y());
        }
        thrust::device_vector<uint2> dEV2 = hEV2;

        auto c_begin_p2 = thrust::make_counting_iterator<int>(0);

        cuda_start = std::chrono::high_resolution_clock::now();

        thrust::for_each(c_begin_p2, c_begin_p2 + int(patchPolys),
            TopoAssembleOp(
                thrust::raw_pointer_cast(d_polyFaces.data()),
                thrust::raw_pointer_cast(d_faceEdges.data()),
                thrust::raw_pointer_cast(dEV2.data()),
                pvOffset, fvOffset, evOffset, vvOffset,
                thrust::raw_pointer_cast(d_topo.data())
            )
        );

        cuda_end = std::chrono::high_resolution_clock::now();
        elapsed = cuda_end - cuda_start;
        patch_elapsed += elapsed.count();
        std::cout << "execution time of current patch: " << patch_elapsed << " ms" << std::endl;
        cuda_elapsed += patch_elapsed;

        thrust::host_vector<uint> h_topo = d_topo;
        polys.reserve(polys.size() + h_topo.size());
        polys.insert(polys.end(), h_topo.begin(), h_topo.end());
    }

    // ---- host/device vec 转换 ----
    static inline dvec3 to_d(const vec3d& v) {
        return dvec3(v.x(), v.y(), v.z());
    }
    static inline vec3d to_h(const dvec3& v) {
        return vec3d(v.x, v.y, v.z);
    }
} // namespace cinolib