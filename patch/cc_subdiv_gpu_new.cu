// cc_subdiv_gpu_new.cu  —— 纯 kernel + 显式 CUDA Graph 版本
#define CINO_STATIC_LIB
#include "patch.h"

#include <cuda_runtime.h>
#include <cuda.h>
#include <vector>
#include <cassert>
#include <chrono>
#include <iostream>

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

    // ==================== 基础类型/工具 ====================

    struct alignas(16) dvec3 {
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

    static inline dvec3 to_d(const vec3d& v) { return dvec3(v.x(), v.y(), v.z()); }
    static inline vec3d to_h(const dvec3& v) { return vec3d(v.x, v.y, v.z); }

    // ==================== 基础 Kernels（替代 Thrust 算子） ====================

    __global__ void EdgeCentroidKernel(
        const uint2* __restrict__ EV,
        const dvec3* __restrict__ Vpos,
        dvec3* __restrict__ Ecentroids,
        int ne)
    {
        int e = blockIdx.x * blockDim.x + threadIdx.x;
        if (e >= ne) return;
        uint2 ev = EV[e];
        dvec3 a = Vpos[ev.x];
        dvec3 b = Vpos[ev.y];
        Ecentroids[e] = (a + b) * 0.5;
    }

    __global__ void FaceCentroidKernel(
        const int* __restrict__ faceEdges, // 4 per face
        const dvec3* __restrict__ Ecentroids,
        dvec3* __restrict__ Fcentroids,
        int nf)
    {
        int f = blockIdx.x * blockDim.x + threadIdx.x;
        if (f >= nf) return;
        dvec3 c(0, 0, 0);
#pragma unroll
        for (int k = 0; k < 4; ++k) {
            int e = idx_abs(faceEdges[f * 4 + k]);
            c += Ecentroids[e];
        }
        Fcentroids[f] = c / 4.0;
    }

    __global__ void PolyCentroidKernel(
        const int* __restrict__ polyFaces, // 6 per poly
        const dvec3* __restrict__ Fcentroids,
        dvec3* __restrict__ Pcentroids,
        int np)
    {
        int p = blockIdx.x * blockDim.x + threadIdx.x;
        if (p >= np) return;
        dvec3 c(0, 0, 0);
#pragma unroll
        for (int k = 0; k < 6; ++k) {
            int f = idx_abs(polyFaces[p * 6 + k]);
            c += Fcentroids[f];
        }
        Pcentroids[p] = c / 6.0;
    }

    __global__ void NewFaceVertKernel(
        const uint8_t* __restrict__ faceOnSurf,
        const uint* __restrict__ facePolysOffset,
        const int* __restrict__ facePolys,
        const dvec3* __restrict__ Fcentroids,
        const dvec3* __restrict__ Pcentroids,
        dvec3* __restrict__ outNewF,
        int patchFaces)
    {
        int f = blockIdx.x * blockDim.x + threadIdx.x;
        if (f >= patchFaces) return;
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

    __global__ void CopyFirstNKernel(
        const dvec3* __restrict__ src,
        dvec3* __restrict__ dst,
        int n)
    {
        int i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i < n) dst[i] = src[i];
    }

    // ==================== 共享内存 Kernels（保留） ====================

    template<int BLOCK_SIZE, int MAX_DEG = 32>
    __launch_bounds__(BLOCK_SIZE, 2)
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
        int* faces = smem + threadIdx.x * MAX_DEG;
        int* polys = smem + BLOCK_SIZE * MAX_DEG + threadIdx.x * MAX_DEG;

        const int e = blockIdx.x * blockDim.x + threadIdx.x;
        if (e >= patchEdges) return;

        const uint begin = edgeFacesOffset[e];
        const uint end = edgeFacesOffset[e + 1];
        const int  N = int(end - begin);

        if (!edgeOnSurf[e]) {
            int nf = 0;
#pragma unroll
            for (int i = 0; i < MAX_DEG; ++i) {
                int idx = (i < N) ? idx_abs(edgeFaces[begin + i]) : -1;
                faces[i] = idx;
                if (idx >= 0) ++nf;
            }
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
    __launch_bounds__(BLOCK_SIZE, 2)
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
        int* edges = smem + threadIdx.x * MAX_DEG;
        int* faces = smem + BLOCK_SIZE * MAX_DEG + threadIdx.x * MAX_DEG;
        int* polys = smem + 2 * BLOCK_SIZE * MAX_DEG + threadIdx.x * MAX_DEG;

        const int v = blockIdx.x * blockDim.x + threadIdx.x;
        if (v >= patchVerts) return;

        const uint begin = vertEdgesOffset[v];
        const uint end = vertEdgesOffset[v + 1];
        const int  N = int(end - begin);

        if (!vertOnSurf[v]) {
            int ne = 0;
            for (int i = 0; i < N && ne < MAX_DEG; ++i) edges[ne++] = idx_abs(vertEdges[begin + i]);

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

            int n = ne;
            dvec3 vnew = faceAvg + (edgeAvg * 2.0) + (Vpos[v] * double(n - 3));
            if (n > 0) vnew /= double(n);
            outNewV[v] = vnew;
        }
    }

    // ==================== 拓扑装配 Kernel ====================

    __device__ __forceinline__ bool contains_u_u32(const uint* arr, int n, uint v) {
        for (int i = 0; i < n; ++i) if (arr[i] == v) return true;
        return false;
    }

    __global__ void TopoAssembleKernel(
        const int* __restrict__ polyFaces,   // 6 per poly
        const int* __restrict__ faceEdges,   // 4 per face
        const uint2* __restrict__ edgeVerts,
        uint pvOffset, uint fvOffset, uint evOffset, uint vvOffset,
        uint* __restrict__ out, // 64 per poly
        int patchPolys)
    {
        int p = blockIdx.x * blockDim.x + threadIdx.x;
        if (p >= patchPolys) return;

        uint VV[8];  int nVV = 0;
        uint EV[12]; int nEV = 0;
        uint FV[6];  int nFV = 0;

        int FACES[6];
#pragma unroll
        for (int fo = 0; fo < 6; ++fo) {
            int v = polyFaces[p * 6 + fo];
            FACES[fo] = (v < 0) ? (-v - 1) : v;
        }

        int raw_f1 = polyFaces[p * 6 + 0];
        bool f1Reverse = (raw_f1 < 0);
        int f1 = (raw_f1 < 0) ? (-raw_f1 - 1) : raw_f1;

        int F1E[4];
#pragma unroll
        for (int eOff = 0; eOff < 4; ++eOff) {
            int v = faceEdges[f1 * 4 + eOff];
            F1E[eOff] = (v < 0) ? (-v - 1) : v;
        }

        int f2 = -1;
        for (int fo = 1; fo < 6; ++fo) {
            int f = FACES[fo];
            bool share = false;
#pragma unroll
            for (int eOff = 0; eOff < 4; ++eOff) {
                int e = faceEdges[f * 4 + eOff];
                e = (e < 0) ? (-e - 1) : e;
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
                uint e = (uint)((faceEdges[f * 4 + eOff] < 0) ? (-faceEdges[f * 4 + eOff] - 1) : faceEdges[f * 4 + eOff]);
                if (!contains_u_u32(EDGES, nEDGES, e) && nEDGES < 12) EDGES[nEDGES++] = e;
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
            for (int t = 0; t < 12; ++t) {
                uint e = EDGES[t];
                bool used = false; for (int q = 0; q < nEV; ++q) { if (EV[q] == e) { used = true; break; } }
                if (!used) {
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
                bool has = false; for (int j = 0; j < nFV; ++j) { if (FV[j] == f) { has = true; break; } }
                if (!has) {
                    bool hit = false;
                    for (int eOff = 0; eOff < 4; ++eOff) {
                        int fe = faceEdges[f * 4 + eOff];
                        fe = (fe < 0) ? (-fe - 1) : fe;
                        if (fe == eNeed) { hit = true; break; }
                    }
                    if (hit) {
                        FV[nFV++] = f;
                        for (int eOff = 0; eOff < 4; ++eOff) {
                            int fe = faceEdges[f * 4 + eOff];
                            fe = (fe < 0) ? (-fe - 1) : fe;
                            bool used = false; for (int q = 0; q < nEV; ++q) { if (EV[q] == (uint)fe) { used = true; break; } }
                            if (!used) { EV[nEV++] = (uint)fe; break; }
                        }
                    }
                }
            }
        }

        uint PV = (uint)p + pvOffset;
        for (int i = 0; i < 6; ++i) FV[i] += fvOffset;
        for (int i = 0; i < 12; ++i) EV[i] += evOffset;
        for (int i = 0; i < 8; ++i) VV[i] += vvOffset;

        uint base = (uint)p * 64u;
        out[base + 0] = VV[0]; out[base + 1] = EV[3];  out[base + 2] = FV[0]; out[base + 3] = EV[0];
        out[base + 4] = EV[4]; out[base + 5] = FV[5];  out[base + 6] = PV;    out[base + 7] = FV[2];
        out[base + 8] = EV[0]; out[base + 9] = FV[0];  out[base + 10] = EV[1]; out[base + 11] = VV[1];
        out[base + 12] = FV[2]; out[base + 13] = PV;     out[base + 14] = FV[3]; out[base + 15] = EV[5];
        out[base + 16] = FV[0]; out[base + 17] = EV[2];  out[base + 18] = VV[2]; out[base + 19] = EV[1];
        out[base + 20] = PV;    out[base + 21] = FV[4];  out[base + 22] = EV[6]; out[base + 23] = FV[3];
        out[base + 24] = EV[3]; out[base + 25] = VV[3];  out[base + 26] = EV[2]; out[base + 27] = FV[0];
        out[base + 28] = FV[5]; out[base + 29] = EV[7];  out[base + 30] = FV[4]; out[base + 31] = PV;
        out[base + 32] = EV[4]; out[base + 33] = FV[5];  out[base + 34] = PV;    out[base + 35] = FV[2];
        out[base + 36] = VV[4]; out[base + 37] = EV[11]; out[base + 38] = FV[1]; out[base + 39] = EV[8];
        out[base + 40] = FV[2]; out[base + 41] = PV;     out[base + 42] = FV[3]; out[base + 43] = EV[5];
        out[base + 44] = EV[8]; out[base + 45] = FV[1];  out[base + 46] = EV[9]; out[base + 47] = VV[5];
        out[base + 48] = PV;    out[base + 49] = FV[4];  out[base + 50] = EV[6]; out[base + 51] = FV[3];
        out[base + 52] = FV[1]; out[base + 53] = EV[10]; out[base + 54] = VV[6]; out[base + 55] = EV[9];
        out[base + 56] = FV[5]; out[base + 57] = EV[7];  out[base + 58] = FV[4]; out[base + 59] = PV;
        out[base + 60] = EV[11]; out[base + 61] = VV[7];  out[base + 62] = EV[10]; out[base + 63] = FV[1];
    }

    // ==================== 主流程：显式 CUDA Graph（无捕获） ====================

    void Patch::singlePatch::subdiv_cuda(std::vector<vec3d>& pos, std::vector<uint>& polys)
    {
        constexpr int BLOCK = 128;
        constexpr int MAX_DEG = 32;

        const uint ne = static_cast<uint>(edgeVerts.size());
        const uint nf = static_cast<uint>(faceEdges.size() / 4);
        const uint np = static_cast<uint>(polyFaces.size() / 6);
        const uint nv = static_cast<uint>(vertsPos.size());

        auto ceil_div = [](int a, int b) { return (a + b - 1) / b; };
        dim3 gE(ceil_div((int)ne, BLOCK)), gF(ceil_div((int)nf, BLOCK)), gP(ceil_div((int)np, BLOCK));
        dim3 gPF(ceil_div((int)patchFaces, BLOCK)), gPE(ceil_div((int)patchEdges, BLOCK)),
            gPV(ceil_div((int)patchVerts, BLOCK)), gPP(ceil_div((int)patchPolys, BLOCK));

        // ---- Host buffers ----
        std::vector<dvec3> hV(nv);
        for (uint i = 0; i < nv; ++i) hV[i] = to_d(vertsPos[i]);

        std::vector<uint2> hEV(ne);
        for (uint e = 0; e < ne; ++e) hEV[e] = make_uint2(edgeVerts[e].x(), edgeVerts[e].y());

        auto pack_bool = [](const std::vector<bool>& v) {
            std::vector<uint8_t> out(v.size());
            for (size_t i = 0; i < v.size(); ++i) out[i] = v[i] ? 1u : 0u;
            return out;
        };
        std::vector<uint8_t> h_vOn = pack_bool(vertOnSurf);
        std::vector<uint8_t> h_eOn = pack_bool(edgeOnSurf);
        std::vector<uint8_t> h_fOn = pack_bool(faceOnSurf);

        // ---- Device allocations ----
        dvec3* dV = nullptr, * d_EC = nullptr, * d_FC = nullptr, * d_PC = nullptr;
        dvec3* d_newPoly = nullptr, * d_newFace = nullptr, * d_newEdge = nullptr, * d_newVert = nullptr;
        uint2* dEV = nullptr, * dEV2 = nullptr;
        int* d_faceEdges = nullptr, * d_polyFaces = nullptr, * d_vertEdges = nullptr, * d_edgeFaces = nullptr, * d_facePolys = nullptr;
        uint* d_veOff = nullptr, * d_efOff = nullptr, * d_fpOff = nullptr;
        uint* d_topo = nullptr;
        uint8_t* d_vOn = nullptr, * d_eOn = nullptr, * d_fOn = nullptr;

        CUDA_CHECK(cudaMalloc(&dV, nv * sizeof(dvec3)));
        CUDA_CHECK(cudaMalloc(&dEV, ne * sizeof(uint2)));
        CUDA_CHECK(cudaMalloc(&dEV2, ne * sizeof(uint2)));
        CUDA_CHECK(cudaMalloc(&d_EC, ne * sizeof(dvec3)));
        CUDA_CHECK(cudaMalloc(&d_FC, nf * sizeof(dvec3)));
        CUDA_CHECK(cudaMalloc(&d_PC, np * sizeof(dvec3)));
        CUDA_CHECK(cudaMalloc(&d_newPoly, patchPolys * sizeof(dvec3)));
        CUDA_CHECK(cudaMalloc(&d_newFace, patchFaces * sizeof(dvec3)));
        CUDA_CHECK(cudaMalloc(&d_newEdge, patchEdges * sizeof(dvec3)));
        CUDA_CHECK(cudaMalloc(&d_newVert, patchVerts * sizeof(dvec3)));
        CUDA_CHECK(cudaMalloc(&d_topo, patchPolys * 64u * sizeof(uint)));

        CUDA_CHECK(cudaMalloc(&d_faceEdges, (size_t)faceEdges.size() * sizeof(int)));
        CUDA_CHECK(cudaMalloc(&d_polyFaces, (size_t)polyFaces.size() * sizeof(int)));
        CUDA_CHECK(cudaMalloc(&d_vertEdges, (size_t)vertEdges.size() * sizeof(int)));
        CUDA_CHECK(cudaMalloc(&d_edgeFaces, (size_t)edgeFaces.size() * sizeof(int)));
        CUDA_CHECK(cudaMalloc(&d_facePolys, (size_t)facePolys.size() * sizeof(int)));

        CUDA_CHECK(cudaMalloc(&d_veOff, (size_t)vertEdgesOffset.size() * sizeof(uint)));
        CUDA_CHECK(cudaMalloc(&d_efOff, (size_t)edgeFacesOffset.size() * sizeof(uint)));
        CUDA_CHECK(cudaMalloc(&d_fpOff, (size_t)facePolysOffset.size() * sizeof(uint)));

        CUDA_CHECK(cudaMalloc(&d_vOn, (size_t)h_vOn.size() * sizeof(uint8_t)));
        CUDA_CHECK(cudaMalloc(&d_eOn, (size_t)h_eOn.size() * sizeof(uint8_t)));
        CUDA_CHECK(cudaMalloc(&d_fOn, (size_t)h_fOn.size() * sizeof(uint8_t)));

        // ---- H2D ----
        CUDA_CHECK(cudaMemcpy(dV, hV.data(), nv * sizeof(dvec3), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(dEV, hEV.data(), ne * sizeof(uint2), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(dEV2, hEV.data(), ne * sizeof(uint2), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_faceEdges, faceEdges.data(), faceEdges.size() * sizeof(int), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_polyFaces, polyFaces.data(), polyFaces.size() * sizeof(int), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_vertEdges, vertEdges.data(), vertEdges.size() * sizeof(int), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_edgeFaces, edgeFaces.data(), edgeFaces.size() * sizeof(int), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_facePolys, facePolys.data(), facePolys.size() * sizeof(int), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_veOff, vertEdgesOffset.data(), vertEdgesOffset.size() * sizeof(uint), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_efOff, edgeFacesOffset.data(), edgeFacesOffset.size() * sizeof(uint), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_fpOff, facePolysOffset.data(), facePolysOffset.size() * sizeof(uint), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_vOn, h_vOn.data(), h_vOn.size() * sizeof(uint8_t), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_eOn, h_eOn.data(), h_eOn.size() * sizeof(uint8_t), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_fOn, h_fOn.data(), h_fOn.size() * sizeof(uint8_t), cudaMemcpyHostToDevice));

        // ---- Topo offsets ----
        uint pvOffset = static_cast<uint>(pos.size());
        uint fvOffset = pvOffset + patchPolys;
        uint evOffset = fvOffset + patchFaces;
        uint vvOffset = evOffset + patchEdges;

        // ===================== 显式 Graph 构建（无 stream capture） =====================
        cudaGraph_t graph; CUDA_CHECK(cudaGraphCreate(&graph, 0));

        // 1) EC
        int ne_i = (int)ne;
        cudaKernelNodeParams nEC{}; void* ecArgs[] = { &dEV, &dV, &d_EC, &ne_i };
        nEC.func = (void*)EdgeCentroidKernel;
        nEC.gridDim = gE; nEC.blockDim = dim3(BLOCK); nEC.sharedMemBytes = 0;
        nEC.kernelParams = ecArgs; nEC.extra = nullptr;
        cudaGraphNode_t EC; CUDA_CHECK(cudaGraphAddKernelNode(&EC, graph, nullptr, 0, &nEC));

        // 2) FC (dep: EC)
        int nf_i = (int)nf;
        cudaKernelNodeParams nFC{}; void* fcArgs[] = { &d_faceEdges, &d_EC, &d_FC, &nf_i };
        nFC.func = (void*)FaceCentroidKernel;
        nFC.gridDim = gF; nFC.blockDim = dim3(BLOCK); nFC.sharedMemBytes = 0;
        nFC.kernelParams = fcArgs; nFC.extra = nullptr;
        cudaGraphNode_t FC; CUDA_CHECK(cudaGraphAddKernelNode(&FC, graph, &EC, 1, &nFC));

        // 3) PC (dep: FC)
        int np_i = (int)np;
        cudaKernelNodeParams nPC{}; void* pcArgs[] = { &d_polyFaces, &d_FC, &d_PC, &np_i };
        nPC.func = (void*)PolyCentroidKernel;
        nPC.gridDim = gP; nPC.blockDim = dim3(BLOCK); nPC.sharedMemBytes = 0;
        nPC.kernelParams = pcArgs; nPC.extra = nullptr;
        cudaGraphNode_t PC; CUDA_CHECK(cudaGraphAddKernelNode(&PC, graph, &FC, 1, &nPC));

        // 4) new poly = first patchPolys of PC (dep: PC)
        int patchPolys_i = (int)patchPolys;
        cudaKernelNodeParams nCopyP{}; void* cpArgs[] = { &d_PC, &d_newPoly, &patchPolys_i };
        nCopyP.func = (void*)CopyFirstNKernel;
        nCopyP.gridDim = gPP; nCopyP.blockDim = dim3(BLOCK); nCopyP.sharedMemBytes = 0;
        nCopyP.kernelParams = cpArgs; nCopyP.extra = nullptr;
        cudaGraphNode_t CP; CUDA_CHECK(cudaGraphAddKernelNode(&CP, graph, &PC, 1, &nCopyP));

        // 5) NewFace (dep: PC)
        int patchFaces_i = (int)patchFaces;
        cudaKernelNodeParams nNF{}; void* nfArgs[] = { &d_fOn,&d_fpOff,&d_facePolys,&d_FC,&d_PC,&d_newFace,&patchFaces_i };
        nNF.func = (void*)NewFaceVertKernel;
        nNF.gridDim = gPF; nNF.blockDim = dim3(BLOCK); nNF.sharedMemBytes = 0;
        nNF.kernelParams = nfArgs; nNF.extra = nullptr;
        cudaGraphNode_t NF; CUDA_CHECK(cudaGraphAddKernelNode(&NF, graph, &PC, 1, &nNF));

        // 6) NewEdge (dep: PC)
        int patchEdges_i = (int)patchEdges;
        size_t shmemEdge = BLOCK * (MAX_DEG + MAX_DEG) * sizeof(int);
        cudaKernelNodeParams nNE{}; void* neArgs[] = { &d_eOn,&d_fOn,&d_efOff,&d_edgeFaces,&d_fpOff,&d_facePolys,
                                                       &d_EC,&d_FC,&d_PC,&d_newEdge,&patchEdges_i };
        nNE.func = (void*)NewEdgeVertKernel<BLOCK, MAX_DEG>;
        nNE.gridDim = gPE; nNE.blockDim = dim3(BLOCK); nNE.sharedMemBytes = (unsigned)shmemEdge;
        nNE.kernelParams = neArgs; nNE.extra = nullptr;
        cudaGraphNode_t NE; CUDA_CHECK(cudaGraphAddKernelNode(&NE, graph, &PC, 1, &nNE));

        // 7) NewVert (dep: PC)
        int patchVerts_i = (int)patchVerts;
        size_t shmemVert = BLOCK * (3 * MAX_DEG) * sizeof(int);
        cudaKernelNodeParams nNV{}; void* nvArgs[] = { &d_vOn,&d_eOn,&d_fOn,&d_veOff,&d_vertEdges,&d_efOff,&d_edgeFaces,
                                                       &d_fpOff,&d_facePolys,&dV,&d_EC,&d_FC,&d_PC,&d_newVert,&patchVerts_i };
        nNV.func = (void*)NewVertVertKernel<BLOCK, MAX_DEG>;
        nNV.gridDim = gPV; nNV.blockDim = dim3(BLOCK); nNV.sharedMemBytes = (unsigned)shmemVert;
        nNV.kernelParams = nvArgs; nNV.extra = nullptr;
        cudaGraphNode_t NV; CUDA_CHECK(cudaGraphAddKernelNode(&NV, graph, &PC, 1, &nNV));

        // 8) Topo（独立节点；若想更稳妥可依赖 PC）
        int patchPolys_i2 = (int)patchPolys;
        cudaKernelNodeParams nTopo{}; void* topoArgs[] = { &d_polyFaces,&d_faceEdges,&dEV2,
                                                           &pvOffset,&fvOffset,&evOffset,&vvOffset,
                                                           &d_topo,&patchPolys_i2 };
        nTopo.func = (void*)TopoAssembleKernel;
        nTopo.gridDim = gPP; nTopo.blockDim = dim3(BLOCK); nTopo.sharedMemBytes = 0;
        nTopo.kernelParams = topoArgs; nTopo.extra = nullptr;
        cudaGraphNode_t TOPO; CUDA_CHECK(cudaGraphAddKernelNode(&TOPO, graph, nullptr, 0, &nTopo));
        // 如果希望 Topo 也等 PC：CUDA_CHECK(cudaGraphAddDependencies(graph, &PC, &TOPO, 1));

        // 实例化 & 执行
        cudaGraphExec_t exec;
        CUDA_CHECK(cudaGraphInstantiate(&exec, graph, nullptr, nullptr, 0));
        auto t0 = std::chrono::high_resolution_clock::now();
        CUDA_CHECK(cudaGraphLaunch(exec, 0)); // 默认流即可
        CUDA_CHECK(cudaStreamSynchronize(0));
        auto t1 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> elapsed = t1 - t0;
        cuda_elapsed += elapsed.count();
        std::cout << "execution time of current patch: " << elapsed.count() << " ms\n";
        outlog << "execution time of current patch: " << elapsed.count() << " ms\n";

        // ---- D2H ----
        std::vector<dvec3> h_newPoly(patchPolys), h_newFace(patchFaces), h_newEdge(patchEdges), h_newVert(patchVerts);
        std::vector<uint>  h_topo(patchPolys * 64u);

        CUDA_CHECK(cudaMemcpy(h_newPoly.data(), d_newPoly, patchPolys * sizeof(dvec3), cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(h_newFace.data(), d_newFace, patchFaces * sizeof(dvec3), cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(h_newEdge.data(), d_newEdge, patchEdges * sizeof(dvec3), cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(h_newVert.data(), d_newVert, patchVerts * sizeof(dvec3), cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(h_topo.data(), d_topo, h_topo.size() * sizeof(uint), cudaMemcpyDeviceToHost));

        // ---- 追加到输出 ----
        pos.reserve(pos.size() + h_newPoly.size() + h_newFace.size() + h_newEdge.size() + h_newVert.size());
        for (const auto& v : h_newPoly) pos.push_back(to_h(v));
        for (const auto& v : h_newFace) pos.push_back(to_h(v));
        for (const auto& v : h_newEdge) pos.push_back(to_h(v));
        for (const auto& v : h_newVert) pos.push_back(to_h(v));

        polys.reserve(polys.size() + h_topo.size());
        polys.insert(polys.end(), h_topo.begin(), h_topo.end());

        // ---- 清理 ----
        CUDA_CHECK(cudaGraphExecDestroy(exec));
        CUDA_CHECK(cudaGraphDestroy(graph));

        CUDA_CHECK(cudaFree(dV));   CUDA_CHECK(cudaFree(dEV));   CUDA_CHECK(cudaFree(dEV2));
        CUDA_CHECK(cudaFree(d_EC)); CUDA_CHECK(cudaFree(d_FC));  CUDA_CHECK(cudaFree(d_PC));
        CUDA_CHECK(cudaFree(d_newPoly)); CUDA_CHECK(cudaFree(d_newFace));
        CUDA_CHECK(cudaFree(d_newEdge)); CUDA_CHECK(cudaFree(d_newVert)); CUDA_CHECK(cudaFree(d_topo));
        CUDA_CHECK(cudaFree(d_faceEdges)); CUDA_CHECK(cudaFree(d_polyFaces));
        CUDA_CHECK(cudaFree(d_vertEdges)); CUDA_CHECK(cudaFree(d_edgeFaces)); CUDA_CHECK(cudaFree(d_facePolys));
        CUDA_CHECK(cudaFree(d_veOff)); CUDA_CHECK(cudaFree(d_efOff)); CUDA_CHECK(cudaFree(d_fpOff));
        CUDA_CHECK(cudaFree(d_vOn)); CUDA_CHECK(cudaFree(d_eOn)); CUDA_CHECK(cudaFree(d_fOn));
    }

} // namespace cinolib
