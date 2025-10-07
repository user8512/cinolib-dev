#include "cc_subdiv_gpu.h"
#include <thrust/device_vector.h>
#include <thrust/for_each.h>
#include <thrust/copy.h>
#include <thrust/iterator/counting_iterator.h>
#include <cinolib/geometry/vec_mat.h>

namespace {
// 可在 device 使用的轻量向量
struct Vec3 {
    double x, y, z;
    __host__ __device__ Vec3(double X=0, double Y=0, double Z=0):x(X),y(Y),z(Z){}
    __host__ __device__ Vec3& operator+=(const Vec3& o){ x+=o.x; y+=o.y; z+=o.z; return *this; }
    __host__ __device__ Vec3  operator+(const Vec3& o) const { return Vec3(x+o.x,y+o.y,z+o.z); }
    __host__ __device__ Vec3  operator*(double s) const { return Vec3(x*s,y*s,z*s); }
    __host__ __device__ Vec3  operator/(double s) const { return Vec3(x/s,y/s,z/s); }
};

struct I2 {
    int a, b;
    __host__ __device__ I2(int A=0, int B=0):a(A),b(B){}
    __host__ __device__ int x() const { return a; }
    __host__ __device__ int y() const { return b; }
};

__host__ __device__ inline int  pos_index(int v){ return (v<0)?(-v-1):v; }
__host__ __device__ inline Vec3 div_safe (const Vec3& v, int d){ return (d>0)? (v / (double)d) : v; }

constexpr int MAX_EDGE_POLYS  = 9999;
constexpr int MAX_VERT_FACES  = 9999;
constexpr int MAX_VERT_POLYS  = 9999;

}

namespace cinolib {

void subdiv_cuda(
    const std::vector<vec3d>&  h_vertsPos,
    int patchVerts, int patchEdges, int patchFaces, int patchPolys,
    const std::vector<vec2i>&  h_edgeVerts,
    const std::vector<int>&    h_faceEdges,
    const std::vector<int>&    h_polyFaces,
    const std::vector<uint8_t>& h_faceOnSurf,
    const std::vector<uint8_t>& h_edgeOnSurf,
    const std::vector<uint8_t>& h_vertOnSurf,
    const std::vector<int>&    h_facePolysOffset,
    const std::vector<int>&    h_facePolys,
    const std::vector<int>&    h_edgeFacesOffset,
    const std::vector<int>&    h_edgeFaces,
    const std::vector<int>&    h_vertEdgesOffset,
    const std::vector<int>&    h_vertEdges,
    std::vector<vec3d>&        newPolyVerts,
    std::vector<vec3d>&        newFaceVerts,
    std::vector<vec3d>&        newEdgeVerts,
    std::vector<vec3d>&        newVertVerts
){
    // host: vec3d/vec2i -> POD
    const int nv = (int)h_vertsPos.size();
    const int ne = (int)h_edgeVerts.size();
    const int nf = (int)h_faceEdges.size()/4;
    const int np = (int)h_polyFaces.size()/6;

    std::vector<Vec3> h_pos_pod(nv);
    for (int i=0;i<nv;++i) h_pos_pod[i] = Vec3(h_vertsPos[i].x(), h_vertsPos[i].y(), h_vertsPos[i].z());

    std::vector<I2> h_edge_pod(ne);
    for (int i=0;i<ne;++i) h_edge_pod[i] = I2(h_edgeVerts[i].x(), h_edgeVerts[i].y());

    // H->D
    thrust::device_vector<Vec3> d_pos(h_pos_pod.begin(), h_pos_pod.end());
    thrust::device_vector<I2>   d_edgeVerts(h_edge_pod.begin(), h_edge_pod.end());
    thrust::device_vector<int>  d_faceEdges(h_faceEdges.begin(), h_faceEdges.end());
    thrust::device_vector<int>  d_polyFaces(h_polyFaces.begin(), h_polyFaces.end());
    thrust::device_vector<uint8_t> d_faceOnSurf(h_faceOnSurf.begin(), h_faceOnSurf.end());
    thrust::device_vector<uint8_t> d_edgeOnSurf(h_edgeOnSurf.begin(), h_edgeOnSurf.end());
    thrust::device_vector<uint8_t> d_vertOnSurf(h_vertOnSurf.begin(), h_vertOnSurf.end());
    thrust::device_vector<int>  d_facePolysOffset(h_facePolysOffset.begin(), h_facePolysOffset.end());
    thrust::device_vector<int>  d_facePolys(h_facePolys.begin(), h_facePolys.end());
    thrust::device_vector<int>  d_edgeFacesOffset(h_edgeFacesOffset.begin(), h_edgeFacesOffset.end());
    thrust::device_vector<int>  d_edgeFaces(h_edgeFaces.begin(), h_edgeFaces.end());
    thrust::device_vector<int>  d_vertEdgesOffset(h_vertEdgesOffset.begin(), h_vertEdgesOffset.end());
    thrust::device_vector<int>  d_vertEdges(h_vertEdges.begin(), h_vertEdges.end());

    thrust::device_vector<Vec3> d_edgeC(ne), d_faceC(nf), d_polyC(np);
    thrust::device_vector<Vec3> d_newFace(patchFaces), d_newEdge(patchEdges), d_newVert(patchVerts), d_newPoly(patchPolys);

    auto c0 = thrust::make_counting_iterator<int>(0);

    // 1) 边中点
    {
        Vec3* edgeC = thrust::raw_pointer_cast(d_edgeC.data());
        const I2* eV= thrust::raw_pointer_cast(d_edgeVerts.data());
        const Vec3* pos = thrust::raw_pointer_cast(d_pos.data());
        thrust::for_each(c0, c0+ne, [=] __host__ __device__ (int e){
            edgeC[e] = (pos[eV[e].x()] + pos[eV[e].y()]) / 2.0;
        });
    }

    // 2) 面中点（4 边）
    {
        Vec3* faceC = thrust::raw_pointer_cast(d_faceC.data());
        const Vec3* edgeC = thrust::raw_pointer_cast(d_edgeC.data());
        const int* fE = thrust::raw_pointer_cast(d_faceEdges.data());
        thrust::for_each(c0, c0+nf, [=] __host__ __device__ (int f){
            Vec3 c(0,0,0);
            #pragma unroll
            for(int k=0;k<4;++k){
                int e = pos_index(fE[f*4+k]);
                c += edgeC[e];
            }
            faceC[f] = c/4.0;
        });
    }

    // 3) 体中点（6 面）
    {
        Vec3* polyC = thrust::raw_pointer_cast(d_polyC.data());
        const Vec3* faceC = thrust::raw_pointer_cast(d_faceC.data());
        const int* pF = thrust::raw_pointer_cast(d_polyFaces.data());
        thrust::for_each(c0, c0+np, [=] __host__ __device__ (int p){
            Vec3 c(0,0,0);
            #pragma unroll
            for(int k=0;k<6;++k){
                int f = pos_index(pF[p*6+k]);
                c += faceC[f];
            }
            polyC[p] = c/6.0;
        });
    }

    // 4) 新体点
    {
        const Vec3* polyC = thrust::raw_pointer_cast(d_polyC.data());
        Vec3* newP = thrust::raw_pointer_cast(d_newPoly.data());
        thrust::for_each(c0, c0+patchPolys, [=] __host__ __device__ (int i){
            newP[i] = polyC[i];
        });
    }

    // 5) 新面点
    {
        const Vec3* polyC = thrust::raw_pointer_cast(d_polyC.data());
        const Vec3* faceC = thrust::raw_pointer_cast(d_faceC.data());
        const uint8_t* faceSurf = thrust::raw_pointer_cast(d_faceOnSurf.data());
        const int* fPOff = thrust::raw_pointer_cast(d_facePolysOffset.data());
        const int* fP = thrust::raw_pointer_cast(d_facePolys.data());
        Vec3* newF = thrust::raw_pointer_cast(d_newFace.data());

        thrust::for_each(c0, c0+patchFaces, [=] __host__ __device__ (int f){
            if(!faceSurf[f]){
                int off = fPOff[f];
                int p0 = pos_index(fP[off]);
                int p1 = pos_index(fP[off+1]);
                newF[f] = (polyC[p0] + polyC[p1] + faceC[f]*2.0) / 4.0;
            }else{
                newF[f] = faceC[f];
            }
        });
    }

    // 6) 新边点
    {
        const Vec3* polyC = thrust::raw_pointer_cast(d_polyC.data());
        const Vec3* faceC = thrust::raw_pointer_cast(d_faceC.data());
        const Vec3* edgeC = thrust::raw_pointer_cast(d_edgeC.data());
        const uint8_t* edgeSurf = thrust::raw_pointer_cast(d_edgeOnSurf.data());
        const uint8_t* faceSurf = thrust::raw_pointer_cast(d_faceOnSurf.data());
        const int* eFO = thrust::raw_pointer_cast(d_edgeFacesOffset.data());
        const int* eF  = thrust::raw_pointer_cast(d_edgeFaces.data());
        const int* fPO = thrust::raw_pointer_cast(d_facePolysOffset.data());
        const int* fP  = thrust::raw_pointer_cast(d_facePolys.data());
        Vec3* newE = thrust::raw_pointer_cast(d_newEdge.data());

        thrust::for_each(c0, c0+patchEdges, [=] __host__ __device__ (int e){
            int b = eFO[e], ee = eFO[e+1], N = ee-b;
            if(!edgeSurf[e]){
                Vec3 faceAvg(0,0,0);
                int  polyBuf[MAX_EDGE_POLYS]; int polyCnt=0;

                for(int i=b;i<ee;++i){
                    int f = pos_index(eF[i]);
                    faceAvg += faceC[f];

                    int fb=fPO[f], fe=fPO[f+1];
                    for(int j=fb;j<fe;++j){
                        int p = pos_index(fP[j]);
                        bool seen=false;
                        for(int t=0;t<polyCnt;++t){ if(polyBuf[t]==p){seen=true;break;} }
                        if(!seen && polyCnt<MAX_EDGE_POLYS) polyBuf[polyCnt++]=p;
                    }
                }
                faceAvg = div_safe(faceAvg,N);

                Vec3 polyAvg(0,0,0);
                for(int t=0;t<polyCnt;++t) polyAvg += polyC[ polyBuf[t] ];
                polyAvg = div_safe(polyAvg, polyCnt);

                newE[e] = (polyAvg + faceAvg*2.0 + edgeC[e]*(double)(N-3)) / (double)N;
            }else{
                if(N==1){
                    newE[e] = edgeC[e];
                }else{
                    Vec3 faceAvg(0,0,0); int cnt=0;
                    for(int i=b;i<ee;++i){
                        int f = pos_index(eF[i]);
                        if(faceSurf[f]){ faceAvg += faceC[f]; ++cnt; }
                    }
                    faceAvg = div_safe(faceAvg, cnt);
                    newE[e] = (faceAvg + edgeC[e]) / 2.0;
                }
            }
        });
    }

    // 7) 新点（原顶点）
    {
        const Vec3* polyC = thrust::raw_pointer_cast(d_polyC.data());
        const Vec3* faceC = thrust::raw_pointer_cast(d_faceC.data());
        const Vec3* edgeC = thrust::raw_pointer_cast(d_edgeC.data());
        const Vec3* pos   = thrust::raw_pointer_cast(d_pos.data());
        const uint8_t* vSurf = thrust::raw_pointer_cast(d_vertOnSurf.data());
        const uint8_t* eSurf = thrust::raw_pointer_cast(d_edgeOnSurf.data());
        const uint8_t* fSurf = thrust::raw_pointer_cast(d_faceOnSurf.data());
        const int* vEO = thrust::raw_pointer_cast(d_vertEdgesOffset.data());
        const int* vE  = thrust::raw_pointer_cast(d_vertEdges.data());
        const int* eFO = thrust::raw_pointer_cast(d_edgeFacesOffset.data());
        const int* eF  = thrust::raw_pointer_cast(d_edgeFaces.data());
        const int* fPO = thrust::raw_pointer_cast(d_facePolysOffset.data());
        const int* fP  = thrust::raw_pointer_cast(d_facePolys.data());
        Vec3* newV = thrust::raw_pointer_cast(d_newVert.data());

        thrust::for_each(c0, c0+patchVerts, [=] __host__ __device__ (int v){
            int b=vEO[v], e=vEO[v+1], N=e-b;

            if(!vSurf[v]){
                Vec3 edgeAvg(0,0,0), faceAvg(0,0,0), polyAvg(0,0,0);
                int faceBuf[MAX_VERT_FACES]; int faceCnt=0;
                int polyBuf[MAX_VERT_POLYS]; int polyCnt=0;

                for(int i=b;i<e;++i){
                    int ed = pos_index(vE[i]);
                    edgeAvg += edgeC[ed];

                    int eb=eFO[ed], ee=eFO[ed+1];
                    for(int j=eb;j<ee;++j){
                        int f = pos_index(eF[j]);
                        bool seen=false; for(int t=0;t<faceCnt;++t){ if(faceBuf[t]==f){seen=true;break;} }
                        if(!seen && faceCnt<MAX_VERT_FACES){
                            faceBuf[faceCnt++]=f;
                            int fb=fPO[f], fe=fPO[f+1];
                            for(int k=fb;k<fe;++k){
                                int p = pos_index(fP[k]);
                                bool s2=false; for(int u=0;u<polyCnt;++u){ if(polyBuf[u]==p){s2=true;break;} }
                                if(!s2 && polyCnt<MAX_VERT_POLYS) polyBuf[polyCnt++]=p;
                            }
                        }
                    }
                }
                edgeAvg = div_safe(edgeAvg, N);
                for(int t=0;t<faceCnt;++t) faceAvg += faceC[ faceBuf[t] ];
                faceAvg = div_safe(faceAvg, faceCnt);
                for(int t=0;t<polyCnt;++t) polyAvg += polyC[ polyBuf[t] ];
                polyAvg = div_safe(polyAvg, polyCnt);

                newV[v] = (polyAvg + faceAvg*3.0 + edgeAvg*3.0 + pos[v]) / 8.0;
            }else{
                Vec3 edgeAvg(0,0,0), faceAvg(0,0,0);
                int edgeCnt=0, faceCnt=0;

                for(int i=b;i<e;++i){
                    int ed = pos_index(vE[i]);
                    if(eSurf[ed]){
                        edgeAvg += edgeC[ed]; ++edgeCnt;
                        int eb=eFO[ed], ee=eFO[ed+1];
                        for(int j=eb;j<ee;++j){
                            int f = pos_index(eF[j]);
                            if(fSurf[f]){ faceAvg += faceC[f]; ++faceCnt; }
                        }
                    }
                }
                edgeAvg = div_safe(edgeAvg, edgeCnt);
                faceAvg = div_safe(faceAvg, faceCnt);
                int n = edgeCnt;
                newV[v] = (faceAvg + edgeAvg*2.0 + pos[v]*(double)(n-3)) / (double)(n>0?n:1);
            }
        });
    }

    // D->H (POD -> vec3d)；避免使用 vec3d(0,0,0) 的 assign，直接 resize + 填值
    std::vector<Vec3> h_newP(patchPolys), h_newF(patchFaces), h_newE(patchEdges), h_newV(patchVerts);
    thrust::copy(d_newPoly.begin(), d_newPoly.end(), h_newP.begin());
    thrust::copy(d_newFace.begin(), d_newFace.end(), h_newF.begin());
    thrust::copy(d_newEdge.begin(), d_newEdge.end(), h_newE.begin());
    thrust::copy(d_newVert.begin(), d_newVert.end(), h_newV.begin());

    newPolyVerts.resize(patchPolys);
    newFaceVerts.resize(patchFaces);
    newEdgeVerts.resize(patchEdges);
    newVertVerts.resize(patchVerts);
    for (int i=0;i<patchPolys;++i) newPolyVerts[i] = vec3d(h_newP[i].x, h_newP[i].y, h_newP[i].z);
    for (int i=0;i<patchFaces;++i) newFaceVerts[i] = vec3d(h_newF[i].x, h_newF[i].y, h_newF[i].z);
    for (int i=0;i<patchEdges;++i) newEdgeVerts[i] = vec3d(h_newE[i].x, h_newE[i].y, h_newE[i].z);
    for (int i=0;i<patchVerts;++i) newVertVerts[i] = vec3d(h_newV[i].x, h_newV[i].y, h_newV[i].z);
}

} // namespace gpu_cc
