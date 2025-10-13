#pragma once
#include <cinolib/meshes/meshes.h>
#include <chrono>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <array>
#include <stack>
#include <cstdint>
#include <random>
#include <iterator>
#include <tuple>
#include <cmath>
#include <limits>
#include <functional>
#include <cassert>
#include <cstddef>

#define USE_CUDA
//#define TEST
//#define DRAW
//#define DEBUG
#define OUTPUT
//#define DETAIL
//#define OUTPUT_DETAIL

#define MAX_POLYS_PER_CLUSTER 32768
#define MAX_CLUSTER 1
#define root std::string(DATA_PATH)

namespace cinolib {
	extern float cuda_elapsed;
	extern float cpu_elapsed;

	class Patch {
	public:
		Patch();
		Patch(std::string meshFileName);
		Patch(const std::vector<vec3d>& verts, const std::vector<uint>& polys);
		Patch(const std::vector<vec3d>& verts, const std::vector<std::vector<uint>>& polys);
		~Patch() = default;

		class singlePatch;
		
		void patching(int num_clusters);
		void subdiv(int subdiv_times);
		// 用于各patch完成操作后，合并容差范围内的点
		void deduplicate_verts(std::vector<vec3d>& verts, std::vector<uint>& polys, double tol = 1e-6);

		Hexmesh<> *getMesh() { return &mesh; }
		bool getFaceEdgeSign(uint eid, uint fid);
		bool getPolyFaceSign(uint fid, uint pid);

	private:
		Hexmesh<> mesh;
		std::vector<singlePatch> patches;
	};

	class Patch::singlePatch {
	public:
		singlePatch(uint pPolys, uint pFaces, uint pEdges, uint pVerts, std::vector<vec3d> v_pos,
			std::vector<vec2u> ev, std::vector<int> fe, std::vector<int> cf, std::vector<int> ve, std::vector<int> ef, std::vector<int> fc,
			std::vector<uint> veoff, std::vector<uint> efoff, std::vector<uint> fpoff, std::vector<bool> vos, std::vector<bool> eos, std::vector<bool> fos) :
			patchPolys(pPolys), patchFaces(pFaces), patchEdges(pEdges), patchVerts(pVerts), vertsPos(std::move(v_pos)),
			edgeVerts(std::move(ev)), faceEdges(std::move(fe)), polyFaces(std::move(cf)), vertEdges(std::move(ve)),
			edgeFaces(std::move(ef)), facePolys(std::move(fc)), vertEdgesOffset(std::move(veoff)), edgeFacesOffset(std::move(efoff)), facePolysOffset(std::move(fpoff)),
			vertOnSurf(std::move(vos)), edgeOnSurf(std::move(eos)), faceOnSurf(std::move(fos)) {}
		~singlePatch() = default;
		void subdiv(std::vector<vec3d>& pos, std::vector<uint>& polys);
		void subdiv_cuda(std::vector<vec3d>& pos, std::vector<uint>& polys);

	private:
		uint patchPolys;
		uint patchFaces;
		uint patchEdges;
		uint patchVerts;
		std::vector<vec3d> vertsPos;

		// patches and ribbons with localVerts index
		std::vector<vec2u> edgeVerts;
		std::vector<int> faceEdges;
		std::vector<int> polyFaces;
		std::vector<int> vertEdges;
		std::vector<int> edgeFaces;
		std::vector<int> facePolys;
		std::vector<uint> vertEdgesOffset;
		std::vector<uint> edgeFacesOffset;
		std::vector<uint> facePolysOffset;
		std::vector<bool> vertOnSurf;
		std::vector<bool> edgeOnSurf;
		std::vector<bool> faceOnSurf;
	};

	inline void PrintVec3d(vec3d& v);
	inline bool Vec3dEqual(vec3d& v1, vec3d& v2, double tolerance = 1e-6);
}