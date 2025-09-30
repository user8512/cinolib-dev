#pragma once
#include <string>
#include <cinolib/meshes/meshes.h>
namespace cinolib
{
	class Patch
	{
	public:
		Patch();
		Patch(std::string meshFileName);
		Patch(const std::vector<vec3d>& verts, const std::vector<uint>& polys);
		Patch(const std::vector<vec3d>& verts, const std::vector<std::vector<uint>>& polys);
		~Patch() = default;

		class singlePatch {
		public:
			singlePatch() {};
			singlePatch(std::vector<uint> polys, uint pPolys, uint pFaces, uint pEdges, uint pVerts, std::vector<vec3d> v_pos, 
					std::vector<vec2u> ev, std::vector<int> fe, std::vector<int> cf,std::vector<int> ve, std::vector<int> ef, std::vector<int> fc, 
					std::vector<uint> veoff, std::vector<uint> efoff, std::vector<uint> fpoff, std::vector<bool> vos, std::vector<bool> eos, std::vector<bool> fos):
				polyGlobalIndex(std::move(polys)), patchPolys(pPolys), patchFaces(pFaces), patchEdges(pEdges), patchVerts(pVerts), vertsPos(std::move(v_pos)),
					edgeVerts(std::move(ev)), faceEdges(std::move(fe)), polyFaces(std::move(cf)), vertEdges(std::move(ve)), 
					edgeFaces(std::move(ef)), facePolys(std::move(fc)), vertEdgesOffset(std::move(veoff)), edgeFacesOffset(std::move(efoff)), facePolysOffset(std::move(fpoff)),
					vertOnSurf(std::move(vos)), edgeOnSurf(std::move(eos)), faceOnSurf(std::move(fos)){}
			~singlePatch() = default;
			void subdiv(std::vector<vec3d>& pos, std::vector<uint>& polys);

		private:
			uint patchPolys;
			uint patchFaces;
			uint patchEdges;
			uint patchVerts;
			std::vector<uint> polyGlobalIndex;
			std::vector<uint> vertGlobalIndex;
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
		
		void patching(std::string patchLabelFileName, int num_clusters);
		void Patch::subdiv();
		Hexmesh<> getMesh() { return mesh; }

	private:
		Hexmesh<> mesh;
		std::vector<singlePatch> patches;

		bool getFaceEdgeSign(uint eid, uint fid);
		bool getPolyFaceSign(uint fid, uint pid);
	};

	void PrintVec3d(vec3d& v);
}