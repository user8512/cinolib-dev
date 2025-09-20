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
		Patch(const std::vector<vec3d>& verts, const std::vector<std::vector<uint>>& polys);
		~Patch() = default;

		class singlePatch {
		public:
			singlePatch() {};
			singlePatch(std::vector<uint> polys, uint pSize, std::vector<vec3d> v_pos, std::vector<vec2d> ev, std::vector<int> fe,
					std::vector<int> cf, std::vector<uint> ve, std::vector<int> ef, std::vector<int> fc):
				polyGlobalIndex(std::move(polys)), patchSize(pSize), vertsPos(std::move(v_pos)), edgeVerts(std::move(ev)), faceEdges(std::move(fe)),
					cellFaces(std::move(cf)), vertEdges(std::move(ve)), edgeFaces(std::move(ef)), faceCells(std::move(fc)){}
			~singlePatch() = default;

		private:
			uint patchSize;
			std::vector<uint> polyGlobalIndex;
			std::vector<uint> vertGlobalIndex;
			std::vector<vec3d> vertsPos;

			// patches and ribbons with localVerts index
			std::vector<vec2d> edgeVerts;
			std::vector<int> faceEdges;
			std::vector<int> cellFaces;
			std::vector<uint> vertEdges;
			std::vector<int> edgeFaces;
			std::vector<int> faceCells;
		};
		
		void patching(std::string patchLabelFileName, int num_clusters);
		Hexmesh<> getMesh() { return mesh; }

	private:
		int maxPatchSize = 0;
		Hexmesh<> mesh;
		std::vector<singlePatch> patches;

		bool getFaceEdgeSign(uint eid, uint fid);
		bool getPolyFaceSign(uint fid, uint pid);

		inline void getBoundaryPolys(std::vector<uint>& polys) {
			std::vector<uint> newpolys;
			for (auto pid : polys) {
				if (mesh.poly_is_on_surf(pid)) {
					newpolys.push_back(pid);
				}
			}
			polys = newpolys;
		}

		void setMaxPatchSize(int size) {
			assert(size > 0);
			maxPatchSize = size;
		}
	};
}