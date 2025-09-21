#include "patch.h"
#include <stack>
#include <algorithm>
#include <random>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdint>
#include <iterator>

#define DEBUG
#define OUTPUT
#define DETAIL
#define OUTPUT_DETAIL

namespace cinolib {
	Patch::Patch() {
	}

	Patch::Patch(std::string fileName) {
		mesh.load(fileName.c_str());
	}

	Patch::Patch(const std::vector<vec3d>& verts, const std::vector<std::vector<uint>>& polys) {
		mesh = Hexmesh<>(verts, polys);
	}

	bool Patch::getFaceEdgeSign(uint fid, uint eid) {
		vec3d start = mesh.edge_vert(eid, 0);
		vec3d end = mesh.edge_vert(eid, 1);
		auto fv = mesh.face_verts(fid);
		if (fv[0] == start && fv[1] == end) {
			return true;
		}
		if (fv[1] == start && fv[2] == end) {
			return true;
		}
		if (fv[2] == start && fv[3] == end) {
			return true;
		}
		if (fv[3] == start && fv[0] == end) {
			return true;
		}
		return false;
	}

	bool Patch::getPolyFaceSign(uint pid, uint fid) {
		return mesh.poly_face_is_CCW(pid, fid);
	}

	void Patch::patching(std::string patchLabelFileName, int num_clusters) {
#ifdef OUTPUT
		std::ofstream outfile("D:/data/clustered_hexa/test/patch.txt");
		if (!outfile) {
			std::cerr << "无法打开文件" << std::endl;
		}
		std::ofstream outfile2("D:/data/clustered_hexa/test/ribbon.txt");
		if (!outfile2) {
			std::cerr << "无法打开文件" << std::endl;
		}
		std::ofstream outfile3("D:/data/clustered_hexa/test/detail.txt");
		if (!outfile3) {
			std::cerr << "无法打开文件" << std::endl;
		}
#endif

		std::ifstream file(patchLabelFileName, std::ios::binary);
		if (!file) {
			std::cerr << "无法打开文件" << std::endl;
		}

		file.seekg(0, std::ios::end);
		std::streamsize size = file.tellg();
		file.seekg(0, std::ios::beg);

		size_t num_elements = size / sizeof(int32_t);
		

		std::vector<int32_t> patchLabel(num_elements);
		file.read(reinterpret_cast<char*>(patchLabel.data()), size);

		// temporary, global idx
		std::vector<uint> patch;
		std::vector<uint> verts;
		std::vector<uint> edges;
		std::vector<uint> faces;

		// index of polys of patch and ribbon in order
		std::vector<uint> polys;
		// cells num of patch(no ribbon)
		uint patchSize;
		// pos of each vert
		std::vector<vec3d> v_pos;

		std::map<uint, uint> vertIdxGlobal2Local;
		std::map<uint, uint> edgeIdxGlobal2Local;
		std::map<uint, uint> faceIdxGlobal2Local;
		std::map<uint, uint> polyIdxGlobal2Local;

		std::vector<vec2u> edgeVerts;
		std::vector<int> faceEdges;
		std::vector<int> polyFaces;
		std::vector<int> vertEdges;
		std::vector<int> edgeFaces;
		std::vector<int> facePolys;
		std::vector<uint> vertEdgesOffset;
		std::vector<uint> edgeFacesOffset;
		std::vector<uint> facePolysOffset;

		for (int cluster = 0; cluster < num_clusters; cluster++) {
			patchSize = 0;
			patch.clear();
			verts.clear();
			edges.clear();
			faces.clear();
			vertIdxGlobal2Local.clear();
			edgeIdxGlobal2Local.clear();
			faceIdxGlobal2Local.clear();
			polyIdxGlobal2Local.clear();

			//1. index of patch cells in increasing numerical order
			for (int pid = 0; pid < patchLabel.size(); pid++) {
				if (patchLabel[pid] == cluster) {
					patch.push_back(pid);
					polys.push_back(pid);
					patchSize++;
					polyIdxGlobal2Local[pid] = polys.size() - 1;
				}
			}

			//2. index of ribbon cells ordered by adjacent patch cell
			for (auto &pid : patch) {
				for (auto &adj : mesh.adj_p2p(pid)) {
					if (find(polys.begin(), polys.end(), adj) == polys.end()) {
						polys.push_back(adj);
						polyIdxGlobal2Local[adj] = polys.size() - 1;
					}
				}
			}

#ifdef DEBUG
			std::cout << "current cluster: " << cluster << std::endl;
			std::cout << "polys num: " << patchSize << std::endl;
			std::cout << "ribbon size: " << polys.size() - patchSize << std::endl << std::endl;
#endif

			//3. index of faces ordered by relevant cell
			for (auto &pid : polys) {
#ifdef DETAIL
				std::cout << "poly: " << pid << " has faces: ";
#endif
				for (auto &fid : mesh.adj_p2f(pid)) {
					if (std::find(faces.begin(), faces.end(), fid) == faces.end()) {
						faceIdxGlobal2Local[fid] = faces.size();
						faces.push_back(fid);
#ifdef DETAIL
						std::cout << fid << ", ";
#endif
					}
				}
			}

			//4. index of edges ordered by relevant face
			for (auto& fid : faces) {
#ifdef DETAIL
				std::cout << "face: " << fid << " has edges: ";
#endif
				for (auto& eid : mesh.adj_f2e(fid)) {
					if (std::find(edges.begin(), edges.end(), eid) == edges.end()) {
						edgeIdxGlobal2Local[eid] = edges.size();
						edges.push_back(eid);
#ifdef DETAIL
						std::cout << eid << ", ";
#endif
					}
				}
			}

			//5. index of verts ordered by relevant edge
			for (auto& eid : edges) {
#ifdef DETAIL
				std::cout << "edge: " << eid << " has verts: ";
#endif
				for (auto& vid : mesh.adj_e2v(eid)) {
					if (std::find(verts.begin(), verts.end(), vid) == verts.end()) {
						vertIdxGlobal2Local[vid] = verts.size();
						verts.push_back(vid);
						v_pos.push_back(mesh.vert(vid));
#ifdef DETAIL
						std::cout << vid << ", ";
#endif
					}
				}
			}

#ifdef DETAIL
			std::cout << std::endl;
			for (auto &vid : verts) {
				uint idx = vertIdxGlobal2Local[vid];
				vec3d pos = v_pos[idx];
				std::cout << "local vert: " << idx << ", pos:(" << pos.x() << "," << pos.y() << "," << pos.z() << "), global idx: " << vid << std::endl;
			}

			std::cout << std::endl;
			for (auto& eid : edges) {
				uint idx = edgeIdxGlobal2Local[eid];
				uint start = vertIdxGlobal2Local[mesh.edge_vert_id(eid, 0)];
				uint end = vertIdxGlobal2Local[mesh.edge_vert_id(eid, 1)];
				std::cout << "local edge: " << idx << ", from vert: " << start << " to: " << end << ", global idx: " << eid << std::endl;
			}

			std::cout << std::endl;
			for (auto& fid : faces) {
				uint idx = faceIdxGlobal2Local[fid];
				std::cout << "local face: " << idx << ", has edges: ";
				for (auto& eid : mesh.adj_f2e(fid)) {
					std::cout << edgeIdxGlobal2Local[eid] << ", ";
				}
				std::cout << "global idx: " << fid << std::endl;
			}

			std::cout << std::endl;
			for (auto& pid : polys) {
				uint idx = polyIdxGlobal2Local[pid];
				std::cout << "local poly: " << idx << ", has faces: ";
				for (auto& fid : mesh.adj_f2e(pid)) {
					std::cout << faceIdxGlobal2Local[fid] << ", ";
				}
				std::cout << "global idx: " << pid << std::endl;
			}
#endif
			// build all top-down relations
			for (auto& pid : polys) {
				for (auto& fid : mesh.adj_p2f(pid)) {
					if (getPolyFaceSign(pid, fid)) {
						polyFaces.push_back(faceIdxGlobal2Local[fid]);
					}
					else {
						polyFaces.push_back(-faceIdxGlobal2Local[fid] - 1);
					}
				}
			}

			for (auto& fid : faces) {
				for (auto& eid : mesh.adj_f2e(fid)) {
					if (getFaceEdgeSign(fid, eid)) {
						faceEdges.push_back(edgeIdxGlobal2Local[eid]);
					}
					else {
						faceEdges.push_back(-edgeIdxGlobal2Local[eid] - 1);
					}
				}
			}

			for (auto& eid : edges) {
				edgeVerts.emplace_back(vertIdxGlobal2Local[mesh.edge_vert_id(eid, 0)], vertIdxGlobal2Local[mesh.edge_vert_id(eid, 1)]);
			}

			// build all down-top relations
			vertEdgesOffset.reserve(verts.size());
			edgeFacesOffset.reserve(edges.size());
			facePolysOffset.reserve(faces.size());
			
			for (auto& vid : verts) {
				vertEdgesOffset.push_back(mesh.adj_v2e(vid).size());
				for (auto& eid : mesh.adj_v2e(vid)) {
					if (vid == mesh.edge_vert_id(eid, 0)) {
						vertEdges.push_back(edgeIdxGlobal2Local[eid]);
					}
					else {
						vertEdges.push_back(-edgeIdxGlobal2Local[eid] - 1);
					}
				}
			}

			for (auto& eid : edges) {
				edgeFacesOffset.push_back(mesh.adj_e2f(eid).size());
				for (auto& fid : mesh.adj_e2f(eid)) {
					if (getFaceEdgeSign(fid, eid)) {
						edgeFaces.push_back(faceIdxGlobal2Local[fid]);
					}
					else {
						edgeFaces.push_back(-faceIdxGlobal2Local[fid] - 1);
					}
				}
			}

			for (auto& fid : faces) {
				facePolysOffset.push_back(mesh.adj_f2p(fid).size());
				for (auto& pid : mesh.adj_f2p(fid)) {
					if (getPolyFaceSign(pid, fid)) {
						facePolys.push_back(polyIdxGlobal2Local[pid]);
					}
					else {
						facePolys.push_back(-polyIdxGlobal2Local[pid] - 1);
					}
				}
			}

#ifdef OUTPUT
			int pid = 0;
			outfile << "cluster" << cluster << ", size: " << patchSize << std::endl;
			for (; pid < patchSize; pid++) {
				outfile << polys[pid] << '\n';
			}
			outfile << '\n';

			outfile2 << "cluster" << cluster << ", size: " << polys.size() - patchSize << std::endl;
			for (; pid < polys.size(); pid++) {
				outfile2 << polys[pid] << '\n';
			}
			outfile2 << '\n';
#endif
#ifdef OUTPUT_DETAIL
			outfile3 << "cluster" << cluster << std::endl;

			for (auto& t : vertIdxGlobal2Local) {
				outfile3 << "global vert: " << t.first << ", v_pos: (" << v_pos[t.second].x() << "," << v_pos[t.second].y() << "," << v_pos[t.second].z() << "), local idx: " << t.second << '\n';
			}
			outfile3 << '\n';

			for (auto& t : edgeIdxGlobal2Local) {
				outfile3 << "global edge: " << t.first << ", local idx: " << t.second << '\n';
			}
			outfile3 << '\n';

			for (auto& t : faceIdxGlobal2Local) {
				outfile3 << "global face: " << t.first << ", local idx: " << t.second << '\n';
			}
			outfile3 << '\n';

			for (auto& t : polyIdxGlobal2Local) {
				outfile3 << "global poly: " << t.first << ", local idx: " << t.second << '\n';
			}
			outfile3 << '\n';

			int tempIdx = 0;
			for (auto& t : edgeVerts) {
				outfile3 << "local edge: " << tempIdx++ << " has local verts: " << t.x() << " and " << t.y() << '\n';
			}

			for (tempIdx = 0; tempIdx < faceEdges.size(); tempIdx++) {
				if (tempIdx % 4 == 0) {
					outfile3 << "\nlocal face: " << tempIdx / 4 << " has local edges: " << '\n';
				}
				outfile3  << faceEdges[tempIdx] << "  ";
			}

			outfile3 << '\n';
			for (tempIdx = 0; tempIdx < polyFaces.size(); tempIdx++) {
				if (tempIdx % 6 == 0) {
					outfile3 << "\nlocal poly: " << tempIdx / 6 << " has local faces: " << '\n';
				}
				outfile3 << polyFaces[tempIdx] << "  ";
			}
			outfile3 << '\n';
			
			tempIdx = 0;
			int idx = 0;
			for (auto offset : vertEdgesOffset) {
				outfile3 << "\noffset of vert: " << idx++ << " is: " << offset << '\n';
				for (int i = 0; i < offset; i++) {
					outfile3 << vertEdges[tempIdx++] << " ";
				}
			}
			outfile3 << '\n';

			tempIdx = 0;
			idx = 0;
			for (auto offset : edgeFacesOffset) {
				outfile3 << "\noffset of edge: " << idx++ << " is: " << offset << '\n';
				for (int i = 0; i < offset; i++) {
					outfile3 << edgeFaces[tempIdx++] << " ";
				}
			}
			outfile3 << '\n';

			tempIdx = 0;
			idx = 0;
			for (auto offset : facePolysOffset) {
				outfile3 << "\noffset of face: " << idx++ << " is: " << offset << '\n';
				for (int i = 0; i < offset; i++) {
					outfile3 << facePolys[tempIdx++] << " ";
				}
			}
			outfile3 << '\n';
#endif

			patches.emplace_back(std::move(polys), patchSize, std::move(v_pos), std::move(edgeVerts), std::move(faceEdges), std::move(polyFaces),
				std::move(vertEdges), std::move(edgeFaces), std::move(facePolys), std::move(vertEdgesOffset), std::move(edgeFacesOffset), std::move(facePolysOffset));
		}
	}
};

int main() {
	cinolib::Patch patch("D:/data/test/mesh.mesh");
	std::cout << "已读取" << patch.getMesh().vector_polys().size() << "单元体网格" << std::endl;
	patch.patching("D:/data/test/clustered_id.txt", 2);
}


