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
//#define DETAIL

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

		// temporary
		std::vector<uint> patch;
		std::vector<uint> verts;

		// index of polys of patch and ribbon in order
		std::vector<uint> polys;
		// cells num of patch(no ribbon)
		uint patchSize;
		// pos of each vert
		std::vector<vec3d> v_pos;

		std::map<uint, uint> vertIdxGlobal2Local;
		std::map<uint, uint> edgeIdxGlobal2Local;
		std::map<uint, uint> faceIdxGlobal2Local;

		std::vector<vec2d> edgeVerts;
		std::vector<int> faceEdges;
		std::vector<int> cellFaces;
		std::vector<uint> vertEdges;
		std::vector<int> edgeFaces;
		std::vector<int> faceCells;

		for (int cluster = 0; cluster < num_clusters; cluster++) {
			patchSize = 0;
			patch.clear();
			verts.clear();
			vertIdxGlobal2Local.clear();

			//1. index of patch cells in increasing numerical order
			for (int pid = 0; pid < patchLabel.size(); pid++) {
				if (patchLabel[pid] == cluster) {
					patch.push_back(pid);
					polys.push_back(pid);
					patchSize++;
				}
			}


			//2. index of ribbon cells ordered by adjacent patch cell
			for (auto &pid : patch) {
				for (auto &adj : mesh.adj_p2p(pid)) {
					if (find(patch.begin(), patch.end(), adj) == patch.end()) {
						polys.push_back(adj);
					}
				}
			}
#ifdef DEBUG
			std::cout << "current cluster: " << cluster << std::endl;
			std::cout << "polys num: " << patchSize << std::endl;
			std::cout << "ribbon size: " << polys.size() - patchSize << std::endl << std::endl;
#endif
			//3. index of verts ordered by relevant cell
			for (auto &pid : polys) {
#ifdef DETAIL
				std::cout << "poly: " << pid << " has verts: ";
#endif
				for (auto &vid : mesh.poly_verts_id(pid)) {
					if (std::find(verts.begin(), verts.end(), vid) == verts.end()) {
						verts.push_back(vid);
						v_pos.push_back(mesh.vert(vid));
						vertIdxGlobal2Local[vid] = verts.size() - 1;
#ifdef DETAIL
						std::cout << vid << ", ";
#endif
					}
				}
			}

#ifdef DETAIL
			std::cout << std::endl;
			for (int vid = 0; vid < v_pos.size(); vid++) {
				vec3d pos = v_pos[vid];
				std::cout << "local vid: " << vid << ", pos:(" << pos.x() << "," << pos.y() << "," << pos.z() << ")" << std::endl;
			}
			for (auto& pid : patch) {
				for (auto& fid : mesh.adj_p2f(pid)) {
					std::cout << "patch: " << pid << " has face: " << fid << " , CCW: " << getPolyFaceSign(pid, fid) << " with verts : ";
					for (auto& vid : mesh.adj_f2v(fid)) {
						std::cout << vid << ", ";
					}
					std::cout << " and edges : " << std::endl;
					for (auto& eid : mesh.adj_f2e(fid)) {
						std::cout << eid << " from: " << mesh.edge_vert_id(eid, 0) << " to: " << mesh.edge_vert_id(eid, 1) << " , positive: " << getFaceEdgeSign(fid, eid) << std::endl;
					}
				}
			}
#endif
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
			singlePatch newPatch(std::move(polys), patchSize, std::move(v_pos), std::move(edgeVerts), std::move(faceEdges), std::move(cellFaces),
				std::move(vertEdges), std::move(edgeFaces), std::move(faceCells));
		}
	}
};

int main() {
	cinolib::Patch patch("D:/data/clustered_hexa/mesh.mesh");
	std::cout << "已读取" << patch.getMesh().vector_polys().size() << "单元体网格" << std::endl;
	patch.patching("D:/data/clustered_hexa/clustered_id.txt", 16);
}


