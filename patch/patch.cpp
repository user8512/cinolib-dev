#include "patch.h"
#include <stack>
#include <algorithm>
#include <random>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdint>
#include <iterator>
#include <cinolib/gl/glcanvas.h>

//#define TEST
//#define DRAW
#define DEBUG
#define OUTPUT
//#define DETAIL
//#define OUTPUT_DETAIL

namespace cinolib {
	Patch::Patch() {
	}

	Patch::Patch(std::string fileName) {
		mesh.load(fileName.c_str());
	}

	Patch::Patch(const std::vector<vec3d>& verts, const std::vector<std::vector<uint>>& polys) {
		mesh = Hexmesh<>(verts, polys);
	}

	Patch::Patch(const std::vector<vec3d>& verts, const std::vector<uint>& polys) {
		mesh = Hexmesh<>(verts, polys);
	}

	bool Patch::getFaceEdgeSign(uint fid, uint eid) {
		vec3d start = mesh.edge_vert(eid, 0);
		vec3d end = mesh.edge_vert(eid, 1);
		std::vector<vec3d> fv = mesh.face_verts(fid);
		if ((fv[0] == start && fv[1] == end) || (fv[1] == start && fv[2] == end) || (fv[2] == start && fv[3] == end) || (fv[3] == start && fv[0] == end)) {
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
		// items num of patch(no ribbon)
		uint patchPolys;
		uint patchFaces;
		uint patchEdges;
		uint patchVerts;
		// pos of each vert
		std::vector<vec3d> v_pos;
		// map global index to local index
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
		std::vector<bool> vertOnSurf;
		std::vector<bool> edgeOnSurf;
		std::vector<bool> faceOnSurf;

		// cluster = patch with 2-ring ribbon
		for (int cluster = 0; cluster < num_clusters; cluster++) {
			// initialize
			patchPolys = 0;
			patchFaces = 0;
			patchEdges = 0;
			patchVerts = 0;
			patch.clear();
			verts.clear();
			edges.clear();
			faces.clear();
			edgeVerts.clear();
			faceEdges.clear();
			polyFaces.clear();
			vertEdges.clear();
			edgeFaces.clear();
			facePolys.clear();
			vertEdgesOffset.clear();
			edgeFacesOffset.clear();
			facePolysOffset.clear();
			vertOnSurf.clear();
			edgeOnSurf.clear();
			faceOnSurf.clear();
			vertIdxGlobal2Local.clear();
			edgeIdxGlobal2Local.clear();
			faceIdxGlobal2Local.clear();
			polyIdxGlobal2Local.clear();

			//1. index of patch cells in increasing numerical order
			for (int pid = 0; pid < patchLabel.size(); pid++) {
				if (patchLabel[pid] == cluster) {
					polyIdxGlobal2Local[pid] = polys.size();
					patch.push_back(pid);
					polys.push_back(pid);
					patchPolys++;
				}
			}

			//2. index of ribbon cells ordered by adjacent patch cell
			std::vector<uint> ring1;
			for (auto& pid : patch) {
				for (auto& adj : mesh.adj_p2p(pid)) {
					if (polyIdxGlobal2Local.find(adj) == polyIdxGlobal2Local.end()) {
						polyIdxGlobal2Local[adj] = polys.size();
						polys.push_back(adj);
						ring1.push_back(adj);
					}
				}
			}
			// 从 1-ring 扩到 2-ring
			for (auto& pid : ring1) {
				for (auto& adj : mesh.adj_p2p(pid)) {
					if (polyIdxGlobal2Local.find(adj) == polyIdxGlobal2Local.end()) {
						polyIdxGlobal2Local[adj] = polys.size();
						polys.push_back(adj);
					}
				}
			}

#ifdef DEBUG
			std::cout << "current cluster: " << cluster << std::endl;
			std::cout << "polys num: " << patchPolys << std::endl;
			std::cout << "ribbon size: " << polys.size() - patchPolys << std::endl;
#endif

			//3. index of faces ordered by relevant cell(first scan to record items num of patch)
			for (auto& pid : patch) {
				for (auto& fid : mesh.adj_p2f(pid)) {
					for (auto& eid : mesh.adj_f2e(fid)) {
						for (auto& vid : mesh.adj_e2v(eid)) {
							if (std::find(verts.begin(), verts.end(), vid) == verts.end()) {
								vertIdxGlobal2Local[vid] = verts.size();
								verts.push_back(vid);
								v_pos.push_back(mesh.vert(vid));
								patchVerts++;
								vertOnSurf.push_back(mesh.vert_is_on_srf(vid));
							}
						}
						if (std::find(edges.begin(), edges.end(), eid) == edges.end()) {
							edgeIdxGlobal2Local[eid] = edges.size();
							edges.push_back(eid);
							patchEdges++;
							edgeOnSurf.push_back(mesh.edge_is_on_srf(eid));
						}
					}
					if (std::find(faces.begin(), faces.end(), fid) == faces.end()) {
						faceIdxGlobal2Local[fid] = faces.size();
						faces.push_back(fid);
						patchFaces++;
						faceOnSurf.push_back(mesh.face_is_on_srf(fid));
					}
				}
			}

			for (auto &pid : polys) {
				for (auto &fid : mesh.adj_p2f(pid)) {
					if (std::find(faces.begin(), faces.end(), fid) == faces.end()) {
						faceIdxGlobal2Local[fid] = faces.size();
						faces.push_back(fid);
						faceOnSurf.push_back(mesh.face_is_on_srf(fid));
					}
				}
			}

			//4. index of edges ordered by relevant face
			for (auto& fid : faces) {
				for (auto& eid : mesh.adj_f2e(fid)) {
					if (std::find(edges.begin(), edges.end(), eid) == edges.end()) {
						edgeIdxGlobal2Local[eid] = edges.size();
						edges.push_back(eid);
						edgeOnSurf.push_back(mesh.edge_is_on_srf(eid));
					}
				}
			}

			//5. index of verts ordered by relevant edge
			for (auto& eid : edges) {
				for (auto& vid : mesh.adj_e2v(eid)) {
					if (std::find(verts.begin(), verts.end(), vid) == verts.end()) {
						vertIdxGlobal2Local[vid] = verts.size();
						verts.push_back(vid);
						v_pos.push_back(mesh.vert(vid));
						vertOnSurf.push_back(mesh.vert_is_on_srf(vid));
					}
				}
			}

#ifdef DETAIL
			/*
			for (auto& pid : polys) {
				std::cout << "\npoly: " << polyIdxGlobal2Local[pid] << " has faces: ";
				for (auto& fid : mesh.adj_p2f(pid)) {
						std::cout << faceIdxGlobal2Local[fid] << ", ";
				}
			}

			for (auto& fid : faces) {
				std::cout << "\nface: " << faceIdxGlobal2Local[fid] << " has edges: ";
				for (auto& eid : mesh.adj_f2e(fid)) {
					std::cout << edgeIdxGlobal2Local[eid] << ", ";
				}
			}

			for (auto& eid : edges) {
				std::cout << "\nedge: " << edgeIdxGlobal2Local[eid] << " has verts: ";
				for (auto& vid : mesh.adj_e2v(eid)) {
					std::cout << vertIdxGlobal2Local[vid] << ", ";
				}
			}
			*/
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
				for (auto& fid : mesh.adj_p2f(pid)) {
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
			vertEdgesOffset.reserve(verts.size() + 1);
			edgeFacesOffset.reserve(edges.size() + 1);
			facePolysOffset.reserve(faces.size() + 1);
			
			uint tempOffset = 0;
			for (auto& vid : verts) {
				vertEdgesOffset.push_back(tempOffset);
				tempOffset += mesh.adj_v2e(vid).size();
				for (auto& eid : mesh.adj_v2e(vid)) {
					if (vid == mesh.edge_vert_id(eid, 0)) {
						vertEdges.push_back(edgeIdxGlobal2Local[eid]);
					}
					else {
						vertEdges.push_back(-edgeIdxGlobal2Local[eid] - 1);
					}
				}
			}
			vertEdgesOffset.push_back(tempOffset);

			tempOffset = 0;
			for (auto& eid : edges) {
				edgeFacesOffset.push_back(tempOffset);
				tempOffset += mesh.adj_e2f(eid).size();
				for (auto& fid : mesh.adj_e2f(eid)) {
					if (getFaceEdgeSign(fid, eid)) {
						edgeFaces.push_back(faceIdxGlobal2Local[fid]);
					}
					else {
						edgeFaces.push_back(-faceIdxGlobal2Local[fid] - 1);
					}
				}
			}
			edgeFacesOffset.push_back(tempOffset);

			tempOffset = 0;
			for (auto& fid : faces) {
				facePolysOffset.push_back(tempOffset);
				tempOffset += mesh.adj_f2p(fid).size();
				for (auto& pid : mesh.adj_f2p(fid)) {
					if (getPolyFaceSign(pid, fid)) {
						facePolys.push_back(polyIdxGlobal2Local[pid]);
					}
					else {
						facePolys.push_back(-polyIdxGlobal2Local[pid] - 1);
					}
				}
			}
			facePolysOffset.push_back(tempOffset);

#ifdef OUTPUT
			int pid = 0;
			outfile << "cluster" << cluster << ", size: " << patchPolys << std::endl;
			for (; pid < patchPolys; pid++) {
				outfile << polys[pid] << '\n';
			}
			outfile << '\n';

			outfile2 << "cluster" << cluster << ", size: " << polys.size() - patchPolys << std::endl;
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
			uint offset = 0;
			for (int ve = 0; ve < vertEdgesOffset.size() - 1; ve++) {
				outfile3 << "\noffset of vert: " << idx++ << " is: " << vertEdgesOffset[ve] << '\n';
				offset = vertEdgesOffset[ve + 1] - vertEdgesOffset[ve];
				for (int i = 0; i < offset; i++) {
					outfile3 << vertEdges[tempIdx++] << " ";
				}
			}
			outfile3 << '\n';

			tempIdx = 0;
			idx = 0;
			for (int ef = 0; ef < edgeFacesOffset.size() - 1; ef++) {
				outfile3 << "\noffset of edge: " << idx++ << " is: " << edgeFacesOffset[ef] << '\n';
				offset = edgeFacesOffset[ef + 1] - edgeFacesOffset[ef];
				for (int i = 0; i < offset; i++) {
					outfile3 << edgeFaces[tempIdx++] << " ";
				}
			}
			outfile3 << '\n';

			tempIdx = 0;
			idx = 0;
			for (int fp = 0; fp < facePolysOffset.size() - 1; fp++) {
				outfile3 << "\noffset of face: " << idx++ << " is: " << facePolysOffset[fp] << '\n';
				offset = facePolysOffset[fp + 1] - facePolysOffset[fp];
				for (int i = 0; i < offset; i++) {
					outfile3 << facePolys[tempIdx++] << " ";
				}
			}
			outfile3 << '\n';
#endif

			patches.emplace_back(std::move(polys), patchPolys, patchFaces, patchEdges, patchVerts, std::move(v_pos), std::move(edgeVerts), std::move(faceEdges), std::move(polyFaces),
				std::move(vertEdges), std::move(edgeFaces), std::move(facePolys), std::move(vertEdgesOffset), std::move(edgeFacesOffset), std::move(facePolysOffset),
				std::move(vertOnSurf), std::move(edgeOnSurf), std::move(faceOnSurf));
		}
	}

	void Patch::subdiv() {
		std::vector<vec3d> pos;
		std::vector<uint> polys;
		int temp = 0;
		for (singlePatch patch : patches) {
			std::cout << "subdiv patch: " << temp++ << std::endl;
			patch.subdiv(pos, polys);
		}
		DrawableHexmesh<> newMesh(pos, polys);
		GLcanvas gui;
		gui.push(&newMesh);
		gui.launch();
		newMesh.save("D:/data/clustered_hexa/result.mesh");
	}

	int count = 0;

	void Patch::singlePatch::subdiv(std::vector<vec3d>& pos, std::vector<uint>& polys) {
#ifdef OUTPUT
		std::ofstream outfile("D:/data/clustered_hexa/test/subdiv.txt", std::ios::app);
		if (!outfile) {
			std::cerr << "无法打开文件" << std::endl;
		}
#endif

		int ne = edgeVerts.size();
		int nf = faceEdges.size() / 4;
		int np = polyFaces.size() / 6;

		//新的几何点
		std::vector<vec3d> newFaceVerts, newEdgeVerts, newVertVerts;

		//各中点（中间量）
		std::vector<vec3d> polyCentroids, faceCentroids, edgeCentroids;

		polyCentroids.reserve(np);
		faceCentroids.reserve(nf);
		edgeCentroids.reserve(ne);

		newFaceVerts.reserve(patchFaces);
		newEdgeVerts.reserve(patchEdges);
		newVertVerts.reserve(patchVerts);

		std::vector<uint> tempVerts;

		//求体的中点
		for (int p = 0; p < np; p++) {
			tempVerts.clear();
			for (int faceOff = 0; faceOff < 6; faceOff++) {
				int f = polyFaces[p * 6 + faceOff];
				if (f < 0) {
					f = -f - 1;
				}
				for (int edgeOff = 0; edgeOff < 4; edgeOff++) {
					int e = faceEdges[f * 4 + edgeOff];
					if (e < 0) {
						e = -e - 1;
					}
					auto start = edgeVerts[e].x();
					auto end = edgeVerts[e].y();
					if (find(tempVerts.begin(), tempVerts.end(), start) == tempVerts.end()) {
						tempVerts.push_back(start);
					}
					if (find(tempVerts.begin(), tempVerts.end(), end) == tempVerts.end()) {
						tempVerts.push_back(end);
					}
				}
			}
			vec3d PolyCentroid(0, 0, 0);
			for (uint v : tempVerts) {
				PolyCentroid += vertsPos[v];
			}
			PolyCentroid /= 8;
			polyCentroids.push_back(PolyCentroid);
		}

		//求面的中点
		for (int f = 0; f < nf; f++) {
			tempVerts.clear();
			for (int edgeOff = 0; edgeOff < 4; edgeOff++) {
				int e = faceEdges[f * 4 + edgeOff];
				if (e < 0) {
					e = -e - 1;
				}
				auto start = edgeVerts[e].x();
				auto end = edgeVerts[e].y();
				if (find(tempVerts.begin(), tempVerts.end(), start) == tempVerts.end()) {
					tempVerts.push_back(start);
				}
				if (find(tempVerts.begin(), tempVerts.end(), end) == tempVerts.end()) {
					tempVerts.push_back(end);
				}
			}
			vec3d faceCentroid(0, 0, 0);
			for (uint v : tempVerts) {
				faceCentroid += vertsPos[v];
			}
			faceCentroid /= 4;
			faceCentroids.push_back(faceCentroid);
		}

		//求边的中点
		for (int e = 0; e < ne; e++) {
			auto start = edgeVerts[e].x();
			auto end = edgeVerts[e].y();
			vec3d edgeCentroid(0, 0, 0);
			edgeCentroid += vertsPos[start];
			edgeCentroid += vertsPos[end];
			edgeCentroid /= 2;
			edgeCentroids.push_back(edgeCentroid);
		}

		//求新体点
		std::vector<vec3d> newPolyVerts(polyCentroids.begin(), polyCentroids.begin() + patchPolys);

		//求新面点
		for (int f = 0; f < patchFaces; f++) {
			if (!faceOnSurf[f]) {
				uint offset = facePolysOffset[f];
				int p0 = facePolys[offset];
				if (p0 < 0) {
					p0 = -p0 - 1;
				}
				int p1 = facePolys[offset + 1];
				if (p1 < 0) {
					p1 = -p1 - 1;
				}
				vec3d newFaceVert(0, 0, 0);
				newFaceVert += polyCentroids[p0];
				newFaceVert += polyCentroids[p1];
				newFaceVert += (faceCentroids[f] * 2);
				newFaceVert /= 4;
				newFaceVerts.push_back(newFaceVert);
			}
			//边界面
			else {
				newFaceVerts.push_back(faceCentroids[f]);
			}
		}

		//求新边点
		std::vector<uint> tempFaces;
		std::vector<uint> tempPolys;

		for (int e = 0; e < patchEdges; e++) {
			tempFaces.clear();
			tempPolys.clear();
			int N = edgeFacesOffset[e + 1] - edgeFacesOffset[e];
			if (!edgeOnSurf[e]) {
				for (int i = 0; i < N; i++) {
					int face = edgeFaces[edgeFacesOffset[e] + i];
					if (face < 0) {
						face = -face - 1;
					}
					tempFaces.push_back(face);
					int faceN = facePolysOffset[face + 1] - facePolysOffset[face];
					for (int j = 0; j < faceN; j++) {
						int poly = facePolys[facePolysOffset[face] + j];
						if (poly < 0) {
							poly = -poly - 1;
						}
						if (find(tempPolys.begin(), tempPolys.end(), poly) == tempPolys.end()) {
							tempPolys.push_back(poly);
						}
					}
				}
				vec3d faceAvg(0, 0, 0);
				for (auto& face : tempFaces) {
					faceAvg += faceCentroids[face];
				}
				faceAvg /= tempFaces.size();
				vec3d polyAvg(0, 0, 0);
				for (auto& poly : tempPolys) {
					polyAvg += polyCentroids[poly];
				}
				polyAvg /= tempPolys.size();
				vec3d newEdgeVert(0, 0, 0);
				newEdgeVert += polyAvg;
				newEdgeVert += (faceAvg * 2);
				newEdgeVert += (edgeCentroids[e] * (N - 3));
				newEdgeVert /= N;
				newEdgeVerts.push_back(newEdgeVert);
			}
			else {
				vec3d faceAvg(0, 0, 0);
				for (int i = 0; i < N; i++) {
					int face = edgeFaces[edgeFacesOffset[e] + i];
					if (face < 0) {
						face = -face - 1;
					}
					if (faceOnSurf[face]) {
						tempFaces.push_back(face);
					}
				}
				for (auto& face : tempFaces) {
					faceAvg += faceCentroids[face];
				}
				faceAvg /= tempFaces.size();
				vec3d newEdgeVert(0, 0, 0);
				newEdgeVert += faceAvg;
				newEdgeVert += edgeCentroids[e];
				newEdgeVert /= 2;
				newEdgeVerts.push_back(newEdgeVert);
			}
		}

		//求新点点
		std::vector<uint> tempEdges;

		for (int v = 0; v < patchVerts; v++) {
			tempEdges.clear();
			tempFaces.clear();
			tempPolys.clear();
			int N = vertEdgesOffset[v + 1] - vertEdgesOffset[v];
			if (!vertOnSurf[v]) {
				for (int i = 0; i < N; i++) {
					int edge = vertEdges[vertEdgesOffset[v] + i];
					if (edge < 0) {
						edge = -edge - 1;
					}
					tempEdges.push_back(edge);
					int edgeN = edgeFacesOffset[edge + 1] - edgeFacesOffset[edge];
					for (int j = 0; j < edgeN; j++) {
						int face = edgeFaces[edgeFacesOffset[edge] + j];
						if (face < 0) {
							face = -face - 1;
						}
						if (find(tempFaces.begin(), tempFaces.end(), face) == tempFaces.end()) {
							tempFaces.push_back(face);
							int faceN = facePolysOffset[face + 1] - facePolysOffset[face];
							for (int k = 0; k < faceN; k++) {
								int poly = facePolys[facePolysOffset[face] + k];
								if (poly < 0) {
									poly = -poly - 1;
								}
								if (find(tempPolys.begin(), tempPolys.end(), poly) == tempPolys.end()) {
									tempPolys.push_back(poly);
								}
							}
						}
					}
				}
				vec3d edgeAvg(0, 0, 0);
				for (auto& edge : tempEdges) {
					edgeAvg += edgeCentroids[edge];
				}
				edgeAvg /= tempEdges.size();
				vec3d faceAvg(0, 0, 0);
				for (auto& face : tempFaces) {
					faceAvg += faceCentroids[face];
				}
				faceAvg /= tempFaces.size();
				vec3d polyAvg(0, 0, 0);
				for (auto& poly : tempPolys) {
					polyAvg += polyCentroids[poly];
				}
				polyAvg /= tempPolys.size();
				vec3d newVertVert(0, 0, 0);
				newVertVert += polyAvg;
				newVertVert += (faceAvg * 3);
				newVertVert += (edgeAvg * 3);
				newVertVert += vertsPos[v];
				newVertVert /= 8;
				newVertVerts.push_back(newVertVert);
			}
			else {
				vec3d edgeAvg(0, 0, 0);
				vec3d faceAvg(0, 0, 0);
				for (int i = 0; i < N; i++) {
					int edge = vertEdges[vertEdgesOffset[v] + i];
					if (edge < 0) {
						edge = -edge - 1;
					}
					if (edgeOnSurf[edge]) {
						tempEdges.push_back(edge);
						int edgeN = edgeFacesOffset[edge + 1] - edgeFacesOffset[edge];
						for (int j = 0; j < edgeN; j++) {
							int face = edgeFaces[edgeFacesOffset[edge] + j];
							if (face < 0) {
								face = -face - 1;
							}
							if (faceOnSurf[face] && find(tempFaces.begin(), tempFaces.end(), face) == tempFaces.end()) {
								tempFaces.push_back(face);
							}
						}
					}
				}
				for (auto& edge : tempEdges) {
					edgeAvg += edgeCentroids[edge];
				}
				edgeAvg /= tempEdges.size();
				for (auto& face : tempFaces) {
					faceAvg += faceCentroids[face];
				}
				faceAvg /= tempFaces.size();
				int n = tempEdges.size();
				vec3d newVertVert(0, 0, 0);
				newVertVert += faceAvg;
				newVertVert += (edgeAvg * 2);
				newVertVert += (vertsPos[v] * (n - 3));
				newVertVert /= n;
				newVertVerts.push_back(newVertVert);
			}
		}

#ifdef OUTPUT
		outfile << "cluster: " << count++ << std::endl;
		outfile << "new poly verts: " << std::endl;
		for (auto& v : newPolyVerts) {
			outfile << "(" << v.x() << ", " << v.y() << ", " << v.z() << ")" << std::endl;
		}
		outfile << std::endl << "new face verts: " << std::endl;
		for (auto& v : newFaceVerts) {
			outfile << "(" << v.x() << ", " << v.y() << ", " << v.z() << ")" << std::endl;
		}
		outfile << std::endl << "new edge verts: " << std::endl;
		for (auto& v : newEdgeVerts) {
			outfile << "(" << v.x() << ", " << v.y() << ", " << v.z() << ")" << std::endl;
		}
		outfile << std::endl << "new vert verts: " << std::endl;
		for (auto& v : newVertVerts) {
			outfile << "(" << v.x() << ", " << v.y() << ", " << v.z() << ")" << std::endl;
		}
#endif

		uint pvOffset = pos.size();
		uint fvOffset = pvOffset + newPolyVerts.size();
		uint evOffset = fvOffset + newFaceVerts.size();
		uint vvOffset = evOffset + newEdgeVerts.size();

		for (auto& v : newPolyVerts) {
			pos.push_back(v);
		}
		for (auto& v : newFaceVerts) {
			pos.push_back(v);
		}
		for (auto& v : newEdgeVerts) {
			pos.push_back(v);
		}
		for (auto& v : newVertVerts) {
			pos.push_back(v);
		}

#ifdef DETAIL
		std::cout << "pvOffset: " << pvOffset << ", fvOffset: " << fvOffset << ", evOffset: " << evOffset << ", vvOffset: " << vvOffset << std::endl;
#endif

		std::vector<uint> FV(6);
		std::vector<uint> EV(12);
		std::vector<uint> VV(8);
		bool f1Reverse = false;
		bool f2Reverse = false;
		// 求所有新的体的拓扑
		for (int p = 0; p < patchPolys; p++) {
			// 1.找到两个相对的面
			int f1, f2;
			f1 = polyFaces[p * 6];
			if (f1 < 0) {
				f1 = -f1 - 1;
				f1Reverse = true;
			}
			tempEdges.clear();
			for (int edgeOff = 0; edgeOff < 4; edgeOff++) {
				int e = faceEdges[f1 * 4 + edgeOff];
				if (e < 0) {
					e = -e - 1;
				}
				tempEdges.push_back(e);
			}

			for (int faceOff = 1; faceOff < 6; faceOff++) {
				bool fReverse = false;
				bool flag = true;
				int f = polyFaces[p * 6 + faceOff];
				if (f < 0) {
					f = -f - 1;
					fReverse = true;
				}
				// 检测是否有重复边
				for (int edgeOff = 0; edgeOff < 4; edgeOff++) {
					int e = faceEdges[f * 4 + edgeOff];
					if (e < 0) {
						e = -e - 1;
					}
					if (find(tempEdges.begin(), tempEdges.end(), e) != tempEdges.end()) {
						flag = false;
					}
					else {

					}
				}
				if (flag) {
					f2 = f;
					f2Reverse = fReverse;
					break;
				}
			}
			// f1,f2是相对的面，利用这两个面构建新的八个细分体

			// 2.记录属于当前体的所有点边面
			tempVerts.clear();
			tempEdges.clear();
			tempFaces.clear();
			for (int faceOff = 0; faceOff < 6; faceOff++) {
				int f = polyFaces[p * 6 + faceOff];
				if (f < 0) {
					f = -f - 1;
				}
				if (find(tempFaces.begin(), tempFaces.end(), f) == tempFaces.end()) {
					tempFaces.push_back(f);
				}
				for (int edgeOff = 0; edgeOff < 4; edgeOff++) {
					int e = faceEdges[f * 4 + edgeOff];
					if (e < 0) {
						e = -e - 1;
					}
					if (find(tempEdges.begin(), tempEdges.end(), e) == tempEdges.end()) {
						tempEdges.push_back(e);
					}
					auto start = edgeVerts[e].x();
					auto end = edgeVerts[e].y();
					if (find(tempVerts.begin(), tempVerts.end(), start) == tempVerts.end()) {
						tempVerts.push_back(start);
					}
					if (find(tempVerts.begin(), tempVerts.end(), end) == tempVerts.end()) {
						tempVerts.push_back(end);
					}
				}
			}

			// 3.对点边面进行排序
			FV.clear();
			EV.clear();
			VV.clear();

			bool edgeReverse = false;
			int vStart, vCurrent;
			int eCurrent = faceEdges[f1 * 4];
			if (eCurrent < 0) {
				eCurrent = -eCurrent - 1;
				edgeReverse = true;
			}

			if (f1Reverse ^ edgeReverse) {
				vStart = edgeVerts[eCurrent].y();
				vCurrent = edgeVerts[eCurrent].x();
				VV.push_back(vStart);
				EV.push_back(eCurrent);
			}
			else {
				vStart = edgeVerts[eCurrent].x();
				vCurrent = edgeVerts[eCurrent].y();
				VV.push_back(vStart);
				EV.push_back(eCurrent);
			}

			while (vCurrent != vStart) {
				VV.push_back(vCurrent);
				for (int eOff = 1; eOff < 4; eOff++) {
					edgeReverse = false;
					eCurrent = faceEdges[f1 * 4 + eOff];
					if (eCurrent < 0) {
						eCurrent = -eCurrent - 1;
						edgeReverse = true;
					}
					if (f1Reverse ^ edgeReverse && edgeVerts[eCurrent].y() == vCurrent) {
						vCurrent = edgeVerts[eCurrent].x();
						EV.push_back(eCurrent);
						break;
					}
					else if (!(f1Reverse ^ edgeReverse) && edgeVerts[eCurrent].x() == vCurrent) {
						vCurrent = edgeVerts[eCurrent].y();
						EV.push_back(eCurrent);
						break;
					}
				}
			}
			// f1的点和边已排序

			for (int i = 0; i < 4; i++) {
				vCurrent = VV[i];
				for (auto& e : tempEdges) {
					if (edgeVerts[e].x() == vCurrent && find(EV.begin(), EV.end(), e) == EV.end()) {
						VV.push_back(edgeVerts[e].y());
						EV.push_back(e);
						break;
					}
					else if (edgeVerts[e].y() == vCurrent && find(EV.begin(), EV.end(), e) == EV.end()) {
						VV.push_back(edgeVerts[e].x());
						EV.push_back(e);
						break;
					}
				}
			}
			// 全部点已排序,侧面边已排序

			FV.push_back(f1);
			FV.push_back(f2);

			for (int i = 0; i < 4; i++) {
				eCurrent = EV[i];
				for (auto& f : tempFaces) {
					if (find(FV.begin(), FV.end(), f) == FV.end()) {
						for (int eOff = 0; eOff < 4; eOff++) {
							int fe = faceEdges[f * 4 + eOff];
							if (fe < 0) {
								fe = -fe - 1;
							}
							if (fe == eCurrent) {
								FV.push_back(f);
								for (int eOff = 0; eOff < 4; eOff++) {
									fe = faceEdges[f * 4 + eOff];
									if (fe < 0) {
										fe = -fe - 1;
									}
									if (find(EV.begin(), EV.end(), fe) == EV.end()) {
										EV.push_back(fe);
										break;
									}
								}
								break;
							}
						}
					}
				}
			}
			// 全部面已排序,全部边已排序

#ifdef DETAIL
			for (auto i : FV)
				std::cout << "FV: " << i << std::endl;
			for (auto i : EV)
				std::cout << "EV: " << i << std::endl;
			for (auto i : VV)
				std::cout << "VV: " << i << std::endl;
#endif

			uint PV = p + pvOffset;
			for (int i = 0; i < FV.size(); i++) {
				FV[i] += fvOffset;
			}
			for (int i = 0; i < EV.size(); i++) {
				EV[i] += evOffset;
			}
			for (int i = 0; i < VV.size(); i++) {
				VV[i] += vvOffset;
			}

			polys.insert(polys.end(), { VV[0], EV[0], FV[0], EV[3], EV[4], FV[2], PV, FV[5] });
			polys.insert(polys.end(), { EV[0], VV[1], EV[1], FV[0], FV[2], EV[5], FV[3], PV });
			polys.insert(polys.end(), { FV[0], EV[1], VV[2], EV[2], PV, FV[3], EV[6], FV[4] });
			polys.insert(polys.end(), { EV[3], FV[0], EV[2], VV[3], FV[5], PV, FV[4], EV[7] });
			polys.insert(polys.end(), { EV[4], FV[2], PV, FV[5], VV[4], EV[8], FV[1], EV[11] });
			polys.insert(polys.end(), { FV[2], EV[5], FV[3], PV, EV[8], VV[5], EV[9], FV[1] });
			polys.insert(polys.end(), { PV, FV[3], EV[6], FV[4], FV[1], EV[9], VV[6], EV[10] });
			polys.insert(polys.end(), { FV[5], PV, FV[4], EV[7], EV[11], FV[1], EV[10], VV[7] });
		}
	}

	inline void PrintVec3d(vec3d& v) {
		std::cout << std::endl << "(" << v.x() << ", " << v.y() << ", " << v.z() << ")" << std::endl;
	}
};

int main() {
#ifdef DRAW
	cinolib::DrawableHexmesh<> newMesh("D:/data/clustered_hexa/fixed_mesh.mesh");
	cinolib::GLcanvas gui;
	gui.push(&newMesh);
	gui.launch();
#else
	#ifdef TEST
		cinolib::Patch patch("D:/data/test/mesh.mesh");
		std::cout << "已读取" << patch.getMesh().vector_polys().size() << "单元体网格" << std::endl;
		patch.patching("D:/data/test/clustered_id.txt", 8);
		patch.subdiv();
	#else
		cinolib::Patch patch("D:/data/clustered_hexa/mesh.mesh");
		std::cout << "已读取" << patch.getMesh().vector_polys().size() << "单元体网格" << std::endl;
		patch.patching("D:/data/clustered_hexa/clustered_id.txt", 16);
		patch.subdiv();
	#endif
#endif
}


