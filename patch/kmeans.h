#pragma once
#include <cinolib/meshes/meshes.h>
#include <vector>
#include <random>
#include <limits>
#include <numeric>
#include <algorithm>
#include <cstdint>
#include <iostream>

namespace cinolib
{
    // ------- 工具函数 -------
    inline double sqr_dist(const vec3d& a, const vec3d& b);

    // k-means++ 初始化：从 points 选 k 个初始质心
    std::vector<vec3d> kmeans_init(const std::vector<vec3d>& points, int k, std::mt19937& rng);

    // ------- 核心：对“单元质心”做 K-Means，返回 labels[pid] -------
    std::vector<int> KMeansOnPolyCenters(const Hexmesh<>& mesh, int k, int max_iters = 100, uint32_t random_seed = 42, bool random_init = true, double rel_move_eps = 1e-6);
}