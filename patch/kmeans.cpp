#include "kmeans.h"

namespace cinolib {
    inline double sqr_dist(const vec3d& a, const vec3d& b) {
        const vec3d d = a - b;
        return d.dot(d);
    }

    std::vector<vec3d> kmeans_init(const std::vector<vec3d>& points, int k, std::mt19937& rng) {
        const size_t n = points.size();
        std::uniform_int_distribution<size_t> uni0(0, n - 1);
        std::vector<vec3d> centers;
        centers.reserve(k);
        centers.push_back(points[uni0(rng)]);

        std::vector<double> min_d2(n, std::numeric_limits<double>::infinity());
        for (int c = 1; c < k; ++c) {
            // 更新每个点到最近已选中心的距离
            for (size_t i = 0; i < n; ++i)
            {
                double d2 = sqr_dist(points[i], centers.back());
                if (d2 < min_d2[i]) min_d2[i] = d2;
            }
            // 按距离平方做加权抽样
            double sum = std::accumulate(min_d2.begin(), min_d2.end(), 0.0);
            if (sum == 0.0) { centers.push_back(points[uni0(rng)]); continue; }
            std::uniform_real_distribution<double> uni(0.0, sum);
            double r = uni(rng);
            size_t idx = 0;
            double acc = 0.0;
            for (; idx < n; ++idx)
            {
                acc += min_d2[idx];
                if (acc >= r) break;
            }
            centers.push_back(points[std::min(idx, n - 1)]);
        }
        return centers;
    }

    std::vector<int> KMeansOnPolyCenters(const Hexmesh<>& mesh, int k, int max_iters, uint32_t random_seed, bool random_init, double rel_move_eps) {
        const size_t P = mesh.num_polys();
        std::vector<int> labels(P, -1);
        if (P == 0 || k <= 0) return labels;
        k = std::min<int>(k, (int)P);

        // 1) 收集所有 poly 的质心
        std::vector<vec3d> points;
        points.reserve(P);
        for (uint pid = 0; pid < P; ++pid)
            points.push_back(mesh.poly_centroid(pid));

        // 2) 初始化质心
        std::mt19937 rng(random_seed);
        std::vector<vec3d> centers;
        if (random_init) centers = kmeans_init(points, k, rng);
        else {
            std::vector<size_t> idx(P);
            std::iota(idx.begin(), idx.end(), 0);
            std::shuffle(idx.begin(), idx.end(), rng);
            centers.resize(k);
            for (int c = 0; c < k; ++c) centers[c] = points[idx[c]];
        }

        // 3) 迭代
        double last_total_shift = std::numeric_limits<double>::infinity();
        for (int it = 0; it < max_iters; ++it) {
            // 3.1 赋值
            bool any_change = false;
            for (size_t i = 0; i < P; ++i)
            {
                double best_d = std::numeric_limits<double>::infinity();
                int best_c = -1;
                for (int c = 0; c < k; ++c)
                {
                    double d = sqr_dist(points[i], centers[c]);
                    if (d < best_d) { best_d = d; best_c = c; }
                }
                if (labels[i] != best_c) { labels[i] = best_c; any_change = true; }
            }

            // 3.2 更新质心
            std::vector<vec3d> new_c(k, vec3d(0, 0, 0));
            std::vector<int>   cnt(k, 0);
            for (size_t i = 0; i < P; ++i)
            {
                int c = labels[i];
                new_c[c] += points[i];
                cnt[c] += 1;
            }

            double total_shift = 0.0;
            for (int c = 0; c < k; ++c)
            {
                vec3d oldc = centers[c];
                if (cnt[c] > 0) centers[c] = new_c[c] / double(cnt[c]);
                else {
                    // 空簇：随机重置为某个样本
                    std::uniform_int_distribution<int> uni(0, (int)P - 1);
                    centers[c] = points[uni(rng)];
                }
                total_shift += sqr_dist(oldc, centers[c]);
            }

            // 收敛判据：标签无变化或质心总移动量很小
            if (!any_change) break;
            if (last_total_shift < std::numeric_limits<double>::infinity())
            {
                double rel = total_shift / (last_total_shift + 1e-16);
                if (rel < rel_move_eps) break;
            }
            last_total_shift = total_shift;
        }
        return labels; // labels[i] 对应 pid == i 的单元
    }
}