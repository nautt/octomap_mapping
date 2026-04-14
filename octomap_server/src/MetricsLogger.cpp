#include "octomap_server/MetricsLogger.hpp"

#include <mesher_roi/MeshPoint.h>
#include <octomap/OcTree.h>

#include <chrono>
#include <cmath>
#include <fstream>
#include <unordered_set>

namespace octomap_server {

//Obtener metricas de mesher y escribirlas en json.
static void writeMesherSection(
  std::ofstream & f,
  const std::vector<Clobscode::LeanOctant> & octants,
  const std::vector<Clobscode::MeshPoint> & mesh_points,
  const std::vector<Clobscode::Point3D> & cloud_points,
  double mesher_time_ms)
{

  size_t occupied_count = 0;
  std::unordered_set<unsigned int> covered_indices;
  std::vector<unsigned int> depth_hist(17, 0u);
  std::vector<double> depth_to_edge(17, 0.0);

  for (const auto & oct : octants) {
    const auto & cloud_pts = oct.getContainedCloudPoints();
    if (cloud_pts.empty()) {
      continue;
    }
    ++occupied_count;
    for (unsigned int idx : cloud_pts) {
      covered_indices.insert(idx);
    }
    unsigned short rl = oct.getRefinementLevel();
    if (rl <= 16) {
      depth_hist[rl]++;
      if (depth_to_edge[rl] == 0.0) {
        const auto & pt_idx = oct.getPoints();
        const Clobscode::Point3D & p0 = mesh_points[pt_idx[0]].getPoint();
        const Clobscode::Point3D & p6 = mesh_points[pt_idx[6]].getPoint();
        depth_to_edge[rl] = std::abs(p6.X() - p0.X());
      }
    }
  }

  double coverage = cloud_points.empty() ? 0.0 : static_cast<double>(covered_indices.size()) / cloud_points.size();
  size_t memory_bytes = sizeof(Clobscode::LeanOctant) * octants.size();

  f << "  \"mesher\": {\n";
  f << "    \"time_ms\": " << mesher_time_ms << ",\n";
  f << "    \"occupied_element_count\": " << occupied_count << ",\n";
  f << "    \"total_octant_count\": " << octants.size() << ",\n";
  f << "    \"memory_bytes\": " << memory_bytes << ",\n";
  f << "    \"point_coverage\": " << coverage << ",\n";
  f << "    \"depth_histogram\": [";
  for (int i = 0; i < 17; ++i) {
    f << depth_hist[i];
    if (i < 16) {f << ",";}
  }
  f << "],\n";
  f << "    \"depth_to_edge_m\": [";
  for (int i = 0; i < 17; ++i) {
    f << depth_to_edge[i];
    if (i < 16) {f << ",";}
  }
  f << "]\n";
  f << "  }";
}

//Armar octree temporal sin raycast, calcular metricas y anotarlas en .json.
static void writeOctomapStructuralSection(
  std::ofstream & f,
  const std::vector<Clobscode::Point3D> & points,
  double resolution,
  const std::string & bt_path)
{
  //generar octree temporal
  octomap::OcTree fresh_tree(resolution);

  auto t0 = std::chrono::high_resolution_clock::now(); //timer start
  for (const auto & pt : points) {
    fresh_tree.updateNode(
      octomap::point3d(
        static_cast<float>(pt.X()),
        static_cast<float>(pt.Y()),
        static_cast<float>(pt.Z())),
      true);
  }
  fresh_tree.updateInnerOccupancy();
  auto t1 = std::chrono::high_resolution_clock::now();//timer stop
  double time_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

  size_t occupied_count = 0;
  size_t free_count = 0;
  std::vector<unsigned int> depth_hist(17, 0u);

  for (auto it = fresh_tree.begin_leafs(), end = fresh_tree.end_leafs(); it != end; ++it) {
    unsigned depth = it.getDepth();
    if (fresh_tree.isNodeOccupied(*it)) {
      ++occupied_count;
    } else {
      ++free_count;
    }
    if (depth <= 16) {
      depth_hist[depth]++;
    }
  }

  size_t covered = 0;
  for (const auto & pt : points) {
    octomap::OcTreeNode * node = fresh_tree.search(
      octomap::point3d(
        static_cast<float>(pt.X()),
        static_cast<float>(pt.Y()),
        static_cast<float>(pt.Z())));
    if (node && fresh_tree.isNodeOccupied(node)) {
      ++covered;
    }
  }
  double coverage = points.empty() ? 0.0 : static_cast<double>(covered) / points.size();

  const unsigned int tree_depth = fresh_tree.getTreeDepth();
  std::vector<double> depth_to_edge(17, 0.0);
  for (unsigned int d = 0; d <= tree_depth; ++d) {
    depth_to_edge[d] = fresh_tree.getNodeSize(d);
  }

  fresh_tree.writeBinary(bt_path);

  f << "  \"octomap_snapshot\": {\n";
  f << "    \"source\": \"fresh_single_scan_no_raycasting\",\n";
  f << "    \"time_ms\": " << time_ms << ",\n";
  f << "    \"occupied_element_count\": " << occupied_count << ",\n";
  f << "    \"free_element_count\": " << free_count << ",\n";
  f << "    \"leaf_count\": " << fresh_tree.getNumLeafNodes() << ",\n";
  f << "    \"total_node_count\": " << fresh_tree.size() << ",\n";
  f << "    \"memory_bytes\": " << fresh_tree.memoryUsage() << ",\n";
  f << "    \"point_coverage\": " << coverage << ",\n";
  f << "    \"depth_histogram\": [";
  for (int i = 0; i < 17; ++i) {
    f << depth_hist[i];
    if (i < 16) {f << ",";}
  }
  f << "],\n";
  f << "    \"depth_to_edge_m\": [";
  for (int i = 0; i < 17; ++i) {
    f << depth_to_edge[i];
    if (i < 16) {f << ",";}
  }
  f << "]\n";
  f << "  }";
}

//Obtener metricas de octree nativo de octomap y escribirlas en .json
static void writeOctomapNativeSection(std::ofstream & f, const octomap::OcTree * tree)
{
  size_t occupied_count = 0;
  size_t free_count = 0;
  std::vector<unsigned int> depth_hist(17, 0u);

  for (auto it = tree->begin_leafs(), end = tree->end_leafs(); it != end; ++it) {
    unsigned depth = it.getDepth();
    if (tree->isNodeOccupied(*it)) {
      ++occupied_count;
    } else {
      ++free_count;
    }
    if (depth <= 16) {
      depth_hist[depth]++;
    }
  }

  f << "  \"octomap_native\": {\n";
  f << "    \"source\": \"accumulated_native_with_raycasting\",\n";
  f << "    \"occupied_element_count\": " << occupied_count << ",\n";
  f << "    \"free_element_count\": " << free_count << ",\n";
  f << "    \"leaf_count\": " << tree->getNumLeafNodes() << ",\n";
  f << "    \"total_node_count\": " << tree->size() << ",\n";
  f << "    \"memory_bytes\": " << tree->memoryUsage() << ",\n";
  f << "    \"depth_histogram\": [";
  for (int i = 0; i < 17; ++i) {
    f << depth_hist[i];
    if (i < 16) {f << ",";}
  }
  f << "]\n";
  f << "  }";
}

//mesher + octomap (no raycast)
void MetricsLogger::compute(
  const std::vector<Clobscode::LeanOctant> & octants,
  const std::vector<Clobscode::MeshPoint> & mesh_points,
  const std::vector<Clobscode::Point3D> & cloud_points,
  double resolution,
  double mesher_time_ms,
  const octomap::OcTree * native_tree,
  const std::string & output_path)
{
  bool has_native = native_tree != nullptr && native_tree->size() > 0;

  std::ofstream f(output_path);
  if (!f.is_open()) {
    return;
  }

  std::string bt_path = output_path;
  auto dot = bt_path.rfind('.');
  if (dot != std::string::npos) {
    bt_path = bt_path.substr(0, dot);
  }
  bt_path += ".bt";

  f << "{\n";
  writeMesherSection(f, octants, mesh_points, cloud_points, mesher_time_ms);
  f << ",\n";
  writeOctomapStructuralSection(f, cloud_points, resolution, bt_path);
  if (has_native) {
    f << ",\n";
    writeOctomapNativeSection(f, native_tree);
  }
  f << "\n}\n";
}

//octomap nativo (raycast y origen)
void MetricsLogger::computeNativeOnly(
  const octomap::OcTree * native_tree,
  const std::string & output_path)
{
  if (native_tree == nullptr || native_tree->size() == 0) {
    return;
  }

  std::ofstream f(output_path);
  if (!f.is_open()) {
    return;
  }

  std::string bt_path = output_path;
  auto dot = bt_path.rfind('.');
  if (dot != std::string::npos) {
    bt_path = bt_path.substr(0, dot);
  }
  bt_path += ".bt";
  const_cast<octomap::OcTree *>(native_tree)->writeBinary(bt_path);

  f << "{\n";
  writeOctomapNativeSection(f, native_tree);
  f << "\n}\n";
}

}
