#ifndef OCTOMAP_SERVER_METRICS_LOGGER_HPP
#define OCTOMAP_SERVER_METRICS_LOGGER_HPP

#include <mesher_roi/Mesher.h>
#include <mesher_roi/Point3D.h>
#include <octomap/octomap.h>

#include <string>
#include <vector>

namespace octomap_server {

/**
 * Computes comparative metrics for the mesher and OctoMap and writes them to a JSON file.
 *
 * Three sections are written:
 *   "mesher"             — metrics measured directly from the Clobscode mesher's octants.
 *   "octomap_structural" — fresh single-scan OcTree built from the same point cloud without
 *                          ray-casting, providing a symmetric structural comparison.
 *   "octomap_native"     — the accumulated OcTree built by the native OctoMap pipeline with
 *                          full ray-casting (free-space modelling included). Omitted when the
 *                          tree is empty (no native scans have been inserted).
 */
class MetricsLogger {
public:
  /**
   * @param mesher          The Clobscode mesher after generateMesh() has returned.
   * @param points          The Point3D cloud passed to generateMesh().
   * @param resolution      OctoMap resolution in metres (res_ from OctomapServer).
   * @param mesher_time_ms  Wall-clock time of the generateMesh() call in milliseconds.
   * @param native_tree     Pointer to the accumulated OcTree (octree_.get()); may be null
   *                        or empty, in which case the "octomap_native" section is omitted.
   * @param output_path     Destination path for metrics.json.
   */
  static void compute(
    const Clobscode::Mesher & mesher,
    const std::vector<Clobscode::Point3D> & points,
    double resolution,
    double mesher_time_ms,
    const octomap::OcTree * native_tree,
    const std::string & output_path = "metrics.json");

  /**
   * Write only the octomap_native section to a JSON file.
   * Intended for OCTOMAP_NATIVE mode where no mesher is involved.
   *
   * @param native_tree  Pointer to the accumulated OcTree (must not be null).
   * @param output_path  Destination path for metrics.json.
   */
  static void computeNativeOnly(
    const octomap::OcTree * native_tree,
    const std::string & output_path = "metrics.json");
};

}  // namespace octomap_server

#endif  // OCTOMAP_SERVER_METRICS_LOGGER_HPP
