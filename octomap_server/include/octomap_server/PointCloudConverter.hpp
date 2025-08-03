#include "mesher_roi/Point3D.h"
#include <vector>
#include <pcl-1.14/pcl/point_cloud.h>
// #include "pcl_conversions/pcl_conversions.h"
// #include <pcl/segmentation/sac_segmentation.h>
// #include <pcl/io/pcd_io.h>
// #include <pcl/filters/extract_indices.h>
// #include <pcl/filters/passthrough.h>
// #include "rclcpp/rclcpp.hpp"
 
namespace PointCloudConverter {

    /**
     * Convierte una PCLPointCloud (pcl::PointCloud<PCLPoint>)
     * a un std::vector de Clobscode::Point3D con las mismas coordenadas.
     * @param pcl_cloud La nube de entrada, cada punto debe tener campos x, y, z.
     * @return Vector de Point3D con una copia de los puntos.
     */

    //using PCLPointCloud = pcl::PointCloud<pcl::PointXYZ>;
    template <typename PointT>
    inline std::vector<Clobscode::Point3D> toPoint3D(const pcl::PointCloud<PointT>& pcl_cloud) {
        std::vector<Clobscode::Point3D> result;
        result.reserve(pcl_cloud.size());
        for (typename pcl::PointCloud<PointT>::const_iterator it = pcl_cloud.begin(); it != pcl_cloud.end(); ++it) {
            result.emplace_back(it->x, it->y, it->z);
        }
        return result;
    }
}


