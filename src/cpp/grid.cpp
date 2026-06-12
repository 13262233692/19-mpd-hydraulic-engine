#include "hydraulic_engine/grid.hpp"
#include <cmath>

namespace hydraulic_engine {

Grid::Grid(double total_depth, int num_cells, double wellbore_diameter)
    : num_cells_(num_cells),
      num_faces_(num_cells + 1),
      total_depth_(total_depth),
      wellbore_radius_(wellbore_diameter / 2.0) {
    
    cells_.resize(num_cells_);
    faces_.resize(num_faces_);
    
    double dz = total_depth_ / num_cells_;
    double cross_sectional_area = M_PI * wellbore_radius_ * wellbore_radius_;
    
    for (int i = 0; i < num_faces_; ++i) {
        faces_[i].z = i * dz;
        faces_[i].area = cross_sectional_area;
    }
    
    for (int i = 0; i < num_cells_; ++i) {
        cells_[i].z_bottom = faces_[i].z;
        cells_[i].z_top = faces_[i + 1].z;
        cells_[i].z_center = (cells_[i].z_bottom + cells_[i].z_top) / 2.0;
        cells_[i].dz = dz;
        cells_[i].volume = cross_sectional_area * dz;
        cells_[i].area_bottom = faces_[i].area;
        cells_[i].area_top = faces_[i + 1].area;
    }
}

}
