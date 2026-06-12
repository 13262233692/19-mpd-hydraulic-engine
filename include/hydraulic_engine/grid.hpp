#pragma once

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include <vector>
#include <cmath>

namespace hydraulic_engine {

struct GridFace {
    double z;
    double area;
};

struct GridCell {
    double z_center;
    double z_top;
    double z_bottom;
    double volume;
    double area_top;
    double area_bottom;
    double dz;
};

class Grid {
public:
    Grid(double total_depth, int num_cells, double wellbore_diameter);
    
    int num_cells() const { return num_cells_; }
    int num_faces() const { return num_faces_; }
    
    const GridCell& cell(int i) const { return cells_[i]; }
    const GridFace& face(int i) const { return faces_[i]; }
    
    double total_depth() const { return total_depth_; }
    double wellbore_radius() const { return wellbore_radius_; }
    
    const std::vector<GridCell>& cells() const { return cells_; }
    const std::vector<GridFace>& faces() const { return faces_; }

private:
    int num_cells_;
    int num_faces_;
    double total_depth_;
    double wellbore_radius_;
    
    std::vector<GridCell> cells_;
    std::vector<GridFace> faces_;
};

}
