#ifndef __SIMULATION_UTILS_MPI_H
#define __SIMULATION_UTILS_MPI_H

#include "atoms.h"
#include "domain.h"
#include <argparse/argparse.hpp>
#include <Eigen/Dense>

class Stretcher {
  private:
    size_t interval_;
    double increase_;
    double strain_ = 0;

  public:
    Stretcher(size_t stretch_interval, double length_increase) : interval_(stretch_interval), increase_(length_increase) {}
    ~Stretcher() {}
    void step(Atoms& atoms, Domain& domain, size_t timestep) {
        if (timestep % interval_ == 0) {
            Eigen::Array3d new_length{
                domain.domain_length(0), domain.domain_length(1),
                domain.domain_length(2) + increase_};
            domain.scale(atoms, new_length);
            strain_ += increase_;
        }
    }
    double strain() { return strain_; }
};

Domain init_domain(Atoms& atoms, argparse::ArgumentParser parser) {
    auto ds = parser.get<std::vector<int>>("--domains");
    auto periodic = parser.get<std::vector<int>>("--periodic");
    double shift = parser.get<double>("--shift_atoms");
    bool verbose = parser.get<bool>("--verbose");

    auto max_pos = atoms.positions.rowwise().maxCoeff();
    auto min_pos = atoms.positions.rowwise().minCoeff();
    // bounding box
    Eigen::Array3d bbox{max_pos(0) - min_pos(0),
                        max_pos(1) - min_pos(1),
                        max_pos(2) - min_pos(2)};
    double offset_x = bbox(0) * shift - min_pos(0);
    double offset_y = bbox(1) * shift - min_pos(1);
    double offset_z = bbox(2) * shift - min_pos(2);
    Eigen::Array3d offset{offset_x, offset_y, offset_z};
    for (size_t i = 0; i < 3; i++) {
        if (periodic[i] == 0) {
            atoms.positions.row(i) += offset(i);
        }
    }

    double box_sz_x = periodic[0] == 0 ? bbox(0) + 2 * bbox(0) * shift : max_pos(0) + min_pos(0);
    double box_sz_y = periodic[1] == 0 ? bbox(1) + 2 * bbox(1) * shift : max_pos(1) + min_pos(1);
    double box_sz_z = periodic[2] == 0 ? bbox(2) + 2 * bbox(2) * shift : max_pos(2) + min_pos(2);

    Domain domain(MPI_COMM_WORLD,
        {box_sz_x, box_sz_y, box_sz_z},
        {ds[0], ds[1], ds[2]},
        {periodic[0], periodic[1], periodic[2]},
        verbose
    );
    return domain;
}

double compute_stress(const Domain& domain, const Atoms& atoms, int dim = 2) {
    double left_domain_boundary{domain.coordinate(dim) * domain.domain_length(dim) / domain.decomposition(dim)};
    auto left_mask{atoms.positions.row(dim) < left_domain_boundary};
    double stress = 0;
    for (::Eigen::Index i{0}; i < left_mask.size(); ++i) {
        if (left_mask[i]) {
            stress += atoms.forces(dim, i);
        }
    }
    return stress;
}


#endif // __SIMULATION_UTILS_MPI_H
