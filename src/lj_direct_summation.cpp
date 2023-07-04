#include "lj_direct_summation.h"
#include <cmath>
#include <Eigen/Dense>

inline constexpr double w(double r, double epsilon, double sigma) {
    return 4 * epsilon * (std::pow(sigma / r, 12) - std::pow(sigma / r, 6));
}

inline constexpr double dw_dr(double r, double epsilon, double sigma) {
    return 4 * epsilon *
           (6 * std::pow(sigma, 6) / std::pow(r, 7) -
            12 * std::pow(sigma, 12) / std::pow(r, 13));
}

double lj_direct_summation(Atoms &atoms, double epsilon, double sigma) {
    double epot = 0;
    for (size_t k = 0; k < atoms.nb_atoms(); k++) {
        Eigen::Vector3d force_k;
        force_k.setZero();
        for (size_t i = 0; i < atoms.nb_atoms(); i++) {
            if (k == i) continue;
            Eigen::Vector3d r_ik_vec = atoms.positions.col(i) - atoms.positions.col(k);
            double r_ik = r_ik_vec.norm();
            force_k += dw_dr(r_ik, epsilon, sigma) * r_ik_vec.normalized();
            epot += w(r_ik, epsilon, sigma);
        }
        atoms.forces.col(k) = force_k;
    }
    return epot / 2;
}
