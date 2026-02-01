#include "factor/pose_so3_local_parametrization.h"

bool PoseSO3LocalParameterization::Plus(
    const double* x,
    const double* delta,
    double* x_plus_delta
) const
{
    // x = [ omega (3), position (3) ]
    Eigen::Map<const Eigen::Vector3d> omega(x);
    Eigen::Map<const Eigen::Vector3d> position(x + 3);

    // delta = [ dtheta (3), dp (3) ]
    Eigen::Map<const Eigen::Vector3d> dtheta(delta);
    Eigen::Map<const Eigen::Vector3d> dp(delta + 3);

    // Convert so3 -> SO3
    std::cout << "omega: " << omega.transpose() << "\n";
    Sophus::SO3d R = Sophus::SO3d::exp(omega);

    // Right-multiplicative update
    Sophus::SO3d R_new = R * Sophus::SO3d::exp(dtheta);

    // Write back
    Eigen::Map<Eigen::Vector3d> omega_new(x_plus_delta);
    Eigen::Map<Eigen::Vector3d> position_new(x_plus_delta + 3);

    omega_new = R_new.log();      // back to so3
    position_new = position + dp;

    return true;
}

bool PoseSO3LocalParameterization::ComputeJacobian(
    const double* /*x*/,
    double* jacobian
) const
{
    // Jacobian of Plus(x, delta) w.r.t delta at delta = 0
    // For so3 ⊕ R3, this is identity
    Eigen::Map<Eigen::Matrix<double, 6, 6, Eigen::RowMajor>> J(jacobian);
    J.setIdentity();
    return true;
}
