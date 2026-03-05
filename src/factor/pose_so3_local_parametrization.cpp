#include "factor/pose_so3_local_parametrization.h"

bool PoseSO3LocalParameterization::Plus(
    const double* x,
    const double* delta,
    double* x_plus_delta
) const
{
    Eigen::Map<const Eigen::Vector3d> trans(x + 3);
    SE3 se3_delta = SE3::exp(Eigen::Map<const Vector6d>(delta));

    Eigen::Quaterniond quaterd_plus = toQuaterniond(Eigen::Map<const Eigen::Vector3d>(x)) * se3_delta.rotation();
    Eigen::Map<Eigen::Vector3d> angles_plus(x_plus_delta);
    angles_plus = toAngleAxis(quaterd_plus);

    Eigen::Map<Eigen::Vector3d> trans_plus(x_plus_delta + 3);
    trans_plus = se3_delta.rotation() * trans + se3_delta.translation();
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
