#ifndef POSE_SO3_LOCAL_PARAMETRIZATION_H
#define POSE_SO3_LOCAL_PARAMETRIZATION_H

#include <ceres/ceres.h>
#include <Eigen/Core>
#include <sophus/se3.hpp>
#include <factor/se3.hpp>

class PoseSO3LocalParameterization : public ceres::LocalParameterization {
public:
    PoseSO3LocalParameterization() = default;
    virtual ~PoseSO3LocalParameterization() = default;

    virtual bool Plus(
        const double* x,
        const double* delta,
        double* x_plus_delta
    ) const override;

    virtual bool ComputeJacobian(
        const double* x,
        double* jacobian
    ) const override;

    virtual int GlobalSize() const override { return 6; }
    virtual int LocalSize()  const override { return 6; }
};

#endif // POSE_SO3_LOCAL_PARAMETRIZATION_H