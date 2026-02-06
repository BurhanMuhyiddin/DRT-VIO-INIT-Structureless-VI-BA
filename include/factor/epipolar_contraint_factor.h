#ifndef EPIPOLAR_CONSTRAINT_FACTOR_H
#define EPIPOLAR_CONSTRAINT_FACTOR_H

#include "utils/eigenUtils.hpp"

#include <ceres/ceres.h>
#include <Eigen/Dense>

namespace vio
{
    class EpipolarConstraintFactor final : public ceres::SizedCostFunction<1, 6, 6>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        EpipolarConstraintFactor(const Eigen::Vector3d &zi,
                                const Eigen::Vector3d &zj,
                                const Eigen::Matrix3d &R_IC,
                                const Eigen::Vector3d &p_IC)
            : zi_(zi),
              zj_(zj),
              R_IC_(R_IC),
              p_IC_(p_IC)
        {
            sqrt_info_ = 458.654 / 1.5;
        }

        bool Evaluate(double const* const* parameters,
                    double* residuals,
                    double** jacobians) const;

    private:
        const Eigen::Vector3d zi_;
        const Eigen::Vector3d zj_;
        Eigen::Matrix3d R_IC_;
        Eigen::Vector3d p_IC_;
        double sqrt_info_;
    };
}

#endif // EPIPOLAR_CONSTRAINT_FACTOR_H