#include "epipolar_contraint_factor.h"

namespace vio
{
    bool EpipolarConstraintFactor::Evaluate(double const* const* params,
                    double* residuals,
                    double** jacobians) const override
    {
        Eigen::Vector3d omegai(parameters[0][0], parameters[0][1], parameters[0][2]);
        Eigen::Vector3d Pi(parameters[0][3], parameters[0][4], parameters[0][5]);
        Sophus::SO3d Ri = Sophus::SO3d::exp(omegai);

        Eigen::Vector3d omegaj(parameters[1][0], parameters[1][1], parameters[1][2]);
        Eigen::Vector3d Pj(parameters[1][3], parameters[1][4], parameters[1][5]);
        Sophus::SO3d Rj = Sophus::SO3d::exp(omegai);

        // Calculate A and B vectors eq(14)
        Eigen::Vector3d B = Ri * R_IC_ * zi_;
        Eigen::Vector3d A = Rj * R_IC_ * zj_;

        // Calculate vector t eq(13)
        Eigen::Vector3d p_Ci = pi + Ri * p_IC_;
        Eigen::Vector3d p_Cj = pj + Rj * p_IC_;
        Eigen::Vector3d t = p_Ci - p_Cj;

        double t_norm = t.norm();
        if (t_norm < 1e-8)
        {
            residuals[0] = 0.0;
            return true;
        }

        // Calculate vector C eq(14)
        Eigen::Vector3d C = t / t_norm;

        // Calculate residual eq(14)
        residuals[0] = A.transpose() * (C.cross(B));

        // Calculate analytical Jacobians from eq(15) to eq(19)
        if (jacobians)
        {
            const Eigen::Matrix3d Cx = skewSymmetric(C);
            coonst Eigen::Matrix3d Bx = skewSymmetric(B);

            Eigen::Matrix3d dC_dt = (Eigen::Matrix3d::Identity() / t_norm) - 
                                    (t * t.transpose()) / (t_norm * t_norm * t_norm);

            if (jacobians[0])
            {
                Eigen::Map<Eigen::Matrix<double, 1, 6, Eigen::RowMajor>> Ji(jacobians[0]);
                Ji.setZero();

                Ji.block<1,3>(0, 0) = A.transpose() * Bx * dC_dt;

                Ji.block<1,3>(0, 3) = A.transpose() * Cx * (-Ri * skewSymmetric(R_IC_ * zi_)) + 
                                    A.transpose() * Bx * dC_dt * (-Ri * skewSymmetric(p_IC_));
            }

            if (jacobians[1])
            {
                Eigen::Map<Eigen::Matrix<double, 1, 6, Eigen::RowMajor>> Jj(jacobians[1]);
                Jj.setZero();

                Jj.block<1,3>(0,0) = -A.transpose() * Bx * dC_dt;

                Jj.block<1,3>(0,3) = (-skewSymmetric(Rj * R_IC_ * zj_).transpose() * (C.cross(B))) +
                    A.transpose() * Bx * dC_dt * (Rj * skewSymmetric(p_IC_));
            }
        }

        return true;
    }
}