#include "structureless_vi_ba/structureless_vi_ba.h"

using namespace vio;

StructurelessVIBA::StructurelessVIBA(DRT::drtVioInit::Ptr drt_vio_init_ptr) :
                                    drt_vio_init_ptr_{drt_vio_init_ptr}
{

}

bool StructurelessVIBA::optimize()
{
    // Create ceres problem and loss function
    ceres::Problem problem;
    // ceres::LossFunction *loss_function = new ceres::HuberLoss(1.0);

    states_to_double_array();

    // Add preintegration factors for each keyframe inside sliding window
    for (size_t i = 0; i < drt_vio_init_ptr_->local_active_frames.size(); i++)
    {
        ceres::LocalParameterization* pose_param = new PoseSO3LocalParameterization();
        problem.AddParameterBlock(para_pose[i], SIZE_PARAMETERIZATION::SIZE_POSE, pose_param);
        problem.AddParameterBlock(para_speed_bias[i], SIZE_PARAMETERIZATION::SIZE_SPEEDBIAS);
    }

    for (size_t i = 0; i < drt_vio_init_ptr_->local_active_frames.size() - 1; i++)
    {
        size_t j = i + 1;

        if (drt_vio_init_ptr_->imu_meas[i].sum_dt_ > 10)
            continue;

        ImuIntegFactor* imu_factor = new ImuIntegFactor(&drt_vio_init_ptr_->imu_meas[i]);
        problem.AddResidualBlock(imu_factor, NULL, para_pose[i], para_speed_bias[i], para_pose[j], para_speed_bias[j]);
    }

    // Add epipolar constraint factors
    for (size_t i = 0; i < drt_vio_init_ptr_->int_frameid2_time_frameid.size() - 1; i++)
    {
        auto target1_tid = drt_vio_init_ptr_->int_frameid2_time_frameid.at(i);
        for (size_t j = i + 1; j < drt_vio_init_ptr_->int_frameid2_time_frameid.size(); j++)
        {
            auto target2_tid = drt_vio_init_ptr_->int_frameid2_time_frameid.at(j);  

            for (const auto &pts : drt_vio_init_ptr_->SFMConstruct)
            {
                // if a point is observed by two keyframes
                if (pts.second.obs.find(target1_tid) != pts.second.obs.end() &&
                    pts.second.obs.find(target2_tid) != pts.second.obs.end())
                {
                    Eigen::Vector3d zi = pts.second.obs.at(target1_tid).normalpoint;
                    Eigen::Vector3d zj = pts.second.obs.at(target2_tid).normalpoint;

                    EpipolarConstraintFactor* epipolar_constraint_factor = new EpipolarConstraintFactor(zi, zj, RIC[0], TIC[0]);
                    problem.AddResidualBlock(epipolar_constraint_factor, new ceres::HuberLoss(1.0), para_pose[i], para_pose[j]);
                }
            }
        }
    }

    // Start optimization
    ceres::Solver::Options options;
    options.max_num_iterations = 1000;
    options.gradient_tolerance = 1e-20;
    options.function_tolerance = 1e-20;
    options.parameter_tolerance = 1e-20;
    options.linear_solver_type = ceres::SPARSE_SCHUR;
    options.trust_region_strategy_type = ceres::DOGLEG;
    options.minimizer_progress_to_stdout = false;
    ceres::Solver::Summary summary;

    try
    {
        ceres::Solve(options, &problem, &summary);
    }
    catch (const std::exception& e)
    {
        std::cerr << "[CERES EXCEPTION] " << e.what() << std::endl;
        return false;
    }
    catch (...)
    {
        std::cerr << "[CERES UNKNOWN EXCEPTION]" << std::endl;
        return false;
    }

    if (summary.termination_type != ceres::TerminationType::CONVERGENCE)
    {
        std::cerr << "structureless-vi-ba did not converge" << std::endl;
        return false;
    }

    // update old states with the optimized ones
    double_array_to_states();

    return true;
}

void StructurelessVIBA::states_to_double_array()
{
    for (size_t i = 0; i < WINDOW_SIZE; i++)
    {
        // Convert SO(3) to so(3) for R
        auto R = drt_vio_init_ptr_->rotation[i];
        Sophus::SO3d R_SO3(R);
        Eigen::Vector3d so3_R = R_SO3.log();
        
        // Add pose
        para_pose[i][0] = so3_R[0];
        para_pose[i][1] = so3_R[1];
        para_pose[i][2] = so3_R[2];
        para_pose[i][3] = drt_vio_init_ptr_->position[i][0];
        para_pose[i][4] = drt_vio_init_ptr_->position[i][1];
        para_pose[i][5] = drt_vio_init_ptr_->position[i][2];

        // Add velocity and bias
        para_speed_bias[i][0] = drt_vio_init_ptr_->velocity[i][0];
        para_speed_bias[i][1] = drt_vio_init_ptr_->velocity[i][1];
        para_speed_bias[i][2] = drt_vio_init_ptr_->velocity[i][2];
        para_speed_bias[i][3] = drt_vio_init_ptr_->biasg[0];
        para_speed_bias[i][4] = drt_vio_init_ptr_->biasg[1];
        para_speed_bias[i][5] = drt_vio_init_ptr_->biasg[2];
        para_speed_bias[i][6] = drt_vio_init_ptr_->biasa[0];
        para_speed_bias[i][7] = drt_vio_init_ptr_->biasa[1];
        para_speed_bias[i][8] = drt_vio_init_ptr_->biasa[2];
    }
}

void StructurelessVIBA::double_array_to_states() const
{
    Eigen::Vector3d origin_R0 = Utility::R2ypr(drt_vio_init_ptr_->rotation[0]);
    Eigen::Vector3d origin_P0 = drt_vio_init_ptr_->position[0];

    Eigen::Vector3d origin_R00 = Utility::R2ypr(Utility::so3d2SO3d(para_pose[0]));

    double y_diff = origin_R0.x() - origin_R00.x();
    Eigen::Matrix3d rot_diff = Utility::ypr2R(Eigen::Vector3d(y_diff, 0, 0));
    if (abs(abs(origin_R0.y()) - 90) < 1.0 || abs(abs(origin_R00.y()) - 90) < 1.0)
    {
        std::cerr << "Euler singularity!!!" << std::endl;
        rot_diff = drt_vio_init_ptr_->rotation[0] * Utility::so3d2SO3d(para_pose[0]).transpose();
    }

    for (size_t i = 0; i < WINDOW_SIZE; i++)
    {
        drt_vio_init_ptr_->rotation[i] = rot_diff * Utility::so3d2SO3d(para_pose[i]);

        drt_vio_init_ptr_->position[i] = rot_diff * Eigen::Vector3d(para_pose[i][3] - para_pose[0][3],
                                                                para_pose[i][4] - para_pose[0][4],
                                                                para_pose[i][5] - para_pose[0][5]) + origin_P0;
        
        drt_vio_init_ptr_->velocity[i] = rot_diff * Eigen::Vector3d(para_speed_bias[i][0],
                                                                para_speed_bias[i][1],
                                                                para_speed_bias[i][2]);
    }
    drt_vio_init_ptr_->biasg = Eigen::Vector3d(
        para_speed_bias[0][3], para_speed_bias[0][4], para_speed_bias[0][5]);

    drt_vio_init_ptr_->biasa = Eigen::Vector3d(
        para_speed_bias[0][6], para_speed_bias[0][7], para_speed_bias[0][8]);
}