#include "structureless_vi_ba.h"

StructurelessVIBA::StructurelessVIBA(DRT::drtVioInit::Ptr drt_vio_init_ptr) :
                                    drt_vio_init_ptr_{drt_vio_init_ptr}
{

}

bool StructurelessVIBA::optimize()
{
    // Create ceres problem and loss function
    ceres::Problem problem;
    ceres::LossFunction *loss_function = new ceres::CauchyLoss(1.0);

    states_to_double_array();

    // Add preintegration factors for each keyframe inside sliding window
    for (size_t i = 0; i < drt_vio_init_ptr_->local_active_frames.size(); i++)
    {
        size_t j = i + 1;

        ImuIntegFactor* imu_factor = new ImuIntegFactor(drt_vio_init_ptr_->imu_meas[j]);
        problem.AddResidualBlock(imu_factor, NULL, para_pose[i], para_speed_bias[i], para_pose[j], para_speed_bias[j]);
    }

    // Add epipolar constraint factors
    for (size_t i = 0; i < drt_vio_init_ptr_->int_frameid2_time_frameid.size() - 1; i++)
    {
        auto target1_tid = int_frameid2_time_frameid.at(i);
        auto target2_tid = int_frameid2_time_frameid.at(i + 1);

        for (const auto &pts : drt_vio_init_ptr_->SFMConstruct)
        {
            // if a point is observed by two keyframes
            if (pts.second.obs.find(target1_tid) != pts.second.obs.end() &&
                pts.second.obs.find(target2_tid) != pts.second.obs.end())
            {
                Eigen::Vector3d zi = pts.second.obs.at(target1_tid).normalpoint;
                Eigen::Vector3d zj = pts.second.obs.at(target2_tid).normalpoint;

                EpipolarConstraintFactor* epipolar_constraint_factor = new EpipolarConstraintFactor(zi, zj, RIC[0], TIC[0]);
                problem.AddResidualBlock(epipolar_constraint_factor, NULL, para_pose[i], para_pose[j]);
            }
        }
    }

    // Start optimization
    ceres::Solver::Options options;
    options.max_num_iterations = 200;
    options.gradient_tolerance = 1e-20;
    options.function_tolerance = 1e-20;
    options.parameter_tolerance = 1e-20;
    options.linear_solver_type = ceres::DENSE_SCHUR;
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
        false;
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

    // get optimized parameters here.
    
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