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

        ImuIntegFactor *imu_factor = new ImuIntegFactor(drt_vio_init_ptr_->imu_meas[j]);
        problem.AddResidualBlock(imu_factor, NULL, para_pose[i], para_speed_bias[i], para_pose[j], para_speed_bias[j]);
    }

    // Add epipolar constraint factors
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