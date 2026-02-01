#ifndef STRUCTURELESS_VI_BA_H
#define STRUCTURELESS_VI_BA_H

#include "initMethod/drtVioInit.h"
#include "factor/imuIntegFactor.h"
#include "factor/epipolar_contraint_factor.h"
#include "factor/pose_so3_local_parametrization.h"
#include "featureTracker/parameters.h"
#include "utils/eigenUtils.hpp"

#include <ceres/ceres.h>

class StructurelessVIBA
{
public:
    explicit StructurelessVIBA(DRT::drtVioInit::Ptr drt_vio_init_ptr);

    bool optimize();

    using Ptr = std::shared_ptr<StructurelessVIBA>;

private:
    void states_to_double_array();

    void double_array_to_states() const;

private:
    DRT::drtVioInit::Ptr drt_vio_init_ptr_;

    double para_pose[WINDOW_SIZE][SIZE_PARAMETERIZATION::SIZE_POSE];
    double para_speed_bias[WINDOW_SIZE][SIZE_PARAMETERIZATION::SIZE_SPEEDBIAS];
};

#endif // STRUCTURELESS_VI_BA_H