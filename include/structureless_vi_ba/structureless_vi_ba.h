#ifndef STRUCTURELESS_VI_BA_H
#define STRUCTURELESS_VI_BA_H

#include "initMethod/drtVioInit.h"
#include "factor/imuIntegFactor.h"

#include <ceres/ceres.h>

class StructurelessVIBA
{
public:
    explicit StructurelessVIBA(DRT::drtVioInit::Ptr drt_vio_init_ptr);

    bool optimize();

private:
    void states_to_double_array();

private:
    DRT::drtVioInit::Ptr drt_vio_init_ptr_;
    
    constexpr int SIZE_POSE {6};
    constexpr int SIZE_SPEEDBIAS {9};

    double para_pose[WINDOW_SIZE][SIZE_POSE];
    double para_speed_bias[WINDOW_SIZE][SIZE_SPEEDBIAS];
};

#endif // STRUCTURELESS_VI_BA_H