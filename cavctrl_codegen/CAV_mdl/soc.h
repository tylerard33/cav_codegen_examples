//
// File: soc.h
//
// MATLAB Coder version            : 24.1
// C/C++ source code generated on  : 08-Jan-2025 15:08:19
//

#ifndef SOC_H
#define SOC_H

// Include Files
#include "rtwtypes.h"
#include <cstddef>
#include <cstdlib>

// Type Declarations
struct i_struct_T;

struct h_struct_T;

struct j_struct_T;

struct e_struct_T;

struct f_struct_T;

struct g_struct_T;

struct l_struct_T;

// Function Declarations
namespace coder {
namespace optim {
namespace coder {
namespace fminconsqp {
namespace step {
boolean_T soc(const double Hessian_data[], const int Hessian_size[2],
              const double grad_data[], i_struct_T &b_TrialState,
              h_struct_T &memspace, j_struct_T &WorkingSet,
              e_struct_T &b_QRManager, f_struct_T &b_CholManager,
              g_struct_T &QPObjective, const l_struct_T &qpoptions);

}
} // namespace fminconsqp
} // namespace coder
} // namespace optim
} // namespace coder

#endif
//
// File trailer for soc.h
//
// [EOF]
//
