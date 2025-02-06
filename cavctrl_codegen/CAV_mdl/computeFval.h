//
// File: computeFval.h
//
// MATLAB Coder version            : 24.1
// C/C++ source code generated on  : 08-Jan-2025 15:08:19
//

#ifndef COMPUTEFVAL_H
#define COMPUTEFVAL_H

// Include Files
#include "rtwtypes.h"
#include <cstddef>
#include <cstdlib>

// Type Declarations
struct g_struct_T;

// Function Declarations
namespace coder {
namespace optim {
namespace coder {
namespace qpactiveset {
namespace Objective {
double computeFval(const g_struct_T &obj, double workspace_data[],
                   const double H_data[], const double f_data[],
                   const double x_data[]);

}
} // namespace qpactiveset
} // namespace coder
} // namespace optim
} // namespace coder

#endif
//
// File trailer for computeFval.h
//
// [EOF]
//
