//
// File: iterate.h
//
// MATLAB Coder version            : 24.1
// C/C++ source code generated on  : 08-Jan-2025 15:08:19
//

#ifndef ITERATE_H
#define ITERATE_H

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

// Function Declarations
namespace coder {
namespace optim {
namespace coder {
namespace qpactiveset {
void iterate(const double H_data[], const int H_size[2], const double f_data[],
             i_struct_T &solution, h_struct_T &memspace, j_struct_T &workingset,
             e_struct_T &qrmanager, f_struct_T &cholmanager,
             g_struct_T &objective, const char options_SolverName[7],
             double options_StepTolerance, double options_ObjectiveLimit,
             int runTimeOptions_MaxIterations);

}
} // namespace coder
} // namespace optim
} // namespace coder

#endif
//
// File trailer for iterate.h
//
// [EOF]
//
