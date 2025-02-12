#include "module_base/matrix.h"
#include "module_basis/module_pw/pw_basis_k.h"
#include "module_hamilt_pw/hamilt_pwdft/parallel_grid.h"
#include "module_cell/unitcell.h"
#include "module_cell/klist.h"

#include <complex>

namespace ModuleIO

{
void cal_ldos(const ModulePW::PW_Basis_K* basis,
              const psi::Psi<std::complex<double>>* psi,
              const Parallel_Grid& pgrid,
              const UnitCell& ucell,
              const ModuleBase::matrix& ekb,
              const std::vector<double>& wk,
              const int nks,
              const double eferm,
              const double e_max);
}
