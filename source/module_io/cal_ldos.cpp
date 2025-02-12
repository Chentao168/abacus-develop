#include "module_io/cal_ldos.h"

#include "module_base/global_variable.h"
#include "module_io/cube_io.h"
#include "module_parameter/parameter.h"

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
              double e_max)
{
    if (PARAM.globalv.dos_setemax)
        e_max = PARAM.inp.dos_emax_ev;
    double emax = eferm + e_max;
    double emin = eferm;
    int nbands = psi->get_nbands();
    int current_spin = 0;
    std::vector<double> ldos(basis->nrxx);
    std::vector<std::complex<double>> wfcr(basis->nrxx);

    if (PARAM.globalv.dos_setemin)
        emin = PARAM.inp.dos_emin_ev;

    for (int ik = 0; ik < nks; ++ik)
    {
        for (int ib = 0; ib < nbands; ib++)
        {
            basis->recip2real(&psi[0](ib, 0), wfcr.data(), ik);
            double eigenval = ekb(ik, ib) * ModuleBase::Ry_to_eV;
            if (eigenval >= emin && eigenval <= emax)
            {
                for (int ir = 0; ir < basis->nrxx; ir++)
                    ldos[ir] += wk[ik] * norm(wfcr[ir]);
            }
        }
    }

    std::stringstream ss;
    ss << PARAM.globalv.global_out_dir << "SPIN" << current_spin + 1 << "_LDOS_" << e_max << "eV" << ".cube";
    ModuleIO::write_vdata_palgrid(pgrid, ldos.data(), current_spin, PARAM.inp.nspin, 0, ss.str(), eferm, &ucell);
}
} // namespace ModuleIO
