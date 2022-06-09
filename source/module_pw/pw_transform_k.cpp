#include "fft.h"
#include <complex>
#include "pw_basis_k.h"
#include <cassert>
#include "../module_base/global_function.h"
#include "../module_base/timer.h"
#include "pw_gatherscatter.h"

namespace ModulePW
{

///
/// transform real space to reciprocal space
/// in: (nplane,ny,nx), complex<double> data
/// out: (nz, ns),  complex<double> data
///
void PW_Basis_K:: real2recip(std::complex<double> * in, std::complex<double> * out, int ik)
{
    ModuleBase::timer::tick("PW_Basis_K", "real2recip");

    assert(this->gamma_only == false);
    for(int ir = 0 ; ir < this->nrxx ; ++ir)
    {
        this->ft.aux1[ir] = in[ir];
    }
    this->ft.fftxyfor(ft.aux1,ft.aux1);

    this->gatherp_scatters(this->ft.aux1, this->ft.aux2);
    
    this->ft.fftzfor(ft.aux2,ft.aux1);

    for(int igl = 0 ; igl < this->npwk[ik] ; ++igl)
    {
        out[igl] = this->ft.aux1[this->igl2isz_k[igl+ik*this->npwk_max]] / double(this->bignxyz);
    }
    return;
    ModuleBase::timer::tick("PW_Basis_K", "real2recip");
}

///
/// transform real space to reciprocal space
/// in: (nplane,ny,nx), double data
/// out: (nz, ns), complex<double> data
///
void PW_Basis_K:: real2recip(double * in, std::complex<double> * out, int ik)
{
    ModuleBase::timer::tick("PW_Basis_K", "real2recip_gamma_only");
    assert(this->gamma_only == true);
    // for(int ir = 0 ; ir < this->nrxx ; ++ir)
    // {
    //     this->ft.r_rspace[ir] = in[ir];
    // }
    // r2c in place
    int npy = this->bigny * this->nplane;
    for(int ix = 0 ; ix < this->nx ; ++ix)
    {
        for(int ipy = 0 ; ipy < npy ; ++ipy)
        {
            this->ft.r_rspace[ix*npy*2 + ipy] = in[ix*npy + ipy];
        }
    }

    this->ft.fftxyr2c(ft.r_rspace,ft.aux1);

    this->gatherp_scatters(this->ft.aux1, this->ft.aux2);
    
    this->ft.fftzfor(ft.aux2,ft.aux1);

    for(int igl = 0 ; igl < this->npwk[ik] ; ++igl)
    {
        out[igl] = this->ft.aux1[this->igl2isz_k[igl+ik*this->npwk_max]] / double(this->bignxyz);
    }
    ModuleBase::timer::tick("PW_Basis_K", "real2recip_gamma_only");
    return;
}

///
/// transform reciprocal space to real space
/// in: (nz,ns), complex<double>
/// out: (nplane, ny, nx), complex<double>
///
void PW_Basis_K:: recip2real(std::complex<double> * in, std::complex<double> * out, int ik)
{
    ModuleBase::timer::tick("PW_Basis_K", "recip2real");
    assert(this->gamma_only == false);
    ModuleBase::GlobalFunc::ZEROS(ft.aux1, this->nst * this->nz);

    for(int igl = 0 ; igl < this->npwk[ik] ; ++igl)
    {
        this->ft.aux1[this->igl2isz_k[igl+ik*this->npwk_max]] = in[igl];
    }
    this->ft.fftzbac(ft.aux1, ft.aux2);

    this->gathers_scatterp(this->ft.aux2,this->ft.aux1);

    this->ft.fftxybac(ft.aux1,ft.aux1);
    
    for(int ir = 0 ; ir < this->nrxx ; ++ir)
    {
        out[ir] = this->ft.aux1[ir];
    }
    ModuleBase::timer::tick("PW_Basis_K", "recip2real");

    return;
}

///
/// transform reciprocal space to real space
/// in: (nz,ns), complex<double>
/// out: (nplane, ny, nx), double
///
void PW_Basis_K:: recip2real(std::complex<double> * in, double * out, int ik)
{
    ModuleBase::timer::tick("PW_Basis_K", "recip2real_gamma_only");
    assert(this->gamma_only == true);
    ModuleBase::GlobalFunc::ZEROS(ft.aux1, this->nst * this->nz);

    for(int igl = 0 ; igl < this->npwk[ik] ; ++igl)
    {
        this->ft.aux1[this->igl2isz_k[igl+ik*this->npwk_max]] = in[igl];
    }
   this->ft.fftzbac(ft.aux1, ft.aux2);
    
    this->gathers_scatterp(this->ft.aux2, this->ft.aux1);

    this->ft.fftxyc2r(ft.aux1,ft.r_rspace);

    // for(int ir = 0 ; ir < this->nrxx ; ++ir)
    // {
    //     out[ir] = this->ft.r_rspace[ir] / this->bignxyz;
    // }

    // r2c in place
    int npy = this->bigny * this->nplane;
    for(int ix = 0 ; ix < this->nx ; ++ix)
    {
        for(int ipy = 0 ; ipy < npy ; ++ipy)
        {
            out[ix*npy + ipy] = this->ft.r_rspace[ix*npy*2 + ipy];
        }
    }
    ModuleBase::timer::tick("PW_Basis_K", "recip2real_gamma_only");
    return;
}

#ifdef __MIX_PRECISION
///
/// transform real space to reciprocal space
/// in: (nplane,ny,nx), complex<float> data
/// out: (nz, ns),  complex<float> data
///
void PW_Basis_K:: real2recip(std::complex<float> * in, std::complex<float> * out, int ik)
{
    assert(this->gamma_only == false);
    for(int ir = 0 ; ir < this->nrxx ; ++ir)
    {
        this->ft.auxf1[ir] = in[ir];
    }
    this->ft.fftfxyfor(ft.auxf1,ft.auxf1);

    this->gatherp_scatters(this->ft.auxf1, this->ft.auxf2);
    
    this->ft.fftfzfor(ft.auxf2,ft.auxf1);

    for(int igl = 0 ; igl < this->npwk[ik] ; ++igl)
    {
        out[igl] = this->ft.auxf1[this->igl2isz_k[igl+ik*this->npwk_max]] / float(this->bignxyz);
    }
    return;
}

///
/// transform real space to reciprocal space
/// in: (nplane,ny,nx), float data
/// out: (nz, ns), complex<float> data
///
void PW_Basis_K:: real2recip(float * in, std::complex<float> * out, int ik)
{
    assert(this->gamma_only == true);
    int npy = this->bigny * this->nplane;
    for(int ix = 0 ; ix < this->nx ; ++ix)
    {
        for(int ipy = 0 ; ipy < npy ; ++ipy)
        {
            this->ft.rf_rspace[ix*npy*2 + ipy] = in[ix*npy + ipy];
        }
    }

    this->ft.fftfxyr2c(ft.rf_rspace,ft.auxf1);

    this->gatherp_scatters(this->ft.auxf1, this->ft.auxf2);
    
    this->ft.fftfzfor(ft.auxf2,ft.auxf1);

    for(int igl = 0 ; igl < this->npwk[ik] ; ++igl)
    {
        out[igl] = this->ft.auxf1[this->igl2isz_k[igl+ik*this->npwk_max]] / float(this->bignxyz);
    }
    return;
}

///
/// transform reciprocal space to real space
/// in: (nz,ns), complex<float>
/// out: (nplane, ny, nx), complex<float>
///
void PW_Basis_K:: recip2real(std::complex<float> * in, std::complex<float> * out, int ik)
{
    assert(this->gamma_only == false);
    ModuleBase::GlobalFunc::ZEROS(ft.auxf1, this->nst * this->nz);

    for(int igl = 0 ; igl < this->npwk[ik] ; ++igl)
    {
        this->ft.auxf1[this->igl2isz_k[igl+ik*this->npwk_max]] = in[igl];
    }
    this->ft.fftfzbac(ft.auxf1, ft.auxf2);

    this->gathers_scatterp(this->ft.auxf2,this->ft.auxf1);

    this->ft.fftfxybac(ft.auxf1,ft.auxf1);
    
    for(int ir = 0 ; ir < this->nrxx ; ++ir)
    {
        out[ir] = this->ft.auxf1[ir];
    }

    return;
}

///
/// transform reciprocal space to real space
/// in: (nz,ns), complex<float>
/// out: (nplane, ny, nx), float
///
void PW_Basis_K:: recip2real(std::complex<float> * in, float * out, int ik)
{
    assert(this->gamma_only == true);
    ModuleBase::GlobalFunc::ZEROS(ft.auxf1, this->nst * this->nz);

    for(int igl = 0 ; igl < this->npwk[ik] ; ++igl)
    {
        this->ft.auxf1[this->igl2isz_k[igl+ik*this->npwk_max]] = in[igl];
    }
   this->ft.fftfzbac(ft.auxf1, ft.auxf2);
    
    this->gathers_scatterp(this->ft.auxf2, this->ft.auxf1);

    this->ft.fftfxyc2r(ft.auxf1,ft.rf_rspace);

    int npy = this->bigny * this->nplane;
    for(int ix = 0 ; ix < this->nx ; ++ix)
    {
        for(int ipy = 0 ; ipy < npy ; ++ipy)
        {
            out[ix*npy + ipy] = this->ft.rf_rspace[ix*npy*2 + ipy];
        }
    }
    return;
}

#endif
}