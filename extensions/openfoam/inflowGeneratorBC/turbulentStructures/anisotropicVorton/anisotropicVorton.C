/*
 * This file is part of Insight CAE, a workbench for Computer-Aided Engineering 
 * Copyright (C) 2014  Hannes Kroeger <hannes@kroegeronline.net>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 */

#include "inflowGeneratorBaseFvPatchVectorField.H"
#include "anisotropicVorton.H"
#include "transformField.H"
#include "gsl/gsl_multimin.h"
#include <armadillo>

namespace Foam
{

AnisotropicVorton_Parameters::AnisotropicVorton_Parameters(const inflowInputDataField& ifp, const dictionary& dict)
: turbulentStructure::Parameters(ifp, dict)
{
}

int AnisotropicVorton_Parameters::nParamFields() const
{
  return turbulentStructure::Parameters::nParamFields()+6;
}

string AnisotropicVorton_Parameters::paramFieldName(int i) const
{
  int i0 = i - turbulentStructure::Parameters::nParamFields();
  switch(i0)
    {
    case 0: return "rx"; break;
    case 1: return "ry"; break;
    case 2: return "rz"; break;
    case 3: return "sx"; break;
    case 4: return "sy"; break;
    case 5: return "sz"; break;
    }
  return turbulentStructure::Parameters::paramFieldName(i);
}

turbulentStructure::Parameters::ParamField AnisotropicVorton_Parameters::paramField(int i) const
{
  int i0 = i - turbulentStructure::Parameters::nParamFields();
  switch(i0)
    {
    case 0: return &rx_; break;
    case 1: return &ry_; break;
    case 2: return &rz_; break;
    case 3: return &sx_; break;
    case 4: return &sy_; break;
    case 5: return &sz_; break;
    }
  return turbulentStructure::Parameters::paramField(i);
}


void AnisotropicVorton_Analytic_Parameters::computeParameters()
{
    rx_=scalarField(inputData().num(), 0.0);
    ry_=scalarField(inputData().num(), 0.0);
    rz_=scalarField(inputData().num(), 0.0);
    sx_=scalarField(inputData().num(), 0.0);
    sy_=scalarField(inputData().num(), 0.0);
    sz_=scalarField(inputData().num(), 0.0);
    Lmod_=inputData().L();
    
    for(label i=0; i<inputData().num(); i++)
    {
        scalar c=inputData().c()[i];
        double 
            L1 = inputData().L()[i][0], 
            L2 = inputData().L()[i][1],
            L3 = inputData().L()[i][2]
            ;
        
        double 
            R1 = std::max(this->Rp()[i][0], 1e-3);
            
        double
            R2 = std::max(1e-3*R1, std::min(this->Rp()[i][1], 0.999*R1)), 
            R3 = std::max(1e-3*R1, std::min(this->Rp()[i][2], 0.999*R1));

        // 1.) modify L3   
            
        L3 = L1*L2*sqrt(R3) / ( L2*sqrt(R1) + L1*sqrt(R2) );
        Lmod_[i][2]=L3;
        
        rx_[i]=1.0;
//         ry_[i]=( sqrt(12.)*L2*L3*sqrt(R1) + L1*(sqrt(12.)*L3*sqrt(R2) + L1*rx_[i]) ) / sqr(L2);
        ry_[i]=L1*( sqrt(12.) *L2*R3/sqrt(c*R3)+L1*rx_[i] ) / sqr(L2);
        rz_[i]=L1*( sqrt(12.) *L3*R2/sqrt(c*R2)+L1*rx_[i] ) / sqr(L3);
            

        // 3.) modify L1
        
        //   L1 = L2*L3*sqrt(R1) / ( L3*sqrt(R2) + L2*sqrt(R3) );
        //   L1_ *= L1/mag(L1_);
        //   rz_=1.0;
        //   rx_=/*Rpmult.y() **/	(L3*(2.*sqrt(3.)*L1*sqrt(R2)+L3*rz_))/sqr(L1);
        //   ry_=/*Rpmult.x() **/	(L3*(2.*sqrt(3.)*L2*sqrt(R1)+L3*rz_))/sqr(L2);

        sx_[i]=L1/sqrt(M_PI);
        sy_[i]=L2/sqrt(M_PI);
        sz_[i]=L3/sqrt(M_PI);
        
        Info
            <<R1<<" "
            <<R2<<" "
            <<R3<<" "
            <<L1<<" "
            <<L2<<" "
            <<L3<<": "
            <<rx_[i]<<" "
            <<ry_[i]<<" "
            <<rz_[i]<<" "
            <<sx_[i]<<" "
            <<sy_[i]<<" "
            <<sz_[i]<<endl;
    }
}

AnisotropicVorton_Analytic_Parameters::AnisotropicVorton_Analytic_Parameters(const inflowInputDataField& ifp, const dictionary& dict)
: AnisotropicVorton_Parameters(ifp, dict)
{
    computeParameters();
}
        
void AnisotropicVorton_Analytic_Parameters::write(Ostream&) const
{
//     std::ofstream f("debug_inflo.txt");
//     forAll(inflowPatch(), i)
//     {
//         f << inflowPatch().patch().Cf()[i].y() 
//             << " " << rx_[i]
//             << " " << ry_[i]
//             << " " << rz_[i]
//             << " " << sx_[i]
//             << " " << sy_[i]
//             << " " << sz_[i] <<std::endl;
//     }
}





void AnisotropicVorton_PseudoInv_Parameters::computeParameters()
{
    rx_=scalarField(inputData().num(), 0.0);
    ry_=scalarField(inputData().num(), 0.0);
    rz_=scalarField(inputData().num(), 0.0);
    sx_=scalarField(inputData().num(), 0.0);
    sy_=scalarField(inputData().num(), 0.0);
    sz_=scalarField(inputData().num(), 0.0);

    for(label i=0; i<inputData().num(); i++)
    {
        scalar c=inputData().c()[i];
        
        double 
            L1 = inputData().L()[i][0], 
            L2 = inputData().L()[i][1],
            L3 = inputData().L()[i][2]
            ;
        
        double 
            R1 = std::max(Rp()[i][0], 1e-3);
            
        double
            R2 = std::max(1e-3*R1, std::min(Rp()[i][1], 0.999*R1)), 
            R3 = std::max(1e-3*R1, std::min(Rp()[i][2], 0.999*R1));
        
        /*
        * determine ri via pseudo-inverse
        */
        arma::mat cM, rhs;
        double f=sqrt(c/12.);
        cM << 0. 	    << f*L2/L3 	    << -f*L3/L2 << arma::endr
           << f*L1/L3	<< 0.		    << -f*L3/L1 << arma::endr
           << f*L1/L2	<< -f*L2/L1 	<< 0.	    << arma::endr;

                    
        rhs << sqrt(R1) << arma::endr << sqrt(R2) << arma::endr << sqrt(R3)  << arma::endr;
            
        arma::mat rxryrz = arma::pinv(cM)*rhs;

        rx_[i]=rxryrz(0);
        ry_[i]=rxryrz(1);
        rz_[i]=rxryrz(2);
        sx_[i]=L1/sqrt(M_PI);
        sy_[i]=L2/sqrt(M_PI);
        sz_[i]=L3/sqrt(M_PI);
    }
}


AnisotropicVorton_PseudoInv_Parameters::AnisotropicVorton_PseudoInv_Parameters(const inflowInputDataField& ifp, const dictionary& dict)
: AnisotropicVorton_Parameters(ifp, dict)
{
    computeParameters();
}
        
void AnisotropicVorton_PseudoInv_Parameters::write(Ostream&) const
{
}





void calcS
(
    double Lx, double Ly, double Lz,
    double &sx, double &sy, double &sz
)
{
  double c=sqrt(M_PI);
  
  sx=Lx/c;
  sy=Ly/c;
  sz=Lz/c;
  
  if (mag(sx)<SMALL) sx=SMALL;
  if (mag(sy)<SMALL) sy=SMALL;
  if (mag(sz)<SMALL) sz=SMALL;
}

double anisovf(const gsl_vector *v, void *params)
{
    
    scalar 
        sx, 
        sy, 
        sz,
        rx=gsl_vector_get(v, 0), 
        ry=gsl_vector_get(v, 1), 
        rz=gsl_vector_get(v, 2)
    ;
    
    
    double *par = static_cast<double*>(params);
    scalar R11=par[0], R22=par[1], R33=par[2];
    scalar Lx=par[3], Ly=par[4], Lz=par[5];

    calcS(Lx, Ly, Lz, sx, sy, sz);
  

    vector R=
        vector
        (
        sqr(sqr(Ly)*ry-sqr(Lz)*rz)/sqr(Ly)/sqr(Lz),
        sqr(sqr(Lx)*rx-sqr(Lz)*rz)/sqr(Lx)/sqr(Lz),
        sqr(sqr(Lx)*rx-sqr(Ly)*ry)/sqr(Lx)/sqr(Ly)
        )/12.;
  
  
  return 
      sqr(R[0]-R11) 	 + sqr(R[1]-R22) 	+ sqr(R[2]-R33) ;
}







void AnisotropicVorton_NumOpt_Parameters::computeParameters()
{
    rx_=scalarField(inputData().num(), 0.0);
    ry_=scalarField(inputData().num(), 0.0);
    rz_=scalarField(inputData().num(), 0.0);
    sx_=scalarField(inputData().num(), 0.0);
    sy_=scalarField(inputData().num(), 0.0);
    sz_=scalarField(inputData().num(), 0.0);

    for(label i=0; i<inputData().num(); i++)
    {
        // Numerical determination
        //reinterpret length scales as aligned with eigensystem of Reynolds Stresses
        double 
            Lx = Foam::max(0.0, mag(inputData().L()[i][0])), 
            Lz = mag(inputData().L()[i][2]),
            Ly = Lx*Lz*sqrt(Rp()[i][1]) / ( Lz*sqrt(Rp()[i][0]) + Lx*sqrt(Rp()[i][2]) );


        double par[] = {Rp()[i][0], Rp()[i][1], Rp()[i][2], Lx, Ly, Lz};

        const gsl_multimin_fminimizer_type *T =
            gsl_multimin_fminimizer_nmsimplex2;
        gsl_multimin_fminimizer *s = NULL;
        gsl_vector *ss, *x;
        gsl_multimin_function minex_func;

        int iter = 0;
        int status;
        double size;
        int np=3;

        /* Starting point */
        x = gsl_vector_alloc (np);
        gsl_vector_set_all (x, 1.0);
        gsl_vector_set(x, 0, rx_[i]);
        gsl_vector_set(x, 1, ry_[i]);
        gsl_vector_set(x, 2, rz_[i]);
    //   gsl_vector_set(x, 3, Lx*0.1);
    //   gsl_vector_set(x, 4, Lx*0.1);

        /* Set initial step sizes to 1 */
        ss = gsl_vector_alloc (np);
        gsl_vector_set_all (ss, 0.01);

        /* Initialize method and iterate */
        minex_func.n = np;
        minex_func.f = anisovf;
        minex_func.params = par;

        s = gsl_multimin_fminimizer_alloc (T, np);
        gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

        do
        {
            iter++;
            status = gsl_multimin_fminimizer_iterate(s);

            if (status)
                break;

            size = gsl_multimin_fminimizer_size (s);
            status = gsl_multimin_test_size (size, 1e-6);

        }
        while (status == GSL_CONTINUE && iter < 2000);

        rx_=gsl_vector_get(s->x, 0);
        ry_=gsl_vector_get(s->x, 1);
        rz_=gsl_vector_get(s->x, 2);
    //   double Ly=Foam::max(minL, gsl_vector_get(s->x, 3));
    //   double Lz=Foam::max(minL, gsl_vector_get(s->x, 4));
    //   sy_=gsl_vector_get(s->x, i++);
    //   sz_=gsl_vector_get(s->x, i++);
    //   C1_=gsl_vector_get(s->x, 3);
    //   k0_=gsl_vector_get(s->x, 4);
    //
    //   C1_=std::max(SMALL, C1_);
    //   k0_=std::max(SMALL, k0_);

//         L1_=e1()[i]*Lx;
//         L2_=e2()[i]*Ly;
//         L3_=e3()[i]*Lz;

        calcS(Lx, Ly, Lz, sx_[i], sy_[i], sz_[i]);

    //   scalar err=anisovf(s->x, par);
    //   if (err>0.1)
    //   Info<<"@"<<Rp_<<";"<<mag(L1_)<<" "<<mag(L2_)<<" "<<mag(L3_)<<": \t"
    //     <<rx_<<" "<<ry_<<" "<<rz_<<" / \t"
    //     <<sx_<<" "<<sy_<<" "<<sz_<<" / \t"
    //     <<k0_<<" "<<C1_<<
    //     " \t err="<<err<<" #"<<iter<<endl;

        gsl_vector_free(x);
        gsl_vector_free(ss);
        gsl_multimin_fminimizer_free (s);

//         if (mag(mag(e1()[i])-1.0)>1e-6)
//         {
//             Info<<Rp_<<" er="<<er1_<<er2_<<er3_<<endl;
//         }
    }
}


AnisotropicVorton_NumOpt_Parameters::AnisotropicVorton_NumOpt_Parameters(const inflowInputDataField& ifp, const dictionary& dict)
: AnisotropicVorton_Parameters(ifp, dict)
{
    computeParameters();
}
        
void AnisotropicVorton_NumOpt_Parameters::write(Ostream&) const
{
}

// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// ************************************************************************* //

}
