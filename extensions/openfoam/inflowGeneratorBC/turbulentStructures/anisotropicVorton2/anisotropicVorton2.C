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

#include "anisotropicVorton2.H"
#include "transformField.H"
#include "gsl/gsl_multimin.h"
#include <armadillo>

namespace Foam
{
  
anisotropicVorton2::StructureParameters::StructureParameters()
{
}

anisotropicVorton2::StructureParameters::StructureParameters(const dictionary&)
{
}

void anisotropicVorton2::StructureParameters::autoMap
(
    const fvPatchFieldMapper&
)
{
}

//- Reverse map the given fvPatchField onto this fvPatchField
void anisotropicVorton2::StructureParameters::rmap
(
    const fvPatchField<vector>&,
    const labelList&
)
{
}

void anisotropicVorton2::StructureParameters::write(Ostream&) const
{
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

anisotropicVorton2::anisotropicVorton2()
: turbulentStructure(),
  epsilon_(1.0),
  rx_(1.0),
  ry_(1.0),
  rz_(1.0),
  sx_(1.0),
  sy_(1.0),
  sz_(1.0),
  k0_(1.0),
  C1_(1.0)
{}

anisotropicVorton2::anisotropicVorton2(Istream& s)
{
  s>>(*this);
}


anisotropicVorton2::anisotropicVorton2
(
  BoostRandomGen& r, 
  const vector& loc,
  const vector& initialDelta, 
  const vector& v, 
  const vector& L, 
  scalar minL,
  label creaface,
  const symmTensor& R,
  vector yw,
  vector Rpmult
)
: turbulentStructure(r, loc, initialDelta, v, L, minL, creaface, R, yw, Rpmult),
  epsilon_(1.0),
  rx_(1.0),
  ry_(1.0),
  rz_(1.0),
  sx_(1.0),
  sy_(1.0),
  sz_(1.0),
  k0_(1.0),
  C1_(1.0)
{
  double 
    L1 = mag(L1_), 
    L2 = mag(L2_),
    L3 = mag(L3_)
    ;
  
  double 
    R1 = Rpmult.x()* std::max(Rp_[0], 1e-3);
    
  double
    R2 = Rpmult.y()* std::max(1e-3*R1, std::min(Rp_[1], 0.999*R1)), 
    R3 = Rpmult.z()* std::max(1e-3*R1, std::min(Rp_[2], 0.999*R1));

  L3 = std::min(L3, L1*L2*sqrt(R3) / ( L2*sqrt(R1) + L1*sqrt(R2) ));
  L3_ *= L3/mag(L3_);

//   Pout<<"L1->"<<L1<<", L2->"<<L2<<", L3->"<<L3<<", R1->"<<R1<<", R2->"<<R2<<", R3->"<<R3<<endl;
  
  rx_=sqrt(3.)*(
     L1*L1*L3*L3*R2 + L2*L2*( -L3*L3*R1+L1*L1*R3)
	       )/(
    /*sqrt(c)**/pow(L1,3)*L2*sqrt(R3)
	      );
  
//   Pout<<"rx="<<rx_<<endl;
  
  ry_=sqrt(3.)*(
    L1*L1*L3*L3*R2 -L2*L2*(L3*L3*R1+L1*L1*R3)
	      )/(
    /*sqrt(c)**/L1*pow(L2,3)*sqrt(R3)
	      );
  
//   Pout<<"ry="<<ry_<<endl;
  
  double radi=-pow(L1,4)*pow(L3,4)*R2*R2 -pow(L2,4)*pow(L3*L3*R1-L1*L1*R3, 2) +2.*L1*L1*L2*L2*L3*L3*R2*(L3*L3*R1+L1*L1*R3);
//   Pout<<"radi="<<radi<<endl;
  radi=std::max(0., radi);
  rz_=sqrt(3./2./M_PI)*sqrt( 
    radi 
	) / (
    /*sqrt(c)**/L1*L2*L3*sqrt(R3)
	);
      
//   Pout<<"rz="<<rz_<<endl;

  sx_=L1/sqrt(M_PI);
  sy_=L2/sqrt(M_PI);
  sz_=L3/sqrt(M_PI);
    
}

anisotropicVorton2::anisotropicVorton2(const anisotropicVorton2& o)
: turbulentStructure(o),
  epsilon_(o.epsilon_),
  rx_(o.rx_),
  ry_(o.ry_),
  rz_(o.rz_),
  sx_(o.sx_),
  sy_(o.sy_),
  sz_(o.sz_),
  k0_(o.k0_),
  C1_(o.C1_)
{ 
}

vector anisotropicVorton2::fluctuation(const StructureParameters& pa, const vector& x) const
{
  vector delta_x = x - location();
  
  double Xx=delta_x&er1_;
  double Yy=delta_x&er2_;
  double Zz=delta_x&er3_;
  
  double sx2=sx_*sx_;
  double sy2=sy_*sy_;
  double sz2=sz_*sz_;
  double xx2=Xx*Xx;
  double yy2=Yy*Yy;
  double zz2=Zz*Zz;
  
  double e=exp( -0.5 * (xx2/sx2 + yy2/sy2 + zz2/sz2) );
  
  tensor trs = tensor(er1_, er2_, er3_).T();
  
//   Info<<er1_<<"\t"<<er2_<<"\t"<<er3_<<"\t"<<trs<<"\t"<<Rp_<<endl;
//   Info<<(trs&trs.T())<<endl;

  return vector
  (
    transform
    ( 
    
      trs, 
      
      epsilon_ * e * vector
      (
	Yy*( Zz*ry_/sz2 - rz_/sy2 ),
	Xx*( rz_/sx2 - Zz*rx_/sz2 ),
	Xx*Yy*( rx_/sy2 - ry_/sx2 )
      )
      
    )
  );
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<anisotropicVorton2> anisotropicVorton2::New(Istream& s)
{
  return autoPtr<anisotropicVorton2>(new anisotropicVorton2(s));
}


void anisotropicVorton2::randomize(BoostRandomGen& rand)
{
  epsilon_ = 2.0*(rand() - 0.5);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

anisotropicVorton2::~anisotropicVorton2()
{}


autoPtr<anisotropicVorton2> anisotropicVorton2::clone() const
{
  return autoPtr<anisotropicVorton2>
  (
    new anisotropicVorton2(*this)
  );
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void anisotropicVorton2::operator=(const anisotropicVorton2& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn("anisotropicVorton2::operator=(const anisotropicVorton2&)")
            << "Attempted assignment to self"
            << abort(FatalError);
    }

    turbulentStructure::operator=(rhs);
    epsilon_=rhs.epsilon_;
    rx_=rhs.rx_;
    ry_=rhs.ry_;
    rz_=rhs.rz_;
    sx_=rhs.sx_;
    sy_=rhs.sy_;
    sz_=rhs.sz_;
    k0_=rhs.k0_;
    C1_=rhs.C1_;
}

// bool anisotropicVorton2::operator!=(const anisotropicVorton2& o) const
// {
//     return 
// //     turbulentStructure::operator!=(o)
// // ||
//         (location()!=o.location())
//         ||
//         (epsilon_!=o.epsilon_)
// 	||
//         (rx_!=o.rx_)
// 	||
//         (ry_!=o.ry_)
// 	||
//         (rz_!=o.rz_)
// 	||
//         (sx_!=o.sx_)
// 	||
//         (sy_!=o.sy_)
// 	||
//         (sz_!=o.sz_)
// 	||
//         (k0_!=o.k0_)
// 	||
//         (C1_!=o.C1_);
// }

Ostream& operator<<(Ostream& s, const anisotropicVorton2& ht)
{
    s << static_cast<const turbulentStructure&>(ht);
    s<<ht.epsilon_<<endl;
    s<<ht.rx_<<endl;
    s<<ht.ry_<<endl;
    s<<ht.rz_<<endl;
    s<<ht.sx_<<endl;
    s<<ht.sy_<<endl;
    s<<ht.sz_<<endl;
    s<<ht.k0_<<endl;
    s<<ht.C1_<<endl;
    return s;
}

Istream& operator>>(Istream& s, anisotropicVorton2& ht)
{
    s >> static_cast<turbulentStructure&>(ht);
    ht.epsilon_=readScalar(s);
    ht.rx_=readScalar(s);
    ht.ry_=readScalar(s);
    ht.rz_=readScalar(s);
    ht.sx_=readScalar(s);
    ht.sy_=readScalar(s);
    ht.sz_=readScalar(s);
    ht.k0_=readScalar(s);
    ht.C1_=readScalar(s);
    return s;
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// ************************************************************************* //

}
