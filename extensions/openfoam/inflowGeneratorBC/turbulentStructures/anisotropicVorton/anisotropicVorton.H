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

#ifndef anisotropicVorton_H
#define anisotropicVorton_H

#include "turbulentStructure.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// class anisotropicVorton;
// class AnisotropicVorton_Analytic_Parameters
// 
// Ostream& operator<<(Ostream&, const anisotropicVorton&);
// Istream& operator>>(Istream&, anisotropicVorton&);

/*---------------------------------------------------------------------------*\
                           Class anisotropicVorton Declaration
\*---------------------------------------------------------------------------*/


class AnisotropicVorton_Parameters
: public turbulentStructure::Parameters
{
public:
    scalarField rx_, ry_, rz_, sx_, sy_, sz_;

    AnisotropicVorton_Parameters(const inflowInputDataField& ifp, const dictionary&);

    virtual int nParamFields() const;
    virtual string paramFieldName(int i) const;
    virtual ParamField paramField(int i) const;
};


class AnisotropicVorton_Analytic_Parameters
: public AnisotropicVorton_Parameters
{
protected:
    vectorField Lmod_;
    
    void computeParameters();
public:
    AnisotropicVorton_Analytic_Parameters(const inflowInputDataField& ifp, const dictionary&);
    
    inline const vectorField& Lmod() const { return Lmod_; }
    
    void write(Ostream&) const;
};


class AnisotropicVorton_PseudoInv_Parameters
: public AnisotropicVorton_Parameters
{
protected:
    
    void computeParameters();
public:
    AnisotropicVorton_PseudoInv_Parameters(const inflowInputDataField& ifp, const dictionary&);
    
    inline const vectorField& Lmod() const { return inputData().L(); }
    
    void write(Ostream&) const;
};


class AnisotropicVorton_NumOpt_Parameters
: public AnisotropicVorton_Parameters
{
protected:
    
    void computeParameters();
public:
    AnisotropicVorton_NumOpt_Parameters(const inflowInputDataField& ifp, const dictionary&);
    
    inline const vectorField& Lmod() const { return inputData().L(); }

    void write(Ostream&) const;
};




template<class ParameterClass>
class anisotropicVorton
    : public turbulentStructure
{

public:
    typedef ParameterClass Parameters;
    
protected:
    const Parameters& p_;
    scalar epsilon_;

public:

    // Constructors
    anisotropicVorton(const ParameterClass& p, Istream& is)
    : turbulentStructure(p, is),
      p_(p),
      epsilon_(readScalar(is))
    {}
      
    anisotropicVorton
    (
      const ParameterClass& p,
      BoostRandomGen& r, 
      const point& footPoint, 
      const vector& initialDelta,
      label creaface
    )
    : turbulentStructure(p, r, footPoint, initialDelta, creaface),
      p_(p),
      epsilon_(1.0)
    {}
    
    anisotropicVorton(const anisotropicVorton& o)
    : turbulentStructure(o),
      p_(o.p_),
      epsilon_(o.epsilon_)
    {}

    ~anisotropicVorton()
    {}

    tensor parameters() const
    {
        return tensor
               (
                   epsilon_, 0, 0,
                   0, 0, 0,
                   0, 0, 0
               );
    }

    vector fluctuation(const vector& x) const
    {
        vector delta_x = x - location();
        
        vector L=p_.Lmod()[creaFace_];
        vector 
            e1=p_.e1()[creaFace_], 
            e2=p_.e2()[creaFace_], 
            e3=p_.e3()[creaFace_];
        
        double Xx=delta_x&e1;
        double Yy=delta_x&e2;
        double Zz=delta_x&e3;
        
        double rx=p_.rx_[creaFace_];
        double ry=p_.ry_[creaFace_];
        double rz=p_.rz_[creaFace_];
        double sx=p_.sx_[creaFace_];
        double sy=p_.sy_[creaFace_];
        double sz=p_.sz_[creaFace_];

        double sx2=sx*sx;
        double sy2=sy*sy;
        double sz2=sz*sz;
        double xx2=Xx*Xx;
        double yy2=Yy*Yy;
        double zz2=Zz*Zz;

        double e=exp( -0.5 * (xx2/sx2 + yy2/sy2 + zz2/sz2) );

        tensor trs = tensor(e1, e2, e3).T();

        return vector
        (
            transform
            (

                trs,

                /*::sqrt(::pow(M_PI,3./2.)) **/ epsilon_ * e * vector
                (
                    Yy*Zz*( ry/sz2 - rz/sy2 ),
                    Xx*Zz*( rz/sx2 - rx/sz2 ),
                    Xx*Yy*( rx/sy2 - ry/sx2 )
                )

            )
        );
        
    }
    
    // Check

    // Edit
    void randomize(BoostRandomGen& rand)
    {
        turbulentStructure::randomize(rand);
        epsilon_ = 2.0*(rand() - 0.5);
    }
    
    virtual scalar vol() const
    {
        vector L = p_.Lmod()[creaFace_];
        return L.x()*L.y()*L.z();
    }
    
    void operator=(const anisotropicVorton& rhs)
    {
        // Check for assignment to self
        if (this == &rhs)
        {
            FatalErrorIn("anisotropicVorton::operator=(const anisotropicVorton&)")
                << "Attempted assignment to self"
                << abort(FatalError);
        }

        turbulentStructure::operator=(rhs);
        epsilon_=rhs.epsilon_;
    }

    // Required for storage in SLList
    autoPtr<anisotropicVorton> clone() const
    {
        return autoPtr<anisotropicVorton>
        (
            new anisotropicVorton(*this)
        );
    }

    virtual vector randomParameter() const
    {
        return vector::one*epsilon_;
    }

    // Friend Functions

    // Friend Operators
    friend Ostream& operator<<(Ostream& s, const anisotropicVorton& ht)
    {
        s << static_cast<const turbulentStructure&>(ht);
        s << ht.epsilon_ << endl;
        return s;
    }

    friend Istream& operator>>(Istream& s, anisotropicVorton& ht)
    {
        s >> static_cast<turbulentStructure&>(ht);
        ht.epsilon_ = readScalar(s);
        return s;
    }
};




typedef anisotropicVorton<AnisotropicVorton_Analytic_Parameters> anisotropicVorton_Analytic;
typedef anisotropicVorton<AnisotropicVorton_PseudoInv_Parameters> anisotropicVorton_PseudoInv;
typedef anisotropicVorton<AnisotropicVorton_NumOpt_Parameters> anisotropicVorton_NumOpt;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
