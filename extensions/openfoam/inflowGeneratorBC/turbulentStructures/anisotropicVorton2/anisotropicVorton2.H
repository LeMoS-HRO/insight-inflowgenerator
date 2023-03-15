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

#ifndef anisotropicVorton2_H
#define anisotropicVorton2_H

#include "turbulentStructure.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class anisotropicVorton2;

Ostream& operator<<(Ostream&, const anisotropicVorton2&);
Istream& operator>>(Istream&, anisotropicVorton2&);

/*---------------------------------------------------------------------------*\
                           Class anisotropicVorton2 Declaration
\*---------------------------------------------------------------------------*/

class anisotropicVorton2
: public turbulentStructure
{

public:

    class StructureParameters
    {
    public:
        StructureParameters();
        StructureParameters(const dictionary&);

	void autoMap
	(
	    const fvPatchFieldMapper&
	);

	//- Reverse map the given fvPatchField onto this fvPatchField
	void rmap
	(
	    const fvPatchField<vector>&,
	    const labelList&
	);
    
        void write(Ostream&) const;
	
	scalar calcInfluenceLength(const vector& L);
    };

protected:
    scalar epsilon_;
//   vector epsilon_;
    
    scalar rx_, ry_, rz_;
    scalar sx_, sy_, sz_;
    scalar k0_, C1_;

public:

    // Constructors
    anisotropicVorton2();
    anisotropicVorton2(Istream&);
    anisotropicVorton2
    (
      BoostRandomGen&, 
      const point& p, 
      const vector& initialDelta, 
      const vector& v, 
      const vector& L, 
      scalar minL,
      label creaface,
      const symmTensor& R,
      vector yw,
      vector Rpmult=vector::one
    );
    
    anisotropicVorton2(const anisotropicVorton2& p);
    
    ~anisotropicVorton2();
    
    tensor parameters() const
    {
      return tensor
      (
	rx_, ry_, rz_,
	sx_, sy_, sz_,
	k0_, C1_, epsilon_
      );
    }

    vector fluctuation(const StructureParameters& pa, const vector&) const;
    
    // Check

    // Edit
    void randomize(BoostRandomGen&);

    void operator=(const anisotropicVorton2&);
//     bool operator!=(const anisotropicVorton2&) const;

    // Required for storage in SLList
    static autoPtr<anisotropicVorton2> New(Istream& s);
    autoPtr<anisotropicVorton2> clone() const;
    
    virtual vector randomParameter() const
    {
      return vector::one*epsilon_;
    }

    // Friend Functions

    // Friend Operators
    friend Ostream& operator<<(Ostream&, const anisotropicVorton2&);
    friend Istream& operator>>(Istream&, anisotropicVorton2&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
