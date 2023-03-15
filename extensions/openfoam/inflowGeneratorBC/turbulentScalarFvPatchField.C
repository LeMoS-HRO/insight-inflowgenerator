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

#include "turbulentScalarFvPatchField.H"

#include "uniof.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


turbulentScalarFvPatchField::turbulentScalarFvPatchField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF),
    UName_("U")
{}



turbulentScalarFvPatchField::turbulentScalarFvPatchField
(
    const turbulentScalarFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper),
    UName_(ptf.UName_),
    mean_(ptf.mean_().clone()),
    variance_(ptf.variance_().clone()),
    turbulentFlux_(ptf.turbulentFlux_().clone())
{}



turbulentScalarFvPatchField::turbulentScalarFvPatchField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF),
    UName_(dict.lookupOrDefault<word>("UName", "U")),
    mean_(FieldDataProvider<scalar>::New(dict.lookup("mean"))),
    variance_(FieldDataProvider<scalar>::New(dict.lookup("variance"))),
    turbulentFlux_(FieldDataProvider<vector>::New(dict.lookup("turbulentFlux")))
{
  if (dict.found("value"))
  {
      fvPatchField<scalar>::operator==(Field<scalar>("value", dict, p.size()));
  }
  else
  {
      // Note: we use evaluate() here to trigger updateCoeffs followed
      //       by re-setting of fvatchfield::updated_ flag. This is
      //       so if first use is in the next time step it retriggers
      //       a new update.
      this->evaluate(
#if OF_VERSION>=060000 //defined(OFdev)||defined(OFesi1806)
            Pstream::commsTypes::blocking
#else
            Pstream::blocking
#endif
            );
  }
}



turbulentScalarFvPatchField::turbulentScalarFvPatchField
(
    const turbulentScalarFvPatchField& ptf
)
:
    fixedValueFvPatchField<scalar>(ptf),
    UName_(ptf.UName_),
    mean_(ptf.mean_().clone()),
    variance_(ptf.variance_().clone()),
    turbulentFlux_(ptf.turbulentFlux_().clone())
{}



turbulentScalarFvPatchField::turbulentScalarFvPatchField
(
    const turbulentScalarFvPatchField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(ptf, iF),
    UName_(ptf.UName_),
    mean_(ptf.mean_().clone()),
    variance_(ptf.variance_().clone()),
    turbulentFlux_(ptf.turbulentFlux_().clone())
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void turbulentScalarFvPatchField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchField<scalar>::autoMap(m);
}


void turbulentScalarFvPatchField::rmap
(
    const fvPatchField<scalar>& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchField<scalar>::rmap(ptf, addr);
/*
    const turbulentScalarFvPatchField& tiptf =
        refCast<const turbulentScalarFvPatchField >(ptf);
*/
}



void turbulentScalarFvPatchField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    
#define ADD_FIELD(NAME, TYPE, TYPE2) \
    {\
        NAME##_.reset\
        (\
            new vol##TYPE##Field\
            (\
                IOobject\
                (\
                    #NAME ,\
                    U.mesh().time().timeName(),\
                    U.mesh(),\
                    IOobject::NO_READ,\
                    IOobject::NO_WRITE\
                ),\
                U.mesh(),\
                dimensioned##TYPE("", dimless, pTraits<TYPE2>::zero),\
                "calculated"\
            )\
        );\
    }
        
    
    const inflowGeneratorBaseFvPatchVectorField* Uinfl = NULL;
    
    if (this->db().foundObject<volVectorField>(UName_))
    {
        const volVectorField& U = this->db().lookupObject<volVectorField>(UName_);
        Uinfl = dynamic_cast< const inflowGeneratorBaseFvPatchVectorField* >( &U.boundaryField()[patch().index()] );
    }
    
    scalarField mean (mean_()(this->db().time().timeOutputValue(), patch().Cf()));

    if (!Uinfl)
    {
        WarningIn("turbulentScalarFvPatchField::updateCoeffs")
        << " cannot access underlying velocity inflow generator! Disabling scalar fluctuations for this time step!" <<endl;
        
        fixedValueFvPatchField<scalar>::operator==( mean );
    }
    else
    {   
        const volVectorField& U = this->db().lookupObject<volVectorField>(UName_);
        
        scalarField variance (variance_()(this->db().time().timeOutputValue(), patch().Cf()));
        vectorField turbulentFlux (turbulentFlux_()(this->db().time().timeOutputValue(), patch().Cf()));
        const symmTensorField& R = Uinfl->R();

        vectorField w123 ( inv(R) & turbulentFlux );
        scalarField w4radic 
        (
            variance - /*( (R&w123) & vector::one )*/ 
            (
            R.component(symmTensor::XX)*sqr(w123.component(vector::X))
            +R.component(symmTensor::YY)*sqr(w123.component(vector::Y))
            +R.component(symmTensor::ZZ)*sqr(w123.component(vector::Z))
            +2.*R.component(symmTensor::XY)*w123.component(vector::X)*w123.component(vector::Y)
            +2.*R.component(symmTensor::XZ)*w123.component(vector::X)*w123.component(vector::Z)
            +2.*R.component(symmTensor::YZ)*w123.component(vector::Y)*w123.component(vector::Z)
            )
        );
        scalarField w4 ( sqrt(max(scalar(0.0), w4radic)) );
        
        scalarField raw_fluctuations ( Uinfl->normalizedScalarFluctuations() );
        
        scalarField scalar_fluctuations ( mean + (w123 & Uinfl->fluctuations()) - w4 * raw_fluctuations );
        
        if (!w123_.valid()) ADD_FIELD(w123, Vector, vector);
        if (!w4radic_.valid()) ADD_FIELD(w4radic, Scalar, scalar);
        if (!w4_.valid()) ADD_FIELD(w4, Scalar, scalar);
        if (!flucraw_.valid()) ADD_FIELD(flucraw, Scalar, scalar);
        if (!flucunlim_.valid()) ADD_FIELD(flucunlim, Scalar, scalar);
        if (!flucmean_.valid()) ADD_FIELD(flucmean, Scalar, scalar);
        if (!flucc1_.valid()) ADD_FIELD(flucc1, Scalar, scalar);
        if (!flucc2_.valid()) ADD_FIELD(flucc2, Scalar, scalar);
        
        UNIOF_BOUNDARY_NONCONST(w123_())[patch().index()]=w123;
        UNIOF_BOUNDARY_NONCONST(w4radic_())[patch().index()]=w4radic;
        UNIOF_BOUNDARY_NONCONST(w4_())[patch().index()]=w4;
        UNIOF_BOUNDARY_NONCONST(flucraw_())[patch().index()]=raw_fluctuations;
        UNIOF_BOUNDARY_NONCONST(flucunlim_())[patch().index()]=scalar_fluctuations;
        UNIOF_BOUNDARY_NONCONST(flucmean_())[patch().index()]=mean;
        UNIOF_BOUNDARY_NONCONST(flucc1_())[patch().index()]=(w123 & Uinfl->fluctuations());
        UNIOF_BOUNDARY_NONCONST(flucc2_())[patch().index()]=w4 * raw_fluctuations;
        
        if (U.mesh().time().outputTime())
        {
            w123_->write();
            w4radic_->write();
            w4_->write();
            flucraw_->write();
            flucunlim_->write();
            flucmean_->write();
            flucc1_->write();
            flucc2_->write();
        }
        
        fixedValueFvPatchField<scalar>::operator==( min(1.0, max(0.0, scalar_fluctuations)) );
    }
    
    fixedValueFvPatchField<scalar>::updateCoeffs();
}




void turbulentScalarFvPatchField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    mean_().writeEntry("mean", os);
    variance_().writeEntry("variance", os);
    turbulentFlux_().writeEntry("turbulentFlux", os);
    this->writeEntry("value", os);
}

makePatchTypeField
(
    fvPatchScalarField,
    turbulentScalarFvPatchField
);
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

