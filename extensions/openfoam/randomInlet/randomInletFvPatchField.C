/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "randomInletFvPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::randomInletFvPatchField<Type>::randomInletFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    ranGen_(label(0)),
    variance_(),
    mean_(),
//     alpha_(0.1),
    curTimeIndex_(-1)
{}


template<class Type>
Foam::randomInletFvPatchField<Type>::randomInletFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF/*, dict, false*/),
    ranGen_(label(0)),
    variance_(FieldDataProvider<Type>::New(dict.lookup("variance"))),
    mean_(FieldDataProvider<Type>::New(dict.lookup("mean"))),
//     alpha_(dict.lookupOrDefault<scalar>("alpha", 0.1)),
    curTimeIndex_(-1)
{
    if (dict.found("value"))
    {
        fixedValueFvPatchField<Type>::operator==
        (
            Field<Type>("value", dict, p.size())
        );
    }
    else
    {
        fixedValueFvPatchField<Type>::operator==( mean_()(this->db().time().timeOutputValue(), this->patch().Cf()) );
    }
}


template<class Type>
Foam::randomInletFvPatchField<Type>::randomInletFvPatchField
(
    const randomInletFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    ranGen_(label(0)),
    variance_(ptf.variance_().clone()),
    mean_(ptf.mean_().clone()),
//     alpha_(ptf.alpha_),
    curTimeIndex_(-1)
{}


template<class Type>
Foam::randomInletFvPatchField<Type>::randomInletFvPatchField
(
    const randomInletFvPatchField<Type>& ptf
)
:
    fixedValueFvPatchField<Type>(ptf),
    ranGen_(ptf.ranGen_),
    variance_(ptf.variance_().clone()),
    mean_(ptf.mean_().clone()),
//     alpha_(ptf.alpha_),
    curTimeIndex_(-1)
{}


template<class Type>
Foam::randomInletFvPatchField<Type>::randomInletFvPatchField
(
    const randomInletFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    ranGen_(ptf.ranGen_),
    variance_(ptf.variance_().clone()),
    mean_(ptf.mean_().clone()),
//     alpha_(ptf.alpha_),
    curTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::randomInletFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchField<Type>::autoMap(m);
//     referenceField_.autoMap(m);
}


template<class Type>
void Foam::randomInletFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchField<Type>::rmap(ptf, addr);

//     const randomInletFvPatchField<Type>& tiptf =
//         refCast<const randomInletFvPatchField<Type>>(ptf);
// 
//     referenceField_.rmap(tiptf.referenceField_, addr);
}


template<class Type>
void Foam::randomInletFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        Field<Type>& patchField = *this;
        
        Field<Type> mean = mean_()(this->db().time().timeOutputValue(), this->patch().Cf());
        Field<Type> variance = variance_()(this->db().time().timeOutputValue(), this->patch().Cf());

        Field<Type> randomField(this->size());
        forAll(mean, facei)
        {
            ranGen_.randomise(randomField[facei]);
        }

        // Correction-factor to compensate for the loss of RMS fluctuation
        // due to the temporal correlation introduced by the alpha parameter.
//         scalar varianceCorr = sqrt(12*(2*alpha_ - sqr(alpha_)))/alpha_;

        Field<Type> sf(this->size(), pTraits<Type>::zero);
        forAll(sf, j)
        {
            for (label k = 0; k<pTraits<Type>::nComponents; k++)
            {
                setComponent(sf[j],k) = sqrt( component(variance[j],k)*3.)*2.;
            }
        }

        patchField =
//             (1 - alpha_)*patchField
//           + alpha_*
            (
                mean
              + cmptMultiply
                (
                    randomField - 0.5*pTraits<Type>::one,
                    sf
                )
            );

        curTimeIndex_ = this->db().time().timeIndex();
    }

    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::randomInletFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    variance_().writeEntry("variance", os);
    mean_().writeEntry("mean", os);
//     os.writeKeyword("alpha") << alpha_ << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}


// ************************************************************************* //
