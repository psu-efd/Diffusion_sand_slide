/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
  \\    /   O peration     |
  \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
  \\/     M anipulation  | Copyright (C) 2015-2018 OpenCFD Ltd.
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

#include "timeVaryingMappedPerturbationFixedValueFvPatchField.H"
#include "Time.H"
#include "AverageField.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::timeVaryingMappedPerturbationFixedValueFvPatchField<Type>::
timeVaryingMappedPerturbationFixedValueFvPatchField
(
 const fvPatch& p,
 const DimensionedField<Type, volMesh>& iF
 )
:
    fixedValueFvPatchField<Type>(p, iF),
    uniformValue_(),
    fluctuationScale_(0.0),
    T_(0.0),
    curTimeIndex_(-1)
{}


template<class Type>
Foam::timeVaryingMappedPerturbationFixedValueFvPatchField<Type>::
timeVaryingMappedPerturbationFixedValueFvPatchField
(
 const fvPatch& p,
 const DimensionedField<Type, volMesh>& iF,
 const dictionary& dict
 )
:
    fixedValueFvPatchField<Type>(p, iF, dict, false),
    uniformValue_
    (
     new PatchFunction1Types::MappedFile<Type>
     (
      p.patch(),
      "uniformValue",
      dict,
      iF.name(),          // field table name
      true                // face values
     )
    ),
    fluctuationScale_(dict.lookupOrDefault<scalar>("fluctuationScale", 0.0)),
    T_(dict.lookupOrDefault<scalar>("T", VGREAT)),  //default is VGREAT->not changing
    curTimeIndex_(-1)
{
    if (dict.found("value"))
    {
        fvPatchField<Type>::operator==(Field<Type>("value", dict, p.size()));
    }
    else
    {
        // Note: we use evaluate() here to trigger updateCoeffs followed
        //       by re-setting of fvatchfield::updated_ flag. This is
        //       so if first use is in the next time step it retriggers
        //       a new update.
        this->evaluate(Pstream::commsTypes::blocking);
    }
}


template<class Type>
Foam::timeVaryingMappedPerturbationFixedValueFvPatchField<Type>::
timeVaryingMappedPerturbationFixedValueFvPatchField
(
 const timeVaryingMappedPerturbationFixedValueFvPatchField<Type>& ptf,
 const fvPatch& p,
 const DimensionedField<Type, volMesh>& iF,
 const fvPatchFieldMapper& mapper
 )
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    uniformValue_
    (
     new PatchFunction1Types::MappedFile<Type>
     (
      ptf.uniformValue_,
      p.patch()
     )
    ),
    fluctuationScale_(ptf.fluctuationScale_),
    T_(ptf.T_),
    curTimeIndex_(-1)
{}


template<class Type>
Foam::timeVaryingMappedPerturbationFixedValueFvPatchField<Type>::
timeVaryingMappedPerturbationFixedValueFvPatchField
(
 const timeVaryingMappedPerturbationFixedValueFvPatchField<Type>& ptf
 )
:
    fixedValueFvPatchField<Type>(ptf),
    uniformValue_
    (
     new PatchFunction1Types::MappedFile<Type>
     (
      ptf.uniformValue_,
      this->patch().patch()
     )
    ),
    fluctuationScale_(ptf.fluctuationScale_),
    T_(ptf.T_),
    curTimeIndex_(-1)
{}


template<class Type>
Foam::timeVaryingMappedPerturbationFixedValueFvPatchField<Type>::
timeVaryingMappedPerturbationFixedValueFvPatchField
(
 const timeVaryingMappedPerturbationFixedValueFvPatchField<Type>& ptf,
 const DimensionedField<Type, volMesh>& iF
 )
:
    fixedValueFvPatchField<Type>(ptf, iF),
    uniformValue_
    (
     new PatchFunction1Types::MappedFile<Type>
     (
      ptf.uniformValue_,
      this->patch().patch()
     )
    ),
    fluctuationScale_(ptf.fluctuationScale_),
    T_(ptf.T_),
    curTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::timeVaryingMappedPerturbationFixedValueFvPatchField<Type>::autoMap
(
 const fvPatchFieldMapper& m
 )
{
    fixedValueFvPatchField<Type>::autoMap(m);
    uniformValue_().autoMap(m);
}


template<class Type>
void Foam::timeVaryingMappedPerturbationFixedValueFvPatchField<Type>::rmap
(
 const fvPatchField<Type>& ptf,
 const labelList& addr
 )
{
    fixedValueFvPatchField<Type>::rmap(ptf, addr);

    const timeVaryingMappedPerturbationFixedValueFvPatchField<Type>& tiptf =
        refCast<const timeVaryingMappedPerturbationFixedValueFvPatchField<Type>>(ptf);

    uniformValue_().rmap(tiptf.uniformValue_(), addr);
}


    template<class Type>
void Foam::timeVaryingMappedPerturbationFixedValueFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        const scalar t = this->db().time().timeOutputValue();
        scalar perturb = fluctuationScale_*(Foam::sin(2.0*3.1415926*t/(T_+VSMALL)));

        fvPatchField<Type>::operator==(uniformValue_->value(t)*(1.0+perturb));

        if (debug)
        {
            Pout<< "updateCoeffs : set fixedValue to min:" << gMin(*this)
                << " max:" << gMax(*this)
                << " avg:" << gAverage(*this)
                << " perturb: " << perturb << endl;
        }

    }

    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::timeVaryingMappedPerturbationFixedValueFvPatchField<Type>::write
(
 Ostream& os
 ) const
{
    //debug
    //Info << this->internalField().name() << endl;
    //Info << (*uniformValue_) << endl;

    fvPatchField<Type>::write(os);
    uniformValue_->writeData(os);
    os.writeEntry("fluctuationScale", fluctuationScale_);
    os.writeEntry("T", T_);

    this->writeEntry("value", os);
}


// ************************************************************************* //
