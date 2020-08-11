/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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

#include "fixedGradientSuspSedFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

fixedGradientSuspSedFvPatchScalarField::fixedGradientSuspSedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchField<scalar>(p, iF)
{}


fixedGradientSuspSedFvPatchScalarField::fixedGradientSuspSedFvPatchScalarField
(
    const fixedGradientSuspSedFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchField<scalar>(ptf, p, iF, mapper)
{}


fixedGradientSuspSedFvPatchScalarField::fixedGradientSuspSedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchField<scalar>(p, iF, dict)
{
    evaluate();
}


fixedGradientSuspSedFvPatchScalarField::fixedGradientSuspSedFvPatchScalarField
(
    const fixedGradientSuspSedFvPatchScalarField& ptf
)
:
    fixedGradientFvPatchField<scalar>(ptf)
{}


fixedGradientSuspSedFvPatchScalarField::fixedGradientSuspSedFvPatchScalarField
(
    const fixedGradientSuspSedFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchField<scalar>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void fixedGradientSuspSedFvPatchScalarField::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    Field<scalar>::operator=
    (
//        this->patchInternalField() + gradient_/this->patch().deltaCoeffs()
        this->patchInternalField() // added by Zheyu Zhou on 04/03/2014
    );

    fvPatchField<scalar>::evaluate();
}

tmp<Field<scalar> > fixedGradientSuspSedFvPatchScalarField::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
//    return gradient()/this->patch().deltaCoeffs();
    return gradient()/this->patch().deltaCoeffs() - gradient()/this->patch().deltaCoeffs(); // modified by Zheyu Zhou on 04/04/2014
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        fixedGradientSuspSedFvPatchScalarField
    );
}


// ************************************************************************* //
