/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2010 OpenCFD Ltd.
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

#include "equilibriumCb_starFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "wallFvPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void equilibriumCb_starFvPatchScalarField::checkType()
{
    if (!isA<wallFvPatch>(patch()))
    {
        FatalErrorIn("equilibriumCb_starFvPatchScalarField::checkType()")
            << "Invalid equilibrium Cb_star specification" << nl
            << "    Patch type for patch " << patch().name()
            << " must be wall" << nl
            << "    Current patch type is " << patch().type() << nl << endl
            << abort(FatalError);
    }
}

tmp<scalarField> equilibriumCb_starFvPatchScalarField::calcCb_star() const
{
    //find the Shields parameter for this patch, which should be updated 
    //before calling into this member function.
    const fvPatchField<scalar>& Shields_wall =
        patch().lookupPatchField<volScalarField, scalar>("wallShieldsNumber");
    
    
    //critical Shields number: 0.045 for Zyserman and Fredsoe (1994) formula
    const scalar Shields_critical = 0.045;

    tmp<scalarField> tCb_starw(new scalarField(*this));
    scalarField& Cb_starw = tCb_starw();

    forAll(Cb_starw, faceI)
    {
        if (Shields_wall[faceI] > Shields_critical)
        { 
            scalar temp = pow(Shields_wall[faceI]-Shields_critical,1.75);
            Cb_starw[faceI] = 0.331*temp/(1+0.331*temp/0.46);
        }
        else
        {
            Cb_starw[faceI] = 0.0;
        }

        if (debug)
        {
            Info<< "Shields = " << Shields_wall[faceI]
                << ", Cb_starw = " << Cb_starw[faceI]
                << endl;
        }
    }

    return tCb_starw;
}


void equilibriumCb_starFvPatchScalarField::writeLocalEntries(Ostream& os) const
{
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

equilibriumCb_starFvPatchScalarField::equilibriumCb_starFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF)
{
    checkType();
}


equilibriumCb_starFvPatchScalarField::equilibriumCb_starFvPatchScalarField
(
    const equilibriumCb_starFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper)
{
    checkType();
}


equilibriumCb_starFvPatchScalarField::equilibriumCb_starFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict)
{
    checkType();
}


equilibriumCb_starFvPatchScalarField::equilibriumCb_starFvPatchScalarField
(
    const equilibriumCb_starFvPatchScalarField& wfpsf
)
:
    fixedValueFvPatchScalarField(wfpsf)
{
    checkType();
}


equilibriumCb_starFvPatchScalarField::equilibriumCb_starFvPatchScalarField
(
    const equilibriumCb_starFvPatchScalarField& wfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(wfpsf, iF)
{
    checkType();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void equilibriumCb_starFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    operator==(calcCb_star());

    fixedValueFvPatchScalarField::updateCoeffs();
}

void equilibriumCb_starFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    writeLocalEntries(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    equilibriumCb_starFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
