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

#include "ZyserFredsoeCb_starFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "wallFvPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //
tmp<scalarField> ZyserFredsoeCb_starFvPatchScalarField::calcCb_star() const
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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

ZyserFredsoeCb_starFvPatchScalarField::ZyserFredsoeCb_starFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    equilibriumCb_starFvPatchScalarField(p, iF)
{}


ZyserFredsoeCb_starFvPatchScalarField::ZyserFredsoeCb_starFvPatchScalarField
(
    const ZyserFredsoeCb_starFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    equilibriumCb_starFvPatchScalarField(ptf, p, iF, mapper)
{}


ZyserFredsoeCb_starFvPatchScalarField::ZyserFredsoeCb_starFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    equilibriumCb_starFvPatchScalarField(p, iF, dict)
{}


ZyserFredsoeCb_starFvPatchScalarField::ZyserFredsoeCb_starFvPatchScalarField
(
    const ZyserFredsoeCb_starFvPatchScalarField& rwfpsf
)
:
    equilibriumCb_starFvPatchScalarField(rwfpsf)
{}


ZyserFredsoeCb_starFvPatchScalarField::ZyserFredsoeCb_starFvPatchScalarField
(
    const ZyserFredsoeCb_starFvPatchScalarField& rwfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    equilibriumCb_starFvPatchScalarField(rwfpsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void ZyserFredsoeCb_starFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    equilibriumCb_starFvPatchScalarField::autoMap(m);
}


void ZyserFredsoeCb_starFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    equilibriumCb_starFvPatchScalarField::rmap(ptf, addr);

    const ZyserFredsoeCb_starFvPatchScalarField& nrwfpsf =
        refCast<const ZyserFredsoeCb_starFvPatchScalarField>(ptf);
}


void ZyserFredsoeCb_starFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    writeLocalEntries(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    ZyserFredsoeCb_starFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
