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

#include "vanRijn1984PickupFunctionCb_starFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "wallFvPatch.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //
tmp<scalarField> vanRijn1984PickupFunctionCb_starFvPatchScalarField::calcCb_star() const
{
    //find the Shields parameter for this patch, which should be updated 
    //before calling into this member function.
    const fvPatchField<scalar>& Shields_wall =
        patch().lookupPatchField<volScalarField, scalar>("wallShieldsNumber");

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            this->db().time().constant(),
            this->db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    dimensionedScalar nu(transportProperties.lookup("nu"));
    dimensionedScalar rhow(transportProperties.lookup("rhow"));
    dimensionedScalar rhos(transportProperties.lookup("rhos"));
    dimensionedScalar R=(rhos-rhow)/rhow;
    dimensionedScalar diam(transportProperties.lookup("diam"));

    word sedBC(transportProperties.lookup("sedBC"));

    if(!sedBC.match("vanRijn1984PickupFunction"))
    {
        FatalErrorIn("vanRijn1984PickupFunctionCw_star")
           << "The option for SedBC in transportDict: " << sedBC
           << " is not vanRijn1984PickupFuction." << nl
           << "You cannot call this BC."
           << abort(FatalError);
    }
/*
    IOdictionary gravitationalProperties
    (
       IOobject
       (
          "gravitationalProperties",
          this->db().time().constant(),
          this->db(),
          IOobject::MUST_READ_IF_MODIFIED,
          IOobject::NO_WRITE
       )
    );

    const dimensionedVector g(gravitationalProperties.lookup("g"));
*/
    uniformDimensionedVectorField g
    (
        IOobject
        (
            "g",
            this->db().time().constant(),
            this->db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    
    const fvPatchField<vector>& Vs_wall =
        patch().lookupPatchField<volVectorField, vector>("Vs");

    //particle parameter Dstar
    tmp<scalarField> tD_star(new scalarField(*this));
    scalarField& D_star = tD_star();

    forAll(D_star, faceI)
    {
       D_star[faceI] = diam.value()*pow((R*mag(g)/pow(nu.value(),2.0)).value(), 1.0/3.0);
    }

    //critical Shields number
    scalar Shields_critical = 0.0;

    tmp<scalarField> tCb_starw(new scalarField(*this));
    scalarField& Cb_starw = tCb_starw();

    forAll(Cb_starw, faceI)
    {
        scalar Vs_temp = mag(Vs_wall[faceI]);
        if(D_star[faceI]<10.0)
        {
           Shields_critical = (16.0*pow(Vs_temp,2.0)
                   /(D_star[faceI]*R*mag(g)*diam.value())).value();
        }
        else
        {
           Shields_critical = (0.16*pow(Vs_temp,2.0)
                    /(R*mag(g)*diam.value())).value();
        }

        if (Shields_wall[faceI] > Shields_critical)
        {
            scalar temp1 = 
                pow(Shields_wall[faceI]/Shields_critical-1.0,1.5);
            scalar temp2 = 
                pow(R.value(),0.6)
               *pow(mag(g).value(),0.6)
               *pow(diam.value(),0.8)
               /pow(nu.value(),0.2);
            Cb_starw[faceI] = 0.00033*temp1*temp2/mag(Vs_wall[faceI]);
        }
        else
        {
            Cb_starw[faceI] = 0.0;
        }

        if (debug)
        {
            Info<< "Shields = " << Shields_wall[faceI]
                << ", Shields_critical = " << Shields_critical
                << ", D_star = " << D_star[faceI]
                << ", Vs_wall = " << Vs_temp
                << ", Cb_starw = " << Cb_starw[faceI]
                << endl;
        }
    }

    return tCb_starw;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

vanRijn1984PickupFunctionCb_starFvPatchScalarField::vanRijn1984PickupFunctionCb_starFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    equilibriumCb_starFvPatchScalarField(p, iF)
{}


vanRijn1984PickupFunctionCb_starFvPatchScalarField::vanRijn1984PickupFunctionCb_starFvPatchScalarField
(
    const vanRijn1984PickupFunctionCb_starFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    equilibriumCb_starFvPatchScalarField(ptf, p, iF, mapper)
{}


vanRijn1984PickupFunctionCb_starFvPatchScalarField::vanRijn1984PickupFunctionCb_starFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    equilibriumCb_starFvPatchScalarField(p, iF, dict)
{}


vanRijn1984PickupFunctionCb_starFvPatchScalarField::vanRijn1984PickupFunctionCb_starFvPatchScalarField
(
    const vanRijn1984PickupFunctionCb_starFvPatchScalarField& rwfpsf
)
:
    equilibriumCb_starFvPatchScalarField(rwfpsf)
{}


vanRijn1984PickupFunctionCb_starFvPatchScalarField::vanRijn1984PickupFunctionCb_starFvPatchScalarField
(
    const vanRijn1984PickupFunctionCb_starFvPatchScalarField& rwfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    equilibriumCb_starFvPatchScalarField(rwfpsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void vanRijn1984PickupFunctionCb_starFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    equilibriumCb_starFvPatchScalarField::autoMap(m);
}


void vanRijn1984PickupFunctionCb_starFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    equilibriumCb_starFvPatchScalarField::rmap(ptf, addr);

    const vanRijn1984PickupFunctionCb_starFvPatchScalarField& nrwfpsf =
        refCast<const vanRijn1984PickupFunctionCb_starFvPatchScalarField>(ptf);
}


void vanRijn1984PickupFunctionCb_starFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    writeLocalEntries(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    vanRijn1984PickupFunctionCb_starFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
