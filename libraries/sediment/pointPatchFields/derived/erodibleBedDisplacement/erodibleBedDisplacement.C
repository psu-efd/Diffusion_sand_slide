/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "erodibleBedDisplacement.H"
#include "pointPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"
#include "polyMesh.H"
#include "mathematicalConstants.H"

#include "volFields.H"
#include "primitivePatchInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

erodibleBedDisplacementPointPatchVectorField::
erodibleBedDisplacementPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(p, iF),
    frozenArea_(false),
    x0_(0.0)
{}


erodibleBedDisplacementPointPatchVectorField::
erodibleBedDisplacementPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchField<vector>(p, iF, dict),
    frozenArea_(dict.lookupOrDefault("frozenArea", false)),
    x0_(dict.lookupOrDefault("x0",0.0))
{
    Info << "In erodibleBedDisplacement BC: ";
    Info << "frozenArea = " << frozenArea_;
    Info << ", x0 = " << x0_ << endl;

    if (!dict.found("value"))
    {
        updateCoeffs();
    }
}


erodibleBedDisplacementPointPatchVectorField::
erodibleBedDisplacementPointPatchVectorField
(
    const erodibleBedDisplacementPointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchField<vector>(ptf, p, iF, mapper),
    frozenArea_(ptf.frozenArea_),
    x0_(ptf.x0_)
{}


erodibleBedDisplacementPointPatchVectorField::
erodibleBedDisplacementPointPatchVectorField
(
    const erodibleBedDisplacementPointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(ptf, iF),
    frozenArea_(ptf.frozenArea_),
    x0_(ptf.x0_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void erodibleBedDisplacementPointPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const polyMesh& mesh = patch().boundaryMesh().mesh()();

    const pointField points(patch().localPoints());

    const pointPatch& pp = patch();
    const word& pname = pp.name();

    //vertical direction is fixed at y direction
    vector yDir(0,1,0);

    //find the fluidFaceCenterTotalDeltaY field
    volScalarField& fluidFaceCenterTotalDeltaY = const_cast<volScalarField&>
       (mesh.objectRegistry::lookupObject<volScalarField>("fluidFaceCenterTotalDeltaY"));

    label movingPatchID = mesh.boundaryMesh().findPatchID(pname);

    primitivePatchInterpolation patchInterpolate(mesh.boundaryMesh()[movingPatchID]);

    scalarField tempPointDisplacementy =
             patchInterpolate.faceToPointInterpolate
             (
                 fluidFaceCenterTotalDeltaY.boundaryField()[movingPatchID]
             );

    //freeze the points in frozenArea
    if(frozenArea_)
    {
        tempPointDisplacementy *= pos0(points.component(0) - x0_);
    }

    vectorField tempPointDisplacement = yDir*tempPointDisplacementy;

    Field<vector>::operator=
    (
       tempPointDisplacement
    );

    fixedValuePointPatchField<vector>::updateCoeffs();
}


void erodibleBedDisplacementPointPatchVectorField::write(Ostream& os) const
{
    pointPatchField<vector>::write(os);
    os.writeKeyword("frozenArea")
        << frozenArea_ << token::END_STATEMENT << nl;
    os.writeKeyword("x0")
        << x0_ << token::END_STATEMENT << nl;

    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePointPatchTypeField
(
    pointPatchVectorField,
    erodibleBedDisplacementPointPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
