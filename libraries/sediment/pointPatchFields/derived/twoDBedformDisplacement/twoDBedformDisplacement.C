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

#include "twoDBedformDisplacement.H"
#include "pointPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"
#include "polyMesh.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

twoDBedformDisplacementPointPatchVectorField::
twoDBedformDisplacementPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(p, iF),
    amplitude_(vector::zero),
    lambda_(0.0)
{}


twoDBedformDisplacementPointPatchVectorField::
twoDBedformDisplacementPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchField<vector>(p, iF, dict),
    amplitude_(dict.lookup("amplitude")),
    lambda_(readScalar(dict.lookup("lambda")))
{
    if (!dict.found("value"))
    {
        updateCoeffs();
    }
}


twoDBedformDisplacementPointPatchVectorField::
twoDBedformDisplacementPointPatchVectorField
(
    const twoDBedformDisplacementPointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchField<vector>(ptf, p, iF, mapper),
    amplitude_(ptf.amplitude_),
    lambda_(ptf.lambda_)
{}


twoDBedformDisplacementPointPatchVectorField::
twoDBedformDisplacementPointPatchVectorField
(
    const twoDBedformDisplacementPointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(ptf, iF),
    amplitude_(ptf.amplitude_),
    lambda_(ptf.lambda_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void twoDBedformDisplacementPointPatchVectorField::updateCoeffs()
{
    Info << "this->updated() = " << this->updated() << endl;
    if (this->updated())
    {
        return;
    }

    const polyMesh& mesh = patch().boundaryMesh().mesh()();
    const Time& t = mesh.time();

    const pointField points(patch().localPoints());

    Info << "twoDBedformDisplacementPointPatchVectorField::updateCoeffs() in " 
        << internalField().name() << endl;

    Field<vector>::operator=
    (
        amplitude_*sin(constant::mathematical::twoPi*t.value()/5.0)*sin(constant::mathematical::twoPi*points.component(0)/lambda_)
    );

    fixedValuePointPatchField<vector>::updateCoeffs();

    Info << "this->updated() = " << this->updated() << endl;
}


void twoDBedformDisplacementPointPatchVectorField::write(Ostream& os) const
{
    pointPatchField<vector>::write(os);
    os.writeKeyword("amplitude")
        << amplitude_ << token::END_STATEMENT << nl;
    os.writeKeyword("lambda")
        << lambda_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePointPatchTypeField
(
    pointPatchVectorField,
    twoDBedformDisplacementPointPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
