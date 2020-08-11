/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2016-2017 Wikki Ltd
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

#include "matchedFlowRateInletVelocityFaPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "areaFields.H"
#include "faPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::matchedFlowRateInletVelocityFaPatchVectorField::
matchedFlowRateInletVelocityFaPatchVectorField
(
    const faPatch& p,
    const DimensionedField<vector, areaMesh>& iF
)
:
    fixedValueFaPatchVectorField(p, iF),
    outletPatchName_()
{}


Foam::matchedFlowRateInletVelocityFaPatchVectorField::
matchedFlowRateInletVelocityFaPatchVectorField
(
    const matchedFlowRateInletVelocityFaPatchVectorField& ptf,
    const faPatch& p,
    const DimensionedField<vector, areaMesh>& iF,
    const faPatchFieldMapper& mapper
)
:
    fixedValueFaPatchVectorField(ptf, p, iF, mapper),
    outletPatchName_(ptf.outletPatchName_)
{}


Foam::matchedFlowRateInletVelocityFaPatchVectorField::
matchedFlowRateInletVelocityFaPatchVectorField
(
    const faPatch& p,
    const DimensionedField<vector, areaMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFaPatchVectorField(p, iF, dict),
    outletPatchName_(dict.lookup("outletPatch"))
{}


Foam::matchedFlowRateInletVelocityFaPatchVectorField::
matchedFlowRateInletVelocityFaPatchVectorField
(
    const matchedFlowRateInletVelocityFaPatchVectorField& ptf
)
:
    fixedValueFaPatchVectorField(ptf),
    outletPatchName_(ptf.outletPatchName_)
{}


Foam::matchedFlowRateInletVelocityFaPatchVectorField::
matchedFlowRateInletVelocityFaPatchVectorField
(
    const matchedFlowRateInletVelocityFaPatchVectorField& ptf,
    const DimensionedField<vector, areaMesh>& iF
)
:
    fixedValueFaPatchVectorField(ptf, iF),
    outletPatchName_(ptf.outletPatchName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::matchedFlowRateInletVelocityFaPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Find corresponding outlet patch
    const label outletPatchID =
        patch().boundaryMesh().findPatchID(outletPatchName_);

    if (outletPatchID < 0)
    {
        FatalErrorInFunction
            << "Unable to find outlet patch " << outletPatchName_
            << exit(FatalError);
    }

    const faPatch& p = patch();
    const faPatch& outletPatch = patch().boundaryMesh()[outletPatchID];

    // Lookup non-const access to velocity field
    areaVectorField& U
    (
        const_cast<areaVectorField&>
        (
            dynamic_cast<const areaVectorField&>(internalField())
        )
    );

    // Get the corresponding outlet velocity patch field
    faPatchVectorField& outletPatchU = U.boundaryFieldRef()[outletPatchID];

    // Ensure that the corresponding outlet velocity patch field is up-to-date
    outletPatchU.updateCoeffs();

    // Calculated the total lenght of the outlet patch
    const scalar outletLength = gSum(outletPatch.magEdgeLengths());
    //Pout << "outlet patch on faMesh has a length of " << outletLength << endl;

    // Calculate the outlet patch flow rate
    const scalar flowRate = gSum(outletPatch.edgeLengths() & outletPatchU);

    // Calculate the mean velocity on the inlet patch
    const scalar inletMeanU = flowRate/(outletLength+SMALL);

    operator==(-inletMeanU*p.edgeNormals());
}


void Foam::matchedFlowRateInletVelocityFaPatchVectorField::write(Ostream& os) const
{
    os.writeEntry("outletPatch", outletPatchName_);
    fixedValueFaPatchVectorField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

makeFaPatchTypeField
(
    faPatchVectorField,
    matchedFlowRateInletVelocityFaPatchVectorField
);

} // End namespace Foam


// ************************************************************************* //
