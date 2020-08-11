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

Class
    Foam::matchedFlowRateInletVelocityFaPatchVectorField

Description
    Edge normal fixed value vector field finite area boundary condition

    Describes a surface normal vector boundary condition by its magnitude.
    Note: The value is positive for outward-pointing vectors

Author
    Zeljko Tukovic, FMENA
    Hrvoje Jasak, Wikki Ltd.

SourceFiles
    matchedFlowRateInletVelocityFaPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef matchedFlowRateInletVelocityFaPatchVectorField_H
#define matchedFlowRateInletVelocityFaPatchVectorField_H

#include "faPatchFields.H"
#include "fixedValueFaPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class matchedFlowRateInletVelocityFaPatch Declaration
\*---------------------------------------------------------------------------*/

class matchedFlowRateInletVelocityFaPatchVectorField
:
    public fixedValueFaPatchVectorField
{
    // Private data

        //- Outlet patch name from which the corresponding flow rate is obtained
        word outletPatchName_;


public:

    //- Runtime type information
    TypeName("matchedFlowRateInletVelocity");


    // Constructors

        //- Construct from patch and internal field
        matchedFlowRateInletVelocityFaPatchVectorField
        (
            const faPatch&,
            const DimensionedField<vector, areaMesh>&
        );

        //- Construct from patch, internal field and dictionary
        matchedFlowRateInletVelocityFaPatchVectorField
        (
            const faPatch&,
            const DimensionedField<vector, areaMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  matchedFlowRateInletVelocityFaPatchVectorField
        //  onto a new patch
        matchedFlowRateInletVelocityFaPatchVectorField
        (
            const matchedFlowRateInletVelocityFaPatchVectorField&,
            const faPatch&,
            const DimensionedField<vector, areaMesh>&,
            const faPatchFieldMapper&
        );

        //- Construct as copy
        matchedFlowRateInletVelocityFaPatchVectorField
        (
            const matchedFlowRateInletVelocityFaPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<faPatchVectorField> clone() const
        {
            return tmp<faPatchVectorField>
            (
                new matchedFlowRateInletVelocityFaPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        matchedFlowRateInletVelocityFaPatchVectorField
        (
            const matchedFlowRateInletVelocityFaPatchVectorField&,
            const DimensionedField<vector, areaMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<faPatchVectorField> clone
        (
            const DimensionedField<vector, areaMesh>& iF
        ) const
        {
            return tmp<faPatchVectorField>
            (
                new matchedFlowRateInletVelocityFaPatchVectorField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
