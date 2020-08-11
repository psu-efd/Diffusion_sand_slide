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

Class
    vanRijn1984PickupFunctionCb_starFvPatchScalarField

Description
    Boundary condition to calculate the equilibrium concentration
 
    It uses the van Rijn (1984) pickup function P. Though it inherits 
    from equilibriumCb_starFvPatchScalarField in this BC implementation, 
    it is meant to be used at the wall. That is the reason with "Cw" in
    its name. This approach is somewhat unusual. 
    It has been used by others, such as Zedler and Street (2001, 2006).

    See notes for details.

SourceFiles
    vanRijn1984PickupFunctionCb_starFvPatchScalarField.C

Author
    Xiaofeng Liu

\*---------------------------------------------------------------------------*/

#ifndef vanRijn1984PickupFunctionCb_starFvPatchScalarField_H
#define vanRijn1984PickupFunctionCb_starFvPatchScalarField_H

#include "equilibriumCb_starFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class vanRijn1984PickupFunctionCb_starFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class vanRijn1984PickupFunctionCb_starFvPatchScalarField
:
    public equilibriumCb_starFvPatchScalarField 
{
protected:

    // Protected data


    // Protected Member Functions

        //- Calculate the equilibrium concentration Cb_star 
        virtual tmp<scalarField> calcCb_star() const;

public:

    //- Runtime type information
    TypeName("vanRijn1984PickupFunction");


    // Constructors

        //- Construct from patch and internal field
        vanRijn1984PickupFunctionCb_starFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        vanRijn1984PickupFunctionCb_starFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  vanRijn1984PickupFunctionCb_starFvPatchScalarField
        //  onto a new patch
        vanRijn1984PickupFunctionCb_starFvPatchScalarField
        (
            const vanRijn1984PickupFunctionCb_starFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        vanRijn1984PickupFunctionCb_starFvPatchScalarField
        (
            const vanRijn1984PickupFunctionCb_starFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new vanRijn1984PickupFunctionCb_starFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        vanRijn1984PickupFunctionCb_starFvPatchScalarField
        (
            const vanRijn1984PickupFunctionCb_starFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new vanRijn1984PickupFunctionCb_starFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Acces functions


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchScalarField&,
                const labelList&
            );


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
