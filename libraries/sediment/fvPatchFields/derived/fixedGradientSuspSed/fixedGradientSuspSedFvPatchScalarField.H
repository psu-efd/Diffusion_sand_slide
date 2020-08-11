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

Class
    Foam::fixedGradientSuspSedFvPatchScalarField

Description
    Foam::fixedGradientSuspSedFvPatchScalarField

SourceFiles
    fixedGradientSuspSedFvPatchScalarField.C

Author
    Zheyu Zhou, University of Delaware, zzhou@udel.edu

Date
    04/03/2014

\*---------------------------------------------------------------------------*/

#ifndef fixedGradientSuspSedFvPatchScalarField_H
#define fixedGradientSuspSedFvPatchScalarField_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class fixedGradientSuspSedFvPatch Declaration
\*---------------------------------------------------------------------------*/

class fixedGradientSuspSedFvPatchScalarField
:
    public fixedGradientFvPatchScalarField
{
    // Private data



public:

    //- Runtime type information
    TypeName("fixedGradientSuspSed");


    // Constructors

        //- Construct from patch and internal field
        fixedGradientSuspSedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedGradientSuspSedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping the given fixedGradientSuspSedFvPatchScalarField
        //  onto a new patch
        fixedGradientSuspSedFvPatchScalarField
        (
            const fixedGradientSuspSedFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fixedGradientSuspSedFvPatchScalarField
        (
            const fixedGradientSuspSedFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<scalar> > clone() const
        {
            return tmp<fvPatchField<scalar> >
            (
                new fixedGradientSuspSedFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        fixedGradientSuspSedFvPatchScalarField
        (
            const fixedGradientSuspSedFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<scalar> > clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<scalar> >
            (
                new fixedGradientSuspSedFvPatchScalarField(*this, iF)
            );
        }


    // Member functions


        // Evaluation functions

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType=Pstream::blocking
            );

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<scalar> > valueBoundaryCoeffs
            (
                const tmp<scalarField>&
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


#endif

// ************************************************************************* //
