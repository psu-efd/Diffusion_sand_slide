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
    Foam::incompressible::RASModels::equilibriumCb_starFvPatchScalarField

Description
    Boundary condition for calculating equilibrium suspended sediment
    concentration Cb_star at a reference height above sediment bed

    Though Cb_star is a concentration at some distance (delta) above the bed, 
    it is store in the bed patch.

SourceFiles
    equilibriumCb_starFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef equilibriumCb_starFvPatchScalarField_H
#define equilibriumCb_starFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class equilibriumCb_starFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class equilibriumCb_starFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
protected:

    // Protected data


    // Protected Member Functions

        //- Check the type of the patch
        virtual void checkType();

        //- Calculate equilibrium concentration at delta
        virtual tmp<scalarField> calcCb_star() const;

        //- Write local wall function variables
        virtual void writeLocalEntries(Ostream&) const;


public:

    //- Runtime type information
    TypeName("equilibriumCb_star");


    // Constructors

        //- Construct from patch and internal field
        equilibriumCb_starFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        equilibriumCb_starFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  equilibriumCb_starFvPatchScalarField
        //  onto a new patch
        equilibriumCb_starFvPatchScalarField
        (
            const equilibriumCb_starFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        equilibriumCb_starFvPatchScalarField
        (
            const equilibriumCb_starFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new equilibriumCb_starFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        equilibriumCb_starFvPatchScalarField
        (
            const equilibriumCb_starFvPatchScalarField&,
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
                new equilibriumCb_starFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
