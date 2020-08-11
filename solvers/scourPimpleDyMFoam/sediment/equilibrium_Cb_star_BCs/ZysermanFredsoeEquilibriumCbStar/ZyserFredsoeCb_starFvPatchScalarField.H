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
    ZysermanFredsoeCb_starFvPatchScalarField

Description
    Boundary condition to calculate the equilibrium concentration at reference
    height (delta=2*diameter) for suspended sediment using Zyserman and Fredsoe (1994)

SourceFiles
    ZyserFredsoeCb_starFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef ZyserFredsoeCb_starFvPatchScalarField_H
#define ZyserFredsoeCb_starFvPatchScalarField_H

#include "equilibriumCb_starFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class ZyserFredsoeCb_starFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class ZyserFredsoeCb_starFvPatchScalarField
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
    TypeName("ZyserFredsoeCb_star");


    // Constructors

        //- Construct from patch and internal field
        ZyserFredsoeCb_starFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        ZyserFredsoeCb_starFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  ZyserFredsoeCb_starFvPatchScalarField
        //  onto a new patch
        ZyserFredsoeCb_starFvPatchScalarField
        (
            const ZyserFredsoeCb_starFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        ZyserFredsoeCb_starFvPatchScalarField
        (
            const ZyserFredsoeCb_starFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new ZyserFredsoeCb_starFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        ZyserFredsoeCb_starFvPatchScalarField
        (
            const ZyserFredsoeCb_starFvPatchScalarField&,
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
                new ZyserFredsoeCb_starFvPatchScalarField(*this, iF)
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
