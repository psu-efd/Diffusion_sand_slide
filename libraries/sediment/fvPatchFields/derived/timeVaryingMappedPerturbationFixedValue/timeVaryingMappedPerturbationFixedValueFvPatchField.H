/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::timeVaryingMappedPerturbationFixedValueFvPatchField

Group
    grpInletBoundaryConditions grpCoupledBoundaryConditions

Description
    This boundary conditions interpolates the values from a set of supplied
    points in space and time.

    Supplied data should be specified in constant/boundaryData/\<patchname\>/
      - points             : pointField of locations
      - \<time\>/\<field\> : field of values at time \<time\>

    The default mode of operation (mapMethod planarInterpolation) is to project
    the points onto a plane (constructed from the first threee points) and
    construct a 2D triangulation and finds for the face centres the triangle it
    is in and the weights to the 3 vertices.

    The optional mapMethod nearest will avoid all projection and triangulation
    and just use the value at the nearest vertex.

    Values are interpolated linearly between times.

Usage
    \table
        Property     | Description                      | Required | Default
        setAverage   | Use average value                | no    | false
        perturb      | Perturb points for regular geometries | no | 1e-5
        points       | Name of points file              | no    | points
        fieldTable   | Alternative field name to sample | no | this field name
        mapMethod    | Type of mapping              | no | planarInterpolation
        offset       | Offset to mapped values      | no | Zero
    \endtable

    \verbatim
    <patchName>
    {
        type            timeVaryingMappedPerturbationFixedValue;
    }
    \endverbatim

See also
    Foam::fixedValueFvPatchField
    Foam::Function1Types

SourceFiles
    timeVaryingMappedPerturbationFixedValueFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef timeVaryingMappedPerturbationFixedValueFvPatchField_H
#define timeVaryingMappedPerturbationFixedValueFvPatchField_H

#include "fixedValueFvPatchFields.H"
#include "MappedFile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class timeVaryingMappedPerturbationFixedValueFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class timeVaryingMappedPerturbationFixedValueFvPatchField
:
    public fixedValueFvPatchField<Type>
{
    // Private data

        autoPtr<PatchFunction1Types::MappedFile<Type>> uniformValue_;

        //- Random number generator
        //Random ranGen_;

        //- Fluctuation scale (amplitude)
        scalar fluctuationScale_;

        //- Fluctuation period
        scalar T_;

        //- Current time index (used for updating)
        label curTimeIndex_;



public:

    //- Runtime type information
    TypeName("timeVaryingMappedPerturbationFixedValue");


    // Constructors

        //- Construct from patch and internal field
        timeVaryingMappedPerturbationFixedValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        timeVaryingMappedPerturbationFixedValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given timeVaryingMappedPerturbationFixedValueFvPatchField
        //  onto a new patch
        timeVaryingMappedPerturbationFixedValueFvPatchField
        (
            const timeVaryingMappedPerturbationFixedValueFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        timeVaryingMappedPerturbationFixedValueFvPatchField
        (
            const timeVaryingMappedPerturbationFixedValueFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new timeVaryingMappedPerturbationFixedValueFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        timeVaryingMappedPerturbationFixedValueFvPatchField
        (
            const timeVaryingMappedPerturbationFixedValueFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>
            (
                new timeVaryingMappedPerturbationFixedValueFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<Type>&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "timeVaryingMappedPerturbationFixedValueFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
