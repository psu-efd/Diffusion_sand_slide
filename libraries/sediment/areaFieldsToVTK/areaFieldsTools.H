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

#ifndef areaFieldsTools_H
#define areaFieldsTools_H

#include "areaFields.H"
#include "OFstream.H"

#include "foamVtkSurfaceMeshWriter.H"
#include "foamVtkSurfaceFieldWriter.H"
#include "foamVtkSeriesWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//write areaLabelField, areaScalarField,
//areaVectorField to vtk (one field per file)

//write the area to VTK
void writeFaMesh2VTK
(
  const faMesh& surf,
//  OFstream& os,
  const Time& runTime
);

void areaScalarField2VTK
(
   const faMesh& surf,
   const areaScalarField& aslf,
   OFstream& os
);

void areaVectorField2VTK
(
   const faMesh& surf,
   const areaVectorField& aslf,
   OFstream& os
);

//write PrimitivePatch to VTK
void writePrimitivePatch2VTK
(
 const PrimitivePatch<face, IndirectList, const pointField&>& pp,
 OFstream& os
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
