/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

Application
    setBedShape 

Description
    set the shape (elevation) of the bed. As of now, it is only used for
    column collapse case.
    

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "readSetBedShapeDict.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    //Assume the mesh 2D in x and y direction. z is the vertical.

    //define the initial sand column on the bed
    //assuming the center is at the origin (0,0).
    forAll(eta.primitiveFieldRef(),I)
    {
       scalar R=Foam::sqrt( Foam::pow(mesh.cellCentres()[I].x(),2)
                           +Foam::pow(mesh.cellCentres()[I].y(),2));
 
       if (R<=R0)
       {
           eta.primitiveFieldRef()[I]=hi.value();
       }
    }

    eta.write();

    Info<< "end" << endl;

    return 0;
}

// ************************************************************************* //
