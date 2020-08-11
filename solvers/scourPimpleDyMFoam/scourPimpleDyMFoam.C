/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
  \\    /   O peration     |
  \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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
  scourPimpleDyMFoam.C

  Description
  Scour model based on pimpleDymFoam

  The Enxer equation is solved with the finite area method.

  Author
  Xiaofeng Liu
  Penn State University

  \*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "faCFD.H"
#include "dynamicFvMesh.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "CorrectPhi.H"
#include "fvOptions.H"
#include "pointPatchFields.H"

#include "primitivePatchInterpolation.H"
#include "wallFvPatch.H"

//#include "areaFieldsTools.H"

#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
        (
         "Transient solver for incompressible, turbulent flow"
         " of Newtonian fluids on a moving mesh."
        );

#   include "postProcess.H"

#   include "addCheckCaseOptions.H"
#   include "setRootCaseLists.H"
#   include "createTime.H"
#   include "createDynamicFvMesh.H"
#   include "createFaMesh.H"          //finite area mesh
#   include "initContinuityErrs.H"
#   include "createDyMControls.H"
#   include "createFields.H"
#   include "createSedimentFields.H"  //fields for sediment
#   include "createFaFields.H"        //finite area fields for sediment
#   include "createUfIfPresent.H"
#   include "CourantNo.H"
#   include "setInitialDeltaT.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    //initialization of the counter for morphological step
    label bedCount = 0;

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
#       include "readDyMControls.H"
#       include "CourantNo.H"
#       include "setDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

#       include "flow.H" 

#       include "sediment.H"

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
