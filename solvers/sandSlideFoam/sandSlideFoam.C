/*---------------------------------------------------------------------------*\
    =========                 |
    \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
    \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    sandSlideFoam 

    Description
    Solves a slope-limited diffusion equation for sand slide.

    Authors
    Xiaofeng Liu and Yalan Song
    Penn State University

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"
#include "unitConversion.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    //This solver can not be used for 3D meshes
    if(mesh.nGeometricD() == 3)
    {
        FatalError
           << "This sand slide solver can not be used with a 3D mesh!"
           << abort(FatalError);
    }

    simpleControl simple(mesh);

    #include "createFields.H"


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
    Info<< "\nCalculating sand slide induced bed elevation evolution\n" << endl;

    scalarField oldEta = eta*0 ;
   
 
    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        scalar friction = Foam::tan(angle);   //static friction angle
     
        label counter = 0;
        
        //fvScalarField oldEta = eta*0 ;

         
        if(gradientScheme.match("LGS"))
        {
            const labelList& faceOwner = mesh.faceOwner();
            const labelList& faceNeighbour = mesh.faceNeighbour();

            forAll(K,I)
            {
                point owner(mesh.cellCentres()[faceOwner[I]]);
                point neigh(mesh.cellCentres()[faceNeighbour[I]]);
                scalar length=sqr(owner.x()-neigh.x())+sqr(owner.y()-neigh.y());
                length=Foam::sqrt(length);
                scalar tanPhi=eta[faceOwner[I]]-eta[faceNeighbour[I]];
                tanPhi=mag(tanPhi/length);

                if(mag(tanPhi)<friction)   //linear gradient
                {
                    K.primitiveFieldRef()[I]=0;
                }
                else
                {
                    K.primitiveFieldRef()[I]=K0.value();
                    counter += 1;
                }
           }
        }
        else if(gradientScheme.match("GNGS") or gradientScheme.match("LGS+GNGS"))
        {
	    volVectorField gradEta=fvc::grad(eta);
            surfaceVectorField sGradEta=linearInterpolate(gradEta);
	    forAll(K,I)
           {
                if(mag(sGradEta[I])<friction)   //Gauss gradient
                {
                    K.primitiveFieldRef()[I]=0;
                }
                else
                {
                    K.primitiveFieldRef()[I]=K0.value();
                    counter += 1;
                }
            }
        }
        else
        {
           FatalErrorIn("gradientScheme")
              << "The specified gradientScheme is not valid." << nl
              << abort(FatalError);

        }
           
        Info<< "Total number of faces with exceeding slopes: "
            << counter << endl;

        scalarField dEta = mag(eta-oldEta);

        Info << "    Maximum bed elevation change in this subiteration = "
             << max(dEta) << endl;
          
        if(gradientScheme.match("LGS+GNGS") && (max(dEta)<1e-5))
        {
            Info << "	Switch to LGS to remove residuals   "<< endl;
            gradientScheme = "LGS";
        }

        if(counter == 0 || (max(dEta)<1e-9) )
        {
            break;
        }
 
        oldEta = eta;
       

             
        fvScalarMatrix newEtaEqn
        (
            fvm::ddt(eta)-fvm::laplacian(K, eta) 
        );
        newEtaEqn.solve();

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
         
    }
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
