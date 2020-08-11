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
    geomSlideFoam 

    Description
    Solve the sand slide using the geometric sand slide method.

    Reference
    Khosronejad, A., Kang, S., Borazjani, I., and Sotiropoulos, F. (2011)
    Curvilinear immersed boundary method for simulating coupled flow and 
    bed morphodynamic interactions due to sediment transport phenomena.
    Advances in Water Resources, 34(7): 829-843

    Authors
    Xiaofeng Liu and Yalan Song
    Penn State University


    \*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"
#include "unitConversion.H"
#include "SortableList.H"

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
    
    Info<< "\nCalculating elevation distribution\n" << endl;
    
    scalarField oldEta = eta*0 ;


    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        scalarField dEta = mag(eta-oldEta);

        Info << "    Maximum bed elevation change in this subiteration = "
             << max(dEta) << endl;


        oldEta = eta;
        

  
        scalar mus = Foam::tan(angle);        //static friction angle

       

        const fvMesh& mesh=eta.mesh();        

        const labelListList& cellCells = mesh.cellCells();  //neigbouring cells
	const pointField& C = mesh.C();                     //cell centers
	const scalarField& Ah = mesh.V();                   //cell volume (projected area A*h, 
                                                            //h is a constant height)
	   
        label slopeCounter=0;                               // record  
        volVectorField gradEta = fvc::grad(eta);   
        volScalarField magGradEta = mag(gradEta);           //gradient of each cell
 

          
        forAll(C,I)             // loop over all cells
        {
            point C_p=C[I];
            C_p.z()=0.0;
            const labelList& myCells = cellCells[I];
            scalarField Spi(myCells.size(),0.0);
            scalarField Lpi(myCells.size(),0.0);
            scalar Ci=0.0;
            scalar Cp=Ah[I];   

            forAll(myCells,II)
            {
                label neighID=myCells[II];
                point C_n=C[neighID];
                C_n.z()=0.0;
                Lpi[II]=mag(C_p-C_n);
                Spi[II]=(eta[I]-eta[neighID])/Lpi[II];

                if(mag(Spi[II])>mus)                    
                {
                    slopeCounter++;
                }

                Spi[II]=max(-mus,Spi[II]);
                Spi[II]=min(mus,Spi[II]);
                Ci+=Ah[neighID]*(eta[I]-eta[neighID]-Lpi[II]*Spi[II]); // Eq(17)
                Cp+=Ah[neighID];
            }

            scalar totalMass=0.0;
            scalar dEtap=0.0;
            scalar oldEdap=eta[I];

            if(mag(Cp)>SMALL)
            {
                dEtap=-Ci/Cp;                                          //Eq(18)
            }
            totalMass+=Ah[I]*dEtap;
               
            eta[I]+=dEtap;
                
            forAll(myCells,II)
            {
                label neighID=myCells[II];
                scalar dEtan=oldEdap+dEtap-eta[neighID]-Lpi[II]*Spi[II];
                eta[neighID]+=dEtan;
                totalMass+=Ah[neighID]*dEtan;
            } 
                
            if(mag(totalMass/Ah[I])>0.000001)
            {
                //check Mass balance for each cell         
                Info<<totalMass/Ah[I]<<tab<<eta[I]<<tab<<dEtap<<endl;           
            }
                
        }

        if(slopeCounter == 0 || (max(dEta)<1e-9) )
        {
            break;
        }
           

        runTime.write();
           
        Info<< "total number of large slope cells: "
            << slopeCounter << endl;
            
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
