Info << "Solving flow ..." << endl;

// --- Pressure-velocity PIMPLE corrector loop
while (pimple.loop())
{
    if (pimple.firstIter() || moveMeshOuterCorrectors)
    {
        mesh.update();
        //check mesh validicity after mesh motion
        if(calculateSed && checkMesh)
        {
            if(mesh.checkMesh(true))
            {
                runTime.writeNow();

                FatalErrorIn("flow.H")
                    << "Mesh is not valid anymore after mesh motion. Simulation can not continue." << nl
                    << abort(FatalError);            
            }
        }

        if (mesh.changing())
        {
            MRF.update();

            if (correctPhi)
            {
                // Calculate absolute flux
                // from the mapped surface velocity
                phi = mesh.Sf() & Uf();

#               include "correctPhi.H"

                // Make the flux relative to the mesh motion
                fvc::makeRelative(phi, U);
            }

            if (checkMeshCourantNo)
            {
#               include "meshCourantNo.H"
            }
        }
    }

#   include "UEqn.H"

    // --- Pressure corrector loop
    while (pimple.correct())
    {
#       include "pEqn.H"
    }

    if (pimple.turbCorr())
    {
        laminarTransport.correct();
        turbulence->correct();
    }
}

//for each of the erodible sediment patch, update the 
//height. Only relevant on the erodiable boundaries and 
//only if sediment transport is turned on.
if(calculateSed) //if sediment transport is turned on
{
    forAll(sedimentPatches, sedPatchI)
    {
        word currentSedPatch = sedimentPatches[sedPatchI];

        label movingPatchID = mesh.boundaryMesh().findPatchID(currentSedPatch);
        if(movingPatchID == -1)
        {
            FatalErrorIn("bed_load.H")
                << "Didn't find the sediment patch named" << currentSedPatch << nl
                << abort(FatalError);
        }

        const fvPatch& erodibleBottomfvPatch = mesh.boundary()[movingPatchID];
        const polyPatch& erodibleBottomPolyPatch = erodibleBottomfvPatch.patch();
        const vectorField::subField faceCentres = erodibleBottomPolyPatch.faceCentres();

        //Info << "sediment patch  = " << currentSedPatch << endl;
        //Info << "faceCentres = " << faceCentres << endl;

        fluidFaceCenterHeigh.boundaryFieldRef()[movingPatchID] ==  faceCentres.component(vector::Y);

        //Info << "tt = " << fluidFaceCenterHeigh.boundaryFieldRef()[movingPatchID] << endl;
    }
}
