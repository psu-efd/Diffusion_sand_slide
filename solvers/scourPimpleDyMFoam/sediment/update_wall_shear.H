// update the wall shear stress. 
{
    //calculate/update wall shear stress and its magnitude
    const fvPatchList& patches = mesh.boundary();

    forAll(patches, patchi)
    {
        const fvPatch& currPatch = patches[patchi];

        if (isA<wallFvPatch>(currPatch))
        {
            const volSymmTensorField Reff = turbulence->devReff();

            const vectorField& Sfp = mesh.Sf().boundaryField()[patchi];
            const scalarField& magSfp = mesh.magSf().boundaryField()[patchi];
            const symmTensorField& Reffp = Reff.boundaryField()[patchi];

            wallShearStress.boundaryFieldRef()[patchi] = ((Sfp/magSfp) & Reffp); 

            //need to remove the component of wall shear stress along the wall face normal
            //vector, such that wall shear stress is always parallel to the wall
            const vectorField& faceUnitNormal =
                bedFaceUnitNormal.boundaryField()[patchi];

            wallShearStress.boundaryFieldRef()[patchi] -=
                (wallShearStress.boundaryFieldRef()[patchi] & faceUnitNormal)*faceUnitNormal;

            //update the wall shear stress magnitude
            wallShearStressMag.boundaryFieldRef()[patchi] =
                mag(wallShearStress.boundaryField()[patchi]);
        }
    }

    //calculate/update Shields paramter for each sediment patch
    forAll(sedimentPatches, sedPatchI)
    {
        word currentSedPatch = sedimentPatches[sedPatchI];
        Info << "calculating Shields paramter on the patch named "
            << currentSedPatch << endl;

        label patchI = mesh.boundaryMesh().findPatchID(currentSedPatch);
        if(patchI == -1)
        {
            FatalErrorIn("update_wall_shear.H")
                << "Didn't find the sediment patch named" << currentSedPatch << nl
                << abort(FatalError);
        }

        //wallShearStressMag at current sediment patch
        const fvPatchField<scalar>& wallShearMag_wall =
            wallShearStressMag.boundaryField()[patchI];

        //calculate/update Shields paramter
        fvPatchField<scalar>& Shields_wall =
            wallShieldsNumber.boundaryFieldRef()[patchI];

        forAll(Shields_wall,faceI)
        {
            Shields_wall[faceI] = (wallShearMag_wall[faceI]/
                    (mag(g)*(rhos-rhow)/rhow*diam.value())).value();
        }
    }
}
