////////////////////////////////////////////////////////////////////////////
// Suspended sediment transport calculations
// by specifying a equilibrium reference concentration at 
// a reference height.
// It is a fixedGradient BC at wall (see notes). 
//
// This option ignores the reference height (delta). If the near-wall
// cell is much larger than the reference height, then we can ignore
// the effect of reference height in the FVM integration. As a result,
// we can just apply the equilibrium reference concentration Cb* to 
// the wall. To use this approach, we also need to use a non-zero
// eddy diffusivity for suspended sediment at the wall. 
//
// More rigorously, the assumption should be that the Peclet number
// Ws*delta/(nu+nut) is << 1. (see Liu 2013, JHE paper).
//
// ATTENTION: Application of this option needs to check the validity
//            of the assumption that the near-wall cell is much larger
//            then the reference height. 
//            As such, this option might be only good if the sediment bed
//            is modeled by a wall function approach.  
// 
// If the fluid part is LES and uses nuSgsUSpaldingWallFunction, then nut_w 
// is readily available since it is derived from a fixedValue type. 
//
// On another note, if the reference height is small, then the eddy 
// diffusiivty at the reference height is not too different from the 
// value at the wall. This also justifies of applying the Cb* right 
// at the wall to estimate the entrainmnet flux.
////////////////////////////////////////////////////////////////////////////

     //Calculate equilibrium concentration Cb_star using any formula.
     //The choice of the formula is specified in the boundary condition 
     //section of Cb_star field.
     //See Garcia (2008): Sedimentation Manual, for a comprehensive list. 
     Cb_star.correctBoundaryConditions();

     //Since we ignore the effect of reference height, 
     //Cw_star equals to Cw_star and 
     //Calculate the equilibrium concentration Cw_star at the wall
     forAll(sedimentPatches, sedPatchI)
     {
        word currentSedPatch = sedimentPatches[sedPatchI];

        label patchI = mesh.boundaryMesh().findPatchID(currentSedPatch);
        if(patchI == -1)
        {
           FatalErrorIn("updateSedimentBC.H")
              << "Didn't find the sediment patch named" << currentSedPatch << nl
              << abort(FatalError);
        }

        //Cw_star at the wall
        fvPatchField<scalar>& Cw_star_wall = Cw_star.boundaryFieldRef()[patchI];

        //Cb_star at the wall
        fvPatchField<scalar>& Cb_star_wall = Cb_star.boundaryFieldRef()[patchI];

        Cw_star_wall = Cb_star_wall;

//        Info << "Cw_star_wall=" << Cw_star_wall << endl;
     }

     //correct the sediment bottom boundary conditon for suspended sediment
     //In this case, the normal gradient of s.s. is specified (see notes)
     forAll(sedimentPatches, sedPatchI)
     {
        word currentSedPatch = sedimentPatches[sedPatchI];

        label patchI = mesh.boundaryMesh().findPatchID(currentSedPatch);
        if(patchI == -1)
        {
           FatalErrorIn("updateSedimentBC.H")
              << "Didn't find the sediment patch named" << currentSedPatch << nl
              << abort(FatalError);
        }

        //settling velocity at current sediment patch
        const fvPatchField<vector>& Ws_wall = Ws.boundaryField()[patchI];
        
        //Ws&n: normal component of the settling velocity to the wall
        scalarField Wsn = Ws_wall & Ws_wall.patch().nf();

        //Cw_star at the wall
        const fvPatchField<scalar>& Cw_star_wall = 
                        Cw_star.boundaryField()[patchI];

        //Cb_star at the wall
        //const fvPatchField<scalar>& Cb_star_wall = 
        //                Cb_star.boundaryField()[patchI];

        //nu and nut at the wall
        volScalarField nu = turbulence->nu();
        volScalarField nut = (turbulence->nuEff() - turbulence->nu());
        const fvPatchField<scalar>& nu_wall = nu.boundaryField()[patchI];
        const fvPatchField<scalar>& nut_wall = nut.boundaryField()[patchI];

        fixedGradientFvPatchScalarField&  TC_wall = 
           refCast<fixedGradientFvPatchScalarField>(Cs.boundaryFieldRef()[patchI]);
         
        scalarField& gradTC = TC_wall.gradient();

        scalarField yp = 1.0/TC_wall.patch().deltaCoeffs();

//        Info << nu << endl;
//        Info << nut_wall << endl;
//        Info << Wsn << endl;
//        Info << Cw_star_wall << endl;

        forAll(gradTC,faceI)
        {
          //we need use nu_total at the wall to entrain the sediment
          scalar nu_total = nu_wall[faceI]+nut_wall[faceI]; 

          gradTC[faceI] = (Wsn[faceI]/nu_total*Cw_star_wall[faceI]); 

//          Info << "gradTC = " << gradTC << endl;

          //calculate the flux through wall
//          scalar Cyp = C.boundaryField()[patchI].internalField()[faceI];

          //update the boundary values
          Cs.boundaryFieldRef().evaluate();
          scalar Cw = Cs.boundaryField()[patchI][faceI];

          suspSedFlux.boundaryFieldRef()[patchI][faceI] =
//             Wsn[faceI]*Cw + nu_total*(Cyp-Cw)/yp[faceI];
             Wsn[faceI]*Cw - nu_total*gradTC[faceI];
        }
     }

