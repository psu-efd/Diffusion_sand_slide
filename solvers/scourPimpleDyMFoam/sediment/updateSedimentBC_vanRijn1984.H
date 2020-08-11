////////////////////////////////////////////////////////////////////////////
// Suspended sediment transport calculations
// using van Rijn (1984) pick up function at the wall.
// It is a fixedGradient BC at wall (see notes). 
//
// This approach, which is simple, has been used in several
// papers, such as Zedler and Street (2001, 2006).
// However, it is not clear to me the rationale of using this "AT" the wall.
// van Rijn himself used this pick up function at certain reference
// height away from the wall. 
////////////////////////////////////////////////////////////////////////////

     //Calculate equilibrium concentration Cb_star
     //Since the van Rijn (1984) pickup function P is used at the 
     //wall, it is a bit unusual. Here Cb_star equals to Cw_star and 
     //it is calculated as Cb_star=Cw_star=P/Ws. 
     Cb_star.correctBoundaryConditions();

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
        //const fvPatchField<scalar>& nut_wall = nut.boundaryField()[patchI];

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
          scalar nu_total = nu_wall[faceI]; //only molecular viscosity

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

