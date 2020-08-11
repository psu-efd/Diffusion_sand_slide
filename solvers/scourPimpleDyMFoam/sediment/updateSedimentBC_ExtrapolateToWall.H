////////////////////////////////////////////////////////////////////////////
// Suspended sediment transport calculations
// by specifying a equilibrium reference concentration at 
// a reference height.
// It is a fixedGradient BC at wall (see notes). 
//
// This option extrapolate the equilibrium reference concentration to the
// wall. So the effect of the reference height is considered. This option 
// should be used if Peclet number, Ws*delta/(nu+nut), is large. This condition
// is usually true for large particles and weak turbulence. 
//
// As a result, we have Cw* extraploated from Cb*. 
// To use this approach, we also need to use a proper
// eddy diffusivity for suspended sediment at the wall. Theoretically,
// eddy viscosity should be zero at the wall (no turbulence) and the only
// contribution is from molecular viscosity. However, if only molecular
// viscosity is used for the diffusive flux, the gradient BC at the wall
// will have a HUGE value. It also somewhat contradicts to the discrete nature
// of the FVM code. In a control volume, we assume the eddy viscosity is 
// constant or linearly distributed. 
//
// If the fluid part is LES and uses nuSgsUSpaldingWallFunction, then nut_w 
// if readily available since it is derived from a fixedValue type. 
//
// Reference: Liu (2013), JHE.
////////////////////////////////////////////////////////////////////////////

     //Calculate equilibrium concentration Cb_star using any formula.
     //The choice of the formula is specified in the boundary condition 
     //section of Cb_star field.
     //See Garcia (2008): Sedimentation Manual, for a comprehensive list. 
     Cb_star.correctBoundaryConditions();

     //Calculate the equilibrium concentration Cw_star at the wall
     //by extrapolation
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

        //nu and nut at the wall
        volScalarField nu = turbulence->nu();
        volScalarField nut = (turbulence->nuEff() - turbulence->nu());
        const fvPatchField<scalar>& nu_wall = nu.boundaryField()[patchI];
        const fvPatchField<scalar>& nut_wall = nut.boundaryField()[patchI];


        //diameter at current sediment patch
        scalarField diam_wall(Ws.boundaryField()[patchI].size(), diam.value());

        //settling velocity at current sediment patch
        const fvPatchField<vector>& Ws_wall = Ws.boundaryField()[patchI];

        //the value resulting from integrating from reference height to wall
        scalarField integral_temp(Ws_wall.size(),0.0);

        //the reference height (assume = 2*diameter, needs to be changed if otherwise)
        scalarField delta = 2*diam_wall;

        scalarField Wsn = Ws_wall & Ws_wall.patch().nf();

        //Cw_star at the wall
        fvPatchField<scalar>& Cw_star_wall = Cw_star.boundaryFieldRef()[patchI];

        //Cb_star at the wall
        fvPatchField<scalar>& Cb_star_wall = Cb_star.boundaryFieldRef()[patchI];

        forAll(integral_temp,faceI)
        {
           scalar nu_total = nu_wall[faceI] + nut_wall[faceI];
           integral_temp[faceI] = ::exp(Wsn[faceI]/nu_total*delta[faceI]);
        }

        Cw_star_wall = Cb_star_wall*integral_temp;

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

        //diameter at current sediment patch
        scalarField diam_wall(Ws_wall.size(), diam.value());

        //the reference height (assume = 2*diameter, needs to be changed otherwise)
        scalarField delta = 2*diam_wall;

        //Cw_star at the wall
        const fvPatchField<scalar>& Cw_star_wall = 
                        Cw_star.boundaryField()[patchI];

        //Cb_star at the wall
        const fvPatchField<scalar>& Cb_star_wall = 
                        Cb_star.boundaryField()[patchI];

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

          //exponential formula 
          if(extrapolationOption.match("exponential"))
          {
              gradTC[faceI] = (Wsn[faceI]/nu_total*Cw_star_wall[faceI]);
          }
          else if(extrapolationOption.match("linearization"))
          {
              //linearization (FVM)
//              scalar temp1 = Wsn[faceI]/nu_total;
//              scalar temp2 = Wsn[faceI]*delta[faceI]/nu_total;
//              gradTC[faceI] = temp1/(1.0-temp2)*Cb_star_wall[faceI];

              scalar temp1 = Wsn[faceI]/nu_total;
              scalar temp2 = Wsn[faceI]*delta[faceI]/nu_total;
              gradTC[faceI] = temp1*Cb_star_wall[faceI]*(1.0+temp2);
          }
          else
          {
            FatalErrorIn("extrapolationOption")
               << "The option for extrapolationOption in transportDict: " 
               << extrapolationOption
               << " is not valid." << nl
               << "Valid options are: exponential and linearization "
              << abort(FatalError);
          }

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

