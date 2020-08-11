//if sand slide model is turned on 
if(sandSlide)
{
    areaScalarField faZb_before_sand_slide = faZb;
   
    areaScalarField oldFaZb = faZb*0 ;

    //subiteration loop for sand slide model to converge 
    for(int sandSlideIter=0;sandSlideIter<maxSandSlideIter;sandSlideIter++)
    {
        Info << nl <<  "Sand slide iteration # " << sandSlideIter << " ..." << endl;

        //update the elevation gradient at edge
#       include "update_edge_gradZb.H"
        areaScalarField faDzb = mag(faZb-oldFaZb);          
         
        Info << "    Maximum bed elevation change in this subiteration = " 
             << gMax(faDzb) << endl;

        //check convergence, i.e., whether all bed angles are below
        //ange of repose
        if(!checkSandSlide || (gMax(faDzb)<1e-8 && sandSlideIter!=0) )
        {
            Info << "    Sand slide step converged after " << sandSlideIter << " iterations." << endl;
           
            break;
        }

        //update the diffusion coefficient for sand slide model
        //If the bed elevation gradient (=tan(bed angle)) is larger than the angle of repose,
        //K=K0; otherwise, K=0.
         
        oldFaZb = faZb;

        faeSandSlideK = K0*pos(mag(faeGradZb) - Foam::tan(angleOfRepose.value()));
   
        scalar rDeltaT = 1.0/runTime.deltaT().value();
        faScalarMatrix DDT(faZb,faZb.dimensions()*dimArea/dimTime);
        DDT.diag()=  rDeltaT*aMesh.S();
        DDT.source() = rDeltaT*oldFaZb.primitiveField()*aMesh.S();

        //define diffusion equation for sand slide
        faScalarMatrix zbEqn
            (
              DDT == fam::laplacian(faeSandSlideK, faZb)
            );
        zbEqn.solve();
        
    }

    //update the bed change due to sand slide
    faDeltaZ_sand_slide = faZb - faZb_before_sand_slide;
}

Info <<faDeltaZFlag_sand_slide.dimensions()<<endl;

//count occurrance of sand slide
forAll(faDeltaZ_sand_slide, I)
{       
    if(mag(faDeltaZ_sand_slide[I])>0)
        {
	      faDeltaZFlag_sand_slide[I] += 1;
        }
}

//update the bed elevation change after the sand slide step
//Now faDeltaZ stores the total change in this sediment step
//due to Exner and sand slide.
faDeltaZ = faZb - faZb_org;
faDeltaZ.correctBoundaryConditions();

Info << nl;
Info << "In this sediment step, the max and min bed elevation changes are = " 
<< gMax(faDeltaZ.internalField()) << ", " 
<< gMin(faDeltaZ.internalField()) << endl;

//after all above steps, we finally need to map te bed elevation change from
//faMesh to fvMesh.
//map faDeltaZ from faMesh to fvMesh
vsm.mapToVolume(faDeltaZ, bedCellHeightChange.boundaryFieldRef());
vsm.mapToVolume(faDeltaZ_sand_slide, bedChange_sand_slide.boundaryFieldRef());
vsm.mapToVolume(faDeltaZFlag_sand_slide, bedChangeFlag_sand_slide.boundaryFieldRef());

//accumulate the bed elevation change
fluidFaceCenterTotalDeltaY += bedCellHeightChange;
