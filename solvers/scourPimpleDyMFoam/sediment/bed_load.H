// calculate the bed load 

//for each of the erodible sediment patch
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

    //steepest slope of the bed face given surface normal unit vector
    const fvPatch& erodibleBottomfvPatch = mesh.boundary()[movingPatchID];

    //bed faces project 2D area on the x-z plane
    //x-z plane normal vector
    vector xzPlaneNormal(0,1,0);

/*  This gives wrong 2D projected area.
    bedFace2DArea.boundaryFieldRef()[movingPatchID] ==
        mag(erodibleBottomfvPatch.Sf() & xzPlaneNormal) + VSMALL;
*/
    pointField pf2D = erodibleBottomfvPatch.patch().localPoints();
    forAll(pf2D, pointI)
    {
        pf2D[pointI].y() = 0.0; //even the elevation, i.e., project the points to x-z plane
    }

    const faceList& localFaces = erodibleBottomfvPatch.patch().localFaces();

    forAll(erodibleBottomfvPatch, faceI)
    {
        bedFace2DArea.boundaryFieldRef()[movingPatchID][faceI] = 
            localFaces[faceI].mag(pf2D);
    }

    //The sign of Sf is defined w.r.t to fvMesh. Thus it points out of domain (downward)
    //Needs to be flipped.
    bedFaceUnitNormal.boundaryFieldRef()[movingPatchID] == 
        -erodibleBottomfvPatch.Sf()/erodibleBottomfvPatch.magSf();

    //normalize bedFaceUnitNormal
    bedFaceUnitNormal.boundaryFieldRef()[movingPatchID] ==
         bedFaceUnitNormal.boundaryFieldRef()[movingPatchID]
        /mag(bedFaceUnitNormal.boundaryFieldRef()[movingPatchID]);

    //y-component (vertical) of bed face unit normal 
    bedFaceUnitNormaly.boundaryFieldRef()[movingPatchID] ==
         bedFaceUnitNormal.boundaryFieldRef()[movingPatchID].component(vector::Y);

    //Calculate the steepest slope unit vector on the erodible patch
    forAll(erodibleBottomfvPatch, faceI)
    {
        const vector& faceUnitNormal = bedFaceUnitNormal.boundaryField()[movingPatchID][faceI];

        //Info << "faceUnitNormal = " << faceUnitNormal << endl;


        //bed slope angle = acos(n & xzPlaneNormal) assuming both vectors are normalized
        bedWallSlopeAngle.boundaryFieldRef()[movingPatchID][faceI] =
            Foam::acos(faceUnitNormal & xzPlaneNormal);

        //Info << "bedWallSlopeAngle.boundaryFieldRef()[movingPatchID][faceI] = "
        //    << bedWallSlopeAngle.boundaryFieldRef()[movingPatchID][faceI] << endl;

        //modulate the angle
        if(bedWallSlopeAngle.boundaryFieldRef()[movingPatchID][faceI] < 0.0)
        {
            bedWallSlopeAngle.boundaryFieldRef()[movingPatchID][faceI] += pi;
        }

        if(bedWallSlopeAngle.boundaryFieldRef()[movingPatchID][faceI] > pi/2.0)
        {
            bedWallSlopeAngle.boundaryFieldRef()[movingPatchID][faceI] -= pi;
        }

        //make sure the bed slope angle is between 0 and pi/2
        if(!((bedWallSlopeAngle.boundaryFieldRef()[movingPatchID][faceI]>=0.0) &&
               (bedWallSlopeAngle.boundaryFieldRef()[movingPatchID][faceI] <= pi/2.0)))
        {
            FatalErrorIn("bed_load.H")
                << "Bed slope out of range" << nl
                << abort(FatalError);
        }

        //also make sure the bed slope angle is less than the angle of repose in order
        //to have a sane critical Shields number correction due to bed angle
        //Ref. the sqrt(...) in Eqn. (11) in Liu and Garcia (2008). 
        if(bedWallSlopeAngle.boundaryFieldRef()[movingPatchID][faceI]>angleOfRepose.value())
        {
            bedWallSlopeAngle.boundaryFieldRef()[movingPatchID][faceI] = angleOfRepose.value();
        }

        //the steepest slope should be equal to gravity vector with its component
        //in the surface normal direction removed, i.e., 
        //g - (g & n)n
        vector steepestSlopeVector = g.value() - (g.value() & faceUnitNormal)*faceUnitNormal;

        //Info << "steepestSlopeVector = " << steepestSlopeVector << endl;
        //Info << "steepestSlopeVector & g.value() " << (steepestSlopeVector & g.value()) << endl;

        //make sure the steepest slope is pointing downward
        if(!pos0(steepestSlopeVector & g.value()))
        {
            FatalErrorIn("bed_load.H")
                << "Bed's steepest slope vector in wrong direction" << nl
                << abort(FatalError);
        }

        //normalize the steepest slope
        if(mag(steepestSlopeVector) < 1e-6)
        {
            //very small mag of steepestSlopeVector happens for horizontal bed
            //In this case, the steepestSlopeVector does not matter because the 
            //bed angle is very small. Just use a simple vector.
            steepestSlopeVector = vector(1, 0, 0);
        }
        else
        {
            steepestSlopeVector /= mag(steepestSlopeVector);
        }

        bedSteepestSlopeVector.boundaryFieldRef()[movingPatchID][faceI] = steepestSlopeVector;

        //calculate the angle between wall shear (velocity) vector and steepestSlopeVector
        const vector& curWallShear = wallShearStress.boundaryFieldRef()[movingPatchID][faceI];

        if(mag(curWallShear) < 1e-6)
        {
            bedWallPhiAngle.boundaryFieldRef()[movingPatchID][faceI] = 0.0;
        }
        else
        {
            bedWallPhiAngle.boundaryFieldRef()[movingPatchID][faceI] = 
                Foam::acos((curWallShear & steepestSlopeVector)
                    /mag(curWallShear)/mag(steepestSlopeVector));
        }
    }

    //correct the slope effect on the critical Shields number
    //Flow direction: 1->downslope, 0->sideways, -1->upslope
    bedWallFlowDir.boundaryFieldRef()[movingPatchID] ==
        sign((wallShearStress.boundaryField()[movingPatchID] &
                    bedSteepestSlopeVector.boundaryField()[movingPatchID]
             )()
            );
//Info << "tt = " << 
//Info << "tt = " << 
//(1.0-Foam::sqr(Foam::sin(bedWallPhiAngle.boundaryField()[movingPatchID])
//           *Foam::tan(bedWallSlopeAngle.boundaryField()[movingPatchID])/mus.value()))<< endl;

    //update the critical Shields number based on calculated slope and angles
    //see Eq. (11) in Liu and Garcia (2008) or Eq. (25) in Roulund et al. (2005)
    scalarField aField = Foam::cos(bedWallSlopeAngle.boundaryField()[movingPatchID])
        *Foam::sqrt(1.0-Foam::sqr(Foam::sin(bedWallPhiAngle.boundaryField()[movingPatchID])  
           *Foam::tan(bedWallSlopeAngle.boundaryField()[movingPatchID])/mus.value()));

    scalarField bField = Foam::cos(bedWallPhiAngle.boundaryField()[movingPatchID])
        *Foam::sin(bedWallSlopeAngle.boundaryField()[movingPatchID])/mus.value();

    bedWallCriticalShieldsNum.boundaryFieldRef()[movingPatchID] ==
        theta_cri_b.value()*(aField - bField);

/*
    bedWallq0 =
        pos(wallShieldsNumber-bedWallCriticalShieldsNum)*12*sqrt(mag(g)*(rhos/rhow-1)*pow(diam,3))
        *(wallShieldsNumber-bedWallCriticalShieldsNum)*sqrt(wallShieldsNumber);
*/

    bedWallq0 =
         Foam::pos(wallShieldsNumber-bedWallCriticalShieldsNum)*Foam::sqrt(mag(g)*(rhos/rhow-1)*pow(diam,3))
        *18.74*(wallShieldsNumber-bedWallCriticalShieldsNum)
        *(Foam::sqrt(wallShieldsNumber) - 0.7*Foam::sqrt(bedWallCriticalShieldsNum));


    //slope: dh/dx, dh/dy(=0), dh/dz
    forAll(bedSlope.boundaryFieldRef()[movingPatchID], faceI)
    {
        fvPatchVectorField& fpv = bedSlope.boundaryFieldRef()[movingPatchID];

        fpv[faceI].x() = 
            (-bedFaceUnitNormal.boundaryField()[movingPatchID].component(vector::X)/
            bedFaceUnitNormal.boundaryField()[movingPatchID].component(vector::Y))()[faceI];
        fpv[faceI].z() =
            (-bedFaceUnitNormal.boundaryField()[movingPatchID].component(vector::Z)/
            bedFaceUnitNormal.boundaryField()[movingPatchID].component(vector::Y))()[faceI];
    }

    //Bedload transport rate vector, see Eq. (12) in Liu and Garcia (2008)
    bedWallqi = bedWallq0*wallShearStress/(wallShearStressMag+smallWallShearMag)-C*mag(bedWallq0)*bedSlope;

    //map the bedload transport rate vector from volume mesh to faMesh
    faQb.primitiveFieldRef() = vsm.mapToSurface<vector>(bedWallqi.boundaryField());
    faQb.correctBoundaryConditions();
   
    //kill the vertical component of bedload transport rate vector
    faQb -= xzPlaneNormal*(faQb & xzPlaneNormal);

    //map the bed face's unit normal vector from fvMesh to faMesh
    faBedFaceUnitNormaly.primitiveFieldRef() = 
        vsm.mapToSurface<scalar>(bedFaceUnitNormaly.boundaryField());

    faQb_copy = faQb;
    phi_qb_sign = sign(fac::interpolate(faQb_copy) & aMesh.Le());
    phi_qb = fac::interpolate(faQb) & aMesh.Le();
/*
Info << "edge centers = " << aMesh.edgeCentres() << endl;
Info << "faQb = " << faQb << endl;
Info << "fac::interpolate(faQb) = " << fac::interpolate(faQb) << endl;
Info << "aMesh.Le() = " << aMesh.Le() << endl;
Info << "phi_qb = " << phi_qb << endl;
*/
}//end of loop over erodible sediment patches
