// solve the Exner equation and update the point displacement
// on the erodible patches.

//bed elevation change at face center
//the orginal fac::div(...) is designed for avalance to consider coordiniate transformation
//which needs to transform the result with (I-n*n) where n is the surface normal vector. This
//is not necessary for our purpose. So we can't use that.
//faDeltaZ = runTime.deltaT()*bedInterval/(1.0-porosity.value())*(-fac::div(faQb)); 
//faDeltaZ = runTime.deltaT()*bedInterval/(1.0-porosity.value())*(-fac::edgeIntegrate(phi_qb))
//            *faBedFaceUnitNormaly; 
faDeltaZ = runTime.deltaT()*bedInterval/(1.0-porosity.value())*(-fac::edgeIntegrate(phi_qb));


Info << "Maximum and minimum bed elevation change rate (m/hr) = " 
    << (Foam::max(faDeltaZ.internalField())/runTime.deltaT()).value()*3600.0 << ", "
    << (Foam::min(faDeltaZ.internalField())/runTime.deltaT()).value()*3600.0
    << endl;

    /*
    //manually assign bed change for testing purpose
    areaVectorField c = aMesh.areaCentres();
    forAll(c.internalField(), i)
    {
    vector curCenter = c[i];
    faDeltaZ[i] = 0.001*Foam::sin(curCenter.x()/2.0*2.0*3.1415926);
    }
    //end of manual bed change
    */

    faDeltaZ.correctBoundaryConditions();

    //original bed elevation at face center before Exner solution and sand slide
    areaScalarField faZb_org = faZb;

    //update bed elevation after the Exner equation solver step
    //However, the solution of Exner equation does not guarantee
    //the angle of repose, whicn needs to be taken care of in the
    //sand slide step.
    faZb += faDeltaZ;
    faZb.correctBoundaryConditions();

