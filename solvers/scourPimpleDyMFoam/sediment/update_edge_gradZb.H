
//The following two options take different path to get the edge normal 
//gradient of bed elevation:
//1. LGS: linearCD, directly uses the snGrad function to calculate the 
//   surface normal gradient with the values of cells sharing a face. The 
//   name "LGS" is slightly misleading because one can specify any snGrad
//   scheme available in OpenFOAM.
//2. GNGS: nonlinearGauss, first call the grad function to get the bed elevation 
//   at the cell centers (which uses the Gauss theorem and interntally
//   calls the interpolation of bed elevation from cell centers to face centers).
//   Then, the bed elevation gradient is interpolated from cell centers to face centers,
//   and dot product with the surface normal. In this way, the bed angle at the faces contains
//   elevation information from much bigger cell stencil. The user has the option to specify 
//   different numerical schemes for different differntial operators and interpolations 
//   along the way.
 bool checkSandSlide = false;
 scalar slopeCounter=0;

if(bedAngleAtSurfaceOption.match("LGS"))
{
    //calculate the line normal gradient of the bed elevation
    
   faeGradZb = mag(fac::lnGrad(faZb, "lnGrad(faZb)"));



    forAll (faeGradZb, I)
    {
         if( (mag(faeGradZb[I]) - Foam::tan(angleOfRepose.value()))>0.01 )
         {
             slopeCounter++;
         }
            
    }
   
    if(returnReduce(slopeCounter, maxOp<label>())>0) checkSandSlide=true;

}
else if(bedAngleAtSurfaceOption.match("GNGS"))
{
    //calculate the bed elevation gradient at cell centers
    areaVectorField faGradZb = fac::grad(faZb);


    //interpolate the bed elevation gradient from cell centers to face centers
    faeGradZb = mag(fac::interpolate(faGradZb));
    

    //back to linear gradient to eliminate residuals
    //For now, use a maximum of 50 iterations of GNGS before
    //switch to LGS to clean up the residuals. 
    if (sandSlideIter>50){ faeGradZb = mag(fac::lnGrad(faZb, "lnGrad(faZb)"));}


    forAll(faeGradZb, I)
    {

         if( (mag(faeGradZb[I]) - Foam::tan(angleOfRepose.value()))>0.01 )
         {
             slopeCounter++;
         }
   }
    if(returnReduce(slopeCounter, maxOp<label>())>0) checkSandSlide=true;

}
else
{
   FatalErrorIn("bedAngleAtSurfaceOption")
              << "The specified bedAngleAtSurfaceOption is not valid." << nl 
              << abort(FatalError);
}

Info << "   Maximum bed elevation gradient magnitude = " << max(mag(faeGradZb)).value() << endl;
Info << "   Tan(angle of repose) = " << Foam::tan(angleOfRepose.value()) << endl;
Info << "   Number of large angle = " << returnReduce(slopeCounter, maxOp<label>()) <<endl;
