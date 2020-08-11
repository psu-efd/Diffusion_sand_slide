//Update the motion of erodible beds by solve sediment transport
//The displacement of points is from the solution of Exner equation
//on the faMesh.        

//Morphological change step is solved every 'bedInterval' fluid steps.
//The bedCount variable is a counter to keep record of when the morphological
//step should be performed.

  //bedCount++;
  //bedCount = bedCount % bedInterval;
  //Info << "\nbedCount = " << bedCount << nl << endl;

  //if(bedCount==0 && runTime.value()!=runTime.startTime().value())


//Morphological factor is bedInterval, i.e., fluid is assumed not changing (frozen)
//for every bedInterval morphological steps
  {
    Info << "\nSolving sediment transport and moving mesh ..." << nl <<  endl;

    //We have following steps:
    //1. update wall shear stress and Shields parameter
    //2. calcuate suspended load
    //3. calculate bed load
    //4. solve Exner equation and update mesh motion
    //5. sand slide model to correct bed angles

    //1. update wall shear stress and Shields parameters
    #include "update_wall_shear.H"

    //2. suspended load
    if(calcSuspLoad) 
    {
       #include "suspended_load.H"
    }

    //3. bed load
    if(calcBedLoad) 
    {
       #include "bed_load.H"
    }

    //4. solve the Exner equation
    #include "solve_Exner.H"

    //5. sand slide
    #include "sand_slide.H"
  }
