//Read sediment transport properties
Info<< "Reading sediment transport properties\n" << endl;

IOdictionary sedimentProperties
(
 IOobject
 (
  "sedimentProperties",
  runTime.constant(),
  mesh,
  IOobject::MUST_READ,
  IOobject::NO_WRITE
 )
 );

//a switch for whether solve sediment transport
Switch calculateSed(sedimentProperties.lookup("calculateSed"));
if(calculateSed)
{
    Info << "Sediment transport is turned on." << endl;
}
else
{
    Info << "Sediment transort is turned off." << endl;
}

//if sediment transport is on, whether to check mesh every time
//the mesh is upaded: default is false.
Switch checkMesh(sedimentProperties.lookupOrDefault("checkMesh", false));
if(checkMesh)
{
    Info << "checkMesh is turned on." << endl;
}
else
{
    Info << "checkMesh is turned off." << endl;
}

//a switch for whether solve bedload
Switch calcBedLoad(sedimentProperties.lookup("calcBedLoad"));
if(calculateSed)
{
    if(calcBedLoad)
    {
        Info << "Bedload is turned on." << endl;
    }
    else
    {
        Info << "Bedload is turned off." << endl;
    }
}

//a switch for whether solve suspended load 
Switch calcSuspLoad(sedimentProperties.lookup("calcSuspLoad"));
if(calculateSed)
{
    if(calcSuspLoad)
    {
        Info << "Suspended load is turned on." << endl;
    }
    else
    {
        Info << "Suspended load is turned off." << endl;
    }
}

//Give a warning if sediment transport is turned on but
//neither bedload nor suspended load is on.
if(calculateSed && (!calcSuspLoad) && (!calcBedLoad))
{
    WarningIn
        (
         "createSedimentFields.H"
        ) << endl << "Sediment transport is turned on but " 
        << "neither bedload nor suspended load is on."<< endl
        << endl;
}

//parameters for suspended sediment transport
word suspSedBC(sedimentProperties.lookup("suspSedBC"));
Info << "Choosing suspended sediment B.C. option: " << suspSedBC << endl;

//only relevant if "suspSedBC" = "ExtrapolateToWall"
word extrapolationOption(sedimentProperties.lookup("extrapolationOption"));

//end of parameters for suspended sediment

dimensionedScalar nu
(
 sedimentProperties.lookup("nu")
 );

dimensionedScalar rhos
(
 sedimentProperties.lookup("rhos")
 );

dimensionedScalar rhow
(
 sedimentProperties.lookup("rhow")
 );
dimensionedScalar diam
(
 sedimentProperties.lookup("diam")
 );
dimensionedScalar porosity
(
 sedimentProperties.lookup("porosity")
 );
dimensionedScalar theta_cri_b
(
 sedimentProperties.lookup("theta_cri_b")
 );
dimensionedScalar theta_cri_s
(
 sedimentProperties.lookup("theta_cri_s")
 );

dimensionedScalar angleOfRepose
(
 sedimentProperties.lookup("angleOfRepose")
 );
dimensionedScalar C
(
 sedimentProperties.lookup("C")
 );
dimensionedScalar sigmac
(
 sedimentProperties.lookup("sigmac")
 );

//static friction coefficient
dimensionedScalar mus
(
 sedimentProperties.lookup("mus")
);

label bedInterval(readInt(sedimentProperties.lookup("nFluidStepsPerBedStep")));

//read parameters for sand slide model
const dictionary& sandSlideDict = sedimentProperties.subDict("sandSlide");

//whether turn on the sand slide model: default is false.
Switch sandSlide(sandSlideDict.lookupOrDefault("sandSlide", false));
if(sandSlide)
{
    Info << "Sand slide model is turned on." << endl;
}
else
{
    Info << "Sand slide model is turned off." << endl;
}


//maximum iterations for the sand slide step
label maxSandSlideIter(readInt(sandSlideDict.lookup("maxSandSlideIter")));

//Option for the calculation of bed angle at faces in fvMesh (edges in faMesh)
//1. linearCD: linear central differencing
//2. nonlinearGauss:
word bedAngleAtSurfaceOption(sandSlideDict.lookup("bedAngleAtSurfaceOption"));
if(bedAngleAtSurfaceOption.match("LGS") || bedAngleAtSurfaceOption.match("GNGS"))
{
    Info << "Choosing the option for the calculation of bed angle at faces: "
         << bedAngleAtSurfaceOption << endl;
}
else
{
   FatalErrorIn("bedAngleAtSurfaceOption")
              << "The specified bedAngleAtSurfaceOption: " << bedAngleAtSurfaceOption
              << " is not valid." << nl
              << "Valid options are: LGS and GNGS"
              << abort(FatalError);
}

Info<< "Reading sand slide diffusivity K0\n" << endl;
dimensionedScalar K0
(
    sandSlideDict.lookup("K0")
);

//sand slide diffusivity defined on fvMesh
volScalarField sandSlideK
(
    IOobject
    (
        "sandSlideK",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("0.0", K0.dimensions(), K0.value())
);

//Create the list of sediment boundary patches.
//There might be more than one sediment patches due to various reasons,
//such as meshing requirement, different sediment sizes, etc.
IOdictionary sedimentPatchesDict
(
 IOobject
 (
  "sedimentPatchesDict",
  runTime.constant(),
  mesh,
  IOobject::MUST_READ,
  IOobject::NO_WRITE
 )
 );

wordList sedimentPatches(sedimentPatchesDict.lookup("sedimentPatches"));


Info<< "Initialising field Ws\n" << endl;
volVectorField Ws
(
 IOobject
 (
  "Ws",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
 ),
 mesh,
 dimensionedVector("Ws", U.dimensions(), vector::zero)
 );

Info<< "Reading sediment concentration field Cs\n" << endl;
volScalarField Cs
(
 IOobject
 (
  "Cs",
  runTime.timeName(),
  mesh,
  IOobject::MUST_READ,
  IOobject::AUTO_WRITE
 ),
 mesh
 );

//surface normal gradient of Cs
volScalarField wallGradCs
(
 IOobject
 (
  "wallGradCs",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
 ),
 mesh,
 dimensionedScalar
 (
  "wallGradCs",
  Cs.dimensions()/dimLength,
  0.0
 )
 );


//deposition rate
volScalarField D
(
 IOobject
 (
  "D",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
 ),
 mesh,
 dimensionedScalar
 (
  "D",
  U.dimensions(),
  0.0
 )
 );

//erosion rate
volScalarField E
(
 IOobject
 (
  "E",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
 ),
 mesh,
 dimensionedScalar
 (
  "E",
  U.dimensions(),
  0.0
 )
 );

//The vertical postion change of the bed cells(it equals the bottom elevation
//of fluid mesh. It's important to let the initial values be equal)
volScalarField bedCellHeightChange
(
 IOobject
 (
  "bedCellHeightChange",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
 ),
 mesh,
 dimensionedScalar("bedCellHeightChange", dimLength, 0.0)
 );

//bed change due to sand slide
volScalarField bedChange_sand_slide
(
 IOobject
 (
 "bedChange_sand_slide",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
  ),
 mesh,
 dimensionedScalar("bedChange_sand_slide", dimLength, 0.0)
 );

 //bed change flag due to sand slide
volScalarField bedChangeFlag_sand_slide
(
 IOobject
 (
 "bedChangeFlag_sand_slide",
  runTime.timeName(),
  mesh,
  IOobject::READ_IF_PRESENT,
  IOobject::AUTO_WRITE
 ),
 mesh,
 dimensionedScalar("bedChangeFlag_sand_slide", dimless, 0.0)
 );



//bedload transport rate vector
volVectorField bedWallqi
(
 IOobject
 (
  "bedWallqi",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
 ),
 mesh,
 dimensionedVector("bedWallqi", dimensionSet(0,2,-1,0,0,0,0), vector(0,0,0))
 );


//flow direction: >0 uphill
//                <0 downhill
volScalarField bedWallFlowDir
(
 IOobject
 (
  "bedWallFlowDir",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
 ),
 mesh,
 dimensionedScalar("bedWallFlowDir", dimless, 0.0)
 );


//critical Shields number field
volScalarField bedWallCriticalShieldsNum
(
 IOobject
 (
  "bedWallCriticalShieldsNum",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
 ),
 mesh,
 dimensionedScalar("bedWallCriticalShieldsNum", dimless , theta_cri_b.value())
 );


volVectorField wallShearStress
(
 IOobject
 (
  "wallShearStress",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
 ),
 mesh,
 dimensionedVector("wallShearStress", dimensionSet(0,2,-2,0,0,0,0),
     vector::zero)
 );

//small number for wall shear magnitude
dimensionedScalar smallWallShearMag
(  
 "smallWallShearMag",
 dimensionSet(0,2,-2,0,0,0,0),
 SMALL
 );

volScalarField wallShearStressMag
(
 IOobject
 (
  "wallShearStressMag",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
 ),
 mag(wallShearStress)+smallWallShearMag
 );

volScalarField wallShieldsNumber
(
 IOobject
 (
  "wallShieldsNumber",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
 ),
 mesh,
 dimensionedScalar("wallShieldsNumber", dimless, 0.0)
 );

//Reference concentration at a reference height
//2*diameter from the bed for example in Zyserman and Fredsoe (1994) formula)
volScalarField Cb_star
(
 IOobject
 (
  "Cb_star",
  runTime.timeName(),
  mesh,
  IOobject::MUST_READ,
  IOobject::AUTO_WRITE
 ),
 mesh
 );

//Equilibrium reference concentration Cw_star right on the bed.
volScalarField Cw_star
(
 IOobject
 (
  "Cw_star",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
 ),
 mesh,
 dimensionedScalar("Cw_star", dimless, 0.0)
 );

//Sediment flux through boundary patches. Only meaningful
//for the bed boundaries.
volScalarField suspSedFlux
(
 IOobject
 (
  "suspSedFlux",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
 ),
 mesh,
 dimensionedScalar("suspSedFlux", dimArea, 0.0)
 );


volVectorField bedFaceUnitNormal
(
 IOobject
 (
  "bedFaceUnitNormal",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
 ),
 mesh,
 dimensionedVector("0.0", dimless, Zero)
 );

volScalarField bedFaceUnitNormaly
(
 IOobject
 (
  "bedFaceUnitNormaly",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
 ),
 mesh,
 dimensionedScalar("0.0", dimless, 0.0)
 );


//faces surface steepest slope unit vector(for the slope effects)
volVectorField bedSteepestSlopeVector
(
 IOobject
 (
  "bedSteepestSlopeVector",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
 ),
 mesh,
 vector(SMALL,SMALL,SMALL)
 );

//wall surface slope angle (with respect to horizontal surface)
volScalarField bedWallSlopeAngle
(
 IOobject
 (
  "bedWallSlopeAngle",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
 ),
 mesh,
 dimensionedScalar("bedWallSlopeAngle", dimless, 0.0)
);

//angle between wallshear (velocity) and bedSteepestSlopeVector
volScalarField bedWallPhiAngle
(
 IOobject
 (
  "bedWallPhiAngle",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
 ),
 mesh,
 dimensionedScalar("bedWallPhiAngle", dimless, 0.0)
);


//bed load sediment transport load
//q0: bedload sediment transport rate per unit width (m2/s) at a horizontal bed
volScalarField bedWallq0
(
 IOobject
 (
  "bedWallq0",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
 ),
 mesh,
 dimensionedScalar("bedWallq0", dimensionSet(0,2,-1,0,0,0,0), 0.0)
 );

volVectorField bedSlope
(
 IOobject
 (
  "bedSlope",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
 ),
 mesh,
 dimensionedVector("bedSlope", dimless, Zero)
 );

//face centers heigh
volScalarField fluidFaceCenterHeigh
(
 IOobject
 (
  "fluidFaceCenterHeigh",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
 ),
 mesh,
 dimensionedScalar("fluidFaceCenterHeigh", dimLength, 0.0)
 );

//accumulated bed elevation change at face centers
volScalarField fluidFaceCenterTotalDeltaY
(
 IOobject
 (
  "fluidFaceCenterTotalDeltaY",
  runTime.timeName(),
  mesh,
  IOobject::READ_IF_PRESENT,
  IOobject::AUTO_WRITE
 ),
 mesh,
 dimensionedScalar("fluidFaceCenterTotalDeltaY", dimLength, 0.0)
 );


//bed face project area on x-z plane
volScalarField bedFace2DArea
(
 IOobject
 (
  "bedFace2DArea",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
 ),
 mesh,
 dimensionedScalar("bedFace2DArea", dimArea, 0.0)
);

#include "calcWs.H"



