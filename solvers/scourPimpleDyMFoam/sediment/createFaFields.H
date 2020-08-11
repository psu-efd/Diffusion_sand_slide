Info<< "Reading field faZb" << endl;
areaScalarField faZb
(
 IOobject
 (
  "faZb",
  runTime.timeName(),
  mesh,
  IOobject::MUST_READ,
  IOobject::AUTO_WRITE
 ),
 aMesh
 );

//magnitude of the elevation gradient 
//at the edge of faMesh
edgeScalarField faeGradZb
(
 IOobject
 (
  "faeGradZb",
  runTime.timeName(),
  mesh,
  IOobject::MUST_READ,
  IOobject::AUTO_WRITE
 ),
 aMesh
 );

//sand slide diffusion coefficient
//at edge of faMesh
edgeScalarField faeSandSlideK
(
 IOobject
 (
  "faeSandSlideK",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
 ),
 aMesh,
 dimensionedScalar("0.0", K0.dimensions(), K0.value())
 );

//bedFaceUnitNormaly mapped to the faMesh
areaScalarField faBedFaceUnitNormaly
(
 IOobject
 (
  "faBedFaceUnitNormaly",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
 ),
 aMesh,
 dimensionedScalar("0.0", dimless, 0.0)
 );


//bedload transport rate vector at face center
areaVectorField faQb
(
 IOobject
 (
  "faQb",
  runTime.timeName(),
  mesh,
  IOobject::MUST_READ,
  IOobject::AUTO_WRITE
 ),
 aMesh
 );

areaVectorField faQb_copy
(
 IOobject
 (
  "faQb_copy",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
 ),
 faQb
 );


/*
edgeVectorField faEdgeQb
(
 IOobject
 (
  "faEdgeQb",
  runTime.timeName(),
  mesh,
  IOobject::MUST_READ,
  IOobject::AUTO_WRITE
 ),
 aMesh
 );
*/
//bed elevation change at face center
areaScalarField faDeltaZ
(
 IOobject
 (
  "faDeltaZ",
  runTime.timeName(),
  mesh,
  IOobject::MUST_READ,
  IOobject::AUTO_WRITE
 ),
 aMesh
 );

//bed elevation change due to sand slide
areaScalarField faDeltaZ_sand_slide
(
 IOobject
 (
  "faDeltaZ_sand_slide",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
 ),
 aMesh,
 dimensionedScalar("0.0", dimLength, 0.0)
 );

areaScalarField faDeltaZFlag_sand_slide
(
 IOobject
 (
  "faDeltaZFlag_sand_slide",
  runTime.timeName(),
  mesh,
  IOobject::READ_IF_PRESENT,
  IOobject::AUTO_WRITE
 ),
 aMesh,
 dimensionedScalar("0.0", dimless, 0.0)
 );



//bedload flux direction at edge center
edgeScalarField phi_qb_sign
(
 IOobject
 (
  "phi_qb_sign",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
 ),
 sign(fac::interpolate(faQb_copy) & aMesh.Le())
 );


//bedload flux at edge center
edgeScalarField phi_qb
(
 IOobject
 (
  "phi_qb",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
 ),
 fac::interpolate(faQb) & aMesh.Le()
 );

// Create volume-surface mapping object
volSurfaceMapping vsm(aMesh);

