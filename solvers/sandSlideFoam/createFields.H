Info<< "Reading elevation field eta\n" << endl;

//bed elevation
volScalarField eta
(
    IOobject
    (
        "eta",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

//static friction angle
dimensionedScalar angle0
(
    transportProperties.lookup("angle0")
);

scalar angle=degToRad(angle0.value());

Info<< "Reading diffusivity K\n" << endl;

dimensionedScalar K0
(
    transportProperties.lookup("K0")
);

word gradientScheme
(
    transportProperties.lookup("gradientScheme")
);


surfaceScalarField K
(
    IOobject
    (
        "K",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("0.0", K0.dimensions(), 0)
);

