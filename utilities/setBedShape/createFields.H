    Info<< "Reading elevation field eta\n" << endl;

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
