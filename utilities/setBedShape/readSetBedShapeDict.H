Info<< "Reading setBedShapeDict\n" << endl;

IOdictionary setBedShapeDict
    (
        IOobject
        (
            "setBedShapeDict",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

//define a sand column
//initial diameter
dimensionedScalar ri
(
    setBedShapeDict.lookup("ri")
);

//initial height
dimensionedScalar hi
(
    setBedShapeDict.lookup("hi")
);

scalar R0=ri.value();


