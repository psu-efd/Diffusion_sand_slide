Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, pimple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());


singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

#include "createMRF.H"
#include "createFvOptions.H"

#include "readGravitationalAcceleration.H"
//A not-so-smart trick to check whether the vertical 
//direction is in y. It checks whether the y-component of
//gravity is the largest. 
if((Foam::mag(g.value().y())<Foam::mag(g.value().x())) ||
   (Foam::mag(g.value().y())<Foam::mag(g.value().z())) )
{
    FatalErrorIn("createFields.H")
        << "Judging by the gravity vector, the vertical direction is not in y." << nl
        << abort(FatalError);
}

