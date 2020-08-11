{
    surfaceScalarField phiC
    (
        IOobject
        (
            "phiC",
            runTime.timeName(),
            mesh
        ),
        phi + (mesh.Sf() & fvc::interpolate(Ws))
    );

    solve
    (
          fvm::ddt(Cs)
        + fvm::div(phiC, Cs)
        - fvm::laplacian(turbulence->nuEff(), Cs)
    );

    Info<< "Suspended sediment volume fraction = "
        << Cs.weightedAverage(mesh.V()).value()
        << "  Min(Cs) = " << min(Cs).value()
        << "  Max(Cs) = " << max(Cs).value()
        << endl;

    //Clipping of Cs in case it goes wild.
    //Not a good idea. The numerical algorithm and mesh
    //should guarantee the boundedness.
//    Cs.min(1.0);
//    Cs.max(0.0);
}
