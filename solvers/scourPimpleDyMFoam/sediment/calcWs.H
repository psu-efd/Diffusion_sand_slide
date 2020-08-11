//Sub routine to calculate particle fall velocity

dimensionedScalar Rep = Foam::sqrt((rhos/rhow-1)*mag(g)*diam)*diam/nu;
dimensionedScalar x = Foam::log10(sqr(Rep));
dimensionedScalar y = -3.76715 + 1.92944*x -0.09815*x*x -0.00575*x*x*x + 0.00056*x*x*x*x;
dimensionedScalar Rf = Foam::pow(pow(10,y)/Rep, 1.0/3.0);

Ws = Rf*Foam::sqrt((rhos/rhow-1)*mag(g)*diam)*g/mag(g);
forAll(Ws.boundaryField(), patchI)
{
    Ws.boundaryFieldRef()[patchI] == 
       (Rf*Foam::sqrt((rhos/rhow-1)*mag(g)*diam)*g/mag(g)).value();
}

//Info << Ws << endl;

