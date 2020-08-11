    Info << "Flattening the faMesh ..." << endl;

    pointField& pF = const_cast<pointField&>(aMesh.patch().localPoints());

    //make the patch flat (2D) in y direction
    forAll(aMesh.patch().localPoints(), pointI)
    {
        pF[pointI].y() = 0.0;
    }

