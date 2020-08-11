//update the motion of erodible patches by solving
//sediment transport on them.
if(calculateSed) //if sediment transport is turned on
{
    Info << "Solving sediment ... " << endl;

    //flatten the bed patch in the faMesh.
    //The Exner equation and sand slide solver assume the 
    //bed mesh is flat (2D). However, faMesh is a 3D surface.
    //Need to level the vertical coordinates of all points in the faMesh
    //only the localPoints are modified
//#   include "flatten_bed_faMesh.H"

    //update the point displacement on the erodible
    //patches (sediment bottom) by solving sediment
    //transport equations on them.
#   include "update_erodiblePatches.H"
}

