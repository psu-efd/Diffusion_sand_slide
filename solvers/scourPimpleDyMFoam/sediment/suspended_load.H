//calculate the suspended load

//update suspended sediment BC
if(suspSedBC.match("vanRijn1984PickupFunction"))
{
   #include "updateSedimentBC_vanRijn1984.H"
}
else if(suspSedBC.match("IgnoreRefHeight"))
{
   #include "updateSedimentBC_IgnoreRefHeight.H"
}
else if(suspSedBC.match("ExtrapolateToWall"))
{
   #include "updateSedimentBC_ExtrapolateToWall.H"
}
else
{
   FatalErrorIn("suspSedBC options")
              << "The option for suspSedBC: " << suspSedBC
              << " is not valid." << nl
              << "Valid options are: vanRijn1984PickupFunction, "
              << "IgnoreRefHeight, "
              << "ExtrapolateToWall."
              << abort(FatalError);
}

//solve the suspended sediment transport equation
#include "CsEqn.H"

