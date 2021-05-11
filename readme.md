
 # README #
 
 This repository has the code and cases for the sand slide model. It is developed with
 OpenFOAM v1812.

 1. solvers: geomSlideFoam for geometric sand slide model 
             sandSlideFoam for slope-limited diffusive sand slide model
             scourPimpleDyMFoam for scour model using slope-limited diffusive
             sand slide method.
 2. libraries: library for sediment transport.
 3. utilities: some tools, for example the setup of initial bed elevation 
               field using the setBedShape tool.
 4. cases: all simulation cases.

## Reference ##

Y. Song, Y. Xu, and X. Liu (2020). A gradient-limited diffusive sand slide method in scour models. Journal of Hydraulic Engineering. 146(11): 04020074

## Disclaimer ##
This offering is not approved or endorsed by OpenCFD Limited, producer and distributor of the OpenFOAM software via <www.openfoam.com>, and owner of the OPENFOAM&reg;  and OpenCFD&reg; trade marks.

OPENFOAM&reg; is a registered trade mark of OpenCFD Limited, producer and distributor of the OpenFOAM software via <www.openfoam.com>.


## Authors and contributors: ##

 Xiaofeng Liu and Yalan Song
 Penn State University

 
 
