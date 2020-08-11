/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
  \\    /   O peration     |
  \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
  \\/     M anipulation  |
  -------------------------------------------------------------------------------
  License
  This file is part of OpenFOAM.

  OpenFOAM is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
  for more details.

  You should have received a copy of the GNU General Public License
  along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

  \*---------------------------------------------------------------------------*/



#include "areaFieldsTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    //write the faMesh to VTK
    void writeFaMesh2VTK
        (
         const faMesh& surf,
         //OFstream& os,
         const Time& runTime
        )
        {
            const auto& pp = surf.patch();

            fileName path(runTime.path()/"morphoSurface");
            mkDir(path);

            vtk::surfaceMeshWriter writer
                (
                 pp,
                 path,
                 Pstream::parRun()
                );

            writer.beginFile(surf.name());

            //writer.writeTimeValue(timeValue);
            writer.writeGeometry();

            fileName outputName(writer.output());

            writer.close();

            /*
               if (Pstream::master())
               {
            // Add to file-series and emit as JSON

            fileName seriesName(vtk::seriesWriter::base(outputName));

            vtk::seriesWriter& series = vtkSeries(seriesName);

            // First time?
            // Load from file, verify against filesystem,
            // prune time >= currentTime
            if (series.empty())
            {
            series.load(seriesName, true, timeValue);
            }

            series.append(timeValue, outputName);
            series.write(seriesName);
            }
            */
        }


    void areaScalarField2VTK
        (
         const faMesh& surf, 
         const areaScalarField& asf, 
         OFstream& os
        )
        {
        }

    void areaVectorField2VTK
        (
         const faMesh& surf,
         const areaVectorField& avf,
         OFstream& os
        )
        {
        }


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
