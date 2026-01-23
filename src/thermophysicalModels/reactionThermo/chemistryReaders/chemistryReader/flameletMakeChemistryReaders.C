/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "flameletMakeReactionThermo.H"
#include "flameletThermoPhysicsTypes.H"

#include "flameletChemistryReader.H"
#include "flameletFoamChemistryReader.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Solid chemistry readers based on sensibleEnthalpy

flameletMakeChemistryReader(constGasHThermoPhysics);
flameletMakeChemistryReader(gasHThermoPhysics);
flameletMakeChemistryReader(constIncompressibleGasHThermoPhysics);
flameletMakeChemistryReader(incompressibleGasHThermoPhysics);
flameletMakeChemistryReader(icoPoly8HThermoPhysics);
flameletMakeChemistryReader(constFluidHThermoPhysics);
flameletMakeChemistryReader(constAdiabaticFluidHThermoPhysics);
flameletMakeChemistryReader(constHThermoPhysics);

flameletMakeChemistryReaderType(flameletFoamChemistryReader, constGasHThermoPhysics);
flameletMakeChemistryReaderType(flameletFoamChemistryReader, gasHThermoPhysics);
flameletMakeChemistryReaderType
(
    flameletFoamChemistryReader,
    constIncompressibleGasHThermoPhysics
);
flameletMakeChemistryReaderType(flameletFoamChemistryReader, incompressibleGasHThermoPhysics);
flameletMakeChemistryReaderType(flameletFoamChemistryReader, icoPoly8HThermoPhysics);
flameletMakeChemistryReaderType(flameletFoamChemistryReader, constFluidHThermoPhysics);
flameletMakeChemistryReaderType(flameletFoamChemistryReader, constAdiabaticFluidHThermoPhysics);
flameletMakeChemistryReaderType(flameletFoamChemistryReader, constHThermoPhysics);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
