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

#include "flameletMakeChemistrySolverTypes.H"

#include "flameletThermoPhysicsTypes.H"
#include "flameletRhoReactionThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // Chemistry solvers based on sensibleEnthalpy
    flameletMakeChemistrySolverTypes(flameletRhoReactionThermo, flameletConstGasHThermoPhysics);
    flameletMakeChemistrySolverTypes(flameletRhoReactionThermo, flameletGasHThermoPhysics);
    flameletMakeChemistrySolverTypes
    (
        flameletRhoReactionThermo,
        flameletConstIncompressibleGasHThermoPhysics
    );
    flameletMakeChemistrySolverTypes
    (
        flameletRhoReactionThermo,
        flameletIncompressibleGasHThermoPhysics
    );
    flameletMakeChemistrySolverTypes(flameletRhoReactionThermo, flameletIcoPoly8HThermoPhysics);
    flameletMakeChemistrySolverTypes(flameletRhoReactionThermo, flameletConstFluidHThermoPhysics);
    flameletMakeChemistrySolverTypes
    (
        flameletRhoReactionThermo,
        flameletConstAdiabaticFluidHThermoPhysics
    );
    flameletMakeChemistrySolverTypes(flameletRhoReactionThermo, flameletConstHThermoPhysics);

    // Chemistry solvers based on sensibleInternalEnergy
    flameletMakeChemistrySolverTypes(flameletRhoReactionThermo, flameletConstGasEThermoPhysics);
    flameletMakeChemistrySolverTypes(flameletRhoReactionThermo, flameletGasEThermoPhysics);
    flameletMakeChemistrySolverTypes
    (
        flameletRhoReactionThermo,
        flameletConstIncompressibleGasEThermoPhysics
    );
    flameletMakeChemistrySolverTypes
    (
        flameletRhoReactionThermo,
        flameletIncompressibleGasEThermoPhysics
    );
    flameletMakeChemistrySolverTypes(flameletRhoReactionThermo, flameletIcoPoly8EThermoPhysics);
    flameletMakeChemistrySolverTypes(flameletRhoReactionThermo, flameletConstFluidEThermoPhysics);
    flameletMakeChemistrySolverTypes
    (
        flameletRhoReactionThermo,
        flameletConstAdiabaticFluidEThermoPhysics
    );
    flameletMakeChemistrySolverTypes(flameletRhoReactionThermo, flameletConstEThermoPhysics);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
