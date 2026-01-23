/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2018 OpenFOAM Foundation
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

#include "flameletRhoReactionThermo.H"
#include "flameletHeRhoThermo.H"

#include "specie.H"
#include "perfectGas.H"
#include "incompressiblePerfectGas.H"
#include "hConstThermo.H"
#include "janafThermo.H"
#include "sensibleEnthalpy.H"
#include "thermo.H"
#include "rhoConst.H"
#include "perfectFluid.H"
#include "adiabaticPerfectFluid.H"

#include "flameletConstTransport.H"
#include "sutherlandTransport.H"
#include "WLFTransport.H"

#include "flameletMultiComponentMixture.H"
#include "flameletReactingMixture.H"

#include "flameletThermoPhysicsTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Reaction thermo for sensible enthalpy

flameletMakeThermoPhysicsReactionThermos
(
    flameletRhoThermo,
    flameletRhoReactionThermo,
    flameletHeRhoThermo,
    flameletReactingMixture,
    constGasHThermoPhysics
);

flameletMakeThermoPhysicsReactionThermos
(
    flameletRhoThermo,
    flameletRhoReactionThermo,
    flameletHeRhoThermo,
    flameletReactingMixture,
    gasHThermoPhysics
);

flameletMakeThermoPhysicsReactionThermos
(
    flameletRhoThermo,
    flameletRhoReactionThermo,
    flameletHeRhoThermo,
    flameletReactingMixture,
    constIncompressibleGasHThermoPhysics
);

flameletMakeThermoPhysicsReactionThermos
(
    flameletRhoThermo,
    flameletRhoReactionThermo,
    flameletHeRhoThermo,
    flameletReactingMixture,
    incompressibleGasHThermoPhysics
);

flameletMakeThermoPhysicsReactionThermos
(
    flameletRhoThermo,
    flameletRhoReactionThermo,
    flameletHeRhoThermo,
    flameletReactingMixture,
    icoPoly8HThermoPhysics
);

flameletMakeThermoPhysicsReactionThermos
(
    flameletRhoThermo,
    flameletRhoReactionThermo,
    flameletHeRhoThermo,
    flameletReactingMixture,
    constFluidHThermoPhysics
);

flameletMakeThermoPhysicsReactionThermos
(
    flameletRhoThermo,
    flameletRhoReactionThermo,
    flameletHeRhoThermo,
    flameletReactingMixture,
    constAdiabaticFluidHThermoPhysics
);

flameletMakeThermoPhysicsReactionThermos
(
    flameletRhoThermo,
    flameletRhoReactionThermo,
    flameletHeRhoThermo,
    flameletReactingMixture,
    constHThermoPhysics
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
