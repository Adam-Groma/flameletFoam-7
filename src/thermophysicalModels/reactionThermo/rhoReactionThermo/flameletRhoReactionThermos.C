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

#include "makeReactionThermo.H"

#include "flameletRhoReactionThermo.H"
#include "heRhoThermo.H"

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

#include "constTransport.H"
#include "sutherlandTransport.H"
#include "WLFTransport.H"

#include "flameletMultiComponentMixture.H"
#include "flameletReactingMixture.H"

#include "thermoPhysicsTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Reaction thermo for sensible enthalpy

makeThermoPhysicsReactionThermos
(
    rhoThermo,
    flameletRhoReactionThermo,
    heRhoThermo,
    flameletReactingMixture,
    constGasHThermoPhysics
);

makeThermoPhysicsReactionThermos
(
    rhoThermo,
    flameletRhoReactionThermo,
    heRhoThermo,
    flameletReactingMixture,
    gasHThermoPhysics
);

makeThermoPhysicsReactionThermos
(
    rhoThermo,
    flameletRhoReactionThermo,
    heRhoThermo,
    flameletReactingMixture,
    constIncompressibleGasHThermoPhysics
);

makeThermoPhysicsReactionThermos
(
    rhoThermo,
    flameletRhoReactionThermo,
    heRhoThermo,
    flameletReactingMixture,
    incompressibleGasHThermoPhysics
);

makeThermoPhysicsReactionThermos
(
    rhoThermo,
    flameletRhoReactionThermo,
    heRhoThermo,
    flameletReactingMixture,
    icoPoly8HThermoPhysics
);

makeThermoPhysicsReactionThermos
(
    rhoThermo,
    flameletRhoReactionThermo,
    heRhoThermo,
    flameletReactingMixture,
    constFluidHThermoPhysics
);

makeThermoPhysicsReactionThermos
(
    rhoThermo,
    flameletRhoReactionThermo,
    heRhoThermo,
    flameletReactingMixture,
    constAdiabaticFluidHThermoPhysics
);

makeThermoPhysicsReactionThermos
(
    rhoThermo,
    flameletRhoReactionThermo,
    heRhoThermo,
    flameletReactingMixture,
    constHThermoPhysics
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
