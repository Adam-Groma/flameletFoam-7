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

#include "flameletSpecie.H"
#include "flameletPerfectGas.H"
#include "flameletIncompressiblePerfectGas.H"
#include "flameletHConstThermo.H"
#include "flameletJanafThermo.H"
#include "flameletSensibleEnthalpy.H"
#include "flameletThermo.H"
#include "flameletRhoConst.H"
#include "flameletPerfectFluid.H"
#include "flameletAdiabaticPerfectFluid.H"

#include "flameletConstTransport.H"
#include "flameletSutherlandTransport.H"
#include "flameletWLFTransport.H"

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
    flameletConstGasHThermoPhysics
);

flameletMakeThermoPhysicsReactionThermos
(
    flameletRhoThermo,
    flameletRhoReactionThermo,
    flameletHeRhoThermo,
    flameletReactingMixture,
    flameletGasHThermoPhysics
);

flameletMakeThermoPhysicsReactionThermos
(
    flameletRhoThermo,
    flameletRhoReactionThermo,
    flameletHeRhoThermo,
    flameletReactingMixture,
    flameletConstIncompressibleGasHThermoPhysics
);

flameletMakeThermoPhysicsReactionThermos
(
    flameletRhoThermo,
    flameletRhoReactionThermo,
    flameletHeRhoThermo,
    flameletReactingMixture,
    flameletIncompressibleGasHThermoPhysics
);

flameletMakeThermoPhysicsReactionThermos
(
    flameletRhoThermo,
    flameletRhoReactionThermo,
    flameletHeRhoThermo,
    flameletReactingMixture,
    flameletIcoPoly8HThermoPhysics
);

flameletMakeThermoPhysicsReactionThermos
(
    flameletRhoThermo,
    flameletRhoReactionThermo,
    flameletHeRhoThermo,
    flameletReactingMixture,
    flameletConstFluidHThermoPhysics
);

flameletMakeThermoPhysicsReactionThermos
(
    flameletRhoThermo,
    flameletRhoReactionThermo,
    flameletHeRhoThermo,
    flameletReactingMixture,
    flameletConstAdiabaticFluidHThermoPhysics
);

flameletMakeThermoPhysicsReactionThermos
(
    flameletRhoThermo,
    flameletRhoReactionThermo,
    flameletHeRhoThermo,
    flameletReactingMixture,
    flameletConstHThermoPhysics
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
