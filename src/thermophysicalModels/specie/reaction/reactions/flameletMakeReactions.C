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

#include "flameletReactionTypes.H"
#include "flameletMakeReaction.H"

#include "flameletArrheniusReactionRate.H"
#include "flameletInfiniteReactionRate.H"
#include "flameletLandauTellerReactionRate.H"
#include "flameletThirdBodyArrheniusReactionRate.H"

#include "flameletChemicallyActivatedReactionRate.H"
#include "flameletJanevReactionRate.H"
#include "flameletPowerSeriesReactionRate.H"

#include "flameletFallOffReactionRate.H"
#include "flameletLindemannFallOffFunction.H"
#include "flameletSRIFallOffFunction.H"
#include "flameletTroeFallOffFunction.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define flameletMakeReactions(Thermo, Reaction)                                        \
                                                                               \
    defineTemplateTypeNameAndDebug(Reaction, 0);                               \
    defineTemplateRunTimeSelectionTable(Reaction, dictionary);                 \
                                                                               \
    flameletMakeIRNReactions(Thermo, flameletArrheniusReactionRate)                            \
    flameletMakeIRNReactions(Thermo, flameletInfiniteReactionRate)                             \
    flameletMakeIRNReactions(Thermo, flameletLandauTellerReactionRate)                         \
    flameletMakeIRNReactions(Thermo, flameletThirdBodyArrheniusReactionRate)                   \
                                                                               \
    flameletMakeIRReactions(Thermo, flameletJanevReactionRate)                                 \
    flameletMakeIRReactions(Thermo, flameletPowerSeriesReactionRate)                           \
                                                                               \
    flameletMakePressureDependentReactions                                             \
    (                                                                          \
       Thermo,                                                                 \
       flameletArrheniusReactionRate,                                                  \
       flameletLindemannFallOffFunction                                                \
    )                                                                          \
                                                                               \
    flameletMakePressureDependentReactions                                             \
    (                                                                          \
       Thermo,                                                                 \
       flameletArrheniusReactionRate,                                                  \
       flameletTroeFallOffFunction                                                     \
    )                                                                          \
                                                                               \
    flameletMakePressureDependentReactions                                             \
    (                                                                          \
       Thermo,                                                                 \
       flameletArrheniusReactionRate,                                                  \
       flameletSRIFallOffFunction                                                      \
    )


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // sensible enthalpy based reactions
    flameletMakeReactions(flameletConstGasHThermoPhysics, flameletConstGasHReaction)
    flameletMakeReactions(flameletGasHThermoPhysics, flameletGasHReaction)
    flameletMakeReactions
    (
        flameletConstIncompressibleGasHThermoPhysics,
        flameletConstIncompressibleGasHReaction
    )
    flameletMakeReactions(flameletIncompressibleGasHThermoPhysics, flameletIncompressibleGasHReaction)
    flameletMakeReactions(flameletIcoPoly8HThermoPhysics, flameletIcoPoly8HReaction)
    flameletMakeReactions(flameletConstFluidHThermoPhysics, flameletConstFluidHReaction)
    flameletMakeReactions
    (
        flameletConstAdiabaticFluidHThermoPhysics,
        flameletConstAdiabaticFluidHReaction
    )
    flameletMakeReactions(flameletConstHThermoPhysics, flameletConstHReaction)

    flameletMakeReactions(flameletConstGasEThermoPhysics, flameletConstGasEReaction)
    flameletMakeReactions(flameletGasEThermoPhysics, flameletGasEReaction)
    flameletMakeReactions
    (
        flameletConstIncompressibleGasEThermoPhysics,
        flameletConstIncompressibleGasEReaction
    )
    flameletMakeReactions(flameletIncompressibleGasEThermoPhysics, flameletIncompressibleGasEReaction)
    flameletMakeReactions(flameletIcoPoly8EThermoPhysics, flameletIcoPoly8EReaction)
    flameletMakeReactions(flameletConstFluidEThermoPhysics, flameletConstFluidEReaction)
    flameletMakeReactions
    (
        flameletConstAdiabaticFluidEThermoPhysics,
        flameletConstAdiabaticFluidEReaction
    )
    flameletMakeReactions(flameletConstEThermoPhysics, flameletConstEReaction)
}

// ************************************************************************* //
