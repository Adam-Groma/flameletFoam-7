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

Typedefs
    Foam::flameletReactionTypes

Description
    Type definitions for reactions

\*---------------------------------------------------------------------------*/

#ifndef flameletReactionTypes_H
#define flameletReactionTypes_H

#include "flameletThermoPhysicsTypes.H"
#include "FlameletReaction.H"

#include "flameletIcoPolynomial.H"
#include "flameletHPolynomialThermo.H"
#include "flameletPolynomialTransport.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // sensible enthalpy based reactions
    typedef FlameletReaction<flameletConstGasHThermoPhysics> flameletConstGasHReaction;

    typedef FlameletReaction<flameletGasHThermoPhysics> flameletGasHReaction;

    typedef FlameletReaction<flameletConstIncompressibleGasHThermoPhysics>
        flameletConstIncompressibleGasHReaction;

    typedef FlameletReaction<flameletIncompressibleGasHThermoPhysics>
        flameletIncompressibleGasHReaction;

    typedef FlameletReaction<flameletIcoPoly8HThermoPhysics> flameletIcoPoly8HReaction;

    typedef FlameletReaction<flameletConstFluidHThermoPhysics> flameletConstFluidHReaction;

    typedef FlameletReaction<flameletConstAdiabaticFluidHThermoPhysics>
        flameletConstAdiabaticFluidHReaction;

    typedef FlameletReaction<flameletConstHThermoPhysics> flameletConstHReaction;


    // internal ennergy based reactions
    typedef FlameletReaction<flameletConstGasEThermoPhysics> flameletConstGasEReaction;

    typedef FlameletReaction<flameletGasEThermoPhysics> flameletGasEReaction;

    typedef FlameletReaction<flameletConstIncompressibleGasEThermoPhysics>
        flameletConstIncompressibleGasEReaction;

    typedef FlameletReaction<flameletIncompressibleGasEThermoPhysics>
        flameletIncompressibleGasEReaction;

    typedef FlameletReaction<flameletIcoPoly8EThermoPhysics> flameletIcoPoly8EReaction;

    typedef FlameletReaction<flameletConstFluidEThermoPhysics> flameletConstFluidEReaction;

    typedef FlameletReaction<flameletConstAdiabaticFluidEThermoPhysics>
        flameletConstAdiabaticFluidEReaction;

    typedef FlameletReaction<flameletConstEThermoPhysics> flameletConstEReaction;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
