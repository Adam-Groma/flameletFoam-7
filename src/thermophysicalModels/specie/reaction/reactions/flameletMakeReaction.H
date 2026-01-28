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

InClass
    Foam::flameletMakeReaction

Description
    Macros for instantiating reactions on given thermo packages

\*---------------------------------------------------------------------------*/

#ifndef flameletMakeReaction_H
#define flameletMakeReaction_H

#include "FlameletReaction.H"

#include "FlameletIrreversibleReaction.H"
#include "FlameletReversibleReaction.H"
#include "FlameletNonEquilibriumReversibleReaction.H"

#include "flameletThermo.H"

#include "flameletSutherlandTransport.H"
#include "flameletJanafThermo.H"
#include "flameletPerfectGas.H"

#include "flameletPolynomialTransport.H"
#include "flameletHPolynomialThermo.H"
#include "flameletIcoPolynomial.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define flameletMakeReaction(Thermo, ReactionType, ReactionRate)                       \
                                                                               \
    typedef FlameletReaction<Thermo> Reaction##Thermo;                                 \
                                                                               \
    typedef ReactionType<FlameletReaction, Thermo, ReactionRate>                       \
        ReactionType##Thermo##ReactionRate;                                    \
                                                                               \
    template<>                                                                 \
    const word ReactionType##Thermo##ReactionRate::typeName                    \
    (                                                                          \
        ReactionType::typeName_()                                              \
      + ReactionRate::type()                                                   \
      + Reaction##Thermo::typeName_()                                          \
    );                                                                         \
                                                                               \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        Reaction##Thermo,                                                      \
        ReactionType##Thermo##ReactionRate,                                    \
        dictionary                                                             \
    );


#define flameletMakePressureDependentReaction\
(                                                                              \
    Thermo,                                                                    \
    FlameletReaction,                                                                  \
    PressureDependentReactionRate,                                             \
    ReactionRate,                                                              \
    FallOffFunction                                                            \
)                                                                              \
                                                                               \
    typedef PressureDependentReactionRate<ReactionRate, FallOffFunction>       \
        PressureDependentReactionRate##ReactionRate##FallOffFunction;          \
                                                                               \
    flameletMakeReaction                                                               \
    (                                                                          \
        Thermo,                                                                \
        FlameletReaction,                                                              \
        PressureDependentReactionRate##ReactionRate##FallOffFunction           \
    )


#define flameletMakeIRReactions(Thermo, ReactionRate)                                  \
                                                                               \
    flameletMakeReaction(Thermo, FlameletIrreversibleReaction, ReactionRate)                   \
                                                                               \
    flameletMakeReaction(Thermo, FlameletReversibleReaction, ReactionRate)


#define flameletMakeIRNReactions(Thermo, ReactionRate)                                 \
                                                                               \
    flameletMakeIRReactions(Thermo, ReactionRate)                                      \
                                                                               \
    flameletMakeReaction(Thermo, FlameletNonEquilibriumReversibleReaction, ReactionRate)


#define flameletMakePressureDependentReactions(Thermo, ReactionRate, FallOffFunction)  \
                                                                               \
    flameletMakePressureDependentReaction                                              \
    (                                                                          \
        Thermo,                                                                \
        FlameletIrreversibleReaction,                                                  \
        flameletFallOffReactionRate,                                                   \
        ReactionRate,                                                          \
        FallOffFunction                                                        \
    )                                                                          \
                                                                               \
    flameletMakePressureDependentReaction                                              \
    (                                                                          \
        Thermo,                                                                \
        FlameletReversibleReaction,                                                    \
        flameletFallOffReactionRate,                                                   \
        ReactionRate,                                                          \
        FallOffFunction                                                        \
    )                                                                          \
                                                                               \
    flameletMakePressureDependentReaction                                              \
    (                                                                          \
        Thermo,                                                                \
        FlameletIrreversibleReaction,                                                  \
        flameletChemicallyActivatedReactionRate,                                       \
        ReactionRate,                                                          \
        FallOffFunction                                                        \
    )                                                                          \
                                                                               \
    flameletMakePressureDependentReaction                                              \
    (                                                                          \
        Thermo,                                                                \
        FlameletReversibleReaction,                                                    \
        flameletChemicallyActivatedReactionRate,                                       \
        ReactionRate,                                                          \
        FallOffFunction                                                        \
    )


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
