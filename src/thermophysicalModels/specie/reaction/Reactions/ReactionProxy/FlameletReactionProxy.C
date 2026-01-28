/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2019 OpenFOAM Foundation
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

#include "FlameletReactionProxy.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ReactionThermo>
Foam::FlameletReactionProxy<ReactionThermo>::FlameletReactionProxy
(
    const flameletSpeciesTable& species,
    const List<flameletSpecieCoeffs>& lhs,
    const List<flameletSpecieCoeffs>& rhs,
    const HashPtrTable<ReactionThermo>& thermoDatabase
)
:
    FlameletReaction<ReactionThermo>
    (
        species,
        lhs,
        rhs,
        thermoDatabase
    )
{}


template<class ReactionThermo>
Foam::FlameletReactionProxy<ReactionThermo>::FlameletReactionProxy
(
    const FlameletReaction<ReactionThermo>& r,
    const flameletSpeciesTable& species
)
:
    FlameletReaction<ReactionThermo>
    (
        r,
        species
    )
{}


template<class ReactionThermo>
Foam::FlameletReactionProxy<ReactionThermo>::FlameletReactionProxy
(
    const flameletSpeciesTable& species,
    const HashPtrTable<ReactionThermo>& thermoDatabase,
    const dictionary& dict
)
:
    FlameletReaction<ReactionThermo>
    (
        species,
        thermoDatabase,
        dict
    )
{}


template<class ReactionThermo>
Foam::autoPtr<Foam::FlameletReaction<ReactionThermo>>
Foam::FlameletReactionProxy<ReactionThermo>::clone() const
{
    NotImplemented;
    return autoPtr<Foam::FlameletReaction<ReactionThermo>>();
}


template<class ReactionThermo>
Foam::autoPtr<Foam::FlameletReaction<ReactionThermo>>
Foam::FlameletReactionProxy<ReactionThermo>::clone
(
    const flameletSpeciesTable& species
) const
{
    NotImplemented;
    return autoPtr<FlameletReaction<ReactionThermo>>();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ReactionThermo>
Foam::scalar Foam::FlameletReactionProxy<ReactionThermo>::kf
(
    const scalar p,
    const scalar T,
    const scalarField& c
) const
{
    NotImplemented;
    return 0;
}


template<class ReactionThermo>
Foam::scalar Foam::FlameletReactionProxy<ReactionThermo>::kr
(
    const scalar kfwd,
    const scalar p,
    const scalar T,
    const scalarField& c
) const
{
    NotImplemented;
    return 0;
}


template<class ReactionThermo>
Foam::scalar Foam::FlameletReactionProxy<ReactionThermo>::kr
(
    const scalar p,
    const scalar T,
    const scalarField& c
) const
{
    NotImplemented;
    return 0;
}


template<class ReactionThermo>
Foam::scalar Foam::FlameletReactionProxy<ReactionThermo>::dkfdT
(
    const scalar p,
    const scalar T,
    const scalarField& c
) const
{
    NotImplemented;
    return 0;
}


template<class ReactionThermo>
Foam::scalar Foam::FlameletReactionProxy<ReactionThermo>::dkrdT
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const scalar dkfdT,
    const scalar kr
) const
{
    NotImplemented;
    return 0;
}


template<class ReactionThermo>
const Foam::List<Foam::Tuple2<Foam::label, Foam::scalar>>&
Foam::FlameletReactionProxy<ReactionThermo>::beta() const
{
    NotImplemented;
    return NullObjectRef<List<Tuple2<label, scalar>>>();
}


template<class ReactionThermo>
void Foam::FlameletReactionProxy<ReactionThermo>::dcidc
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    scalarField& dcidc
) const
{
    NotImplemented;
}


template<class ReactionThermo>
Foam::scalar Foam::FlameletReactionProxy<ReactionThermo>::dcidT
(
    const scalar p,
    const scalar T,
    const scalarField& c
) const
{
    NotImplemented;
    return 0;
}


// ************************************************************************* //
