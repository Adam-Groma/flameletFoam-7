/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "FlameletReactionList.H"
#include "IFstream.H"
#include "SLPtrList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::FlameletReactionList<ThermoType>::FlameletReactionList
(
    const flameletSpeciesTable& species,
    const HashPtrTable<ThermoType>& thermoDb
)
:
    SLPtrList<FlameletReaction<ThermoType>>(),
    species_(species),
    thermoDb_(thermoDb),
    dict_(dictionary::null)
{}


template<class ThermoType>
Foam::FlameletReactionList<ThermoType>::FlameletReactionList
(
    const flameletSpeciesTable& species,
    const HashPtrTable<ThermoType>& thermoDb,
    const dictionary& dict
)
:
    SLPtrList<FlameletReaction<ThermoType>>(),
    species_(species),
    thermoDb_(thermoDb),
    dict_(dict)
{
    readReactionDict();
}


template<class ThermoType>
Foam::FlameletReactionList<ThermoType>::FlameletReactionList(const FlameletReactionList& reactions)
:
    SLPtrList<FlameletReaction<ThermoType>>(reactions),
    species_(reactions.species_),
    thermoDb_(reactions.thermoDb_),
    dict_(reactions.dict_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::FlameletReactionList<ThermoType>::~FlameletReactionList()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
bool Foam::FlameletReactionList<ThermoType>::readReactionDict()
{
    const dictionary& reactions(dict_.subDict("reactions"));

    // Set general temperature limits from the dictionary
    FlameletReaction<ThermoType>::TlowDefault =
        dict_.lookupOrDefault<scalar>("Tlow", 0);

    FlameletReaction<ThermoType>::ThighDefault =
        dict_.lookupOrDefault<scalar>("Thigh", great);

    forAllConstIter(dictionary, reactions, iter)
    {
        const word reactionName = iter().keyword();

        this->append
        (
            FlameletReaction<ThermoType>::New
            (
                species_,
                thermoDb_,
                reactions.subDict(reactionName)
            ).ptr()
        );
    }

    return true;
}


template<class ThermoType>
void Foam::FlameletReactionList<ThermoType>::write(Ostream& os) const
{
    os  << "reactions" << nl;
    os  << token::BEGIN_BLOCK << incrIndent << nl;

    forAllConstIter(typename SLPtrList<FlameletReaction<ThermoType>>, *this, iter)
    {
        const FlameletReaction<ThermoType>& r = iter();
        os  << indent << r.name() << nl
            << indent << token::BEGIN_BLOCK << incrIndent << nl;
        writeEntry(os, "type", r.type());
        r.write(os);
        os  << decrIndent << indent << token::END_BLOCK << nl;
    }

    os << decrIndent << token::END_BLOCK << nl;
}


// ************************************************************************* //
