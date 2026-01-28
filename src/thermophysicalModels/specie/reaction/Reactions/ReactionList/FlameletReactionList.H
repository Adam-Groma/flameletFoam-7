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

Class
    Foam::FlameletReactionList

Description
    List of templated reactions

SourceFiles
    FlameletReactionList.C

\*---------------------------------------------------------------------------*/

#ifndef FlameletReactionList_H
#define FlameletReactionList_H

#include "PtrList.H"
#include "SLPtrList.H"
#include "flameletSpeciesTable.H"
#include "HashPtrTable.H"
#include "FlameletReaction.H"
#include "fileName.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class FlameletReactionList Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class FlameletReactionList
:
    public SLPtrList<FlameletReaction<ThermoType>>
{
    // Private Data

        //- Reference to the table of species
        const flameletSpeciesTable& species_;

        //- Reference to the thermo database
        const HashPtrTable<ThermoType>& thermoDb_;

        //- The dictionary used for construction
        const dictionary dict_;


public:

    // Constructors

        //- Construct null
        FlameletReactionList
        (
            const flameletSpeciesTable& species,
            const HashPtrTable<ThermoType>& thermoDatabase
        );

        //- Construct from dictionary
        FlameletReactionList
        (
            const flameletSpeciesTable& species,
            const HashPtrTable<ThermoType>& thermoDatabase,
            const dictionary& dict
        );

        //- Construct copy
        FlameletReactionList(const FlameletReactionList<ThermoType>& reactions);


    //- Destructor
    ~FlameletReactionList();


    // Public Member Functions

        //- Read reactions from dictionary
        bool readReactionDict();

        //- Write
        void write(Ostream& os) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const FlameletReactionList&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "FlameletReactionList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
