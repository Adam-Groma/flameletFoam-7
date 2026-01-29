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
    Foam::flameletChemistryReader

Description
    Abstract class for reading chemistry

SourceFiles
    flameletChemistryReader.C

\*---------------------------------------------------------------------------*/

#ifndef flameletChemistryReader_H
#define flameletChemistryReader_H

#include "typeInfo.H"
#include "flameletSpecieElement.H"
#include "FlameletReaction.H"
#include "FlameletReactionList.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

typedef HashTable<List<flameletSpecieElement>> speciesCompositionTable;


/*---------------------------------------------------------------------------*\
                      Class flameletChemistryReader Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class flameletChemistryReader
{
public:

    //- Runtime type information
    TypeName("flameletChemistryReader");

    //- The type of thermo package the reader was instantiated for
    typedef ThermoType thermoType;


    // Constructors

        //- Construct null
        flameletChemistryReader()
        {}

        //- Disallow default bitwise copy construction
        flameletChemistryReader(const flameletChemistryReader&) = delete;


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            flameletChemistryReader,
            dictionary,
            (
                const dictionary& thermoDict,
                flameletSpeciesTable& species
            ),
            (thermoDict, species)
        );


    // Selectors

        //- Select constructed from dictionary
        static autoPtr<flameletChemistryReader> New
        (
            const dictionary& thermoDict,
            flameletSpeciesTable& species
        );


    //- Destructor
    virtual ~flameletChemistryReader()
    {}


    // Member Functions

        //- Return access to the list of species
        virtual const flameletSpeciesTable& species() const = 0;

        //- Table of species composition
        virtual const speciesCompositionTable& specieComposition() const = 0;

        //- Return access to the thermo packages
        virtual const HashPtrTable<ThermoType>& speciesThermo() const = 0;

        //- Return access to the list of reactions
        virtual const FlameletReactionList<ThermoType>& reactions() const = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const flameletChemistryReader&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "flameletChemistryReader.C"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#define flameletMakeChemistryReader(Thermo)                                            \
    defineTemplateTypeNameAndDebug(flameletChemistryReader<Thermo>, 0);                \
    defineTemplateRunTimeSelectionTable(flameletChemistryReader<Thermo>, dictionary)


#define flameletMakeChemistryReaderType(Reader, Thermo)                                \
    defineNamedTemplateTypeNameAndDebug(Reader<Thermo>, 0);                    \
    flameletChemistryReader<Thermo>::adddictionaryConstructorToTable<Reader<Thermo>> \
        add##Reader##Thermo##ConstructorToTable_


// for non-templated chemistry readers
#define addChemistryReaderType(Reader, Thermo)                                 \
    defineTypeNameAndDebug(Reader, 0);                                         \
    flameletChemistryReader<Thermo>::adddictionaryConstructorToTable<Reader>           \
        add##Reader##Thermo##ConstructorToTable_


// for templated chemistry readers
#define addTemplateChemistryReaderType(Reader, Thermo)                         \
    defineNamedTemplateTypeNameAndDebug(Reader, 0);                            \
    flameletChemistryReader<Thermo>::adddictionaryConstructorToTable<Reader>           \
        add##Reader##Thermo##ConstructorToTable_


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
