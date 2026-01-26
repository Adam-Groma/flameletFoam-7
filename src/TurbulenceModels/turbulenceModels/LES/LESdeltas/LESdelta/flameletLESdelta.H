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
    Foam::flameletLESdelta

Description
    Abstract base class for LES deltas

SourceFiles
    flameletLESdelta.C

\*---------------------------------------------------------------------------*/

#ifndef flameletLESdelta_H
#define flameletLESdelta_H

#include "flameletTurbulenceModel.H"
#include "volFields.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class flameletLESdelta Declaration
\*---------------------------------------------------------------------------*/

class flameletLESdelta
{
protected:

    // Protected data

        const flameletTurbulenceModel& turbulenceModel_;

        volScalarField delta_;


public:

    //- Runtime type information
    TypeName("flameletLESdelta");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            flameletLESdelta,
            dictionary,
            (
                const word& name,
                const flameletTurbulenceModel& turbulence,
                const dictionary& dict
            ),
            (name, turbulence, dict)
        );


    // Constructors

        //- Construct from name, flameletTurbulenceModel and dictionary
        flameletLESdelta
        (
            const word& name,
            const flameletTurbulenceModel& turbulence
        );

        //- Disallow default bitwise copy construction
        flameletLESdelta(const flameletLESdelta&) = delete;


    // Selectors

        //- Return a reference to the selected LES delta
        static autoPtr<flameletLESdelta> New
        (
            const word& name,
            const flameletTurbulenceModel& turbulence,
            const dictionary& dict
        );

        //- Return a reference to the selected LES delta
        static autoPtr<flameletLESdelta> New
        (
            const word& name,
            const flameletTurbulenceModel& turbulence,
            const dictionary& dict,
            const dictionaryConstructorTable&
        );


    //- Destructor
    virtual ~flameletLESdelta()
    {}


    // Member Functions

        //- Return flameletTurbulenceModel reference
        const flameletTurbulenceModel& turbulence() const
        {
            return turbulenceModel_;
        }

        //- Read the flameletLESdelta dictionary
        virtual void read(const dictionary&) = 0;

        // Correct values
        virtual void correct() = 0;


    // Member Operators

        void operator=(const flameletLESdelta&) = delete;

        virtual operator const volScalarField&() const
        {
            return delta_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
