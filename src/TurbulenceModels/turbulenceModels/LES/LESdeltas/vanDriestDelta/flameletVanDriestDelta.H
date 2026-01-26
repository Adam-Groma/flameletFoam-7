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
    Foam::LESModels::flameletVanDriestDelta

Description
    Apply van Driest damping function to the specified geometric delta to
    improve near-wall behavior or LES SGS models.

SourceFiles
    flameletVanDriestDelta.C

\*---------------------------------------------------------------------------*/

#ifndef flameletVanDriestDelta_H
#define flameletVanDriestDelta_H

#include "flameletLESdelta.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                           Class flameletVanDriestDelta Declaration
\*---------------------------------------------------------------------------*/

class flameletVanDriestDelta
:
    public flameletLESdelta
{
    // Private Data

        autoPtr<flameletLESdelta> geometricDelta_;
        scalar kappa_;
        scalar Aplus_;
        scalar Cdelta_;
        label calcInterval_;


    // Private Member Functions

        // Calculate the delta values
        void calcDelta();


public:

    //- Runtime type information
    TypeName("flameletVanDriest");


    // Constructors

        //- Construct from name, flameletTurbulenceModel and dictionary
        flameletVanDriestDelta
        (
            const word& name,
            const flameletTurbulenceModel& turbulence,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        flameletVanDriestDelta(const flameletVanDriestDelta&) = delete;


    //- Destructor
    virtual ~flameletVanDriestDelta()
    {}


    // Member Functions

        //- Read the flameletLESdelta dictionary
        virtual void read(const dictionary&);

        // Correct values
        virtual void correct();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const flameletVanDriestDelta&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
