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
    Foam::LESModels::flameletIDDESDelta

Description
    flameletIDDESDelta used by the IDDES (improved low Re Spalart-Allmaras DES model)
    The min and max delta are calculated using the double distance of the min or
    max from the face centre to the cell centre.

SourceFiles
    flameletIDDESDelta.C

\*---------------------------------------------------------------------------*/

#ifndef flameletIDDESDelta_H
#define flameletIDDESDelta_H

#include "flameletMaxDeltaxyz.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                         Class flameletIDDESDelta Declaration
\*---------------------------------------------------------------------------*/

class flameletIDDESDelta
:
    public flameletLESdelta
{
    // Private Data

        flameletMaxDeltaxyz hmax_;

        scalar Cw_;


    // Private Member Functions

        //- Calculate the delta values
        void calcDelta();


public:

    //- Runtime type information
    TypeName("flameletIDDESDelta");


    // Constructors

        //- Construct from name, mesh and IOdictionary
        flameletIDDESDelta
        (
            const word& name,
            const flameletTurbulenceModel& turbulence,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        flameletIDDESDelta(const flameletIDDESDelta&) = delete;


    // Destructor

        ~flameletIDDESDelta()
        {}


    // Member Functions

        //- Read the flameletLESdelta dictionary
        void read(const dictionary&);

        //- Return the hmax delta field
        const volScalarField& hmax() const
        {
            return hmax_;
        }

        // Correct values
        void correct();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const flameletIDDESDelta&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
