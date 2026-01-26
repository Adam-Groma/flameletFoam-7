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
    Foam::LESModels::flameletPrandtlDelta

Description
    Apply Prandtl mixing-length based damping function to the specified
    geometric delta to improve near-wall behavior or LES SGS models.

    \verbatim
        delta = min(geometricDelta, (kappa/Cdelta)*y)
    \endverbatim

    Example specification in the turbulenceProperties dictionary:
    \verbatim
    delta           flameletPrandtl;

    kappa           0.41;

    flameletPrandtlCoeffs
    {
        delta   flameletCubeRootVol;

        flameletCubeRootVolCoeffs
        {
            deltaCoeff      1;
        }

        Cdelta          0.158;
    }
    \endverbatim

SourceFiles
    flameletPrandtlDelta.C

\*---------------------------------------------------------------------------*/

#ifndef flameletPrandtlDelta_H
#define flameletPrandtlDelta_H

#include "flameletLESdelta.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                           Class flameletPrandtlDelta Declaration
\*---------------------------------------------------------------------------*/

class flameletPrandtlDelta
:
    public flameletLESdelta
{
    // Private Data

        autoPtr<flameletLESdelta> geometricDelta_;
        scalar kappa_;
        scalar Cdelta_;


    // Private Member Functions

        // Calculate the delta values
        void calcDelta();


public:

    //- Runtime type information
    TypeName("flameletPrandtl");


    // Constructors

        //- Construct from name, flameletTurbulenceModel and dictionary
        flameletPrandtlDelta
        (
            const word& name,
            const flameletTurbulenceModel& turbulence,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        flameletPrandtlDelta(const flameletPrandtlDelta&) = delete;


    //- Destructor
    virtual ~flameletPrandtlDelta()
    {}


    // Member Functions

        //- Read the flameletLESdelta dictionary
        virtual void read(const dictionary&);

        // Correct values
        virtual void correct();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const flameletPrandtlDelta&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
