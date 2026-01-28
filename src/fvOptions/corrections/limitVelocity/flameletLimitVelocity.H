/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2019 OpenFOAM Foundation
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
    Foam::fv::flameletLimitVelocity

Description
    Limits the maximum velocity magnitude to the specified \c max value.

Usage
    Example usage:
    \verbatim
    limitU
    {
        type            flameletLimitVelocity;
        active          yes;

        selectionMode   all;
        max             100;
    }
    \endverbatim

SourceFiles
    flameletLimitVelocity.C

\*---------------------------------------------------------------------------*/

#ifndef flameletLimitVelocity_H
#define flameletLimitVelocity_H

#include "flameletCellSetOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                 Class flameletLimitVelocity Declaration
\*---------------------------------------------------------------------------*/

class flameletLimitVelocity
:
    public flameletCellSetOption
{
protected:

    // Protected data

        //- Velocity field name, default = U
        word UName_;

        //- Maximum velocity magnitude
        scalar max_;


public:

    //- Runtime type information
    TypeName("flameletLimitVelocity");


    // Constructors

        //- Construct from components
        flameletLimitVelocity
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Disallow default bitwise copy construction
        flameletLimitVelocity(const flameletLimitVelocity&) = delete;


    //- Destructor
    virtual ~flameletLimitVelocity()
    {}


    // Member Functions

        //- Read dictionary
        virtual bool read(const dictionary& dict);

        //- Correct the energy field
        virtual void correct(volVectorField& U);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const flameletLimitVelocity&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
