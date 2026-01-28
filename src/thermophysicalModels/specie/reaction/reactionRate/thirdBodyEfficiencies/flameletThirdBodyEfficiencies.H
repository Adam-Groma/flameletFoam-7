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
    Foam::flameletThirdBodyEfficiencies

Description
    Third body efficiencies

SourceFiles
    flameletThirdBodyEfficienciesI.H

\*---------------------------------------------------------------------------*/

#ifndef flameletThirdBodyEfficiencies_H
#define flameletThirdBodyEfficiencies_H

#include "scalarList.H"
#include "flameletSpeciesTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class flameletThirdBodyEfficiencies;
Ostream& operator<<(Ostream&, const flameletThirdBodyEfficiencies&);


/*---------------------------------------------------------------------------*\
                  Class flameletThirdBodyEfficiencies Declaration
\*---------------------------------------------------------------------------*/

class flameletThirdBodyEfficiencies
:
    public scalarList
{
    // Private Data

        const flameletSpeciesTable& species_;


public:

    // Constructors

        //- Construct from components
        inline flameletThirdBodyEfficiencies
        (
            const flameletSpeciesTable& species,
            const scalarList& efficiencies
        );

        //- Construct from dictionary
        inline flameletThirdBodyEfficiencies
        (
            const flameletSpeciesTable& species,
            const dictionary& dict
        );


    // Member Functions

        //- Calculate and return M, the concentration of the third-bodies
        inline scalar M(const scalarList& c) const;

        //- Write to stream
        inline void write(Ostream& os) const;



    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const flameletThirdBodyEfficiencies&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "flameletThirdBodyEfficienciesI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
