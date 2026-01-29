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
    Foam::flameletEulerImplicit

Description
    An Euler implicit solver for chemistry

SourceFiles
    flameletEulerImplicit.C

\*---------------------------------------------------------------------------*/

#ifndef flameletEulerImplicit_H
#define flameletEulerImplicit_H

#include "flameletChemistrySolver.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template <class Type>
class simpleMatrix;

/*---------------------------------------------------------------------------*\
                       Class flameletEulerImplicit Declaration
\*---------------------------------------------------------------------------*/

template<class ChemistryModel>
class flameletEulerImplicit
:
    public flameletChemistrySolver<ChemistryModel>
{
    // Private Data

        //- Coefficients dictionary
        dictionary coeffsDict_;


        // Model constants

            //- Chemistry timescale
            scalar cTauChem_;

            //- Equilibrium rate limiter flag (on/off)
            Switch eqRateLimiter_;

        // Solver data
        mutable scalarField cTp_;


public:

    //- Runtime type information
    TypeName("flameletEulerImplicit");


    // Constructors

        //- Construct from thermo
        flameletEulerImplicit(typename ChemistryModel::reactionThermo& thermo);


    //- Destructor
    virtual ~flameletEulerImplicit();


    // Member Functions

    void updateRRInReactionI
    (
        const label index,
        const scalar pr,
        const scalar pf,
        const scalar corr,
        const label lRef,
        const label rRef,
        const scalar p,
        const scalar T,
        simpleMatrix<scalar>& RR
    ) const;

        //- Update the concentrations and return the chemical time
        virtual void solve
        (
            scalarField& c,
            scalar& T,
            scalar& p,
            scalar& deltaT,
            scalar& subDeltaT
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "flameletEulerImplicit.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
