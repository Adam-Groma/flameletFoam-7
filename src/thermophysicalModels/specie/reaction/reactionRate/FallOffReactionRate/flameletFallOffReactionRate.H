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
    Foam::flameletFallOffReactionRate

Description
    General class for handling unimolecular/recombination fall-off reactions.

SourceFiles
    flameletFallOffReactionRateI.H

\*---------------------------------------------------------------------------*/

#ifndef flameletFallOffReactionRate_H
#define flameletFallOffReactionRate_H

#include "flameletThirdBodyEfficiencies.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class ReactionRate, class FallOffFunction>
class flameletFallOffReactionRate;

template<class ReactionRate, class FallOffFunction>
inline Ostream& operator<<
(
    Ostream&,
    const flameletFallOffReactionRate<ReactionRate, FallOffFunction>&
);


/*---------------------------------------------------------------------------*\
                    Class flameletFallOffReactionRate Declaration
\*---------------------------------------------------------------------------*/

template<class ReactionRate, class FallOffFunction>
class flameletFallOffReactionRate
{
    // Private Data

        ReactionRate k0_;
        ReactionRate kInf_;
        FallOffFunction F_;
        flameletThirdBodyEfficiencies flameletThirdBodyEfficiencies_;
        List<Tuple2<label, scalar>> beta_;


public:

    // Constructors

        //- Construct from components
        inline flameletFallOffReactionRate
        (
            const ReactionRate& k0,
            const ReactionRate& kInf,
            const FallOffFunction& F,
            const flameletThirdBodyEfficiencies& tbes
        );

        //- Construct from dictionary
        inline flameletFallOffReactionRate
        (
            const flameletSpeciesTable& species,
            const dictionary& dict
        );


    // Member Functions

        //- Return the type name
        static word type()
        {
            return ReactionRate::type() + FallOffFunction::type() + "FallOff";
        }

        inline scalar operator()
        (
            const scalar p,
            const scalar T,
            const scalarField& c
        ) const;

        inline scalar ddT
        (
            const scalar p,
            const scalar T,
            const scalarField& c
        ) const;

        //- Third-body efficiencies (beta = 1-alpha)
        //  non-empty only for third-body reactions
        //  with enhanced molecularity (alpha != 1)
        inline const List<Tuple2<label, scalar>>& beta() const
        {
            return beta_;
        }

        //- Species concentration derivative of the pressure dependent term
        inline void dcidc
        (
            const scalar p,
            const scalar T,
            const scalarField& c,
            scalarField& dcidc
        ) const;

        //- Temperature derivative of the pressure dependent term
        inline scalar dcidT
        (
            const scalar p,
            const scalar T,
            const scalarField& c
        ) const;

        //- Write to stream
        inline void write(Ostream& os) const;


    // Ostream Operator

        friend Ostream& operator<< <ReactionRate, FallOffFunction>
        (
            Ostream&,
            const flameletFallOffReactionRate<ReactionRate, FallOffFunction>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "flameletFallOffReactionRateI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
