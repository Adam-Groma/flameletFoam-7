/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2019 OpenFOAM Foundation
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
    Foam::flameletMichaelisMentenReactionRate

Description
    Michaelis-Menten reaction rate for enzymatic reactions.

    Reference:
    \verbatim
        Michaelis, L., & Menten, M. L. (1913).
        Die Kinetik der InwertinWirkung.
        Biochem., (49), 333-369.
    \endverbatim

SourceFiles
    flameletMichaelisMentenReactionRateI.H

\*---------------------------------------------------------------------------*/

#ifndef flameletMichaelisMentenReactionRate_H
#define flameletMichaelisMentenReactionRate_H

#include "FlameletReaction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class flameletMichaelisMentenReactionRate;

Ostream& operator<<(Ostream&, const flameletMichaelisMentenReactionRate&);


/*---------------------------------------------------------------------------*\
               Class flameletMichaelisMentenReactionRate Declaration
\*---------------------------------------------------------------------------*/

class flameletMichaelisMentenReactionRate
{
    // Private Data

        //- List of species present in reaction system
        const flameletSpeciesTable& species_;

        //- The maximum reaction rate at saturating substrate concentration
        scalar Vmax_;

        //- The Michaelis constant
        //  the substrate concentration at which the reaction rate is half Vmax_
        scalar Km_;

        //- The substrate specie index
        label s_;

        List<Tuple2<label, scalar>> beta_;


public:

    // Constructors

        //- Construct from dictionary
        inline flameletMichaelisMentenReactionRate
        (
            const flameletSpeciesTable& species,
            const dictionary& dict
        );


    // Member Functions

        //- Return the type name
        static word type()
        {
            return "flameletMichaelisMenten";
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
        inline const List<Tuple2<label, scalar>>& beta() const;

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

        inline friend Ostream& operator<<
        (
            Ostream&,
            const flameletMichaelisMentenReactionRate&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "flameletMichaelisMentenReactionRateI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
