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
    Foam::FlameletReactionProxy

Description
    Proxy version of FlameletReaction which provides dummy implementations of the
    abstract virtual functions.

    Used for read-construction and format conversion.

SourceFiles
    FlameletReactionProxy.C

\*---------------------------------------------------------------------------*/

#ifndef FlameletReactionProxy_H
#define FlameletReactionProxy_H

#include "FlameletReaction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class FlameletReactionProxy Declaration
\*---------------------------------------------------------------------------*/

template<class ReactionThermo>
class FlameletReactionProxy
:
    public FlameletReaction<ReactionThermo>
{

public:

    // Constructors

        //- Construct from components
        FlameletReactionProxy
        (
            const flameletSpeciesTable& species,
            const List<flameletSpecieCoeffs>& lhs,
            const List<flameletSpecieCoeffs>& rhs,
            const HashPtrTable<ReactionThermo>& thermoDatabase
        );

        //- Construct as copy given new flameletSpeciesTable
        FlameletReactionProxy
        (
            const FlameletReaction<ReactionThermo>&,
            const flameletSpeciesTable& species
        );

        //- Construct from dictionary
        FlameletReactionProxy
        (
            const flameletSpeciesTable& species,
            const HashPtrTable<ReactionThermo>& thermoDatabase,
            const dictionary& dict
        );

        //- Construct and return a clone
        virtual autoPtr<FlameletReaction<ReactionThermo>> clone() const;

        //- Construct and return a clone with new flameletSpeciesTable
        virtual autoPtr<FlameletReaction<ReactionThermo>> clone
        (
            const flameletSpeciesTable& species
        ) const;


    //- Destructor
    virtual ~FlameletReactionProxy()
    {}


    // Member Functions

        // FlameletReaction rate coefficients

            //- Forward rate constant
            virtual scalar kf
            (
                const scalar p,
                const scalar T,
                const scalarField& c
            ) const;

            //- Reverse rate constant from the given forward rate constant
            virtual scalar kr
            (
                const scalar kfwd,
                const scalar p,
                const scalar T,
                const scalarField& c
            ) const;

            //- Reverse rate constant
            virtual scalar kr
            (
                const scalar p,
                const scalar T,
                const scalarField& c
            ) const;


        // Jacobian coefficients

            //- Temperature derivative of forward rate
            virtual scalar dkfdT
            (
                const scalar p,
                const scalar T,
                const scalarField& c
            ) const;

            //- Temperature derivative of reverse rate
            virtual scalar dkrdT
            (
                const scalar p,
                const scalar T,
                const scalarField& c,
                const scalar dkfdT,
                const scalar kr
            ) const;

            //- Third-body efficiencies (beta = 1-alpha)
            //  non-empty only for third-body reactions
            //  with enhanced molecularity (alpha != 1)
            virtual const List<Tuple2<label, scalar>>& beta() const;

            //- Species concentration derivative of the pressure dependent term
            virtual void dcidc
            (
                const scalar p,
                const scalar T,
                const scalarField& c,
                scalarField& dcidc
            ) const;

            //- Temperature derivative of the pressure dependent term
            virtual scalar dcidT
            (
                const scalar p,
                const scalar T,
                const scalarField& c
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "FlameletReactionProxy.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
