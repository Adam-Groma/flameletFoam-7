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
    Foam::FlameletReaction

Description
    Simple extension of ReactionThermo to handle reaction kinetics in addition
    to the equilibrium thermodynamics already handled.

SourceFiles
    FlameletReactionI.H
    FlameletReaction.C

\*---------------------------------------------------------------------------*/

#ifndef FlameletReaction_H
#define FlameletReaction_H

#include "flameletSpeciesTable.H"
#include "flameletSpecieCoeffs.H"
#include "HashPtrTable.H"
#include "scalarField.H"
#include "simpleMatrix.H"
#include "Tuple2.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class ReactionThermo>
class FlameletReaction;

template<class ReactionThermo>
inline Ostream& operator<<(Ostream&, const FlameletReaction<ReactionThermo>&);


/*---------------------------------------------------------------------------*\
                           Class FlameletReaction Declaration
\*---------------------------------------------------------------------------*/

template<class ReactionThermo>
class FlameletReaction
:
    public ReactionThermo::thermoType
{

public:

    // Static data

        //- Number of un-named reactions
        static label nUnNamedReactions;

        //- Default temperature limits of applicability of reaction rates
        static scalar TlowDefault, ThighDefault;

private:

    // Private Data

        //- Name of reaction
        const word name_;

        //- List of specie names present in reaction system
        const flameletSpeciesTable& species_;

        //- Temperature limits of applicability of reaction rates
        scalar Tlow_, Thigh_;

        //- Specie info for the left-hand-side of the reaction
        List<flameletSpecieCoeffs> lhs_;

        //- Specie info for the right-hand-side of the reaction
        List<flameletSpecieCoeffs> rhs_;


    // Private Member Functions

        //- Construct reaction thermo
        void setThermo(const HashPtrTable<ReactionThermo>& thermoDatabase);

        //- Return new reaction ID for un-named reactions
        label getNewReactionID();


public:

    //- Runtime type information
    TypeName("FlameletReaction");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            FlameletReaction,
            dictionary,
            (
                const flameletSpeciesTable& species,
                const HashPtrTable<ReactionThermo>& thermoDatabase,
                const dictionary& dict
            ),
            (species, thermoDatabase, dict)
        );


    // Constructors

        //- Construct from components
        FlameletReaction
        (
            const flameletSpeciesTable& species,
            const List<flameletSpecieCoeffs>& lhs,
            const List<flameletSpecieCoeffs>& rhs,
            const HashPtrTable<ReactionThermo>& thermoDatabase
        );

        //- Construct as copy given new flameletSpeciesTable
        FlameletReaction(const FlameletReaction<ReactionThermo>&, const flameletSpeciesTable& species);

        //- Construct from dictionary
        FlameletReaction
        (
            const flameletSpeciesTable& species,
            const HashPtrTable<ReactionThermo>& thermoDatabase,
            const dictionary& dict
        );

        //- Construct and return a clone
        virtual autoPtr<FlameletReaction<ReactionThermo>> clone() const = 0;

        //- Construct and return a clone with new flameletSpeciesTable
        virtual autoPtr<FlameletReaction<ReactionThermo>> clone
        (
            const flameletSpeciesTable& species
        ) const = 0;


    // Selectors

        //- Return a pointer to new patchField created on freestore from dict
        static autoPtr<FlameletReaction<ReactionThermo>> New
        (
            const flameletSpeciesTable& species,
            const HashPtrTable<ReactionThermo>& thermoDatabase,
            const dictionary& dict
        );


    //- Destructor
    virtual ~FlameletReaction()
    {}


    // Member Functions

        // Access

            //- Return the name of the reaction
            inline const word& name() const;

            //- Return the lower temperature limit for the reaction
            inline scalar Tlow() const;

            //- Return the upper temperature limit for the reaction
            inline scalar Thigh() const;

            //- Return the components of the left hand side
            inline const List<flameletSpecieCoeffs>& lhs() const;

            //- Return the components of the right hand side
            inline const List<flameletSpecieCoeffs>& rhs() const;

            //- Return the specie list
            const flameletSpeciesTable& species() const;


        // FlameletReaction rate coefficients

            //- Forward reaction rate
            void ddot
            (
                const scalar p,
                const scalar T,
                const scalarField& c,
                scalarField& d
            ) const;

            //- Backward reaction rate
            void fdot
            (
                const scalar p,
                const scalar T,
                const scalarField& c,
                scalarField& f
            ) const;

            //- Net reaction rate for individual species
            void omega
            (
                const scalar p,
                const scalar T,
                const scalarField& c,
                scalarField& dcdt
            ) const;

            //- Net reaction rate
            scalar omega
            (
                const scalar p,
                const scalar T,
                const scalarField& c,
                scalar& pf,
                scalar& cf,
                label& lRef,
                scalar& pr,
                scalar& cr,
                label& rRef
            ) const;

        // FlameletReaction rate coefficients

            //- Forward rate constant
            virtual scalar kf
            (
                const scalar p,
                const scalar T,
                const scalarField& c
            ) const = 0;

            //- Reverse rate constant from the given forward rate constant
            virtual scalar kr
            (
                const scalar kfwd,
                const scalar p,
                const scalar T,
                const scalarField& c
            ) const = 0;

            //- Reverse rate constant
            virtual scalar kr
            (
                const scalar p,
                const scalar T,
                const scalarField& c
            ) const = 0;


        // Jacobian coefficients

            //- Derivative of the net reaction rate for each species involved
            //  w.r.t. the species concentration
            void dwdc
            (
                const scalar p,
                const scalar T,
                const scalarField& c,
                scalarSquareMatrix& J,
                scalarField& dcdt,
                scalar& omegaI,
                scalar& kfwd,
                scalar& kbwd,
                const bool reduced,
                const List<label>& c2s
            ) const;

            //- Derivative of the net reaction rate for each species involved
            //  w.r.t. the temperature
            void dwdT
            (
                const scalar p,
                const scalar T,
                const scalarField& c,
                const scalar omegaI,
                const scalar kfwd,
                const scalar kbwd,
                scalarSquareMatrix& J,
                const bool reduced,
                const List<label>& c2s,
                const label indexT
            ) const;

            //- Temperature derivative of forward rate
            virtual scalar dkfdT
            (
                const scalar p,
                const scalar T,
                const scalarField& c
            ) const = 0;

            //- Temperature derivative of reverse rate
            virtual scalar dkrdT
            (
                const scalar p,
                const scalar T,
                const scalarField& c,
                const scalar dkfdT,
                const scalar kr
            ) const = 0;

            //- Third-body efficiencies (beta = 1-alpha)
            //  non-empty only for third-body reactions
            //  with enhanced molecularity (alpha != 1)
            virtual const List<Tuple2<label, scalar>>& beta() const = 0;

            //- Species concentration derivative of the pressure dependent term
            virtual void dcidc
            (
                const scalar p,
                const scalar T,
                const scalarField& c,
                scalarField& dcidc
            ) const = 0;

            //- Temperature derivative of the pressure dependent term
            virtual scalar dcidT
            (
                const scalar p,
                const scalar T,
                const scalarField& c
            ) const = 0;


        //- Write
        virtual void write(Ostream&) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const FlameletReaction<ReactionThermo>&) = delete;


    // Ostream Operator

        friend Ostream& operator<< <ReactionThermo>
        (
            Ostream&,
            const FlameletReaction<ReactionThermo>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "FlameletReactionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "FlameletReaction.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
