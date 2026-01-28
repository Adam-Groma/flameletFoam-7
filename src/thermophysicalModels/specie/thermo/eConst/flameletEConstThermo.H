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
    Foam::flameletEConstThermo

Description
    Constant properties thermodynamics package templated on an equation of
    state

SourceFiles
    flameletEConstThermoI.H
    flameletEConstThermo.C

\*---------------------------------------------------------------------------*/

#ifndef flameletEConstThermo_H
#define flameletEConstThermo_H

#include "flameletThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class EquationOfState> class flameletEConstThermo;

template<class EquationOfState>
inline flameletEConstThermo<EquationOfState> operator+
(
    const flameletEConstThermo<EquationOfState>&,
    const flameletEConstThermo<EquationOfState>&
);

template<class EquationOfState>
inline flameletEConstThermo<EquationOfState> operator*
(
    const scalar,
    const flameletEConstThermo<EquationOfState>&
);

template<class EquationOfState>
inline flameletEConstThermo<EquationOfState> operator==
(
    const flameletEConstThermo<EquationOfState>&,
    const flameletEConstThermo<EquationOfState>&
);

template<class EquationOfState>
Ostream& operator<<
(
    Ostream&,
    const flameletEConstThermo<EquationOfState>&
);


/*---------------------------------------------------------------------------*\
                        Class flameletEConstThermo Declaration
\*---------------------------------------------------------------------------*/

template<class EquationOfState>
class flameletEConstThermo
:
    public EquationOfState
{
    // Private Data

        //- Heat capacity at constant volume
        //  Note: input in [J/kg/K], but internally uses [J/kmol/K]
        scalar Cv_;

        //- Heat of formation
        //  Note: input in [J/kg], but internally uses [J/kmol]
        scalar Hf_;


    // Private Member Functions

        //- Construct from components
        inline flameletEConstThermo
        (
            const EquationOfState& st,
            const scalar cv,
            const scalar hf
        );


public:

    // Constructors

        //- Construct from dictionary
        flameletEConstThermo(const dictionary& dict);

        //- Construct as named copy
        inline flameletEConstThermo(const word&, const flameletEConstThermo&);

        //- Construct and return a clone
        inline autoPtr<flameletEConstThermo> clone() const;

        // Selector from dictionary
        inline static autoPtr<flameletEConstThermo> New(const dictionary& dict);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "eConst<" + EquationOfState::typeName() + '>';
        }

        //- Limit the temperature to be in the range Tlow_ to Thigh_
        inline scalar limit(const scalar T) const;


        // Fundamental properties

            //- Heat capacity at constant volume [J/kg/K]
            inline scalar Cv(const scalar p, const scalar T) const;

            //- Sensible internal energy [J/kg]
            inline scalar Es(const scalar p, const scalar T) const;

            //- Absolute internal energy [J/kg]
            inline scalar Ea(const scalar p, const scalar T) const;

            //- Chemical enthalpy [J/kg]
            inline scalar Hc() const;

            //- Entropy [J/kg/K]
            inline scalar S(const scalar p, const scalar T) const;

            #include "flameletEtoHthermo.H"


        // Derivative term used for Jacobian

            //- Derivative of Gibbs free energy w.r.t. temperature
            inline scalar dGdT(const scalar p, const scalar T) const;

            //- Temperature derivative of heat capacity at constant pressure
            inline scalar dCpdT(const scalar p, const scalar T) const;


        // I-O

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member Operators

        inline void operator+=(const flameletEConstThermo&);


    // Friend operators

        friend flameletEConstThermo operator+ <EquationOfState>
        (
            const flameletEConstThermo&,
            const flameletEConstThermo&
        );

        friend flameletEConstThermo operator* <EquationOfState>
        (
            const scalar,
            const flameletEConstThermo&
        );

        friend flameletEConstThermo operator== <EquationOfState>
        (
            const flameletEConstThermo&,
            const flameletEConstThermo&
        );


    // IOstream Operators

        friend Ostream& operator<< <EquationOfState>
        (
            Ostream&,
            const flameletEConstThermo&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "flameletEConstThermoI.H"

#ifdef NoRepository
    #include "flameletEConstThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
