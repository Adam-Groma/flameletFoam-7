/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2019 OpenFOAM Foundation
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
    Foam::flameletHRefConstThermo

Description
    Constant properties thermodynamics package
    templated into the EquationOfState.

SourceFiles
    flameletHRefConstThermoI.H
    flameletHRefConstThermo.C

\*---------------------------------------------------------------------------*/

#ifndef flameletHRefConstThermo_H
#define flameletHRefConstThermo_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class EquationOfState> class flameletHRefConstThermo;

template<class EquationOfState>
inline flameletHRefConstThermo<EquationOfState> operator+
(
    const flameletHRefConstThermo<EquationOfState>&,
    const flameletHRefConstThermo<EquationOfState>&
);

template<class EquationOfState>
inline flameletHRefConstThermo<EquationOfState> operator*
(
    const scalar,
    const flameletHRefConstThermo<EquationOfState>&
);

template<class EquationOfState>
inline flameletHRefConstThermo<EquationOfState> operator==
(
    const flameletHRefConstThermo<EquationOfState>&,
    const flameletHRefConstThermo<EquationOfState>&
);

template<class EquationOfState>
Ostream& operator<<
(
    Ostream&,
    const flameletHRefConstThermo<EquationOfState>&
);


/*---------------------------------------------------------------------------*\
                           Class flameletHRefConstThermo Declaration
\*---------------------------------------------------------------------------*/

template<class EquationOfState>
class flameletHRefConstThermo
:
    public EquationOfState
{
    // Private Data

        scalar Cp_;
        scalar Hf_;
        scalar Tref_;
        scalar Href_;


    // Private Member Functions

        //- Construct from components
        inline flameletHRefConstThermo
        (
            const EquationOfState& st,
            const scalar cp,
            const scalar hf,
            const scalar tref,
            const scalar href
        );


public:

    // Constructors

        //- Construct from dictionary
        flameletHRefConstThermo(const dictionary& dict);

        //- Construct as named copy
        inline flameletHRefConstThermo(const word&, const flameletHRefConstThermo&);

        //- Construct and return a clone
        inline autoPtr<flameletHRefConstThermo> clone() const;

        //- Selector from dictionary
        inline static autoPtr<flameletHRefConstThermo> New(const dictionary& dict);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "hRefConst<" + EquationOfState::typeName() + '>';
        }

        //- Limit the temperature to be in the range Tlow_ to Thigh_
        inline scalar limit(const scalar T) const;


        // Fundamental properties

            //- Heat capacity at constant pressure [J/kg/K]
            inline scalar Cp(const scalar p, const scalar T) const;

            //- Absolute Enthalpy [J/kg]
            inline scalar Ha(const scalar p, const scalar T) const;

            //- Sensible enthalpy [J/kg]
            inline scalar Hs(const scalar p, const scalar T) const;

            //- Chemical enthalpy [J/kg]
            inline scalar Hc() const;

            //- Entropy [J/kg/K]
            inline scalar S(const scalar p, const scalar T) const;

            #include "flameletHtoEthermo.H"


        // Derivative term used for Jacobian

            //- Derivative of Gibbs free energy w.r.t. temperature
            inline scalar dGdT(const scalar p, const scalar T) const;

            //- Temperature derivative of heat capacity at constant pressure
            inline scalar dCpdT(const scalar p, const scalar T) const;


        // I-O

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member Operators

        inline void operator+=(const flameletHRefConstThermo&);


    // Friend operators

        friend flameletHRefConstThermo operator+ <EquationOfState>
        (
            const flameletHRefConstThermo&,
            const flameletHRefConstThermo&
        );

        friend flameletHRefConstThermo operator* <EquationOfState>
        (
            const scalar,
            const flameletHRefConstThermo&
        );

        friend flameletHRefConstThermo operator== <EquationOfState>
        (
            const flameletHRefConstThermo&,
            const flameletHRefConstThermo&
        );


    // IOstream Operators

        friend Ostream& operator<< <EquationOfState>
        (
            Ostream&,
            const flameletHRefConstThermo&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "flameletHRefConstThermoI.H"

#ifdef NoRepository
    #include "flameletHRefConstThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
