/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2019 OpenFOAM Foundation
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
    Foam::flameletHPowerThermo

Description
    Power-function based thermodynamics package templated on EquationOfState.

    In this thermodynamics package the heat capacity is a simple power of
    temperature:

        Cp(T) = c0*(T/Tref)^n0;

    which is particularly suitable for solids.

SourceFiles
    flameletHPowerThermoI.H
    flameletHPowerThermo.C

\*---------------------------------------------------------------------------*/

#ifndef flameletHPowerThermo_H
#define flameletHPowerThermo_H

#include "scalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class EquationOfState> class flameletHPowerThermo;

template<class EquationOfState>
inline flameletHPowerThermo<EquationOfState> operator+
(
    const flameletHPowerThermo<EquationOfState>&,
    const flameletHPowerThermo<EquationOfState>&
);

template<class EquationOfState>
inline flameletHPowerThermo<EquationOfState> operator*
(
    const scalar,
    const flameletHPowerThermo<EquationOfState>&
);


template<class EquationOfState>
inline flameletHPowerThermo<EquationOfState> operator==
(
    const flameletHPowerThermo<EquationOfState>&,
    const flameletHPowerThermo<EquationOfState>&
);


template<class EquationOfState>
Ostream& operator<<
(
    Ostream&,
    const flameletHPowerThermo<EquationOfState>&
);


/*---------------------------------------------------------------------------*\
                         Class flameletHPowerThermo Declaration
\*---------------------------------------------------------------------------*/

template<class EquationOfState>
class flameletHPowerThermo
:
    public EquationOfState
{
    // Private Data

        scalar c0_;
        scalar n0_;
        scalar Tref_;
        scalar Hf_;


    // Private Member Functions

        //- Check given temperature is within the range of the fitted coeffs
        inline void checkT(const scalar T) const;

        //- Construct from components
        inline flameletHPowerThermo
        (
            const EquationOfState& st,
            const scalar c0,
            const scalar n0,
            const scalar Tref,
            const scalar Hf
        );


public:

    // Constructors

        //- Construct from dictionary
        flameletHPowerThermo(const dictionary&);

        //- Construct as a named copy
        inline flameletHPowerThermo
        (
            const word&,
            const flameletHPowerThermo&
        );

         //- Construct and return a clone
        inline autoPtr<flameletHPowerThermo> clone() const;

        //- Selector from dictionary
        inline static autoPtr<flameletHPowerThermo> New(const dictionary& dict);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "hPower<" + EquationOfState::typeName() + '>';
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


    // Member Operators

        inline void operator+=(const flameletHPowerThermo&);


    // Friend operators

        friend flameletHPowerThermo operator+ <EquationOfState>
        (
            const flameletHPowerThermo&,
            const flameletHPowerThermo&
        );

        friend flameletHPowerThermo operator* <EquationOfState>
        (
            const scalar,
            const flameletHPowerThermo&
        );


        friend flameletHPowerThermo operator== <EquationOfState>
        (
            const flameletHPowerThermo&,
            const flameletHPowerThermo&
        );


    // Ostream Operator

        friend Ostream& operator<< <EquationOfState>
        (
            Ostream&,
            const flameletHPowerThermo&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "flameletHPowerThermoI.H"
    #include "flameletHPowerThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
