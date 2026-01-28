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
    Foam::flameletHPolynomialThermo

Description
    Thermodynamics package templated on the equation of state, using polynomial
    functions for \c cp, \c h and \c s.

    Polynomials for \c h and \c s derived from \c cp.

Usage

    \table
        Property     | Description
        Hf           | Heat of formation
        Sf           | Standard entropy
        CpCoeffs<8>  | Specific heat at constant pressure polynomial coeffs
    \endtable

    Example of the specification of the thermodynamic properties:
    \verbatim
    thermodynamics
    {
        Hf              0;
        Sf              0;
        CpCoeffs<8>     ( 1000 -0.05 0.003 0 0 0 0 0 );
    }
    \endverbatim

    The polynomial expression is evaluated as so:

        \f[
            Cp = 1000 - 0.05 T + 0.003 T^2
        \f]

Note
    - Heat of formation is inputted in [J/kg], but internally uses [J/kmol]
    - Standard entropy is inputted in [J/kg/K], but internally uses [J/kmol/K]
    - Specific heat at constant pressure polynomial coefficients evaluate to an
      expression in [J/kg/K].

SourceFiles
    flameletHPolynomialThermoI.H
    flameletHPolynomialThermo.C

See also
    Foam::Polynomial

\*---------------------------------------------------------------------------*/

#ifndef flameletHPolynomialThermo_H
#define flameletHPolynomialThermo_H

#include "scalar.H"
#include "Polynomial.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class EquationOfState, int PolySize>
class flameletHPolynomialThermo;

template<class EquationOfState, int PolySize>
inline flameletHPolynomialThermo<EquationOfState, PolySize> operator+
(
    const flameletHPolynomialThermo<EquationOfState, PolySize>&,
    const flameletHPolynomialThermo<EquationOfState, PolySize>&
);

template<class EquationOfState, int PolySize>
inline flameletHPolynomialThermo<EquationOfState, PolySize> operator*
(
    const scalar,
    const flameletHPolynomialThermo<EquationOfState, PolySize>&
);

template<class EquationOfState, int PolySize>
inline flameletHPolynomialThermo<EquationOfState, PolySize> operator==
(
    const flameletHPolynomialThermo<EquationOfState, PolySize>&,
    const flameletHPolynomialThermo<EquationOfState, PolySize>&
);

template<class EquationOfState, int PolySize>
Ostream& operator<<
(
    Ostream&,
    const flameletHPolynomialThermo<EquationOfState, PolySize>&
);


/*---------------------------------------------------------------------------*\
                      Class flameletHPolynomialThermo Declaration
\*---------------------------------------------------------------------------*/

template<class EquationOfState, int PolySize=8>
class flameletHPolynomialThermo
:
    public EquationOfState
{
    // Private Data

        //- Heat of formation
        scalar Hf_;

        //- Standard entropy
        scalar Sf_;

        //- Specific heat at constant pressure polynomial coeffs
        Polynomial<PolySize> CpCoeffs_;

        //- Enthalpy polynomial coeffs - derived from cp [J/kg]
        //  NOTE: relative to Tstd
        typename Polynomial<PolySize>::intPolyType hCoeffs_;

        //- Entropy - derived from Cp [J/kg/K] - relative to Tstd
        Polynomial<PolySize> sCoeffs_;


    // Private Member Functions

        //- Construct from components
        inline flameletHPolynomialThermo
        (
            const EquationOfState& pt,
            const scalar Hf,
            const scalar Sf,
            const Polynomial<PolySize>& CpCoeffs,
            const typename Polynomial<PolySize>::intPolyType& hCoeffs,
            const Polynomial<PolySize>& sCoeffs
        );


public:

    // Constructors

        //- Construct from dictionary
        flameletHPolynomialThermo(const dictionary& dict);

        //- Construct as a named copy
        inline flameletHPolynomialThermo(const word&, const flameletHPolynomialThermo&);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "hPolynomial<" + EquationOfState::typeName() + '>';
        }

        //- Limit the temperature to be in the range Tlow_ to Thigh_
        inline scalar limit(const scalar) const;


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

        inline void operator+=(const flameletHPolynomialThermo&);
        inline void operator*=(const scalar);


    // Friend operators

        friend flameletHPolynomialThermo operator+ <EquationOfState, PolySize>
        (
            const flameletHPolynomialThermo&,
            const flameletHPolynomialThermo&
        );

        friend flameletHPolynomialThermo operator* <EquationOfState, PolySize>
        (
            const scalar,
            const flameletHPolynomialThermo&
        );

        friend flameletHPolynomialThermo operator== <EquationOfState, PolySize>
        (
            const flameletHPolynomialThermo&,
            const flameletHPolynomialThermo&
        );


    // Ostream Operator

        friend Ostream& operator<< <EquationOfState, PolySize>
        (
            Ostream&,
            const flameletHPolynomialThermo&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "flameletHPolynomialThermoI.H"

#ifdef NoRepository
    #include "flameletHPolynomialThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
