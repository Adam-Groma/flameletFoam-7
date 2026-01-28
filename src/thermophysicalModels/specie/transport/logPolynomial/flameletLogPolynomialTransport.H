/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2019 OpenFOAM Foundation
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
    Foam::flameletLogPolynomialTransport

Description
    Transport package using polynomial functions of \c ln(T) for \c mu and
    \c kappa:

        \f[
            ln(mu)    = \sum_{i=1}^N \left( a[i] * ln(T)^{i-1} \right)
        \f]

        \f[
            ln(kappa) = \sum_{i=1}^N \left( b[i] * ln(T)^{i-1} \right)
        \f]

Usage

    \table
        Property        | Description
        muCoeffs<8>     | Dynamic viscosity polynomial coefficients
        kappaCoeffs<8>  | Thermal conductivity polynomial coefficients
    \endtable

    Example of the specification of the transport properties:
    \verbatim
    transport
    {
        muCoeffs<8>     ( 1000 -0.05 0.003 0 0 0 0 0 );
        kappaCoeffs<8>  ( 2000 -0.15 0.023 0 0 0 0 0 );
    }
    \endverbatim

    The polynomial expressions are evaluated as so:

        \f[
            \mu    = 1000 - 0.05 ln(T) + 0.003 ln(T)^2
        \f]

        \f[
            \kappa = 2000 - 0.15 ln(T) + 0.023 ln(T)^2
        \f]

Note
    - Dynamic viscosity polynomial coefficients evaluate to an expression in
      [Pa.s], but internally uses [Pa.s/kmol].
    - Thermal conductivity polynomial coefficients evaluate to an expression in
      [W/m/K], but internally uses [W/m/K/kmol].

SourceFiles
    flameletLogPolynomialTransportI.H
    flameletLogPolynomialTransport.C

See also
    Foam::Polynomial

\*---------------------------------------------------------------------------*/

#ifndef flameletLogPolynomialTransport_H
#define flameletLogPolynomialTransport_H

#include "Polynomial.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Thermo, int PolySize> class flameletLogPolynomialTransport;

template<class Thermo, int PolySize>
inline flameletLogPolynomialTransport<Thermo, PolySize> operator+
(
    const flameletLogPolynomialTransport<Thermo, PolySize>&,
    const flameletLogPolynomialTransport<Thermo, PolySize>&
);

template<class Thermo, int PolySize>
inline flameletLogPolynomialTransport<Thermo, PolySize> operator*
(
    const scalar,
    const flameletLogPolynomialTransport<Thermo, PolySize>&
);

template<class Thermo, int PolySize>
Ostream& operator<<
(
    Ostream&,
    const flameletLogPolynomialTransport<Thermo, PolySize>&
);


/*---------------------------------------------------------------------------*\
                     Class flameletLogPolynomialTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo, int PolySize=8>
class flameletLogPolynomialTransport
:
    public Thermo
{
    // Private Data

        //- Dynamic viscosity polynomial coefficients
        //  Note: input in [Pa.s], but internally uses [Pa.s/kmol]
        Polynomial<PolySize> muCoeffs_;

        //- Thermal conductivity polynomial coefficients
        //  Note: input in [W/m/K], but internally uses [W/m/K/kmol]
        Polynomial<PolySize> kappaCoeffs_;


    // Private Member Functions

        //- Construct from components
        inline flameletLogPolynomialTransport
        (
            const Thermo& t,
            const Polynomial<PolySize>& muPoly,
            const Polynomial<PolySize>& kappaPoly
        );


public:

    // Constructors

        //- Construct as named copy
        inline flameletLogPolynomialTransport
        (
            const word&,
            const flameletLogPolynomialTransport&
        );

        //- Construct from dictionary
        flameletLogPolynomialTransport(const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<flameletLogPolynomialTransport> clone() const;

        // Selector from dictionary
        inline static autoPtr<flameletLogPolynomialTransport> New
        (
            const dictionary& dict
        );


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "logPolynomial<" + Thermo::typeName() + '>';
        }

        //- Dynamic viscosity [kg/m/s]
        inline scalar mu(const scalar p, const scalar T) const;

        //- Thermal conductivity [W/m/K]
        inline scalar kappa(const scalar p, const scalar T) const;

        //- Thermal diffusivity of enthalpy [kg/m/s]
        inline scalar alphah(const scalar p, const scalar T) const;

        // Species diffusivity
        // inline scalar D(const scalar p, const scalar T) const;

        //- Write to Ostream
        void write(Ostream& os) const;


    // Member Operators

        inline void operator+=(const flameletLogPolynomialTransport&);
        inline void operator*=(const scalar);


    // Friend operators

        friend flameletLogPolynomialTransport operator+ <Thermo, PolySize>
        (
            const flameletLogPolynomialTransport&,
            const flameletLogPolynomialTransport&
        );

        friend flameletLogPolynomialTransport operator* <Thermo, PolySize>
        (
            const scalar,
            const flameletLogPolynomialTransport&
        );


    // Ostream Operator

        friend Ostream& operator<< <Thermo, PolySize>
        (
            Ostream&,
            const flameletLogPolynomialTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "flameletLogPolynomialTransportI.H"

#ifdef NoRepository
    #include "flameletLogPolynomialTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
