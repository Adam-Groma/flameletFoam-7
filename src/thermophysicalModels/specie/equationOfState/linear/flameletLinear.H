/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2019 OpenFOAM Foundation
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
    Foam::flameletLinear

Description
    Linear equation of state with constant compressibility

    \verbatim
        rho = rho0 + psi*p
    \endverbatim

SourceFiles
    flameletLinearI.H
    flameletLinear.C

\*---------------------------------------------------------------------------*/

#ifndef flameletLinear_H
#define flameletLinear_H

#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie> class flameletLinear;

template<class Specie>
inline flameletLinear<Specie> operator+
(
    const flameletLinear<Specie>&,
    const flameletLinear<Specie>&
);

template<class Specie>
inline flameletLinear<Specie> operator*
(
    const scalar,
    const flameletLinear<Specie>&
);

template<class Specie>
inline flameletLinear<Specie> operator==
(
    const flameletLinear<Specie>&,
    const flameletLinear<Specie>&
);

template<class Specie>
Ostream& operator<<
(
    Ostream&,
    const flameletLinear<Specie>&
);


/*---------------------------------------------------------------------------*\
                           Class flameletLinear Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class flameletLinear
:
    public Specie
{
    // Private Data

        //- Compressibility
        scalar psi_;

        //- The reference density
        scalar rho0_;


public:

    // Constructors

        //- Construct from components
        inline flameletLinear
        (
            const Specie& sp,
            const scalar psi,
            const scalar rho0
        );

        //- Construct from dictionary
        flameletLinear(const dictionary& dict);

        //- Construct as named copy
        inline flameletLinear(const word& name, const flameletLinear&);

        //- Construct and return a clone
        inline autoPtr<flameletLinear> clone() const;

        // Selector from dictionary
        inline static autoPtr<flameletLinear> New(const dictionary& dict);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "flameletLinear<" + word(Specie::typeName_()) + '>';
        }


        // Fundamental properties

            //- Is the equation of state is incompressible i.e. rho != f(p)
            static const bool incompressible = false;

            //- Is the equation of state is isochoric i.e. rho = const
            static const bool isochoric = false;

            //- Return density [kg/m^3]
            inline scalar rho(scalar p, scalar T) const;

            //- Return enthalpy departure [J/kg]
            inline scalar H(const scalar p, const scalar T) const;

            //- Return Cp departure [J/(kg K]
            inline scalar Cp(scalar p, scalar T) const;

            //- Return internal energy departure [J/kg]
            inline scalar E(const scalar p, const scalar T) const;

            //- Return Cv departure [J/(kg K]
            inline scalar Cv(scalar p, scalar T) const;

            //- Return entropy [J/kg/K]
            inline scalar S(const scalar p, const scalar T) const;

            //- Return compressibility rho/p [s^2/m^2]
            inline scalar psi(scalar p, scalar T) const;

            //- Return compression factor []
            inline scalar Z(scalar p, scalar T) const;

            //- Return (Cp - Cv) [J/(kg K]
            inline scalar CpMCv(scalar p, scalar T) const;


        // IO

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member Operators

        inline void operator+=(const flameletLinear&);
        inline void operator*=(const scalar);


    // Friend operators

        friend flameletLinear operator+ <Specie>
        (
            const flameletLinear&,
            const flameletLinear&
        );

        friend flameletLinear operator* <Specie>
        (
            const scalar s,
            const flameletLinear&
        );

        friend flameletLinear operator== <Specie>
        (
            const flameletLinear&,
            const flameletLinear&
        );


    // Ostream Operator

        friend Ostream& operator<< <Specie>
        (
            Ostream&,
            const flameletLinear&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "flameletLinearI.H"

#ifdef NoRepository
    #include "flameletLinear.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
