/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2019 OpenFOAM Foundation
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
    Foam::flameletPengRobinsonGas

Description
    flameletPengRobinsonGas gas equation of state.

SourceFiles
    flameletPengRobinsonGasI.H
    flameletPengRobinsonGas.C

\*---------------------------------------------------------------------------*/

#ifndef flameletPengRobinsonGas_H
#define flameletPengRobinsonGas_H

#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie> class flameletPengRobinsonGas;

template<class Specie>
inline flameletPengRobinsonGas<Specie> operator+
(
    const flameletPengRobinsonGas<Specie>&,
    const flameletPengRobinsonGas<Specie>&
);

template<class Specie>
inline flameletPengRobinsonGas<Specie> operator*
(
    const scalar,
    const flameletPengRobinsonGas<Specie>&
);

template<class Specie>
inline flameletPengRobinsonGas<Specie> operator==
(
    const flameletPengRobinsonGas<Specie>&,
    const flameletPengRobinsonGas<Specie>&
);

template<class Specie>
Ostream& operator<<
(
    Ostream&,
    const flameletPengRobinsonGas<Specie>&
);



/*---------------------------------------------------------------------------*\
                           Class flameletPengRobinsonGas Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class flameletPengRobinsonGas
:
    public Specie
{
    // Private Data

        //- Critical Temperature [K]
        scalar Tc_;

        //- Critical volume [m^3/kmol]
        scalar Vc_;

        //- Critical compression factor [-]
        scalar Zc_;

        //- Critical Pressure [Pa]
        scalar Pc_;

        //- Acentric factor [-]
        scalar omega_;


public:

    // Constructors

        //- Construct from components
        inline flameletPengRobinsonGas
        (
            const Specie& sp,
            const scalar& Tc,
            const scalar& Vc,
            const scalar& Zc,
            const scalar& Pc,
            const scalar& omega
        );

        //- Construct from dictionary
        flameletPengRobinsonGas(const dictionary& dict);

        //- Construct as named copy
        inline flameletPengRobinsonGas(const word& name, const flameletPengRobinsonGas&);

        //- Construct and return a clone
        inline autoPtr<flameletPengRobinsonGas> clone() const;

        // Selector from dictionary
        inline static autoPtr<flameletPengRobinsonGas> New
        (
            const dictionary& dict
        );


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "flameletPengRobinsonGas<" + word(Specie::typeName_()) + '>';
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

        inline void operator+=(const flameletPengRobinsonGas&);
        inline void operator*=(const scalar);


    // Friend operators

        friend flameletPengRobinsonGas operator+ <Specie>
        (
            const flameletPengRobinsonGas&,
            const flameletPengRobinsonGas&
        );

        friend flameletPengRobinsonGas operator* <Specie>
        (
            const scalar s,
            const flameletPengRobinsonGas&
        );

        friend flameletPengRobinsonGas operator== <Specie>
        (
            const flameletPengRobinsonGas&,
            const flameletPengRobinsonGas&
        );


    // Ostream Operator

        friend Ostream& operator<< <Specie>
        (
            Ostream&,
            const flameletPengRobinsonGas&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "flameletPengRobinsonGasI.H"

#ifdef NoRepository
    #include "flameletPengRobinsonGas.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
