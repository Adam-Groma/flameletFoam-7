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
    Foam::flameletIncompressiblePerfectGas

Description
    Incompressible gas equation of state using a constant reference pressure in
    the perfect gas equation of state rather than the local pressure so that the
    density only varies with temperature and composition.

SourceFiles
    flameletIncompressiblePerfectGasI.H
    flameletIncompressiblePerfectGas.C

\*---------------------------------------------------------------------------*/

#ifndef flameletIncompressiblePerfectGas_H
#define flameletIncompressiblePerfectGas_H

#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie> class flameletIncompressiblePerfectGas;

template<class Specie>
inline flameletIncompressiblePerfectGas<Specie> operator+
(
    const flameletIncompressiblePerfectGas<Specie>&,
    const flameletIncompressiblePerfectGas<Specie>&
);

template<class Specie>
inline flameletIncompressiblePerfectGas<Specie> operator*
(
    const scalar,
    const flameletIncompressiblePerfectGas<Specie>&
);

template<class Specie>
inline flameletIncompressiblePerfectGas<Specie> operator==
(
    const flameletIncompressiblePerfectGas<Specie>&,
    const flameletIncompressiblePerfectGas<Specie>&
);

template<class Specie>
Ostream& operator<<
(
    Ostream&,
    const flameletIncompressiblePerfectGas<Specie>&
);


/*---------------------------------------------------------------------------*\
                    Class flameletIncompressiblePerfectGas Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class flameletIncompressiblePerfectGas
:
    public Specie
{
    // Private Data

        //- Reference pressure
        scalar pRef_;


public:

    // Constructors

        //- Construct from components
        inline flameletIncompressiblePerfectGas(const Specie& sp, const scalar pRef);

        //- Construct from dictionary
        flameletIncompressiblePerfectGas(const dictionary& dict);

        //- Construct as named copy
        inline flameletIncompressiblePerfectGas
        (
            const word& name,
            const flameletIncompressiblePerfectGas&
        );

        //- Construct and return a clone
        inline autoPtr<flameletIncompressiblePerfectGas> clone() const;

        // Selector from dictionary
        inline static autoPtr<flameletIncompressiblePerfectGas> New
        (
            const dictionary& dict
        );


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return
                "flameletIncompressiblePerfectGas<"
              + word(Specie::typeName_()) + '>';
        }


        // Fundamental properties

            //- Is the equation of state is incompressible i.e. rho != f(p)
            static const bool incompressible = true;

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

        inline void operator+=(const flameletIncompressiblePerfectGas&);
        inline void operator*=(const scalar);


    // Friend operators

        friend flameletIncompressiblePerfectGas operator+ <Specie>
        (
            const flameletIncompressiblePerfectGas&,
            const flameletIncompressiblePerfectGas&
        );

        friend flameletIncompressiblePerfectGas operator* <Specie>
        (
            const scalar s,
            const flameletIncompressiblePerfectGas&
        );

        friend flameletIncompressiblePerfectGas operator== <Specie>
        (
            const flameletIncompressiblePerfectGas&,
            const flameletIncompressiblePerfectGas&
        );


    // Ostream Operator

        friend Ostream& operator<< <Specie>
        (
            Ostream&,
            const flameletIncompressiblePerfectGas&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "flameletIncompressiblePerfectGasI.H"

#ifdef NoRepository
    #include "flameletIncompressiblePerfectGas.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
