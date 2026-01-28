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
    Foam::flameletPerfectFluid

Description
    Perfect gas equation of state.

SourceFiles
    flameletPerfectFluidI.H
    flameletPerfectFluid.C

\*---------------------------------------------------------------------------*/

#ifndef flameletPerfectFluid_H
#define flameletPerfectFluid_H

#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie> class flameletPerfectFluid;

template<class Specie>
inline flameletPerfectFluid<Specie> operator+
(
    const flameletPerfectFluid<Specie>&,
    const flameletPerfectFluid<Specie>&
);

template<class Specie>
inline flameletPerfectFluid<Specie> operator*
(
    const scalar,
    const flameletPerfectFluid<Specie>&
);

template<class Specie>
inline flameletPerfectFluid<Specie> operator==
(
    const flameletPerfectFluid<Specie>&,
    const flameletPerfectFluid<Specie>&
);

template<class Specie>
Ostream& operator<<
(
    Ostream&,
    const flameletPerfectFluid<Specie>&
);


/*---------------------------------------------------------------------------*\
                           Class flameletPerfectFluid Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class flameletPerfectFluid
:
    public Specie
{
    // Private Data

        //- Fluid constant
        scalar R_;

        //- The reference density
        scalar rho0_;


public:

    // Constructors

        //- Construct from components
        inline flameletPerfectFluid
        (
            const Specie& sp,
            const scalar R,
            const scalar rho0
        );

        //- Construct from dictionary
        flameletPerfectFluid(const dictionary& dict);

        //- Construct as named copy
        inline flameletPerfectFluid(const word& name, const flameletPerfectFluid&);

        //- Construct and return a clone
        inline autoPtr<flameletPerfectFluid> clone() const;

        // Selector from dictionary
        inline static autoPtr<flameletPerfectFluid> New(const dictionary& dict);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "flameletPerfectFluid<" + word(Specie::typeName_()) + '>';
        }


        // Fundamental properties

            //- Is the equation of state is incompressible i.e. rho != f(p)
            static const bool incompressible = false;

            //- Is the equation of state is isochoric i.e. rho = const
            static const bool isochoric = false;

            //- Return fluid constant [J/kg/K]
            inline scalar R() const;

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

        inline void operator+=(const flameletPerfectFluid&);
        inline void operator*=(const scalar);


    // Friend operators

        friend flameletPerfectFluid operator+ <Specie>
        (
            const flameletPerfectFluid&,
            const flameletPerfectFluid&
        );

        friend flameletPerfectFluid operator* <Specie>
        (
            const scalar s,
            const flameletPerfectFluid&
        );

        friend flameletPerfectFluid operator== <Specie>
        (
            const flameletPerfectFluid&,
            const flameletPerfectFluid&
        );


    // Ostream Operator

        friend Ostream& operator<< <Specie>
        (
            Ostream&,
            const flameletPerfectFluid&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "flameletPerfectFluidI.H"

#ifdef NoRepository
    #include "flameletPerfectFluid.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
