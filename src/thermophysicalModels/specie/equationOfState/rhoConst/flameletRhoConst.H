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
    Foam::flameletRhoConst

Description
    RhoConst (rho = const) of state.

SourceFiles
    flameletRhoConstI.H
    flameletRhoConst.C

\*---------------------------------------------------------------------------*/

#ifndef flameletRhoConst_H
#define flameletRhoConst_H

#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie> class flameletRhoConst;

template<class Specie>
inline flameletRhoConst<Specie> operator+
(
    const flameletRhoConst<Specie>&,
    const flameletRhoConst<Specie>&
);

template<class Specie>
inline flameletRhoConst<Specie> operator*
(
    const scalar,
    const flameletRhoConst<Specie>&
);

template<class Specie>
inline flameletRhoConst<Specie> operator==
(
    const flameletRhoConst<Specie>&,
    const flameletRhoConst<Specie>&
);

template<class Specie>
Ostream& operator<<
(
    Ostream&,
    const flameletRhoConst<Specie>&
);


/*---------------------------------------------------------------------------*\
                        Class flameletRhoConst Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class flameletRhoConst
:
    public Specie
{
    // Private Data

        //- Density
        scalar rho_;


public:

    // Constructors

        //- Construct from components
        inline flameletRhoConst(const Specie& sp, const scalar rho);

        //- Construct from dictionary
        flameletRhoConst(const dictionary& dict);

        //- Construct as named copy
        inline flameletRhoConst(const word& name, const flameletRhoConst&);

        //- Construct and return a clone
        inline autoPtr<flameletRhoConst> clone() const;


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "flameletRhoConst<" + word(Specie::typeName_()) + '>';
        }


        // Fundamental properties

            //- Is the equation of state is incompressible i.e. rho != f(p)
            static const bool incompressible = true;

            //- Is the equation of state is isochoric i.e. rho = const
            static const bool isochoric = true;

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

        inline void operator+=(const flameletRhoConst&);
        inline void operator*=(const scalar);


    // Friend operators

        friend flameletRhoConst operator+ <Specie>
        (
            const flameletRhoConst&,
            const flameletRhoConst&
        );

        friend flameletRhoConst operator* <Specie>
        (
            const scalar s,
            const flameletRhoConst&
        );

        friend flameletRhoConst operator== <Specie>
        (
            const flameletRhoConst&,
            const flameletRhoConst&
        );


    // Ostream Operator

        friend Ostream& operator<< <Specie>
        (
            Ostream&,
            const flameletRhoConst&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "flameletRhoConstI.H"

#ifdef NoRepository
    #include "flameletRhoConst.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
