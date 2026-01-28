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
    Foam::flameletSutherlandTransport

Description
    Transport package using Sutherland's formula.

    Templated into a given thermodynamics package (needed for thermal
    conductivity).

    Dynamic viscosity [kg/m.s]
    \f[
        \mu = A_s \frac{\sqrt{T}}{1 + T_s / T}
    \f]

SourceFiles
    flameletSutherlandTransportI.H
    flameletSutherlandTransport.C

\*---------------------------------------------------------------------------*/

#ifndef flameletSutherlandTransport_H
#define flameletSutherlandTransport_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Thermo> class flameletSutherlandTransport;

template<class Thermo>
inline flameletSutherlandTransport<Thermo> operator+
(
    const flameletSutherlandTransport<Thermo>&,
    const flameletSutherlandTransport<Thermo>&
);

template<class Thermo>
inline flameletSutherlandTransport<Thermo> operator*
(
    const scalar,
    const flameletSutherlandTransport<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const flameletSutherlandTransport<Thermo>&
);


/*---------------------------------------------------------------------------*\
                    Class flameletSutherlandTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class flameletSutherlandTransport
:
    public Thermo
{
    // Private Data

        // Sutherland's coefficients
        scalar As_, Ts_;


    // Private Member Functions

        //- Calculate the Sutherland coefficients
        //  given two viscosities and temperatures
        inline void calcCoeffs
        (
            const scalar mu1, const scalar T1,
            const scalar mu2, const scalar T2
        );

        //- Read coefficient from dictionary
        scalar readCoeff(const word& coeffName, const dictionary& dict);


public:

    // Constructors

        //- Construct from components
        inline flameletSutherlandTransport
        (
            const Thermo& t,
            const scalar As,
            const scalar Ts
        );

        //- Construct from two viscosities
        inline flameletSutherlandTransport
        (
            const Thermo& t,
            const scalar mu1, const scalar T1,
            const scalar mu2, const scalar T2
        );

        //- Construct as named copy
        inline flameletSutherlandTransport(const word&, const flameletSutherlandTransport&);

        //- Construct from dictionary
        flameletSutherlandTransport(const dictionary& dict);

        //- Construct from base thermo and dictionary
        flameletSutherlandTransport(const Thermo& t,const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<flameletSutherlandTransport> clone() const;

        // Selector from dictionary
        inline static autoPtr<flameletSutherlandTransport> New(const dictionary& dict);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "sutherland<" + Thermo::typeName() + '>';
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

        inline void operator+=(const flameletSutherlandTransport&);
        inline void operator*=(const scalar);


    // Friend operators

        friend flameletSutherlandTransport operator+ <Thermo>
        (
            const flameletSutherlandTransport&,
            const flameletSutherlandTransport&
        );

        friend flameletSutherlandTransport operator* <Thermo>
        (
            const scalar,
            const flameletSutherlandTransport&
        );


    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const flameletSutherlandTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "flameletSutherlandTransportI.H"

#ifdef NoRepository
    #include "flameletSutherlandTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
