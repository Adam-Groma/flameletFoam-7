/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2019 OpenFOAM Foundation
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
    Foam::flameletWLFTransport

Description
    Transport package using the Williams-Landel-Ferry model.

    Templated into a given thermodynamics package (needed for thermal
    conductivity).

    Dynamic viscosity [kg/m.s]
    \f[
        \mu = \mu_0 \exp \left(\frac{-C_1 ( T - T_r )}{C_2 + T - T_r}\right)
    \f]

    References:
    \verbatim
        Williams, M. L., Landel, R. F., & Ferry, J. D. (1955).
        The temperature dependence of relaxation mechanisms
        in amorphous polymers and other glass-forming liquids.
        Journal of the American Chemical society, 77(14), 3701-3707.
    \endverbatim

SourceFiles
    flameletWLFTransportI.H
    flameletWLFTransport.C

\*---------------------------------------------------------------------------*/

#ifndef flameletWLFTransport_H
#define flameletWLFTransport_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Thermo> class flameletWLFTransport;

template<class Thermo>
inline flameletWLFTransport<Thermo> operator+
(
    const flameletWLFTransport<Thermo>&,
    const flameletWLFTransport<Thermo>&
);

template<class Thermo>
inline flameletWLFTransport<Thermo> operator*
(
    const scalar,
    const flameletWLFTransport<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const flameletWLFTransport<Thermo>&
);


/*---------------------------------------------------------------------------*\
                    Class flameletWLFTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class flameletWLFTransport
:
    public Thermo
{
    // Private Data

        //- Dynamic viscosity at the reference temperature [Pa.s]
        scalar mu0_;

        //- Reference temperature [T]
        scalar Tr_;

        //- WLF coefficient 1 []
        scalar C1_;

        //- WLF coefficient 2 [T]
        scalar C2_;

        //- Reciprocal Prandtl Number []
        scalar rPr_;


    // Private Member Functions

        //- Read coefficient from dictionary
        scalar readCoeff(const word& coeffName, const dictionary& dict);


public:

    // Constructors

        //- Construct as named copy
        inline flameletWLFTransport(const word&, const flameletWLFTransport&);

        //- Construct from dictionary
        flameletWLFTransport(const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<flameletWLFTransport> clone() const;

        // Selector from dictionary
        inline static autoPtr<flameletWLFTransport> New(const dictionary& dict);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "WLF<" + Thermo::typeName() + '>';
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

        inline void operator+=(const flameletWLFTransport&);
        inline void operator*=(const scalar);


    // Friend operators

        friend flameletWLFTransport operator+ <Thermo>
        (
            const flameletWLFTransport&,
            const flameletWLFTransport&
        );

        friend flameletWLFTransport operator* <Thermo>
        (
            const scalar,
            const flameletWLFTransport&
        );


    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const flameletWLFTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "flameletWLFTransportI.H"

#ifdef NoRepository
    #include "flameletWLFTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
