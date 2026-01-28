/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "flameletPerfectGas.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Specie>
inline Foam::flameletPerfectGas<Specie>::flameletPerfectGas(const Specie& sp)
:
    Specie(sp)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::flameletPerfectGas<Specie>::flameletPerfectGas
(
    const word& name,
    const flameletPerfectGas<Specie>& pg
)
:
    Specie(name, pg)
{}


template<class Specie>
inline Foam::autoPtr<Foam::flameletPerfectGas<Specie>>
Foam::flameletPerfectGas<Specie>::clone() const
{
    return autoPtr<flameletPerfectGas<Specie>>(new flameletPerfectGas<Specie>(*this));
}


template<class Specie>
inline Foam::autoPtr<Foam::flameletPerfectGas<Specie>> Foam::flameletPerfectGas<Specie>::New
(
    const dictionary& dict
)
{
    return autoPtr<flameletPerfectGas<Specie>>(new flameletPerfectGas<Specie>(dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::scalar Foam::flameletPerfectGas<Specie>::rho(scalar p, scalar T) const
{
    return p/(this->R()*T);
}


template<class Specie>
inline Foam::scalar Foam::flameletPerfectGas<Specie>::H(scalar p, scalar T) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::flameletPerfectGas<Specie>::Cp(scalar p, scalar T) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::flameletPerfectGas<Specie>::E(scalar p, scalar T) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::flameletPerfectGas<Specie>::Cv(scalar p, scalar T) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::flameletPerfectGas<Specie>::S(scalar p, scalar T) const
{
    return -this->R()*log(p/Pstd);
}


template<class Specie>
inline Foam::scalar Foam::flameletPerfectGas<Specie>::psi(scalar p, scalar T) const
{
    return 1.0/(this->R()*T);
}


template<class Specie>
inline Foam::scalar Foam::flameletPerfectGas<Specie>::Z(scalar p, scalar T) const
{
    return 1;
}


template<class Specie>
inline Foam::scalar Foam::flameletPerfectGas<Specie>::CpMCv(scalar p, scalar T) const
{
    return this->R();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie>
inline void Foam::flameletPerfectGas<Specie>::operator+=(const flameletPerfectGas<Specie>& pg)
{
    Specie::operator+=(pg);
}


template<class Specie>
inline void Foam::flameletPerfectGas<Specie>::operator*=(const scalar s)
{
    Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::flameletPerfectGas<Specie> Foam::operator+
(
    const flameletPerfectGas<Specie>& pg1,
    const flameletPerfectGas<Specie>& pg2
)
{
    return flameletPerfectGas<Specie>
    (
        static_cast<const Specie&>(pg1) + static_cast<const Specie&>(pg2)
    );
}


template<class Specie>
inline Foam::flameletPerfectGas<Specie> Foam::operator*
(
    const scalar s,
    const flameletPerfectGas<Specie>& pg
)
{
    return flameletPerfectGas<Specie>(s*static_cast<const Specie&>(pg));
}


template<class Specie>
inline Foam::flameletPerfectGas<Specie> Foam::operator==
(
    const flameletPerfectGas<Specie>& pg1,
    const flameletPerfectGas<Specie>& pg2
)
{
    return flameletPerfectGas<Specie>
    (
        static_cast<const Specie&>(pg1) == static_cast<const Specie&>(pg2)
    );
}


// ************************************************************************* //
