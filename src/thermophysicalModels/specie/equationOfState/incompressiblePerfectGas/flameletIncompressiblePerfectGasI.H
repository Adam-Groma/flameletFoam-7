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

\*---------------------------------------------------------------------------*/

#include "flameletIncompressiblePerfectGas.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::flameletIncompressiblePerfectGas<Specie>::flameletIncompressiblePerfectGas
(
    const Specie& sp,  const scalar pRef
)
:
    Specie(sp),
    pRef_(pRef)
{}


template<class Specie>
inline Foam::flameletIncompressiblePerfectGas<Specie>::flameletIncompressiblePerfectGas
(
    const word& name,
    const flameletIncompressiblePerfectGas<Specie>& ipg
)
:
    Specie(name, ipg),
    pRef_(ipg.pRef_)
{}


template<class Specie>
inline Foam::autoPtr<Foam::flameletIncompressiblePerfectGas<Specie>>
Foam::flameletIncompressiblePerfectGas<Specie>::clone() const
{
    return autoPtr<flameletIncompressiblePerfectGas<Specie>>
    (
        new flameletIncompressiblePerfectGas<Specie>(*this)
    );
}


template<class Specie>
inline Foam::autoPtr<Foam::flameletIncompressiblePerfectGas<Specie>>
Foam::flameletIncompressiblePerfectGas<Specie>::New
(
    const dictionary& dict
)
{
    return autoPtr<flameletIncompressiblePerfectGas<Specie>>
    (
        new flameletIncompressiblePerfectGas<Specie>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::scalar Foam::flameletIncompressiblePerfectGas<Specie>::rho
(
    scalar p,
    scalar T
) const
{
    return pRef_/(this->R()*T);
}


template<class Specie>
inline Foam::scalar Foam::flameletIncompressiblePerfectGas<Specie>::H
(
    scalar p,
    scalar T
) const
{
    return p/this->rho(p, T);
}


template<class Specie>
inline Foam::scalar Foam::flameletIncompressiblePerfectGas<Specie>::Cp
(
    scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::flameletIncompressiblePerfectGas<Specie>::E
(
    scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::flameletIncompressiblePerfectGas<Specie>::Cv
(
    scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::flameletIncompressiblePerfectGas<Specie>::S
(
    scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::flameletIncompressiblePerfectGas<Specie>::psi
(
    scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::flameletIncompressiblePerfectGas<Specie>::Z
(
    scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::flameletIncompressiblePerfectGas<Specie>::CpMCv
(
    scalar p,
    scalar T
) const
{
    return 0;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie>
inline void Foam::flameletIncompressiblePerfectGas<Specie>::operator+=
(
    const flameletIncompressiblePerfectGas<Specie>& ipg
)
{
    scalar Y1 = this->Y();
    Specie::operator+=(ipg);

    if (mag(this->Y()) > small)
    {
        Y1 /= this->Y();
        const scalar Y2 = ipg.Y()/this->Y();

        pRef_ = Y1*pRef_ + Y2*ipg.pRef_;
    }
}


template<class Specie>
inline void Foam::flameletIncompressiblePerfectGas<Specie>::operator*=(const scalar s)
{
    Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::flameletIncompressiblePerfectGas<Specie> Foam::operator+
(
    const flameletIncompressiblePerfectGas<Specie>& ipg1,
    const flameletIncompressiblePerfectGas<Specie>& ipg2
)
{
    Specie sp
    (
        static_cast<const Specie&>(ipg1)
      + static_cast<const Specie&>(ipg2)
    );

    if (mag(sp.Y()) < small)
    {
        return flameletIncompressiblePerfectGas<Specie>
        (
            sp,
            ipg1.pRef_
        );
    }
    else
    {
        const scalar Y1 = ipg1.Y()/sp.Y();
        const scalar Y2 = ipg2.Y()/sp.Y();

        return flameletIncompressiblePerfectGas<Specie>
        (
            sp,
            Y1*ipg1.pRef_ + Y2*ipg2.pRef_
        );
    }
}


template<class Specie>
inline Foam::flameletIncompressiblePerfectGas<Specie> Foam::operator*
(
    const scalar s,
    const flameletIncompressiblePerfectGas<Specie>& ipg
)
{
    return flameletIncompressiblePerfectGas<Specie>
    (
        s*static_cast<const Specie&>(ipg),
        ipg.pRef_
    );
}


template<class Specie>
inline Foam::flameletIncompressiblePerfectGas<Specie> Foam::operator==
(
    const flameletIncompressiblePerfectGas<Specie>& ipg1,
    const flameletIncompressiblePerfectGas<Specie>& ipg2
)
{
    Specie sp
    (
        static_cast<const Specie&>(ipg1)
     == static_cast<const Specie&>(ipg2)
    );

    const scalar Y1 = ipg1.Y()/sp.Y();
    const scalar Y2 = ipg2.Y()/sp.Y();

    return flameletIncompressiblePerfectGas<Specie>
    (
        sp,
        Y2*ipg2.pRef_ - Y1*ipg1.pRef_
    );
}


// ************************************************************************* //
