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

\*---------------------------------------------------------------------------*/

#include "flameletIcoPolynomial.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Specie, int PolySize>
inline Foam::flameletIcoPolynomial<Specie, PolySize>::flameletIcoPolynomial
(
    const Specie& sp,
    const Polynomial<PolySize>& rhoCoeffs
)
:
    Specie(sp),
    rhoCoeffs_(rhoCoeffs)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie, int PolySize>
inline Foam::flameletIcoPolynomial<Specie, PolySize>::flameletIcoPolynomial
(
    const word& name,
    const flameletIcoPolynomial<Specie, PolySize>& ip
)
:
    Specie(name, ip),
    rhoCoeffs_(ip.rhoCoeffs_)
{}


template<class Specie, int PolySize>
inline Foam::autoPtr<Foam::flameletIcoPolynomial<Specie, PolySize>>
Foam::flameletIcoPolynomial<Specie, PolySize>::clone() const
{
    return autoPtr<flameletIcoPolynomial<Specie, PolySize>>
    (
        new flameletIcoPolynomial<Specie, PolySize>(*this)
    );
}


template<class Specie, int PolySize>
inline Foam::autoPtr<Foam::flameletIcoPolynomial<Specie, PolySize>>
Foam::flameletIcoPolynomial<Specie, PolySize>::New(const dictionary& dict)
{
    return autoPtr<flameletIcoPolynomial<Specie, PolySize>>
    (
        new flameletIcoPolynomial<Specie, PolySize>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie, int PolySize>
inline Foam::scalar Foam::flameletIcoPolynomial<Specie, PolySize>::rho
(
    scalar p,
    scalar T
) const
{
    return rhoCoeffs_.value(T);
}


template<class Specie, int PolySize>
inline Foam::scalar Foam::flameletIcoPolynomial<Specie, PolySize>::H
(
    scalar p,
    scalar T
) const
{
    return p/this->rho(p, T);
}


template<class Specie, int PolySize>
inline Foam::scalar Foam::flameletIcoPolynomial<Specie, PolySize>::Cp
(
    scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie, int PolySize>
inline Foam::scalar Foam::flameletIcoPolynomial<Specie, PolySize>::E
(
    scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie, int PolySize>
inline Foam::scalar Foam::flameletIcoPolynomial<Specie, PolySize>::Cv
(
    scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie, int PolySize>
inline Foam::scalar Foam::flameletIcoPolynomial<Specie, PolySize>::S
(
    scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie, int PolySize>
inline Foam::scalar Foam::flameletIcoPolynomial<Specie, PolySize>::psi
(
    scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie, int PolySize>
inline Foam::scalar Foam::flameletIcoPolynomial<Specie, PolySize>::Z
(
    scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie, int PolySize>
inline Foam::scalar Foam::flameletIcoPolynomial<Specie, PolySize>::CpMCv
(
    scalar p,
    scalar T
) const
{
    return 0;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie, int PolySize>
inline void Foam::flameletIcoPolynomial<Specie, PolySize>::operator+=
(
    const flameletIcoPolynomial<Specie, PolySize>& ip
)
{
    scalar Y1 = this->Y();
    Specie::operator+=(ip);

    if (mag(this->Y()) > small)
    {
        Y1 /= this->Y();
        const scalar Y2 = ip.Y()/this->Y();

        rhoCoeffs_ = Y1*rhoCoeffs_ + Y2*ip.rhoCoeffs_;
    }
}


template<class Specie, int PolySize>
inline void Foam::flameletIcoPolynomial<Specie, PolySize>::operator*=(const scalar s)
{
    Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Specie, int PolySize>
Foam::flameletIcoPolynomial<Specie, PolySize> Foam::operator+
(
    const flameletIcoPolynomial<Specie, PolySize>& ip1,
    const flameletIcoPolynomial<Specie, PolySize>& ip2
)
{
    Specie sp
    (
        static_cast<const Specie&>(ip1)
      + static_cast<const Specie&>(ip2)
    );

    if (mag(sp.Y()) < small)
    {
        return flameletIcoPolynomial<Specie, PolySize>
        (
            sp,
            ip1.rhoCoeffs_
        );
    }
    else
    {
        const scalar Y1 = ip1.Y()/sp.Y();
        const scalar Y2 = ip2.Y()/sp.Y();

        return flameletIcoPolynomial<Specie, PolySize>
        (
            sp,
            Y1*ip1.rhoCoeffs_ + Y2*ip2.rhoCoeffs_
        );
    }
}


template<class Specie, int PolySize>
Foam::flameletIcoPolynomial<Specie, PolySize> Foam::operator*
(
    const scalar s,
    const flameletIcoPolynomial<Specie, PolySize>& ip
)
{
    return flameletIcoPolynomial<Specie, PolySize>
    (
        s*static_cast<const Specie&>(ip),
        ip.rhoCoeffs_
    );
}


template<class Specie, int PolySize>
Foam::flameletIcoPolynomial<Specie, PolySize> Foam::operator==
(
    const flameletIcoPolynomial<Specie, PolySize>& ip1,
    const flameletIcoPolynomial<Specie, PolySize>& ip2
)
{
    Specie sp
    (
        static_cast<const Specie&>(ip1)
     == static_cast<const Specie&>(ip2)
    );

    const scalar Y1 = ip1.Y()/sp.Y();
    const scalar Y2 = ip2.Y()/sp.Y();

    return flameletIcoPolynomial<Specie, PolySize>
    (
        sp,
        Y2*ip2.rhoCoeffs_ - Y1*ip1.rhoCoeffs_
    );
}


// ************************************************************************* //
