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

#include "flameletSpecie.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo, int PolySize>
inline Foam::flameletPolynomialTransport<Thermo, PolySize>::flameletPolynomialTransport
(
    const Thermo& t,
    const Polynomial<PolySize>& muCoeffs,
    const Polynomial<PolySize>& kappaCoeffs
)
:
    Thermo(t),
    muCoeffs_(muCoeffs),
    kappaCoeffs_(kappaCoeffs)
{}


template<class Thermo, int PolySize>
inline Foam::flameletPolynomialTransport<Thermo, PolySize>::flameletPolynomialTransport
(
    const word& name,
    const flameletPolynomialTransport& pt
)
:
    Thermo(name, pt),
    muCoeffs_(pt.muCoeffs_),
    kappaCoeffs_(pt.kappaCoeffs_)
{}


template<class Thermo, int PolySize>
inline Foam::autoPtr<Foam::flameletPolynomialTransport<Thermo, PolySize>>
Foam::flameletPolynomialTransport<Thermo, PolySize>::clone() const
{
    return autoPtr<flameletPolynomialTransport<Thermo, PolySize>>
    (
        new flameletPolynomialTransport<Thermo, PolySize>(*this)
    );
}


template<class Thermo, int PolySize>
inline Foam::autoPtr<Foam::flameletPolynomialTransport<Thermo, PolySize>>
Foam::flameletPolynomialTransport<Thermo, PolySize>::New(const dictionary& dict)
{
    return autoPtr<flameletPolynomialTransport<Thermo, PolySize>>
    (
        new flameletPolynomialTransport<Thermo, PolySize>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo, int PolySize>
inline Foam::scalar Foam::flameletPolynomialTransport<Thermo, PolySize>::mu
(
    const scalar p,
    const scalar T
) const
{
    return muCoeffs_.value(T);
}


template<class Thermo, int PolySize>
inline Foam::scalar Foam::flameletPolynomialTransport<Thermo, PolySize>::kappa
(
    const scalar p,
    const scalar T
) const
{
    return kappaCoeffs_.value(T);
}


template<class Thermo, int PolySize>
inline Foam::scalar Foam::flameletPolynomialTransport<Thermo, PolySize>::alphah
(
    const scalar p, const scalar T
) const
{
    return kappa(p, T)/this->Cp(p, T);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo, int PolySize>
inline void Foam::flameletPolynomialTransport<Thermo, PolySize>::operator+=
(
    const flameletPolynomialTransport<Thermo, PolySize>& pt
)
{
    scalar Y1 = this->Y();

    Thermo::operator+=(pt);

    if (mag(this->Y()) > small)
    {
        Y1 /= this->Y();
        scalar Y2 = pt.Y()/this->Y();

        muCoeffs_ = Y1*muCoeffs_ + Y2*pt.muCoeffs_;
        kappaCoeffs_ = Y1*kappaCoeffs_ + Y2*pt.kappaCoeffs_;
    }
}


template<class Thermo, int PolySize>
inline void Foam::flameletPolynomialTransport<Thermo, PolySize>::operator*=
(
    const scalar s
)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo, int PolySize>
inline Foam::flameletPolynomialTransport<Thermo, PolySize> Foam::operator+
(
    const flameletPolynomialTransport<Thermo, PolySize>& pt1,
    const flameletPolynomialTransport<Thermo, PolySize>& pt2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(pt1) + static_cast<const Thermo&>(pt2)
    );

    if (mag(t.Y()) < small)
    {
        return flameletPolynomialTransport<Thermo>
        (
            t,
            0,
            pt1.muCoeffs_,
            pt1.kappaCoeffs_
        );
    }
    else
    {
        scalar Y1 = pt1.Y()/t.Y();
        scalar Y2 = pt2.Y()/t.Y();

        return flameletPolynomialTransport<Thermo, PolySize>
        (
            t,
            Y1*pt1.muCoeffs_ + Y2*pt2.muCoeffs_,
            Y1*pt1.kappaCoeffs_ + Y2*pt2.kappaCoeffs_
        );
    }
}


template<class Thermo, int PolySize>
inline Foam::flameletPolynomialTransport<Thermo, PolySize> Foam::operator*
(
    const scalar s,
    const flameletPolynomialTransport<Thermo, PolySize>& pt
)
{
    return flameletPolynomialTransport<Thermo, PolySize>
    (
        s*static_cast<const Thermo&>(pt),
        pt.muCoeffs_,
        pt.kappaCoeffs_
    );
}


// ************************************************************************* //
