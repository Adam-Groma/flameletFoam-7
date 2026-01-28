/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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

#include "flameletBoussinesq.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::flameletBoussinesq<Specie>::flameletBoussinesq
(
    const Specie& sp, const scalar rho0, const scalar T0, const scalar beta
)
:
    Specie(sp),
    rho0_(rho0),
    T0_(T0),
    beta_(beta)
{}


template<class Specie>
inline Foam::flameletBoussinesq<Specie>::flameletBoussinesq
(
    const word& name,
    const flameletBoussinesq<Specie>& b
)
:
    Specie(name, b),
    rho0_(b.rho0_),
    T0_(b.T0_),
    beta_(b.beta_)
{}


template<class Specie>
inline Foam::autoPtr<Foam::flameletBoussinesq<Specie>>
Foam::flameletBoussinesq<Specie>::clone() const
{
    return autoPtr<flameletBoussinesq<Specie>>
    (
        new flameletBoussinesq<Specie>(*this)
    );
}


template<class Specie>
inline Foam::autoPtr<Foam::flameletBoussinesq<Specie>>
Foam::flameletBoussinesq<Specie>::New
(
    const dictionary& dict
)
{
    return autoPtr<flameletBoussinesq<Specie>>
    (
        new flameletBoussinesq<Specie>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::scalar Foam::flameletBoussinesq<Specie>::rho
(
    scalar p,
    scalar T
) const
{
    return rho0_*(1.0 - beta_*(T - T0_));
}


template<class Specie>
inline Foam::scalar Foam::flameletBoussinesq<Specie>::H(scalar p, scalar T) const
{
    return p/this->rho(p, T);
}


template<class Specie>
inline Foam::scalar Foam::flameletBoussinesq<Specie>::Cp(scalar p, scalar T) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::flameletBoussinesq<Specie>::E(scalar p, scalar T) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::flameletBoussinesq<Specie>::Cv(scalar p, scalar T) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::flameletBoussinesq<Specie>::S
(
    scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::flameletBoussinesq<Specie>::psi
(
    scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::flameletBoussinesq<Specie>::Z
(
    scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::flameletBoussinesq<Specie>::CpMCv
(
    scalar p,
    scalar T
) const
{
    return 0;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie>
inline void Foam::flameletBoussinesq<Specie>::operator=
(
    const flameletBoussinesq<Specie>& b
)
{
    Specie::operator=(b);

    rho0_ = b.rho0_;
    T0_ = b.T0_;
    beta_ = b.beta_;
}


template<class Specie>
inline void Foam::flameletBoussinesq<Specie>::operator+=
(
    const flameletBoussinesq<Specie>& b
)
{
    scalar Y1 = this->Y();
    Specie::operator+=(b);

    if (mag(this->Y()) > small)
    {
        Y1 /= this->Y();
        const scalar Y2 = b.Y()/this->Y();

        rho0_ = Y1*rho0_ + Y2*b.rho0_;
        T0_ = Y1*T0_ + Y2*b.T0_;
        beta_ = Y1*beta_ + Y2*b.beta_;
    }
}


template<class Specie>
inline void Foam::flameletBoussinesq<Specie>::operator*=(const scalar s)
{
    Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::flameletBoussinesq<Specie> Foam::operator+
(
    const flameletBoussinesq<Specie>& b1,
    const flameletBoussinesq<Specie>& b2
)
{
    Specie sp(static_cast<const Specie&>(b1) + static_cast<const Specie&>(b2));

    if (mag(sp.Y()) < small)
    {
        return flameletBoussinesq<Specie>
        (
            sp,
            b1.rho0_,
            b1.T0_,
            b1.beta_
        );
    }
    else
    {
        const scalar Y1 = b1.Y()/sp.Y();
        const scalar Y2 = b2.Y()/sp.Y();

        return flameletBoussinesq<Specie>
        (
            sp,
            Y1*b1.rho0_ + Y2*b2.rho0_,
            Y1*b1.T0_ + Y2*b2.T0_,
            Y1*b1.beta_ + Y2*b2.beta_
        );
    }
}


template<class Specie>
inline Foam::flameletBoussinesq<Specie> Foam::operator*
(
    const scalar s,
    const flameletBoussinesq<Specie>& b
)
{
    return flameletBoussinesq<Specie>
    (
        s*static_cast<const Specie&>(b),
        b.rho0_,
        b.T0_,
        b.beta_
    );
}


template<class Specie>
inline Foam::flameletBoussinesq<Specie> Foam::operator==
(
    const flameletBoussinesq<Specie>& b1,
    const flameletBoussinesq<Specie>& b2
)
{
    Specie sp(static_cast<const Specie&>(b1) == static_cast<const Specie&>(b2));

    const scalar Y1 = b1.Y()/sp.Y();
    const scalar Y2 = b2.Y()/sp.Y();

    return flameletBoussinesq<Specie>
    (
        sp,
        Y2*b2.rho0_ - Y1*b1.rho0_,
        Y2*b2.T0_   - Y1*b1.T0_,
        Y2*b2.beta_ - Y1*b1.beta_
    );
}


// ************************************************************************* //
