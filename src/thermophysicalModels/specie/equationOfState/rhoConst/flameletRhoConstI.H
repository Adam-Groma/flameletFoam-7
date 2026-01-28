/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2018 OpenFOAM Foundation
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

#include "flameletRhoConst.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Specie>
inline Foam::flameletRhoConst<Specie>::flameletRhoConst
(
    const Specie& sp,
    const scalar rho
)
:
    Specie(sp),
    rho_(rho)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::flameletRhoConst<Specie>::flameletRhoConst
(
    const word& name,
    const flameletRhoConst<Specie>& ico
)
:
    Specie(name, ico),
    rho_(ico.rho_)
{}


template<class Specie>
inline Foam::autoPtr<Foam::flameletRhoConst<Specie>>
Foam::flameletRhoConst<Specie>::clone() const
{
    return autoPtr<flameletRhoConst<Specie>>(new flameletRhoConst<Specie>(*this));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::scalar Foam::flameletRhoConst<Specie>::rho(scalar p, scalar T) const
{
    return rho_;
}


template<class Specie>
inline Foam::scalar Foam::flameletRhoConst<Specie>::H(scalar p, scalar T) const
{
    return p/rho_;
}


template<class Specie>
inline Foam::scalar Foam::flameletRhoConst<Specie>::Cp(scalar p, scalar T) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::flameletRhoConst<Specie>::E(scalar p, scalar T) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::flameletRhoConst<Specie>::Cv(scalar p, scalar T) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::flameletRhoConst<Specie>::S(scalar p, scalar T) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::flameletRhoConst<Specie>::psi(scalar p, scalar T) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::flameletRhoConst<Specie>::Z(scalar p, scalar T) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::flameletRhoConst<Specie>::CpMCv(scalar p, scalar T) const
{
    return 0;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie>
inline void Foam::flameletRhoConst<Specie>::operator+=(const flameletRhoConst<Specie>& ico)
{
    scalar Y1 = this->Y();
    Specie::operator+=(ico);

    if (mag(this->Y()) > small)
    {
        Y1 /= this->Y();
        const scalar Y2 = ico.Y()/this->Y();

        rho_ = Y1*rho_ + Y2*ico.rho_;
    }
}


template<class Specie>
inline void Foam::flameletRhoConst<Specie>::operator*=(const scalar s)
{
    Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::flameletRhoConst<Specie> Foam::operator+
(
    const flameletRhoConst<Specie>& ico1,
    const flameletRhoConst<Specie>& ico2
)
{
    Specie sp
    (
        static_cast<const Specie&>(ico1)
      + static_cast<const Specie&>(ico2)
    );

    if (mag(sp.Y()) < small)
    {
        return flameletRhoConst<Specie>
        (
            sp,
            ico1.rho_
        );
    }
    else
    {
        const scalar Y1 = ico1.Y()/sp.Y();
        const scalar Y2 = ico2.Y()/sp.Y();

        return flameletRhoConst<Specie>
        (
            sp,
            Y1*ico1.rho_ + Y2*ico2.rho_
        );
    }
}


template<class Specie>
inline Foam::flameletRhoConst<Specie> Foam::operator*
(
    const scalar s,
    const flameletRhoConst<Specie>& ico
)
{
    return flameletRhoConst<Specie>(s*static_cast<const Specie&>(ico), ico.rho_);
}


template<class Specie>
inline Foam::flameletRhoConst<Specie> Foam::operator==
(
    const flameletRhoConst<Specie>& ico1,
    const flameletRhoConst<Specie>& ico2
)
{
    Specie sp
    (
        static_cast<const Specie&>(ico1)
     == static_cast<const Specie&>(ico2)
    );

    const scalar Y1 = ico1.Y()/sp.Y();
    const scalar Y2 = ico2.Y()/sp.Y();

    return flameletRhoConst<Specie>
    (
        sp,
        Y2*ico2.rho_ - Y1*ico1.rho_
    );
}


// ************************************************************************* //
