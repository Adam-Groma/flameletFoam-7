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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::flameletConstTransport<Thermo>::flameletConstTransport
(
    const Thermo& t,
    const scalar mu,
    const scalar Pr
)
:
    Thermo(t),
    mu_(mu),
    rPr_(1.0/Pr)
{}


template<class Thermo>
inline Foam::flameletConstTransport<Thermo>::flameletConstTransport
(
    const word& name,
    const flameletConstTransport& ct
)
:
    Thermo(name, ct),
    mu_(ct.mu_),
    rPr_(ct.rPr_)
{}


template<class Thermo>
inline Foam::autoPtr<Foam::flameletConstTransport<Thermo>>
Foam::flameletConstTransport<Thermo>::clone() const
{
    return autoPtr<flameletConstTransport<Thermo>>
    (
        new flameletConstTransport<Thermo>(*this)
    );
}


template<class Thermo>
inline Foam::autoPtr<Foam::flameletConstTransport<Thermo>>
Foam::flameletConstTransport<Thermo>::New
(
    const dictionary& dict
)
{
    return autoPtr<flameletConstTransport<Thermo>>
    (
        new flameletConstTransport<Thermo>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::flameletConstTransport<Thermo>::mu
(
    const scalar p,
    const scalar T
) const
{
    return mu_;
}


template<class Thermo>
inline Foam::scalar Foam::flameletConstTransport<Thermo>::kappa
(
    const scalar p,
    const scalar T
) const
{
    return this->Cpv(p, T)*mu(p, T)*rPr_;
}


template<class Thermo>
inline Foam::scalar Foam::flameletConstTransport<Thermo>::alphah
(
    const scalar p,
    const scalar T
) const
{
    return mu(p, T)*rPr_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline void Foam::flameletConstTransport<Thermo>::operator+=
(
    const flameletConstTransport<Thermo>& st
)
{
    scalar Y1 = this->Y();

    Thermo::operator+=(st);

    Y1 /= this->Y();
    scalar Y2 = st.Y()/this->Y();

    mu_ = Y1*mu_ + Y2*st.mu_;
    rPr_ = 1.0/(Y1/rPr_ + Y2/st.rPr_);
}


template<class Thermo>
inline void Foam::flameletConstTransport<Thermo>::operator*=
(
    const scalar s
)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::flameletConstTransport<Thermo> Foam::operator+
(
    const flameletConstTransport<Thermo>& ct1,
    const flameletConstTransport<Thermo>& ct2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(ct1) + static_cast<const Thermo&>(ct2)
    );

    scalar Y1 = ct1.Y()/t.Y();
    scalar Y2 = ct2.Y()/t.Y();

    return flameletConstTransport<Thermo>
    (
        t,
        Y1*ct1.mu_ + Y2*ct2.mu_,
        1.0/(Y1/ct1.rPr_ + Y2/ct2.rPr_)
    );
}


template<class Thermo>
inline Foam::flameletConstTransport<Thermo> Foam::operator*
(
    const scalar s,
    const flameletConstTransport<Thermo>& ct
)
{
    return flameletConstTransport<Thermo>
    (
        s*static_cast<const Thermo&>(ct),
        ct.mu_,
        1.0/ct.rPr_
    );
}


// ************************************************************************* //
