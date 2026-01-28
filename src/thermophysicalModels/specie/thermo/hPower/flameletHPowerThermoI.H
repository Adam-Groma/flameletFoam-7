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

#include "flameletHPowerThermo.H"
#include "flameletSpecie.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class EquationOfState>
inline void Foam::flameletHPowerThermo<EquationOfState>::checkT
(
    const scalar T
) const
{
    if (T < 0)
    {
        FatalErrorInFunction
            << "attempt to evaluate flameletHPowerThermo<EquationOfState>"
               " for negative temperature " << T
            << abort(FatalError);
    }
}


template<class EquationOfState>
inline Foam::flameletHPowerThermo<EquationOfState>::flameletHPowerThermo
(
    const word& name,
    const flameletHPowerThermo& jt
)
:
    EquationOfState(name, jt),
    c0_(jt.c0_),
    n0_(jt.n0_),
    Tref_(jt.Tref_),
    Hf_(jt.Hf_)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::flameletHPowerThermo<EquationOfState>::flameletHPowerThermo
(
    const EquationOfState& st,
    const scalar c0,
    const scalar n0,
    const scalar Tref,
    const scalar Hf
)
:
    EquationOfState(st),
    c0_(c0),
    n0_(n0),
    Tref_(Tref),
    Hf_(Hf)
{}


template<class EquationOfState>
inline Foam::autoPtr<Foam::flameletHPowerThermo<EquationOfState>>
Foam::flameletHPowerThermo<EquationOfState>::clone() const
{
    return autoPtr<flameletHPowerThermo<EquationOfState>>
    (
        new flameletHPowerThermo<EquationOfState>(*this)
    );
}


template<class EquationOfState>
inline Foam::autoPtr<Foam::flameletHPowerThermo<EquationOfState>>
Foam::flameletHPowerThermo<EquationOfState>::New(const dictionary& dict)
{
    return autoPtr<flameletHPowerThermo<EquationOfState>>
    (
        new flameletHPowerThermo<EquationOfState>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::scalar Foam::flameletHPowerThermo<EquationOfState>::limit
(
    const scalar T
) const
{
    return T;
}


template<class EquationOfState>
inline Foam::scalar Foam::flameletHPowerThermo<EquationOfState>::Cp
(
    const scalar p, const scalar T
) const
{
    return c0_*pow(T/Tref_, n0_) + EquationOfState::Cp(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::flameletHPowerThermo<EquationOfState>::Ha
(
    const scalar p, const scalar T
) const
{
    return Hs(p, T) + Hc();
}


template<class EquationOfState>
inline Foam::scalar Foam::flameletHPowerThermo<EquationOfState>::Hs
(
    const scalar p, const scalar T
) const
{
    return
        c0_*(pow(T, n0_ + 1) - pow(Tstd, n0_ + 1))/(pow(Tref_, n0_)*(n0_ + 1))
      + EquationOfState::H(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::flameletHPowerThermo<EquationOfState>::Hc() const
{
    return Hf_;
}


template<class EquationOfState>
inline Foam::scalar Foam::flameletHPowerThermo<EquationOfState>::S
(
    const scalar p, const scalar T
) const
{
    return
        c0_*(pow(T, n0_) - pow(Tstd, n0_))/(pow(Tref_, n0_)*n0_)
      + EquationOfState::S(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::flameletHPowerThermo<EquationOfState>::dGdT
(
    const scalar p, const scalar T
) const
{
    // To be implemented
    NotImplemented;
    return 0;
}


template<class EquationOfState>
inline Foam::scalar Foam::flameletHPowerThermo<EquationOfState>::dCpdT
(
    const scalar p, const scalar T
) const
{
    // To be implemented
    NotImplemented;
    return 0;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline void Foam::flameletHPowerThermo<EquationOfState>::operator+=
(
    const flameletHPowerThermo<EquationOfState>& ct
)
{
    scalar Y1 = this->Y();

    EquationOfState::operator+=(ct);

    if (mag(this->Y()) > small)
    {
        Y1 /= this->Y();
        const scalar Y2 = ct.Y()/this->Y();

        Hf_ = Y1*Hf_ + Y2*ct.Hf_;
        c0_ = Y1*c0_ + Y2*ct.c0_;
        n0_ = Y1*n0_ + Y2*ct.n0_;
        Tref_ = Y1*Tref_ + Y2*ct.Tref_;
    }
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::flameletHPowerThermo<EquationOfState> Foam::operator+
(
    const flameletHPowerThermo<EquationOfState>& ct1,
    const flameletHPowerThermo<EquationOfState>& ct2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(ct1)
      + static_cast<const EquationOfState&>(ct2)
    );

    if (mag(eofs.Y()) < small)
    {
        return flameletHPowerThermo<EquationOfState>
        (
            eofs,
            ct1.c0_,
            ct1.n0_,
            ct1.Tref_,
            ct1.Hf_
        );
    }
    else
    {
        return flameletHPowerThermo<EquationOfState>
        (
            eofs,
            ct1.Y()/eofs.Y()*ct1.c0_
          + ct2.Y()/eofs.Y()*ct2.c0_,
            ct1.Y()/eofs.Y()*ct1.n0_
          + ct2.Y()/eofs.Y()*ct2.n0_,
            ct1.Y()/eofs.Y()*ct1.Tref_
          + ct2.Y()/eofs.Y()*ct2.Tref_,
            ct1.Y()/eofs.Y()*ct1.Hf_
          + ct2.Y()/eofs.Y()*ct2.Hf_
        );
    }
}


template<class EquationOfState>
inline Foam::flameletHPowerThermo<EquationOfState> Foam::operator*
(
    const scalar s,
    const flameletHPowerThermo<EquationOfState>& ct
)
{
    return flameletHPowerThermo<EquationOfState>
    (
        s*static_cast<const EquationOfState&>(ct),
        ct.c0_,
        ct.n0_,
        ct.Tref_,
        ct.Hf_
    );
}


template<class EquationOfState>
inline Foam::flameletHPowerThermo<EquationOfState> Foam::operator==
(
    const flameletHPowerThermo<EquationOfState>& ct1,
    const flameletHPowerThermo<EquationOfState>& ct2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(ct1)
     == static_cast<const EquationOfState&>(ct2)
    );

    return flameletHPowerThermo<EquationOfState>
    (
        eofs,
        ct2.Y()/eofs.Y()*ct2.c0_
      - ct1.Y()/eofs.Y()*ct1.c0_,
        ct2.Y()/eofs.Y()*ct2.n0_
      - ct1.Y()/eofs.Y()*ct1.n0_,
        ct2.Y()/eofs.Y()*ct2.Tref_
      - ct1.Y()/eofs.Y()*ct1.Tref_,
        ct2.Y()/eofs.Y()*ct2.Hf_
      - ct1.Y()/eofs.Y()*ct1.Hf_
    );
}


// ************************************************************************* //
