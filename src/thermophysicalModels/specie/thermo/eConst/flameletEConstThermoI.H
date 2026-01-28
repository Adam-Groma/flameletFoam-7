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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::flameletEConstThermo<EquationOfState>::flameletEConstThermo
(
    const EquationOfState& st,
    const scalar cv,
    const scalar hf
)
:
    EquationOfState(st),
    Cv_(cv),
    Hf_(hf)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::flameletEConstThermo<EquationOfState>::flameletEConstThermo
(
    const word& name,
    const flameletEConstThermo<EquationOfState>& ct
)
:
    EquationOfState(name, ct),
    Cv_(ct.Cv_),
    Hf_(ct.Hf_)
{}


template<class EquationOfState>
inline Foam::autoPtr<Foam::flameletEConstThermo<EquationOfState>>
Foam::flameletEConstThermo<EquationOfState>::clone() const
{
    return autoPtr<flameletEConstThermo<EquationOfState>>
    (
        new flameletEConstThermo<EquationOfState>(*this)
    );
}


template<class EquationOfState>
inline Foam::autoPtr<Foam::flameletEConstThermo<EquationOfState>>
Foam::flameletEConstThermo<EquationOfState>::New(const dictionary& dict)
{
    return autoPtr<flameletEConstThermo<EquationOfState>>
    (
        new flameletEConstThermo<EquationOfState>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::scalar Foam::flameletEConstThermo<EquationOfState>::limit
(
    const scalar T
) const
{
    return T;
}


template<class EquationOfState>
inline Foam::scalar Foam::flameletEConstThermo<EquationOfState>::Cv
(
    const scalar p,
    const scalar T
) const
{
    return Cv_ + EquationOfState::Cv(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::flameletEConstThermo<EquationOfState>::Es
(
    const scalar p,
    const scalar T
) const
{
    return Cv_*T + EquationOfState::E(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::flameletEConstThermo<EquationOfState>::Hc() const
{
    return Hf_;
}


template<class EquationOfState>
inline Foam::scalar Foam::flameletEConstThermo<EquationOfState>::Ea
(
    const scalar p,
    const scalar T
) const
{
    return Es(p, T) + Hc();
}


template<class EquationOfState>
inline Foam::scalar Foam::flameletEConstThermo<EquationOfState>::S
(
    const scalar p,
    const scalar T
) const
{
    return Cp(p, T)*log(T/Tstd) + EquationOfState::S(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::flameletEConstThermo<EquationOfState>::dGdT
(
    const scalar p,
    const scalar T
) const
{
    return 0;
}


template<class EquationOfState>
inline Foam::scalar Foam::flameletEConstThermo<EquationOfState>::dCpdT
(
    const scalar p,
    const scalar T
) const
{
    return 0;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline void Foam::flameletEConstThermo<EquationOfState>::operator+=
(
    const flameletEConstThermo<EquationOfState>& ct
)
{
    scalar Y1 = this->Y();

    EquationOfState::operator+=(ct);

    if (mag(this->Y()) > small)
    {
        Y1 /= this->Y();
        const scalar Y2 = ct.Y()/this->Y();

        Cv_ = Y1*Cv_ + Y2*ct.Cv_;
        Hf_ = Y1*Hf_ + Y2*ct.Hf_;
    }
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::flameletEConstThermo<EquationOfState> Foam::operator+
(
    const flameletEConstThermo<EquationOfState>& ct1,
    const flameletEConstThermo<EquationOfState>& ct2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(ct1)
      + static_cast<const EquationOfState&>(ct2)
    );

    if (mag(eofs.Y()) < small)
    {
        return flameletEConstThermo<EquationOfState>
        (
            eofs,
            ct1.Cv_,
            ct1.Hf_
        );
    }
    else
    {
        return flameletEConstThermo<EquationOfState>
        (
            eofs,
            ct1.Y()/eofs.Y()*ct1.Cv_
          + ct2.Y()/eofs.Y()*ct2.Cv_,
            ct1.Y()/eofs.Y()*ct1.Hf_
          + ct2.Y()/eofs.Y()*ct2.Hf_
        );
    }
}


template<class EquationOfState>
inline Foam::flameletEConstThermo<EquationOfState> Foam::operator*
(
    const scalar s,
    const flameletEConstThermo<EquationOfState>& ct
)
{
    return flameletEConstThermo<EquationOfState>
    (
        s*static_cast<const EquationOfState&>(ct),
        ct.Cv_,
        ct.Hf_
    );
}


template<class EquationOfState>
inline Foam::flameletEConstThermo<EquationOfState> Foam::operator==
(
    const flameletEConstThermo<EquationOfState>& ct1,
    const flameletEConstThermo<EquationOfState>& ct2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(ct1)
     == static_cast<const EquationOfState&>(ct2)
    );

    return flameletEConstThermo<EquationOfState>
    (
        eofs,
        ct2.Y()/eofs.Y()*ct2.Cv_
      - ct1.Y()/eofs.Y()*ct1.Cv_,
        ct2.Y()/eofs.Y()*ct2.Hf_
      - ct1.Y()/eofs.Y()*ct1.Hf_
    );
}


// ************************************************************************* //
