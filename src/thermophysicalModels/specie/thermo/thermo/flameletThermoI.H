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

#include "flameletThermo.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Thermo, template<class> class Type>
inline Foam::species::flameletThermo<Thermo, Type>::flameletThermo
(
    const Thermo& sp
)
:
    Thermo(sp)
{}


template<class Thermo, template<class> class Type>
inline Foam::scalar Foam::species::flameletThermo<Thermo, Type>::T
(
    const scalar f,
    const scalar p,
    const scalar T0,
    scalar (flameletThermo<Thermo, Type>::*F)(const scalar, const scalar) const,
    scalar (flameletThermo<Thermo, Type>::*dFdT)(const scalar, const scalar)
        const,
    scalar (flameletThermo<Thermo, Type>::*limit)(const scalar) const,
    const bool diagnostics
) const
{
    if (T0 < 0)
    {
        FatalErrorInFunction
            << "Negative initial temperature T0: " << T0
            << abort(FatalError);
    }

    scalar Test = T0;
    scalar Tnew = T0;
    scalar Ttol = T0*tol_;
    int    iter = 0;

    if (diagnostics)
    {
        const unsigned int width = IOstream::defaultPrecision() + 8;

        InfoInFunction
            << "Energy -> temperature conversion failed to converge:" << endl;
        Pout<< setw(width) << "iter"
            << setw(width) << "Test"
            << setw(width) << "e/h"
            << setw(width) << "Cv/p"
            << setw(width) << "Tnew"
            << endl;
    }

    do
    {
        Test = Tnew;
        Tnew =
            (this->*limit)
            (Test - ((this->*F)(p, Test) - f)/(this->*dFdT)(p, Test));

        if (diagnostics)
        {
            const unsigned int width = IOstream::defaultPrecision() + 8;

            Pout<< setw(width) << iter
                << setw(width) << Test
                << setw(width) << ((this->*F)(p, Test))
                << setw(width) << ((this->*dFdT)(p, Test))
                << setw(width) << Tnew
                << endl;
        }

        if (iter++ > maxIter_)
        {
            if (!diagnostics)
            {
                T(f, p, T0, F, dFdT, limit, true);
            }

            FatalErrorInFunction
                << "Maximum number of iterations exceeded: " << maxIter_
                << abort(FatalError);
        }

    } while (mag(Tnew - Test) > Ttol);

    return Tnew;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo, template<class> class Type>
inline Foam::species::flameletThermo<Thermo, Type>::flameletThermo
(
    const word& name,
    const flameletThermo& st
)
:
    Thermo(name, st)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo, template<class> class Type>
inline Foam::word
Foam::species::flameletThermo<Thermo, Type>::heName()
{
    return Type<flameletThermo<Thermo, Type>>::energyName();
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::flameletThermo<Thermo, Type>::Cpv(const scalar p, const scalar T) const
{
    return Type<flameletThermo<Thermo, Type>>::Cpv(*this, p, T);
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::flameletThermo<Thermo, Type>::gamma(const scalar p, const scalar T) const
{
    const scalar Cp = this->Cp(p, T);
    return Cp/(Cp - this->CpMCv(p, T));
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::flameletThermo<Thermo, Type>::CpByCpv
(
    const scalar p,
    const scalar T
) const
{
    return Type<flameletThermo<Thermo, Type>>::CpByCpv(*this, p, T);
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::flameletThermo<Thermo, Type>::HE(const scalar p, const scalar T) const
{
    return Type<flameletThermo<Thermo, Type>>::HE(*this, p, T);
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::flameletThermo<Thermo, Type>::G(const scalar p, const scalar T) const
{
    return this->Ha(p, T) - T*this->S(p, T);
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::flameletThermo<Thermo, Type>::A(const scalar p, const scalar T) const
{
    return this->Ea(p, T) - T*this->S(p, T);
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::flameletThermo<Thermo, Type>::cp(const scalar p, const scalar T) const
{
    return this->Cp(p, T)*this->W();
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::flameletThermo<Thermo, Type>::ha(const scalar p, const scalar T) const
{
    return this->Ha(p, T)*this->W();
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::flameletThermo<Thermo, Type>::hs(const scalar p, const scalar T) const
{
    return this->Hs(p, T)*this->W();
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::flameletThermo<Thermo, Type>::hc() const
{
    return this->Hc()*this->W();
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::flameletThermo<Thermo, Type>::s(const scalar p, const scalar T) const
{
    return this->S(p, T)*this->W();
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::flameletThermo<Thermo, Type>::he(const scalar p, const scalar T) const
{
    return this->HE(p, T)*this->W();
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::flameletThermo<Thermo, Type>::cv(const scalar p, const scalar T) const
{
    return this->Cv(p, T)*this->W();
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::flameletThermo<Thermo, Type>::es(const scalar p, const scalar T) const
{
    return this->Es(p, T)*this->W();
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::flameletThermo<Thermo, Type>::ea(const scalar p, const scalar T) const
{
    return this->Ea(p, T)*this->W();
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::flameletThermo<Thermo, Type>::g(const scalar p, const scalar T) const
{
    return this->G(p, T)*this->W();
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::flameletThermo<Thermo, Type>::a(const scalar p, const scalar T) const
{
    return this->A(p, T)*this->W();
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::flameletThermo<Thermo, Type>::K(const scalar p, const scalar T) const
{
    scalar arg = -this->Y()*this->G(Pstd, T)/(RR*T);

    if (arg < 600)
    {
        return exp(arg);
    }
    else
    {
        return rootVGreat;
    }
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::flameletThermo<Thermo, Type>::Kp(const scalar p, const scalar T) const
{
    return K(p, T);
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::flameletThermo<Thermo, Type>::Kc(const scalar p, const scalar T) const
{
    const scalar nm = this->Y()/this->W();

    if (equal(nm, small))
    {
        return Kp(p, T);
    }
    else
    {
        return Kp(p, T)*pow(Pstd/(RR*T), nm);
    }
}


template<class Thermo, template<class> class Type>
inline Foam::scalar Foam::species::flameletThermo<Thermo, Type>::Kx
(
    const scalar p,
    const scalar T
) const
{
    const scalar nm = this->Y()/this->W();

    if (equal(nm, small))
    {
        return Kp(p, T);
    }
    else
    {
        return Kp(p, T)*pow(Pstd/p, nm);
    }
}


template<class Thermo, template<class> class Type>
inline Foam::scalar Foam::species::flameletThermo<Thermo, Type>::Kn
(
    const scalar p,
    const scalar T,
    const scalar n
) const
{
    const scalar nm = this->Y()/this->W();

    if (equal(nm, small))
    {
        return Kp(p, T);
    }
    else
    {
        return Kp(p, T)*pow(n*Pstd/p, nm);
    }
}


template<class Thermo, template<class> class Type>
inline Foam::scalar Foam::species::flameletThermo<Thermo, Type>::THE
(
    const scalar he,
    const scalar p,
    const scalar T0
) const
{
    return Type<flameletThermo<Thermo, Type>>::THE(*this, he, p, T0);
}


template<class Thermo, template<class> class Type>
inline Foam::scalar Foam::species::flameletThermo<Thermo, Type>::THs
(
    const scalar hs,
    const scalar p,
    const scalar T0
) const
{
    return T
    (
        hs,
        p,
        T0,
        &flameletThermo<Thermo, Type>::Hs,
        &flameletThermo<Thermo, Type>::Cp,
        &flameletThermo<Thermo, Type>::limit
    );
}


template<class Thermo, template<class> class Type>
inline Foam::scalar Foam::species::flameletThermo<Thermo, Type>::THa
(
    const scalar ha,
    const scalar p,
    const scalar T0
) const
{
    return T
    (
        ha,
        p,
        T0,
        &flameletThermo<Thermo, Type>::Ha,
        &flameletThermo<Thermo, Type>::Cp,
        &flameletThermo<Thermo, Type>::limit
    );
}


template<class Thermo, template<class> class Type>
inline Foam::scalar Foam::species::flameletThermo<Thermo, Type>::TEs
(
    const scalar es,
    const scalar p,
    const scalar T0
) const
{
    return T
    (
        es,
        p,
        T0,
        &flameletThermo<Thermo, Type>::Es,
        &flameletThermo<Thermo, Type>::Cv,
        &flameletThermo<Thermo, Type>::limit
    );
}


template<class Thermo, template<class> class Type>
inline Foam::scalar Foam::species::flameletThermo<Thermo, Type>::TEa
(
    const scalar ea,
    const scalar p,
    const scalar T0
) const
{
    return T
    (
        ea,
        p,
        T0,
        &flameletThermo<Thermo, Type>::Ea,
        &flameletThermo<Thermo, Type>::Cv,
        &flameletThermo<Thermo, Type>::limit
    );
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::flameletThermo<Thermo, Type>::dKcdTbyKc
(
    const scalar p,
    const scalar T
) const
{
    const scalar nm = this->Y()/this->W();

    if (equal(nm, small))
    {
        return -this->dGdT(Pstd, T)*this->Y()/RR;
    }
    else
    {
        return -(nm/T + this->dGdT(Pstd, T)*this->Y()/RR);
    }
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::flameletThermo<Thermo, Type>::dcpdT(const scalar p, const scalar T) const
{
    return this->dCpdT(p, T)*this->W();;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo, template<class> class Type>
inline void Foam::species::flameletThermo<Thermo, Type>::operator+=
(
    const flameletThermo<Thermo, Type>& st
)
{
    Thermo::operator+=(st);
}


template<class Thermo, template<class> class Type>
inline void Foam::species::flameletThermo<Thermo, Type>::operator*=(const scalar s)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo, template<class> class Type>
inline Foam::species::flameletThermo<Thermo, Type> Foam::species::operator+
(
    const flameletThermo<Thermo, Type>& st1,
    const flameletThermo<Thermo, Type>& st2
)
{
    return flameletThermo<Thermo, Type>
    (
        static_cast<const Thermo&>(st1) + static_cast<const Thermo&>(st2)
    );
}


template<class Thermo, template<class> class Type>
inline Foam::species::flameletThermo<Thermo, Type> Foam::species::operator*
(
    const scalar s,
    const flameletThermo<Thermo, Type>& st
)
{
    return flameletThermo<Thermo, Type>
    (
        s*static_cast<const Thermo&>(st)
    );
}


template<class Thermo, template<class> class Type>
inline Foam::species::flameletThermo<Thermo, Type> Foam::species::operator==
(
    const flameletThermo<Thermo, Type>& st1,
    const flameletThermo<Thermo, Type>& st2
)
{
    return flameletThermo<Thermo, Type>
    (
        static_cast<const Thermo&>(st1) == static_cast<const Thermo&>(st2)
    );
}


// ************************************************************************* //
