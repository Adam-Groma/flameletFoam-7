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

inline Foam::flameletPowerSeriesReactionRate::flameletPowerSeriesReactionRate
(
    const scalar A,
    const scalar beta,
    const scalar Ta,
    const FixedList<scalar, nCoeff_> coeffs
)
:
    A_(A),
    beta_(beta),
    Ta_(Ta),
    coeffs_(coeffs)
{}


inline Foam::flameletPowerSeriesReactionRate::flameletPowerSeriesReactionRate
(
    const flameletSpeciesTable&,
    const dictionary& dict
)
:
    A_(readScalar(dict.lookup("A"))),
    beta_(readScalar(dict.lookup("beta"))),
    Ta_(readScalar(dict.lookup("Ta"))),
    coeffs_(dict.lookup("coeffs"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::flameletPowerSeriesReactionRate::operator()
(
    const scalar p,
    const scalar T,
    const scalarField&
) const
{
    scalar lta = A_;

    if (mag(beta_) > vSmall)
    {
        lta *= pow(T, beta_);
    }

    scalar expArg = 0;

    forAll(coeffs_, n)
    {
        expArg += coeffs_[n]/pow(T, n + 1);
    }

    lta *= exp(expArg);

    return lta;
}


inline Foam::scalar Foam::flameletPowerSeriesReactionRate::ddT
(
    const scalar p,
    const scalar T,
    const scalarField&
) const
{
    scalar lta = A_;

    if (mag(beta_) > vSmall)
    {
        lta *= pow(T, beta_);
    }

    scalar expArg = 0;
    scalar deriv = 0;

    forAll(coeffs_, n)
    {
        scalar cT = coeffs_[n]/pow(T, n + 1);
        expArg += cT;
        deriv -= (n + 1)*cT;
    }

    lta *= exp(expArg);

    return lta*(beta_+deriv)/T;
}


inline const Foam::List<Foam::Tuple2<Foam::label, Foam::scalar>>&
Foam::flameletPowerSeriesReactionRate::beta() const
{
    return NullObjectRef<List<Tuple2<label, scalar>>>();
}


inline void Foam::flameletPowerSeriesReactionRate::dcidc
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    scalarField& dcidc
) const
{}


inline Foam::scalar Foam::flameletPowerSeriesReactionRate::dcidT
(
    const scalar p,
    const scalar T,
    const scalarField& c
) const
{
    return 0;
}


inline void Foam::flameletPowerSeriesReactionRate::write(Ostream& os) const
{
    writeEntry(os, "A", A_);
    writeEntry(os, "beta", beta_);
    writeEntry(os, "Ta", Ta_);
    writeEntry(os, "coeffs", coeffs_);
}


inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const flameletPowerSeriesReactionRate& psrr
)
{
    psrr.write(os);
    return os;
}


// ************************************************************************* //
