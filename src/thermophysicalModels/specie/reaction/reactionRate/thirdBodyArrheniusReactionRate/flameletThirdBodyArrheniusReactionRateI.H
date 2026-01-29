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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::flameletThirdBodyArrheniusReactionRate::flameletThirdBodyArrheniusReactionRate
(
    const scalar A,
    const scalar beta,
    const scalar Ta,
    const flameletThirdBodyEfficiencies& tbes
)
:
    flameletArrheniusReactionRate(A, beta, Ta),
    thirdBodyEfficiencies_(tbes)
{
    forAll(tbes, i)
    {
        beta_.append(Tuple2<label, scalar>(i, tbes[i]));
    }
}


inline Foam::flameletThirdBodyArrheniusReactionRate::flameletThirdBodyArrheniusReactionRate
(
    const flameletSpeciesTable& species,
    const dictionary& dict
)
:
    flameletArrheniusReactionRate
    (
        species,
        dict
    ),
    thirdBodyEfficiencies_(species, dict)
{
    forAll(thirdBodyEfficiencies_, i)
    {
        beta_.append
        (
            Tuple2<label, scalar>
            (
                i,
                thirdBodyEfficiencies_[i]
            )
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::flameletThirdBodyArrheniusReactionRate::operator()
(
    const scalar p,
    const scalar T,
    const scalarField& c
) const
{
    return
        thirdBodyEfficiencies_.M(c)
       *flameletArrheniusReactionRate::operator()(p, T, c);
}


inline Foam::scalar Foam::flameletThirdBodyArrheniusReactionRate::ddT
(
    const scalar p,
    const scalar T,
    const scalarField& c
) const
{
    return
        thirdBodyEfficiencies_.M(c)
       *flameletArrheniusReactionRate::ddT(p, T, c);
}


inline void Foam::flameletThirdBodyArrheniusReactionRate::dcidc
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    scalarField& dcidc
) const
{
    scalar M = thirdBodyEfficiencies_.M(c);
    forAll(beta_, i)
    {
        dcidc[i] = beta_[i].second()/max(M, small);
    }
}


inline Foam::scalar Foam::flameletThirdBodyArrheniusReactionRate::dcidT
(
    const scalar p,
    const scalar T,
    const scalarField& c
) const
{
    return -1.0/T;
}


inline void Foam::flameletThirdBodyArrheniusReactionRate::write(Ostream& os) const
{
    flameletArrheniusReactionRate::write(os);
    thirdBodyEfficiencies_.write(os);
}


inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const flameletThirdBodyArrheniusReactionRate& arr
)
{
    arr.write(os);
    return os;
}


// ************************************************************************* //
