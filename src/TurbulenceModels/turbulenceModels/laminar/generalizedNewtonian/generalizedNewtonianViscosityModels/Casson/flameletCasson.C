/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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

#include "flameletCasson.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarModels
{
namespace generalizedNewtonianViscosityModels
{
    defineTypeNameAndDebug(flameletCasson, 0);
    addToRunTimeSelectionTable
    (
        flameletGeneralizedNewtonianViscosityModel,
        flameletCasson,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::laminarModels::generalizedNewtonianViscosityModels::flameletCasson::flameletCasson
(
    const dictionary& viscosityProperties
)
:
    flameletGeneralizedNewtonianViscosityModel(viscosityProperties),
    m_("m", dimViscosity, 0),
    tau0_("tau0", dimViscosity/dimTime, 0),
    nuMin_("nuMin", dimViscosity, 0),
    nuMax_("nuMax", dimViscosity, 0)
{
    read(viscosityProperties);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::laminarModels::generalizedNewtonianViscosityModels::flameletCasson::read
(
    const dictionary& viscosityProperties
)
{
    flameletGeneralizedNewtonianViscosityModel::read(viscosityProperties);

    const dictionary& coeffs =
        viscosityProperties.optionalSubDict(typeName + "Coeffs");

    coeffs.lookup("m") >> m_;
    coeffs.lookup("tau0") >> tau0_;
    coeffs.lookup("nuMin_") >> nuMin_;
    coeffs.lookup("nuMax_") >> nuMax_;

    return true;
}


Foam::tmp<Foam::volScalarField>
Foam::laminarModels::generalizedNewtonianViscosityModels::flameletCasson::
nu
(
    const volScalarField& nu0,
    const volScalarField& strainRate
) const
{
    return max
    (
        nuMin_,
        min
        (
            nuMax_,
            sqr
            (
                sqrt
                (
                    tau0_
                   /max
                    (
                        strainRate,
                        dimensionedScalar(dimless/dimTime, vSmall)
                    )
                ) + sqrt(m_)
            )
        )
    );
}


// ************************************************************************* //
