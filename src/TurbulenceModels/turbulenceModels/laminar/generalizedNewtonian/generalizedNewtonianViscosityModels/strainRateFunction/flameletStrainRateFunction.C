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

#include "flameletStrainRateFunction.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarModels
{
namespace generalizedNewtonianViscosityModels
{
    defineTypeNameAndDebug(flameletStrainRateFunction, 0);

    addToRunTimeSelectionTable
    (
        flameletGeneralizedNewtonianViscosityModel,
        flameletStrainRateFunction,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::laminarModels::generalizedNewtonianViscosityModels::flameletStrainRateFunction::
flameletStrainRateFunction
(
    const dictionary& viscosityProperties
)
:
    flameletGeneralizedNewtonianViscosityModel(viscosityProperties),
    strainRateFunction_
    (
        Function1<scalar>::New
        (
            "function",
            viscosityProperties.optionalSubDict(typeName + "Coeffs")
        )
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::laminarModels::generalizedNewtonianViscosityModels::
flameletStrainRateFunction::read
(
    const dictionary& viscosityProperties
)
{
    flameletGeneralizedNewtonianViscosityModel::read(viscosityProperties);

    strainRateFunction_.clear();
    strainRateFunction_ = Function1<scalar>::New
    (
        "function",
        viscosityProperties.optionalSubDict
        (
            typeName + "Coeffs"
        )
    );

    return true;
}


Foam::tmp<Foam::volScalarField>
Foam::laminarModels::generalizedNewtonianViscosityModels::flameletStrainRateFunction::
nu
(
    const volScalarField& nu0,
    const volScalarField& strainRate
) const
{
    tmp<volScalarField> tnu
    (
        volScalarField::New
        (
            IOobject::groupName(type() + ":nu", nu0.group()),
            nu0.mesh(),
            dimensionedScalar(dimViscosity, 0)
        )
    );

    tnu.ref().primitiveFieldRef() = strainRateFunction_->value(strainRate);

    volScalarField::Boundary& nuBf = tnu.ref().boundaryFieldRef();
    const volScalarField::Boundary& sigmaBf = strainRate.boundaryField();

    forAll(nuBf, patchi)
    {
        nuBf[patchi] = strainRateFunction_->value(sigmaBf[patchi]);
    }

    return tnu;
}


// ************************************************************************* //
