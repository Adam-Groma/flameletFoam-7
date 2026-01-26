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

#include "flameletPrandtlDelta.H"
#include "wallDist.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{
    defineTypeNameAndDebug(flameletPrandtlDelta, 0);
    addToRunTimeSelectionTable(flameletLESdelta, flameletPrandtlDelta, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::LESModels::flameletPrandtlDelta::calcDelta()
{
    delta_ = min
    (
        static_cast<const volScalarField&>(geometricDelta_()),
        (kappa_/Cdelta_)*wallDist::New(turbulenceModel_.mesh()).y()
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LESModels::flameletPrandtlDelta::flameletPrandtlDelta
(
    const word& name,
    const flameletTurbulenceModel& turbulence,
    const dictionary& dict
)
:
    flameletLESdelta(name, turbulence),
    geometricDelta_
    (
        flameletLESdelta::New
        (
            name,
            turbulence,
            dict.optionalSubDict(type() + "Coeffs")
        )
    ),
    kappa_(dict.lookupOrDefault<scalar>("kappa", 0.41)),
    Cdelta_
    (
        dict.optionalSubDict(type() + "Coeffs").lookupOrDefault<scalar>
        (
            "Cdelta",
            0.158
        )
    )
{
    calcDelta();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::LESModels::flameletPrandtlDelta::read(const dictionary& dict)
{
    const dictionary& coeffDict(dict.optionalSubDict(type() + "Coeffs"));

    geometricDelta_().read(coeffDict);
    dict.readIfPresent<scalar>("kappa", kappa_);
    coeffDict.readIfPresent<scalar>("Cdelta", Cdelta_);
    calcDelta();
}


void Foam::LESModels::flameletPrandtlDelta::correct()
{
    geometricDelta_().correct();

    if (turbulenceModel_.mesh().changing())
    {
        calcDelta();
    }
}


// ************************************************************************* //
