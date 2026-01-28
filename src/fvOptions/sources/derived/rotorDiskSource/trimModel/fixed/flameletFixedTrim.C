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

#include "flameletFixedTrim.H"
#include "addToRunTimeSelectionTable.H"
#include "unitConversion.H"
#include "mathematicalConstants.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(flameletFixedTrim, 0);

    addToRunTimeSelectionTable(flameletTrimModel, flameletFixedTrim, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::flameletFixedTrim::flameletFixedTrim
(
    const fv::flameletRotorDiskSource& rotor,
    const dictionary& dict
)
:
    flameletTrimModel(rotor, dict, typeName),
    thetag_(rotor.cells().size(), 0.0)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

Foam::flameletFixedTrim::~flameletFixedTrim()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::flameletFixedTrim::read(const dictionary& dict)
{
    flameletTrimModel::read(dict);

    scalar theta0 = degToRad(readScalar(coeffs_.lookup("theta0")));
    scalar theta1c = degToRad(readScalar(coeffs_.lookup("theta1c")));
    scalar theta1s = degToRad(readScalar(coeffs_.lookup("theta1s")));

    const List<point>& x = rotor_.x();
    forAll(thetag_, i)
    {
        scalar psi = x[i].y();
        thetag_[i] = theta0 + theta1c*cos(psi) + theta1s*sin(psi);
    }
}


Foam::tmp<Foam::scalarField> Foam::flameletFixedTrim::thetag() const
{
    return tmp<scalarField>(thetag_);
}


void Foam::flameletFixedTrim::correct
(
    const vectorField& U,
    vectorField& force
)
{}


void Foam::flameletFixedTrim::correct
(
    const volScalarField rho,
    const vectorField& U,
    vectorField& force)
{}


// ************************************************************************* //
