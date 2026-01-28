/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019 OpenFOAM Foundation
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

#include "flameletIsotropicDamping.H"
#include "fvMatrix.H"
#include "fvmSup.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(flameletIsotropicDamping, 0);
    addToRunTimeSelectionTable(option, flameletIsotropicDamping, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::flameletIsotropicDamping::add
(
    const volScalarField::Internal& forceCoeff,
    fvMatrix<vector>& eqn
)
{
    eqn -= fvm::Sp(forceCoeff, eqn.psi());
    eqn += forceCoeff*value_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::flameletIsotropicDamping::flameletIsotropicDamping
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    flameletDamping(name, modelType, dict, mesh),
    value_("value", dimVelocity, coeffs_.lookup("value"))
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::flameletIsotropicDamping::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    add(this->forceCoeff(), eqn);
}


void Foam::fv::flameletIsotropicDamping::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    add(rho*forceCoeff(), eqn);
}


void Foam::fv::flameletIsotropicDamping::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    add(alpha()*rho()*this->forceCoeff(), eqn);
}


bool Foam::fv::flameletIsotropicDamping::read(const dictionary& dict)
{
    if (flameletDamping::read(dict))
    {
        value_ =
            dimensionedVector
            (
                value_.name(),
                value_.dimensions(),
                coeffs_.lookup(value_.name())
            );

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
