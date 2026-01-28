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

#include "flameletExplicitPorositySource.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "porosityModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(flameletExplicitPorositySource, 0);
    addToRunTimeSelectionTable
    (
        option,
        flameletExplicitPorositySource,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::flameletExplicitPorositySource::flameletExplicitPorositySource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    flameletCellSetOption(name, modelType, dict, mesh),
    porosityPtr_(nullptr)
{
    read(dict);

    porosityPtr_.reset
    (
        porosityModel::New
        (
            name_,
            mesh_,
            coeffs_,
            cellSetName_
        ).ptr()
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::flameletExplicitPorositySource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    fvMatrix<vector> porosityEqn(eqn.psi(), eqn.dimensions());
    porosityPtr_->addResistance(porosityEqn);
    eqn -= porosityEqn;
}


void Foam::fv::flameletExplicitPorositySource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    fvMatrix<vector> porosityEqn(eqn.psi(), eqn.dimensions());
    porosityPtr_->addResistance(porosityEqn);
    eqn -= porosityEqn;
}


void Foam::fv::flameletExplicitPorositySource::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    fvMatrix<vector> porosityEqn(eqn.psi(), eqn.dimensions());
    porosityPtr_->addResistance(porosityEqn);
    eqn -= alpha*porosityEqn;
}


bool Foam::fv::flameletExplicitPorositySource::read(const dictionary& dict)
{
    if (flameletCellSetOption::read(dict))
    {
        if (coeffs_.found("UNames"))
        {
            coeffs_.lookup("UNames") >> fieldNames_;
        }
        else if (coeffs_.found("U"))
        {
            word UName(coeffs_.lookup("U"));
            fieldNames_ = wordList(1, UName);
        }
        else
        {
            fieldNames_ = wordList(1, "U");
        }

        applied_.setSize(fieldNames_.size(), false);

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
