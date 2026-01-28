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

#include "flameletActuationDiskSource.H"
#include "fvMesh.H"
#include "fvMatrix.H"
#include "geometricOneField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(flameletActuationDiskSource, 0);
    addToRunTimeSelectionTable
    (
        option,
        flameletActuationDiskSource,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::flameletActuationDiskSource::checkData() const
{
    if (magSqr(diskArea_) <= vSmall)
    {
        FatalErrorInFunction
           << "diskArea is approximately zero"
           << exit(FatalIOError);
    }
    if (Cp_ <= vSmall || Ct_ <= vSmall)
    {
        FatalErrorInFunction
           << "Cp and Ct must be greater than zero"
           << exit(FatalIOError);
    }
    if (mag(diskDir_) < vSmall)
    {
        FatalErrorInFunction
           << "disk direction vector is approximately zero"
           << exit(FatalIOError);
    }
    if (returnReduce(upstreamCellId_, maxOp<label>()) == -1)
    {
        FatalErrorInFunction
           << "upstream location " << upstreamPoint_  << " not found in mesh"
           << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::flameletActuationDiskSource::flameletActuationDiskSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    flameletCellSetOption(name, modelType, dict, mesh),
    diskDir_(coeffs_.lookup("diskDir")),
    Cp_(readScalar(coeffs_.lookup("Cp"))),
    Ct_(readScalar(coeffs_.lookup("Ct"))),
    diskArea_(readScalar(coeffs_.lookup("diskArea"))),
    upstreamPoint_(coeffs_.lookup("upstreamPoint")),
    upstreamCellId_(-1)
{
    coeffs_.lookup("fields") >> fieldNames_;
    applied_.setSize(fieldNames_.size(), false);

    Info<< "    - creating actuation disk zone: "
        << this->name() << endl;

    upstreamCellId_ = mesh.findCell(upstreamPoint_);

    checkData();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::flameletActuationDiskSource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    const scalarField& cellsV = mesh_.V();
    vectorField& Usource = eqn.source();
    const vectorField& U = eqn.psi();

    if (V() > vSmall)
    {
        addActuationDiskAxialInertialResistance
        (
            Usource,
            cells_,
            cellsV,
            geometricOneField(),
            U
        );
    }
}


void Foam::fv::flameletActuationDiskSource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    const scalarField& cellsV = mesh_.V();
    vectorField& Usource = eqn.source();
    const vectorField& U = eqn.psi();

    if (V() > vSmall)
    {
        addActuationDiskAxialInertialResistance
        (
            Usource,
            cells_,
            cellsV,
            rho,
            U
        );
    }
}


bool Foam::fv::flameletActuationDiskSource::read(const dictionary& dict)
{
    if (flameletCellSetOption::read(dict))
    {
        coeffs_.readIfPresent("diskDir", diskDir_);
        coeffs_.readIfPresent("Cp", Cp_);
        coeffs_.readIfPresent("Ct", Ct_);
        coeffs_.readIfPresent("diskArea", diskArea_);

        checkData();

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
