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

#include "flameletNutLowReWallFunctionFvPatchScalarField.H"
#include "flameletTurbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<scalarField> flameletNutLowReWallFunctionFvPatchScalarField::nut() const
{
    return tmp<scalarField>(new scalarField(patch().size(), 0.0));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

flameletNutLowReWallFunctionFvPatchScalarField::flameletNutLowReWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    flameletNutWallFunctionFvPatchScalarField(p, iF)
{}


flameletNutLowReWallFunctionFvPatchScalarField::flameletNutLowReWallFunctionFvPatchScalarField
(
    const flameletNutLowReWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    flameletNutWallFunctionFvPatchScalarField(ptf, p, iF, mapper)
{}


flameletNutLowReWallFunctionFvPatchScalarField::flameletNutLowReWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    flameletNutWallFunctionFvPatchScalarField(p, iF, dict)
{}


flameletNutLowReWallFunctionFvPatchScalarField::flameletNutLowReWallFunctionFvPatchScalarField
(
    const flameletNutLowReWallFunctionFvPatchScalarField& nlrwfpsf
)
:
    flameletNutWallFunctionFvPatchScalarField(nlrwfpsf)
{}


flameletNutLowReWallFunctionFvPatchScalarField::flameletNutLowReWallFunctionFvPatchScalarField
(
    const flameletNutLowReWallFunctionFvPatchScalarField& nlrwfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    flameletNutWallFunctionFvPatchScalarField(nlrwfpsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<scalarField> flameletNutLowReWallFunctionFvPatchScalarField::yPlus() const
{
    const label patchi = patch().index();
    const flameletTurbulenceModel& turbModel = db().lookupObject<flameletTurbulenceModel>
    (
        IOobject::groupName
        (
            flameletTurbulenceModel::propertiesName,
            internalField().group()
        )
    );
    const scalarField& y = turbModel.y()[patchi];
    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();
    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];

    return y*sqrt(nuw*mag(Uw.snGrad()))/nuw;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    flameletNutLowReWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
