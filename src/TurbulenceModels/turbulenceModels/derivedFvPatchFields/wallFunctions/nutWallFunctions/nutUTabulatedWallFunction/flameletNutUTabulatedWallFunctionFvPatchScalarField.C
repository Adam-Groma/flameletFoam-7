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

#include "flameletNutUTabulatedWallFunctionFvPatchScalarField.H"
#include "flameletTurbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<scalarField> flameletNutUTabulatedWallFunctionFvPatchScalarField::nut() const
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
    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];
    const scalarField magUp(mag(Uw.patchInternalField() - Uw));
    const scalarField magGradU(mag(Uw.snGrad()));
    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    return
        max
        (
            scalar(0),
            sqr(magUp/(calcUPlus(magUp*y/nuw) + rootVSmall))
           /(magGradU + rootVSmall)
          - nuw
        );
}


tmp<scalarField> flameletNutUTabulatedWallFunctionFvPatchScalarField::calcUPlus
(
    const scalarField& Rey
) const
{
    tmp<scalarField> tuPlus(new scalarField(patch().size(), 0.0));
    scalarField& uPlus = tuPlus.ref();

    forAll(uPlus, facei)
    {
        uPlus[facei] = uPlusTable_.interpolateLog10(Rey[facei]);
    }

    return tuPlus;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

flameletNutUTabulatedWallFunctionFvPatchScalarField::
flameletNutUTabulatedWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    flameletNutWallFunctionFvPatchScalarField(p, iF),
    uPlusTableName_("undefined-uPlusTableName"),
    uPlusTable_
    (
        IOobject
        (
            uPlusTableName_,
            patch().boundaryMesh().mesh().time().constant(),
            patch().boundaryMesh().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        false
    )
{}


flameletNutUTabulatedWallFunctionFvPatchScalarField::
flameletNutUTabulatedWallFunctionFvPatchScalarField
(
    const flameletNutUTabulatedWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    flameletNutWallFunctionFvPatchScalarField(ptf, p, iF, mapper),
    uPlusTableName_(ptf.uPlusTableName_),
    uPlusTable_(ptf.uPlusTable_)
{}


flameletNutUTabulatedWallFunctionFvPatchScalarField::
flameletNutUTabulatedWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    flameletNutWallFunctionFvPatchScalarField(p, iF, dict),
    uPlusTableName_(dict.lookup("uPlusTable")),
    uPlusTable_
    (
        IOobject
        (
            uPlusTableName_,
            patch().boundaryMesh().mesh().time().constant(),
            patch().boundaryMesh().mesh(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false
        ),
        true
    )
{}


flameletNutUTabulatedWallFunctionFvPatchScalarField::
flameletNutUTabulatedWallFunctionFvPatchScalarField
(
    const flameletNutUTabulatedWallFunctionFvPatchScalarField& wfpsf
)
:
    flameletNutWallFunctionFvPatchScalarField(wfpsf),
    uPlusTableName_(wfpsf.uPlusTableName_),
    uPlusTable_(wfpsf.uPlusTable_)
{}


flameletNutUTabulatedWallFunctionFvPatchScalarField::
flameletNutUTabulatedWallFunctionFvPatchScalarField
(
    const flameletNutUTabulatedWallFunctionFvPatchScalarField& wfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    flameletNutWallFunctionFvPatchScalarField(wfpsf, iF),
    uPlusTableName_(wfpsf.uPlusTableName_),
    uPlusTable_(wfpsf.uPlusTable_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<scalarField> flameletNutUTabulatedWallFunctionFvPatchScalarField::yPlus() const
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
    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];
    const scalarField magUp(mag(Uw.patchInternalField() - Uw));
    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();
    const scalarField Rey(magUp*y/nuw);

    return Rey/(calcUPlus(Rey) + rootVSmall);
}


void flameletNutUTabulatedWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    writeEntry(os, "uPlusTable", uPlusTableName_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    flameletNutUTabulatedWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
