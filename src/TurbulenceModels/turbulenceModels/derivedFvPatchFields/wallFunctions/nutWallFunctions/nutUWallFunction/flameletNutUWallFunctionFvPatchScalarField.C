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

#include "flameletNutUWallFunctionFvPatchScalarField.H"
#include "flameletTurbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<scalarField> flameletNutUWallFunctionFvPatchScalarField::nut() const
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
    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];
    const scalarField magUp(mag(Uw.patchInternalField() - Uw));

    const scalarField yPlus(this->yPlus(magUp));

    tmp<scalarField> tnutw(new scalarField(patch().size(), 0.0));
    scalarField& nutw = tnutw.ref();

    forAll(yPlus, facei)
    {
        if (yPlus[facei] > yPlusLam_)
        {
            nutw[facei] =
                nuw[facei]*(yPlus[facei]*kappa_/log(E_*yPlus[facei]) - 1);
        }
    }

    return tnutw;
}


tmp<scalarField> flameletNutUWallFunctionFvPatchScalarField::yPlus
(
    const scalarField& magUp
) const
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

    tmp<scalarField> tyPlus(new scalarField(patch().size(), 0.0));
    scalarField& yPlus = tyPlus.ref();

    forAll(yPlus, facei)
    {
        const scalar Re = magUp[facei]*y[facei]/nuw[facei];
        const scalar ryPlusLam = 1/yPlusLam_;

        int iter = 0;
        scalar yp = yPlusLam_;
        scalar yPlusLast = yp;

        do
        {
            yPlusLast = yp;
            if (yp > yPlusLam_)
            {
                yp = (kappa_*Re + yp)/(1 + log(E_*yp));
            }
            else
            {
                yp = sqrt(Re);
            }
        } while(mag(ryPlusLam*(yp - yPlusLast)) > 0.0001 && ++iter < 20);

        yPlus[facei] = yp;
    }

    return tyPlus;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

flameletNutUWallFunctionFvPatchScalarField::flameletNutUWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    flameletNutWallFunctionFvPatchScalarField(p, iF)
{}


flameletNutUWallFunctionFvPatchScalarField::flameletNutUWallFunctionFvPatchScalarField
(
    const flameletNutUWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    flameletNutWallFunctionFvPatchScalarField(ptf, p, iF, mapper)
{}


flameletNutUWallFunctionFvPatchScalarField::flameletNutUWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    flameletNutWallFunctionFvPatchScalarField(p, iF, dict)
{}


flameletNutUWallFunctionFvPatchScalarField::flameletNutUWallFunctionFvPatchScalarField
(
    const flameletNutUWallFunctionFvPatchScalarField& sawfpsf
)
:
    flameletNutWallFunctionFvPatchScalarField(sawfpsf)
{}


flameletNutUWallFunctionFvPatchScalarField::flameletNutUWallFunctionFvPatchScalarField
(
    const flameletNutUWallFunctionFvPatchScalarField& sawfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    flameletNutWallFunctionFvPatchScalarField(sawfpsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<scalarField> flameletNutUWallFunctionFvPatchScalarField::yPlus() const
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
    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];
    const scalarField magUp(mag(Uw.patchInternalField() - Uw));

    return yPlus(magUp);
}


void flameletNutUWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    writeLocalEntries(os);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    flameletNutUWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
