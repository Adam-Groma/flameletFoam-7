/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2019 OpenFOAM Foundation
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

#include "addToRunTimeSelectionTable.H"
#include "flameletEnergyJumpAMIFvPatchScalarField.H"
#include "fixedJumpAMIFvPatchFields.H"
#include "flameletBasicThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::flameletEnergyJumpAMIFvPatchScalarField::flameletEnergyJumpAMIFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedJumpAMIFvPatchField<scalar>(p, iF)
{}


Foam::flameletEnergyJumpAMIFvPatchScalarField::flameletEnergyJumpAMIFvPatchScalarField
(
    const flameletEnergyJumpAMIFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedJumpAMIFvPatchField<scalar>(ptf, p, iF, mapper)
{}


Foam::flameletEnergyJumpAMIFvPatchScalarField::flameletEnergyJumpAMIFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedJumpAMIFvPatchField<scalar>(p, iF)
{
    if (dict.found("value"))
    {
        fvPatchScalarField::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        evaluate(Pstream::commsTypes::blocking);
    }
}


Foam::flameletEnergyJumpAMIFvPatchScalarField::flameletEnergyJumpAMIFvPatchScalarField
(
    const flameletEnergyJumpAMIFvPatchScalarField& ptf
)
:
    fixedJumpAMIFvPatchField<scalar>(ptf)
{}


Foam::flameletEnergyJumpAMIFvPatchScalarField::flameletEnergyJumpAMIFvPatchScalarField
(
    const flameletEnergyJumpAMIFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedJumpAMIFvPatchField<scalar>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::flameletEnergyJumpAMIFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (this->cyclicAMIPatch().owner())
    {
        const flameletBasicThermo& thermo = flameletBasicThermo::lookupThermo(*this);
        label patchID = patch().index();

        const scalarField& pp = thermo.p().boundaryField()[patchID];
        const fixedJumpAMIFvPatchScalarField& TbPatch =
            refCast<const fixedJumpAMIFvPatchScalarField>
            (
                thermo.T().boundaryField()[patchID]
            );

        fixedJumpAMIFvPatchScalarField& Tbp =
            const_cast<fixedJumpAMIFvPatchScalarField&>(TbPatch);

        // force update of jump
        Tbp.updateCoeffs();

        const labelUList& faceCells = this->patch().faceCells();

        jump_ = thermo.he(pp, Tbp.jump(), faceCells);
    }

    fixedJumpAMIFvPatchField<scalar>::updateCoeffs();
}


void Foam::flameletEnergyJumpAMIFvPatchScalarField::write(Ostream& os) const
{
    fixedJumpAMIFvPatchField<scalar>::write(os);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchScalarField,
       flameletEnergyJumpAMIFvPatchScalarField
   );
}

// ************************************************************************* //
