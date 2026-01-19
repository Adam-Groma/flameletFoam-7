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

#include "externalWallHeatFluxTemperatureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "physicoChemicalConstants.H"

using Foam::constant::physicoChemical::sigma;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    template<>
    const char*
    NamedEnum
    <
        externalWallHeatFluxTemperatureFvPatchScalarField::operationMode,
        3
    >::names[] =
    {
        "power",
        "flux",
        "coefficient"
    };
}

const Foam::NamedEnum
<
    Foam::externalWallHeatFluxTemperatureFvPatchScalarField::operationMode,
    3
> Foam::externalWallHeatFluxTemperatureFvPatchScalarField::operationModeNames;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::externalWallHeatFluxTemperatureFvPatchScalarField::
externalWallHeatFluxTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    temperatureCoupledBase(patch(), "undefined", "undefined-K"),
    mode_(fixedHeatFlux),
    Q_(0),
    Ta_(),
    thicknessLayers_(),
    kappaLayers_()
{
    this->refValue() = 0;
    this->refGrad() = 0;
    this->valueFraction() = 1;
}


Foam::externalWallHeatFluxTemperatureFvPatchScalarField::
externalWallHeatFluxTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    temperatureCoupledBase(patch(), dict),
    mode_(operationModeNames.read(dict.lookup("mode"))),
    Q_(0),
    Ta_(),
    thicknessLayers_(),
    kappaLayers_()
{
    switch (mode_)
    {
        case fixedPower:
        {
            dict.lookup("Q") >> Q_;

            break;
        }
        case fixedHeatFlux:
        {
            q_ = scalarField("q", dict, p.size());

            break;
        }
        case fixedHeatTransferCoeff:
        {
            h_ = scalarField("h", dict, p.size());
            Ta_ = Function1<scalar>::New("Ta", dict);

            if (dict.found("thicknessLayers"))
            {
                dict.lookup("thicknessLayers") >> thicknessLayers_;
                dict.lookup("kappaLayers") >> kappaLayers_;
            }

            break;
        }
    }

    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    if (dict.found("refValue"))
    {
        // Full restart
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        // Start from user entered data. Assume fixedValue.
        refValue() = *this;
        refGrad() = 0;
        valueFraction() = 1;
    }
}


Foam::externalWallHeatFluxTemperatureFvPatchScalarField::
externalWallHeatFluxTemperatureFvPatchScalarField
(
    const externalWallHeatFluxTemperatureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    temperatureCoupledBase(patch(), ptf.KMethod(), ptf.kappaName()),
    mode_(ptf.mode_),
    Q_(ptf.Q_),
    Ta_(ptf.Ta_, false),
    thicknessLayers_(ptf.thicknessLayers_),
    kappaLayers_(ptf.kappaLayers_)
{
    switch (mode_)
    {
        case fixedPower:
        {
            break;
        }
        case fixedHeatFlux:
        {
            mapper(q_, ptf.q_);
            break;
        }
        case fixedHeatTransferCoeff:
        {
            mapper(h_, ptf.h_);
            break;
        }
    }
}


Foam::externalWallHeatFluxTemperatureFvPatchScalarField::
externalWallHeatFluxTemperatureFvPatchScalarField
(
    const externalWallHeatFluxTemperatureFvPatchScalarField& tppsf
)
:
    mixedFvPatchScalarField(tppsf),
    temperatureCoupledBase(tppsf),
    mode_(tppsf.mode_),
    Q_(tppsf.Q_),
    q_(tppsf.q_),
    h_(tppsf.h_),
    Ta_(tppsf.Ta_, false),
    thicknessLayers_(tppsf.thicknessLayers_),
    kappaLayers_(tppsf.kappaLayers_)
{}


Foam::externalWallHeatFluxTemperatureFvPatchScalarField::
externalWallHeatFluxTemperatureFvPatchScalarField
(
    const externalWallHeatFluxTemperatureFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(tppsf, iF),
    temperatureCoupledBase(patch(), tppsf.KMethod(), tppsf.kappaName()),
    mode_(tppsf.mode_),
    Q_(tppsf.Q_),
    q_(tppsf.q_),
    h_(tppsf.h_),
    Ta_(tppsf.Ta_, false),
    thicknessLayers_(tppsf.thicknessLayers_),
    kappaLayers_(tppsf.kappaLayers_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::externalWallHeatFluxTemperatureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchScalarField::autoMap(m);

    switch (mode_)
    {
        case fixedPower:
        {
            break;
        }
        case fixedHeatFlux:
        {
            m(q_, q_);

            break;
        }
        case fixedHeatTransferCoeff:
        {
            m(h_, h_);

            break;
        }
    }
}


void Foam::externalWallHeatFluxTemperatureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);

    const externalWallHeatFluxTemperatureFvPatchScalarField& tiptf =
        refCast<const externalWallHeatFluxTemperatureFvPatchScalarField>(ptf);

    switch (mode_)
    {
        case fixedPower:
        {
            break;
        }
        case fixedHeatFlux:
        {
            q_.rmap(tiptf.q_, addr);

            break;
        }
        case fixedHeatTransferCoeff:
        {
            h_.rmap(tiptf.h_, addr);

            break;
        }
    }
}


void Foam::externalWallHeatFluxTemperatureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    scalarField q(size(), 0.0);
    const scalarField Tc(patchInternalField());
    const scalarField& Tp(*this);
    const scalarField KWall(kappa(Tp));
    const scalarField KDelta(KWall*patch().deltaCoeffs());

    switch (mode_)
    {
        case fixedPower:
        {
            q = Q_/gSum(patch().magSf());
            break;
        }
        case fixedHeatFlux:
        {
            q = q_;
            break;
        }
        case fixedHeatTransferCoeff:
        {
            scalar totalSolidRes = 0;
            if (thicknessLayers_.size())
            {
                forAll(thicknessLayers_, iLayer)
                {
                    const scalar l = thicknessLayers_[iLayer];
                    if (kappaLayers_[iLayer] > 0)
                    {
                        totalSolidRes += l/kappaLayers_[iLayer];
                    }
                }
            }
            const scalar Ta = Ta_->value(this->db().time().timeOutputValue());

            q = (Ta - Tp)*(1.0/h_ + totalSolidRes);
            break;
        }
    }

    forAll(*this, i)
    {
        if (q[i] > 0) //in
        {
            this->refGrad()[i] = q[i]/KWall[i];
            this->refValue()[i] = 0.0;
            this->valueFraction()[i] = 0.0;
        }
        else //out
        {
            this->refGrad()[i] = 0.0;
            this->refValue()[i] = q[i]/KDelta[i] + Tc[i];
            this->valueFraction()[i] = 1.0;
        }
    }

    mixedFvPatchScalarField::updateCoeffs();

    if (debug)
    {
        const scalar Q = gSum(KWall*patch().magSf()*snGrad());

        Info<< patch().boundaryMesh().mesh().name() << ':'
            << patch().name() << ':'
            << this->internalField().name() << " :"
            << " heat transfer rate:" << Q
            << " walltemperature "
            << " min:" << gMin(*this)
            << " max:" << gMax(*this)
            << " avg:" << gAverage(*this)
            << endl;
    }
}


void Foam::externalWallHeatFluxTemperatureFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);

    writeEntry(os, "mode", operationModeNames[mode_]);
    temperatureCoupledBase::write(os);

    switch (mode_)
    {
        case fixedPower:
        {
            writeEntry(os, "Q", Q_);

            break;
        }
        case fixedHeatFlux:
        {
            writeEntry(os, "q", q_);

            break;
        }
        case fixedHeatTransferCoeff:
        {
            writeEntry(os, "h", h_);
            writeEntry(os, Ta_());

            if (thicknessLayers_.size())
            {
                writeEntry(os, "thicknessLayers", thicknessLayers_);
                writeEntry(os, "kappaLayers", kappaLayers_);
            }

            break;
        }
    }

    writeEntry(os, "refValue", refValue());
    writeEntry(os, "refGradient", refGrad());
    writeEntry(os, "valueFraction", valueFraction());
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        externalWallHeatFluxTemperatureFvPatchScalarField
    );
}

// ************************************************************************* //
