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

#include "flameletVolumeFractionSource.H"
#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvcDiv.H"
#include "fvmLaplacian.H"
#include "surfaceInterpolate.H"
#include "flameletTurbulentFluidThermoModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(flameletVolumeFractionSource, 0);
    addToRunTimeSelectionTable
    (
        option,
        flameletVolumeFractionSource,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::volScalarField& Foam::fv::flameletVolumeFractionSource::alpha() const
{
    const word alphaName = IOobject::groupName("alpha", phaseName_);

    if (!mesh_.foundObject<volScalarField>(alphaName))
    {
        volScalarField* alphaPtr =
            new volScalarField
            (
                IOobject
                (
                    alphaName,
                    mesh_.time().constant(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh_
            );

        alphaPtr->store();
    }

    return mesh_.lookupObject<volScalarField>(alphaName);
}


Foam::tmp<Foam::volScalarField> Foam::fv::flameletVolumeFractionSource::D
(
    const label fieldi
) const
{
    const surfaceScalarField& phi =
        mesh().lookupObject<surfaceScalarField>(phiName_);

    if (phi.dimensions() == dimVolume/dimTime)
    {
        const flameletTurbulenceModel& turbulence =
            mesh().lookupObject<flameletTurbulenceModel>
            (
                flameletTurbulenceModel::propertiesName
            );

        return turbulence.nuEff();
    }
    else if (phi.dimensions() == dimMass/dimTime)
    {
        const compressible::flameletTurbulenceModel& turbulence =
            mesh().lookupObject<compressible::flameletTurbulenceModel>
            (
                flameletTurbulenceModel::propertiesName
            );

        return
            fieldNames_[fieldi] == turbulence.transport().T().name()
          ? turbulence.kappaEff()
          : fieldNames_[fieldi] == turbulence.transport().he().name()
          ? turbulence.alphaEff()
          : turbulence.muEff();
    }
    else
    {
        FatalErrorInFunction
            << "Dimensions of " << phi.name() << " not recognised"
            << exit(FatalError);
        return tmp<volScalarField>(nullptr);
    }
}


template <class Type>
void Foam::fv::flameletVolumeFractionSource::addDivSup
(
    fvMatrix<Type>& eqn,
    const label fieldi
)
{
    const surfaceScalarField& phi =
        mesh().lookupObject<surfaceScalarField>(phiName_);

    const volScalarField AByB(this->alpha()/(1 - this->alpha()));

    eqn -= AByB*fvm::div(phi, eqn.psi());
}


void Foam::fv::flameletVolumeFractionSource::addUDivSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    const surfaceScalarField& phi =
        mesh().lookupObject<surfaceScalarField>(phiName_);

    const volScalarField AByB(this->alpha()/(1 - this->alpha()));

    const word scheme("div(" + phiName_ + "," + eqn.psi().name() + ")");

    eqn -= fvm::div(fvc::interpolate(AByB)*phi, eqn.psi(), scheme);
}


void Foam::fv::flameletVolumeFractionSource::addRhoDivSup
(
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    const surfaceScalarField& phi =
        mesh().lookupObject<surfaceScalarField>(phiName_);

    const volScalarField AByB(this->alpha()/(1 - this->alpha()));

    eqn -= AByB*fvc::div(phi);
}


template <class Type, class AlphaFieldType>
void Foam::fv::flameletVolumeFractionSource::addLaplacianSup
(
    const AlphaFieldType& alpha,
    fvMatrix<Type>& eqn,
    const label fieldi
)
{
    const volScalarField B(1 - this->alpha());

    const volScalarField D(this->D(fieldi));

    const word scheme("laplacian(" + D.name() + "," + eqn.psi().name() + ")");

    eqn +=
        fvm::laplacian(D, eqn.psi())
      - 1/B*fvm::laplacian(B*D, eqn.psi(), scheme);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::flameletVolumeFractionSource::flameletVolumeFractionSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    option(name, modelType, dict, mesh),
    phaseName_(dict.lookupType<word>("phase")),
    phiName_("phi"),
    rhoName_("rho"),
    UName_("U")
{
    read(dict);
    alpha();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::flameletVolumeFractionSource::~flameletVolumeFractionSource()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::flameletVolumeFractionSource::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    if (fieldNames_[fieldi] == rhoName_)
    {
        addRhoDivSup(eqn, fieldi);
    }
    else
    {
        addDivSup(eqn, fieldi);
        addLaplacianSup(geometricOneField(), eqn, fieldi);
    }
}


void Foam::fv::flameletVolumeFractionSource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    if (fieldNames_[fieldi] == UName_)
    {
        addUDivSup(eqn, fieldi);
    }
    else
    {
        addDivSup(eqn, fieldi);
        addLaplacianSup(geometricOneField(), eqn, fieldi);
    }
}


void Foam::fv::flameletVolumeFractionSource::addSup
(
    fvMatrix<sphericalTensor>& eqn,
    const label fieldi
)
{
    addDivSup(eqn, fieldi);
    addLaplacianSup(geometricOneField(), eqn, fieldi);
}


void Foam::fv::flameletVolumeFractionSource::addSup
(
    fvMatrix<symmTensor>& eqn,
    const label fieldi
)
{
    addDivSup(eqn, fieldi);
    addLaplacianSup(geometricOneField(), eqn, fieldi);
}


void Foam::fv::flameletVolumeFractionSource::addSup
(
    fvMatrix<tensor>& eqn,
    const label fieldi
)
{
    addDivSup(eqn, fieldi);
    addLaplacianSup(geometricOneField(), eqn, fieldi);
}


void Foam::fv::flameletVolumeFractionSource::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    if (fieldNames_[fieldi] == rhoName_)
    {
        addRhoDivSup(eqn, fieldi);
    }
    else
    {
        addDivSup(eqn, fieldi);
        addLaplacianSup(geometricOneField(), eqn, fieldi);
    }
}


void Foam::fv::flameletVolumeFractionSource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    if (fieldNames_[fieldi] == UName_)
    {
        addUDivSup(eqn, fieldi);
    }
    else
    {
        addDivSup(eqn, fieldi);
        addLaplacianSup(geometricOneField(), eqn, fieldi);
    }
}


void Foam::fv::flameletVolumeFractionSource::addSup
(
    const volScalarField& rho,
    fvMatrix<sphericalTensor>& eqn,
    const label fieldi
)
{
    addDivSup(eqn, fieldi);
    addLaplacianSup(geometricOneField(), eqn, fieldi);
}


void Foam::fv::flameletVolumeFractionSource::addSup
(
    const volScalarField& rho,
    fvMatrix<symmTensor>& eqn,
    const label fieldi
)
{
    addDivSup(eqn, fieldi);
    addLaplacianSup(geometricOneField(), eqn, fieldi);
}


void Foam::fv::flameletVolumeFractionSource::addSup
(
    const volScalarField& rho,
    fvMatrix<tensor>& eqn,
    const label fieldi
)
{
    addDivSup(eqn, fieldi);
    addLaplacianSup(geometricOneField(), eqn, fieldi);
}


void Foam::fv::flameletVolumeFractionSource::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    if (fieldNames_[fieldi] == rhoName_)
    {
        addRhoDivSup(eqn, fieldi);
    }
    else
    {
        addDivSup(eqn, fieldi);
        addLaplacianSup(alpha, eqn, fieldi);
    }
}


void Foam::fv::flameletVolumeFractionSource::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    if (fieldNames_[fieldi] == UName_)
    {
        addUDivSup(eqn, fieldi);
    }
    else
    {
        addDivSup(eqn, fieldi);
        addLaplacianSup(alpha, eqn, fieldi);
    }
}


void Foam::fv::flameletVolumeFractionSource::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<sphericalTensor>& eqn,
    const label fieldi
)
{
    addDivSup(eqn, fieldi);
    addLaplacianSup(alpha, eqn, fieldi);
}


void Foam::fv::flameletVolumeFractionSource::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<symmTensor>& eqn,
    const label fieldi
)
{
    addDivSup(eqn, fieldi);
    addLaplacianSup(alpha, eqn, fieldi);
}


void Foam::fv::flameletVolumeFractionSource::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<tensor>& eqn,
    const label fieldi
)
{
    addDivSup(eqn, fieldi);
    addLaplacianSup(alpha, eqn, fieldi);
}


bool Foam::fv::flameletVolumeFractionSource::read(const dictionary& dict)
{
    if (option::read(dict))
    {
        if (coeffs_.found("fields"))
        {
            coeffs_.lookup("fields") >> fieldNames_;
        }

        applied_.setSize(fieldNames_.size(), false);

        dict.readIfPresent("phi", phiName_);

        dict.readIfPresent("rho", rhoName_);

        dict.readIfPresent("U", UName_);

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
