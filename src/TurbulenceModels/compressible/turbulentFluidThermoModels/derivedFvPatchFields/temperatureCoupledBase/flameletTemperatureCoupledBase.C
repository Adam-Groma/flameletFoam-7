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

#include "flameletTemperatureCoupledBase.H"
#include "volFields.H"
#include "flameletFluidThermo.H"
//#include "solidThermo.H"
#include "flameletTurbulentFluidThermoModel.H"

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

namespace Foam
{
    template<>
    const char* Foam::NamedEnum
    <
        Foam::flameletTemperatureCoupledBase::KMethodType,
        4
    >::names[] =
    {
        "flameletFluidThermo",
        "solidThermo",
        "directionalSolidThermo",
        "lookup"
    };
}


const Foam::NamedEnum<Foam::flameletTemperatureCoupledBase::KMethodType, 4>
    Foam::flameletTemperatureCoupledBase::KMethodTypeNames_;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::flameletTemperatureCoupledBase::flameletTemperatureCoupledBase
(
    const fvPatch& patch,
    const word& calculationType,
    const word& kappaName
)
:
    patch_(patch),
    method_(KMethodTypeNames_[calculationType]),
    kappaName_(kappaName)
{}


Foam::flameletTemperatureCoupledBase::flameletTemperatureCoupledBase
(
    const fvPatch& patch,
    const dictionary& dict
)
:
    patch_(patch),
    method_(KMethodTypeNames_.read(dict.lookup("kappaMethod"))),
    kappaName_(dict.lookupOrDefault<word>("kappa", "none"))
{
    switch (method_)
    {
/*----------------------------------------------------------\
        case mtDirectionalSolidThermo:
        {
            if (!dict.found("alphaAni"))
            {
                FatalIOErrorInFunction(dict)
                    << "Did not find entry 'alphaAni'"
                       " required for 'kappaMethod' "
                    << KMethodTypeNames_[method_]
                    << exit(FatalIOError);
            }

            break;
        }
\---------------------------------------------------------*/

        case mtLookup:
        {
            if (!dict.found("kappa"))
            {
                FatalIOErrorInFunction(dict)
                    << "Did not find entry 'kappa'"
                       " required for 'kappaMethod' "
                    <<  KMethodTypeNames_[method_] << nl
                    << "    Please set 'kappa' to the name of a volScalarField"
                       " or volSymmTensorField"
                    << exit(FatalIOError);
            }

            break;
        }

        default:
            break;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::flameletTemperatureCoupledBase::kappa
(
    const scalarField& Tp
) const
{
    const fvMesh& mesh = patch_.boundaryMesh().mesh();
    const label patchi = patch_.index();

    switch (method_)
    {
        case mtFluidThermo:
        {
            typedef compressible::flameletTurbulenceModel flameletTurbulenceModel;

            word turbName(flameletTurbulenceModel::propertiesName);

            if
            (
                mesh.foundObject<flameletTurbulenceModel>(turbName)
            )
            {
                const flameletTurbulenceModel& turbModel =
                    mesh.lookupObject<flameletTurbulenceModel>(turbName);

                return turbModel.kappaEff(patchi);
            }
            else if (mesh.foundObject<flameletFluidThermo>(flameletBasicThermo::dictName))
            {
                const flameletFluidThermo& thermo =
                    mesh.lookupObject<flameletFluidThermo>(flameletBasicThermo::dictName);

                return thermo.kappa(patchi);
            }
            else
            {
                FatalErrorInFunction
                    << "kappaMethod defined to employ "
                    << KMethodTypeNames_[method_]
                    << " method, but thermo package not available"
                    << exit(FatalError);
            }

            break;
        }
/*-----------------------------------------------------------------------\
        case mtSolidThermo:
        {
            const solidThermo& thermo =
                mesh.lookupObject<solidThermo>(flameletBasicThermo::dictName);

            return thermo.kappa(patchi);
            break;
        }

        case mtDirectionalSolidThermo:
        {
            const solidThermo& thermo =
                mesh.lookupObject<solidThermo>(flameletBasicThermo::dictName);

            const vectorField kappa(thermo.Kappa(patchi));

            tmp<scalarField> tmeanKappa(Tp);
            scalarField& meanKappa = tmeanKappa();
            forAll(meanKappa, i)
            {
                meanKappa[i] = (kappa[i].x() + kappa[i].y() + kappa[i].z())/3.0;
            }

            return meanKappa;
            break;
        }
\----------------------------------------------------------------------------*/

        case mtLookup:
        {
            if (mesh.foundObject<volScalarField>(kappaName_))
            {
                return patch_.lookupPatchField<volScalarField, scalar>
                (
                    kappaName_
                );
            }
            else if (mesh.foundObject<volSymmTensorField>(kappaName_))
            {
                const symmTensorField& KWall =
                    patch_.lookupPatchField<volSymmTensorField, scalar>
                    (
                        kappaName_
                    );

                const vectorField n(patch_.nf());

                return n & KWall & n;
            }
            else
            {
                FatalErrorInFunction
                    << "Did not find field " << kappaName_
                    << " on mesh " << mesh.name() << " patch " << patch_.name()
                    << nl
                    << "    Please set 'kappa' to the name of a volScalarField"
                       " or volSymmTensorField."
                    << exit(FatalError);
            }

            break;
        }

        default:
        {
            FatalErrorInFunction
                << "Unimplemented method " << KMethodTypeNames_[method_] << nl
                << "    Please set 'kappaMethod' to one of "
                << KMethodTypeNames_.toc()
                << " and 'kappa' to the name of the volScalar"
                << " or volSymmTensor field (if kappa=lookup)"
                << exit(FatalError);
        }
    }

    return scalarField(0);
}


void Foam::flameletTemperatureCoupledBase::write(Ostream& os) const
{
    writeEntry(os, "kappaMethod", KMethodTypeNames_[method_]);
    writeEntry(os, "kappa", kappaName_);
}


// ************************************************************************* //
