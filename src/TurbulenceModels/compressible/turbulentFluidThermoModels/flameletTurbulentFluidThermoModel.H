/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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

Typedef
    Foam::compressible::flameletTurbulenceModel

Typedef
    Foam::compressible::flameletRASModel

Typedef
    Foam::compressible::flameletLESModel

Description
    Typedefs for turbulence, RAS and LES models for compressible flow
    based on the standard laminar transport package.

SourceFiles
    flameletTurbulentFluidThermoModel.C
    flameletTurbulentFluidThermoModels.C

\*---------------------------------------------------------------------------*/

#ifndef flameletTurbulentFluidThermoModel_H
#define flameletTurbulentFluidThermoModel_H

#include "FlameletCompressibleTurbulenceModel.H"
#include "FlameletThermalDiffusivity.H"
#include "FlameletEddyDiffusivity.H"
#include "flameletLaminarModel.H"
#include "flameletRASModel.H"
#include "flameletLESModel.H"
#include "flameletFluidThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    namespace compressible
    {
        typedef FlameletThermalDiffusivity<FlameletCompressibleTurbulenceModel<flameletFluidThermo>>
            flameletTurbulenceModel;

        typedef flameletLaminarModel<flameletTurbulenceModel> laminarModel;
        typedef flameletRASModel<FlameletEddyDiffusivity<flameletTurbulenceModel>> flameletRASModel;
        typedef flameletLESModel<FlameletEddyDiffusivity<flameletTurbulenceModel>> flameletLESModel;

        template<class BasicCompressibleTurbulenceModel>
        autoPtr<BasicCompressibleTurbulenceModel> New
        (
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const typename BasicCompressibleTurbulenceModel::transportModel&
                transport,
            const word& propertiesName = flameletTurbulenceModel::propertiesName
        );
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "flameletTurbulentFluidThermoModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
