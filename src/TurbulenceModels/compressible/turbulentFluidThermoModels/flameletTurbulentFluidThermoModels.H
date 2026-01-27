/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
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

#include "FlameletCompressibleTurbulenceModel.H"
#include "compressibleTransportModel.H"
#include "flameletFluidThermo.H"
#include "addToRunTimeSelectionTable.H"
#include "flameletMakeTurbulenceModel.H"

#include "FlameletThermalDiffusivity.H"
#include "FlameletEddyDiffusivity.H"

#include "flameletLaminarModel.H"
#include "flameletRASModel.H"
#include "flameletLESModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

flameletMakeTurbulenceModelTypes
(
    geometricOneField,
    volScalarField,
    flameletCompressibleTurbulenceModel,
    FlameletCompressibleTurbulenceModel,
    FlameletThermalDiffusivity,
    flameletFluidThermo
);

#define flameletMakeLaminarModel(Type)                                                 \
    flameletMakeTemplatedLaminarModel                                                  \
    (flameletFluidThermoFlameletCompressibleTurbulenceModel, laminar, Type)

#define flameletMakeRASModel(Type)                                                     \
    flameletMakeTemplatedTurbulenceModel                                               \
    (flameletFluidThermoFlameletCompressibleTurbulenceModel, RAS, Type)

#define flameletMakeLESModel(Type)                                                     \
    flameletMakeTemplatedTurbulenceModel                                               \
    (flameletFluidThermoFlameletCompressibleTurbulenceModel, LES, Type)


// ************************************************************************* //
