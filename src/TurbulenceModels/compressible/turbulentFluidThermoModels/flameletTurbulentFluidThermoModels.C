/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2019 OpenFOAM Foundation
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

#include "flameletTurbulentFluidThermoModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

flameletMakeBaseTurbulenceModel
(
    geometricOneField,
    volScalarField,
    flameletCompressibleTurbulenceModel,
    FlameletCompressibleTurbulenceModel,
    FlameletThermalDiffusivity,
    flameletFluidThermo
);


// -------------------------------------------------------------------------- //
// Laminar models
// -------------------------------------------------------------------------- //

//#include "flameletStokes.H"
//flameletMakeLaminarModel(flameletStokes);

///#include "flameletGeneralizedNewtonian.H"
///flameletMakeLaminarModel(flameletGeneralizedNewtonian);
///
///#include "flameletMaxwell.H"
///flameletMakeLaminarModel(flameletMaxwell);
///
///#include "flameletGiesekus.H"
///flameletMakeLaminarModel(flameletGiesekus);


// -------------------------------------------------------------------------- //
// RAS models
// -------------------------------------------------------------------------- //

#include "flameletSpalartAllmaras.H"
flameletMakeRASModel(flameletSpalartAllmaras);

#include "flameletKEpsilon.H"
flameletMakeRASModel(flameletKEpsilon);

#include "flameletRNGkEpsilon.H"
flameletMakeRASModel(flameletRNGkEpsilon);

#include "flameletRealizableKE.H"
flameletMakeRASModel(flameletRealizableKE);

#include "flameletBuoyantKEpsilon.H"
flameletMakeRASModel(flameletBuoyantKEpsilon);

#include "flameletLaunderSharmaKE.H"
flameletMakeRASModel(flameletLaunderSharmaKE);

#include "flameletKOmega.H"
flameletMakeRASModel(flameletKOmega);

#include "flameletKOmegaSST.H"
flameletMakeRASModel(flameletKOmegaSST);

#include "flameletKOmegaSSTSAS.H"
flameletMakeRASModel(flameletKOmegaSSTSAS);

#include "flameletKOmegaSSTLM.H"
flameletMakeRASModel(flameletKOmegaSSTLM);

#include "flameletV2f.H"
flameletMakeRASModel(flameletV2f);

#include "flameletLRR.H"
flameletMakeRASModel(flameletLRR);

#include "flameletSSG.H"
flameletMakeRASModel(flameletSSG);


// -------------------------------------------------------------------------- //
// LES models
// -------------------------------------------------------------------------- //

#include "flameletSmagorinsky.H"
flameletMakeLESModel(flameletSmagorinsky);

#include "flameletWALE.H"
flameletMakeLESModel(flameletWALE);

#include "flameletKEqn.H"
flameletMakeLESModel(flameletKEqn);

#include "flameletDynamicKEqn.H"
flameletMakeLESModel(flameletDynamicKEqn);

#include "flameletDynamicLagrangian.H"
flameletMakeLESModel(flameletDynamicLagrangian);

#include "flameletKOmegaSSTDES.H"
flameletMakeLESModel(flameletKOmegaSSTDES);

#include "flameletSpalartAllmarasDES.H"
flameletMakeLESModel(flameletSpalartAllmarasDES);

#include "flameletSpalartAllmarasDDES.H"
flameletMakeLESModel(flameletSpalartAllmarasDDES);

#include "flameletSpalartAllmarasIDDES.H"
flameletMakeLESModel(flameletSpalartAllmarasIDDES);

#include "flameletDeardorffDiffStress.H"
flameletMakeLESModel(flameletDeardorffDiffStress);


// ************************************************************************* //
