/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

InClass
    Foam::psiChemistryModel

Description
    Creates chemistry model instances templated on the type of thermodynamics

\*---------------------------------------------------------------------------*/

#include "flameletMakeChemistryModel.H"

#include "flameletRhoReactionThermo.H"

#include "FlameletStandardChemistryModel.H"
#include "flameletThermoPhysicsTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // Make base types
    flameletMakeChemistryModel(flameletRhoReactionThermo);

    // Chemistry moldels based on sensibleEnthalpy
    flameletMakeChemistryModelType
    (
        FlameletStandardChemistryModel,
        flameletRhoReactionThermo,
        constGasHThermoPhysics
    );

    flameletMakeChemistryModelType
    (
        FlameletStandardChemistryModel,
        flameletRhoReactionThermo,
        gasHThermoPhysics
    );

    flameletMakeChemistryModelType
    (
        FlameletStandardChemistryModel,
        flameletRhoReactionThermo,
        constIncompressibleGasHThermoPhysics
    );

    flameletMakeChemistryModelType
    (
        FlameletStandardChemistryModel,
        flameletRhoReactionThermo,
        incompressibleGasHThermoPhysics
    );

    flameletMakeChemistryModelType
    (
        FlameletStandardChemistryModel,
        flameletRhoReactionThermo,
        icoPoly8HThermoPhysics
    );

    flameletMakeChemistryModelType
    (
        FlameletStandardChemistryModel,
        flameletRhoReactionThermo,
        constFluidHThermoPhysics
    );

    flameletMakeChemistryModelType
    (
        FlameletStandardChemistryModel,
        flameletRhoReactionThermo,
        constAdiabaticFluidHThermoPhysics
    );

    flameletMakeChemistryModelType
    (
        FlameletStandardChemistryModel,
        flameletRhoReactionThermo,
        constHThermoPhysics
    );



    // Chemistry moldels based on sensibleInternalEnergy
    flameletMakeChemistryModelType
    (
        FlameletStandardChemistryModel,
        flameletRhoReactionThermo,
        constGasEThermoPhysics
    );

    flameletMakeChemistryModelType
    (
        FlameletStandardChemistryModel,
        flameletRhoReactionThermo,
        gasEThermoPhysics
    );

    flameletMakeChemistryModelType
    (
        FlameletStandardChemistryModel,
        flameletRhoReactionThermo,
        constIncompressibleGasEThermoPhysics
    );

    flameletMakeChemistryModelType
    (
        FlameletStandardChemistryModel,
        flameletRhoReactionThermo,
        incompressibleGasEThermoPhysics
    );

    flameletMakeChemistryModelType
    (
        FlameletStandardChemistryModel,
        flameletRhoReactionThermo,
        icoPoly8EThermoPhysics
    );

    flameletMakeChemistryModelType
    (
        FlameletStandardChemistryModel,
        flameletRhoReactionThermo,
        constFluidEThermoPhysics
    );

    flameletMakeChemistryModelType
    (
        FlameletStandardChemistryModel,
        flameletRhoReactionThermo,
        constAdiabaticFluidEThermoPhysics
    );

    flameletMakeChemistryModelType
    (
        FlameletStandardChemistryModel,
        flameletRhoReactionThermo,
        constEThermoPhysics
    );
}

// ************************************************************************* //
