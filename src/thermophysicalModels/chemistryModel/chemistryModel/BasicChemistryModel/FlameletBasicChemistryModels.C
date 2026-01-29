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
        flameletConstGasHThermoPhysics
    );

    flameletMakeChemistryModelType
    (
        FlameletStandardChemistryModel,
        flameletRhoReactionThermo,
        flameletGasHThermoPhysics
    );

    flameletMakeChemistryModelType
    (
        FlameletStandardChemistryModel,
        flameletRhoReactionThermo,
        flameletConstIncompressibleGasHThermoPhysics
    );

    flameletMakeChemistryModelType
    (
        FlameletStandardChemistryModel,
        flameletRhoReactionThermo,
        flameletIncompressibleGasHThermoPhysics
    );

    flameletMakeChemistryModelType
    (
        FlameletStandardChemistryModel,
        flameletRhoReactionThermo,
        flameletIcoPoly8HThermoPhysics
    );

    flameletMakeChemistryModelType
    (
        FlameletStandardChemistryModel,
        flameletRhoReactionThermo,
        flameletConstFluidHThermoPhysics
    );

    flameletMakeChemistryModelType
    (
        FlameletStandardChemistryModel,
        flameletRhoReactionThermo,
        flameletConstAdiabaticFluidHThermoPhysics
    );

    flameletMakeChemistryModelType
    (
        FlameletStandardChemistryModel,
        flameletRhoReactionThermo,
        flameletConstHThermoPhysics
    );



    // Chemistry moldels based on sensibleInternalEnergy
    flameletMakeChemistryModelType
    (
        FlameletStandardChemistryModel,
        flameletRhoReactionThermo,
        flameletConstGasEThermoPhysics
    );

    flameletMakeChemistryModelType
    (
        FlameletStandardChemistryModel,
        flameletRhoReactionThermo,
        flameletGasEThermoPhysics
    );

    flameletMakeChemistryModelType
    (
        FlameletStandardChemistryModel,
        flameletRhoReactionThermo,
        flameletConstIncompressibleGasEThermoPhysics
    );

    flameletMakeChemistryModelType
    (
        FlameletStandardChemistryModel,
        flameletRhoReactionThermo,
        flameletIncompressibleGasEThermoPhysics
    );

    flameletMakeChemistryModelType
    (
        FlameletStandardChemistryModel,
        flameletRhoReactionThermo,
        flameletIcoPoly8EThermoPhysics
    );

    flameletMakeChemistryModelType
    (
        FlameletStandardChemistryModel,
        flameletRhoReactionThermo,
        flameletConstFluidEThermoPhysics
    );

    flameletMakeChemistryModelType
    (
        FlameletStandardChemistryModel,
        flameletRhoReactionThermo,
        flameletConstAdiabaticFluidEThermoPhysics
    );

    flameletMakeChemistryModelType
    (
        FlameletStandardChemistryModel,
        flameletRhoReactionThermo,
        flameletConstEThermoPhysics
    );
}

// ************************************************************************* //
