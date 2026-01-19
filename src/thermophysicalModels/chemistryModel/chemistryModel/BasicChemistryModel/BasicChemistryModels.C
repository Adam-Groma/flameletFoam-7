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

#include "makeChemistryModel.H"

#include "rhoReactionThermo.H"

#include "StandardChemistryModel.H"
#include "thermoPhysicsTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // Make base types
    makeChemistryModel(rhoReactionThermo);

    // Chemistry moldels based on sensibleEnthalpy
    makeChemistryModelType
    (
        StandardChemistryModel,
        rhoReactionThermo,
        constGasHThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModel,
        rhoReactionThermo,
        gasHThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModel,
        rhoReactionThermo,
        constIncompressibleGasHThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModel,
        rhoReactionThermo,
        incompressibleGasHThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModel,
        rhoReactionThermo,
        icoPoly8HThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModel,
        rhoReactionThermo,
        constFluidHThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModel,
        rhoReactionThermo,
        constAdiabaticFluidHThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModel,
        rhoReactionThermo,
        constHThermoPhysics
    );



    // Chemistry moldels based on sensibleInternalEnergy
    makeChemistryModelType
    (
        StandardChemistryModel,
        rhoReactionThermo,
        constGasEThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModel,
        rhoReactionThermo,
        gasEThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModel,
        rhoReactionThermo,
        constIncompressibleGasEThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModel,
        rhoReactionThermo,
        incompressibleGasEThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModel,
        rhoReactionThermo,
        icoPoly8EThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModel,
        rhoReactionThermo,
        constFluidEThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModel,
        rhoReactionThermo,
        constAdiabaticFluidEThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModel,
        rhoReactionThermo,
        constEThermoPhysics
    );
}

// ************************************************************************* //
