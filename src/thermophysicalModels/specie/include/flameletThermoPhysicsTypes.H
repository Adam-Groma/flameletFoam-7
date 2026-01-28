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

Typedefs
    Foam::flameletThermoPhysicsTypes

Description
    Type definitions for thermo-physics models

\*---------------------------------------------------------------------------*/

#ifndef flameletThermoPhysicsTypes_H
#define flameletThermoPhysicsTypes_H

#include "flameletSpecie.H"
#include "flameletPerfectGas.H"
#include "flameletIncompressiblePerfectGas.H"
#include "flameletPerfectFluid.H"
#include "flameletAdiabaticPerfectFluid.H"
#include "flameletRhoConst.H"
#include "flameletHConstThermo.H"
#include "flameletEConstThermo.H"
#include "flameletJanafThermo.H"

#include "flameletSensibleEnthalpy.H"
#include "flameletSensibleInternalEnergy.H"
#include "flameletThermo.H"
#include "flameletSutherlandTransport.H"
#include "flameletConstTransport.H"

#include "flameletIcoPolynomial.H"
#include "flameletHPolynomialThermo.H"
#include "flameletPolynomialTransport.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // thermo physics types based on flameletSensibleEnthalpy
    typedef
    flameletConstTransport
    <
        species::flameletThermo
        <
            flameletHConstThermo
            <
                flameletPerfectGas<flameletSpecie>
            >,
            flameletSensibleEnthalpy
        >
    > flameletConstGasHThermoPhysics;

    typedef
    flameletSutherlandTransport
    <
        species::flameletThermo
        <
            flameletJanafThermo
            <
                flameletPerfectGas<flameletSpecie>
            >,
            flameletSensibleEnthalpy
        >
    > flameletGasHThermoPhysics;

    typedef
    flameletConstTransport
    <
        species::flameletThermo
        <
            flameletHConstThermo
            <
                flameletIncompressiblePerfectGas<flameletSpecie>
            >,
            flameletSensibleEnthalpy
        >
    > flameletConstIncompressibleGasHThermoPhysics;

    typedef
    flameletSutherlandTransport
    <
        species::flameletThermo
        <
            flameletJanafThermo
            <
                flameletIncompressiblePerfectGas<flameletSpecie>
            >,
            flameletSensibleEnthalpy
        >
    > flameletIncompressibleGasHThermoPhysics;

    typedef
    flameletPolynomialTransport
    <
        species::flameletThermo
        <
            flameletHPolynomialThermo
            <
                flameletIcoPolynomial<flameletSpecie, 8>,
                8
            >,
            flameletSensibleEnthalpy
        >,
        8
    > flameletIcoPoly8HThermoPhysics;

    typedef
    flameletConstTransport
    <
        species::flameletThermo
        <
            flameletHConstThermo
            <
                flameletPerfectFluid<flameletSpecie>
            >,
            flameletSensibleEnthalpy
        >
    >
    flameletConstFluidHThermoPhysics;

    typedef
    flameletConstTransport
    <
        species::flameletThermo
        <
            flameletHConstThermo
            <
                flameletAdiabaticPerfectFluid<flameletSpecie>
            >,
            flameletSensibleEnthalpy
        >
    >
    flameletConstAdiabaticFluidHThermoPhysics;

    typedef
    flameletConstTransport
    <
        species::flameletThermo
        <
            flameletHConstThermo
            <
                flameletRhoConst<flameletSpecie>
            >,
            flameletSensibleEnthalpy
        >
    >
    flameletConstHThermoPhysics;


    // thermo physics types based on flameletSensibleInternalEnergy
    typedef
    flameletConstTransport
    <
        species::flameletThermo
        <
            flameletEConstThermo
            <
                flameletPerfectGas<flameletSpecie>
            >,
            flameletSensibleInternalEnergy
        >
    > flameletConstGasEThermoPhysics;

    typedef
    flameletSutherlandTransport
    <
        species::flameletThermo
        <
            flameletJanafThermo
            <
                flameletPerfectGas<flameletSpecie>
            >,
            flameletSensibleInternalEnergy
        >
    > flameletGasEThermoPhysics;

    typedef
    flameletConstTransport
    <
        species::flameletThermo
        <
            flameletEConstThermo
            <
                flameletIncompressiblePerfectGas<flameletSpecie>
            >,
            flameletSensibleInternalEnergy
        >
    > flameletConstIncompressibleGasEThermoPhysics;

    typedef
    flameletSutherlandTransport
    <
        species::flameletThermo
        <
            flameletJanafThermo
            <
                flameletIncompressiblePerfectGas<flameletSpecie>
            >,
            flameletSensibleInternalEnergy
        >
    > flameletIncompressibleGasEThermoPhysics;

    typedef
    flameletPolynomialTransport
    <
        species::flameletThermo
        <
            flameletHPolynomialThermo
            <
                flameletIcoPolynomial<flameletSpecie, 8>,
                8
            >,
            flameletSensibleInternalEnergy
        >,
        8
    > flameletIcoPoly8EThermoPhysics;

    typedef
    flameletConstTransport
    <
        species::flameletThermo
        <
            flameletEConstThermo
            <
                flameletPerfectFluid<flameletSpecie>
            >,
            flameletSensibleInternalEnergy
        >
    >
    flameletConstFluidEThermoPhysics;

    typedef
    flameletConstTransport
    <
        species::flameletThermo
        <
            flameletEConstThermo
            <
                flameletAdiabaticPerfectFluid<flameletSpecie>
            >,
            flameletSensibleInternalEnergy
        >
    >
    flameletConstAdiabaticFluidEThermoPhysics;

    typedef
    flameletConstTransport
    <
        species::flameletThermo
        <
            flameletEConstThermo
            <
                flameletRhoConst<flameletSpecie>
            >,
            flameletSensibleInternalEnergy
        >
    >
    flameletConstEThermoPhysics;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
