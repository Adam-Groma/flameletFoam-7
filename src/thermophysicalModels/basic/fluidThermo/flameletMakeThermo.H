/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2018 OpenFOAM Foundation
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
    Foam::flameletFluidThermo

Description
    Macros for creating basic fluid thermo packages

\*---------------------------------------------------------------------------*/

#ifndef flameletMakeThermo_H
#define flameletMakeThermo_H

#include "flameletFluidThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define flameletTypedefThermoPhysics(Transport,Type,Thermo,EqnOfState,Specie)          \
                                                                               \
    typedef                                                                    \
        Transport                                                              \
        <                                                                      \
            species::flameletThermo                                                    \
            <                                                                  \
                Thermo                                                         \
                <                                                              \
                    EqnOfState                                                 \
                    <                                                          \
                        Specie                                                 \
                    >                                                          \
                >,                                                             \
                Type                                                           \
            >                                                                  \
        >                                                                      \
        Transport##Type##Thermo##EqnOfState##Specie


#define flameletDefineThermoPhysicsThermo(BaseThermo,CThermo,Mixture,ThermoPhys)       \
                                                                               \
    typedef                                                                    \
        CThermo                                                                \
        <                                                                      \
            BaseThermo,                                                        \
            Mixture<ThermoPhys>                                                \
        >                                                                      \
        CThermo##Mixture##ThermoPhys;                                          \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        CThermo##Mixture##ThermoPhys,                                          \
        (#CThermo"<" + Mixture<ThermoPhys>::typeName() + ">").c_str(),         \
        0                                                                      \
    )


#define flameletAddThermoPhysicsThermo(BaseThermo,CThermoMixtureThermoPhys)            \
                                                                               \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        BaseThermo,                                                            \
        CThermoMixtureThermoPhys,                                              \
        fvMesh                                                                 \
    );                                                                         \


#define flameletMakeThermoPhysicsThermo(BaseThermo,CThermo,Mixture,ThermoPhys)         \
                                                                               \
    flameletDefineThermoPhysicsThermo(BaseThermo, CThermo, Mixture, ThermoPhys);       \
                                                                               \
    flameletAddThermoPhysicsThermo(BaseThermo, CThermo##Mixture##ThermoPhys)


#define flameletMakeThermoPhysicsThermos(BaseThermo,CThermo,Mixture,ThermoPhys)        \
                                                                               \
    flameletDefineThermoPhysicsThermo(BaseThermo, CThermo, Mixture, ThermoPhys);       \
                                                                               \
    flameletAddThermoPhysicsThermo(flameletBasicThermo, CThermo##Mixture##ThermoPhys);         \
    flameletAddThermoPhysicsThermo(flameletFluidThermo, CThermo##Mixture##ThermoPhys);         \
    flameletAddThermoPhysicsThermo(BaseThermo, CThermo##Mixture##ThermoPhys)


#define flameletMakeThermo(BaseThermo,CThermo,Mixture,Transport,Type,Thermo,EqnOfState,Specie) \
                                                                               \
    flameletTypedefThermoPhysics(Transport,Type,Thermo,EqnOfState,Specie);             \
                                                                               \
    flameletMakeThermoPhysicsThermo                                                    \
    (                                                                          \
        BaseThermo,                                                            \
        CThermo,                                                               \
        Mixture,                                                               \
        Transport##Type##Thermo##EqnOfState##Specie                            \
    )


#define flameletMakeThermos(BaseThermo,CThermo,Mixture,Transport,Type,Thermo,EqnOfState,Specie) \
                                                                               \
    flameletTypedefThermoPhysics(Transport,Type,Thermo,EqnOfState,Specie);             \
                                                                               \
    flameletMakeThermoPhysicsThermos                                                   \
    (                                                                          \
        BaseThermo,                                                            \
        CThermo,                                                               \
        Mixture,                                                               \
        Transport##Type##Thermo##EqnOfState##Specie                            \
    )

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
