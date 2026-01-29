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

\*---------------------------------------------------------------------------*/

#ifndef flameletMakeChemistrySolverTypes_H
#define flameletMakeChemistrySolverTypes_H

#include "flameletChemistrySolver.H"

#include "FlameletStandardChemistryModel.H"

#include "flameletNoChemistrySolver.H"
#include "flameletEulerImplicit.H"
#include "flameletOde.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define flameletMakeChemistrySolverType(SS, Comp, Thermo)                              \
                                                                               \
    typedef SS<FlameletStandardChemistryModel<Comp, Thermo>> SS##Comp##Thermo;         \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        SS##Comp##Thermo,                                                      \
        (#SS"<" + word(FlameletStandardChemistryModel<Comp, Thermo>::typeName_()) + "<"\
        + word(Comp::typeName_()) + "," + Thermo::typeName() + ">>").c_str(),  \
        0                                                                      \
    );                                                                         \
                                                                               \
    FlameletBasicChemistryModel<Comp>::                                                \
        add##thermo##ConstructorToTable<SS##Comp##Thermo>                      \
        add##SS##Comp##Thermo##thermo##ConstructorTo##FlameletBasicChemistryModel##Comp\
##Table_;


#define flameletMakeChemistrySolverTypes(Comp, Thermo)                                 \
                                                                               \
    flameletMakeChemistrySolverType                                                    \
    (                                                                          \
        flameletNoChemistrySolver,                                                     \
        Comp,                                                                  \
        Thermo                                                                 \
    );                                                                         \
                                                                               \
    flameletMakeChemistrySolverType                                                    \
    (                                                                          \
        flameletEulerImplicit,                                                         \
        Comp,                                                                  \
        Thermo                                                                 \
    );                                                                         \
                                                                               \
    flameletMakeChemistrySolverType                                                    \
    (                                                                          \
        flameletOde,                                                                   \
        Comp,                                                                  \
        Thermo                                                                 \
    );                                                                         \


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
