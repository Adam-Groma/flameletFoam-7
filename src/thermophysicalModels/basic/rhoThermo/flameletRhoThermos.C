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

#include "flameletRhoThermo.H"
#include "flameletMakeThermo.H"

#include "flameletSpecie.H"
#include "flameletPerfectGas.H"
#include "flameletIncompressiblePerfectGas.H"
#include "flameletBoussinesq.H"
#include "flameletRhoConst.H"
#include "flameletPerfectFluid.H"

#include "flameletHConstThermo.H"
#include "flameletEConstThermo.H"
#include "flameletJanafThermo.H"
#include "flameletSensibleEnthalpy.H"
#include "flameletSensibleInternalEnergy.H"
#include "flameletThermo.H"

#include "flameletConstTransport.H"
#include "flameletSutherlandTransport.H"
#include "flameletWLFTransport.H"

#include "flameletIcoPolynomial.H"
#include "flameletHPolynomialThermo.H"
#include "flameletPolynomialTransport.H"

#include "flameletHeRhoThermo.H"
#include "flameletPureMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

flameletMakeThermos
(
    flameletRhoThermo,
    flameletHeRhoThermo,
    flameletPureMixture,
    flameletConstTransport,
    flameletSensibleEnthalpy,
    flameletHConstThermo,
    flameletPerfectGas,
    flameletSpecie
);

flameletMakeThermos
(
    flameletRhoThermo,
    flameletHeRhoThermo,
    flameletPureMixture,
    flameletSutherlandTransport,
    flameletSensibleEnthalpy,
    flameletHConstThermo,
    flameletPerfectGas,
    flameletSpecie
);

flameletMakeThermos
(
    flameletRhoThermo,
    flameletHeRhoThermo,
    flameletPureMixture,
    flameletSutherlandTransport,
    flameletSensibleEnthalpy,
    flameletJanafThermo,
    flameletPerfectGas,
    flameletSpecie
);

flameletMakeThermos
(
    flameletRhoThermo,
    flameletHeRhoThermo,
    flameletPureMixture,
    flameletConstTransport,
    flameletSensibleEnthalpy,
    flameletHConstThermo,
    flameletRhoConst,
    flameletSpecie
);

flameletMakeThermos
(
    flameletRhoThermo,
    flameletHeRhoThermo,
    flameletPureMixture,
    flameletConstTransport,
    flameletSensibleEnthalpy,
    flameletHConstThermo,
    flameletPerfectFluid,
    flameletSpecie
);

flameletMakeThermos
(
    flameletRhoThermo,
    flameletHeRhoThermo,
    flameletPureMixture,
    flameletPolynomialTransport,
    flameletSensibleEnthalpy,
    flameletHPolynomialThermo,
    flameletIcoPolynomial,
    flameletSpecie
);

flameletMakeThermos
(
    flameletRhoThermo,
    flameletHeRhoThermo,
    flameletPureMixture,
    flameletConstTransport,
    flameletSensibleEnthalpy,
    flameletHConstThermo,
    flameletIncompressiblePerfectGas,
    flameletSpecie
);

flameletMakeThermos
(
    flameletRhoThermo,
    flameletHeRhoThermo,
    flameletPureMixture,
    flameletSutherlandTransport,
    flameletSensibleEnthalpy,
    flameletHConstThermo,
    flameletIncompressiblePerfectGas,
    flameletSpecie
);

flameletMakeThermos
(
    flameletRhoThermo,
    flameletHeRhoThermo,
    flameletPureMixture,
    flameletSutherlandTransport,
    flameletSensibleEnthalpy,
    flameletJanafThermo,
    flameletIncompressiblePerfectGas,
    flameletSpecie
);

flameletMakeThermos
(
    flameletRhoThermo,
    flameletHeRhoThermo,
    flameletPureMixture,
    flameletConstTransport,
    flameletSensibleEnthalpy,
    flameletHConstThermo,
    flameletBoussinesq,
    flameletSpecie
);

flameletMakeThermos
(
    flameletRhoThermo,
    flameletHeRhoThermo,
    flameletPureMixture,
    flameletSutherlandTransport,
    flameletSensibleEnthalpy,
    flameletHConstThermo,
    flameletBoussinesq,
    flameletSpecie
);

flameletMakeThermos
(
    flameletRhoThermo,
    flameletHeRhoThermo,
    flameletPureMixture,
    flameletSutherlandTransport,
    flameletSensibleEnthalpy,
    flameletJanafThermo,
    flameletBoussinesq,
    flameletSpecie
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

flameletMakeThermos
(
    flameletRhoThermo,
    flameletHeRhoThermo,
    flameletPureMixture,
    flameletConstTransport,
    flameletSensibleInternalEnergy,
    flameletHConstThermo,
    flameletPerfectGas,
    flameletSpecie
);

flameletMakeThermos
(
    flameletRhoThermo,
    flameletHeRhoThermo,
    flameletPureMixture,
    flameletSutherlandTransport,
    flameletSensibleInternalEnergy,
    flameletHConstThermo,
    flameletPerfectGas,
    flameletSpecie
);

flameletMakeThermos
(
    flameletRhoThermo,
    flameletHeRhoThermo,
    flameletPureMixture,
    flameletSutherlandTransport,
    flameletSensibleInternalEnergy,
    flameletJanafThermo,
    flameletPerfectGas,
    flameletSpecie
);

flameletMakeThermos
(
    flameletRhoThermo,
    flameletHeRhoThermo,
    flameletPureMixture,
    flameletConstTransport,
    flameletSensibleInternalEnergy,
    flameletHConstThermo,
    flameletRhoConst,
    flameletSpecie
);

flameletMakeThermos
(
    flameletRhoThermo,
    flameletHeRhoThermo,
    flameletPureMixture,
    flameletConstTransport,
    flameletSensibleInternalEnergy,
    flameletEConstThermo,
    flameletRhoConst,
    flameletSpecie
);

flameletMakeThermos
(
    flameletRhoThermo,
    flameletHeRhoThermo,
    flameletPureMixture,
    flameletConstTransport,
    flameletSensibleInternalEnergy,
    flameletHConstThermo,
    flameletPerfectFluid,
    flameletSpecie
);

flameletMakeThermos
(
    flameletRhoThermo,
    flameletHeRhoThermo,
    flameletPureMixture,
    flameletConstTransport,
    flameletSensibleInternalEnergy,
    flameletEConstThermo,
    flameletPerfectFluid,
    flameletSpecie
);

flameletMakeThermos
(
    flameletRhoThermo,
    flameletHeRhoThermo,
    flameletPureMixture,
    flameletPolynomialTransport,
    flameletSensibleInternalEnergy,
    flameletHPolynomialThermo,
    flameletIcoPolynomial,
    flameletSpecie
);

flameletMakeThermos
(
    flameletRhoThermo,
    flameletHeRhoThermo,
    flameletPureMixture,
    flameletConstTransport,
    flameletSensibleInternalEnergy,
    flameletHConstThermo,
    flameletIncompressiblePerfectGas,
    flameletSpecie
);

flameletMakeThermos
(
    flameletRhoThermo,
    flameletHeRhoThermo,
    flameletPureMixture,
    flameletSutherlandTransport,
    flameletSensibleInternalEnergy,
    flameletHConstThermo,
    flameletIncompressiblePerfectGas,
    flameletSpecie
);

flameletMakeThermos
(
    flameletRhoThermo,
    flameletHeRhoThermo,
    flameletPureMixture,
    flameletSutherlandTransport,
    flameletSensibleInternalEnergy,
    flameletJanafThermo,
    flameletIncompressiblePerfectGas,
    flameletSpecie
);

flameletMakeThermos
(
    flameletRhoThermo,
    flameletHeRhoThermo,
    flameletPureMixture,
    flameletConstTransport,
    flameletSensibleInternalEnergy,
    flameletEConstThermo,
    flameletBoussinesq,
    flameletSpecie
);

flameletMakeThermos
(
    flameletRhoThermo,
    flameletHeRhoThermo,
    flameletPureMixture,
    flameletConstTransport,
    flameletSensibleInternalEnergy,
    flameletHConstThermo,
    flameletBoussinesq,
    flameletSpecie
);

flameletMakeThermos
(
    flameletRhoThermo,
    flameletHeRhoThermo,
    flameletPureMixture,
    flameletSutherlandTransport,
    flameletSensibleInternalEnergy,
    flameletHConstThermo,
    flameletBoussinesq,
    flameletSpecie
);

flameletMakeThermos
(
    flameletRhoThermo,
    flameletHeRhoThermo,
    flameletPureMixture,
    flameletSutherlandTransport,
    flameletSensibleInternalEnergy,
    flameletJanafThermo,
    flameletBoussinesq,
    flameletSpecie
);

flameletMakeThermos
(
    flameletRhoThermo,
    flameletHeRhoThermo,
    flameletPureMixture,
    flameletWLFTransport,
    flameletSensibleInternalEnergy,
    flameletEConstThermo,
    flameletRhoConst,
    flameletSpecie
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
