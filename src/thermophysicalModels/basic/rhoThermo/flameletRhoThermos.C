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

#include "specie.H"
#include "perfectGas.H"
#include "incompressiblePerfectGas.H"
#include "Boussinesq.H"
#include "rhoConst.H"
#include "perfectFluid.H"

#include "hConstThermo.H"
#include "eConstThermo.H"
#include "janafThermo.H"
#include "sensibleEnthalpy.H"
#include "sensibleInternalEnergy.H"
#include "thermo.H"

#include "flameletConstTransport.H"
#include "sutherlandTransport.H"
#include "WLFTransport.H"

#include "icoPolynomial.H"
#include "hPolynomialThermo.H"
#include "polynomialTransport.H"

#include "flameletHeRhoThermo.H"
#include "pureMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

flameletMakeThermos
(
    flameletRhoThermo,
    flameletHeRhoThermo,
    pureMixture,
    flameletConstTransport,
    sensibleEnthalpy,
    hConstThermo,
    perfectGas,
    specie
);

flameletMakeThermos
(
    flameletRhoThermo,
    flameletHeRhoThermo,
    pureMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    hConstThermo,
    perfectGas,
    specie
);

flameletMakeThermos
(
    flameletRhoThermo,
    flameletHeRhoThermo,
    pureMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    janafThermo,
    perfectGas,
    specie
);

flameletMakeThermos
(
    flameletRhoThermo,
    flameletHeRhoThermo,
    pureMixture,
    flameletConstTransport,
    sensibleEnthalpy,
    hConstThermo,
    rhoConst,
    specie
);

flameletMakeThermos
(
    flameletRhoThermo,
    flameletHeRhoThermo,
    pureMixture,
    flameletConstTransport,
    sensibleEnthalpy,
    hConstThermo,
    perfectFluid,
    specie
);

flameletMakeThermos
(
    flameletRhoThermo,
    flameletHeRhoThermo,
    pureMixture,
    polynomialTransport,
    sensibleEnthalpy,
    hPolynomialThermo,
    icoPolynomial,
    specie
);

flameletMakeThermos
(
    flameletRhoThermo,
    flameletHeRhoThermo,
    pureMixture,
    flameletConstTransport,
    sensibleEnthalpy,
    hConstThermo,
    incompressiblePerfectGas,
    specie
);

flameletMakeThermos
(
    flameletRhoThermo,
    flameletHeRhoThermo,
    pureMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    hConstThermo,
    incompressiblePerfectGas,
    specie
);

flameletMakeThermos
(
    flameletRhoThermo,
    flameletHeRhoThermo,
    pureMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    janafThermo,
    incompressiblePerfectGas,
    specie
);

flameletMakeThermos
(
    flameletRhoThermo,
    flameletHeRhoThermo,
    pureMixture,
    flameletConstTransport,
    sensibleEnthalpy,
    hConstThermo,
    Boussinesq,
    specie
);

flameletMakeThermos
(
    flameletRhoThermo,
    flameletHeRhoThermo,
    pureMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    hConstThermo,
    Boussinesq,
    specie
);

flameletMakeThermos
(
    flameletRhoThermo,
    flameletHeRhoThermo,
    pureMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    janafThermo,
    Boussinesq,
    specie
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

flameletMakeThermos
(
    flameletRhoThermo,
    flameletHeRhoThermo,
    pureMixture,
    flameletConstTransport,
    sensibleInternalEnergy,
    hConstThermo,
    perfectGas,
    specie
);

flameletMakeThermos
(
    flameletRhoThermo,
    flameletHeRhoThermo,
    pureMixture,
    sutherlandTransport,
    sensibleInternalEnergy,
    hConstThermo,
    perfectGas,
    specie
);

flameletMakeThermos
(
    flameletRhoThermo,
    flameletHeRhoThermo,
    pureMixture,
    sutherlandTransport,
    sensibleInternalEnergy,
    janafThermo,
    perfectGas,
    specie
);

flameletMakeThermos
(
    flameletRhoThermo,
    flameletHeRhoThermo,
    pureMixture,
    flameletConstTransport,
    sensibleInternalEnergy,
    hConstThermo,
    rhoConst,
    specie
);

flameletMakeThermos
(
    flameletRhoThermo,
    flameletHeRhoThermo,
    pureMixture,
    flameletConstTransport,
    sensibleInternalEnergy,
    eConstThermo,
    rhoConst,
    specie
);

flameletMakeThermos
(
    flameletRhoThermo,
    flameletHeRhoThermo,
    pureMixture,
    flameletConstTransport,
    sensibleInternalEnergy,
    hConstThermo,
    perfectFluid,
    specie
);

flameletMakeThermos
(
    flameletRhoThermo,
    flameletHeRhoThermo,
    pureMixture,
    flameletConstTransport,
    sensibleInternalEnergy,
    eConstThermo,
    perfectFluid,
    specie
);

flameletMakeThermos
(
    flameletRhoThermo,
    flameletHeRhoThermo,
    pureMixture,
    polynomialTransport,
    sensibleInternalEnergy,
    hPolynomialThermo,
    icoPolynomial,
    specie
);

flameletMakeThermos
(
    flameletRhoThermo,
    flameletHeRhoThermo,
    pureMixture,
    flameletConstTransport,
    sensibleInternalEnergy,
    hConstThermo,
    incompressiblePerfectGas,
    specie
);

flameletMakeThermos
(
    flameletRhoThermo,
    flameletHeRhoThermo,
    pureMixture,
    sutherlandTransport,
    sensibleInternalEnergy,
    hConstThermo,
    incompressiblePerfectGas,
    specie
);

flameletMakeThermos
(
    flameletRhoThermo,
    flameletHeRhoThermo,
    pureMixture,
    sutherlandTransport,
    sensibleInternalEnergy,
    janafThermo,
    incompressiblePerfectGas,
    specie
);

flameletMakeThermos
(
    flameletRhoThermo,
    flameletHeRhoThermo,
    pureMixture,
    flameletConstTransport,
    sensibleInternalEnergy,
    eConstThermo,
    Boussinesq,
    specie
);

flameletMakeThermos
(
    flameletRhoThermo,
    flameletHeRhoThermo,
    pureMixture,
    flameletConstTransport,
    sensibleInternalEnergy,
    hConstThermo,
    Boussinesq,
    specie
);

flameletMakeThermos
(
    flameletRhoThermo,
    flameletHeRhoThermo,
    pureMixture,
    sutherlandTransport,
    sensibleInternalEnergy,
    hConstThermo,
    Boussinesq,
    specie
);

flameletMakeThermos
(
    flameletRhoThermo,
    flameletHeRhoThermo,
    pureMixture,
    sutherlandTransport,
    sensibleInternalEnergy,
    janafThermo,
    Boussinesq,
    specie
);

flameletMakeThermos
(
    flameletRhoThermo,
    flameletHeRhoThermo,
    pureMixture,
    WLFTransport,
    sensibleInternalEnergy,
    eConstThermo,
    rhoConst,
    specie
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
