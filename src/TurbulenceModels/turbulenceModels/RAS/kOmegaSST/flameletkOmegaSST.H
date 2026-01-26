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

Class
    Foam::RASModels::flameletkOmegaSST

Description
    Specialisation for RAS of the generic flameletkOmegaSSTBase base class.
    For more information, see Description of flameletkOmegaSSTBase.H

See also
    Foam::flameletkOmegaSST

SourceFiles
    flameletkOmegaSST.C

\*---------------------------------------------------------------------------*/

#ifndef flameletkOmegaSST_H
#define flameletkOmegaSST_H

#include "flameletkOmegaSSTBase.H"
#include "flameletRASModel.H"
#include "flameletEddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                          Class flameletkOmegaSST Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class flameletkOmegaSST
:
    public Foam::flameletkOmegaSST
    <
        flameletEddyViscosity<flameletRASModel<BasicTurbulenceModel>>,
        BasicTurbulenceModel
    >
{

public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("flameletkOmegaSST");


    // Constructors

        //- Construct from components
        flameletkOmegaSST
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            transportModel& transport,
            const word& propertiesName = flameletTurbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~flameletkOmegaSST()
    {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
    #include "flameletkOmegaSST.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //
