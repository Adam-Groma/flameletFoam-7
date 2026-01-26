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

Class
    Foam::LESModels::flameletLESflameletLESModel

Description
    Eddy viscosity LES SGS model base class

SourceFiles
    flameletLESflameletLESModel.C

\*---------------------------------------------------------------------------*/

#ifndef flameletLESflameletLESModel_H
#define flameletLESflameletLESModel_H

#include "flameletLESModel.H"
#include "flameletEddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                           Class flameletLESflameletLESModel Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class flameletLESflameletLESModel
:
    public flameletEddyViscosity<flameletLESModel<BasicTurbulenceModel>>
{
protected:

    // Protected data

        dimensionedScalar Ce_;


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    // Constructors

        //- Construct from components
        flameletLESflameletLESModel
        (
            const word& type,
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName
        );

        //- Disallow default bitwise copy construction
        flameletLESflameletLESModel(const flameletLESflameletLESModel&) = delete;


    //- Destructor
    virtual ~flameletLESflameletLESModel()
    {}


    // Member Functions

        //- Read model coefficients if they have changed
        virtual bool read();

        //- Return sub-grid disipation rate
        virtual tmp<volScalarField> epsilon() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const flameletLESflameletLESModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "flameletLESflameletLESModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
