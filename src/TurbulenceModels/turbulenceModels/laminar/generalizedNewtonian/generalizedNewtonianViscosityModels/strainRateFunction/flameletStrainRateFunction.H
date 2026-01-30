/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2019 OpenFOAM Foundation
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
    Foam::laminarModels::generalizedNewtonianViscosityModels::flameletStrainRateFunction

Description
    Run-time selected strain-rate function generalized Newtonian viscosity model

    Example linear function of strain-rate:
    \verbatim
        generalizedNewtonianModel flameletStrainRateFunction;

        function polynomial ((0 0.1) (1 1.3));
    \endverbatim

See also
    Foam::flameletGeneralizedNewtonianViscosityModel
    Foam::Function1

SourceFiles
    flameletStrainRateFunction.C

\*---------------------------------------------------------------------------*/

#ifndef flameletStrainRateFunction_H
#define flameletStrainRateFunction_H

#include "flameletGeneralizedNewtonianViscosityModel.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarModels
{
namespace generalizedNewtonianViscosityModels
{

/*---------------------------------------------------------------------------*\
                           Class flameletStrainRateFunction Declaration
\*---------------------------------------------------------------------------*/

class flameletStrainRateFunction
:
    public flameletGeneralizedNewtonianViscosityModel
{
    // Private Data

        //- Strain-rate function
        autoPtr<Function1<scalar>> strainRateFunction_;


public:

    //- Runtime type information
    TypeName("flameletStrainRateFunction");


    // Constructors

        //- Construct from components
        flameletStrainRateFunction(const dictionary& viscosityProperties);


    //- Destructor
    virtual ~flameletStrainRateFunction()
    {}


    // Member Functions

        //- Read transportProperties dictionary
        virtual bool read(const dictionary& viscosityProperties);

        //- Return the laminar viscosity
        virtual tmp<volScalarField> nu
        (
            const volScalarField& nu0,
            const volScalarField& strainRate
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace generalizedNewtonianViscosityModels
} // End namespace laminarModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
