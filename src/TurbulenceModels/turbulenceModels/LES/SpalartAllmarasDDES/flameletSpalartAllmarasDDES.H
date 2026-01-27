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

Class
    Foam::LESModels::flameletSpalartAllmarasDDES

Description
    flameletSpalartAllmaras DDES turbulence model for incompressible and compressible
    flows

    Reference:
    \verbatim
        Spalart, P. R., Deck, S., Shur, M. L., Squires, K. D., Strelets, M. K.,
        & Travin, A. (2006).
        A new version of detached-eddy simulation, resistant to ambiguous grid
        densities.
        Theoretical and computational fluid dynamics, 20(3), 181-195.
    \endverbatim

SourceFiles
    flameletSpalartAllmarasDDES.C

\*---------------------------------------------------------------------------*/

#ifndef flameletSpalartAllmarasDDES_H
#define flameletSpalartAllmarasDDES_H

#include "flameletSpalartAllmarasDES.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                     Class flameletSpalartAllmarasDDES Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class flameletSpalartAllmarasDDES
:
    public flameletSpalartAllmarasDES<BasicTurbulenceModel>
{
    // Private Member Functions

        tmp<volScalarField> fd(const volScalarField& magGradU) const;

        tmp<volScalarField> rd(const volScalarField& magGradU) const;


protected:

    // Protected Member Functions

        //- Length scale
        virtual tmp<volScalarField> dTilda
        (
            const volScalarField& chi,
            const volScalarField& fv1,
            const volTensorField& gradU
        ) const;


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("flameletSpalartAllmarasDDES");


    // Constructors

        //- Construct from components
        flameletSpalartAllmarasDDES
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

        //- Disallow default bitwise copy construction
        flameletSpalartAllmarasDDES(const flameletSpalartAllmarasDDES&) = delete;


    //- Destructor
    virtual ~flameletSpalartAllmarasDDES()
    {}


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const flameletSpalartAllmarasDDES&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "flameletSpalartAllmarasDDES.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
