/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019 OpenFOAM Foundation
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
    Foam::laminarModels::flameletGiesekus

Description
    flameletGiesekus model for viscoelasticity.

    Reference:
    \verbatim
        Giesekus, H., 1982.
        A simple constitutive equation for polymer fluids based on the
        concept of deformation-dependent tensional mobility.
        J. Non-Newton. Fluid. 11, 69–109.
    \endverbatim

See also
    Foam::laminarModels::flameletMaxwell

SourceFiles
    flameletGiesekus.C

\*---------------------------------------------------------------------------*/

#ifndef flameletGiesekus_H
#define flameletGiesekus_H

#include "flameletMaxwell.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarModels
{

/*---------------------------------------------------------------------------*\
                           Class flameletGiesekus Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class flameletGiesekus
:
    public flameletMaxwell<BasicTurbulenceModel>
{
protected:

    // Protected data

        // Model coefficients

            dimensionedScalar alphaG_;


    // Protected Member Functions

        virtual tmp<fvSymmTensorMatrix> sigmaSource() const;


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("flameletGiesekus");


    // Constructors

        //- Construct from components
        flameletGiesekus
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
        flameletGiesekus(const flameletGiesekus&) = delete;



    //- Destructor
    virtual ~flameletGiesekus()
    {}


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const flameletGiesekus&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace laminarModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "flameletGiesekus.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
