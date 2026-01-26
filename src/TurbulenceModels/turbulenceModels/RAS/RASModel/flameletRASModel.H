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
    Foam::flameletRASModel

Description
    Templated abstract base class for RAS turbulence models

SourceFiles
    flameletRASModel.C

\*---------------------------------------------------------------------------*/

#ifndef flameletRASModel_H
#define flameletRASModel_H

#include "FlameletTurbulenceModel.H"
#include "fvm.H"
#include "bound.H"
#include "flameletRhoReactionThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class flameletRASModel Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class flameletRASModel
:
    public BasicTurbulenceModel
{

protected:

    // Protected data

        //- RAS coefficients dictionary
        dictionary RASDict_;

        //- Turbulence on/off flag
        Switch turbulence_;

        //- Flag to print the model coeffs at run-time
        Switch printCoeffs_;

        //- Model coefficients dictionary
        dictionary coeffDict_;
	
        //- scalar dissipation rate model constant
  	    dimensionedScalar Cchi_;

	      // laminar Schmidt number
	      dimensionedScalar Sc_;

        // turbulent Schmidt number
	      dimensionedScalar Sct_;

        //- Lower limit of k
        dimensionedScalar kMin_;

        //- Lower limit of epsilon
        dimensionedScalar epsilonMin_;

        //- Lower limit for omega
        dimensionedScalar omegaMin_;

        //- Mixture fraction variance
        volScalarField& varZ_;

        //- Scalar Dissipation Rate
        volScalarField& Chi_;


    // Protected Member Functions

        //- Print model coefficients
        virtual void printCoeffs(const word& type);


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("flameletRAS");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            flameletRASModel,
            dictionary,
            (
                const alphaField& alpha,
                const rhoField& rho,
                const volVectorField& U,
                const surfaceScalarField& alphaRhoPhi,
                const surfaceScalarField& phi,
                transportModel& transport,
                const word& propertiesName
            ),
            (alpha, rho, U, alphaRhoPhi, phi, transport, propertiesName)
        );


    // Constructors

        //- Construct from components
        flameletRASModel
        (
            const word& type,
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            transportModel& transport,
            const word& propertiesName
        );

        //- Disallow default bitwise copy construction
        flameletRASModel(const flameletRASModel&) = delete;


    // Selectors

        //- Return a reference to the selected RAS model
        static autoPtr<flameletRASModel> New
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            transportModel& transport,
            const word& propertiesName = flameletTurbulenceModel::propertiesName
        );


    //- Destructor
    virtual ~flameletRASModel()
    {}


    // Member Functions

        //- Read model coefficients if they have changed
        virtual bool read();


        // Access

            //- Return the lower allowable limit for k (default: small)
            const dimensionedScalar& kMin() const
            {
                return kMin_;
            }

            //- Return the lower allowable limit for epsilon (default: small)
            const dimensionedScalar& epsilonMin() const
            {
                return epsilonMin_;
            }

            //- Return the lower allowable limit for omega (default: small)
            const dimensionedScalar& omegaMin() const
            {
                return omegaMin_;
            }

            //- Allow kMin to be changed
            dimensionedScalar& kMin()
            {
                return kMin_;
            }

            //- Allow epsilonMin to be changed
            dimensionedScalar& epsilonMin()
            {
                return epsilonMin_;
            }

            //- Allow omegaMin to be changed
            dimensionedScalar& omegaMin()
            {
                return omegaMin_;
            }

            //- Const access to the coefficients dictionary
            virtual const dictionary& coeffDict() const
            {
                return coeffDict_;
            }

            //- Return mixture Fraction variance
            virtual tmp<volScalarField> varZ() const
            {
                return varZ_;
            }

            //- Return thermal diffusivity
            virtual tmp<volScalarField> Chi() const
            {
                return Chi_;
            }

            //- Return the turbulent mixture fraction diffusivity
            virtual tmp<volScalarField> DZt() const
            {
                return tmp<volScalarField>
                (
                    new volScalarField("DZt", this->mut()/Sct_)
                );
            }

            //- Return the effective mixture fraction diffusivity
            virtual tmp<volScalarField> DZEff() const
            {
                return tmp<volScalarField>
                (
                    new volScalarField("DZEff", this->mut()/Sct_ + this->mu()/Sc_)
                );
            }


        //- Return the effective viscosity
        virtual tmp<volScalarField> nuEff() const
        {
            return volScalarField::New
            (
                IOobject::groupName("nuEff", this->alphaRhoPhi_.group()),
                this->nut() + this->nu()
            );
        }

        //- Return the effective viscosity on patch
        virtual tmp<scalarField> nuEff(const label patchi) const
        {
            return this->nut(patchi) + this->nu(patchi);
        }

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();

        // Added for calculation of varZ
        virtual void correctVarZ();

        // Added for calculation of chi
        virtual void correctChi();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const flameletRASModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "flameletRASModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
