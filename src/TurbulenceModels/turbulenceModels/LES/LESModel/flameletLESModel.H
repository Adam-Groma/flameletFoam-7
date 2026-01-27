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

Namespace
    Foam::LESModels

Description
    Namespace for LES SGS models.

Class
    Foam::flameletLESModel

Description
    Templated abstract base class for LES SGS models

SourceFiles
    flameletLESModel.C

\*---------------------------------------------------------------------------*/

#ifndef flameletLESModel_H
#define flameletLESModel_H

#include "FlameletTurbulenceModel.H"
#include "flameletLESdelta.H"
#include "flameletRhoReactionThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class flameletLESModel Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class flameletLESModel
:
    public BasicTurbulenceModel
{

protected:

    // Protected data

        //- LES coefficients dictionary
        dictionary LESDict_;

        //- Turbulence on/off flag
        Switch turbulence_;

        //- Flag to print the model coeffs at run-time
        Switch printCoeffs_;

        //- Flag to print transport of Z variance at run-time
        Switch transportVarZ_;

        //- Model coefficients dictionary
        dictionary coeffDict_;
        
        //- 
        dimensionedScalar CvarZ_;

        //- 
        dimensionedScalar Cchi_;

        //- 
        dimensionedScalar Sc_;

        //- 
        dimensionedScalar Sct_;

        //- 
        const Time& runTime_;

        //- Lower limit of k
        dimensionedScalar kMin_;

        //- Lower limit of epsilon
        dimensionedScalar epsilonMin_;

        //- Lower limit for omega
        dimensionedScalar omegaMin_;

        //- Run-time selectable delta model
        autoPtr<Foam::flameletLESdelta> delta_;

        //- Mixture fraction variance
        volScalarField& varZ_;

        //- Scalar dissipation rate
        volScalarField& Chi_;


    // Protected Member Functions

        //- Print model coefficients
        virtual void printCoeffs(const word& type);


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("LES");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            flameletLESModel,
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
        flameletLESModel
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
        flameletLESModel(const flameletLESModel&) = delete;


    // Selectors

        //- Return a reference to the selected LES model
        static autoPtr<flameletLESModel> New
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
    virtual ~flameletLESModel()
    {}


    // Member Functions

        //- Read model coefficients if they have changed
        virtual bool read();


        // Access

            //- Const access to the coefficients dictionary
            virtual const dictionary& coeffDict() const
            {
                return coeffDict_;
            }

            //- Return the lower allowable limit for k (default: small)
            const dimensionedScalar& kMin() const
            {
                return kMin_;
            }

            //- Allow kMin to be changed
            dimensionedScalar& kMin()
            {
                return kMin_;
            }

            //- Access function to filter width
            virtual const volScalarField& delta() const
            {
                return delta_();
            }

            //- Return mixture Fraction variance
            virtual tmp<volScalarField> varZ() const
            {
                return varZ_;
            }

            //- Return scalar dissipation rate
            virtual tmp<volScalarField> Chi() const
            {
                return Chi_;
            }

            //- Return SGS mixture fraction diffusivity
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

        //- Calculate mixture fraction variance
        virtual void correctVarZ();

        //- Calculate scalar dissipation rate
        virtual void correctChi();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const flameletLESModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "flameletLESModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
