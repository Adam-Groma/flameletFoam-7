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
    Foam::fv::flameletInterRegionHeatTransferModel

Description
    Base class for inter region heat exchange. The derived classes must
    provide the heat transfer coeffisine (htc) which is used as follows
    in the energy equation.

     \f[
        -htc*Tmapped + Sp(htc, T)
     \f]

\*---------------------------------------------------------------------------*/

#ifndef flameletInterRegionHeatTransferModel_H
#define flameletInterRegionHeatTransferModel_H

#include "flameletInterRegionOption.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                Class flameletInterRegionHeatTransferModel Declaration
\*---------------------------------------------------------------------------*/

class flameletInterRegionHeatTransferModel
:
    public flameletInterRegionOption
{
protected:

    // Protected data

        //- Name of the model in the neighbour mesh
        word nbrModelName_;

        //- Pointer to neighbour flameletInterRegionHeatTransferModel
        flameletInterRegionHeatTransferModel* nbrModel_;

        //- First iteration
        bool firstIter_;

        //- Time index - used for updating htc
        label timeIndex_;

        //- Heat transfer coefficient [W/m^2/k] times area/volume [1/m]
        volScalarField htc_;

        //- Flag to activate semi-implicit coupling
        bool semiImplicit_;

        //- Name of temperature field; default = "T"
        word TName_;

        //- Name of neighbour temperature field; default = "T"
        word TNbrName_;


    // Protected member functions

        //- Set the neighbour flameletInterRegionHeatTransferModel
        void setNbrModel();

        //- Inherit correct from flameletInterRegionOption
        using flameletInterRegionOption::correct;

        //- Correct to calculate the inter-region heat transfer coefficient
        void correct();

        //- Interpolate field with nbrModel specified
        template<class Type>
        tmp<Field<Type>> interpolate
        (
            const flameletInterRegionHeatTransferModel& nbrModel,
            const Field<Type>& field
        ) const;

        //- Interpolate field without nbrModel specified
        template<class Type>
        tmp<Field<Type>> interpolate
        (
            const Field<Type>& field
        ) const;

        //- Interpolate field with nbrModel specified
        template<class Type>
        void interpolate
        (
            const flameletInterRegionHeatTransferModel& nbrModel,
            const Field<Type>& field,
            Field<Type>& result
        ) const;

        //- Interpolate field without nbrModel specified
        template<class Type>
        void interpolate
        (
            const Field<Type>& field,
            Field<Type>& result
        ) const;


public:

    //- Runtime type information
    TypeName("flameletInterRegionHeatTransferModel");


    // Constructors

        //- Construct from dictionary
        flameletInterRegionHeatTransferModel
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~flameletInterRegionHeatTransferModel();


    // Member Functions

        // Access

            //- Return const access to the neighbour region name
            inline const word& nbrRegionName() const;

            //- Return const access to the mapToMap pointer
            inline const meshToMesh& meshInterp() const;

            //- Return the heat transfer coefficient
            inline const volScalarField& htc() const;

            //- Return const access to the neighbour model
            inline const flameletInterRegionHeatTransferModel& nbrModel() const;

            //- Return access to the neighbour model
            inline flameletInterRegionHeatTransferModel& nbrModel();

            //- Source term to energy equation
            virtual void addSup
            (
                fvMatrix<scalar>& eqn,
                const label fieldi
            );

            //- Source term to compressible energy equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<scalar>& eqn,
                const label fieldi
            );

            //- Calculate heat transfer coefficient
            virtual void calculateHtc() = 0;


        // IO

            //- Read dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "flameletInterRegionHeatTransferModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "flameletInterRegionHeatTransferModelTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
