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
    Foam::flameletGradientUnburntEnthalpyFvPatchScalarField

Description
    gradient boundary condition for unburnt

SourceFiles
    flameletGradientUnburntEnthalpyFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef flameletGradientUnburntEnthalpyFvPatchScalarField_H
#define flameletGradientUnburntEnthalpyFvPatchScalarField_H

#include "fixedGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
         Class flameletGradientUnburntEnthalpyFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class flameletGradientUnburntEnthalpyFvPatchScalarField
:
    public fixedGradientFvPatchScalarField
{

public:

    //- Runtime type information
    TypeName("flameletGradientUnburntEnthalpy");


    // Constructors

        //- Construct from patch and internal field
        flameletGradientUnburntEnthalpyFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        flameletGradientUnburntEnthalpyFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given flameletGradientUnburntEnthalpyFvPatchScalarField
        // onto a new patch
        flameletGradientUnburntEnthalpyFvPatchScalarField
        (
            const flameletGradientUnburntEnthalpyFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        flameletGradientUnburntEnthalpyFvPatchScalarField
        (
            const flameletGradientUnburntEnthalpyFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new flameletGradientUnburntEnthalpyFvPatchScalarField(*this)
            );
        }

        //- Copy constructor setting internal field reference
        flameletGradientUnburntEnthalpyFvPatchScalarField
        (
            const flameletGradientUnburntEnthalpyFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new flameletGradientUnburntEnthalpyFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
