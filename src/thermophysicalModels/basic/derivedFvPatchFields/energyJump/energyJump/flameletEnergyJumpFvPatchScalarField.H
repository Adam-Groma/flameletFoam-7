/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2019 OpenFOAM Foundation
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
    Foam::flameletEnergyJumpFvPatchScalarField

Description
    This boundary condition provides an energy jump condition across a pair
    of coupled patches.  It is not applied directly, but is employed on-the-fly
    when converting temperature boundary conditions into energy.

See also
    Foam::fixedJumpFvPatchField

SourceFiles
    flameletEnergyJumpFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef flameletEnergyJumpFvPatchScalarField_H
#define flameletEnergyJumpFvPatchScalarField_H

#include "fixedJumpFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class flameletEnergyJumpFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class flameletEnergyJumpFvPatchScalarField
:
    public fixedJumpFvPatchField<scalar>
{

public:

    //- Runtime type information
    TypeName("flameletEnergyJump");

    // Constructors

        //- Construct from patch and internal field
        flameletEnergyJumpFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        flameletEnergyJumpFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given flameletEnergyJumpFvPatchScalarField onto a
        //  new patch
        flameletEnergyJumpFvPatchScalarField
        (
            const flameletEnergyJumpFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        flameletEnergyJumpFvPatchScalarField
        (
            const flameletEnergyJumpFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<scalar>> clone() const
        {
            return tmp<fvPatchField<scalar>>
            (
                new flameletEnergyJumpFvPatchScalarField(*this)
            );
        }

        //- Copy constructor setting internal field reference
        flameletEnergyJumpFvPatchScalarField
        (
            const flameletEnergyJumpFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<scalar>> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<scalar>>
            (
                new flameletEnergyJumpFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        // Evaluation functions

            //- Update the coefficients
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
