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
    Foam::flameletWallHeatTransferFvPatchScalarField

Description
    This boundary condition provides an enthalpy condition for wall heat
    transfer

Usage
    \table
        Property     | Description             | Required    | Default value
        Tinf         | wall temperature        | yes         |
        alphaWall    | thermal diffusivity     | yes         |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            flameletWallHeatTransfer;
        Tinf            uniform 500;
        alphaWall       uniform 1;
    }
    \endverbatim

SourceFiles
    flameletWallHeatTransferFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef flameletWallHeatTransferFvPatchScalarField_H
#define flameletWallHeatTransferFvPatchScalarField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class flameletWallHeatTransferFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class flameletWallHeatTransferFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private Data

        //- Temperature at the wall
        scalarField Tinf_;

        //- Thermal diffusivity at the wall
        scalarField alphaWall_;


public:

    //- Runtime type information
    TypeName("flameletWallHeatTransfer");


    // Constructors

        //- Construct from patch and internal field
        flameletWallHeatTransferFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        flameletWallHeatTransferFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given flameletWallHeatTransferFvPatchScalarField
        // onto a new patch
        flameletWallHeatTransferFvPatchScalarField
        (
            const flameletWallHeatTransferFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        flameletWallHeatTransferFvPatchScalarField
        (
            const flameletWallHeatTransferFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new flameletWallHeatTransferFvPatchScalarField(*this)
            );
        }

        //- Copy constructor setting internal field reference
        flameletWallHeatTransferFvPatchScalarField
        (
            const flameletWallHeatTransferFvPatchScalarField&,
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
                new flameletWallHeatTransferFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        // Access

            //- Return Tinf
            const scalarField& Tinf() const
            {
                return Tinf_;
            }

            //- Return reference to Tinf to allow adjustment
            scalarField& Tinf()
            {
                return Tinf_;
            }

            //- Return alphaWall
            const scalarField& alphaWall() const
            {
                return alphaWall_;
            }

            //- Return reference to alphaWall to allow adjustment
            scalarField& alphaWall()
            {
                return alphaWall_;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            //  Used to update fields following mesh topology change
            virtual void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            //  Used to reconstruct fields
            virtual void rmap(const fvPatchScalarField&, const labelList&);


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
