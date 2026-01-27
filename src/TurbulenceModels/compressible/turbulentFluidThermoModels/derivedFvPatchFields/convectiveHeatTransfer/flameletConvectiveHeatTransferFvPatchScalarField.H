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
    Foam::compressible::flameletConvectiveHeatTransferFvPatchScalarField

Description
    This boundary condition provides a convective heat transfer coefficient
    condition

    if Re > 500000
    \f[
        htc_p = \frac{0.664 Re^{0.5} Pr^{0.333} \kappa_p}{L}
    \f]
    else
    \f[
        htc_p = \frac{0.037 Re^{0.8} Pr^{0.333} \kappa_p}{L}
    \f]

    where

    \vartable
        htc_p   | patch convective heat transfer coefficient
        Re      | Reynolds number
        Pr      | Prandtl number
        \kappa_p | thermal conductivity
        L       | length scale
    \endvartable

Usage
    \table
        Property     | Description             | Required    | Default value
        L            | Length scale [m]        | yes |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            flameletConvectiveHeatTransfer;
        L               0.1;
    }
    \endverbatim

See also
    Foam::fixedValueFvPatchField

SourceFiles
    flameletConvectiveHeatTransferFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef flameletConvectiveHeatTransferFvPatchScalarField_H
#define flameletConvectiveHeatTransferFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

/*---------------------------------------------------------------------------*\
           Class flameletConvectiveHeatTransferFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class flameletConvectiveHeatTransferFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
protected:

    // Protected data

        //- L Length scale [m]
        const scalar L_;


public:

    //- Runtime type information
    TypeName("flameletConvectiveHeatTransfer");


    // Constructors

        //- Construct from patch and internal field
        flameletConvectiveHeatTransferFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        flameletConvectiveHeatTransferFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  flameletConvectiveHeatTransferFvPatchScalarField
        //  onto a new patch
        flameletConvectiveHeatTransferFvPatchScalarField
        (
            const flameletConvectiveHeatTransferFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        flameletConvectiveHeatTransferFvPatchScalarField
        (
            const flameletConvectiveHeatTransferFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new flameletConvectiveHeatTransferFvPatchScalarField(*this)
            );
        }

        //- Copy constructor setting internal field reference
        flameletConvectiveHeatTransferFvPatchScalarField
        (
            const flameletConvectiveHeatTransferFvPatchScalarField&,
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
                new flameletConvectiveHeatTransferFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
