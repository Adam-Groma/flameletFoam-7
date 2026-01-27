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
    Foam::flameletTotalFlowRateAdvectiveDiffusiveFvPatchScalarField

Description
    This BC is used for species inlets. The diffusion and advection fluxes are
    considered to calculate the inlet value for the species
    The massFluxFraction sets the fraction of the flux of each particular
    species.

SourceFiles
    flameletTotalFlowRateAdvectiveDiffusiveFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef flameletTotalFlowRateAdvectiveDiffusiveFvPatchScalarField_H
#define flameletTotalFlowRateAdvectiveDiffusiveFvPatchScalarField_H

#include "mixedFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
     Class flameletTotalFlowRateAdvectiveDiffusiveFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class flameletTotalFlowRateAdvectiveDiffusiveFvPatchScalarField
:
    public mixedFvPatchField<scalar>
{
    // Private Data

        //- Name of the flux transporting the field
        word phiName_;

        //- Name of the density field used to normalise the mass flux
        //  if necessary
        word rhoName_;

        //- Mass flux fraction
        scalar massFluxFraction_;

public:

    //- Runtime type information
    TypeName("flameletTotalFlowRateAdvectiveDiffusive");


    // Constructors

        //- Construct from patch and internal field
        flameletTotalFlowRateAdvectiveDiffusiveFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        flameletTotalFlowRateAdvectiveDiffusiveFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  flameletTotalFlowRateAdvectiveDiffusiveFvPatchScalarField
        //  onto a new patch
        flameletTotalFlowRateAdvectiveDiffusiveFvPatchScalarField
        (
            const flameletTotalFlowRateAdvectiveDiffusiveFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        flameletTotalFlowRateAdvectiveDiffusiveFvPatchScalarField
        (
            const flameletTotalFlowRateAdvectiveDiffusiveFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<scalar>> clone() const
        {
            return tmp<fvPatchField<scalar>>
            (
                new
                flameletTotalFlowRateAdvectiveDiffusiveFvPatchScalarField(*this)
            );
        }

        //- Copy constructor setting internal field reference
        flameletTotalFlowRateAdvectiveDiffusiveFvPatchScalarField
        (
            const flameletTotalFlowRateAdvectiveDiffusiveFvPatchScalarField&,
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
                new
                flameletTotalFlowRateAdvectiveDiffusiveFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member Functions

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
