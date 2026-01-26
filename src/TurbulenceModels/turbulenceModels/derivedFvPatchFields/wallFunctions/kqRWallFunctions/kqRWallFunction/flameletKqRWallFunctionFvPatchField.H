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
    Foam::flameletKqRWallFunctionFvPatchField

Description
    This boundary condition provides a suitable condition for turbulence
    \c k, \c q, and \c R fields for the case of high Reynolds number flow using
    wall functions.

    It is a simple wrapper around the zero-gradient condition.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            flameletKqRWallFunction;
    }
    \endverbatim

See also
    Foam::zeroGradientFvPatchField

SourceFiles
    flameletKqRWallFunctionFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef flameletKqRWallFunctionFvPatchField_H
#define flameletKqRWallFunctionFvPatchField_H

#include "zeroGradientFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class flameletKqRWallFunctionFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class flameletKqRWallFunctionFvPatchField
:
    public zeroGradientFvPatchField<Type>
{

public:

    //- Runtime type information
    TypeName("flameletKqRWallFunction");


    // Constructors

        //- Construct from patch and internal field
        flameletKqRWallFunctionFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        flameletKqRWallFunctionFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        // flameletKqRWallFunctionFvPatchField
        //  onto a new patch
        flameletKqRWallFunctionFvPatchField
        (
            const flameletKqRWallFunctionFvPatchField&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        flameletKqRWallFunctionFvPatchField
        (
            const flameletKqRWallFunctionFvPatchField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new flameletKqRWallFunctionFvPatchField(*this)
            );
        }

        //- Copy constructor setting internal field reference
        flameletKqRWallFunctionFvPatchField
        (
            const flameletKqRWallFunctionFvPatchField&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>
            (
                new flameletKqRWallFunctionFvPatchField(*this, iF)
            );
        }


    // Member Functions

        //- Evaluate the patchField
        virtual void evaluate
        (
            const Pstream::commsTypes commsType = Pstream::commsTypes::blocking
        );

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "flameletKqRWallFunctionFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
