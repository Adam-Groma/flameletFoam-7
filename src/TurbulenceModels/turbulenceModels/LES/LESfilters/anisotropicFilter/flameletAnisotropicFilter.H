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
    Foam::flameletAnisotropicFilter

Description
    anisotropic filter

    \verbatim
    Kernel                 as filter          as Test filter with ratio 2
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Box filter:            g = delta2/24  ->  g = delta2/6
    Spherical box filter:  g = delta2/64  ->  g = delta2/16
    Gaussian filter:       g = delta2/24  ->  g = delta2/6
    \endverbatim

SourceFiles
    flameletAnisotropicFilter.C

\*---------------------------------------------------------------------------*/

#ifndef flameletAnisotropicFilter_H
#define flameletAnisotropicFilter_H

#include "flameletLESfilter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class flameletAnisotropicFilter Declaration
\*---------------------------------------------------------------------------*/

class flameletAnisotropicFilter
:
    public flameletLESfilter
{
    // Private Data

        scalar widthCoeff_;
        volVectorField coeff_;


public:

    //- Runtime type information
    TypeName("anisotropic");

    // Constructors

        //- Construct from components
        flameletAnisotropicFilter(const fvMesh& mesh, scalar widthCoeff);

        //- Construct from IOdictionary
        flameletAnisotropicFilter(const fvMesh& mesh, const dictionary&);

        //- Disallow default bitwise copy construction
        flameletAnisotropicFilter(const flameletAnisotropicFilter&) = delete;


    //- Destructor
    virtual ~flameletAnisotropicFilter()
    {}


    // Member Functions

        //- Read the LESfilter dictionary
        virtual void read(const dictionary&);


    // Member Operators

        void operator=(const flameletAnisotropicFilter&) = delete;

        virtual tmp<volScalarField> operator()
        (
            const tmp<volScalarField>&
        ) const;

        virtual tmp<volVectorField> operator()
        (
            const tmp<volVectorField>&
        ) const;

        virtual tmp<volSymmTensorField> operator()
        (
            const tmp<volSymmTensorField>&
        ) const;

        virtual tmp<volTensorField> operator()
        (
            const tmp<volTensorField>&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
