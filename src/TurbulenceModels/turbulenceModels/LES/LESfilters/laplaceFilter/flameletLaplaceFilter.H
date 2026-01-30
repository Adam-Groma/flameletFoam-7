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
    Foam::flameletLaplaceFilter

Description
    Laplace filter for LES

    \verbatim
    Kernel                 as filter          as Test filter with ratio 2
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Box filter:            g = delta2/24  ->  g = delta2/6
    Spherical box filter:  g = delta2/64  ->  g = delta2/16
    Gaussian filter:       g = delta2/24  ->  g = delta2/6
    \endverbatim

SourceFiles
    flameletLaplaceFilter.C

\*---------------------------------------------------------------------------*/

#ifndef flameletLaplaceFilter_H
#define flameletLaplaceFilter_H

#include "flameletLESfilter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class flameletLaplaceFilter Declaration
\*---------------------------------------------------------------------------*/

class flameletLaplaceFilter
:
    public flameletLESfilter
{
    // Private Data

        scalar widthCoeff_;
        volScalarField coeff_;


public:

    //- Runtime type information
    TypeName("laplace");


    // Constructors

        //- Construct from components
        flameletLaplaceFilter(const fvMesh& mesh, scalar widthCoeff);

        //- Construct from IOdictionary
        flameletLaplaceFilter(const fvMesh& mesh, const dictionary&);

        //- Disallow default bitwise copy construction
        flameletLaplaceFilter(const flameletLaplaceFilter&) = delete;


    //- Destructor
    virtual ~flameletLaplaceFilter()
    {}


    // Member Functions

        //- Read the LESfilter dictionary
        virtual void read(const dictionary&);


    // Member Operators

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


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const flameletLaplaceFilter&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
