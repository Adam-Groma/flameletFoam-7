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
    Foam::fv::flameletInterRegionExplicitPorositySource

Description
    Inter-region explicit porosity source.

    Sources described by, for example using the DarcyForchheimer model:

    \verbatim
    flameletInterRegionExplicitPorositySourceCoeffs
    {
        type            DarcyForchheimer;
        DarcyForchheimerCoeffs
        {
            d   d [0 -2 0 0 0 0 0] (5e7 -1000 -1000);
            f   f [0 -1 0 0 0 0 0] (0 0 0);

            coordinateSystem
            {
                e1  (0.70710678 0.70710678 0);
                e2  (0 0 1);
            }
        }
    }
    \endverbatim

Note
    The porous region must be selected as a cellZone.

SourceFiles
    flameletInterRegionExplicitPorositySource.C

\*---------------------------------------------------------------------------*/

#ifndef flameletInterRegionExplicitPorositySource_H
#define flameletInterRegionExplicitPorositySource_H

#include "flameletInterRegionOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class porosityModel;

namespace fv
{


/*---------------------------------------------------------------------------*\
                   Class flameletInterRegionExplicitPorositySource Declaration
\*---------------------------------------------------------------------------*/

class flameletInterRegionExplicitPorositySource
:
    public flameletInterRegionOption
{

protected:

    // Protected data

        //- Run-time selectable porosity model
        autoPtr<porosityModel> porosityPtr_;

        //- First iteration
        bool firstIter_;

        //- Velocity field name, default = U
        word UName_;

        //- Dynamic viscosity field name (compressible case only)
        //  default = thermo:mu
        word muName_;


    // Protected Member Functions

        //- Initialise
        void initialise();


public:

    //- Runtime type information
    TypeName("flameletInterRegionExplicitPorositySource");


    // Constructors

        //- Construct from components
        flameletInterRegionExplicitPorositySource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Disallow default bitwise copy construction
        flameletInterRegionExplicitPorositySource
        (
            const flameletInterRegionExplicitPorositySource&
        );


    //- Destructor
    virtual ~flameletInterRegionExplicitPorositySource()
    {}


    // Member Functions

        // Add explicit and implicit contributions

            //- Vector
            virtual void addSup
            (
                fvMatrix<vector>& eqn,
                const label fieldi
            );


        // Add explicit and implicit contributions to compressible equation

            //- Vector
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const label fieldi
            );


        // IO

            //- Read dictionary
            virtual bool read(const dictionary& dict);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const flameletInterRegionExplicitPorositySource&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
