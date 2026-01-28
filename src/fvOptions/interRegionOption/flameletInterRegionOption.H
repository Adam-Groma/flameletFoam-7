/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::fv::flameletInterRegionOption

Description
    Base class for inter-region exchange.

\*---------------------------------------------------------------------------*/

#ifndef flameletInterRegionOption_H
#define flameletInterRegionOption_H

#include "fvOption.H"
#include "volFields.H"
#include "autoPtr.H"
#include "meshToMesh.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                Class flameletInterRegionOption Declaration
\*---------------------------------------------------------------------------*/

class flameletInterRegionOption
:
    public option
{
protected:

    // Protected data

        //- Master or slave region
        bool master_;

        //- Name of the neighbour region to map
        word nbrRegionName_;

        //- Mesh to mesh interpolation object
        autoPtr<meshToMesh> meshInterpPtr_;


    // Protected member functions

        //- Set the mesh to mesh interpolation object
        void setMapper();


public:

    //- Runtime type information
    TypeName("flameletInterRegionOption");


    // Constructors

        //- Construct from dictionary
        flameletInterRegionOption
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~flameletInterRegionOption();


    // Member Functions

        // Access

            //- Return const access to the neighbour region name
            inline const word& nbrRegionName() const;

            //- Return const access to the mapToMap pointer
            inline const meshToMesh& meshInterp() const;


        // IO

            //- Read dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "flameletInterRegionOptionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
