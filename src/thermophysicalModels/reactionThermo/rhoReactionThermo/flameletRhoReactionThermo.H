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
    Foam::flameletRhoReactionThermo

Description
    Foam::flameletRhoReactionThermo

SourceFiles
    flameletRhoReactionThermo.C

\*---------------------------------------------------------------------------*/

#ifndef flameletRhoReactionThermo_H
#define flameletRhoReactionThermo_H

#include "flameletRhoThermo.H"
#include "flameletBasicSpecieMixture.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class flameletRhoReactionThermo Declaration
\*---------------------------------------------------------------------------*/

class flameletRhoReactionThermo
:
    public flameletRhoThermo
{

public:

    //- Runtime type information
    TypeName("flameletRhoReactionThermo");


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        flameletRhoReactionThermo,
        fvMesh,
        (const fvMesh& mesh, const word& phaseName),
        (mesh, phaseName)
    );


    // Constructors

        //- Construct from mesh and phase name
        flameletRhoReactionThermo
        (
            const fvMesh&,
            const word& phaseName
        );


    // Selectors

        //- Standard selection based on fvMesh
        static autoPtr<flameletRhoReactionThermo> New
        (
            const fvMesh&,
            const word& phaseName=word::null
        );


    //- Destructor
    virtual ~flameletRhoReactionThermo();


    // Member Functions

        //- Return the composition of the multi-component mixture
        virtual flameletBasicSpecieMixture& composition() = 0;

        //- Return the composition of the multi-component mixture
        virtual const flameletBasicSpecieMixture& composition() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
