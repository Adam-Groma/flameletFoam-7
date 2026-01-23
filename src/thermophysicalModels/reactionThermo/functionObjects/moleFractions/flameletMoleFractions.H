/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2019 OpenFOAM Foundation
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
    Foam::flameletMoleFractions

Description
    This function object calculates mole-fraction fields from the mass-fraction
    fields of the psi/flameletRhoReactionThermo and caches them for output and further
    post-processing.

    The names of the mole-fraction fields are obtained from the corresponding
    mass-fraction fields prepended by "X_"

    Example of function object specification:
    \verbatim
    flameletMoleFractions
    {
        type psiReactionThermoMoleFractions;
    }
    \endverbatim
    or
    \verbatim
    flameletMoleFractions
    {
        type flameletRhoReactionThermoMoleFractions;
    }
    \endverbatim
    depending on the thermodynamics package used in the solver.

See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    flameletMoleFractions.C

\*---------------------------------------------------------------------------*/

#ifndef flameletMoleFractions_H
#define flameletMoleFractions_H

#include "fvMeshFunctionObject.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class flameletMoleFractions Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class flameletMoleFractions
:
    public functionObjects::fvMeshFunctionObject
{
    // Private Data

        //- Species mole fractions
        PtrList<volScalarField> X_;


    // Private Member Functions

        //- Calculate the mole fraction fields
        virtual void calculateMoleFractions();


public:

    //- Runtime type information
    TypeName("flameletMoleFractions");


    // Constructors

        //- Construct from Time and dictionary
        flameletMoleFractions
        (
            const word& name,
            const Time& t,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        flameletMoleFractions(const flameletMoleFractions&) = delete;


    //- Destructor
    virtual ~flameletMoleFractions();


    // Member Functions

        //- Read the flameletMoleFractions data
        virtual bool read(const dictionary&);

        //- Calculate the mole-fraction fields
        virtual bool execute();

        //- The mole-fraction fields auto-write
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const flameletMoleFractions&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "flameletMoleFractions.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
