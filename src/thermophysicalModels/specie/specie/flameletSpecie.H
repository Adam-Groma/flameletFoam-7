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
    Foam::flameletSpecie

Description
    Base class of the thermophysical property types.

SourceFiles
    flameletSpecieI.H
    flameletSpecie.C

\*---------------------------------------------------------------------------*/

#ifndef flameletSpecie_H
#define flameletSpecie_H

#include "word.H"
#include "scalar.H"
#include "dictionary.H"

#include "thermodynamicConstants.H"
using namespace Foam::constant::thermodynamic;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class flameletSpecie;

inline flameletSpecie operator+(const flameletSpecie&, const flameletSpecie&);
inline flameletSpecie operator*(const scalar, const flameletSpecie&);
inline flameletSpecie operator==(const flameletSpecie&, const flameletSpecie&);

Ostream& operator<<(Ostream&, const flameletSpecie&);


/*---------------------------------------------------------------------------*\
                           Class flameletSpecie Declaration
\*---------------------------------------------------------------------------*/

class flameletSpecie
{
    // Private Data

        //- Name of flameletSpecie
        word name_;

        //- Number of moles of this component in the mixture
        scalar Y_;

        //- Molecular weight of flameletSpecie [kg/kmol]
        scalar molWeight_;


public:

    //- Runtime type information
    ClassName("flameletSpecie");


    // Constructors

        //- Construct from components without name
        inline flameletSpecie(const scalar Y, const scalar molWeight);

        //- Construct from components with name
        inline flameletSpecie
        (
            const word& name,
            const scalar Y,
            const scalar molWeight
        );

        //- Construct as named copy
        inline flameletSpecie(const word& name, const flameletSpecie&);

        //- Construct from dictionary
        flameletSpecie(const dictionary& dict);

        //- Copy constructor
        flameletSpecie(const flameletSpecie&) = default;


    // Member Functions

        // Access

            //- Name
            inline const word& name() const;

            //- Molecular weight [kg/kmol]
            inline scalar W() const;

            //- No of moles of this species in mixture
            inline scalar Y() const;

            //- Gas constant [J/kg/K]
            inline scalar R() const;


        // I-O

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member Operators

        inline void operator=(const flameletSpecie&);
        inline void operator+=(const flameletSpecie&);
        inline void operator*=(const scalar);


    // Friend operators

        inline friend flameletSpecie operator+(const flameletSpecie&, const flameletSpecie&);
        inline friend flameletSpecie operator*(const scalar, const flameletSpecie&);
        inline friend flameletSpecie operator==(const flameletSpecie&, const flameletSpecie&);


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const flameletSpecie&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "flameletSpecieI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
