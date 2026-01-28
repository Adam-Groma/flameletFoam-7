/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019 OpenFOAM Foundation
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
    Foam::flameletSpecieCoeffs

Description
    Class to hold the specie index and its coefficients in the
    reaction rate expression.

SourceFiles
    flameletSpecieCoeffs.C

\*---------------------------------------------------------------------------*/

#ifndef flameletSpecieCoeffs_H
#define flameletSpecieCoeffs_H

#include "flameletSpeciesTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
class flameletSpecieCoeffs;
Ostream& operator<<(Ostream&, const flameletSpecieCoeffs&);


/*---------------------------------------------------------------------------*\
                         Class flameletSpecieCoeffs Declaration
\*---------------------------------------------------------------------------*/

class flameletSpecieCoeffs
{
public:

    label index;
    scalar stoichCoeff;
    scalar exponent;

    flameletSpecieCoeffs()
    :
        index(-1),
        stoichCoeff(0),
        exponent(1)
    {}

    flameletSpecieCoeffs(const flameletSpeciesTable& species, Istream& is);

    bool operator==(const flameletSpecieCoeffs& sc) const
    {
        return index == sc.index;
    }

    bool operator!=(const flameletSpecieCoeffs& sc) const
    {
        return index != sc.index;
    }

    friend Ostream& operator<<(Ostream& os, const flameletSpecieCoeffs& sc)
    {
        os  << sc.index << token::SPACE
            << sc.stoichCoeff << token::SPACE
            << sc.exponent;
        return os;
    }

    //- Construct the left- and right-hand-side reaction coefficients
    static void setLRhs
    (
        Istream&,
        const flameletSpeciesTable&,
        List<flameletSpecieCoeffs>& lhs,
        List<flameletSpecieCoeffs>& rhs
    );

    //- Write the string representation of the flameletSpecieCoeffs list
    static void reactionStr
    (
        OStringStream& reaction,
        const flameletSpeciesTable&,
        const List<flameletSpecieCoeffs>& scs
    );

    //- Return string representation of reaction
    static string reactionStr
    (
        OStringStream& reaction,
        const flameletSpeciesTable&,
        const List<flameletSpecieCoeffs>& lhs,
        const List<flameletSpecieCoeffs>& rhs
    );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
