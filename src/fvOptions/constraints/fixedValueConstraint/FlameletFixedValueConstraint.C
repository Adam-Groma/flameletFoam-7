/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "FlameletFixedValueConstraint.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "DimensionedField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::fv::FlameletFixedValueConstraint<Type>::FlameletFixedValueConstraint
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    flameletCellSetOption(name, modelType, dict, mesh)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::fv::FlameletFixedValueConstraint<Type>::read(const dictionary& dict)
{
    if (flameletCellSetOption::read(dict))
    {
        const dictionary& fieldValuesDict = coeffs_.subDict("fieldValues");

        fieldNames_.setSize(fieldValuesDict.size());
        fieldValues_.setSize(fieldNames_.size());

        label i = 0;
        forAllConstIter(dictionary, fieldValuesDict, iter)
        {
            fieldNames_[i] = iter().keyword();
            fieldValuesDict.lookup(iter().keyword()) >> fieldValues_[i];
            i++;
        }

        applied_.setSize(fieldNames_.size(), false);

        return true;
    }
    else
    {
        return false;
    }
}


template<class Type>
void Foam::fv::FlameletFixedValueConstraint<Type>::constrain
(
    fvMatrix<Type>& eqn,
    const label fieldi
)
{
    DebugInfo
        << "FlameletFixedValueConstraint<"
        << pTraits<Type>::typeName
        << ">::constrain for source " << name_ << endl;

    eqn.setValues(cells_, List<Type>(cells_.size(), fieldValues_[fieldi]));
}


// ************************************************************************* //
