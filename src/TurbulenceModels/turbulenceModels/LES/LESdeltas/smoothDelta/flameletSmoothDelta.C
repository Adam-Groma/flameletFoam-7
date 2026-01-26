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

\*---------------------------------------------------------------------------*/

#include "flameletSmoothDelta.H"
#include "addToRunTimeSelectionTable.H"
#include "FaceCellWave.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{
    defineTypeNameAndDebug(flameletSmoothDelta, 0);
    addToRunTimeSelectionTable(flameletLESdelta, flameletSmoothDelta, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::LESModels::flameletSmoothDelta::setChangedFaces
(
    const polyMesh& mesh,
    const volScalarField& delta,
    DynamicList<label>& changedFaces,
    DynamicList<deltaData>& changedFacesInfo
)
{
    for (label facei = 0; facei < mesh.nInternalFaces(); facei++)
    {
        scalar ownDelta = delta[mesh.faceOwner()[facei]];

        scalar neiDelta = delta[mesh.faceNeighbour()[facei]];

        // Check if owner delta much larger than neighbour delta or vice versa

        if (ownDelta > maxDeltaRatio_ * neiDelta)
        {
            changedFaces.append(facei);
            changedFacesInfo.append(deltaData(ownDelta));
        }
        else if (neiDelta > maxDeltaRatio_ * ownDelta)
        {
            changedFaces.append(facei);
            changedFacesInfo.append(deltaData(neiDelta));
        }
    }

    // Insert all faces of coupled patches no matter what. Let FaceCellWave
    // sort it out.
    forAll(mesh.boundaryMesh(), patchi)
    {
        const polyPatch& patch = mesh.boundaryMesh()[patchi];

        if (patch.coupled())
        {
            forAll(patch, patchFacei)
            {
                label meshFacei = patch.start() + patchFacei;

                scalar ownDelta = delta[mesh.faceOwner()[meshFacei]];

                changedFaces.append(meshFacei);
                changedFacesInfo.append(deltaData(ownDelta));
            }
        }
    }

    changedFaces.shrink();
    changedFacesInfo.shrink();
}


void Foam::LESModels::flameletSmoothDelta::calcDelta()
{
    const fvMesh& mesh = turbulenceModel_.mesh();

    const volScalarField& geometricDelta = geometricDelta_();

    // Fill changed faces with info
    DynamicList<label> changedFaces(mesh.nFaces()/100 + 100);
    DynamicList<deltaData> changedFacesInfo(changedFaces.size());

    setChangedFaces(mesh, geometricDelta, changedFaces, changedFacesInfo);

    // Set initial field on cells.
    List<deltaData> cellDeltaData(mesh.nCells());

    forAll(geometricDelta, celli)
    {
        cellDeltaData[celli] = geometricDelta[celli];
    }

    // Set initial field on faces.
    List<deltaData> faceDeltaData(mesh.nFaces());


    // Propagate information over whole domain.
    FaceCellWave<deltaData, scalar> deltaCalc
    (
        mesh,
        changedFaces,
        changedFacesInfo,
        faceDeltaData,
        cellDeltaData,
        mesh.globalData().nTotalCells()+1,  // max iterations
        maxDeltaRatio_
    );

    forAll(delta_, celli)
    {
        delta_[celli] = cellDeltaData[celli].delta();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LESModels::flameletSmoothDelta::flameletSmoothDelta
(
    const word& name,
    const flameletTurbulenceModel& turbulence,
    const dictionary& dict
)
:
    flameletLESdelta(name, turbulence),
    geometricDelta_
    (
        flameletLESdelta::New
        (
            "geometricDelta",
            turbulence,
            dict.optionalSubDict(type() + "Coeffs")
        )
    ),
    maxDeltaRatio_
    (
        readScalar
        (
            dict.optionalSubDict(type() + "Coeffs").lookup("maxDeltaRatio")
        )
    )
{
    calcDelta();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::LESModels::flameletSmoothDelta::read(const dictionary& dict)
{
    const dictionary& coeffsDict(dict.optionalSubDict(type() + "Coeffs"));

    geometricDelta_().read(coeffsDict);
    coeffsDict.lookup("maxDeltaRatio") >> maxDeltaRatio_;
    calcDelta();
}


void Foam::LESModels::flameletSmoothDelta::correct()
{
    geometricDelta_().correct();

    if (turbulenceModel_.mesh().changing())
    {
        calcDelta();
    }
}


// ************************************************************************* //
