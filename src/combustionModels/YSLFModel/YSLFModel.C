/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License

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

Contributors/Copyright
    2014 Hagen Müller <hagen.mueller@unibw.de> Universität der Bundeswehr München
    2014 Likun Ma <L.Ma@tudelft.nl> TU Delft

\*---------------------------------------------------------------------------*/

#include "YSLFModel.H"
#include "fvmSup.H"
#include "hashedWordList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ReactionThermo>
Foam::combustionModels::YSLFModel<ReactionThermo>::YSLFModel
(
    const word& modelType,
    ReactionThermo& thermo,
    const compressibleTurbulenceModel& turb,
    const word& combustionProperties
)
:
    ThermoCombustion<ReactionThermo>(modelType, thermo, turb),
    solver_(tableSolver(this->mesh(), tables())),
    Y_(this->thermo().composition().Y()),
    he_(this->thermo().he()),             
    Z_(this->thermo().Z()),
    varZ_(this->thermo().varZ()),
    Chi_(this->thermo().Chi()),
    ubIF_(this->mesh().cells().size()),
    ubP_(),
    posIF_(this->mesh().cells().size()),
    posP_(),
    useScalarDissipation_(this->coeffs().lookup("useScalarDissipation")),
    useMixtureFractionVariance_(this->coeffs().lookup("useMixtureFractionVariance"))
{
  const polyBoundaryMesh& patches = this->mesh().boundaryMesh();
  int patchSize = 0;
  forAll(patches, patchI)
  {
    const polyPatch& pp = patches[patchI];
    if (pp.size() > patchSize) patchSize = pp.size();

    ubP_.setSize(patchSize);
    posP_.setSize(patchSize);
  }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ReactionThermo>
Foam::combustionModels::YSLFModel<ReactionThermo>::~YSLFModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ReactionThermo>
Foam::hashedWordList
Foam::combustionModels::YSLFModel<ReactionThermo>::tables()
{
  hashedWordList tableNames = this->thermo().composition().species();
  tableNames.append("he");

  return tableNames;
}

template<class ReactionThermo>
void Foam::combustionModels::YSLFModel<ReactionThermo>::correct()
{
    // Limit the scalar dissipation rate to avoid instabilities at extiction
    scalar chiLimiter = solver_.maxChi();

    const scalarField& ZCells = Z_.internalField();
    const scalarField& varZCells = varZ_.internalField();
    const scalarField& chiCells = Chi_.internalField();

    scalarField& heCells = he_.primitiveFieldRef();

    //- Update the species and enthalpy field
    scalarList x(3, 0.0);
    double Zeta;

    // Interpolate for internal Field
    forAll(Y_, i)
    {
        scalarField& YCells = Y_[i].primitiveFieldRef();

        forAll(ZCells, cellI)
        {
    	      if (i == 0)
    	      {
    	  	      Zeta = sqrt(varZCells[cellI]/max(ZCells[cellI]*(1 - ZCells[cellI]), SMALL));
                if (useScalarDissipation_)   x[0] = min(chiCells[cellI], chiLimiter);
                if (useMixtureFractionVariance_) x[1] = min(Zeta, 0.99);
                x[2] = ZCells[cellI];

                ubIF_[cellI] = solver_.upperBounds(x);
                posIF_[cellI] = solver_.position(ubIF_[cellI], x);

          	    heCells[cellI] = solver_.interpolate(ubIF_[cellI], posIF_[cellI], (solver_.sizeTableNames() - 1));
    	      }

    	      YCells[cellI] = solver_.interpolate(ubIF_[cellI], posIF_[cellI], i);
        }
    }

    // Interpolate for patches
    forAll(he_.boundaryField(), patchi)   
    {
       const fvPatchScalarField& pChi = Chi_.boundaryField()[patchi];
       const fvPatchScalarField& pvarZ = varZ_.boundaryField()[patchi];
       const fvPatchScalarField& pZ = Z_.boundaryField()[patchi];

       fvPatchScalarField& pHe = he_.boundaryFieldRef()[patchi];

       forAll(Y_, i)
       {
     	  fvPatchScalarField& pY = Y_[i].boundaryFieldRef()[patchi];

           forAll(pY , facei)
           {
          	 if (i == 0)
          	 {
                  Zeta = sqrt(pvarZ[facei]/max(pZ[facei]*(1 - pZ[facei]), SMALL));

                  if (useScalarDissipation_) x[0] = min(pChi[facei], chiLimiter);
                  if (useMixtureFractionVariance_) x[1] = min(Zeta, 0.99);
                  x[2] = pZ[facei];

                  ubP_[facei] = solver_.upperBounds(x);
                  posP_[facei] = solver_.position(ubP_[facei], x);

                  pHe[facei] = solver_.interpolate(ubP_[facei], posP_[facei], (solver_.sizeTableNames() - 1));
          	 }

         	 pY[facei] = solver_.interpolate(ubP_[facei], posP_[facei], i);
          }
       }
    }

    // Calculate thermodynamic Properties
    this->thermo().correct();
}


template<class ReactionThermo>
Foam::Switch
Foam::combustionModels::YSLFModel<ReactionThermo>::correctDensity()
{
    return true;
}


template<class ReactionThermo>
Foam::tmp<Foam::fvScalarMatrix>
Foam::combustionModels::YSLFModel<ReactionThermo>::R
(
    volScalarField& Y
) const
{
    tmp<fvScalarMatrix> tSu
    (
        new fvScalarMatrix(Y, dimMass/dimTime)
    );

    return tSu;
}


template<class ReactionThermo>
Foam::tmp<Foam::volScalarField>
Foam::combustionModels::YSLFModel<ReactionThermo>::Qdot() const
{
    return volScalarField::New
    (
        this->thermo().phasePropertyName(typeName + ":Qdot"),
        this->mesh(),
        dimensionedScalar(dimEnergy/dimVolume/dimTime, 0)
    );
}


template<class ReactionThermo>
bool Foam::combustionModels::YSLFModel<ReactionThermo>::read()
{
    if (ThermoCombustion<ReactionThermo>::read())
    {
        this->coeffs().lookup("useScalarDissipation") >> useScalarDissipation_;
        this->coeffs().lookup("useMixtureFractionVariance") >> useMixtureFractionVariance_;
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //

