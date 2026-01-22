/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    translates cantera table Data to OpenFOAM

Contributors/Copyright
    2014 Hagen Müller <hagen.mueller@unibw.de> Universität der Bundeswehr München
    2014 Gabriele Frank <gabriele.frank@unibw.de> Universität der Bundeswehr München

\*---------------------------------------------------------------------------*/

#include "OFstream.H"
#include "fvCFD.H"
//#include "rhoCombustionModel.H"
//#include "rhoChemistryModel.H"  
#include "turbulenceModel.H"
#include "multivariateScheme.H"
#include "pimpleControl.H"
#include "canteraReader.H"
#include "flameletBasicSpecieMixture.H"
#include "flameletRhoReactionThermo.H"
#include <stdio.h>
#include <stdlib.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
   #include "setRootCase.H"
   #include "createTime.H"
   #include "createMesh.H"

   IOdictionary tableDict
   (
       IOobject
       (
          "tableProperties",
          runTime.constant(),
          runTime,
          IOobject::MUST_READ,
          IOobject::NO_WRITE,
          false
       )
   );

   Info<< "Reading thermophysical properties\n" << endl;

   autoPtr<flameletRhoReactionThermo> pThermo(flameletRhoReactionThermo::New(mesh));

   flameletRhoReactionThermo& thermo = pThermo();
   flameletBasicSpecieMixture& composition = thermo.composition();

   Info << "Creating dummy tables\n" << endl;
   //create dummy tables
   hashedWordList dummytable(thermo.composition().species());
//   hashedWordList dummytable(List<word>{
//		   "T", "H2", "H", "O", "O2", "OH", "H2O", "HO2", "H2O2",
//		   "C", "CH", "CH2", "CH2(S)", "CH3", "CH4", "CO", "CO2",
//		   "HCO", "CH2O", "CH2OH", "CH3O", "CH3OH", "C2H", "C2H2",
//		   "C2H3", "C2H4", "C2H5", "C2H6", "HCCO", "CH2CO", "HCCOH",
//		   "N", "NH", "NH2", "NH3", "NNH", "NO", "NO2", "N2O", "HNO",
//		   "CN", "HCN", "H2CN", "HCNN", "HCNO", "HOCN", "HNCO", "NCO",
//		   "N2", "AR", "C3H7", "C3H8", "CH2CHO", "CH3CHO", "Z"
//		   });
   List<scalar> dummy(0);

   for (int i=0; i<dummytable.size(); i++)
   {
      IOdictionary dictionary
      (
         IOobject
         (
            dummytable[i]+"_table",
            runTime.constant(),
            runTime,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
         )
      );

      word dummyName=dummytable[i]+"_table";
      dictionary.set(dummyName, dummy);
      OFstream output("constant/"+dummytable[i]+"_table");
      dictionary.writeHeader(output);
      dictionary.writeData(output);
   }
   //end creating dummy tables

   //read the cantera-data
   canteraReader canteraRead(tableDict, thermo, composition);

   Info<<"lines = "<<canteraRead.numberOfLines()<<endl;
   Info<<"columns = "<<canteraRead.numberOfColumns()<<endl;

   /*Write the tables, constant/tables/ must exist*/
   for (int i=0; i<canteraRead.getNames().size(); i++)
   {
      IOdictionary dictionary
      (
         IOobject
         (
            canteraRead.getNames()[i],
            runTime.constant(),
            runTime,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
         )
      );

      OFstream output("constant/"+canteraRead.getNames()[i]+"_table");
      canteraRead.write(i, dictionary, output);
   }
   
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
Info<< "End\n" << endl;

return 0;
}

