/*
 * This file is part of Insight CAE, a workbench for Computer-Aided Engineering 
 * Copyright (C) 2014  Hannes Kroeger <hannes@kroegeronline.net>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 */

#include "fvCFD.H"
#include "OFstream.H"
#include "transformField.H"
#include "transformGeometricField.H"
#include "fixedGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "inflowGeneratorBaseFvPatchVectorField.H"
#include "wallDist.H"
#include "interpolationTable.H"

#include <boost/concept_check.hpp>
#include <boost/assign.hpp>

#include "pipe.h"
#include "channel.h"
// #include "refdata.h"

#include "boostRandomGen.H"
#include "TurbulentSpotVolumeInserter.H"

#include "uniof.h"
// #include "turbulentStructures/anisotropicVorton/anisotropicVorton.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;
using namespace std;
using namespace boost;
using namespace boost::assign;
using namespace insight;


int main(int argc, char *argv[])
{

    argList::validArgs.append("structure type");
    argList::validArgs.append("p0");
    argList::validArgs.append("R");
    argList::validArgs.append("L");
    argList::validOptions.insert("scale", "fluctuation scale factor");
    argList::validOptions.insert("scaleToe", "scale fluctuation to energy");
    
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
    
    wallDist y(mesh);
    volVectorField gradyw=fvc::grad(y.y());
    gradyw/=mag(gradyw)+dimensionedScalar("", gradyw.dimensions(), SMALL);
    gradyw*=y.y();
  
    Foam::word structureType( UNIOF_ADDARG(args,0) );
    Foam::vector p0(IStringStream( UNIOF_ADDARG(args,1) )());
    Foam::symmTensor R(IStringStream( UNIOF_ADDARG(args,2) )());
    Foam::vector L(IStringStream( UNIOF_ADDARG(args,3) )());

    Info << "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    
    autoPtr<turbulentSpotVolumeInserter> ins = turbulentSpotVolumeInserter::New(structureType, mesh, "");
    
    ins->insertSpot(p0, L, 0., R, gradyw[mesh.findCell(p0)]); 
    
    U=ins->uPrime();
    
    if (args.optionFound("scale"))
    {
      U*=readScalar(IStringStream(args.options()["scale"])());
    }
    
    dimensionedScalar e=fvc::domainIntegrate( magSqr(U) );
    Info<<"total kinetic energy e = "<<e<<endl;

    if (args.optionFound("scaleToe"))
    {
      scalar target_e=readScalar(IStringStream(args.options()["scaleToe"])());
      U *= Foam::sqrt(target_e/e.value());
      
      dimensionedScalar e_new=fvc::domainIntegrate( magSqr(U) );
      Info<<"after rescale: total kinetic energy e = "<<e_new<<endl;
    }
    
/*
    anisotropicVorton::StructureParameters sp;
    BoostRandomGen rg;
    anisotropicVorton aiv(rg, p0, Foam::vector::zero, Foam::vector::zero, L, 0.0, -1, R);

    forAll(mesh.C(), ci)
    {
      U[ci]+=aiv.fluctuation(sp, mesh.C()[ci]);
    }
    */
    
    U.write();

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
