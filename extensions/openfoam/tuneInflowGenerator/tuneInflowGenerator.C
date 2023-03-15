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

#include "uniof.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;
using namespace std;
using namespace boost;
using namespace boost::assign;
using namespace insight;


int main(int argc, char *argv[])
{
#ifdef OF16ext
    argList::validOptions.insert("nSteps", "number of time steps to run");
    argList::validOptions.insert("patches", "patches to test");
    argList::validOptions.insert("doFullTimeLoop", "");
#else
    argList::addOption("nSteps", "number of time steps to run");
    argList::addOption("patches", "patches to test");
    argList::addOption("doFullTimeLoop", "");
#endif

//   argList::validArgs.append("patches");

#   include "setRootCase.H"
#   include "createTime.H"

    label nsteps=10;
    bool do_fulltimeloop=false;

#ifdef OF16ext
#define PATCHNAMELIST wordList
#else
#define PATCHNAMELIST wordReList
#endif

    autoPtr<PATCHNAMELIST> patchNames;

    IOobject dicthead
    (
        "checkInflowGeneratorStatisticsDict",
        runTime.system(),
        runTime,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );
    if (UNIOF_HEADEROK(dicthead,IOdictionary))
    {
        IOdictionary dict(dicthead);
        if (dict.found("nSteps")) nsteps=readLabel(dict.lookup("nSteps"));
        if (dict.found("doFullTimeLoop")) do_fulltimeloop=Switch(dict.lookup("doFullTimeLoop"));
        patchNames.reset(new PATCHNAMELIST(dict.lookup("patches")));
    }

    if (args.optionFound("nSteps"))
    {
        nsteps=readLabel(IStringStream(args.option("nSteps"))());
    }

    if (args.optionFound("doFullTimeLoop"))
    {
        do_fulltimeloop=Switch(IStringStream(args.option("doFullTimeLoop"))());
    }

    if (args.optionFound("patches"))
    {
        patchNames.reset(new PATCHNAMELIST(IStringStream(args.option("patches"))()));
    }

    if (!patchNames.valid())
    {
        FatalErrorIn("main")
                <<"patches to test are not specified!\n"
                <<"Either provide a list on the command line or in the system/checkInflowGeneratorStatisticsDict!\n"
                <<abort(FatalError);
    }

#   include "createMesh.H"

    labelHashSet patches =
        mesh.boundaryMesh().patchSet(patchNames());

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

    forAllConstIter(labelHashSet, patches, iter)
    {
        label patchI=iter.key();
        if (isA<inflowGeneratorBaseFvPatchVectorField>(U.boundaryField()[patchI]))
        {
            inflowGeneratorBaseFvPatchVectorField &ifpf =
                refCast<inflowGeneratorBaseFvPatchVectorField>(UNIOF_BOUNDARY_NONCONST(U)[patchI]);

            ifpf.computeConditioningFactor(1, nsteps);
        }
    }

    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
