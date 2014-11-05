/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Application
    Co

Description
    Calculates and writes the Co number as a volScalarField obtained
    from field phi.

    The -noWrite option just outputs the max values without writing the
    field.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void calcCompressibleTotalKineticEnergy
(
    const fvMesh& mesh,
    const Time& runTime,
    const volVectorField& U
)
{

    IOobject rhoHeader
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (!rhoHeader.headerOk())
    {
	Info << "    no " << rhoHeader.name() << " field\n" << endl;
	return;
    }

    Info << "Reading field rho\n" << endl;
    volScalarField rho(rhoHeader, mesh);

    dimensionedScalar totalKE("totalKE", dimensionSet(1,2,-2,0,0,0,0), gSum(0.5*(U&U)*mesh.V()*rho));	
    Info << "    Total kinetic energy: " << totalKE.value() << " [J]\n" << endl;
}

void calcIncompressibleTotalKineticEnergy
(
    const fvMesh& mesh,
    const Time& runTime,
    const volVectorField& U
)
{
    Info<< "Reading transportProperties dictionary \n" << endl;

    IOdictionary transportProperties
    (
	IOobject
	(
	    "transportProperties",
	    runTime.constant(),
	    mesh,
	    IOobject::MUST_READ,
	    IOobject::NO_WRITE
	)
    );

    dimensionedScalar rhoConst (transportProperties.lookup("rho"));

    dimensionedScalar totalKE("totalKE", dimensionSet(1,2,-2,0,0,0,0), gSum(0.5*(U&U)*mesh.V()*rhoConst.value()));
    
    if (totalKE.dimensions() == ( U.dimensions()* U.dimensions() * rhoConst.dimensions() * mesh.V().dimensions()))
    {
	Info << "    Total kinetic energy: " << totalKE.value() << " [J]\n" << endl;
    }
    else
    {
	FatalError
		<<"Incorrect dimensions of totalKE: " << U.dimensions()* U.dimensions() * rhoConst.dimensions() * mesh.V().dimensions()
		<< " should be " << totalKE.dimensions() << endl << endl
		<< "Dimensions in calculation are:" << endl
		<< "U   " << U.dimensions() << endl
		<< "rho " << rhoConst.dimensions() << endl
                << "V   " << mesh.V().dimensions() << endl << exit(FatalError);
    }
}

int main(int argc, char *argv[])
{
    timeSelector::addOptions();

    #include "addRegionOption.H"

    argList::addBoolOption
    (
	"compressible",
	"calculate kinetic energy for compressible cases"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createNamedMesh.H"

    const bool compressible = args.optionFound("compressible");

    forAll(timeDirs, timeI)
    {
 	runTime.setTime(timeDirs[timeI], timeI);
	Info<< "Time = " << runTime.timeName() << endl;

    	IOobject Uheader 
	(
            "U",
	    runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
	    IOobject::NO_WRITE
	);

        if (Uheader.headerOk())
        {
	    Info << "Reading field U" << endl;
	    volVectorField U(Uheader, mesh);

    	    if (compressible)
	    {
		calcCompressibleTotalKineticEnergy(mesh, runTime, U);
	    }
	    else
	    {
		calcIncompressibleTotalKineticEnergy(mesh, runTime, U);
	    }
	}
    }

    Info << "End\n" << endl;
    return 0;
}

// ************************************************************************* //
