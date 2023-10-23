/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "MachTurb.H"
#include "fluidThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace functionObjects
    {
        defineTypeNameAndDebug(MachTurb, 0);
        addToRunTimeSelectionTable(functionObject, MachTurb, dictionary);
    }
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::MachTurb::calc()
{
    if (
        foundObject<volScalarField>(fieldName_) && foundObject<fluidThermo>(fluidThermo::dictName))
    {
        const fluidThermo &thermo =
            lookupObject<fluidThermo>(fluidThermo::dictName);

        const volScalarField &k = lookupObject<volScalarField>(fieldName_);

        return store(
            resultName_,
            sqrt(2 * k) / sqrt(thermo.gamma() * thermo.p() / thermo.rho()));
    }

    return false;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::MachTurb::MachTurb(
    const word &name,
    const Time &runTime,
    const dictionary &dict)
    : fieldExpression(name, runTime, dict, "k")
{
    setResultName("MaT", "k");
}

// ************************************************************************* //
