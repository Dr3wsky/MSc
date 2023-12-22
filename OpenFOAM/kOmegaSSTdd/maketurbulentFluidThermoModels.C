/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
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

#include "turbulentFluidThermoModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// -------------------------------------------------------------------------- //
// Compressibility Correction  models
// -------------------------------------------------------------------------- //

// Dilatation Dissipation - modified dissipation term with MachT func, Sarkar
#include "kOmegaSSTdd.H"
makeRASModel(kOmegaSSTdd);

// Pressur Dissipation - additional term in both k and omega
#include "kOmegaSSTpd.H"
makeRASModel(kOmegaSSTpd);

// Compressible Dissipation - combination of both DD and PD
#include "kOmegaSSTcd.H"
makeRASModel(kOmegaSSTcd);

// ************************************************************************* //
