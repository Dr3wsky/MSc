/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "kOmegaSSTdd.H"
#include "MachNo.H" //
#include "fluidThermo.H"
#include "addToRunTimeSelectionTable.H" //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    namespace RASModels
    {

        // * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

        template <class BasicTurbulenceModel>
        void kOmegaSSTdd<BasicTurbulenceModel>::correctNut(const volScalarField &S2)
        {
            // Correct the turbulence viscosity
            kOmegaSSTBase<eddyViscosity<RASModel<BasicTurbulenceModel>>>::correctNut(
                S2);

            // Correct the turbulence thermal diffusivity
            BasicTurbulenceModel::correctNut();
        }

        template <class BasicTurbulenceModel>
        void kOmegaSSTdd<BasicTurbulenceModel>::correctNut()
        {
            correctNut(2 * magSqr(symm(fvc::grad(this->U_))));
        }

        // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

        template <class BasicTurbulenceModel>
        kOmegaSSTdd<BasicTurbulenceModel>::kOmegaSSTdd(
            const alphaField &alpha,
            const rhoField &rho,
            const volVectorField &U,
            const surfaceScalarField &alphaRhoPhi,
            const surfaceScalarField &phi,
            const transportModel &transport,
            const word &propertiesName,
            const word &type)
            : kOmegaSSTBase<eddyViscosity<RASModel<BasicTurbulenceModel>>>(
                  type,
                  alpha,
                  rho,
                  U,
                  alphaRhoPhi,
                  phi,
                  transport,
                  propertiesName),
              MachTurb_(
                  IOobject(
                      "MachTurb",
                      this->mesh().time().timeName(),
                      this->mesh(),
                      IOobject::NO_READ,
                      IOobject::AUTO_WRITE),
                  this->mesh(),
                  dimensionedScalar(dimless, Zero)),
              gammaThermo_(
                  IOobject(
                      "gammaThermo",
                      this->mesh().time().timeName(),
                      this->mesh(),
                      IOobject::NO_READ,
                      IOobject::AUTO_WRITE),
                  this->mesh(),
                  dimensionedScalar(dimless, Zero))
        {
            if (type == typeName)
            {
                this->printCoeffs(type);
            }
        }

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
        // Used to update turbulent Mach number
        template <class BasicTurbulenceModel>
        void kOmegaSSTdd<BasicTurbulenceModel>::correctMachTurb()
        {
            const fluidThermo &thermo = this->mesh().objectRegistry::lookupObject<fluidThermo>(fluidThermo::dictName);
            // const fluidThermo& thermo = Foam::functionObjects::MachNo::lookupObject<fluidThermo>(fluidThermo::dictName);
            // const fluidThermo& thermo =  Foam::functionObjects::MachNo::lookupObject<fluidThermo>(fluidThermo::dictName);
            gammaThermo_ = thermo.gamma();
            MachTurb_ = sqrt(2 * this->k_) / sqrt(gammaThermo_ * thermo.p() / this->rho_);
        }

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
        template <class BasicTurbulenceModel>
        void kOmegaSSTdd<BasicTurbulenceModel>::correct()
        {
            if (!this->turbulence_)
            {
                return;
            }

            // Output statement to show turbulence model being used
            Info << "-----------------------------------------------------------------------------" << endl;
            Info << "This is a modified version of kOmegaSST to account for dilatation dissipation" << endl;
            Info << "-----------------------------------------------------------------------------" << endl;

            // Local references
            const alphaField &alpha = this->alpha_;
            const rhoField &rho = this->rho_;
            const surfaceScalarField &alphaRhoPhi = this->alphaRhoPhi_;
            const volVectorField &U = this->U_;
            volScalarField &nut = this->nut_;
            fv::options &fvOptions(fv::options::New(this->mesh_));

            // Updating turbulent Mach No
            correctMachTurb();

            BasicTurbulenceModel::correct();

            volScalarField::Internal divU(fvc::div(fvc::absolute(this->phi(), U)));

            tmp<volTensorField> tgradU = fvc::grad(U);
            volScalarField S2(2 * magSqr(symm(tgradU())));
            volScalarField::Internal GbyNu0(
                this->type() + ":GbyNu",
                (tgradU() && dev(twoSymm(tgradU()))));
            volScalarField::Internal G(this->GName(), nut * GbyNu0);

            // Update omega and G at the wall
            this->omega_.boundaryFieldRef().updateCoeffs();

            volScalarField CDkOmega(
                (2 * this->alphaOmega2_) * (fvc::grad(this->k_) & fvc::grad(this->omega_)) / this->omega_);

            volScalarField F1(this->F1(CDkOmega));
            volScalarField F23(this->F23());

            {
                volScalarField::Internal gamma(this->gamma(F1));
                volScalarField::Internal beta(this->beta(F1));

                GbyNu0 = this->GbyNu(GbyNu0, F23(), S2());

                // Turbulent frequency equation
                tmp<fvScalarMatrix> omegaEqn(
                    fvm::ddt(alpha, rho, this->omega_) + fvm::div(alphaRhoPhi, this->omega_) - fvm::laplacian(alpha * rho * this->DomegaEff(F1), this->omega_) ==
                    alpha() * rho() * gamma * GbyNu0 - fvm::SuSp((2.0 / 3.0) * alpha() * rho() * gamma * divU, this->omega_) - fvm::Sp(alpha() * rho() * beta * this->omega_(), this->omega_) - fvm::SuSp(alpha() * rho() * (F1() - scalar(1)) * CDkOmega() / this->omega_(), this->omega_) + alpha() * rho() * beta * sqr(this->omegaInf_) + this->Qsas(S2(), gamma, beta) + this->omegaSource() + fvOptions(alpha, rho, this->omega_));

                omegaEqn.ref().relax();
                fvOptions.constrain(omegaEqn.ref());
                omegaEqn.ref().boundaryManipulate(this->omega_.boundaryFieldRef());
                solve(omegaEqn);
                fvOptions.correct(this->omega_);
                bound(this->omega_, this->omegaMin_);
            }

            // Turbulent kinetic energy equation
            tmp<fvScalarMatrix> kEqn(
                fvm::ddt(alpha, rho, this->k_) + fvm::div(alphaRhoPhi, this->k_) - fvm::laplacian(alpha * rho * this->DkEff(F1), this->k_) ==
                alpha() * rho() * this->Pk(G) - fvm::SuSp((2.0 / 3.0) * alpha() * rho() * divU, this->k_) - fvm::Sp(alpha() * rho() * this->epsilonByk(F1, tgradU()), this->k_) + alpha() * rho() * this->betaStar_ * this->omegaInf_ * this->kInf_ + this->kSource() + fvOptions(alpha, rho, this->k_));

            tgradU.clear();

            kEqn.ref().relax();
            fvOptions.constrain(kEqn.ref());
            solve(kEqn);
            fvOptions.correct(this->k_);
            bound(this->k_, this->kMin_);

            correctNut(S2);
        }
    } // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
