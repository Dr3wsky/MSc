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

#include "kOmegaSSTdilationDissipation.H"
#include "MachNo.H"
#include "fluidThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    namespace functionObjects
    {
        defineTypeNameAndDebug(MachNo, 0);
        addToRunTimeSelectionTable(functionObject, MachNo, dictionary);

        bool MachNo::calc()
        {
            if (
                foundObject<volVectorField>(fieldName_) && foundObject<fluidThermo>(fluidThermo::dictName))
            {
                const fluidThermo &thermo =
                    lookupObject<fluidThermo>(fluidThermo::dictName);

                const volVectorField &U = lookupObject<volVectorField>(fieldName_);

                return store(
                    resultName_,
                    mag(U) / sqrt(thermo.gamma() * thermo.p() / thermo.rho()));
            }

            return false;
        }

        MachNo::MachNo(
            const word &name,
            const Time &runTime,
            const dictionary &dict)
            : fieldExpression(name, runTime, dict, "U")
        {
            setResultName("Ma", "U");
        }
    } // End namespace functionObjects

    // * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //w

    // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
    namespace RASModels
    {

        // * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

        template <class BasicTurbulenceModel>
        void kOmegaSSTdilationDissipation<BasicTurbulenceModel>::correctNut(const volScalarField &S2)
        {
            // Correct the turbulence viscosity
            kOmegaSSTBase<eddyViscosity<RASModel<BasicTurbulenceModel>>>::correctNut(
                S2);

            // Correct the turbulence thermal diffusivity
            BasicTurbulenceModel::correctNut();
        }

        template <class BasicTurbulenceModel>
        void kOmegaSSTdilationDissipation<BasicTurbulenceModel>::correctNut()
        {
            correctNut(2 * magSqr(symm(fvc::grad(this->U_))));
        }

        // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

        template <class BasicTurbulenceModel>
        kOmegaSSTdilationDissipation<BasicTurbulenceModel>::kOmegaSSTdilationDissipation(
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
                  propertiesName)
        {
            if (type == typeName)
            {
                this->printCoeffs(type);
            }
        }

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
        template <class BasicTurbulenceModel>
        void kOmegaSSTdilationDissipation<BasicTurbulenceModel>::correct()
        {
            if (!this->turbulence_)
            {
                return;
            }

            // Local references
            const alphaField &alpha = this->alpha_;
            const rhoField &rho = this->rho_;
            const surfaceScalarField &alphaRhoPhi = this->alphaRhoPhi_;
            const volVectorField &U = this->U_;
            // volScalarField &MachTest = this->MachTest_;
            volScalarField &omega = this->omega_;
            volScalarField &k = this->k_;
            volScalarField &nut = this->nut_;

            fv::options &
                fvOptions(fv::options::New(this->mesh_));

            // Output statements to confirm use
            Info << "----------------------------------------------" << endl;
            Info << "RUNNING kOmegaSSTdilationDissipation" << endl;
            Info << "----------------------------------------------" << endl;

            BasicTurbulenceModel::correct();

            volScalarField::Internal divU(fvc::div(fvc::absolute(this->phi(), U)));

            tmp<volTensorField> tgradU = fvc::grad(U);
            volScalarField S2(2 * magSqr(symm(tgradU())));
            volScalarField::Internal GbyNu0(
                this->type() + ":GbyNu",
                (tgradU() && dev(twoSymm(tgradU()))));
            volScalarField::Internal G(this->GName(), nut * GbyNu0);

            // Update omega and G at the wall
            omega.boundaryFieldRef().updateCoeffs();

            volScalarField CDkOmega(
                (2 * this->alphaOmega2_) * (fvc::grad(k) & fvc::grad(omega)) / omega);

            volScalarField F1(this->F1(CDkOmega));
            volScalarField F23(this->F23());

            {
                volScalarField::Internal gamma(this->gamma(F1));
                volScalarField::Internal beta(this->beta(F1));

                GbyNu0 = this->GbyNu(GbyNu0, F23(), S2());

                // Turbulent frequency equation
                tmp<fvScalarMatrix> omegaEqn(
                    fvm::ddt(alpha, rho, omega) + fvm::div(alphaRhoPhi, omega) - fvm::laplacian(alpha * rho * this->DomegaEff(F1), omega) ==
                    alpha() * rho() * gamma * GbyNu0 - fvm::SuSp((2.0 / 3.0) * alpha() * rho() * gamma * divU, omega) - fvm::Sp(alpha() * rho() * beta * omega(), omega) - fvm::SuSp(alpha() * rho() * (F1() - scalar(1)) * CDkOmega() / omega(), omega) + alpha() * rho() * beta * sqr(this->omegaInf_) + this->Qsas(S2(), gamma, beta) + this->omegaSource() + fvOptions(alpha, rho, omega));

                omegaEqn.ref().relax();
                fvOptions.constrain(omegaEqn.ref());
                omegaEqn.ref().boundaryManipulate(omega.boundaryFieldRef());
                solve(omegaEqn);
                fvOptions.correct(omega);
                bound(omega, this->omegaMin_);
            }

            // Turbulent kinetic energy equation
            tmp<fvScalarMatrix> kEqn(
                fvm::ddt(alpha, rho, k) + fvm::div(alphaRhoPhi, k) - fvm::laplacian(alpha * rho * this->DkEff(F1), k) ==
                alpha() * rho() * this->Pk(G) - fvm::SuSp((2.0 / 3.0) * alpha() * rho() * divU, k) - fvm::Sp(alpha() * rho() * this->epsilonByk(F1, tgradU()), k) + alpha() * rho() * this->betaStar_ * this->omegaInf_ * this->kInf_ + this->kSource() + fvOptions(alpha, rho, k));

            tgradU.clear();

            kEqn.ref().relax();
            fvOptions.constrain(kEqn.ref());
            solve(kEqn);
            fvOptions.correct(k);
            bound(k, this->kMin_);

            correctNut(S2);
        }
    } // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
