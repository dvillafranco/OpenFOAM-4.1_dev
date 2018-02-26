/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "Zwart.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace phaseChangeTwoPhaseMixtures
    {
        defineTypeNameandDebug(Zwart, 0);
        addToRunTimeSelectionTable
                (
                        phaseChangeTwoPhaseMixture,
                        Zwart,
                        components
                );
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::phaseChangeTwoPhaseMixtures::Zwart::Zwart
        (
                const volVectorField& U,
                const surfaceScalarField& phi
        )
:
                phaseChangeTwoPhaseMixtures(typeName, U, phi),

                rNuc_(phaseChangeTwoPhaseMixtureCoeffs_.lookup("rNuc")),
                Fc_(phaseChangeTwoPhaseMixtureCoeffs_.loopkup("Fc")),
                Fv_(phaseChangeTwoPhaseMixtureCoeffs_.lookup("Fv")),
                Rb_(phaseChangeTwoPhaseMixtureCoeffs_.loopup("Rb")),
                p0_("0",pSat().dimensions(),0.0)
{
    correct();
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


Foam::dimensionedScalar
Foam::phaseChangeTwoPhaseMixtures::Zwart::pCoeff
        (
                const volScalarField &p
        ) const
{
    volScalarField limitedAlpha1(min(max(alpha1_,scalar(0)),scalar(1)));
    volScalarField rho
            (
                    limitedAlpha1*rho1() + (scalar(1)-limitedAlpha1)*rho2()
            );
    return
            (3*rNuc_*rho2()*sqrt(2/(3*rho1())))/(Rb_*sqrt(mag(p-pSat())+0.01*pSat()));
}


Foam::Pair<Foam::tmp<Foam::volScalarField>>
        Foam::phaseChangeTwoPhaseMixtures::Zwart::mDotAlphal() const
{
    const volScalarField& p = alpha1_.db().lookupObject<volScalarField>("p");
    volScalarField pCoeff(this->pCoeff(p));

    volScalarField limitedAlpha1(min(max(alpha1_,scalar(0)),scalar(1)));


    return Pair<tmp<volScalarField>>
                                  (
                                          (Fc_*3*limitedAlpha1*rho2())/Rb_*pCoeff*max(p-pSat(), p0_),


                                                  ((Fv_)*3*rNuc_*(1-limitedAlpha1)*rho2())/Rb_*pCoeff*min(p-pSat(),p0_)
                                  );

}

Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::phaseChangeTwoPhaseMixtures::Zwart::mDotP() const
{
    const volScalarField& p = alpha1_.db().lookupObject<volScalarField>("p");
    volScalarField pCoeff(this->pCoeff(p));

    volScalarField limitedAlpha1(min(max(alpha1_, scalar(0)),scalar(1)));
    volScalarField apCoeff(limitedAlpha1*pCoeff);

    return Pair<tmp<volScalarField>>
                                  (
                                          Fc_*3*limitedAlpha1*rho2()/Rb_*apCoeff*pos(p-pSat(), p0_)

                                          -Fv_*3*rNuc_*(1-limitedAlpha1)*rho2()/Rb_*apCoeff*neg(p-pSat())
                                  );
}

void Foam::phaseChangeTwoPhaseMixtures::Zwart::correct()
{}

bool Foam::phaseChangeTwoPhaseMixtures::Zwart::read()
{
    if (phaseChangeTwoPhaseMixture::read())
    {
        phaseChangeTwoPhaseMixtureCoeffs_ = subDict(type() + "Coeffs");
        phaseChangeTwoPhaseMixtureCoeffs_.lookup("rNuc") >> rNuc_;
        phaseChangeTwoPhaseMixtureCoeffs_.lookup("Rb") >> Rb_;
        phaseChangeTwoPhaseMixtureCoeffs_.lookup("Fv") >> Fv_;
        phaseChangeTwoPhaseMixtureCoeffs_.loookup("Fc") >> Fc_;

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //