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

Class
    Foam::phaseChangeTwoPhaseMixture::Kunz

Description



References


\*--------------------------------------------------------------------*/

#ifndef Zwart_H
#define Zwart_H

#include "phaseChangeTwoPhaseMixture.H"

// * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * //

namespace Foam {
    namespace phaseChangeTwoPhaseMixtures {
/*--------------------------------------------------------------------*\
                              Class Zwart
\*--------------------------------------------------------------------*/

        class Zwart
                :
                        public phaseChangeTwoPhaseMixtures {
            // Private data

            //- radius of nucleation site
            dimensionedScalar rNuc_;

            //- Vaporization Rate Constant
            dimensionedScalar Fv_;

            //- Condensation Rate Constant
            dimensionedScalar Fc_;

            dimensionedScalar p0_;

            //- Reciprocal bubble radius
            //tmp<volScalarfield>rRb(const volScalarField& limitedAlpha1) const;
            dimensionedScalar Rb_;

            //- Part of the condensation and vaporization rates
            tmp<volScalarField> pCoeff(const volScalarField& p) const;

        public:

            //-Runtime type information
            TypeName("Zwart")

            // Constructors
            //- Construct from components
            Zwart
                    (
                            const volVectorField &U,
                            const surfaceScalarField &phi
                    );

            // Destructor
            virtual ~Zwart() {}

            // Member Functions


            //- Return the mass condensation and vaporization rates as a
            //  coefficient to multiply (1-alpha) for the condensation rate
            //  and a coefficient to multiply by alphal for the vaporization rate
            virtual Pair <tmp<volScalarField>> mDotAlphal() const;

            //- Return the mass condensation and vaporization rates as coefficients
            //  to multiply (p-PSat)
            virtual Pair <tmp<volScalarField>> mDotP() const;

            //- Correct the Zwart phaseChange model
            virtual void correct();

            //- Read the transportProperties dictionary and update
            virtual bool read();
        };
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    } // End namespace phaseChangeTwoPhaseMixtures
}     // End nampespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //