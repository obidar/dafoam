/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v3

    Description:
    Child class for the SpalartAllmarasFv3 model with a betaFieldInversion field multiplying 
    to the production term. This betaFieldInversion term can then be trained to improve the
    SA turbulence model.

    This file is modified from OpenFOAM's source code
    src/TurbulenceModels/turbulenceModels/RAS/SpalartAllmaras/SpalartAllmaras.H

    OpenFOAM: The Open Source CFD Toolbox

    Copyright (C): 2011-2016 OpenFOAM Foundation

    OpenFOAM License:

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

#ifndef DASpalartAllmarasFv3FIML_H
#define DASpalartAllmarasFv3FIML_H

#include "DATurbulenceModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class DASpalartAllmarasFv3FIML Declaration
\*---------------------------------------------------------------------------*/

class DASpalartAllmarasFv3FIML
    : public DATurbulenceModel
{

protected:
    /// \name SA parameters
    //@{
    dimensionedScalar sigmaNut_;
    dimensionedScalar kappa_;
    dimensionedScalar Cb1_;
    dimensionedScalar Cb2_;
    dimensionedScalar Cw1_;
    dimensionedScalar Cw2_;
    dimensionedScalar Cw3_;
    dimensionedScalar Cv1_;
    dimensionedScalar Cv2_;

    /// inputs and outputs for the beta calculation
    scalar* inputs_ = new scalar[mesh_.nCells() * 7];
    scalar* outputs_ = new scalar[mesh_.nCells()];

#ifdef CODI_AD_FORWARD
    double* inputsDouble_ = new double[mesh_.nCells() * 9];
    double* outputsDouble_ = new double[mesh_.nCells()];
#endif

    /// calculate the beta field using the trained model
    void calcBetaField();

    //@}

    /// \name SA functions
    //@{
    tmp<volScalarField> chi() const;

    tmp<volScalarField> fv1(const volScalarField& chi) const;

    tmp<volScalarField> fv2(
        const volScalarField& chi,
        const volScalarField& fv1) const;

    tmp<volScalarField> fv3(
        const volScalarField& chi,
        const volScalarField& fv1) const;

    tmp<volScalarField> fw(const volScalarField& Stilda) const;

    //@}

    /// \name Augmented variables for residual computation
    //@{
    volScalarField& nuTilda_;
    volScalarField nuTildaRes_;
    volScalarField nuTildaResPartDeriv_;
    //@}

    /// A beta field multiplying to the production term
    volScalarField betaFieldInversion_; 

    volScalarField betaFieldInversionML_; 

    /// 3D wall distance
    const volScalarField& y_;

    /// whether to solve for turb states
    label solveTurbState_ = 0;

    /// time step interval to print residual
    label printInterval_;

public:
    TypeName("SpalartAllmarasFv3FIML");
    // Constructors

    //- Construct from components
    DASpalartAllmarasFv3FIML(
        const word modelType,
        const fvMesh& mesh,
        const DAOption& daOption);

    //- Destructor
    virtual ~DASpalartAllmarasFv3FIML()
    {
        delete[] inputs_;
        delete[] outputs_;
#ifdef CODI_AD_FORWARD
        delete[] inputsDouble_;
        delete[] outputsDouble_;
#endif
    }

    // Member functions

    /// Return the effective diffusivity for nuTilda
    tmp<volScalarField> DnuTildaEff() const;

    /// update the turbulence state for DAStateInfo::regStates_
    virtual void correctModelStates(wordList& modelStates) const;

    /// update nut based on other turbulence variables and update the BCs
    virtual void correctNut();

    /// update turbulence variable boundary values
    virtual void correctBoundaryConditions();

    /// update any intermediate variables that are dependent on state variables and are used in calcResiduals
    virtual void updateIntermediateVariables();

    /// update the original variable connectivity for the adjoint state residuals in stateCon
    virtual void correctStateResidualModelCon(List<List<word>>& stateCon) const;

    /// add the model residual connectivity to stateCon
    virtual void addModelResidualCon(HashTable<List<List<word>>>& allCon) const;

    /// compute the turbulence residuals
    virtual void calcResiduals(const dictionary& options);

    /// solve the residual equations and update the state
    virtual void correct();

    /// return the value of the production term from the turbulence model 
    virtual void getTurbProdTerm(scalarList& prodTerm) const;

#ifdef CODI_AD_REVERSE

    /// these two functions are for AD external functions
    static void betaCompute(
        const double* x,
        size_t n,
        double* y,
        size_t m,
        codi::ExternalFunctionUserData* d)
    {
        DAUtility::pyCalcBetaInterface(x, n, y, m, DAUtility::pyCalcBeta);
    }

    static void betaJacVecProd(
        const double* x,
        double* x_b,
        size_t n,
        const double* y,
        const double* y_b,
        size_t m,
        codi::ExternalFunctionUserData* d)
    {
        DAUtility::pyCalcBetaJacVecProdInterface(x, x_b, n, y, y_b, m, DAUtility::pyCalcBetaJacVecProd);
    }
#endifs
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
