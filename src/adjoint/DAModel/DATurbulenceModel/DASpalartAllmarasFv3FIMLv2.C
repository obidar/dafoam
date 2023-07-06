/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v3

    This file is modified from OpenFOAM's source code
    src/TurbulenceModels/turbulenceModels/RAS/SpalartAllmaras/SpalartAllmaras.C

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

#include "DASpalartAllmarasFv3FIMLv2.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(DASpalartAllmarasFv3FIMLv2, 0);
addToRunTimeSelectionTable(DATurbulenceModel, DASpalartAllmarasFv3FIMLv2, dictionary);
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

DASpalartAllmarasFv3FIMLv2::DASpalartAllmarasFv3FIMLv2(
    const word modelType,
    const fvMesh& mesh,
    const DAOption& daOption)
    : DATurbulenceModel(modelType, mesh, daOption),
      // SA parameters
      sigmaNut_(dimensioned<scalar>::lookupOrAddToDict(
          "sigmaNut",
          this->coeffDict_,
          0.66666)),
      kappa_(dimensioned<scalar>::lookupOrAddToDict(
          "kappa",
          this->coeffDict_,
          0.41)),

      Cb1_(dimensioned<scalar>::lookupOrAddToDict(
          "Cb1",
          this->coeffDict_,
          0.1355)),
      Cb2_(dimensioned<scalar>::lookupOrAddToDict(
          "Cb2",
          this->coeffDict_,
          0.622)),
      Cw1_(Cb1_ / sqr(kappa_) + (1.0 + Cb2_) / sigmaNut_),
      Cw2_(dimensioned<scalar>::lookupOrAddToDict(
          "Cw2",
          this->coeffDict_,
          0.3)),
      Cw3_(dimensioned<scalar>::lookupOrAddToDict(
          "Cw3",
          this->coeffDict_,
          2.0)),
      Cv1_(dimensioned<scalar>::lookupOrAddToDict(
          "Cv1",
          this->coeffDict_,
          7.1)),
      Cv2_(dimensioned<scalar>::lookupOrAddToDict(
          "Cv2",
          this->coeffDict_,
          5.0)),
      // Augmented variables
      nuTilda_(const_cast<volScalarField&>(
          mesh.thisDb().lookupObject<volScalarField>("nuTilda"))),
      nuTildaRes_(
          IOobject(
              "nuTildaRes",
              mesh.time().timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::NO_WRITE),
          mesh,
#ifdef CompressibleFlow
          dimensionedScalar("nuTildaRes", dimensionSet(1, -1, -2, 0, 0, 0, 0), 0.0),
#endif
#ifdef IncompressibleFlow
          dimensionedScalar("nuTildaRes", dimensionSet(0, 2, -2, 0, 0, 0, 0), 0.0),
#endif
          zeroGradientFvPatchField<scalar>::typeName),
      nuTildaResPartDeriv_(
          IOobject(
              "nuTildaResPartDeriv",
              mesh.time().timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::NO_WRITE),
          nuTildaRes_),
      betaFieldInversion_(
          IOobject(
              "betaFieldInversion",
              mesh.time().timeName(),
              mesh_,
              IOobject::READ_IF_PRESENT,
              IOobject::AUTO_WRITE),
          mesh_,
          dimensionedScalar("betaFieldInversion", dimensionSet(0, 0, 0, 0, 0, 0, 0), 1.0),
          zeroGradientFvPatchScalarField::typeName),
      betaFieldInversionML_(
          IOobject(
              "betaFieldInversionML",
              mesh.time().timeName(),
              mesh_,
              IOobject::NO_READ,
              IOobject::AUTO_WRITE),
          mesh_,
          dimensionedScalar("betaFieldInversionML", dimensionSet(0, 0, 0, 0, 0, 0, 0), 1.0),
          zeroGradientFvPatchScalarField::typeName),
      QCriterion_(
          IOobject(
              "QCriterion",
              mesh.time().timeName(),
              mesh_,
              IOobject::NO_READ,
              IOobject::AUTO_WRITE),
          mesh_,
          dimensionedScalar("QCriterion", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.0),
          zeroGradientFvPatchScalarField::typeName),
    p_(const_cast<volScalarField&>(
          mesh_.thisDb().lookupObject<volScalarField>("p"))),
    pGradAlongStream_(
          IOobject(
              "pGradAlongStream",
              mesh.time().timeName(),
              mesh_,
              IOobject::NO_READ,
              IOobject::AUTO_WRITE),
          mesh_,
          dimensionedScalar("pGradAlongStream", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.0),
          zeroGradientFvPatchScalarField::typeName),
    pressureStress_(
          IOobject(
              "pressureStress",
              mesh.time().timeName(),
              mesh_,
              IOobject::NO_READ,
              IOobject::AUTO_WRITE),
          mesh_,
          dimensionedScalar("pressureStress", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.0),
          zeroGradientFvPatchScalarField::typeName),
    curvature_(
          IOobject(
              "curvature",
              mesh.time().timeName(),
              mesh_,
              IOobject::NO_READ,
              IOobject::AUTO_WRITE),
          mesh_,
          dimensionedScalar("curvature", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.0),
          zeroGradientFvPatchScalarField::typeName),
    UGradMisalignment_(
          IOobject(
              "UGradMisalignment",
              mesh.time().timeName(),
              mesh_,
              IOobject::NO_READ,
              IOobject::AUTO_WRITE),
          mesh_,
          dimensionedScalar("UGradMisalignment", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.0),
          zeroGradientFvPatchScalarField::typeName),
    viscosityRatio_(
          IOobject(
              "viscosityRatio",
              mesh.time().timeName(),
              mesh_,
              IOobject::NO_READ,
              IOobject::AUTO_WRITE),
          mesh_,
          dimensionedScalar("viscosityRatio", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.0),
          zeroGradientFvPatchScalarField::typeName),
    wallInfluence_(
          IOobject(
              "wallInfluence",
              mesh.time().timeName(),
              mesh_,
              IOobject::NO_READ,
              IOobject::AUTO_WRITE),
          mesh_,
          dimensionedScalar("wallInfluence", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.0),
          zeroGradientFvPatchScalarField::typeName),
    ratioProductionToDiffusion_(
          IOobject(
              "ratioProductionToDiffusion",
              mesh.time().timeName(),
              mesh_,
              IOobject::NO_READ,
              IOobject::AUTO_WRITE),
          mesh_,
          dimensionedScalar("ratioProductionToDiffusion", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.0),
          zeroGradientFvPatchScalarField::typeName),
    ratioDestructionToDiffusion_(
          IOobject(
              "ratioDestructionToDiffusion",
              mesh.time().timeName(),
              mesh_,
              IOobject::NO_READ,
              IOobject::AUTO_WRITE),
          mesh_,
          dimensionedScalar("ratioDestructionToDiffusion", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.0),
          zeroGradientFvPatchScalarField::typeName),
     refViscosity_(transportProperties.lookup("nu")), 
      y_(mesh.thisDb().lookupObject<volScalarField>("yWall"))
{

    // initialize printInterval_ we need to check whether it is a steady state
    // or unsteady primal solver
    IOdictionary fvSchemes(
        IOobject(
            "fvSchemes",
            mesh.time().system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false));
    word ddtScheme = word(fvSchemes.subDict("ddtSchemes").lookup("default"));
    if (ddtScheme == "steadyState")
    {
        printInterval_ =
            daOption.getAllOptions().lookupOrDefault<label>("printInterval", 100);
    }
    else
    {
        printInterval_ =
            daOption.getAllOptions().lookupOrDefault<label>("printIntervalUnsteady", 500);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// SA member functions. these functions are copied from
tmp<volScalarField> DASpalartAllmarasFv3FIMLv2::chi() const
{
    return nuTilda_ / this->nu();
}

tmp<volScalarField> DASpalartAllmarasFv3FIMLv2::fv1(
    const volScalarField& chi) const
{
    const volScalarField chi3(pow3(chi));
    return chi3 / (chi3 + pow3(Cv1_));
}

tmp<volScalarField> DASpalartAllmarasFv3FIMLv2::fv2(
    const volScalarField& chi,
    const volScalarField& fv1) const
{
    return 1.0 / pow3(scalar(1) + chi / Cv2_);
}

tmp<volScalarField> DASpalartAllmarasFv3FIMLv2::fv3(
    const volScalarField& chi,
    const volScalarField& fv1) const
{

    const volScalarField chiByCv2((1 / Cv2_) * chi);

    return (scalar(1) + chi * fv1)
        * (1 / Cv2_)
        * (3 * (scalar(1) + chiByCv2) + sqr(chiByCv2))
        / pow3(scalar(1) + chiByCv2);
}

tmp<volScalarField> DASpalartAllmarasFv3FIMLv2::fw(
    const volScalarField& Stilda) const
{
    volScalarField r(
        min(
            nuTilda_
                / (max(
                       Stilda,
                       dimensionedScalar("SMALL", Stilda.dimensions(), SMALL))
                   * sqr(kappa_ * y_)),
            scalar(10.0)));
    r.boundaryFieldRef() == 0.0;

    const volScalarField g(r + Cw2_ * (pow6(r) - r));

    return g * pow((1.0 + pow6(Cw3_)) / (pow6(g) + pow6(Cw3_)), 1.0 / 6.0);
}

tmp<volScalarField> DASpalartAllmarasFv3FIMLv2::DnuTildaEff() const
{
    return tmp<volScalarField>(
        new volScalarField("DnuTildaEff", (nuTilda_ + this->nu()) / sigmaNut_));
}

// Augmented functions
void DASpalartAllmarasFv3FIMLv2::correctModelStates(wordList& modelStates) const
{
    /*
    Description:
        Update the name in modelStates based on the selected physical model at runtime

    Example:
        In DAStateInfo, if the modelStates reads:
        
        modelStates = {"nut"}
        
        then for the SA model, calling correctModelStates(modelStates) will give:
    
        modelStates={"nuTilda"}
        
        while calling correctModelStates(modelStates) for the SST model will give 
        
        modelStates={"k","omega"}
        
        We don't udpate the names for the radiation model becasue users are 
        supposed to set modelStates={"G"}
    */

    // replace nut with nuTilda
    forAll(modelStates, idxI)
    {
        word stateName = modelStates[idxI];
        if (stateName == "nut")
        {
            modelStates[idxI] = "nuTilda";
        }
    }
}

void DASpalartAllmarasFv3FIMLv2::correctNut()
{
    /*
    Description:
        Update nut based on other turbulence variables and update the BCs
        Also update alphat if is present
    */

    const volScalarField chi(this->chi());
    const volScalarField fv1(this->fv1(chi));
    nut_ = nuTilda_ * fv1;

    nut_.correctBoundaryConditions();

    // this is basically BasicTurbulenceModel::correctNut();
    this->correctAlphat();

    return;
}

void DASpalartAllmarasFv3FIMLv2::correctBoundaryConditions()
{
    /*
    Description:
        Update turbulence variable boundary values
    */

    // correct the BCs for the perturbed fields
    nuTilda_.correctBoundaryConditions();
}

void DASpalartAllmarasFv3FIMLv2::updateIntermediateVariables()
{
    /*
    Description:
        Update nut based on nuTilda. Note: we need to update nut and its BC since we 
        may have perturbed other turbulence vars that affect the nut values
    */

    this->correctNut();
}

void DASpalartAllmarasFv3FIMLv2::correctStateResidualModelCon(List<List<word>>& stateCon) const
{
    /*
    Description:
        Update the original variable connectivity for the adjoint state 
        residuals in stateCon. Basically, we modify/add state variables based on the
        original model variables defined in stateCon.

    Input:
    
        stateResCon: the connectivity levels for a state residual, defined in Foam::DAJacCon

    Example:
        If stateCon reads:
        stateCon=
        {
            {"U", "p", "nut"},
            {"p"}
        }
    
        For the SA turbulence model, calling this function for will get a new stateCon
        stateCon=
        {
            {"U", "p", "nuTilda"},
            {"p"}
        }
    
        For the SST turbulence model, calling this function will give
        stateCon=
        {
            {"U", "p", "k", "omega"},
            {"p", "U"}
        }
        ***NOTE***: we add a extra level of U connectivity because nut is 
        related to grad(U), k, and omega in SST!
    */

    forAll(stateCon, idxI)
    {
        forAll(stateCon[idxI], idxJ)
        {
            word conStateName = stateCon[idxI][idxJ];
            if (conStateName == "nut")
            {
                stateCon[idxI][idxJ] = "nuTilda";
            }
        }
    }
}

void DASpalartAllmarasFv3FIMLv2::addModelResidualCon(HashTable<List<List<word>>>& allCon) const
{
    /*
    Description:
        Add the connectivity levels for all physical model residuals to allCon

    Input:
        allCon: the connectivity levels for all state residual, defined in DAJacCon

    Example:
        If stateCon reads:
        allCon=
        {
            "URes":
            {
               {"U", "p", "nut"},
               {"p"}
            }
        }
    
        For the SA turbulence model, calling this function for will get a new stateCon,
        something like this:
        allCon=
        {
            "URes":
            {
               {"U", "p", "nuTilda"},
               {"p"}
            },
            "nuTildaRes": 
            {
                {"U", "phi", "nuTilda"},
                {"U"}
            }
        }

    */

    word pName;

    if (mesh_.thisDb().foundObject<volScalarField>("p"))
    {
        pName = "p";
    }
    else if (mesh_.thisDb().foundObject<volScalarField>("p_rgh"))
    {
        pName = "p_rgh";
    }
    else
    {
        FatalErrorIn(
            "Neither p nor p_rgh was found in mesh.thisDb()!"
            "addModelResidualCon failed to setup turbulence residuals!")
            << exit(FatalError);
    }

    // NOTE: for compressible flow, it depends on rho so we need to add T and p
#ifdef IncompressibleFlow
    allCon.set(
        "nuTildaRes",
        {
            {"U", "nuTilda", "phi"}, // lv0
            {"U", "nuTilda"}, // lv1
            {"nuTilda"} // lv2
        });
#endif

#ifdef CompressibleFlow
    allCon.set(
        "nuTildaRes",
        {
            {"U", "T", pName, "nuTilda", "phi"}, // lv0
            {"U", "T", pName, "nuTilda"}, // lv1
            {"T", pName, "nuTilda"} // lv2
        });
#endif
}

void DASpalartAllmarasFv3FIMLv2::correct()
{
    /*
    Descroption:
        Solve the residual equations and update the state. This function will be called 
        by the DASolver. It is needed because we want to control the output frequency
        of the residual convergence every 100 steps. If using the correct from turbulence
        it will output residual every step which will be too much of information.
    */

    // We set the flag solveTurbState_ to 1 such that in the calcResiduals function
    // we will solve and update nuTilda
    solveTurbState_ = 1;
    dictionary dummyOptions;
    this->calcResiduals(dummyOptions);
    // after it, we reset solveTurbState_ = 0 such that calcResiduals will not
    // update nuTilda when calling from the adjoint class, i.e., solveAdjoint from DASolver.
    solveTurbState_ = 0;
}

void DASpalartAllmarasFv3FIMLv2::calcBetaField()
{

    // COMPUTE MACHINE LEARNING FEATURES
    volTensorField UGrad(fvc::grad(U_));
    volTensorField Omega("Omega",skew(UGrad));
    volScalarField magOmegaSqr(magSqr(Omega));
    volSymmTensorField S("S",symm(UGrad));
    volScalarField magS(mag(S));
    volScalarField magSSqr(magSqr(S));
    QCriterion_ = (magOmegaSqr -  magSSqr)/(magOmegaSqr + magSSqr); 

    volVectorField pGrad("gradP",fvc::grad(p_));
    volScalarField pG_denominator (mag(U_) * mag(pGrad) + mag(U_ & pGrad));
    pGradAlongStream_ =  (U_ & pGrad) / Foam::max(pG_denominator, dimensionedScalar("minpG",dimensionSet(0,2,-3,0,0,0,0),SMALL)); 
    volVectorField diagUGrad
    (IOobject("diagUGrad",runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh,
        dimensionedVector("diagUGrad", dimensionSet(0,0,0,0,0,0,0),  Foam::vector(0,0,0)),
        zeroGradientFvPatchScalarField::typeName
    );
    forAll(mesh_.cells(), cI)
    {
        diagUGrad[cI].component(0) = UGrad[cI].xx(); 
        diagUGrad[cI].component(1) = UGrad[cI].yy(); 
        diagUGrad[cI].component(2) =  UGrad[cI].zz(); 
        pressureStress_[cI] = mag(pGrad[cI]) / (mag(pGrad[cI]) + mag(3.0*cmptAv(U_[cI] & diagUGrad[cI])));
    }

    forAll(mesh_.cells(), cI)
    {
        curvature_[cI] = mag(U_[cI] & UGrad[cI]) 
                        /
                        (
                            mag(U_[cI] &  U_[cI])
                          + mag(U_[cI] & UGrad[cI])
                        );
    }

    forAll(mesh_.cells(), cI)
    {
        UGradMisalignment_[cI] = mag(U_[cI] & UGrad[cI] & U_[cI]) 
                        / 
                        ( 
                            mag(U_[cI]) * mag(UGrad[cI] & U_[cI])
                            + 
                            mag(U_[cI] & UGrad[cI] & U_[cI]) 
                        );
    }

    forAll(mesh_.cells(), cI)
    {
        viscosityRatio_[cI] = nut_[cI] / (nut_[cI] + 100 * refViscosity_.value()); 
    }   

    volScalarField d(wallDist::New(mesh).y_());
    volScalarField chi(nuTilda_ / refViscosity_);
    volScalarField fv1(pow3(chi) / (pow3(chi) + pow3(Cv1_))); 
    volScalarField fv2(1 / pow3(1 + chi / Cv2_)); 
    volScalarField fv3(((1 + chi * fv1_) * (1 - fv2_)) /chi);
    volScalarField STilda(fv3_ * Foam::sqrt(2.0) * mag(skew(fvc::grad(U))) + (fv2 * nuTilda_ / (Foam::sqr(kappa_ * d))));
    
    volScalarField r(
    Foam::min
    (
            nuTilda_
           /(
               Foam::max
               (
                   STilda,
                   dimensionedScalar("SMALL", STilda.dimensions(), SMALL)
               )
              *sqr(kappa_*d)
            ),
            scalar(10.0)
        )
    ); 
    
    //volScalarField r(Foam::min((nuTilda/(Foam::sqr(kappa * d) * STilda),10));
    forAll(mesh_.cells(), cI)
    {
        wallInfluence_[cI] = 1 / (1 + r[cI]); 
    }

    volScalarField numTemp(Cb1_ * STilda * nuTilda_);
    volScalarField magNumTemp(mag(numTemp)); 
    volScalarField epsilonTemp(Cb2_/sigmaNut * magSqr(fvc::grad(nuTilda_))); 
    forAll(mesh_.cells(), cI)
    {
        ratioProductionToDiffusion[cI] = numTemp[cI] / (magNumTemp[cI] + epsilonTemp[cI]); 
    }

    volScalarField g(r + Cw2_*(pow6(r) - r));
    volScalarField fw(g*pow((1.0 + pow6(Cw3_))/(pow6(g) + pow6(Cw3_)), 1.0/6.0)); 
    volScalarField num2Temp(Cw1*_fw*Foam::sqr(nuTilda_/d)); 
    forAll(mesh_.cells(), cI)
    {
        ratioDestructionToDiffusion_[cI] = num2Temp[cI] / (mag(num2Temp[cI]) + epsilonTemp[cI]);
    }

    forAll(mesh_.cells(), cI)
    {
        inputs_[cI * 9 + 0] = QCriterion_[cI];
        inputs_[cI * 9 + 1] = pGradAlongStream_[cI];
        inputs_[cI * 9 + 2] = pressureStress_[cI];
        inputs_[cI * 9 + 3] = curvature_[cI];
        inputs_[cI * 9 + 4] = UGradMisalignment_[cI];
        inputs_[cI * 9 + 5] = viscosityRatio_[cI];
        inputs_[cI * 9 + 6] = wallInfluence_[cI];
        inputs_[cI * 9 + 7] = ratioProductionToDiffusion_[cI];
        inputs_[cI * 9 + 8] = ratioDestructionToDiffusion_[cI];
    }

    // NOTE: forward mode not supported..
#if defined(CODI_AD_REVERSE)

    // we need to use the external function helper from CoDiPack to propagate the AD

    codi::ExternalFunctionHelper<codi::RealReverse> externalFunc;
    for (label i = 0; i < mesh_.nCells() * 7; i++)
    {
        externalFunc.addInput(inputs_[i]);
    }

    for (label i = 0; i < mesh_.nCells(); i++)
    {
        externalFunc.addOutput(outputs_[i]);
    }

    externalFunc.callPrimalFunc(DASpalartAllmarasFv3FIMLv2::betaCompute);

    codi::RealReverse::Tape& tape = codi::RealReverse::getTape();

    if (tape.isActive())
    {
        externalFunc.addToTape(DASpalartAllmarasFv3FIMLv2::betaJacVecProd);
    }

    forAll(betaFieldInversionML_, cellI)
    {
        betaFieldInversionML_[cellI] = outputs_[cellI];
    }

#elif defined(CODI_AD_FORWARD)

    for (label i = 0; i < n; i++)
    {
        inputsDouble_[i] = inputs_[i].value();
    }

    for (label i = 0; i < m; i++)
    {
        outputsDouble_[i] = outputs_[i].value();
    }

    // python callback function
    DAUtility::pyCalcBetaInterface(inputsDouble_, n, outputsDouble_, m, DAUtility::pyCalcBeta);

    forAll(betaFieldInversionML_, cellI)
    {
        betaFieldInversionML_[cellI] = outputsDouble_[cellI];
    }

#else

    // python callback function
    DAUtility::pyCalcBetaInterface(inputs_, n, outputs_, m, DAUtility::pyCalcBeta);

    forAll(betaFieldInversionML_, cellI)
    {
        betaFieldInversionML_[cellI] = outputs_[cellI];
    }

#endif
}

void DASpalartAllmarasFv3FIMLv2::calcResiduals(const dictionary& options)
{
    /*
    Descroption:
        If solveTurbState_ == 1, this function solve and update nuTilda, and 
        is the same as calling turbulence.correct(). If solveTurbState_ == 0,
        this function compute residuals for turbulence variables, e.g., nuTildaRes_

    Input:
        options.isPC: 1 means computing residuals for preconditioner matrix.
        This essentially use the first order scheme for div(phi,nuTilda)

        p_, U_, phi_, etc: State variables in OpenFOAM
    
    Output:
        nuTildaRes_: If solveTurbState_ == 0, update the residual field variable

        nuTilda_: If solveTurbState_ == 1, update nuTilda
    */

    // Copy and modify based on the "correct" function

    label printToScreen = this->isPrintTime(mesh_.time(), printInterval_);

    word divNuTildaScheme = "div(phi,nuTilda)";

    label isPC = 0;

    if (!solveTurbState_)
    {
        isPC = options.getLabel("isPC");

        if (isPC)
        {
            divNuTildaScheme = "div(pc)";
        }
    }

    this->calcBetaField();

    const volScalarField chi(this->chi());
    const volScalarField fv1(this->fv1(chi));

    const volScalarField Stilda(
        this->fv3(chi, fv1) * ::sqrt(2.0) * mag(skew(fvc::grad(U_)))
        + this->fv2(chi, fv1) * nuTilda_ / sqr(kappa_ * y_));

    tmp<fvScalarMatrix> nuTildaEqn(
        fvm::ddt(phase_, rho_, nuTilda_)
            + fvm::div(phaseRhoPhi_, nuTilda_, divNuTildaScheme)
            - fvm::laplacian(phase_ * rho_ * DnuTildaEff(), nuTilda_)
            - Cb2_ / sigmaNut_ * phase_ * rho_ * magSqr(fvc::grad(nuTilda_))
        == Cb1_ * phase_ * rho_ * Stilda * nuTilda_ * betaFieldInversionML_
            - fvm::Sp(Cw1_ * phase_ * rho_ * fw(Stilda) * nuTilda_ / sqr(y_), nuTilda_));

    nuTildaEqn.ref().relax();

    if (solveTurbState_)
    {

        // get the solver performance info such as initial
        // and final residuals
        SolverPerformance<scalar> solverNuTilda = solve(nuTildaEqn);
        if (printToScreen)
        {
            Info << "nuTilda Initial residual: " << solverNuTilda.initialResidual() << endl
                 << "          Final residual: " << solverNuTilda.finalResidual() << endl;
        }

        DAUtility::boundVar(allOptions_, nuTilda_, printToScreen);
        nuTilda_.correctBoundaryConditions();

        // ***************** NOTE*****************
        // In the original SA, it is correctNut(fv1) and fv1 is not
        // updated based on the latest nuTilda. We use correctNut which
        // recompute fv1 with the latest nuTilda
        this->correctNut();
    }
    else
    {
        // calculate residuals
        nuTildaRes_ = nuTildaEqn.ref() & nuTilda_;
        // need to normalize residuals
        normalizeResiduals(nuTildaRes);
    }

    return;
}

void DASpalartAllmarasFv3FIMLv2::getTurbProdTerm(scalarList& prodTerm) const
{
    /*
    Description:
        Return the value of the production term from the Spalart Allmaras model 
    */

    const volScalarField chi(this->chi());
    const volScalarField fv1(this->fv1(chi));

    const volScalarField Stilda(
        this->fv3(chi, fv1) * ::sqrt(2.0) * mag(skew(fvc::grad(U_)))
        + this->fv2(chi, fv1) * nuTilda_ / sqr(kappa_ * y_));
    
    volScalarField SAProdTerm = Cb1_ * phase_ * rho_ * Stilda * nuTilda_;

    forAll(SAProdTerm, cellI)
    {
        prodTerm[cellI] = SAProdTerm[cellI];
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
