/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v3

\*---------------------------------------------------------------------------*/

#include "SpalartAllmarasFv3FIMLv2.H"
#include "fvOptions.H"
#include "bound.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void SpalartAllmarasFv3FIMLv2<BasicTurbulenceModel>::correctNut()
{
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
SpalartAllmarasFv3FIMLv2<BasicTurbulenceModel>::SpalartAllmarasFv3FIMLv2(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type)
    : eddyViscosity<RASModel<BasicTurbulenceModel>>(
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName),
      nuTilda_(
          IOobject(
              "nuTilda",
              this->runTime_.timeName(),
              this->mesh_,
              IOobject::MUST_READ,
              IOobject::AUTO_WRITE),
          this->mesh_),
      /*betaFieldInversion_(
          IOobject(
              "betaFieldInversion",
              this->runTime_.timeName(),
              this->mesh_,
              IOobject::MUST_READ,
              IOobject::AUTO_WRITE),
          this->mesh_),
      betaFieldInversionML_(
          IOobject(
              "betaFieldInversionML",
              this->runTime_.timeName(),
              this->mesh_,
              IOobject::NO_READ,
              IOobject::AUTO_WRITE),
          this->mesh_),*/
      y_(wallDist::New(this->mesh_).y())
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool SpalartAllmarasFv3FIMLv2<BasicTurbulenceModel>::read()
{
    return true;
}

template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmarasFv3FIMLv2<BasicTurbulenceModel>::k() const
{
    return this->nut_;
}

template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmarasFv3FIMLv2<BasicTurbulenceModel>::epsilon() const
{

    return this->nut_;
}

template<class BasicTurbulenceModel>
void SpalartAllmarasFv3FIMLv2<BasicTurbulenceModel>::correct()
{
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
