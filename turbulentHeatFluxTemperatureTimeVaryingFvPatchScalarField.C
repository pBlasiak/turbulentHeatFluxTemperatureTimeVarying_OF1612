/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2016 OpenFOAM Foundation
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

#include "turbulentHeatFluxTemperatureTimeVaryingFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "turbulenceModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    // declare specialization within 'Foam' namespace
    template<>
    const char* NamedEnum
    <
        Foam::incompressible::
        turbulentHeatFluxTemperatureTimeVaryingFvPatchScalarField::heatSourceType,
        2
    >::names[] =
    {
        "power",
        "flux"
    };
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


namespace Foam
{

namespace incompressible
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const NamedEnum
<
    turbulentHeatFluxTemperatureTimeVaryingFvPatchScalarField::heatSourceType,
    2
> turbulentHeatFluxTemperatureTimeVaryingFvPatchScalarField::heatSourceTypeNames_;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

turbulentHeatFluxTemperatureTimeVaryingFvPatchScalarField::
turbulentHeatFluxTemperatureTimeVaryingFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    uniformFixedGradientFvPatchField<scalar>(p, iF),
    heatSource_(hsPower),
    //q_(p.size(), 0.0),
    alphaEffName_("undefinedAlphaEff"),
	q_()
{

	Info<< "JESTEM3" << endl;
}

turbulentHeatFluxTemperatureTimeVaryingFvPatchScalarField::turbulentHeatFluxTemperatureTimeVaryingFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    uniformFixedGradientFvPatchField<scalar>(p, iF),
    heatSource_(heatSourceTypeNames_.read(dict.lookup("heatSource"))),
    alphaEffName_(dict.lookup("alphaEff")),
    q_(Function1<scalar>::New("q", dict))
{
	Info<< "DDDDDDDDDD" << endl;
    //this->evaluate();
    //if (dict.found("value") && dict.found("gradient"))
    //{
    //    fvPatchField<scalar>::operator=(Field<scalar>("value", dict, p.size()));
    //    gradient() = Field<scalar>("gradient", dict, p.size());
    //}
    //else
    //{
    //    fvPatchField<scalar>::operator=(patchInternalField());
    //    gradient() = 0.0;
    //}
    const scalar t = this->db().time().timeOutputValue();
	Info<< "q = " << q_->value(t) << endl;
}

//template<class Type>
//Foam::turbulentHeatFluxTemperatureTimeVaryingFvPatchScalarField<Type>::turbulentHeatFluxTemperatureTimeVaryingFvPatchScalarField
//(
//    const fvPatch& p,
//    const DimensionedField<Type, volMesh>& iF,
//    const Field<Type>& fld
//)
//:
//    uniformFixedGradientFvPatchField<Type>(p, iF, fld)
//{}




turbulentHeatFluxTemperatureTimeVaryingFvPatchScalarField::turbulentHeatFluxTemperatureTimeVaryingFvPatchScalarField
(
    const turbulentHeatFluxTemperatureTimeVaryingFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    uniformFixedGradientFvPatchField<scalar>(ptf, p, iF, mapper),
    heatSource_(ptf.heatSource_),
    //q_(ptf.q_, mapper),
    alphaEffName_(ptf.alphaEffName_),
    q_(ptf.q_, false)
{

	Info<< "JESTEM2" << endl;
}


turbulentHeatFluxTemperatureTimeVaryingFvPatchScalarField::turbulentHeatFluxTemperatureTimeVaryingFvPatchScalarField
(
    const turbulentHeatFluxTemperatureTimeVaryingFvPatchScalarField& ptf
)
:
    uniformFixedGradientFvPatchField<scalar>(ptf),
    heatSource_(ptf.heatSource_),
    //q_(ptf.q_),
    alphaEffName_(ptf.alphaEffName_),
    q_(ptf.q_, false)
{

	Info<< "JESTEM1" << endl;
}


turbulentHeatFluxTemperatureTimeVaryingFvPatchScalarField::turbulentHeatFluxTemperatureTimeVaryingFvPatchScalarField
(
    const turbulentHeatFluxTemperatureTimeVaryingFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    uniformFixedGradientFvPatchField<scalar>(ptf, iF),
    heatSource_(ptf.heatSource_),
    //q_(ptf.q_),
    alphaEffName_(ptf.alphaEffName_),
    q_(ptf.q_, false)
{
	Info<< "JESTEM" << endl;
    // Evaluate the profile if defined
    if (ptf.q_.valid())
    {
        this->evaluate();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void turbulentHeatFluxTemperatureTimeVaryingFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const scalar t = this->db().time().timeOutputValue();

    const scalarField& alphaEffp =
        patch().lookupPatchField<volScalarField, scalar>(alphaEffName_);
	
	if (!db().foundObject<volScalarField>("cp"))
	{
		// retrieve (constant) specific heat capacity from transport dictionary
    	const IOdictionary& transportProperties =
    	    db().lookupObject<IOdictionary>("transportProperties");
    	const scalar rho(readScalar(transportProperties.lookup("rho")));
    	const scalar cp(readScalar(transportProperties.lookup("cp")));

		switch (heatSource_)
    	{
    	    case hsPower:
    	    {
    	        const scalar Ap = gSum(patch().magSf());
    	        gradient() = q_->value(t)/(Ap*rho*cp*alphaEffp);
    	        break;
    	    }
    	    case hsFlux:
    	    {
    	        gradient() = q_->value(t)/(rho*cp*alphaEffp);
    	        break;
    	    }
    	    default:
    	    {
    	        FatalErrorIn
    	        (
    	            "turbulentHeatFluxTemperatureFvPatchScalarField"
    	            "("
    	                "const fvPatch&, "
    	                "const DimensionedField<scalar, volMesh>&, "
    	                "const dictionary&"
    	            ")"
    	        )   << "Unknown heat source type. Valid types are: "
    	            << heatSourceTypeNames_ << nl << exit(FatalError);
    	    }
    	}

	}
	else
	{
		const scalarField& cp0 =
			patch().lookupPatchField<volScalarField, scalar>("cp");
		const scalarField& rho =
			patch().lookupPatchField<volScalarField, scalar>("rho");

		switch (heatSource_)
    	{
    	    case hsPower:
    	    {
    	        const scalar Ap = gSum(patch().magSf());
    	        gradient() = q_->value(t)/(Ap*rho*cp0*alphaEffp);
    	        break;
    	    }
    	    case hsFlux:
    	    {
    	        gradient() = q_->value(t)/(rho*cp0*alphaEffp);
    	        break;
    	    }
    	    default:
    	    {
    	        FatalErrorIn
    	        (
    	            "turbulentHeatFluxTemperatureFvPatchScalarField"
    	            "("
    	                "const fvPatch&, "
    	                "const DimensionedField<scalar, volMesh>&, "
    	                "const dictionary&"
    	            ")"
    	        )   << "Unknown heat source type. Valid types are: "
    	            << heatSourceTypeNames_ << nl << exit(FatalError);
    	    }
    	}
	}

    //uniformFixedGradientFvPatchField<scalar>::updateCoeffs();
    fixedGradientFvPatchField<scalar>::updateCoeffs();
}


void turbulentHeatFluxTemperatureTimeVaryingFvPatchScalarField::write(Ostream& os) const
{
    uniformFixedGradientFvPatchField<scalar>::write(os);
    os.writeKeyword("heatSource") << heatSourceTypeNames_[heatSource_]
        << token::END_STATEMENT << nl;
    q_->writeData(os);
    os.writeKeyword("alphaEff") << alphaEffName_ << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    turbulentHeatFluxTemperatureTimeVaryingFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
