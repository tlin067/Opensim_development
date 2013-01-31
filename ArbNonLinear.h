#ifndef __ArbNonLinear_h__
#define __ArbNonLinear_h__
/* -------------------------------------------------------------------------- *
 *                              OpenSim:  ArbNonLinear.h                              *
 * -------------------------------------------------------------------------- *
 * The OpenSim API is a toolkit for musculoskeletal modeling and simulation.  *
 * See http://opensim.stanford.edu and the NOTICE file for more information.  *
 * OpenSim is developed at Stanford University and supported by the US        *
 * National Institutes of Health (U54 GM072970, R24 HD065690) and by DARPA    *
 * through the Warrior Web program.                                           *
 *                                                                            *
 * Copyright (c) 2005-2012 Stanford University and the Authors                *
 * Author(s): Ajay Seth                                                       *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */


// INCLUDES
#include <string>
#include "OpenSim\Common\Function.h"
#include "OpenSim\Common\FunctionAdapter.h"
#include "OpenSim\Common\PropertyDbl.h"

namespace OpenSim {

//=============================================================================
//=============================================================================
/**
 * A class for representing a ArbNonLinear function.
 *
 * This class inherits from Function and so can be used as input to
 * any class requiring a Fuction as input. Implements f(x) = Ax^N + B
 *
 * @author Ajay Seth
 * @version 1.0  
 */

// Was in front of ArbNonLinear "OSIMCOMMON_API".... but only works without it
class ArbNonLinear : public Function {
OpenSim_DECLARE_CONCRETE_OBJECT(ArbNonLinear, Function);

//=============================================================================
// MEMBER VARIABLES
//=============================================================================
protected:

	PropertyDbl _AProp;
	double &_A;

	PropertyDbl _BProp;
	double &_B;

	PropertyDbl _NProp;
	double &_N;

//=============================================================================
// METHODS
//=============================================================================
public:
	//--------------------------------------------------------------------------
	// CONSTRUCTION
	//--------------------------------------------------------------------------
	ArbNonLinear() : _A(_AProp.getValueDbl()), _B(_BProp.getValueDbl()), _N(_NProp.getValueDbl()) { setupProperties();}
	// Convenience Constructor
	ArbNonLinear(double A, double B, double N) : _A(_AProp.getValueDbl()), _B(_BProp.getValueDbl()), _N(_NProp.getValueDbl()) {
		setupProperties();
		_A = A;  _B = B;  _N = N; 
	}
	// Copy Constructor
	ArbNonLinear(const ArbNonLinear &aFunc): _A(_AProp.getValueDbl()), _B(_BProp.getValueDbl()), _N(_NProp.getValueDbl()) {
			setupProperties();
			_A = aFunc._A;  _B = aFunc._B;  _N = aFunc._N; 
	};
	virtual ~ArbNonLinear() {};

private:
	void setupProperties() {
		_AProp.setName("A");
		_AProp.setComment("Coeff of the ArbNonLinear function");
		_AProp.setValue(1);
		_propertySet.append(&_AProp);

		_BProp.setName("B");
		_BProp.setComment("Offset of the ArbNonLinear function");
		_BProp.setValue(0);
		_propertySet.append(&_BProp);

		_NProp.setName("N");
		_NProp.setComment("Power of the ArbNonLinear function");
		_NProp.setValue(1);
		_propertySet.append(&_NProp);
	}

	//--------------------------------------------------------------------------
	// OPERATORS
	//--------------------------------------------------------------------------
public:
	ArbNonLinear& operator=(const ArbNonLinear &func)
	{
		Function::operator=(func);
		_A = func._A;  _B = func._B;  _N = func._N; 
		return(*this);
	}
	//--------------------------------------------------------------------------
	// SET AND GET
	//--------------------------------------------------------------------------
public:
	void setValue(double aValue);

	//--------------------------------------------------------------------------
	// EVALUATION
	//--------------------------------------------------------------------------
    virtual double calcValue(const SimTK::Vector& x) const
	{
		double returnValue;

		if (x[0]<0){
			returnValue = -(pow(-x[0],_N)*_A) + _B;
		}
		else if (x[0]>=0){
			returnValue = pow(x[0],_N)*_A + _B;
		}
		
		return returnValue;
	}
	
	double calcDerivative(const std::vector<int>& derivComponents, const SimTK::Vector& x) const
	{
		int n = derivComponents.size();

		return _N*_A*pow(fabs(x[0]),(_N-1));
	}

	SimTK::Function* createSimTKFunction() const {
		return new FunctionAdapter(*this);
	}
   
	int getArgumentSize() const {return 1;}
	int getMaxDerivativeOrder() const {return 1;}

//=============================================================================
};	// END class ArbNonLinear;
//=============================================================================
//=============================================================================

} // end of namespace OpenSim

#endif  // __ArbNonLinear_h__