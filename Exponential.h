#ifndef __Exponential_h__
#define __Exponential_h__
/* -------------------------------------------------------------------------- *
 *                              OpenSim:  Exponential.h                              *
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
 * A class for representing a Exponential function.
 *
 * This class inherits from Function and so can be used as input to
 * any class requiring a Fuction as input. Implements f(x) = Aexp^Bx + C
 *
 * @author Ajay Seth
 * @version 1.0  
 */

// Was in front of Exponential "OSIMCOMMON_API".... but only works without it
class Exponential : public Function {
OpenSim_DECLARE_CONCRETE_OBJECT(Exponential, Function);

//=============================================================================
// MEMBER VARIABLES
//=============================================================================
protected:

	PropertyDbl _AProp;
	double &_A;

	PropertyDbl _BProp;
	double &_B;

	PropertyDbl _CProp;
	double &_C;

//=============================================================================
// METHODS
//=============================================================================
public:
	//--------------------------------------------------------------------------
	// CONSTRUCTION
	//--------------------------------------------------------------------------
	Exponential() : _A(_AProp.getValueDbl()), _B(_BProp.getValueDbl()), _C(_CProp.getValueDbl()) { setupProperties();}
	// Convenience Constructor
	Exponential(double A, double B, double C) : _A(_AProp.getValueDbl()), _B(_BProp.getValueDbl()), _C(_CProp.getValueDbl()) {
		setupProperties();
		_A = A;  _B = B;  _C = C; 
	}
	// Copy Constructor
	Exponential(const Exponential &aFunc): _A(_AProp.getValueDbl()), _B(_BProp.getValueDbl()), _C(_CProp.getValueDbl()) {
			setupProperties();
			_A = aFunc._A;  _B = aFunc._B;  _C = aFunc._C; 
	};
	virtual ~Exponential() {};

private:
	void setupProperties() {
		_AProp.setName("A");
		_AProp.setComment("A of the Exponential function");
		_AProp.setValue(1);
		_propertySet.append(&_AProp);

		_BProp.setName("B");
		_BProp.setComment("B of the Exponential function");
		_BProp.setValue(1);
		_propertySet.append(&_BProp);

		_CProp.setName("C");
		_CProp.setComment("Offset of the exponential function");
		_CProp.setValue(0);
		_propertySet.append(&_CProp);
	}

	//--------------------------------------------------------------------------
	// OPERATORS
	//--------------------------------------------------------------------------
public:
	Exponential& operator=(const Exponential &func)
	{
		Function::operator=(func);
		_A = func._A;  _B = func._B;  _C = func._C; 
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
			returnValue = -(_A*exp(_B*-x[0]) -_A + _C);
		}
		else if (x[0]>=0){
			returnValue = _A*exp(_B*x[0]) -_A + _C;
		}
		
		return returnValue;
	}
	
	double calcDerivative(const std::vector<int>& derivComponents, const SimTK::Vector& x) const
	{
		int n = derivComponents.size();

		return _A*_B*exp(_B*x[0]);
	}

	SimTK::Function* createSimTKFunction() const {
		return new FunctionAdapter(*this);
	}
   
	int getArgumentSize() const {return 1;}
	int getMaxDerivativeOrder() const {return 1;}

//=============================================================================
};	// END class Exponential;
//=============================================================================
//=============================================================================

} // end of namespace OpenSim

#endif  // __Exponential_h__