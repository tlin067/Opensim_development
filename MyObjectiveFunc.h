#include <OpenSim/OpenSim.h>
#include "SimTKsimbody.h"
#include "OpenSim/Tools/InverseKinematicsTool.h"
#include <OpenSim/Simulation/CoordinateReference.h>
#include <OpenSim/Simulation/InverseKinematicsSolver.h>


#include "MyOptimizerSystem.h"


using namespace OpenSim;
using namespace SimTK;
using namespace std;

// This is a single scalar function of a vector of parameters.
class MyObjectiveFunc : public Differentiator::GradientFunction {
public:
	// constructor
    MyObjectiveFunc(int ny, const Vector &newGuess, const MyOptimizerSystem *myOptSys) ;
    
	// existing functions from template 
    void setTime(Real t);
    Real getTime() const;

	// Calls objective function from within an optimizersystem class: 
	Real Call_obj(Vector GUESS) const;

    // Must provide this pure virtual function.
	// evalulates optimizer class objective function so can numerically differentiate
    int f(const Vector& y, Real& fy) const;

private:
    Real time;
    const Vector &newGuess;
	const MyOptimizerSystem *myOptSys;
};

//************************************//
//************************************//
// scalar function from a scalar parameter: For testing:
class SinOmegaX : public Differentiator::ScalarFunction {
public:
	//SinOmegaX();
	SinOmegaX(const Vector &newGuess, const MyOptimizerSystem *myOptSys); //: w(omega) { }
	Real Call_obj() const;
    // Must provide this virtual function.
    int f(Real x, Real& fx) const ;//{
	

	//	fx = w;//std::sin(w*x);
    //    return 0; // success
    //}
protected:
    const Vector &newGuess;
	const MyOptimizerSystem *myOptSys;
};