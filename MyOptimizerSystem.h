#include <OpenSim/OpenSim.h>
#include "SimTKsimbody.h"
#include "OpenSim/Tools/InverseKinematicsTool.h"
#include <OpenSim/Simulation/CoordinateReference.h>
#include <OpenSim/Simulation/InverseKinematicsSolver.h>
//#include "my_classes.h"


using namespace OpenSim;
using namespace SimTK;
using namespace std;


// optimizaer system class
class MyOptimizerSystem : public OptimizerSystem{
public:
	
	///* Constructor class. Parameters passed are accessed in the objectiveFunc() class. */CoordinateCouplerConstraint& FE21_c,  FE21_c(FE21_c),  Storage& data_trc
	MyOptimizerSystem(int& numParams, Model& osimModel, Storage& data_trc, double& ti, double& tf, Array<double>& ICs, string& output_fd);
	int objectiveFunc(  const Vector &newGuess, const bool new_coefficients, Real& f ) const; 
	int gradientFunc( const Vector &newGuess, bool new_coefficients, Vector &gradient ) const;
	
	//int constraintFunc( const Vector &newGuess, const bool new_coefficients, Vector &constraints) const;
	//int constraintJacobian( const Vector& newGuess, const bool new_coefficients, Matrix& jac)  const;



protected: 
	//CoordinateCouplerConstraint& FE21_c;
	//SinOmegaX& sin_func;
	int& numParams;
	//State& si;
	OpenSim::Model& osimModel;
	OpenSim::Storage& data_trc;
	//Storage& data_trc;
	double& ti;
	double& tf;
	OpenSim::Array<double>& ICs;
	/*OpenSim::AnalysisSet& AnalysisSet;*/
	OpenSim::CoordinateLimitForce& LIMIT_FORCE;
	string& output_fd;
	

};