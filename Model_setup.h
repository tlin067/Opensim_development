#include <OpenSim/OpenSim.h>
#include "SimTKsimbody.h"

using namespace OpenSim;
using namespace SimTK;
using namespace std;

// my simulation tools
class Model_setup {
public:
	Model_setup();
	
	void addExtForce(OpenSim::Model& osimModel, string& force_filename); 
	void addInitialCoupledBushing(OpenSim::Model& osimModel); 
	void addPKAnalyses(OpenSim::Model& osimModel); 


		
private:
};