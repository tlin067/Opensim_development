// test.cpp   

/* Copyright (c)  2009 Stanford University
* Use of the OpenSim software in source form is permitted provided that the following
* conditions are met:
*   1. The software is used only for non-commercial research and education. It may not
*     be used in relation to any commercial activity.
*   2. The software is not distributed or redistributed.  Software distribution is allowed 
*     only through https://simtk.org/home/opensim.
*   3. Use of the OpenSim software or derivatives must be acknowledged in all publications,
*      presentations, or documents describing work in which OpenSim or derivatives are used.
*   4. Credits to developers may not be removed from executables
*     created from modifications of the source.
*   5. Modifications of source code must retain the above copyright notice, this list of
*     conditions and the following disclaimer. 
* 
*  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
*  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
*  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
*  SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
*  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
*  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
*  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
*  OR BUSINESS INTERRUPTION) OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY
*  WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/* 
*  Below is an example of an OpenSim application that provides its own 
*  main() routine.  
*/

//==============================================================================
//==============================================================================
#include <OpenSim/OpenSim.h>
#include "SimTKsimbody.h"
#include "CoupledBushingForceEDIT.h"
#include "my_classes.h"
#include "MyObjectiveFunc.h"
#include "Model_setup.h"
#include "Force_Plugin\FunctionBasedBushingForce.h"

//#include <Cfsqp\cfsqpusr.h>
#include <time.h>    // for clock()
#include <iostream>

#ifdef _WIN32
	#include <direct.h> 
	#define GetCurrentDir _getcwd 
#else 
	#include <unistd.h> 
	#define GetCurrentDir getcwd 
#endif 

double bestSoFar = Infinity;
int stepCount = 0;
//______________________________________________________________________________
/**
* Run an fwd simulation optimization
*/
int main()
{
	try {

		//////////////////////////////////////////////////////////////////////////////////////////////
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
		//////////////////////////////////////////////////////////////////////////////////////////////

		string osim_filename = "";
		string force_filename = "";
		string marker_filename = "";
		string ik_setup_filename = "";

		

		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
		//// Set up file directories for LINUX
		//      //Current directory... 
		//      char cCurrentPath[100]; 
		//      if (!GetCurrentDir(cCurrentPath, sizeof(cCurrentPath))) 
		//           { 
		//           return 0; 
		//           } 
		//      string currentPath = cCurrentPath;         
		//      // Define osim file location 
		//      string osim_fd = currentPath.append("/Opensim_development/"); 
		//      currentPath = cCurrentPath; 

		//       
		//      // Specify force and trc mocap file 
		//      string fd = currentPath.append("/Lamb2 data files/"); 
		//      currentPath = cCurrentPath; 
		//      //FLEXION 
		//      string expt_file = "l2flexv2"; 
		//      osim_filename = osim_filename.append(osim_fd).append("OSIM_contrained.osim"); // for flex     
		//      //Open existing XML model 
		//      Model osimModel(osim_filename);  
		//      osimModel.printBasicInfo(cout);         
		//      osimModel.updCoordinateSet().get("t2ToGND_FE").setDefaultValue(-Pi/2); // to flip up model so initial ik guess is close 
		//      osimModel.updCoordinateSet().get("t2ToGND_LB").setDefaultValue(-Pi/2); 
		//      string output_fd = currentPath.append("/OpenSim Output/flex_output.sto"); 

        //force_filename = force_filename.append(fd).append(expt_file).append(".mot"); 
        //marker_filename = marker_filename.append(fd).append(expt_file).append(".trc"); 
        //ik_setup_filename = ik_setup_filename.append(osim_fd).append(expt_file).append("_initial_ik.xml"); 

		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//

		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
		// Set up file directories for WINDOWS...

		// Define osim file location
		string osim_fd = "T:/Lamb expts 2012/Lamb experiments 2012/23rd - 24th July Expts/Opensim/OpenSim new versions (from Version10)/Opensim_development/";

		// Specify force and trc mocap file
		string fd = "T:/Lamb expts 2012/Lamb experiments 2012/23rd - 24th July Expts/Opensim/Lamb2 data files/";

		//FLEXION
		string expt_file = "l2flexv2";
		osim_filename = osim_filename.append(osim_fd).append("OSIM_contrained.osim"); // for flex	
		//Open existing XML model
		Model osimModel(osim_filename); 
		osimModel.printBasicInfo(cout);		
		osimModel.updCoordinateSet().get("t2ToGND_FE").setDefaultValue(-Pi/2); // to flip up model so initial ik guess is close
		osimModel.updCoordinateSet().get("t2ToGND_LB").setDefaultValue(-Pi/2);
		string output_fd = "T:/Lamb expts 2012/Lamb experiments 2012/23rd - 24th July Expts/Opensim/OpenSim new versions (from Version10)/OpenSim Output/flex_output.sto";
		double ti,tf;
		ti = 2;
		tf = 7;	

		//// EXTENSION
		//string expt_file = "l2extv2a";
		//osim_filename = osim_filename.append(osim_fd).append("OSIM_contrained.osim"); // for extension
		////Open existing XML model
		//Model osimModel(osim_filename); 
		//osimModel.printBasicInfo(cout);		
		//osimModel.updCoordinateSet().get("t2ToGND_FE").setDefaultValue(0); // to flip up model so initial ik guess is close
		//osimModel.updCoordinateSet().get("t2ToGND_LB").setDefaultValue(0);
		//string output_fd = "T:/Lamb expts 2012/Lamb experiments 2012/23rd - 24th July Expts/Opensim/OpenSim new versions (from Version10)/OpenSim Output/ext_output.sto";
		//double ti,tf;
		//ti = 5;
		//tf = 18;

		force_filename = force_filename.append(fd).append(expt_file).append(".mot");
		marker_filename = marker_filename.append(fd).append(expt_file).append(".trc");
		ik_setup_filename = ik_setup_filename.append(osim_fd).append(expt_file).append("_initial_ik.xml");

		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//	

		OpenSim::CoordinateSet& CS = osimModel.updCoordinateSet();
		BodySet& bodyset = osimModel.updBodySet();

		// Create initial set up class
		Model_setup* modelSetup = new Model_setup();
		
		//add external force measurements
		modelSetup->addExtForce(osimModel, force_filename);
		//add default coupled bushing so properties can be changed later.
		modelSetup->addInitialCoupledBushing(osimModel);
		modelSetup->addPKAnalyses(osimModel);

		//osimModel.print("T:/Lamb expts 2012/Lamb experiments 2012/23rd - 24th July Expts/Opensim/OpenSim new versions (from Version10)/OpenSim Output/XXXmodel output.osim");

		cout<<"\n\n"<<osimModel.updAnalysisSet();
		cout<<"\n\n";


		//////////////////////////////////////////////////////////////////////////////////////////////
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
		//////////////////////////////////////////////////////////////////////////////////////////////


		// Run initial ik step to find Initial conditions for forward simulation

		SimTools* simtools = new SimTools();

		Array<double> ICs;

		// Returns state values at initial time -> ICs has length of coordinate set (computed by single ik step)
		ICs = simtools->ik_constrain_torso(osimModel,marker_filename,ik_setup_filename,ti,tf);

		// Display initial conditions for forward simulation
		cout<<"\n\n Initial Conditions = "<<ICs<<"\n\n";


		// Sets torso coordinates to be set @ value specified by ik analysis (assuming torso does not move....)
		for(int i = 0; i<6; i++)
		{
			OpenSim::Function* function_ik = new OpenSim::Constant(ICs[i]);
			// First 6 coordinates in cooridnate set correspond to rotx,roty,rotz,tx,ty,tz
			CS[i].setPrescribedFunction(*function_ik);
			CS[i].setDefaultValue(ICs[i]);
			CS[i].setDefaultIsPrescribed(true);
		} 

		// Prints out ammended osim file -> Torso coords will be prescribed to some constant value
		//osimModel.print("T:/Lamb expts 2012/Lamb experiments 2012/23rd - 24th July Expts/Opensim/Version4/test.osim");


		//////////////////////////////////////////////////////////////////////////////////////////////
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
		//////////////////////////////////////////////////////////////////////////////////////////////

		//////////////////////////////////////////////////////////////////////////////////////////////
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
		//////////////////////////////////////////////////////////////////////////////////////////////

		//// RUN SIMULATION for a given set of bushing force properties

		cout<<"\n\n Force set: "<<osimModel.getForceSet()<<"\n\n";

		//// Define bushing properties
		//// Define translational stiffness (t_stiff - x,y,z), rotational stiffness (r_stiff - x,y,z) and corresponding damping
		//Vec3 t_stiff(0),r_stiff(0.01),t_damp(0),r_damp(10);

		// Open trc file and store in class MarkerData
		MarkerData md = OpenSim::MarkerData(marker_filename);

		// Create storage object with marker data (md) in it
		Storage data_trc;
		md.makeRdStorage(data_trc);
		// crop storage object to only include data in time frame of simulation
		data_trc.crop(ti,tf); // -> data_trc is used to calculate rms error between model marker and motion capture markers

		// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //

		// Initialise and run optimization
		int numParams = 22;
		// Parameters:
		// - Theta_STAR
		// - k1		
		// - k2
		// - bushing offset

		MyOptimizerSystem sys(numParams,osimModel,data_trc,ti,tf,ICs,output_fd);

//(OpenSim::PointKinematics&)osimModel.getAnalysisSet().get(0),(OpenSim::PointKinematics&)osimModel.getAnalysisSet().get(1),(OpenSim::PointKinematics&)osimModel.getAnalysisSet().get(2),(OpenSim::PointKinematics&)osimModel.getAnalysisSet().get(3),(OpenSim::PointKinematics&)osimModel.getAnalysisSet().get(4),(OpenSim::PointKinematics&)osimModel.getAnalysisSet().get(5),output_fd);

		Real f = NaN;

		//Vector lower_bounds(numParams);
		//Vector upper_bounds(numParams);
		Vector guess(numParams);
		//// set bounds
		for (int i=0;i<(numParams-1);i++){
			guess[i] = 0;//0.0;//*Pi/180;
			//upper_bounds[i] = 100;//1.0;//*Pi/180;
		}

		guess[21] = 0.8;
		// set parameter limits
		//sys.setParameterLimits(lower_bounds,upper_bounds);

		// initialise vector of initial parameter guesses
		

		//guess[0] = 5;
		//guess[1] = 0; 
		//guess[2] = 0; 
		//guess[3] = 0;
		//guess[4] = 40;
		//guess[5] = 0; 
		//guess[6] = 0;
		//guess[7] = 0; 
		//guess[8] = 40; 

		//guess[9] = 0.65;


		//// Try optimisation
		clock_t t1,t2,t3,t4;
		t1 = clock();
		//try{

		//	// intialise optimizer
		//	Optimizer opt(sys, SimTK::InteriorPoint);
		//	opt.setDiagnosticsLevel(5);

		//	// Optimisation settings					
		//	opt.setConvergenceTolerance(1e-5);
		//	opt.setMaxIterations(1000);
		//	opt.useNumericalGradient(true);
		//	opt.setLimitedMemoryHistory(500);

		//	// return optimum solution
		//	f = opt.optimize(guess);

		//	cout<<"\nf = "<<f;
		//	cout<<"\nguess = "<<guess;
		//}
		//catch(const std::exception& e) {
		//	std::cout << "OptimizationExample.cpp Caught exception :"  << std::endl;
		//	std::cout << e.what() << std::endl;
		//}
		t2 = clock();
		 //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //	

		cout<<"\n\nOptimiszation time: "<<((float)t2-(float)t1)/ CLOCKS_PER_SEC << " seconds";

		// Run simulation and save point kinematic reporter in data_sim storage object. Print simulation results to file.
		Array<Array<double> > pk_data;

		Vector PARAMS(numParams);
		for (int i = 0; i<numParams; i++){
			PARAMS[i] = guess[i];
		}

		OpenSim::PointKinematics *m1h = (OpenSim::PointKinematics*)&osimModel.updAnalysisSet().get("m1");
		OpenSim::PointKinematics *m2h = (OpenSim::PointKinematics*)&osimModel.updAnalysisSet().get("m2");
		OpenSim::PointKinematics *m3h = (OpenSim::PointKinematics*)&osimModel.updAnalysisSet().get("m3");
		OpenSim::PointKinematics *m4h = (OpenSim::PointKinematics*)&osimModel.updAnalysisSet().get("m4");
		OpenSim::PointKinematics *m5h = (OpenSim::PointKinematics*)&osimModel.updAnalysisSet().get("m5");
		OpenSim::PointKinematics *m6h = (OpenSim::PointKinematics*)&osimModel.updAnalysisSet().get("m6");

		bool toWrite = 1;
		//pk_data = simtools->RunSimulation_LIMITSTOP(osimModel,PARAMS,ti,tf,ICs,toWrite,output_fd,*m1h,*m2h,*m3h,*m4h,*m5h,*m6h);
		double rms = simtools->RunSimulation_wRMS(data_trc,osimModel,PARAMS,ti,tf,ICs,false,output_fd,*m1h,*m2h,*m3h,*m4h,*m5h,*m6h);
		cout<<"\nRMS: "<<rms<<endl;

		////cout<<"\n\npk_data: "<<pk_data;
		//t3 = clock();
		//cout<<"\n\nSimulation time: "<<((float)t3-(float)t2)/ CLOCKS_PER_SEC << " seconds";
		//double rms = simtools->Calc_rms_VERSION2(data_trc,ti,tf,pk_data);
		//t4 = clock();
		//cout<<"\n\nRMS time: "<<((float)t4-(float)t3)/ CLOCKS_PER_SEC << " seconds";

		//cout<<"\n\nRMS ERROR: "<<rms;
		//
		//osimModel.print("T:/Lamb expts 2012/Lamb experiments 2012/23rd - 24th July Expts/Opensim/Version5/test.osim");

		//ofstream myfile;
		//myfile.open("T:/Lamb expts 2012/Lamb experiments 2012/23rd - 24th July Expts/Opensim/Version6/Model/fileConstrainedwBushingAndLimit.txt");
		//double rms = 0;
		//int n = 25;
		//double val = 0;
		//Vector k(n,(double)0);
		//Vector obj(n,(double)0);
		//for (int i=0; i<n; i++){
		//	t1 = clock();
		//	val = i*0.5 + 4;
		//	t_slackL_opt = val;
		//	//Vec3 r_stiff_opt(val);
		//	pk_data = simtools->RunSimulation_LIMITSTOP(osimModel,t_stiff,r_stiff_opt,t_damp,r_damp,ti,tf,t_slackL_opt,vert_mass_opt,head_mass_opt,ICs,false,osim_fd,*m1h,*m2h,*m3h,*m4h,*m5h,*m6h,fibreL_opt);
		//	rms = simtools->Calc_rms_VERSION2(data_trc,ti,tf,pk_data);
		//	k[i] = val;
		//	obj[i] = rms;
		//	cout<<i<<"\n";
		//	myfile<<k[i]<<"\t"<<obj[i]<<"\n";
		//	t2 = clock();
		//	cout<<k[i]<<"\t"<<obj[i]<<"\n";
		//	cout<<"\n\nLoop time: "<<((float)t2-(float)t1)/ CLOCKS_PER_SEC << " seconds";
		//}

		//myfile.close();
		//cout<<"\n\nk = "<<k;
		//cout<<"\n\nobj = "<<obj;


		//ofstream myfile_2param;
		//myfile_2param.open("T:/Lamb expts 2012/Lamb experiments 2012/23rd - 24th July Expts/Opensim/Version6/Model/fileConstrainedwBushingAndLimit_2param_v2.txt");
		//double rms = 0;
		//int n =10;
		//int nn = 20;
		//double val_k = 0, val_m = 0;
		//Vector k(n*nn,(double)0);
		//Vector m(n*nn,(double)0);
		//Vector obj(n*nn,(double)0);
		//for (int i=0; i<n; i++){
		//	for (int j=0;j<nn; j++){
		//		val_k = i*0.05 + 0.3;
		//		val_m = j*0.5 + 5;
		//		t_slackL_opt = val_m;
		//		head_mass_opt = val_k;
		//		//Vec3 r_stiff_opt(val_k);
		//		t_slackL_opt = val_m;
		//		pk_data = simtools->RunSimulation_LIMITSTOP(osimModel,t_stiff,r_stiff_opt,t_damp,r_damp,ti,tf,t_slackL_opt,vert_mass_opt,head_mass_opt,ICs,false,osim_fd,*m1h,*m2h,*m3h,*m4h,*m5h,*m6h,fibreL_opt);
		//		rms = simtools->Calc_rms_VERSION2(data_trc,ti,tf,pk_data);
		//		k[i*n + j] = val_k;
		//		m[i*n + j] = val_m;
		//		obj[i*n + j] = rms;
		//		cout<<(i)<<"\t"<<j<<"\n";
		//		cout<<k[i*n + j]<<"\n";
		//		cout<<m[i*n + j]<<"\n";
		//		cout<<rms<<"\n";
		//		myfile_2param<<k[i*n + j]<<"\t"<<m[i*n + j]<<"\t"<<obj[i*n + j]<<"\n";
		//	}
		//}

		//myfile_2param.close();
		//cout<<"\n\nk = "<<k;
		//cout<<"\n\nm = "<<m;
		//cout<<"\n\nobj = "<<obj;
		//////////////////////////////////////////////////////////////////////////////////////////////
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
		//////////////////////////////////////////////////////////////////////////////////////////////

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////		
	}
	catch (OpenSim::Exception ex)
	{
		std::cout << ex.getMessage() << std::endl;
		return 1;
	}
	catch (std::exception ex)
	{
		std::cout << ex.what() << std::endl;
		return 1;
	}
	catch (...)
	{
		std::cout << "UNRECOGNIZED EXCEPTION" << std::endl;
		return 1;
	}

	//std::cout << "main() routine time = " << 1.e3*(std::clock()-startTime)/CLOCKS_PER_SEC << "ms\n";

	std::cout << "\n\nOpenSim example completed successfully.\n";
	std::cin.get();
	return 0;
}
