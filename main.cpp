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

#include "my_classes.h"
#include "MyObjectiveFunc.h"

#include <Cfsqp\cfsqpusr.h>
#include <time.h>    // for clock()
#include <iostream>

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
		SimTools* simtools = new SimTools(); // Instance of simtools class to use my functions

		string osim_filename = "";
		string force_filename_flex = "";
		string marker_filename_flex = "";
		string ik_setup_filename_flex = "";

		string force_filename_ext = "";
		string marker_filename_ext = "";
		string ik_setup_filename_ext = "";

		// Set up file directories
		
		// Define osim file location
		string osim_fd = "T:/Lamb expts 2012/Lamb experiments 2012/23rd - 24th July Expts/Opensim/Version9/";

		// Specify force and trc mocap file
		string fd = "T:/Lamb expts 2012/Lamb experiments 2012/23rd - 24th July Expts/Opensim/Lamb2 data files/";
		

		osim_filename = osim_filename.append(osim_fd).append("Version9_contrained.osim"); // for flex		
		//Open existing XML model
		Model osimModel(osim_filename); 
		Model osimModelE(osim_filename); 
		osimModel.printBasicInfo(cout);		
		
		//FLEXION
		string expt_file_flex = "l2flexv2";


		string output_fd_flex = "T:/Lamb expts 2012/Lamb experiments 2012/23rd - 24th July Expts/Opensim/Version9/flex_output.sto";
		double tiFlex,tfFlex;
		tiFlex = 2; // Static equilibrium
		tfFlex = 15;		

		// EXTENSION
		string expt_file_ext = "l2extv2a";

		string output_fd_ext = "T:/Lamb expts 2012/Lamb experiments 2012/23rd - 24th July Expts/Opensim/Version9/ext_output.sto";
		double tiExt,tfExt;
		tiExt = 5;
		tfExt = 18;

		//flex filenames
		force_filename_flex = force_filename_flex.append(fd).append(expt_file_flex).append(".mot");
		marker_filename_flex = marker_filename_flex.append(fd).append(expt_file_flex).append(".trc");
		ik_setup_filename_flex = ik_setup_filename_flex.append(osim_fd).append(expt_file_flex).append("_initial_ik.xml");

		//ext filenames
		force_filename_ext = force_filename_ext.append(fd).append(expt_file_ext).append(".mot");
		marker_filename_ext = marker_filename_ext.append(fd).append(expt_file_ext).append(".trc");
		ik_setup_filename_ext = ik_setup_filename_ext.append(osim_fd).append(expt_file_ext).append("_initial_ik.xml");

		// Create storage object for force file
		OpenSim::Storage* force_storage_ext = new Storage(force_filename_ext);
		OpenSim::Storage* force_storage_flex = new Storage(force_filename_flex);

		// smooths out all data....
		force_storage_flex->smoothSpline(3,1); // order = 3, cutoff freq = 1Hz
		force_storage_flex->resampleLinear(0.1); // resample to 10Hz

		force_storage_ext->smoothSpline(3,1); // order = 3, cutoff freq = 1Hz
		force_storage_ext->resampleLinear(0.1); // resample to 10Hz
		//force_storage->print("T:/Lamb expts 2012/Lamb experiments 2012/23rd - 24th July Expts/Opensim/Version4/force_flex_smooth.sto");
		//////////////////////////////////////////////////////////////////////////////////////////////
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
		//////////////////////////////////////////////////////////////////////////////////////////////

	

		//OpenSim::CoordinateSet& CS = osimModel.updCoordinateSet();
		// Run initial ik step to find Initial conditions for forward simulation
		Array<double> ICs_Flex, ICs_Ext;

		osimModel.updCoordinateSet().get("t2ToGND_FE").setDefaultValue(0); // to flip up model so initial ik guess is close
		osimModel.updCoordinateSet().get("t2ToGND_LB").setDefaultValue(0);
		ICs_Ext = simtools->ik_constrain_torso(osimModel,marker_filename_ext,ik_setup_filename_ext,tiExt,tfExt);
		cout<<"\n\next ic: :"<<ICs_Ext<<endl;
        osimModel.print("T:/Lamb expts 2012/Lamb experiments 2012/23rd - 24th July Expts/Opensim/Version9/test1.osim");		
		

		OpenSim::CoordinateSet& CS = osimModel.updCoordinateSet();//.get("t2ToGND_FE").setDefaultIsPrescribed(false);

		for(int i = 0; i<6; i++)
		{
			// First 6 coordinates in cooridnate set correspond to rotx,roty,rotz,tx,ty,tz
			CS[i].setDefaultIsPrescribed(false);
			CS[i].setDefaultValue(0.0);

		} 
		// For Flex
		osimModel.updCoordinateSet().get("t2ToGND_FE").setDefaultValue(-Pi/2); // to flip up model so initial ik guess is close
		osimModel.updCoordinateSet().get("t2ToGND_LB").setDefaultValue(-Pi/2);
		osimModel.print("T:/Lamb expts 2012/Lamb experiments 2012/23rd - 24th July Expts/Opensim/Version9/testa.osim");
		// Returns state values at initial time -> ICs has length of coordinate set (computed by single ik step)
		ICs_Flex = simtools->ik_constrain_torso(osimModel,marker_filename_flex,ik_setup_filename_flex,tiFlex,tfFlex);	
		osimModel.print("T:/Lamb expts 2012/Lamb experiments 2012/23rd - 24th July Expts/Opensim/Version9/test.osim");
		cout<<"\n\nflex ic: :"<<ICs_Flex<<endl;
		
		// Set up point kinematic analyses -> to track marker positions for use in objective function
		//AnalysisSet pk_set = 
		//simtools->AddPKAtest(osimModel);

		//cout<<"name : "<<osimModel.getNumAnalyses();

		//cout<<"\n\nname : "<<osimModel.updAnalysisSet().get("m1");

		////// USE force file to prescribe forces to head_markers body
		////simtools->ApplyForce(*force_storage_flex,osimModel);


		// Create measurement data structure for each simulation // 
		// Open trc file and store in class MarkerData
		MarkerData md_flex = OpenSim::MarkerData(marker_filename_flex);
		MarkerData md_ext = OpenSim::MarkerData(marker_filename_ext);

		// Create storage object with marker data (md) in it
		Storage data_trc_flex, data_trc_ext;
		md_flex.makeRdStorage(data_trc_flex);
		md_ext.makeRdStorage(data_trc_ext);
		// crop storage object to only include data in time frame of simulation
		data_trc_flex.crop(tiFlex,tfFlex); // -> data_trc is used to calculate rms error between model marker and motion capture markers
		data_trc_ext.crop(tiExt,tfExt); // -> data_trc is used to calculate rms error between model marker and motion capture markers



		//osimModel.print("T:/Lamb expts 2012/Lamb experiments 2012/23rd - 24th July Expts/Opensim/Version9/test.osim");

		////////////////////////////////////////////////////////////////////////////////////////////////
		////%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
		////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////
		//		 Create list of point kinematics reporters which match marker positions from osim model
		PointKinematics* m1h = new PointKinematics(&osimModel);
		PointKinematics* m2h = new PointKinematics(&osimModel);
		PointKinematics* m3h = new PointKinematics(&osimModel);
		PointKinematics* m4h = new PointKinematics(&osimModel);
		PointKinematics* m5h = new PointKinematics(&osimModel);
		PointKinematics* m6h = new PointKinematics(&osimModel);


		// points in head_marker ref frame specified by inspection
		Vec3 p1(0.071700000,0.00000000,-0.054772000);
		Vec3 p2(-0.071700000,0.00000000,-0.038524000);
		Vec3 p3(0.071700000,0.00000000,0.005228);
		Vec3 p4(-0.0300000,0.00000000,0.026928);
		Vec3 p5(0.0300000,0.00000000,0.026928);
		Vec3 p6(-0.071700000,0.00000000,-0.008524000);

		m1h->setBodyPoint("head_markers",p1);
		m2h->setBodyPoint("head_markers",p2);
		m3h->setBodyPoint("head_markers",p3);
		m4h->setBodyPoint("head_markers",p4);
		m5h->setBodyPoint("head_markers",p5);
		m6h->setBodyPoint("head_markers",p6);
		
		m1h->setRelativeToBody(&osimModel.updBodySet().get("ground"));
		m2h->setRelativeToBody(&osimModel.updBodySet().get("ground"));
		m3h->setRelativeToBody(&osimModel.updBodySet().get("ground"));
		m4h->setRelativeToBody(&osimModel.updBodySet().get("ground"));
		m5h->setRelativeToBody(&osimModel.updBodySet().get("ground"));
		m6h->setRelativeToBody(&osimModel.updBodySet().get("ground"));

		m1h->setName("m1");
		m2h->setName("m2");
		m3h->setName("m3");
		m4h->setName("m4");
		m5h->setName("m5");
		m6h->setName("m6");

		osimModel.addAnalysis(m1h);
		osimModel.addAnalysis(m2h);
		osimModel.addAnalysis(m3h);
		osimModel.addAnalysis(m4h);
		osimModel.addAnalysis(m5h);
		osimModel.addAnalysis(m6h);

		//////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////

		//// RUN SIMULATION for a given set of bushing force properties

		//// Define bushing properties
		//// Define translational stiffness (t_stiff - x,y,z), rotational stiffness (r_stiff - x,y,z) and corresponding damping
		Vec3 t_stiff(0),r_stiff(0.01),t_damp(0),r_damp(10);

		//// initialise vector of initial parameter guesses
		//int numParams = 4;
		//Vector guess(numParams);
	
		//guess[0] = 6.82;//6.00553; // theta star (degrees)
		//guess[1] = 1.22; // k1 (N*m/degree)
		//guess[2] = 7.29; // k2 (N*m/degree)
		//guess[3] = 0.22;//0.409055; // damping (Nm/(deg/sec))

		//Vector PARAMS(numParams);
		//for (int i = 0; i<numParams; i++){
		//	PARAMS[i] = guess[i];
		//}
		////string fd = "";
		//simtools->RunSimulation_FlexExt(osimModel,PARAMS,tiFlex,tiExt,tfFlex,tfExt,ICs_Flex,ICs_Ext,false,fd,*force_storage_flex,*force_storage_ext,*m1h,*m2h,*m3h,*m4h,*m5h,*m6h,data_trc_flex,data_trc_ext);
		//
		
		// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
		//
		//// Initialise and run optimization
		int numParams = 5;
		//// Parameters:
		//// - Theta_STAR
		//// - k1		
		//// - k2
		//// - damping
		//MyOptimizerSystem sys(numParams,osimModel,data_trc,ti,tf,ICs,*m1h,*m2h,*m3h,*m4h,*m5h,*m6h,output_fd);

		MyOptimizerSystem sys(numParams,osimModel,tiFlex,tiExt,tfFlex,tfExt,ICs_Flex,ICs_Ext,*force_storage_flex,*force_storage_ext,*m1h,*m2h,*m3h,*m4h,*m5h,*m6h,data_trc_flex,data_trc_ext);

		Real f = NaN;

		Vector lower_bounds(numParams);
		Vector upper_bounds(numParams);

		// theta star flex
		lower_bounds[0] = 2;//0.0;//*Pi/180;
		upper_bounds[0] = 12;//1.0;//*Pi/180;

		// theta star ext 
		lower_bounds[1] = 2;//0.0;//*Pi/180;
		upper_bounds[1] = 12;//1.0;//*Pi/180;

		// k1
		lower_bounds[2] = 1e-6;//0.0;//*Pi/180;
		upper_bounds[2] = 15;//1.0;//*Pi/180;

		// k2
		lower_bounds[3] = 0.1;
		upper_bounds[3] = 100;

		// damping
		lower_bounds[4] = -2;
		upper_bounds[4] = 2;

		//head mass
		//lower_bounds[3] = -2.0;
		//upper_bounds[3] = 2.0;

		//// theta star sk
		//lower_bounds[4] = 6;//0.0;//*Pi/180;
		//upper_bounds[4] = 30;//1.0;//*Pi/180;

		//// k1 sk
		//lower_bounds[5] = 1e-6;//0.0;//*Pi/180;
		//upper_bounds[5] = 15;//1.0;//*Pi/180;

		//// k2 sk
		//lower_bounds[6] = 0.1;
		//upper_bounds[6] = 100;

		// set parameter limits
		sys.setParameterLimits(lower_bounds,upper_bounds);

		// initialise vector of initial parameter guesses
		Vector guess(numParams);
	
		guess[0] = 6.82; //flexion theta star (degrees)
		guess[1] = 6.82; // extension theta star (degrees)
		guess[2] = 1.22; // k1 (N*m/degree)
		guess[3] = 7.29; // k2 (N*m/degree)
		guess[4] = 0.3; // bushing offset

		//guess[4] = 0.0; // head mass
		//		guess[4] = 45; // theta star sk (degrees)
		//guess[5] = guess[1]; // k1 sk (N*m/degree)
		//guess[6] = guess[2]; // k2 sk (N*m/degree)
		
		// Try optimisation
		clock_t t1,t2,t3,t4;
		t1 = clock();
		try{
		
			// intialise optimizer
			Optimizer opt(sys, SimTK::InteriorPoint);
			opt.setDiagnosticsLevel(5);
			
			// Optimisation settings					
			opt.setConvergenceTolerance(1e-3);
			opt.setMaxIterations(1000);
			opt.useNumericalGradient(true);
			opt.setLimitedMemoryHistory(500);
			
			// return optimum solution
			f = opt.optimize(guess);
			
			cout<<"\nf = "<<f;
			cout<<"\nguess = "<<guess;
		}
		catch(const std::exception& e) {
		std::cout << "OptimizationExample.cpp Caught exception :"  << std::endl;
		std::cout << e.what() << std::endl;
		}
		t2 = clock();
		//// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //	

		//cout<<"\n\nOptimiszation time: "<<((float)t2-(float)t1)/ CLOCKS_PER_SEC << " seconds";
		//
		//// Run simulation and save point kinematic reporter in data_sim storage object. Print simulation results to file.
		//Array<Array<double>> pk_data;

		Vector PARAMS(numParams);
		for (int i = 0; i<numParams; i++){
			PARAMS[i] = guess[i];
		}
		//string fd = "";
		simtools->RunSimulation_FlexExt(osimModel,PARAMS,tiFlex,tiExt,tfFlex,tfExt,ICs_Flex,ICs_Ext,false,fd,*force_storage_flex,*force_storage_ext,*m1h,*m2h,*m3h,*m4h,*m5h,*m6h,data_trc_flex,data_trc_ext);
			
		////
		//////cout<<"\n\npk_data: "<<pk_data;
		////t3 = clock();
		////cout<<"\n\nSimulation time: "<<((float)t3-(float)t2)/ CLOCKS_PER_SEC << " seconds";
		////double rms = simtools->Calc_rms_VERSION2(data_trc,ti,tf,pk_data);
		////t4 = clock();
		////cout<<"\n\nRMS time: "<<((float)t4-(float)t3)/ CLOCKS_PER_SEC << " seconds";

		////cout<<"\n\nRMS ERROR: "<<rms;
		////
		////osimModel.print("T:/Lamb expts 2012/Lamb experiments 2012/23rd - 24th July Expts/Opensim/Version5/test.osim");

		////ofstream myfile;
		////myfile.open("T:/Lamb expts 2012/Lamb experiments 2012/23rd - 24th July Expts/Opensim/Version6/Model/fileConstrainedwBushingAndLimit.txt");
		////double rms = 0;
		////int n = 25;
		////double val = 0;
		////Vector k(n,(double)0);
		////Vector obj(n,(double)0);
		////for (int i=0; i<n; i++){
		////	t1 = clock();
		////	val = i*0.5 + 4;
		////	t_slackL_opt = val;
		////	//Vec3 r_stiff_opt(val);
		////	pk_data = simtools->RunSimulation_LIMITSTOP(osimModel,t_stiff,r_stiff_opt,t_damp,r_damp,ti,tf,t_slackL_opt,vert_mass_opt,head_mass_opt,ICs,false,osim_fd,*m1h,*m2h,*m3h,*m4h,*m5h,*m6h,fibreL_opt);
		////	rms = simtools->Calc_rms_VERSION2(data_trc,ti,tf,pk_data);
		////	k[i] = val;
		////	obj[i] = rms;
		////	cout<<i<<"\n";
		////	myfile<<k[i]<<"\t"<<obj[i]<<"\n";
		////	t2 = clock();
		////	cout<<k[i]<<"\t"<<obj[i]<<"\n";
		////	cout<<"\n\nLoop time: "<<((float)t2-(float)t1)/ CLOCKS_PER_SEC << " seconds";
		////}

		////myfile.close();
		////cout<<"\n\nk = "<<k;
		////cout<<"\n\nobj = "<<obj;
		//

		////ofstream myfile_2param;
		////myfile_2param.open("T:/Lamb expts 2012/Lamb experiments 2012/23rd - 24th July Expts/Opensim/Version6/Model/fileConstrainedwBushingAndLimit_2param_v2.txt");
		////double rms = 0;
		////int n =10;
		////int nn = 20;
		////double val_k = 0, val_m = 0;
		////Vector k(n*nn,(double)0);
		////Vector m(n*nn,(double)0);
		////Vector obj(n*nn,(double)0);
		////for (int i=0; i<n; i++){
		////	for (int j=0;j<nn; j++){
		////		val_k = i*0.05 + 0.3;
		////		val_m = j*0.5 + 5;
		////		t_slackL_opt = val_m;
		////		head_mass_opt = val_k;
		////		//Vec3 r_stiff_opt(val_k);
		////		t_slackL_opt = val_m;
		////		pk_data = simtools->RunSimulation_LIMITSTOP(osimModel,t_stiff,r_stiff_opt,t_damp,r_damp,ti,tf,t_slackL_opt,vert_mass_opt,head_mass_opt,ICs,false,osim_fd,*m1h,*m2h,*m3h,*m4h,*m5h,*m6h,fibreL_opt);
		////		rms = simtools->Calc_rms_VERSION2(data_trc,ti,tf,pk_data);
		////		k[i*n + j] = val_k;
		////		m[i*n + j] = val_m;
		////		obj[i*n + j] = rms;
		////		cout<<(i)<<"\t"<<j<<"\n";
		////		cout<<k[i*n + j]<<"\n";
		////		cout<<m[i*n + j]<<"\n";
		////		cout<<rms<<"\n";
		////		myfile_2param<<k[i*n + j]<<"\t"<<m[i*n + j]<<"\t"<<obj[i*n + j]<<"\n";
		////	}
		////}

		////myfile_2param.close();
		////cout<<"\n\nk = "<<k;
		////cout<<"\n\nm = "<<m;
		////cout<<"\n\nobj = "<<obj;
		////////////////////////////////////////////////////////////////////////////////////////////////
		////%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
		////////////////////////////////////////////////////////////////////////////////////////////////

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
