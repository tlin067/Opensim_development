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

		string osim_filename = "";
		string force_filename = "";
		string marker_filename = "";
		string ik_setup_filename = "";

		// Set up file directories
		
		// Define osim file location
		string osim_fd = "T:/Lamb expts 2012/Lamb experiments 2012/23rd - 24th July Expts/Opensim/OpenSim new versions (from Version10)/Opensim_development/";
		
		// Specify force and trc mocap file
		string fd = "T:/Lamb expts 2012/Lamb experiments 2012/23rd - 24th July Expts/Opensim/Lamb2 data files/";
		
		//FLEXION
		string expt_file = "l2flexv2";
		osim_filename = osim_filename.append(osim_fd).append("OSIM_constrainedwFORCES.osim"); // for flex	
		//Open existing XML model
		Model osimModel(osim_filename); 
		osimModel.printBasicInfo(cout);		
		osimModel.updCoordinateSet().get("t2ToGND_FE").setDefaultValue(-Pi/2); // to flip up model so initial ik guess is close
		osimModel.updCoordinateSet().get("t2ToGND_LB").setDefaultValue(-Pi/2);
		string output_fd = "T:/Lamb expts 2012/Lamb experiments 2012/23rd - 24th July Expts/Opensim/OpenSim new versions (from Version10)/OpenSim Output/flex_output.sto";
		
		double ti,tf;
		ti = 2;
		tf = 15;		

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

		
		// Create storage object for force file
		OpenSim::Storage* force_storage = new Storage(force_filename);
		

		// smooths out all data....
		force_storage->smoothSpline(3,1); // order = 3, cutoff freq = 1Hz
		force_storage->resampleLinear(0.1); // resample to 10Hz
		//force_storage->print("T:/Lamb expts 2012/Lamb experiments 2012/23rd - 24th July Expts/Opensim/Version4/force_flex_smooth.sto");
		//////////////////////////////////////////////////////////////////////////////////////////////
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
		//////////////////////////////////////////////////////////////////////////////////////////////

		// USE force file to prescribe forces to head_markers body

		/////////////////////////////////////
		// NOTE: There seemed to be problem invoking "setdatacolumn" and saving result directly to double. Done a work around where save result as array, ...
		// and then convert to double for input into GCVSpline
		Array<double> time = Array<double>();
		Array<double> fx = Array<double>();
		Array<double> fy = Array<double>();
		Array<double> fz = Array<double>();
		Array<double> tx = Array<double>();
		Array<double> ty = Array<double>();
		Array<double> tz = Array<double>();

		force_storage->getTimeColumn(time,-1);
		//cout<<"\n\n"<<time;

		force_storage->getDataColumn("f1_vx",fx);
		force_storage->getDataColumn("f1_vy",fy);
		force_storage->getDataColumn("f1_vz",fz);
		force_storage->getDataColumn("t1_x",tx);
		force_storage->getDataColumn("t1_y",ty);
		force_storage->getDataColumn("t1_z",tz);
		//cout<<"\n\n"<<fx;

		

		vector<double> fx_v,fy_v,fz_v,tx_v,ty_v,tz_v;
		vector<double> time_v;

		for(int j = 0; j<force_storage->getSize(); j++)
		{
			time_v.push_back(time.get(j));
			fx_v.push_back(fx.get(j));
			fy_v.push_back(fy.get(j));
			fz_v.push_back(fz.get(j));
			tx_v.push_back(tx.get(j));
			ty_v.push_back(ty.get(j));
			tz_v.push_back(tz.get(j));
		} 

		double* t = &time_v[0];
		double* f1x = &fx_v[0];
		double* f1y = &fy_v[0];
		double* f1z = &fz_v[0];
		double* t1x = &tx_v[0];
		double* t1y = &ty_v[0];
		double* t1z = &tz_v[0];
		//cout<<"\n\n"<<*(t+10);
		//cout<<"\n\n"<<*(f1z+10);

		// Create function for all forces and torques
		OpenSim::Function* spline_fx = new OpenSim::GCVSpline(3,force_storage->getSize(),t,f1x,"fx",0);
		OpenSim::Function* spline_fy = new OpenSim::GCVSpline(3,force_storage->getSize(),t,f1y,"fy",0);
		OpenSim::Function* spline_fz = new OpenSim::GCVSpline(3,force_storage->getSize(),t,f1z,"fz",0);
		OpenSim::Function* spline_tx = new OpenSim::GCVSpline(3,force_storage->getSize(),t,t1x,"tx",0);
		OpenSim::Function* spline_ty = new OpenSim::GCVSpline(3,force_storage->getSize(),t,t1y,"ty",0);
		OpenSim::Function* spline_tz = new OpenSim::GCVSpline(3,force_storage->getSize(),t,t1z,"tz",0);
		/////////////////////////////////////

		OpenSim::CoordinateSet& CS = osimModel.updCoordinateSet();
		OpenSim::Function* function1 = new OpenSim::Constant(0);

		// Create pointer to BodySet inside osimModel
		BodySet& bodyset = osimModel.updBodySet();

		// Add external prescribed force to head_markers body using spline functions created. 
		OpenSim::Body * head_markers = &bodyset.get("head_markers");
		PrescribedForce *force1 = new PrescribedForce(head_markers); 

		// Prescribe functions for specifying the point that force is applied to (function 1 = constant value of 0)
		force1->setPointFunctions(function1,function1,function1);

		// Create linear function
		// Array<double> function_coeffs(0,2);
		// function_coeffs[0] = 5;
		// function_coeffs[1] = 0;
		// OpenSim::Function* function2 = new OpenSim::LinearFunction(function_coeffs);
		// Create constant function
		// OpenSim::Function* function3 = new OpenSim::Constant(8);

		// Prescribe force and torque functions using splines from force measurement file
		force1->setForceFunctions(spline_fx,spline_fy,spline_fz);
		force1->setTorqueFunctions(spline_tx,spline_ty,spline_tz);

		// Specify that force and point are NOT in global frame (i.e. expressed in frame of head_markers body)
		force1->setForceIsInGlobalFrame(0);
		force1->setPointIsInGlobalFrame(0);
		force1->setName("External_Load");
		// add force to osim model
		osimModel.addForce(force1); 



		//////////////////////////////////////////////////////////////////////////////////////////////
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
		//////////////////////////////////////////////////////////////////////////////////////////////

		//////////////////////////////////////////////////////////////////////////////////////////////
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
		//////////////////////////////////////////////////////////////////////////////////////////////

		// Set up point kinematic analyses -> to track marker positions for use in objective function

		//cout<<"\n\nINITAL ANALYSIS SET: "<<osimModel.getAnalysisSet()<<"\n\n";

		// Create list of point kinematics reporters which match marker positions from osim model
		PointKinematics* m1h = new PointKinematics(&osimModel);
		PointKinematics* m2h = new PointKinematics(&osimModel);
		PointKinematics* m3h = new PointKinematics(&osimModel);
		PointKinematics* m4h = new PointKinematics(&osimModel);
		PointKinematics* m5h = new PointKinematics(&osimModel);
		PointKinematics* m6h = new PointKinematics(&osimModel);
		
		// Set to return coordinates of point (specified wrt head_markers) but expressed in ground ref frame
		//m1h->setBody(&osimModel.updBodySet().get("head_markers"));
		//m2h->setBody(&osimModel.updBodySet().get("head_markers"));
		//m3h->setBody(&osimModel.updBodySet().get("head_markers"));
		//m4h->setBody(&osimModel.updBodySet().get("head_markers"));
		//m5h->setBody(&osimModel.updBodySet().get("head_markers"));
		//m6h->setBody(&osimModel.updBodySet().get("head_markers"));

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

		//cout<<"\n\nINITAL ANALYSIS SET: "<<osimModel.getAnalysisSet()<<"\n\n";

		//cout<<"\n\nMARKER SET: "<<osimModel.getMarkerSet()<<"\n\n";

		//////////////////////////////////////////////////////////////////////////////////////////////
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
		//////////////////////////////////////////////////////////////////////////////////////////////

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
		Vec3 t_stiff(0),r_stiff(0.01),t_damp(0),r_damp(10);
		
		// Open trc file and store in class MarkerData
		MarkerData md = OpenSim::MarkerData(marker_filename);

		// Create storage object with marker data (md) in it
		Storage data_trc;
		md.makeRdStorage(data_trc);
		// crop storage object to only include data in time frame of simulation
		data_trc.crop(ti,tf); // -> data_trc is used to calculate rms error between model marker and motion capture markers

		// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
		
		// Initialise and run optimization
		int numParams = 4;
		// Parameters:
		// - Theta_STAR
		// - k1		
		// - k2
		// - bushing offset
		
		MyOptimizerSystem sys(numParams,osimModel,data_trc,ti,tf,ICs,*m1h,*m2h,*m3h,*m4h,*m5h,*m6h,output_fd);

		Real f = NaN;

		Vector lower_bounds(numParams);
		Vector upper_bounds(numParams);

		// theta star
		lower_bounds[0] = 6;//0.0;//*Pi/180;
		upper_bounds[0] = 15;//1.0;//*Pi/180;

		// k1
		lower_bounds[1] = 1e-6;//0.0;//*Pi/180;
		upper_bounds[1] = 15;//1.0;//*Pi/180;

		// k2
		lower_bounds[2] = 0.1;
		upper_bounds[2] = 100;

		//// damping
		//lower_bounds[3] = 0.5;
		//upper_bounds[3] = 30;

		//bushing offset
		lower_bounds[3] = -2.0;
		upper_bounds[3] = 2.0;

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
	
		guess[0] = 6.82;//6.00553; // theta star (degrees)
		guess[1] = 1.22; // k1 (N*m/degree)
		guess[2] = 7.29; // k2 (N*m/degree)
		guess[3] = 0.22;//0.409055; // bushing offset
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
			opt.setConvergenceTolerance(1e-5);
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
		// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //	

		cout<<"\n\nOptimiszation time: "<<((float)t2-(float)t1)/ CLOCKS_PER_SEC << " seconds";
		
		//// Run simulation and save point kinematic reporter in data_sim storage object. Print simulation results to file.
		//Array<Array<double>> pk_data;

		//Vector PARAMS(numParams);
		//for (int i = 0; i<numParams; i++){
		//	PARAMS[i] = guess[i];
		//}

		//bool toWrite = 1;
		////pk_data = simtools->RunSimulation_LIMITSTOP(osimModel,PARAMS,ti,tf,ICs,toWrite,output_fd,*m1h,*m2h,*m3h,*m4h,*m5h,*m6h);
		//double rms = simtools->RunSimulation_wRMS(data_trc,osimModel,PARAMS,ti,tf,ICs,false,output_fd,*m1h,*m2h,*m3h,*m4h,*m5h,*m6h);
		//cout<<"\nRMS: "<<rms<<endl;
		//
		//////cout<<"\n\npk_data: "<<pk_data;
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
