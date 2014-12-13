#include "harmonic_pes_main.h"

bool harmonic_pes_main (const char* xmlFileName)
{
  simpleXMLparser xmlF;
  xmlF.assignFile(xmlFileName);

  //======= read "global" job variables  =====================================================
  xmlF.reset().node("input").node("job_parameters");

  // read
  double temperature;
  temperature=xmlF.getDoubleValue("temperature");

  // fcf threshold (from the <job_parameters> tag)
  double fcf_threshold=sqrt( xmlF.getDoubleValue("spectrum_intensity_threshold") );

  xmlF.reset().node("input");

  // check if print normal modes after transformations & overlap matrix
  bool if_print_normal_modes=false;
  if ( xmlF.CheckSubNode("print_normal_modes") )
    {
      xmlF.node("print_normal_modes");
      if_print_normal_modes = xmlF.getBoolValue("flag"); 
      xmlF.stepBack();
    }

  // check if the web version format of the output (do not print the input file & create a ".nmoverlap" file)
  bool if_web_version=false;
  if ( xmlF.CheckSubNode("if_web_version") )
    {
      xmlF.node("if_web_version");
      if_web_version = xmlF.getBoolValue("flag"); 
      xmlF.stepBack();
    }



  //===========================================================================================
  //Read initial state and N target states; i.e. (N+1) electronic states total
  std::vector <MolState> elStates;
  int state_i=-1;
  bool ifElStateExist=true;
  bool ifAnyNormalModesReordered=false;

  do {
    state_i++;
    
    MolState elSt(xmlFileName, ATOMIC_MASSES_FILE);
    
    if (state_i==0)
      elSt.reset().node("input").node("initial_state");
    else
      elSt.reset().node("input").node("target_state",state_i);
    
    if (elSt.Check()) 	
      {
	std::cout << '\n';
	if (state_i==0)
	  std::cout << "====== Reading the initial state ======\n";
	else
	  std::cout << "===== Reading the target state #" << state_i << " =====\n";
    	elSt.Read();
    	elStates.push_back(elSt);
	
	// index of the last state:
	int state=elStates.size()-1;

	if ( elStates[state].ifNMReorderedManually() )
	  {
	    ifAnyNormalModesReordered=true;
	    if ( state==0)
	      {
		std::cout<<"\nError! Manual reordering of the normal modes is not allowed for the initial state\n\n";
		exit(2);
	      }
	  }


	// check that states are similar (same number of atoms, same order of the atomic names, same "liniarity")
	if (state>0)
	  if ( not(elStates[state].ifSimilar(elStates[0])) )
	    {
	      std::cout << "Error: excited state #" << state << " is different from the initial state\n\n";
	      exit(2);
	    }


	// apply automatic transformations to the last loaded state (if no manual were requested):

	if ( not(elStates[state].ifAlignedManually()) )
	  {
	    // align each state: center of mass in the coordinates origin, moment of ineretia principal axes along the coordinate axes
	    elStates[state].align();

	    // align each target state with the initial one
	    if (state>0)
	      elStates[state].align(elStates[0]);

	    // get clean zeros
	    elStates[state].applyCoordinateThreshold(COORDINATE_THRESHOLD);
	  }

	std::cout << "\nNew molecular geometry:\n";
	elStates[state].printGeometry();
	// centerOfMass=elStates[i].getCenterOfMass().applyThreshold(COORDINATE_THRESHOLD).print("Center of mass: ");
	elStates[state].getMomentOfInertiaTensor().Print("\nMOI tensor:");

	if (if_print_normal_modes) 
	  {
	    std::cout << "Normal modes after the geometry transformations:\n\n";
	    elStates[state].printNormalModes();
	  }

      } // end: one more state in the input
    else
      ifElStateExist=false;
  } while (ifElStateExist);

  if (state_i==0) std::cout << "No initial (e.g. neutral) state info was found in \"" << xmlFileName << "\"\n";
  if (state_i==1) std::cout << "No target (e.g. cation) state(s) info was found in \"" << xmlFileName << "\"\n";
  if (state_i<=1) exit(-1);
  
  std::cout << '\n';
  std::cout << "Initial state and " << state_i-1 << " target state" << ((state_i==2)?" was":"s were") <<" loaded from \"" << xmlFileName << "\"\n\n";
  std::cout << "------------------------------------------------------------------------------\n";


  // total number of the normal modes (in the initial state)
  int n_norm_modes = elStates[0].NNormModes();

  // if parallel or dushinsky
  bool if_something_to_do=false;

  //======================================================================
  // Parallel approximation section

  xmlF.reset().node("input");
  if (xmlF.CheckSubNode("parallel_approximation"))
    {
      if_something_to_do=true;

      std::cout << "\n=== Reading the parallel approximation job parameters ===\n"<< std::flush;

      xmlF.reset().node("input").node("parallel_approximation");
      // Maximum number of vibrational levels to take:
      int max_n_initial, max_n_target;  //maximum number of vibrational levels for initial and target state
      max_n_initial = xmlF.getIntValue("max_vibr_excitations_in_initial_el_state"); // i.e. =2 <=> three vibr. states GS and two excited states
      max_n_target  = xmlF.getIntValue("max_vibr_excitations_in_target_el_state");
      bool if_comb_bands = xmlF.getBoolValue("combination_bands");
      bool if_use_target_nm = xmlF.getBoolValue("use_normal_coordinates_of_target_states");
      
      // check if print normal modes after transformations & overlap matrix
      bool if_print_fcfs=false;
      if ( xmlF.CheckSubNode("print_franck_condon_matrices") )
	{
	  xmlF.node("print_franck_condon_matrices");
	  if_print_fcfs = xmlF.getBoolValue("flag"); 
	  xmlF.stepBack();
	}
      
      // read energy thresholds (if provided)
      double energy_threshold_initial = DBL_MAX;//eV
      double energy_threshold_target = DBL_MAX; //eV
      if ( xmlF.CheckSubNode("energy_thresholds") )
	{
	  xmlF.node("energy_thresholds");
	  if ( xmlF.CheckSubNode("initial_state") )
	    {
	      xmlF.node("initial_state");
	      std::string units=xmlF.value("units");
	      energy_threshold_initial=xmlF.getDoubleValue();
	      if (units=="cm-1")
		energy_threshold_initial*=WAVENUMBERS2EV;
	      else if (units=="K")
		energy_threshold_initial*=KELVINS2EV;
	      else if (units!="eV")
		{
		  std::cout << "\nError! Unknow units of the initial state threshold: \"" << units <<"\"\n  (should be equal to \"eV\", \"K\", or \"cm-1\")\n\n";
		  exit(1);
		}
	      xmlF.stepBack();
	    }
	  if ( xmlF.CheckSubNode("target_state") )
	    {
	      xmlF.node("target_state");
	      std::string units=xmlF.value("units");
	      energy_threshold_target=xmlF.getDoubleValue();
	      if (units=="cm-1")
		energy_threshold_target*=WAVENUMBERS2EV;
	      else if (units=="K")
		energy_threshold_target*=KELVINS2EV;
	      else if (units!="eV")
		{
		  std::cout << "\nError! Unknow units of the target state threshold: \"" << units <<"\"\n  (should be equal to \"eV\", \"K\", or \"cm-1\")\n\n";
		  exit(1);
		}
	      xmlF.stepBack();
	    }
	  xmlF.stepBack();
	}
      
      // read normal modes do_not_excite subspace (for parallel approximation only):
      std::set<int> do_not_excite_subspace;
      bool if_use_do_not_excite_subspace=false;
      int do_not_excite_subspace_size=0;
      int do_not_excite_subspace_max=0; // maximum value in lists -- error check later (dirty)
      std::istringstream tmp_iStr;
      int tmpInt;
      xmlF.reset().node("input").node("parallel_approximation");
      if ( xmlF.CheckSubNode("do_not_excite_subspace") )
	{
	  if_use_do_not_excite_subspace=true;
	  xmlF.node("do_not_excite_subspace");
	  tmp_iStr.clear();
	  tmp_iStr.str(xmlF.value("normal_modes")); 
	  do_not_excite_subspace_size=xmlF.getIntValue("size");
	  for (int nm=0; nm<do_not_excite_subspace_size; nm++)
	    {
	      tmp_iStr >> tmpInt;
	      //input error check:
	      if (tmp_iStr.fail())
		{
		  std::cout << "\nFormat error in \"input\"->\"do_not_excite_subspace\"->\"normal_modes\"\n"
			    << "(non numeric symbol or less entries then specified by the \"size\" value)\n\n";
		  xmlF.exitOnFormatError(true);
		}
	      if (tmpInt<0) 
		{
		  std::cout << "\nFormat error in \"input\"->\"do_not_excite_subspace\"->\"normal_modes\"\n"
			    << "Entry ["<< tmpInt<<"] is negative.\n\n";
		  xmlF.exitOnFormatError(true);
		}
	      // keep the maximum value of the list
	      if (tmpInt>do_not_excite_subspace_max) do_not_excite_subspace_max=tmpInt;
	      
	      //check if tmpInt is already in the set:
	      std::set<int>::const_iterator intSet_iter;
	      intSet_iter = do_not_excite_subspace.find(tmpInt);
	      if ( intSet_iter != do_not_excite_subspace.end( ) )
		{
		  std::cout << "\nFormat error in \"input\"->\"do_not_excite_subspace\"->\"normal_modes\"\n"
			    << "Entry ["<< tmpInt<<"] is not unique.\n\n";
		  xmlF.exitOnFormatError(true);
		}
	      do_not_excite_subspace.insert(tmpInt);
	    }

	  if (do_not_excite_subspace.size()!=0)
	    {

	      if(ifAnyNormalModesReordered)
		std::cout <<"WARNING! The normal modes of the target state were reordered!\n"
			  <<"         New order is used for the \"do_not_excite_subspace\".\n\n";
	      
	      std::cout << "The following normal modes will have no vibrational excitations:\n";

	      for (std::set<int>::const_iterator intSet_iter=do_not_excite_subspace.begin(); intSet_iter!=do_not_excite_subspace.end(); intSet_iter++)
		std::cout << *intSet_iter << ' ';
	      std::cout<<"\n";
	      xmlF.stepBack();
	    }

	}

      //check that numbers in do_not_excite_subspace are less than the number_of_normal_modes
      if ((do_not_excite_subspace_max >= n_norm_modes)and(if_use_do_not_excite_subspace==true))
	{
	  std::cout << "\nError! Maximum normal mode number in \"do_not_excite_subspace\" is ["<< do_not_excite_subspace_max<<"],\n"
		    << "  which is greater than (number_of_normal_modes-1)="<< n_norm_modes-1 <<"\n\n";
	  exit(2);
	}
      
      //create nms_parallel -- "excite subspace" (full_space-do_not_excite_subspace)
      std::vector<int> nms_parallel;
      for (int nm=0; nm<n_norm_modes; nm++)
	{
	  std::set<int>::const_iterator intSet_iter;
	  intSet_iter = do_not_excite_subspace.find(nm);
	  if ( intSet_iter == do_not_excite_subspace.end( ) )
	    nms_parallel.push_back(nm);
	}
      
      //================================================================================
      // print the overlap matrix with the initial state for each target states:
      
      for (int state=1; state<elStates.size();state++)
	{
	  std::cout << "\n===== Overlap matrix of the target state #" << state << " with the initial state =====\n";
	  
	  std::vector <int> normal_modes_list;
	  KMatrix NMoverlap;  //normal modes overlap matric (for each target state the same matrix is used)
	  bool if_overlap_diagonal;
	  
	  
	  // select nondiagonal submatrix of the overlap matrix:
	  if_overlap_diagonal=elStates[state].getNormalModeOverlapWithOtherState(elStates[0], NMoverlap, normal_modes_list);
	  // rows -- norm modes of the target state; colums norm modes of the initial state;

	  // remove normal modes from normal_modes_list that are in the do_not_excite_subspace:
	  std::set<int>::iterator iter_set;
	  std::vector<int> new_normal_modes_list; 
	  
	  for (int nm=0; nm<normal_modes_list.size(); nm++)
	    {
	      // if nm is not in the do_not_excite set:
	      iter_set = do_not_excite_subspace.find(normal_modes_list[nm]);
	      if ( iter_set == do_not_excite_subspace.end( ) )
		//than copy it to the new list:
		new_normal_modes_list.push_back(normal_modes_list[nm]);
	    }
	  
	  //Create an overlap submatrix:
	  if ((if_overlap_diagonal) or (new_normal_modes_list.size()<=1))
	    {
	      std::cout << "The normal modes overlap matrix whit the initial state is diagonal\n";
	      if (new_normal_modes_list.size()<=1)
		std::cout<<"  (do_not_excite_subspace is excluded)\n";
	      std::cout<<"\n";
	    }
	  else
	    {        
	      std::cout << "WARNING! The normal mmodes overlap matrix whith the initial state\n"
			<< "         is non-diagonal! Consider normal modes reordering.\n\n";
	      // create a normal mode submatrix:
	      KMatrix overlap_submatrix(new_normal_modes_list.size(),new_normal_modes_list.size());
	      overlap_submatrix.Set(0.0);
	      for (int nm1=0; nm1<new_normal_modes_list.size(); nm1++)
		for (int nm2=0; nm2<new_normal_modes_list.size(); nm2++)
		  overlap_submatrix.Elem2(nm1,nm2)=NMoverlap.Elem2(new_normal_modes_list[nm1],new_normal_modes_list[nm2]);
	      
	      //print the overlap_submatrix (with correct column/row labbels):
	      std::cout << "  The non diagonal part of the normal modes overlap matrix (do_not_excite_subspace is excluded):";
	      std::cout << "\n     ";
	      
	      for (int j=0; j<new_normal_modes_list.size(); j++)
		std::cout << std::fixed << std::setprecision(0) << std::setw(8) << new_normal_modes_list[j];
	      for (int i=0; i<new_normal_modes_list.size(); i++)
		{
		  std::cout << "\n  "<< std::fixed << std::setprecision(0) << std::setw(3) << new_normal_modes_list[i];
		  for (int j=0; j<new_normal_modes_list.size(); j++)
		    if (fabs(overlap_submatrix.Elem2(i,j)) >= 0.001)
		      std::cout << std::fixed << std::setprecision(3) << std::setw(8) << overlap_submatrix.Elem2(i,j);
		    else
		      std::cout << "      --";
		}
	      std::cout <<"\n\n";
	    }

	  // print in a "fit 80 chars wide termial" form
	  if(if_print_normal_modes)
	    NMoverlap.Print("Normal modes overlap matrix with the initial state \n(if significantly non diagonal, please consider normal modes reordering)");
	}
      
      // for the web version: save the overlap matrix (with displacements) in an xml file
      std::stringstream nmoverlapFName; 
      nmoverlapFName << xmlFileName << ".nmoverlap";
      

      std::cout << "------------------------------------------------------------------------------\n\n";
      std::cout << "Photoelectron spectrum in the parallel approximation will be evaluated\n\n"<< std::flush;

      //================================================================================
      //================================================================================
      //================================================================================
      //================================================================================
      
      // create a new parallel approximation object (evaluates and stores FCFs in the harmonic approximation)
      Parallel* parallel_ptr = new Parallel(elStates, nms_parallel, 
					    fcf_threshold, temperature, 
					    max_n_initial, max_n_target, 
					    if_comb_bands, if_use_target_nm, if_print_fcfs, if_web_version, nmoverlapFName.str().c_str(),  
					    energy_threshold_initial,  energy_threshold_target);
      
      //================================================================================
      //================================================================================
      //================================================================================
      //================================================================================
      
      
      //--------------------------------------------------------------------------------
      // Print the updated spectrum:
      (*parallel_ptr).getSpectrum().Sort();
      std::cout << "------------------------------------------------------------------------------\n";
      std::cout << "           Stick photoelectron spectrum (parallel approximation)\n";
      std::cout << "------------------------------------------------------------------------------\n";
      if (ifAnyNormalModesReordered)
	std::cout <<"\nWARNING! The normal modes of one of the target states were reordered!\n"
		  <<"         New order is used for the target state assignment.\n";
      if( nms_parallel.size()!=n_norm_modes )
	{
	  std::cout << "\nNOTE: only the following normal modes were excited: (\"excite subspace\"):\n  ";
	  for (int nm=0; nm< nms_parallel.size(); nm++)
	    std::cout << nms_parallel[nm] << ' ';
	  std::cout << "\n";
	  if (ifAnyNormalModesReordered)
	    std::cout <<"\nWARNING! The normal modes of one of the target states were reordered!\n"
		      <<"         New order is used for the \"excite subspace\"\n";
 	    
	}
      std::cout << "\n";
      if ( (*parallel_ptr).getSpectrum().getNSpectralPoints()>0)
	(*parallel_ptr).getSpectrum().PrintStickTable();
      else                
	std::cout << "\n\n\n"
		  <<"WARNING! The spectrum is empty.\n\n"
		  <<"         Plese refer to \"The spectrum is empty!\" in the\n"
		  <<"         \"Common problems\" section of the manual\n\n\n\n";

      std::cout << "------------------------------------------------------------------------------\n";
  
      // save this spectrum to the file
      std::stringstream spectrumFName; 
      spectrumFName << xmlFileName << ".spectrum_parallel";
      (*parallel_ptr).getSpectrum().PrintStickTable(spectrumFName.str().c_str());
      std::cout << "\nStick spectrum was also saved in \"" << spectrumFName.str() << "\" file \n";
      if(if_use_do_not_excite_subspace)
	std::cout << " (Full list of the normal modes was used for assining transitions)\n";

      std::cout << "------------------------------------------------------------------------------\n\n";
  
      delete parallel_ptr;
    }



  //=========================================================================
  // Dushinski rotation (nonparallel approximation)
  //=========================================================================
  // the notation and equations are from [Berger et al. JPCA 102:7157(1998)]
  //

  bool ifDushinsky=false;
  xmlF.reset().node("input");
  if ( xmlF.CheckSubNode("dushinsky_rotations") )
    ifDushinsky=true;
 
  if (ifDushinsky)
    {
      if_something_to_do=true;

      std::cout << "\n\n=== Reading the Dushinsky rotations job parameters ===\n\n"<< std::flush;

      //----------------------------------------------------------------------
      // load parameters 

      // indexes of the initial and target electronic states:
      xmlF.reset().node("input").node("dushinsky_rotations");
      int iniN=0;
      int targN=xmlF.getIntValue("target_state"); 
      if ((targN>elStates.size()-1)or(targN<1))
	{
	  std::cout<<"\nFormat error: \"target_state\" value should be positive and not greater than "<<elStates.size()-1<< "\n\n";
	  xmlF.exitOnFormatError(true);
	}
      std::cout<<"Target state number "<<targN<<" from the input will be used\n\n";
	  
      // maximum number of quanta to store:     
      int max_quanta_ini = xmlF.getIntValue("max_vibr_excitations_in_initial_el_state");
      int max_quanta_targ = xmlF.getIntValue("max_vibr_excitations_in_target_el_state");

      // read "do not excite subspace"
      std::set<int> do_not_excite_subspace;
      std::istringstream tmp_iStr;
      int tmpInt;

      int Kp_max_to_save=32000;
      if ( xmlF.CheckSubNode("max_vibr_to_store") )
	{
	  xmlF.node("max_vibr_to_store");
	  Kp_max_to_save=xmlF.getIntValue("target_el_state");
      	  xmlF.stepBack();
	}

      if ( xmlF.CheckSubNode("do_not_excite_subspace") )
	{
	  xmlF.node("do_not_excite_subspace");

	  int ss_size=xmlF.getIntValue("size");
	  if ((ss_size<0) or (ss_size>n_norm_modes))
	    {
	      std::cout << "\nError: subspace size is out of range\n\n";
	      xmlF.exitOnFormatError(true);
	    }
	  if (ss_size>0)
	    {

	      if (elStates[targN].ifNMReorderedManually())
	      std::cout <<"WARNING! The normal modes of the target state were reordered!\n"
			<<"         New order is used for the \"do_not_excite_subspace\".\n\n";

	      tmp_iStr.clear();
	      tmp_iStr.str(xmlF.value("normal_modes")); 
	      for (int nm=0; nm<ss_size; nm++)
		{
		  tmp_iStr >> tmpInt;
		  //input error check:
		  if (tmp_iStr.fail())
		    {
		      std::cout << "\nError: non numeric symbol or less entries then specified by the \"size\" value\n\n";
		      xmlF.exitOnFormatError(true);
		    }
		  if ((tmpInt<0)or(tmpInt>n_norm_modes-1))
		    {
		      std::cout << "\nError: Entry ["<< tmpInt<<"] is out of range.\n\n";
		      xmlF.exitOnFormatError(true);
		    }

		  //check if tmpInt is already in the set:
		  std::set<int>::const_iterator intSet_iter;
		  intSet_iter = do_not_excite_subspace.find(tmpInt);
		  if ( intSet_iter != do_not_excite_subspace.end( ) )
		    {
		      std::cout << "\nEntry ["<< tmpInt<<"] is not unique.\n\n";
		      xmlF.exitOnFormatError(true);
		    }
		  do_not_excite_subspace.insert(tmpInt);
		}
	    }
      	  xmlF.stepBack();
	}

      //----------------------------------------------------------------------
      // create the "excite subspace"; only normal modes from this subspace will be excited
      std::vector<int> nms_dushinsky;
      for (int nm=0; nm<n_norm_modes;nm++)
	{
	  std::set<int>::const_iterator intSet_iter;
	  intSet_iter = do_not_excite_subspace.find(nm);
	  if ( intSet_iter == do_not_excite_subspace.end( ) )
	    nms_dushinsky.push_back(nm);
	}
      
      if (nms_dushinsky.size()<n_norm_modes)
	{
	  std::cout << "The following normal will be excited\n(for both states the order is same as in in the input):\n ";
	  for (int nm=0; nm<nms_dushinsky.size(); nm++)
	    std::cout << nms_dushinsky[nm] << ' ';
	  std::cout<<"\n\n";
	}
      else
	std::cout << "All normal modes modes will be excited\n\n";



      std::cout << "=== Photoelectron spectrum with the Dushinsky rotations will be evaluated ==== \n\n"<< std::flush;

      //================================================================================
      //================================================================================
      //================================================================================
      //================================================================================

      // create a new dushinsky object for a given set of normal modes
      // All matrices and zero-zero integral are evaluated for the full space;
      // than excitations are only added to the normal modes from "nms_dushinsky"
      Dushinsky* dushinsky_ptr = new Dushinsky(elStates, nms_dushinsky, fcf_threshold, targN, max_quanta_targ, max_quanta_ini);

      //================================================================================
      //================================================================================
      //================================================================================
      //================================================================================

      // print estimated size of each layer up to K'_max:
      std::cout << "Number of normal modes to excite: " << nms_dushinsky.size() << "\n\n";
      std::cout << "Size of layers with exactly K' excitations in the target state (in bytes):\n";
      (*dushinsky_ptr).printLayersSizes( (max_quanta_targ<Kp_max_to_save)? max_quanta_targ : Kp_max_to_save );
      std::cout << "\n";
      //--------------------------------------------------------------------------------
      
      // go over all layers and add points to the spectrum:
      for (int Kp=1; Kp<=max_quanta_targ; Kp++)
	{
	  if ( Kp<=Kp_max_to_save )
	    std:: cout << "Layer K'="<< Kp <<" is being evaluated (will be saved in memory)... " << std::flush;
	  else
	    std:: cout << "Layer K'="<< Kp <<" is being evaluated... " << std::flush;

	  int n_fresh_points=(*dushinsky_ptr).evalNextLayer( Kp<=Kp_max_to_save );

	  std:: cout << "Done\n";
	  if (n_fresh_points>0)            
	    std:: cout << n_fresh_points <<" points above the intensity threhold were added to the spectrum\n\n" << std::flush;
	  else
	    std:: cout << "No points above the intensity threhold were found in this layer\n\n" << std::flush;
	}

      //----------------------------------------------------------------------
      // add the hot bands if requested

      double energy_threshold_initial = DBL_MAX;//eV
      double energy_threshold_target = DBL_MAX; //eV
      if (max_quanta_ini!=0)
	{
	  // read the energy thresholds (if provided)
	  if ( xmlF.CheckSubNode("energy_thresholds") )
	    {
	      xmlF.node("energy_thresholds");
	      if ( xmlF.CheckSubNode("initial_state") )
		{
		  xmlF.node("initial_state");
		  std::string units=xmlF.value("units");
		  energy_threshold_initial=xmlF.getDoubleValue();
		  if (units=="cm-1")
		    energy_threshold_initial*=WAVENUMBERS2EV;
		  else if (units=="K")
		    energy_threshold_initial*=KELVINS2EV;
		  else if (units!="eV")
		    {
		      std::cout << "\nError! Unknow units of the initial state threshold: \"" << units <<"\"\n  (should be equal to \"eV\", \"K\", or \"cm-1\")\n\n";
		      exit(1);
		    }
		  xmlF.stepBack();
		}
	      if ( xmlF.CheckSubNode("target_state") )
		{
		  xmlF.node("target_state");
		  std::string units=xmlF.value("units");
		  energy_threshold_target=xmlF.getDoubleValue();
		  if (units=="cm-1")
		    energy_threshold_target*=WAVENUMBERS2EV;
		  else if (units=="K")
		    energy_threshold_target*=KELVINS2EV;
		  else if (units!="eV")
		    {
		      std::cout << "\nError! Unknow units of the target state threshold: \"" << units <<"\"\n  (should be equal to \"eV\", \"K\", or \"cm-1\")\n\n";
		      exit(1);
		    }
		  xmlF.stepBack();
		}
	      xmlF.stepBack();
	    }
	  
	  
	  int n_hot_bands=(*dushinsky_ptr).addHotBands(elStates, nms_dushinsky, fcf_threshold, temperature, 
						       max_quanta_ini, max_quanta_targ, 
						       energy_threshold_initial,energy_threshold_target);
	  std:: cout << n_hot_bands <<" hot bands were added to the spectrum\n"
		     <<"Note: the Boltzmann distribution will be applied later\n\n" << std::flush;
	}

      //----------------------------------------------------------------------
      // now load the list of single transitions to evaluate FCFs recursively
      // single_transition is in the "full space"; do_not_excite_subspace does not applies;
      SpectralPoint single_transition;
      for (int nm=0; nm<elStates[iniN].NNormModes(); nm++)
	{
	  single_transition.getVibrState1().addVibrQuanta(0, nm);
	  single_transition.getVibrState2().addVibrQuanta(0, nm);
	}
      single_transition.getVibrState1().setElStateIndex(iniN);
      single_transition.getVibrState2().setElStateIndex(targN);

      
      int ex_n=1;
      xmlF.node("single_excitation", ex_n);
      if ( xmlF.Check() )
	{
	  if (elStates[targN].ifNMReorderedManually())
	    std::cout <<"WARNING! The normal modes of the target state were reordered!\n"
		      <<"         New order is used for the single transitions.\n\n";

	  std::cout << "The following single transitions were added to the spectrum:\n" << std::flush;
	}

      while ( xmlF.Check() )
	{
	  My_istringstream ini_str(xmlF.value("ini"));
	  fillVibrState(ini_str, single_transition.getVibrState1(), xmlF, n_norm_modes);
  
	  My_istringstream targ_str(xmlF.value("targ"));
	  fillVibrState(targ_str, single_transition.getVibrState2(), xmlF, n_norm_modes);

	  // evaluate FCF for each transition and add to the spectrum:
	  int K =single_transition.getVibrState1().getTotalQuantaCount();
	  int Kp=single_transition.getVibrState2().getTotalQuantaCount();
	  double s_fcf=(*dushinsky_ptr).evalSingleFCF_full_space(single_transition.getVibrState1(), K, single_transition.getVibrState2(),Kp);
	  (*dushinsky_ptr).addSpectralPoint(s_fcf, single_transition.getVibrState1(), single_transition.getVibrState2()); 

	  std::cout << "FCF=" << std::scientific << std::setprecision(6) << s_fcf << " ";
	  single_transition.getVibrState1().print();
	  std::cout << "->";
	  single_transition.getVibrState2().print();
	  std::cout << "\n" << std::flush;

	  ex_n++;
	  xmlF.stepBack();
	  xmlF.node("single_excitation", ex_n);
 	}
                     
      std:: cout << "\nUpdating the energies and applying the Boltzmann distribution..." << std::flush;
      //--------------------------------------------------------------------------------
      //update (fill) energies for every point in the spectrum and add the Boltzmann distribution:
      int points_removed=0;
      for (int pt=0; pt<(*dushinsky_ptr).getSpectrum().getNSpectralPoints(); pt++)
	{
	  double energy = -elStates[targN].Energy();
	  double E_prime_prime = 0; // no hot bands

	  // run it over the full space, if nm not in the nms_dushinsky subspace, getV_full_dim() returns zero (no excitations):
	  for (int nm=0; nm <elStates[iniN].NNormModes(); nm++)
	    {
	      energy += elStates[iniN].getNormMode(nm).getFreq() * WAVENUMBERS2EV * (*dushinsky_ptr).getSpectrum().getSpectralPoint(pt).getVibrState1().getV_full_dim(nm);
	      energy -= elStates[targN].getNormMode(nm).getFreq() * WAVENUMBERS2EV * (*dushinsky_ptr).getSpectrum().getSpectralPoint(pt).getVibrState2().getV_full_dim(nm);

	      E_prime_prime += 
		elStates[iniN].getNormMode(nm).getFreq() * WAVENUMBERS2EV * (*dushinsky_ptr).getSpectrum().getSpectralPoint(pt).getVibrState1().getV_full_dim(nm);
	    }

	  // add the Boltzmann distribution to the initial state population
	  double IExponent;
	  if (temperature==0)
	    if (E_prime_prime==0)
	      IExponent=0;   // intensity unchanged 
	    else
	      IExponent=100; //(intensity < 10e-44 for non ground states
	  else
	    {
	      IExponent= E_prime_prime / (temperature * KELVINS2EV );
	      if (IExponent > 100) 
		IExponent=100; // keep the intensity >= 10e-44 == exp(-100)
	    }
	  (*dushinsky_ptr).getSpectrum().getSpectralPoint(pt).getIntensity() *= exp ( -IExponent); 
	  

	  (*dushinsky_ptr).getSpectrum().getSpectralPoint(pt).getEnergy()=energy;
	  (*dushinsky_ptr).getSpectrum().getSpectralPoint(pt).getE_prime_prime()=E_prime_prime;

	  // if intensity below the fcf_threshold^2 or energy above the threshold -- do not print
	  if (
	      ( (*dushinsky_ptr).getSpectrum().getSpectralPoint(pt).getIntensity() < fcf_threshold*fcf_threshold ) 
	      or 
	      ( -(energy-E_prime_prime+elStates[targN].Energy()) > energy_threshold_target ) 
	      or 
	      ( E_prime_prime > energy_threshold_initial )
	      )
	    {
	      (*dushinsky_ptr).getSpectrum().getSpectralPoint(pt).setIfPrint(false);
	      points_removed++;
	    }
	}
      std:: cout << "Done\n" << std::flush;

      if (max_quanta_ini!=0)
	{
	  if (points_removed>0)                 
	    std::cout << "  " << points_removed << " hotbands were removed from the spectrum\n";
	  else
	    std::cout << "All hotbands are above the intensity threshold\n";
	  std::cout << "\n" << std::flush;
	      
	}

      //--------------------------------------------------------------------------------
      // Print the updated spectrum:
      (*dushinsky_ptr).getSpectrum().Sort();
      std::cout << "------------------------------------------------------------------------------\n";
      std::cout << "        Stick photoelectron spectrum (with Dushinsky rotations) \n";
      std::cout << "------------------------------------------------------------------------------\n";
      if (elStates[targN].ifNMReorderedManually())
	std::cout <<"\nWARNING! The normal modes of the target state were reordered!\n"
		  <<"         New order is used for the target state assignment.\n";
      if(nms_dushinsky.size()!=n_norm_modes)
	{
	  std::cout << "\nNOTE: only the following normal modes were excited (\"excite subspace\"):\n  ";
	  for (int nm=0; nm<nms_dushinsky.size();nm++)
	    std::cout <<nms_dushinsky[nm]<< ' ';
	  std::cout << "\n";
	  if (elStates[targN].ifNMReorderedManually())
	    std::cout <<"\nWARNING! The normal modes of the target state were reordered!\n"
		      <<"         New order is used for the \"excite subspace\"\n";
 	}
      std::cout << "\n";
      (*dushinsky_ptr).getSpectrum().PrintStickTable();
      std::cout << "------------------------------------------------------------------------------\n";
  
      // save the spectrum to the file
      std::stringstream spectrumFName; 
      spectrumFName << xmlFileName << ".spectrum_dushinsky";
      (*dushinsky_ptr).getSpectrum().PrintStickTable(spectrumFName.str().c_str());
      std::cout << "\nStick spectrum was also saved in \"" << spectrumFName.str() << "\" file \n";
      if(nms_dushinsky.size()!=n_norm_modes)
	std::cout << " (Full list of the normal modes was used for assining transitions)\n";
      std::cout << "\n\n";

      delete dushinsky_ptr;
    }

  if (!if_something_to_do)
    {                 
      std::cout << "\nError! No \"parallel_approximation\" or \"dushinsky_rotations\" section\n       was found in the input. Nothing to do.\n\n";
      exit(2);
    }
    

  return true;
};









//! splits string of type "3v21" into two integers 3 and 21
void get_qnt_nm(std::string& ex_str, simpleXMLparser& xmlF, int& qnt, int& nm )
{
  if (ex_str.find("v")==std::string::npos)
    {
      std::cout << "\nFormat error in [" << ex_str << "] excitation: should contain symbol \'v\'\n\n";
      xmlF.exitOnFormatError(true);
    }
  ex_str.replace( ex_str.find("v"), 1, " " );
  std::istringstream ex_strs(ex_str);
  ex_strs>>qnt>>nm;

  if (ex_strs.fail())
    {
      ex_str.replace( ex_str.find(" "), 1, "v" );
      std::cout << "\nFormat error in [" << ex_str << "] excitation. Should be two integers separated by the symbol 'v' \n\n";
      xmlF.exitOnFormatError(true);
    }
}






//! converts string of type "1v1,1v2,1v3,3v19" into a vibrational state (i.e. vector of integers)
void fillVibrState(My_istringstream& vibr_str, VibronicState& v_state, simpleXMLparser& xmlF, const int nm_max )
{
  // quanta & normal mode number (for parcing strings like "3v19", where qnt=3 and nm=19)
  int qnt=0, nm=0;
  // string like 3v19"
  std::string ex_str;
  // string like "1v1,1v2,1v3,3v19"
  vibr_str.getNextWord(ex_str);
  // reset vibrational state
  for (int i=0; i<v_state.getVibrQuantaSize(); i++)
    v_state.setVibrQuanta(i,0);
  // fill vibrational state (if ==0 -- nothing to do)
  if (ex_str!="0")
    {
      get_qnt_nm(ex_str, xmlF, qnt, nm );

      if (nm>nm_max)
	{
	  std::cout << "\nError! Normal mode "<< nm <<" (in ["<< qnt<<'v'<<nm <<"] excitation) is out of range.\n\n";
	  xmlF.exitOnFormatError(true);
	}

      v_state.setVibrQuanta(nm,qnt);
      while (not(vibr_str.fail()))
	{
	  vibr_str.getNextWord(ex_str);
	  get_qnt_nm(ex_str, xmlF, qnt, nm );
	  v_state.setVibrQuanta(nm,qnt);
	}
    }
}


