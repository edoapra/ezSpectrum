#include "dushinsky.h"
Dushinsky::Dushinsky(std::vector <MolState>& molStates, std::vector<int>& nm_list, double in_fcf_threshold, int in_targN, int max_quanta_target, int max_quanta_initial)
{
  //--------------------------------------------------------------------------------
  // creates the Dushinsky object. Constractor evaluates <0|0> integral and all required matrices for a given set of normal modes.
  // i.e. the layer L=0 and is ready for iterative evaluation of the next layers

  // subspace defined by the nm_list (list of normal modes to include)
  // molStates -- vector of molecular stated loaded from the input file;
  // in_targN -- which target state to take from the vector;

  //IMPORTANT! All matrices and zero-zero integral are evaluated for the FULL space of the normal modes of the molecule
  //than space is shrinked to space described by nm_list, so excitations will be added only to those normal modes.

  Kp_max=0;
  Kp_max_saved=0;
  fcf_threshold=in_fcf_threshold;

  // index of the inital electronic state
  int iniN = 0;
  targN=in_targN;

  // number of normal modes in the subset
  N = molStates[iniN].NNormModes();
  
  //----------------------------------------------------------------------
  // get frequencies, displacements and rotation matrix for the current subspace:

  // Load frequencies:
  std::vector<double> selected_nm_freq_ini;
  selected_nm_freq_ini.clear();
  for (int nm=0; nm<N; nm++)
    selected_nm_freq_ini.push_back( molStates[iniN].getNormMode(nm).getFreq() );
  std::vector<double> selected_nm_freq_targ;
  selected_nm_freq_targ.clear();
  for (int nm=0; nm<N; nm++)
    selected_nm_freq_targ.push_back( molStates[targN].getNormMode(nm).getFreq() );

  // get w & w'; frequencies; p=='==target electronic state; [N]; frequency==w_i=nu_i[cm-1]*BOHR_2_CM*2*PI*SPEEDOFLIGHT_AU; 
  KMatrix w(N,1), wp(N,1);
  for (int nm=0; nm<N; nm++)
    w[nm]=selected_nm_freq_ini[nm]*BOHR_2_CM*2*PI*SPEEDOFLIGHT_AU;
  //w.Print("w");
  for (int nm=0; nm<N; nm++)
    wp[nm]=selected_nm_freq_targ[nm]*BOHR_2_CM*2*PI*SPEEDOFLIGHT_AU;
  //wp.Print("wp");

  // get normal modes: L and Lp (for the selected subspace):
  KMatrix L( CARTDIM*(molStates[iniN].NAtoms()), N );
  L.Set(0.0);  
  for (int a=0; a < molStates[iniN].NAtoms(); a++) 
    for (int nm = 0; nm < N; nm++) 
      for (int k=0; k < CARTDIM; k++)
	L.Elem2(a*CARTDIM+k, nm) = molStates[iniN].getNormMode(nm).getDisplacement()[a*CARTDIM+k];
  KMatrix Lp( CARTDIM*(molStates[targN].NAtoms()), N );
  Lp.Set(0.0);  
  for (int a=0; a < molStates[targN].NAtoms(); a++) 
    for (int nm = 0; nm < N; nm++) 
      for (int k=0; k < CARTDIM; k++)
	Lp.Elem2(a*CARTDIM+k, nm) = molStates[targN].getNormMode(nm).getDisplacement()[a*CARTDIM+k];

  // get S; [NxN]; S=(Lp^T)*L; n.m.rotation matrix; if S==I, than norm. modes are parallel, det(S) sould be close to 1
  KMatrix S;
  Lp.Transpose();
  S=Lp;
  S*=L;
  //S.Print("S");

  // get d displacements ; [N]; d= Lp^T*aqrt(T)*(x-xp); x & xp -- cartesian geometries of to states;
  // displacement of the target geometry in normal coordinates  d[angstr*sqrt(amu)] *ANGSTROM2AU*sqrt(AMU_2_ELECTRONMASS) ); 
  // mass weighted -- so m=1amu everywhere;
  KMatrix d;
  KMatrix tmpX(CARTDIM*molStates[iniN].NAtoms(),1), tmpXp(CARTDIM*molStates[iniN].NAtoms(),1);
  for(int a=0; a<molStates[iniN].NAtoms(); a++)
    for (int k=0; k <CARTDIM; k++)
      {
	tmpX[a*CARTDIM+k] =molStates[iniN].getAtom(a).Coord(k);
	tmpXp[a*CARTDIM+k]=molStates[targN].getAtom(a).Coord(k);
      }
  //(x-xp)
  tmpX-=tmpXp;
  //d= Lp^T*aqrt(T)*(x-xp); Lp is already transposed;
  d=Lp;
  //d.Print("Lp^T");
  //Make squrt(T)-matrix (diagonal matrix with sqrt(atomic masses) in cartesian coordinates); the matrix is the same for all states;
  KMatrix SqrtT( CARTDIM*(molStates[iniN].NAtoms()), CARTDIM*(molStates[iniN].NAtoms()), true);
  SqrtT.Set(0.0);
  for(int a=0; a<molStates[iniN].NAtoms(); a++)
    SqrtT.Elem2(a*CARTDIM,a*CARTDIM) = SqrtT.Elem2(a*CARTDIM+1,a*CARTDIM+1)= SqrtT.Elem2(a*CARTDIM+2,a*CARTDIM+2)=sqrt(molStates[iniN].getAtom(a).Mass());
  d*=SqrtT;
  d*=tmpX;
  //d.Print("d tmp");
  // switch to atomic units
  d*=ANGSTROM2AU*sqrt(AMU_2_ELECTRONMASS);

  //----------------------------------------------------------------------
  // now evaluate all the matrices:  

  // temporary NxN matrices:
  KMatrix tmpM(N,N), tmpM2(N,N);

  // diag(1) NxN matrix:
  KMatrix I(N,N);
  I.Set();
  I.SetDiagonal(1,true);

  // get Lm & Lmp; lamda & lamda'; [NxN] diag.; sqrt.freq.; \lamda=diag(sqrt(w_i)); 
  KMatrix Lm(N,N), Lmp(N,N);
  Lm.Set(0);
  Lmp.Set(0);
  for (int nm=0; nm<N; nm++)
    Lm.Elem2(nm,nm)=sqrt(w[nm]);
  //Lm.Print("diag. Lm");
  for (int nm=0; nm<N; nm++)
    Lmp.Elem2(nm,nm)=sqrt(wp[nm]);
  //Lmp.Print("diag. Lm\'");
  
  // get Dt; [N]; \delta=Lmp*d;
  KMatrix Dt(N,1);
  tmpM=Lmp;
  tmpM*=d;
  Dt=tmpM;
  //Dt.Print("Dt=Lmp*d");
      
  // get J;  [NxN]; J=Lmp*S*Lm^{-1};  ^{-1} -- inverse;
  KMatrix J(N,N);
  J=Lmp;
  J*=S;
  tmpM=Lm;
  tmpM.Inverse();
  J*=tmpM;
  //J.Print("J=Lmp*S*Lm^{-1}");

  // get Q;  [NxN] symm. pos.; Q = (1 + J^T * J)^{-1}; ^{T} -- transposed {KMatrix::Transpose()};  ^{-1} -- inverse;
  KMatrix Q(N,N);
  tmpM=J;
  tmpM.Transpose();
  tmpM*=J;
  tmpM+=I;
  Q=tmpM.Inverse();
  //Q.Print("Q = (1 + J^T * J)^{-1}; symm. pos. ");
  double detQ=Q.Determinant();
      
  // get P;  [NxN] symm.; P = J * Q * J^T
  KMatrix P(N,N);
  P=J;
  P*=Q;
  tmpM=J;
  tmpM.Transpose();
  P*=tmpM;
  //P.Print("P = J * Q * J^T; symm.");

  // get R;  [NxN]; R = Q * J^T
  KMatrix R(N,N);
  tmpM=J;
  tmpM.Transpose();
  R=Q;
  R*=tmpM;
  //R.Print("R = Q * J^T");


  // get Det(S)
  double detS=S.Determinant();
  std::cout << "Determinant of the normal modes rotation matrix: |Det(S)| =" << std::fixed << std::setw(12) << std::setprecision(8) << fabs(detS)<<"\n";
  if (fabs(detS)<0.5)
    {   
      std::cout << "\nError: |Det(S)| is too small (<0.5). Please see \"|Det(S)| is less then one\n   (or even zero)\" in the \"Common problems\" section of the manual\n\n";
      exit(2);
    }
  

  //--------------------------------------------------------------------------------
  // zero_zero, K'=0: one element = <0|0>
  double zero_zero=1;
  for (int nm=0; nm<N; nm++)
    zero_zero*=w[nm]/wp[nm];
  // ORIGINAL: pow 0.25; MODIFIED TO: pow -0.25
  zero_zero=pow(2.0,N*0.5)*pow(zero_zero,-0.25)*sqrt(detQ)/sqrt(fabs(detS));
  tmpM=Dt;
  tmpM.Transpose();
  tmpM2=I;
  tmpM2-=P;
  tmpM*=tmpM2;
  tmpM*=Dt;
  tmpM*=-0.5;
  zero_zero*=exp(tmpM[0]);

  //--------------------------------------------------------------------------------
  // get the the frequently used vectors/matrices:

  // for target state excitations (no hot bands):

  // ompd=sqrt(2)*(1-P)*Dt; [N]
  ompd_full=I;
  ompd_full-=P;
  ompd_full*=Dt;
  ompd_full*=sqrt(2);
  //ompd_full.Print("sqrt(2)*(1-P)*Dt");
  
  // tpmo=(2P-1); [NxN]
  tpmo_full=P;
  tpmo_full*=2;
  tpmo_full-=I;
  //tpmo_full.Print("2P-1");


  // for the hot bands:

  // rd=sqrt(2)*R*Dt; [N]
  rd_full = R;
  rd_full *= Dt;
  rd_full *= sqrt(2);
  //rd_full.Print("rd=sqrt(2)*R*Dt");


  // tqmo=(2Q-1); [NxN]
  tqmo_full=Q;
  tqmo_full*=2;
  tqmo_full-=I;
  //tqmo_full.Print("tqmo=(2Q-1)");

  // tr=2*R; [NxN]
  tr_full=R;
  tr_full*=2;
  //tr_full.Print("tr=2*R");

  
  //----------------------------------------------------------------------
  // shrink the space to the nm_list space

  if (N!=nm_list.size())
    {
      N = nm_list.size();
      
      // shrink "frequently used matrices"
      ompd.Adjust(N,1);
      tpmo.Adjust(N,N);
      rd.Adjust(N,1);
      tqmo.Adjust(N,N);
      tr.Adjust(N,N);

      for (int nm1=0; nm1<N; nm1++)
	{
	  ompd[nm1]=ompd_full[ nm_list[nm1] ];
	  rd[nm1]=rd[ nm_list[nm1] ];
	  for (int nm2=0; nm2<N; nm2++)
	    {
	      tpmo.Elem2(nm1,nm2)=tpmo_full.Elem2(nm_list[nm1],nm_list[nm2]);
	      tqmo.Elem2(nm1,nm2)=tqmo_full.Elem2(nm_list[nm1],nm_list[nm2]);
	      tr.Elem2(nm1,nm2)=tr_full.Elem2(nm_list[nm1],nm_list[nm2]);
	    }
	}
    }
  else
    {
      ompd=ompd_full;
      tpmo=tpmo_full;
      rd=rd_full;
      tqmo=tqmo_full;
      tr=tr_full;
    }
 

  //----------------------------------------------------------------------
  // create the initial ground vibr. state <0000...000|
  state0.clear();
  state0.setElStateIndex(iniN);
  for (int nm=0; nm<nm_list.size(); nm++)
    state0.addVibrQuanta(0, nm_list[nm]);

  // create a target ground vibr. state |0000...000>
  state.clear();
  state.setElStateIndex(targN);
  for (int nm=0; nm<nm_list.size(); nm++)
    state.addVibrQuanta(0, nm_list[nm]);


  //----------------------------------------------------------------------
  //add the zero-th layer, K=0:
  std::vector<double>* layer_ptr= new std::vector<double>;
  (*layer_ptr).push_back(zero_zero);
  layers.push_back(layer_ptr);
  if (fabs(zero_zero) > fcf_threshold)
    // set energies later (for now 0.0):
    addSpectralPoint(zero_zero, state0, state  );


  // fill matrix C with combinations C(n,k)=C_n^k=n!/(k!*(n-k)!); 
  // n=0..N+K-1; k=0..N-1=0..min(n,N-1); matrix (N+K)x(N) -- N is known everywhere;
  // Combination(n,k)=C(n,k)=C[n*N+k]
  // N-total number of normal modes; K -- max number of quanta in target state
  
  int K = max_quanta_target;
  int size=(N)*(N+K);
  
  /* // ZZZ 4/11/2012 removed, and the combinations are calculated now on the fly
  C=new unsigned long[size];
  memset(C,0,sizeof(unsigned long)*size);
  for (int n=0; n<N+K; n++)
    for (int k=0; k<MIN(n+1,N); k++)
	C[n*N+k]=Combination(n,k);
  */

  // Create an array of sqrt() 0..K+1
  K = MAX(max_quanta_target, max_quanta_initial);
  sqrtArray=new double[K+2];
  for (int i=0; i<K+1; i++)
    sqrtArray[i]=sqrt(i);
  
}


Dushinsky::~Dushinsky()
{
  for (int Kp=0; Kp<layers.size(); Kp++)
    delete layers[Kp];
  layers.clear();
  
  // ZZZ 4/11/2012 removed, and the combinations are calculated now on the fly
  // delete [] C;
  delete [] sqrtArray;
}




int Dushinsky::evalNextLayer(const bool if_save)
{
  // debug check:
  if (if_save and (Kp_max)!=(Kp_max_saved))
    {
      std::cerr << "\nDebug Error!: Layer " << Kp_max+1 << " can not be saved. The number of already saved layers is " 
		<< Kp_max_saved << " but should be " << Kp_max << "\n\n";
      exit(2);
    }

  // points above the threshold
  int points_added=0;

  // reset the target state
  for (int i=0; i<N; i++)
    state.getV()[i]=-1;
  
  //create a new layer
  std::vector<double>* layer_ptr= new std::vector<double>;

  //check the memory avaliability (dirty):
  if (if_save)
    {
      unsigned long total_combs = Combination(  (Kp_max+1 + state.getVibrQuantaSize() - 1), (state.getVibrQuantaSize() - 1)  );
      double * buffer = (double*) malloc (total_combs*sizeof(double));
      if (buffer==NULL)
	{
	  std::cout << "\n\nError: not enough memory available to store layer K'="<<Kp_max_saved+1<<"\n\n"
		    << "Please reduce \"max_vibr_to_store\" value to limit the memory use.\n"
		    <<" You also can reduce \"max_vibr_excitations_in_target_el_state\" value,\n"
		    << "and/or use \"single_excitation\" elements to add higher excitations.\n";
	  exit(2);
	}
      free (buffer);
    }
  
  unsigned long index_counter=0;
  while (enumerateVibrStates(state.getVibrQuantaSize(), Kp_max+1, state.getV(), true))
    {

      // ZZZ 4/11/2012 removed, and the combinations are calculated now on the fly
      // unsigned long index_rev=convVibrState2Index(state.getV(), N, C, Kp_max+1);
      unsigned long index_rev=convVibrState2Index(state.getV(), N, Kp_max+1);
      
      // check if the reverse index is ok
      if (index_rev!=index_counter) 
	{
	  // if numbers are too large, factorials in Combination() fanction will be out of "unsigned long" range ...
	  std::cout << "\n Error!\n[Debug info: reverse index function convVibrState2Index(state.getV()) for the state:\n";
	  state.print();
	  std::cout <<"\n"<<"returns index=" << index_rev<< "; should be index="<< index_counter<<"]\n\n"
		    << "Layer #" << Kp_max << " is the maximum layer which can be handled for this system.\n"
		    << "Please reduce \"max_vibr_excitations_in_target_el_state\" value to "<<Kp_max<<".\n" 
		    << "You can also use \"single_excitation\" elements to add higher excitations manually.\n";
	  exit(2);
	}
      

      double fcf=evalSingleFCF(state0, 0, state, Kp_max+1);

      if (if_save)
	(*layer_ptr).push_back(fcf);

      if (fabs(fcf) > fcf_threshold)
	{
	  points_added++;
	  addSpectralPoint(fcf,state0, state );
	} 

      index_counter++;
    }

  Kp_max++;
  // save the layer
  if (if_save)
    {
      layers.push_back(layer_ptr);
      Kp_max_saved++;
    }

  return points_added;
}





double Dushinsky::evalSingleFCF(VibronicState& state_ini, int K, VibronicState& state_targ, int Kp)
{
  double fcf;
 
  if (K==0) // if there are no excitations in the initial state left

    if (Kp<=Kp_max_saved) // if there are less than K'max excitations in the target state left (i.e. saved)
      // ZZZ 4/11/2012 removed, and the combinations are calculated now on the fly
      // fcf=(*layers[Kp])[convVibrState2Index(state_targ.getV(), N, C, Kp)];
      fcf=(*layers[Kp])[convVibrState2Index(state_targ.getV(), N, Kp)];
  
    else
      {
	// find first non zero quanta
	int ksi=0;
	while (state_targ.getV()[ksi]==0)
	  ksi++;
	
	// get the first term
	state_targ.getV()[ksi]--;
	fcf = ompd[ksi]*evalSingleFCF(state_ini, K, state_targ, Kp-1);
	
	// add the second term for K'>=2
	if (Kp>1)
	  for (int theta=ksi; theta<N; theta++)
	    if (state_targ.getV()[theta]>0)
	      {
		double tmp_dbl = tpmo.Elem2(ksi,theta) * sqrtArray[ state_targ.getV()[theta] ];
		state_targ.getV()[theta]--;
		tmp_dbl*= evalSingleFCF(state_ini, K, state_targ, Kp-2);
		fcf += tmp_dbl;
		state_targ.getV()[theta]++;
	      }
	fcf/=sqrtArray[ state_targ.getV()[ksi]+1 ];
	
	// get back to the original state ("to be calculated" state)
	state_targ.getV()[ksi]++;
      }
  else 
    {
      // find first non zero quanta
      int ksi=0;
      while (state_ini.getV()[ksi]==0)
	ksi++;
      
      // get the first term
      state_ini.getV()[ksi]--;
      fcf = -rd[ksi]*evalSingleFCF(state_ini, K-1, state_targ, Kp);

      // add the second term
      if (K>1)
	for (int theta=ksi; theta<N; theta++)
	  if (state_ini.getV()[theta]>0)
	    {
	      double tmp_dbl = tqmo.Elem2(ksi,theta) * sqrtArray[state_ini.getV()[theta] ];
	      state_ini.getV()[theta]--;
	      tmp_dbl*= evalSingleFCF(state_ini, K-2, state_targ, Kp);
	      fcf += tmp_dbl;
	      state_ini.getV()[theta]++;
	    }

      // add the third term
      if (Kp>0)
	for (int theta=0; theta<N; theta++)
	  if (state_targ.getV()[theta]>0)
	    {
	      double tmp_dbl = tr.Elem2(ksi,theta) * sqrtArray[ state_targ.getV()[theta] ];
	      state_targ.getV()[theta]--;
	      tmp_dbl*= evalSingleFCF(state_ini, K-1, state_targ, Kp-1);
	      fcf += tmp_dbl;
	      state_targ.getV()[theta]++;
	    }
      
      // normolize
      fcf/=sqrtArray[ state_ini.getV()[ksi]+1 ];
     
      // get back to the original state ("to be calculated" state)
      state_ini.getV()[ksi]++;
    }

  return fcf;
}




double Dushinsky::evalSingleFCF_full_space(VibronicState& state_ini, int K, VibronicState& state_targ, int Kp)
{
  double fcf;
 
  if (K==0) // if there are no excitations in the initial state left
    if (Kp==0) // if there are no excitations in the target state left
      fcf=(*layers[0])[0]; // return <0|0> integral
    else
      {
	// find first non zero quanta
	int ksi=0;
	while (state_targ.getV()[ksi]==0)
	  ksi++;
	
	// get the first term
	state_targ.getV()[ksi]--;
	fcf = ompd_full[ksi]*evalSingleFCF_full_space(state_ini, K, state_targ, Kp-1);
	
	// add the second term for K'>=2
	if (Kp>1)
	  for (int theta=ksi; theta<N; theta++)
	    if (state_targ.getV()[theta]>0)
	      {
		double tmp_dbl = tpmo_full.Elem2(ksi,theta) * sqrt( state_targ.getV()[theta] );
		state_targ.getV()[theta]--;
		tmp_dbl*= evalSingleFCF_full_space(state_ini, K, state_targ, Kp-2);
		fcf += tmp_dbl;
		state_targ.getV()[theta]++;
	      }
	fcf/=sqrt( state_targ.getV()[ksi]+1 );
	
	// get back to the original state ("to be calculated" state)
	state_targ.getV()[ksi]++;
      }
  else 
    {
      // find first non zero quanta
      int ksi=0;
      while (state_ini.getV()[ksi]==0)
	ksi++;
      
      // get the first term
      state_ini.getV()[ksi]--;
      fcf = -rd_full[ksi]*evalSingleFCF_full_space(state_ini, K-1, state_targ, Kp);

      // add the second term
      if (K>1)
	for (int theta=ksi; theta<N; theta++)
	  if (state_ini.getV()[theta]>0)
	    {
	      double tmp_dbl = tqmo_full.Elem2(ksi,theta) * sqrt( state_ini.getV()[theta] );
	      state_ini.getV()[theta]--;
	      tmp_dbl*= evalSingleFCF_full_space(state_ini, K-2, state_targ, Kp);
	      fcf += tmp_dbl;
	      state_ini.getV()[theta]++;
	    }

      // add the third term
      if (Kp>0)
	for (int theta=0; theta<N; theta++)
	  if (state_targ.getV()[theta]>0)
	    {
	      double tmp_dbl = tr_full.Elem2(ksi,theta) * sqrt( state_targ.getV()[theta] );
	      state_targ.getV()[theta]--;
	      tmp_dbl*= evalSingleFCF_full_space(state_ini, K-1, state_targ, Kp-1);
	      fcf += tmp_dbl;
	      state_targ.getV()[theta]++;
	    }
      fcf/=sqrt( state_ini.getV()[ksi]+1 );
     
      // get back to the original state ("to be calculated" state)
      state_ini.getV()[ksi]++;
    }

  return fcf;
}



int Dushinsky::addHotBands(std::vector <MolState>& molStates, std::vector<int>& nm_list, 
			    double fcf_threshold, double temperature, 
			    int max_n_initial, int max_n_target, 
			    double energy_threshold_initial,  double energy_threshold_target)
{

  int points_added=0;

  // vector of states below the energy thresholds and with total number of excitations up to requested number
  std::vector < VibronicState > selected_states_ini, selected_states_targ;

  // reset the initial state
  for (int i=0; i<N; i++)
    state0.getV()[i]=-1;
  // reset the target state
  for (int i=0; i<N; i++)
    state.getV()[i]=-1;

  unsigned long total_combs_ini=0, total_combs_targ=0;
  for (int curr_max_ini=0; curr_max_ini<=max_n_initial; curr_max_ini++)
    total_combs_ini += Combination(  (curr_max_ini + nm_list.size() - 1), (nm_list.size() - 1)  );
  for (int curr_max_targ=0; curr_max_targ<=max_n_target; curr_max_targ++)
    total_combs_targ += Combination(  (curr_max_targ + nm_list.size() -1), (nm_list.size() - 1)  );
  std::cout << "Maximum number of combination bands = " << total_combs_ini*total_combs_targ  
	    << "\n   = " << total_combs_ini << " (# of vibrational states in the initial electronic state)"
	    << "\n   * " << total_combs_targ  << " (# of vibrational states in the target electronic state)\n\n" << std::flush;

  // find INITIAL states with up to 'max_n_initial' vibrational quanta and with energy below 'energy_threshold_initial':
  std::cout << "A set of initial vibrational states is being created...\n";
  if (energy_threshold_initial < DBL_MAX)    
    std::cout << "  energy threshold = " << std::fixed << energy_threshold_initial <<" eV ("<< energy_threshold_initial/KELVINS2EV <<" K)\n" << std::flush;
  else
    {
      std::cout << "  energy threshold is not specified in the input (Please consider adding\n"
		<< "  the \"energy_thresholds\" tag for a faster calculation.)\n\n";
    }

  // std::cout <<"NOTE: ezSpectrum may crash at this point if memory is insufficient to store\n"
  //	    <<"      all vibrational states. If so, please reduce the initial state's energy\n"
  //	    <<"      threshold or max_vibr_excitations_in_initial_el_state\n\n" << std::flush;

  // start with one quanta in the initial state (hot bands)
  for (int curr_max_ini=1; curr_max_ini<=max_n_initial; curr_max_ini++)
    {
      while ( enumerateVibrStates(state0.getVibrQuantaSize(), curr_max_ini, state0.getV(), true) )
	{
	  double energy = 0;

	  for (int nm=0; nm <molStates[state0.getElStateIndex()].NNormModes(); nm++)
	    energy += 
	      molStates[state0.getElStateIndex()].getNormMode(nm).getFreq() * WAVENUMBERS2EV * state0.getV_full_dim(nm);

	  if (energy < energy_threshold_initial)
	    {
	      // check memory available (dirty, but anyway one copying of the state is requared...)
	      //VibronicState * state_tmp = (VibronicState*) malloc ( sizeof(VibronicState)*2 );
	      //if (state_tmp==NULL)
	      //{
	      //  std::cout << "\nError: not enough memory available to store all initial vibrational states\n\n";
	      //  exit(2);
	      //}
	      //free (state_tmp);

	      selected_states_ini.push_back(state0);
	    }
	}
	  //reset the initial state's vibration "population"
      state0.getV()[0]=-1;
    }
  std::cout << "  " << selected_states_ini.size() << " vibrational states below the energy threshold\n\n"<<std::flush;


  // find TARGET states with up to 'max_n_target' vibrational quanta and with energy below 'energy_threshold_target':
  std::cout << "A set of target vibrational states is being created...\n";
  if (energy_threshold_target < DBL_MAX)
    std::cout << "  energy threshold = " << std::fixed <<energy_threshold_target <<" eV ("<< energy_threshold_target/WAVENUMBERS2EV <<" cm-1)\n"<<std::flush;
  else
    std::cout << "  energy threshold is not specified in the input (Please consider adding\n"
	      << "  the \"energy_thresholds\" tag for a faster calculation.)\n\n";

  //std::cout <<"NOTE: ezSpectrum may crash at this point if memory is insufficient to store\n"
  //	    <<"      all vibrational states. If so, please reduce the target state's energy\n"
  //	    <<"      threshold or max_vibr_excitations_in_target_el_state\n\n" << std::flush;

  for (int curr_max_targ=0; curr_max_targ<=max_n_target; curr_max_targ++)
    {
      while ( enumerateVibrStates(state.getVibrQuantaSize(), curr_max_targ, state.getV(), true) )
	{
	  double energy = 0;

	  for (int nm=0; nm < molStates[state.getElStateIndex()].NNormModes(); nm++)
	    energy += molStates[state.getElStateIndex()].getNormMode(nm).getFreq() * WAVENUMBERS2EV * state.getV_full_dim(nm);

	  if ( energy < energy_threshold_target )
	    {
	      // check memory available (dirty, but anyway one copying of the state is requared...)
	      //VibronicState * state_tmp = (VibronicState*) malloc ( sizeof(VibronicState)*2 );
	      //if (state_tmp==NULL)
	      //	{
	      //	  std::cout << "\nError: not enough memory available to store all target vibrational states\n\n";
	      //	  exit(2);
	      //	}
	      //free (state_tmp);

	      selected_states_targ.push_back(state);
	    }
	}
      //reset the target state's vibration "population"
      state.getV()[0]=-1;
    }
  std::cout << "  " << selected_states_targ.size() << " vibrational states below the energy threshold\n\n";

  std::cout << "Total number of combination bands with thresholds applied: " << selected_states_ini.size()*selected_states_targ.size() <<"\n\n"<<std::flush;

  std::cout << "Hot bands are being calculated..." <<std::flush;


  // evaluate all possible cross integrals 
  for (int curr_ini=0; curr_ini<selected_states_ini.size(); curr_ini++)
    for (int curr_targ=0; curr_targ<selected_states_targ.size(); curr_targ++)
      {
	int K = selected_states_ini[curr_ini].getTotalQuantaCount();
	int Kp= selected_states_targ[curr_targ].getTotalQuantaCount();

	double s_fcf=evalSingleFCF(selected_states_ini[curr_ini], K, selected_states_targ[curr_targ], Kp);

	if (fabs(s_fcf) > fcf_threshold)
	  {
	    points_added++;
	    addSpectralPoint(s_fcf,selected_states_ini[curr_ini], selected_states_targ[curr_targ] );
	  } 
      }

  std::cout << "Done\n\n" << std::flush;

  return points_added;
}




void  Dushinsky::addSpectralPoint(const double fcf, VibronicState state_ini, VibronicState state_targ )
{
  // set energies later
  spectrum.AddSpectralPoint( 0.0, fcf*fcf, fcf, 0.0, state_ini, state_targ);
};


void Dushinsky::printLayersSizes(const int uptoKp)
{
  unsigned long elements_per_layer, size_per_layer=0, size_per_layer_prev=0;

  for (int Kp=0; Kp<=uptoKp; Kp++)
    {
      elements_per_layer = Combination(Kp+N-1, N-1 );
      size_per_layer_prev = size_per_layer;
      size_per_layer = elements_per_layer * sizeof(double);

      if (size_per_layer<size_per_layer_prev)
	{  
	  std::cout<< "\nError! Layer #" << Kp-1 << " is the maximum layer which can be saved "
		   <<"for this molecule. Please use \"max_vibr_to_store\" or/and"
		   <<"\"max_vibr_excitations_in_target_el_state\" value to "<<Kp-1<<".\n";
	  exit(2);
	}

      std::cout << "   layer K'=" << Kp << ": " ;
      std::cout.precision(2);
      if (size_per_layer<1024) 
	std::cout << size_per_layer;
      else if (size_per_layer<1024*1024)
	std::cout << size_per_layer/1024.0 << " K";
      else if (size_per_layer<1024*1024*1024)
	std::cout << size_per_layer/(1024.0*1024.0) << " M";
      else 
	std::cout << size_per_layer/(1024.0*1024.0*1024.0) << " G";
      std::cout << "\n";
    }          

  std::cout << "Please be sure that you have enough memory to store all this layers,\n" 
	    << "otherwise use \"max_vibr_to_store\" tag to limit the memory use.\n"
	    << "You also can reduce \"max_vibr_excitations_in_target_el_state\" value,\n"
	    << "and/or use \"single_excitation\" elements to add higher excitations.\n";
}

