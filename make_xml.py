#! /usr/bin/env python

"""
This script creates an input XML file for the 'ezSpectrum 3.0' programm

Script supports the following ab initio outputs:

Q-Chem 
ACES II (versions 0.3 [old] and 2.5.0 [new])
Molpro
GAMESS (linear molecules are not supported)
Gaussian (we can not guarantee it works, since we are not allowed to touch it)

Script requires a set of ab initio outputs with frequency jobs in
the optimized geometries. First file in the argument's list is
the initial state, others are the target states.
From each file this script copies a geometry, frequencies and
normal modes into the XML file.
Default job parameters and energy difference between the electronic
states (IP, EA, etc) should be adjusted as necessary in the created XML file.
"""

import sys


def main(xml_filename, ai_filenames, run_type):

    if len(ai_filenames)==0:
        if run_type!="web":
            print '\nTo create an XML file from ab initio outputs type:\n\
            "make_xml.py  <filename.xml> <initial_state.out> <target_state_1.out> <target_state_2.out> etc..."\n\n'
            sys.exit(2)
        else:
            return "Error. make_xml.py: no ab-initio file names found"

    # ##################################################
    # # Function ReadAndWriteState
    # ##################################################
    def ReadAndWriteState(FileName,which_state, run_type):
        """Reads one state from the qchem output file. Writes to the XML file"""

        #check which ab initio program created it
        #input_type="q-chem", "aces", or "molpro"
        input_type=""
        StateF=open(FileName,'r')
        Line=StateF.readline()
        while Line and (input_type==""):
            if (Line.find('Welcome to Q-Chem')>=0) and (input_type==""):
                input_type="q-chem"
            if (Line.find('* ACES2:')>=0) and (input_type==""):
                input_type="aces_old"
            if (Line.find('* ACES :')>=0) and (input_type==""):
                input_type="aces_new"
            if (Line.find('PROGRAM SYSTEM MOLPRO')>=0) and (input_type==""):
                input_type="molpro"
            if (Line.find('GAMESS VERSION = ')>=0) and (input_type==""):
                input_type="gamess"
            if (Line.find('RESTRICTED RIGHTS')>=0) and (input_type==""):
                input_type="other"

            Line=StateF.readline()
        if (input_type==""):
            if run_type!="web":
                print "Error: File \""+FileName+"\" has an unknown format. Q-Chem, Molpro, ACES II, GAMESS, or Gaussian frequency jobs are supported"
                sys.exit(2)
            else:
                return "Error: File \""+FileName+"\" has an unknown format. Q-Chem, Molpro, ACES II, GAMESS, or Gaussian frequency jobs are supported"
        
        StateF.close

        #read all the info from file:
        NAtoms=0
        NLinesWithFrequencies=0
        ifLinear = "false"
        Geometry=""
        atoms_list=""
        NormalModes=""
        Frequencies=""
        ifGeometryIsLoaded='false'
        ifFrequenciesLoaded='false'
        ifNormalModesLoaded='false'

        StateF=open(FileName,'r')


        if (input_type=="q-chem"):
            Line=StateF.readline()
            while Line:
                if (Line.find('Standard Nuclear Orientation')>=0) and (ifGeometryIsLoaded=='false'):
                    StateF.readline()
                    StateF.readline()
                    Line=StateF.readline()
                    while Line.find('----') == -1:
                        NAtoms+=1
                        Geometry+=Line[5:]
                        atoms_list+=Line[8:14]
                        Line=StateF.readline()
                    ifGeometryIsLoaded='true'

                if Line.find('Raman Active: ')>=0:
                    StateF.readline()
                    for i in range(NAtoms):
                        Line=StateF.readline()
                        NormalModes=NormalModes+Line[2:]
                    NormalModes+='\n'

                # remove end of the line symbols!!!
                if Line.find('Frequency: ') >= 0:
                    Frequencies+=Line.replace('Frequency:','')
                    NLinesWithFrequencies+=1

                Line=StateF.readline()

            if NLinesWithFrequencies==NAtoms-1:
                ifLinear="true"

            if_normal_modes_weighted="true"
            geometry_units="angstr"

            # END Q-CHEM
            # ================================================================================

               
        if (input_type=="aces_old"):
            Line=StateF.readline()
            while Line:
                # in FREQ job only one such line. Check it if it will be from OPT job....
                if (Line.find('Coordinates (in bohr)')>=0) and (ifGeometryIsLoaded=='false'):
                    StateF.readline()
                    StateF.readline()
                    Line=StateF.readline()
                    while Line.find('----') == -1:
                        NAtoms+=1
                        Geometry=Geometry+"      "+Line[5:9]+Line[20:]
                        atoms_list+=Line[5:9]
                        Line=StateF.readline()
                        ifGeometryIsLoaded='true'

                if Line.find('Normal Coordinates')>=0:
                    Line=StateF.readline()
                    while (Line=="\n"):
                        StateF.readline()
                        Line=StateF.readline()
                        NLinesWithFrequencies+=1
                        Frequencies+=Line
                        StateF.readline()
                        StateF.readline()
                        for i in range(NAtoms):
                            Line=StateF.readline()
                            NormalModes=NormalModes+Line[4:]
                        NormalModes+='\n'
                        Line=StateF.readline()
                Line=StateF.readline()
            if NLinesWithFrequencies==NAtoms-1:
                ifLinear="true"

            if_normal_modes_weighted="false"
            geometry_units="au"

           # END ACES_OLD
           # ================================================================================
           
        if (input_type=="aces_new"):
            Line=StateF.readline()
            while Line:
                # in FREQ job only one such line. Check it if it will be from OPT job....
                if (Line.find('C o o r d i n a t e s')>=0) and (ifGeometryIsLoaded=='false'):
                    StateF.readline()
                    StateF.readline()
                    Line=StateF.readline()
                    while Line.find('----') == -1:
                        NAtoms+=1
                        Geometry=Geometry+"      "+Line[5:9]+Line[20:]
                        atoms_list+=Line[5:9]
                        Line=StateF.readline()
                    ifGeometryIsLoaded='true'

                if Line.find('Normal Coordinates')>=0:
                    Line=StateF.readline()
                    Line=StateF.readline()
                    while (Line=="\n"):
                        StateF.readline()
                        Line=StateF.readline()
                        NLinesWithFrequencies+=1
                        Frequencies+=Line
                        StateF.readline()
                        StateF.readline()
                        for i in range(NAtoms):
                            Line=StateF.readline()
                            NormalModes=NormalModes+Line[4:]
                        NormalModes+='\n'
                        Line=StateF.readline()
                Line=StateF.readline()
            if NLinesWithFrequencies==NAtoms-1:
                ifLinear="true"

            if_normal_modes_weighted="false"
            geometry_units="au"

           # END ACES_NEW
           # ================================================================================

        if (input_type=="molpro"):
            normal_coordinates=[]
            n_normal_modes=0
        
            Line=StateF.readline()
            while Line:
                # geometry
                if (Line.find('ATOMIC COORDINATES')>=0) and (ifGeometryIsLoaded=='false'):
                    StateF.readline()
                    StateF.readline()
                    StateF.readline()
                    Line=StateF.readline()
                    while Line.find('Bond lengths in Bohr (Angstrom)') == -1:
                        NAtoms+=1
                        Geometry=Geometry+"      "+Line[5:8]+Line[19:]
                        atoms_list+=Line[5:8]
                        Line=StateF.readline()
                    ifGeometryIsLoaded='true'
                    NAtoms-=1
                    # create an empty list of coordinates for every atom in form "X_nm1 Y_nm1 Z_nm1 X_nm2 Y_nm2 Z_nm2 X_nm3..."
                    for i in range(NAtoms):
                        normal_coordinates.append([])
 
                if (Line.find('Intensities [relative]')>=0) and (ifNormalModesLoaded=='false'):
                    for i in range(NAtoms):
                        LineX=StateF.readline()
                        LineY=StateF.readline()
                        LineZ=StateF.readline()
                        
                        # analyse the length of the line 20+12*j:
                        nm_per_line= (len(LineX)-23)/12
                        
                        n_normal_modes+=nm_per_line
                        
                        for j in range(nm_per_line):
                            normal_coordinates[i].append( LineX[ 23+12*j : 23+12*(j+1) ]+LineY[ 23+12*j : 23+12*(j+1) ]+LineZ[ 23+12*j : 23+12*(j+1)] )

                    n_normal_modes/=NAtoms # repeated "NAtoms"-times
                    
                if (Line.find('Wavenumbers [cm-1]') >= 0) and (ifFrequenciesLoaded=='false'):
                    Frequencies+=Line.replace('Wavenumbers [cm-1]','')
                    NLinesWithFrequencies+=1

                if Line.find('Normal Modes of low/zero frequencies') >= 0:
                    ifFrequenciesLoaded='true'
                    ifNormalModesLoaded='true'

                Line=StateF.readline()

            # now create normal modes in q-chem format:
            printed_normal_modes=0
            for j in range( len(normal_coordinates[0])/3 ):
                for i in range(NAtoms):
                    NormalModes=NormalModes+normal_coordinates[i][j*3]+"   "+normal_coordinates[i][j*3+1]+"   "+normal_coordinates[i][j*3+2]+'\n'
                printed_normal_modes+=3
                NormalModes+='\n'

            # add the reminder of 3
            if ( (len(normal_coordinates[0]) % 3) > 0 ):
                for i in range(NAtoms):
                    for j in range( len(normal_coordinates[0]) % 3 ):
                        NormalModes=NormalModes+normal_coordinates[i][printed_normal_modes+j]+"   "
                    NormalModes+='\n'
                NormalModes+='\n'

            if (n_normal_modes==(3*NAtoms-5)):
                ifLinear="true"

            if_normal_modes_weighted="true"
            geometry_units="au"

            # END MOLPRO
            # ================================================================================

        if (input_type=="gamess"):
            
            if run_type!="web":
                print "\nWarning! All geometries from GAMESS' outputs treated as non-linear."
    
            ifLinear="false"

            normal_coordinates=[]
            n_normal_modes=0
        
            Line=StateF.readline()
            while Line:
                # geometry
                if (Line.find('COORDINATES (BOHR)')>=0) and (ifGeometryIsLoaded=='false'):
                    StateF.readline()
                    Line=StateF.readline()
                    while Line.find('INTERNUCLEAR DISTANCES (ANGS.)') == -1:
                        NAtoms+=1
                        Geometry=Geometry+"     "+Line[0:5]+Line[19:]
                        atoms_list+=Line[0:5]
                        Line=StateF.readline()
                    ifGeometryIsLoaded='true'
                    NAtoms-=1
                    # create an empty list of coordinates for every atom in form "X_nm1 Y_nm1 Z_nm1 X_nm2 Y_nm2 Z_nm2 X_nm3..."
                    for i in range(NAtoms):
                        normal_coordinates.append([])
 
                if (Line.find(' IR INTENSITY:')>=0) and (ifNormalModesLoaded=='false'):
                    Line=StateF.readline()
                   
                    for i in range(NAtoms):
                        LineX=StateF.readline()
                        LineY=StateF.readline()
                        LineZ=StateF.readline()
                        
                        # analyse the length of the line 20+12*j:
                        nm_per_line= (len(LineX)-20)/12
                        
                        for j in range(nm_per_line):
                            normal_coordinates[i].append( LineX[ 20+12*j : 20+12*(j+1) ]+LineY[ 20+12*j : 20+12*(j+1) ]+LineZ[ 20+12*j : 20+12*(j+1)] )

                    n_normal_modes/=NAtoms # repeated "NAtoms"-times

                if (Line.find('   FREQUENCY:') >= 0) and (ifFrequenciesLoaded=='false'):
                    Frequencies+=Line.replace('   FREQUENCY:','')
                    NLinesWithFrequencies+=1

                if Line.find('REFERENCE ON SAYVETZ CONDITIONS') >= 0:
                    ifFrequenciesLoaded='true'
                    ifNormalModesLoaded='true'

                Line=StateF.readline()

            # now create normal modes in q-chem format:
            printed_normal_modes=0
            for j in range( len(normal_coordinates[0])/3 ):
                for i in range(NAtoms):
                    NormalModes=NormalModes+normal_coordinates[i][j*3]+"   "+normal_coordinates[i][j*3+1]+"   "+normal_coordinates[i][j*3+2]+'\n'
                printed_normal_modes+=3
                NormalModes+='\n'

            # add the reminder of 3
            if ( (len(normal_coordinates[0]) % 3) > 0 ):
                for i in range(NAtoms):
                    for j in range( len(normal_coordinates[0]) % 3 ):
                        NormalModes=NormalModes+normal_coordinates[i][printed_normal_modes+j]+"   "
                    NormalModes+='\n'
                NormalModes+='\n'

            # remove first 6 normal modes (two lines of three normal modes)
            # i.e. rotation+translation (works only for non-linear molecules!):
            for i in range((NAtoms+1)*2):
                NormalModes=NormalModes[NormalModes.find('\n')+1:]
            # remove first 6 frequencies:
            Frequencies_set=Frequencies.split()
            Frequencies_set=Frequencies_set[6:]
            Frequencies = " ".join(["%s" % (f) for f in Frequencies_set])
            Frequencies = "       " + Frequencies + "\n"
            
            if_normal_modes_weighted="true"
            geometry_units="au"

            # END GAMESS
            # ================================================================================

        if (input_type=="other"):
            ifAtomsLoaded='false'
            
            Line=StateF.readline()
            while Line:
                if (Line.find('Distance matrix')>=0) and (ifAtomsLoaded=='false'):
                    StateF.readline()
                    while Line.find('Stoichiometry') == -1:
                        atoms_list+=Line[8:12]
                        Line=StateF.readline()
                    atoms_set=atoms_list.split()
                    ifAtomsLoaded='true'

                if (Line.find('Standard orientation')>=0) and (ifGeometryIsLoaded=='false'):
                    for L in range(4):
                        StateF.readline()
                    Line=StateF.readline()
                    while Line.find('----') == -1:
                        Geometry=Geometry+"      "+atoms_set[NAtoms]+Line[32:]
                        Line=StateF.readline()
                        NAtoms+=1
                    ifGeometryIsLoaded='true'

                if (Line.find('Atom AN')>=0) or (Line.find('Atom  AN')>=0):
                    for i in range(NAtoms):
                        Line=StateF.readline()
                        NormalModes=NormalModes+Line[10:]
                    NormalModes+='\n'

                # remove end of the line symbols!!!
                if Line.find('Frequencies --') >= 0:
                    Frequencies+=Line.replace('Frequencies --','')
                    NLinesWithFrequencies+=1

                Line=StateF.readline()

            if NLinesWithFrequencies==NAtoms-1:
                ifLinear="true"

            if_normal_modes_weighted="true"
            geometry_units="angstr"

            # END OTHER
            # ================================================================================

        # Write the state to the xml file
        xmlF.write('  <geometry\n'+'    number_of_atoms = "'+str(NAtoms)+'"\n')
        # linear?
        if ifLinear=='true':
            xmlF.write('    linear = "true"\n')
        else:
            xmlF.write('    linear = "false"\n')

        xmlF.write('    units   ="'+geometry_units+'"\n')
            
        xmlF.write('    text   = "\n')
        xmlF.write(Geometry)
        xmlF.write('             "\n        />\n\n')

        atoms_order = " ".join(["%s" % (nm) for nm in range(NAtoms)])

        xmlF.write('  <OPT_manual_atoms_reordering\n     new_order="'+atoms_order+'" />\n\n')

        xmlF.write('  <normal_modes\n    if_mass_weighted="'+if_normal_modes_weighted+'"\n    text = "\n')
        xmlF.write(NormalModes)
        xmlF.write('             "\n   atoms = "')
        xmlF.write(atoms_list)
        xmlF.write('"\n        />\n\n')

        if ifLinear=='true':
            normal_modes = " ".join(["%s" % (nm) for nm in range(3*NAtoms-5)])
        else:
            normal_modes = " ".join(["%s" % (nm) for nm in range(3*NAtoms-6)])
            
        if which_state=="target":
            if run_type!="web":
                xmlF.write('  <OPT_manual_normal_modes_reordering\n     new_order="' +normal_modes+ '" />\n\n')
            else: #web version:
                xmlF.write('  <manual_normal_modes_reordering\n     new_order="'+normal_modes+'" />\n\n')

        xmlF.write('  <frequencies\n    text = "\n')
        xmlF.write(Frequencies)
        xmlF.write('             "\n        />\n\n')

        return ""

###################################################
## END
###################################################

    
    xmlF=open(xml_filename,'w')

    # Write default job parameters
    xmlF.write("""<input
  job = "harmonic_pes">

<job_parameters 
        temperature                              = "300"
        spectrum_intensity_threshold             = "0.001"
/>

<!-- 
  ______________________________________________________________________ 

    Tags which start with "OPT_" will be ignored.
    To include these optional keywords please "uncomment" by removing 
    "OPT_" from the start and the corresponding end tag (if present)
  ______________________________________________________________________ 

 -->

<parallel_approximation
        max_vibr_excitations_in_initial_el_state = "1"
        max_vibr_excitations_in_target_el_state  = "4"
        combination_bands                        = "true"
        use_normal_coordinates_of_target_states  = "true"
 >

  <OPT_do_not_excite_subspace size = "0" normal_modes = " " />

  <OPT_energy_thresholds  units="eV, K, cm-1">
    <initial_state   units="K">      1000    </initial_state>
    <target_state   units="eV">      0.25    </target_state>
  </OPT_energy_thresholds> 

  <OPT_print_franck_condon_matrices flag="true"/>

</parallel_approximation>

<!-- 
  ______________________________________________________________________ 

 -->

<OPT_dushinsky_rotations target_state="1"
      max_vibr_excitations_in_initial_el_state = "1"
      max_vibr_excitations_in_target_el_state  = "4"
      > 
  <OPT_max_vibr_to_store  target_el_state="4"/>

  <OPT_do_not_excite_subspace size = "2" normal_modes = "0 1" />

  <OPT_energy_thresholds  units="eV, K, cm-1">
    <initial_state   units="K">      1000    </initial_state>
    <target_state   units="eV">      0.25    </target_state>
  </OPT_energy_thresholds> 

  <OPT_single_excitation 
       ini="0" 
       targ="1v1"/>

</OPT_dushinsky_rotations>

<!-- 
  ______________________________________________________________________ 

 -->\n\n""")
  
    # web version
    if run_type=="web": 
        xmlF.write('<if_web_version flag="true"/>\n')
 
    xmlF.write('<initial_state>\n')
    xmlF.write('  <!-- THIS INITIAL STATE IS FROM "'+ai_filenames[0]+'" FILE -->\n\n')

    return_string=ReadAndWriteState(ai_filenames[0],"initial", run_type)
    if not(return_string==""):
        return return_string
    
    xmlF.write('</initial_state>\n\n')
    xmlF.write("""<!-- 
  ______________________________________________________________________ 

 -->\n\n""")

    state_n=0
    for targetStateFileName in ai_filenames[1:]:
        state_n+=1
        xmlF.write('<target_state>\n\n')
        # xmlF.write('  <ip units="eV"> '+str(state_n)+' </ip>\n\n')
        xmlF.write('  <ip units="eV"> '+str(state_n)+' </ip>\n\n')
        xmlF.write('  <!-- THIS TARGET STATE IS FROM "'+targetStateFileName+'" FILE -->\n')
        return_string=ReadAndWriteState(targetStateFileName,"target", run_type)
        if not(return_string==""):
            return return_string
        xmlF.write('</target_state>\n\n')
        xmlF.write("""<!-- 
  ______________________________________________________________________ 

 -->\n\n""")
    xmlF.write('</input>\n')

    return ""

if __name__=="__main__":
    
    if len(sys.argv) < 2:
        print '\nTo create an XML file from ab initio outputs type:\n\
        "make_xml.py  <filename.xml> <initial_state.out> <target_state_1.out> <target_state_2.out> etc..."\n\n'
        sys.exit(2)
    main(sys.argv[1],sys.argv[2:],"command_line")
