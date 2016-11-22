#Simulation script for a bead-spring polymer with attached (charged) side beads
#usage: script.tcl [simulation_identifier] [simulation_max_time]
# 0 simulation_identifier: id for the current simulation run, will be part of the name of the generated output files, simulation configuration will be loaded accordingly, also important for checkpointing
# 1 simulation_max_time: simulation time after which the simulation is terminated, if the given time starts with a + this is the time that will be added to the simulation time


package require Tclx

# Filenames & paths
#############################################################

#simulation identifier here set to microtime but could also be passed as argument - easy checkpointing possible
set simulation_identifier [lindex $argv 0]

if { $simulation_identifier=="" } {
    set simulation_identifier "default_[clock microseconds]"
    set input_simulation_config_file "default_simulation_config.txt"
} else {
    set input_simulation_config_file "[lindex $argv 0]_simulation_config.txt"
}

puts "Simulation Identifier: $simulation_identifier"


# set simulation time (needed here for filenames)
set simulation_max_time_default 1000

set simulation_max_time [lindex $argv 1]

set add_sim_time 0

if { ![string is double $simulation_max_time] || $simulation_max_time == "" } {
    set simulation_max_time $simulation_max_time_default
    puts "No value for maximum simulation time found. Set 'simulation_max_time' to '$simulation_max_time'"
} else {
    if { [string index $simulation_max_time 0] == "+" } {
        set add_sim_time 1
        puts "Simulation will be run for additional '$simulation_max_time' time units."
    } else {
        puts "Set 'simulation_max_time' to '$simulation_max_time'"
    }
}


set host [info hostname]

puts "Hostname: $host"
#t_random seed 20160414
puts "t_random seed: [t_random seed]"

# Define paths containing no simulation_max_time
#############################################################

set output_vtf_path  ./output_vtf/
set output_vtf_file_warmup "${simulation_identifier}_polymer_sidechains_warmup.vtf"
set output_vtf_file_equi "${simulation_identifier}_polymer_sidechains_equi.vtf"

set output_analyze_path ./output_analyze/
set output_analyze_file "${simulation_identifier}_analyze.txt"

set output_npt_path ./output_npt/
set output_npt_file "${simulation_identifier}_npt_info.txt"

set output_checkpointing_path ./checkpointing/
set output_checkpointing_file "${simulation_identifier}_checkpointing.txt"
set output_sub_checkpointing_path "${output_checkpointing_path}subcheckpoints/"

set input_simulation_config_path ./simulation_config/

# also see output paths defined after main simulation loop



# Check in/out directories
#############################################################

proc makedirexist {path} {
    if { ![file exists $path] } {
	file mkdir $path
    } {
	if { ![file isdir $path]} {
	    puts "$path has to be a directory. Remove it or change the path setup"
	    exit
	}
    }
}

makedirexist $output_vtf_path
makedirexist $output_analyze_path
makedirexist $output_npt_path
makedirexist $output_checkpointing_path
makedirexist $output_sub_checkpointing_path
#makedirexist $input_simulation_config_path
  #the last line is probably not necessary...


#try to read in simulation configuration file
if { [file exists $input_simulation_config_path$input_simulation_config_file] } {
    puts "Read in Simulation Configuration file '$input_simulation_config_file'..."
    
    set simulation_config [open "$input_simulation_config_path$input_simulation_config_file" "r"]
    while { [blockfile $simulation_config read auto] != "eof" } {}
    close $simulation_config
} else {
    puts "Cannot find Simulation Configuration file '$input_simulation_config_file'. Quit Simulation with Error."
    exit 1
}


if { [info exists enableCoulombInt] } {
    set allowCoulomb $enableCoulombInt
} else {
    set allowCoulomb 1
}


if { ![info exists NPTEqui] } {
    set NPTEqui 0
}

if { ![info exists vmd_writeCoordinatesInterval] } {
    set vmd_writeCoordinatesInterval 1
}


#############################################################
#  Parameters                                               #
#############################################################

# System parameters
#############################################################
set sigma 1.0

#if this variable is true, a neutral backbone with charged side chains will be created. otherwise the backbone is charged and no side chains are attached
#(read from simulation config)
#set createSidechains [lindex $argv 2]

#if { ![string is digit $createSidechains] || $createSidechains == "" } {
#    puts "Set 'createSidechains' to 0"
#    set createSidechains 0
#}

#if this variable is true, the integrated counterion distribution P(r) will be calculated in each analyze step and averaged before writing after every checkpoint interval
#(read from simulation config)
#set calculateIntegratedCIDistribution 0

#set N_P 1
   # number of polymers (read from simulation config)
#set MPC 300
   # number of monomers per polymer chain (= its length) (read from simulation config)
set bond_l [expr 1.0*$sigma]
   # length of the bonds between monomers
set bond_l_side [expr 1.0*$sigma]
   # length of the bonds between monomers of the side beads

set val_cM 1.0
   # valency of charged monomers

set val_CI -1.0
   # valency of counterions

#set density 0.00005
   # particle density of the system (from which the box_l will be derived), this will be ignored if monomer_density is given
   # (read from simulation config)
#set monomer_density 0.0005
  # monomer density of the system (from which the box_l will be derived)

set start_particle_id 0
    # defines the id of the first monomer of the backbone that is connected with a side bead

#set particle_id_shift 1000000
    # defines the number of which the particle ids of the side monomers are shifted (first side monomer on the first sidechain of the first polymer has the id particle_id_shift, second side monomer on the first sidechain of the first polymer has the id 2*particle_id_shift and so on, this implicitly defines the maximum number of polymer beads+counterion beads

#set side_bead_distance 3
    # defines the distance between monomers of the backbone that are connected with a side bead
    # (read from simulation config)

set N_side_chains_per_polymer [expr int(ceil(double($MPC - $start_particle_id)/$side_bead_distance))]


# System setup options
#############################################################

set poly_mode SAW
set counterion_mode SAW
#set shield [expr 0.1*$sigma]
#set shield [expr 0.6*$sigma]
#set shield [expr 0.5*$sigma]
set shield [expr 0.27*$sigma]
   # $poly_mode=SAW uses a Self-Avoiding Random Walk to place the monomers, observing an excluded volume $shield around each particle 
   # $poly_mode=RW  uses a simple random walk without any additional fancy stuff
#set max_try 30000
#set max_try 1000000
set max_try 1000000
   # how many trials to attempt to place a particle, crosslink the network, etc.


# Interaction parameters
#############################################################

# Set the WCA interaction parameters
set LJ_sig $sigma
#set LJ_eps 1.0
#(read from simulation config)
#set LJ_cut_factor 1.122462048
#(read from simulation config)
set LJ_cut [expr ($LJ_cut_factor*$LJ_sig)]
set LJ_shift  [calc_lj_shift $LJ_sig $LJ_cut]
#set lj2_shift  [calc_lj_shift 1 $lj2_cut]


# Set the interaction parameters for the attractive LJ potential
set LJ2_sig $sigma
#set LJ2_eps 1.75
#(read from simulation config)
#set LJ2_cut_factor 2.5
#(read from simulation config)
set LJ2_cut [expr ($LJ2_cut_factor*$LJ2_sig)]
set LJ2_shift  [calc_lj_shift $LJ2_sig $LJ2_cut]



# Set the FENE interaction parameters (read from simulation config)
#set fene_k       7.0
#set fene_r       2.0

# Set the Coulomb interaction parameters
#set l_bjerrum      3.0
# (read from simulation config)

set p3m_accuracy 1e-4


# Integration parameters
#############################################################
# parameters for initial warm-up process with capped LJ-interactions
#set warmup_dt 0.00005
set warmup_dt 0.0001
   # time step to be used during initial warm-up
#set warmup_steps_per_iteration 20000
set warmup_steps_per_iteration 10000
#set warmup_steps_per_iteration 1000
set warmup_forcecap_start 10
set warmup_forcecap_incr 2
set warmup_forcecap_stop 50
set warmup_temp 1.0
set warmup_gamma 5.0
#set warmup_gamma 100000.0
   # temperature friction coefficient to be used during the warm-up (may be even up to 1e4 if desired)


#parameters for equilibration process (read from simulation config)
#set equi_dt 0.01
set equi_steps_per_subiteration 10
set equi_steps_iteration_factor 10
set equi_steps_per_iteration [expr $equi_steps_per_subiteration*$equi_steps_iteration_factor]

#set equi_iterations 1000 #not needed since simulation_max_time is used
#set equi_temp 1.0
#set equi_gamma 1.0

#set the checkpoint iteration interval (after how many iterations the checkpoint file is rewritten)
#set checkpoint_iteration_interval 500 (read from simulation config)

#set the interval after which the integrated counterion distribution is calculated and stored for later averaging
#set distribution_iteration_interval 10 (read from simulation config)

#set the interval after which the structure factor S(q) is calculated and stored for later averaging
#set structurefactor_iteration_interval 10 (read from simulation config)

#set the interval after which the spherically averaged form factor of a single chain S_1(q) is calculated and stored for later averaging
#set formfactor_iteration_interval 10 (read from simulation config)




# Espresso parameters
#############################################################

#better tune skin with tune_skin min max tol steps (see ug, sec 7.4), see equilibration section below
setmd skin 0.3

# bond types
set type_FENE 0


# particle types
set type_nM 1
   # neutral monomer
set type_cM 2
   # charged monomer
set type_nSM 3
   # neutral side monomer
set type_cSM 4
   # charged side monomer
set type_CI 5
   #counterion
set types [list $type_nM $type_cM $type_nSM $type_cSM $type_CI]
set types_excl_CI [list $type_nM $type_cM $type_nSM $type_cSM]

set types_rdf_backbone [list $type_nM]
set types_rdf_charged_polymer_beads [list $type_cM $type_cSM]
set types_rdf_CI [list $type_CI]


# Other parameters
#############################################################

#set tcl_precision 5
#set tcl_precision 7
   # tells tcl which precision to use for formating the output
set pi [expr acos(-1.0)]

#ranges for the angles to attach monomers in spherical coordinates
set phi_min 0.0
set phi_max [expr 2*$pi]
set theta_min 0.0
set theta_max $pi

### variables that are possibly overwritten by checkpoint ##############################################
# loop counter for npt runs
set j_start 0
# npt history variables for determining pressure convergence
#set npt_history_avg_rel_pressure_err 0.0
set npt_history_avg_pressure 0.0
if { [info exists NPTInitPistonMass] } {
    set piston_mass $NPTInitPistonMass
} else {
    set piston_mass 1e-5
}

#nvt subiteration counter
set subiteration_counter_start 0
########################################################################################################





# Complete parameter set
#############################################################

if { [info exists backboneFixed] && $backboneFixed } {
    set box_l_xy [expr ($bond_l*$monomer_density)**(-0.5)]
    #also define box_l variable for CIdistribution range
    set box_l $box_l_xy
    set box_l_z [expr $MPC*$bond_l]
    setmd box_l $box_l_xy $box_l_xy $box_l_z
    #puts "Box length was set to ($box_l_xy, $box_l_xy, $box_l_z) (fixed backbone)."
    puts "Box length was set to ([setmd box_l]) (fixed backbone)."
} else {
    set box_l  [expr pow(($N_P * $MPC / $monomer_density),(1.0/3.0))]
    setmd box_l $box_l $box_l $box_l
    #puts "Box length was set to ($box_l, $box_l, $box_l)."
    puts "Box length was set to ([setmd box_l])."
}


#set periodic boundary conditions (already default value)
setmd periodic 1 1 1



# Analyze Variables
#############################################################
set list_simulation_time ""
set list_end_to_end_distance ""
set list_end_to_end_distance_error ""
set list_radius_of_gyration ""
set list_radius_of_gyration_error ""
set list_hydrodynamic_radius ""
set list_hydrodynamic_radius_error ""
set list_energy_total ""
set list_energy_kinetic ""
set list_energy_coulomb ""
set list_boxlength ""

# integrated counterion distribution
set distr_r_min [expr (0.1*$sigma)]
#set distr_r_max $box_l ##################################### SETTING FROM SINGLE CHAIN SIMULATIONS #############################################################################
#adjusted to box_l/2
set distr_r_max [expr $box_l/2.0]
#set distr_bins 100 (set in configuration file)
if { ![info exists distr_log_flag] } {
    #distr_log_flag is now set in simulation config file
    #in case distr_log_flag is not set (old format) do it here
    set distr_log_flag 1
}
if { ![info exists distr_int_flag] } {
    #distr_int_flag is now set in simulation config file
    #in case distr_int_flag is not set (old format) do it here
    set distr_int_flag 1
}

if { $calculateIntegratedCIDistribution } {
    #create initial list for average integrated counterion distribution
    for {set i 0} {$i < $distr_bins} {incr i} {
        lappend avg_int_CI_distr 0.0
    }
}

set rdf_r_min [expr (0.1*$sigma)]
set rdf_r_max [expr $box_l/2.0]

if { [info exists calculateRDF] && $calculateRDF } {
    #create initial lists for average rdf values
    #6 lists in total:
    #1. neutral backbone - neutral backbone
    #2. charged polymer beads - charged polymer beads
    #3. charged polymer beads - neutral backbone
    #4. neutral backbone - counterions
    #5. charged polymer beads - counterions
    #6. counterions - counterions
    
    
    #set rdf_rlist ""
    #set rdf_list1_nb_nb ""
    #set rdf_list2_cpb_cpb ""
    #set rdf_list3_cpb_nb ""
    #set rdf_list4_nb_ci ""
    #set rdf_list5_cpb_ci ""
    #set rdf_list6_ci_ci ""
    
    #varibale rdf_bins exists in all configuration files containing the variable calculateRDF (no compatibility problems with old config files)
    
    for {set i 0} {$i < $rdf_bins} {incr i} {
        lappend avg_rdf_list1_nb_nb 0.0
        lappend avg_rdf_list2_cpb_cpb 0.0
        lappend avg_rdf_list3_cpb_nb 0.0
        lappend avg_rdf_list4_nb_ci 0.0
        lappend avg_rdf_list5_cpb_ci 0.0
        lappend avg_rdf_list6_ci_ci 0.0
    }
}
if { ![info exists calculateRDF] } {
    #to avoid compatibility problems for the code below the integration
    set calculateRDF 0
}




# structure factor variables are set below after system is set up

# spherically averaged form factor of a single chain
#set formfac_q_min_fac 0.01 (set in configuration file)
set formfac_q_min [expr ($formfac_q_min_fac*$sigma**(-1))]
#set formfac_q_max_fac 100. (set in configuration file)
set formfac_q_max [expr ($formfac_q_max_fac*$sigma**(-1))]
#set formfac_bins 100 (set in configuration file)
if { $calculateFormFactor } {
    #create initial list for average form factor
    for {set i 0} {$i <= $formfac_bins} {incr i} {
        lappend avg_formfac 0.0
        lappend avg_formfac_polymer 0.0        
    }
}

# Interaction setup
#############################################################
puts "    Creating interactions..."

# fene r0 is identical for polymer and side bead - need to be adjusted if bond_l != bond_l_side
#inter $type_FENE fene $fene_k $fene_r $bond_l
inter $type_FENE fene $fene_k $fene_r
puts "        [inter $type_FENE]"


# pairwise lennard_jones (WCA) for all particles
foreach i $types {
    foreach j $types {
        inter $i $j lennard-jones $LJ_eps $LJ_sig $LJ_cut $LJ_shift 0.0
        puts "        [inter $i $j]"
    }
}

# pairwise WCA potential between counterions and all other particles
#foreach i $types_excl_CI {
#    inter $type_CI $i lennard-jones $LJ_eps $LJ_sig $LJ_cut $LJ_shift 0.0
#    puts "        [inter $type_CI $i]"
#}

# overwrite WCA with pairwise attractive LJ potential for all particles except the counterions
foreach i $types_excl_CI {
    foreach j $types_excl_CI {
        inter $i $j lennard-jones $LJ2_eps $LJ2_sig $LJ2_cut $LJ2_shift 0.0
        puts "        [inter $i $j]"
    }
}

puts "    Interaction setup complete."


#############################################################
#  Checkpointing                                            #
#############################################################

if { [file exists "$output_checkpointing_path$output_checkpointing_file"] } {
    #checkpointing file found - load configurations
    puts "Found simulation checkpoint '$output_checkpointing_file'."

    set checkpoint [open "$output_checkpointing_path$output_checkpointing_file" "r"]
    while { [blockfile $checkpoint read auto] != "eof" } {}
    close $checkpoint

    # Analysis options & I/0-options
    #############################################################
    set vtf_file_equi [open "$output_vtf_path$output_vtf_file_equi" "a"]
    set analyze_file_equi [open "$output_analyze_path$output_analyze_file" "a"]

    #change the buffersize to 1mb and write only to vtf file when the buffer is full - avoiding checkpoint problems
    fconfigure $vtf_file_equi -buffersize 1048576
    fconfigure $vtf_file_equi -buffering full
    
} else {
    #start a new simulation
    puts "No checkpoint found - start a new simulation."

    # Analysis options & I/0-options
    #############################################################
    set vtf_file_warmup [open "$output_vtf_path$output_vtf_file_warmup" "w"]
    set vtf_file_equi [open "$output_vtf_path$output_vtf_file_equi" "w"]

    #change the buffersize to 1mb and write only to vtf file when the buffer is full - avoiding checkpoint problems
    fconfigure $vtf_file_equi -buffersize 1048576
    fconfigure $vtf_file_equi -buffering full


    # Particle setup
    #############################################################
    #check which particle setup should be created
    if { $createSidechains } {
        #sidechains.
        if { [info exists backboneFixed] && $backboneFixed } {
            puts "Create single neutral and fixed polymer chain with charged side chains of length $createSidechains attached."
            #source src_sidechains_fixed.tcl
            source src_var_sidechains_fixed.tcl
        } else {
            puts "Create neutral polymer chain with charged side chains of length $createSidechains attached."
            #source src_sidechains.tcl
            source src_var_sidechains.tcl
        }
    } else {
        #no sidechains
        if { [info exists backboneFixed] && $backboneFixed } {
            puts "Create single charged and fixed polymer chain with no side chains attached."
            #source src_nosidechains_fixed.tcl
            source src_var_sidechains_fixed.tcl
        } else {
            puts "Create charged polymer with no side chains attached."
            #source src_nosidechains.tcl
            source src_var_sidechains.tcl
        }
    }



    #############################################################
    #  Warm-up Integration                                      #
    #############################################################

    puts "\n Warm-up integration"

    thermostat langevin $warmup_temp $warmup_gamma

    setmd time_step $warmup_dt
    
    #puts "Particle properties:\n[part]"
    #set all_particle_ids ""
    #foreach value [part] {
    #    #puts [lindex $value 0]
    #    lappend all_particle_ids [lindex $value 0]
    #}
    #puts $all_particle_ids

    #write structe in vtf file - structures are writed regardless of vmd_writeCoordinates flag - this enables to export coordinates from checkpoint
    writevsf $vtf_file_warmup short
    #puts "writevsf executed"
    #writevsf $vtf_file_warmup verbose

    #write out initial coordinates in vtf file
    if { $vmd_writeCoordinates } {
        writevcf $vtf_file_warmup short
    }

    for {set cap $warmup_forcecap_start} {$cap < $warmup_forcecap_stop} {incr cap $warmup_forcecap_incr} {
        inter forcecap $cap
    
        integrate $warmup_steps_per_iteration
        puts "forcecap=$cap"
        if { $vmd_writeCoordinates } {
            writevcf $vtf_file_warmup short
        }
    
        #puts $out_ecoul "$cap [analyze energy coulomb]"
    }
    inter forcecap 0
    
    if { [info exists coulombWarmup] && $coulombWarmup && $allowCoulomb } {
        puts "Warm-up for force caped coulomb interactions..."
        #turn on coulomb interaction
        puts "turn on coulomb interactions..."
        #puts "[inter coulomb $l_bjerrum p3m tune accuracy $p3m_accuracy mesh 16]"
        puts "[inter coulomb $l_bjerrum p3m tunev2 accuracy $p3m_accuracy]"
        #puts "[inter coulomb $l_bjerrum p3m tune accuracy $p3m_accuracy]"
        puts "Interactions are now: {[inter]}"
        
        for {set cap $warmup_forcecap_start} {$cap < $warmup_forcecap_stop} {incr cap $warmup_forcecap_incr} {
            inter forcecap $cap
        
            integrate $warmup_steps_per_iteration
            if { $vmd_writeCoordinates } {
                writevcf $vtf_file_warmup short
            }
        
        
        }
        inter forcecap 0        
    }
    
    puts "Warm-up complete."

    flush $vtf_file_warmup
    close $vtf_file_warmup

    #reset the simulation time to 0 (warmup time should not be counted)
    setmd time 0.0


    #write first vtf lines for equilibration iteration
    #write structe in vtf file - structures are writed regardless of vmd_writeCoordinates flag - this enables to export coordinates from checkpoint
    writevsf $vtf_file_equi short

    #write first line in analyze file
    set analyze_file_equi [open "$output_analyze_path$output_analyze_file" "w"]
    puts $analyze_file_equi "\#time\tr_e\terr_r_e\tr_g\terr_r_g\tr_h\terr_r_h\tE_tot\tE_kin\tE_coul\tbox_l_x"
    flush $analyze_file_equi
}

### end of checkpointing ###############################################################################


# structure factor (better be done in postprocessing with python?)
if { $createSidechains } {
    set strufac_type $type_nM
} else {
    set strufac_type $type_cM
      #this might cause unwanted artefacts in the produced structure factor data (other spacing when $side_bead_distance>1)
}
#set strufac_order 3 (set in configuration file)

if { $calculateStructureFactor } {
    #create initial list for average structure factor
    set length_strufac_list [llength [analyze structurefactor $strufac_type $strufac_order]]
    for {set i 0} {$i < $length_strufac_list} {incr i} {
        lappend avg_strufac 0.0
    }
}





#############################################################
#  Equilibration                                            #
#############################################################
#do not turn on coulomb interactions again in case of a new simulation (no checkpoint) and coulombWarmup activated, allowCoulomb must also be true to turn on coulomb interactions
if { ( [file exists "$output_checkpointing_path$output_checkpointing_file"] || !([info exists coulombWarmup] && $coulombWarmup ) ) && $allowCoulomb } {
    #turn on coulomb interaction after warmup
    puts "turn on coulomb interactions..."
    #puts "[inter coulomb $l_bjerrum p3m tune accuracy $p3m_accuracy mesh 16]"
    puts "[inter coulomb $l_bjerrum p3m tunev2 accuracy $p3m_accuracy]"
    #puts "[inter coulomb $l_bjerrum p3m tune accuracy $p3m_accuracy]"
    puts "Interactions are now: {[inter]}"
    
    #find an optimized skin value first
    #tune_skin (vgl. user-guide 7.4)
}

#write coordinates in vtf file
#writevcf $vtf_file_equi short

proc writeCheckpoint {} {
    global list_simulation_time list_end_to_end_distance list_end_to_end_distance_error list_radius_of_gyration list_radius_of_gyration_error list_hydrodynamic_radius list_hydrodynamic_radius_error list_energy_total list_energy_kinetic list_energy_coulomb list_boxlength analyze_file_equi vtf_file_equi output_checkpointing_path output_checkpointing_file save_all_checkpoints output_sub_checkpointing_path simulation_identifier j_start npt_history_avg_pressure piston_mass subiteration_counter_start
    
    puts "Writing checkpoint file..."
        
    #update analyze file
    foreach t $list_simulation_time re $list_end_to_end_distance err_re $list_end_to_end_distance_error rg $list_radius_of_gyration err_rg $list_radius_of_gyration_error rh $list_hydrodynamic_radius err_rh $list_hydrodynamic_radius_error etot $list_energy_total ekin $list_energy_kinetic ecoul $list_energy_coulomb boxl $list_boxlength {
        puts $analyze_file_equi "$t\t$re\t$err_re\t$rg\t$err_rg\t$rh\t$err_rh\t$etot\t$ekin\t$ecoul\t$boxl"
    }
    flush $analyze_file_equi

    #write vcf data from buffer
    flush $vtf_file_equi

    #reset analyze list variables
    set list_simulation_time ""
    set list_end_to_end_distance ""
    set list_end_to_end_distance_error ""
    set list_radius_of_gyration ""
    set list_radius_of_gyration_error ""
    set list_hydrodynamic_radius ""
    set list_hydrodynamic_radius_error ""
    set list_energy_total ""
    set list_energy_kinetic ""
    set list_energy_coulomb ""
    set list_boxlength ""

    set checkpoint [open "$output_checkpointing_path$output_checkpointing_file" "w"]
    blockfile $checkpoint write variable all
    #blockfile $checkpoint write tclvariable all
    blockfile $checkpoint write tclvariable { j_start npt_history_avg_pressure piston_mass subiteration_counter_start }
    blockfile $checkpoint write particles "id type pos q v f fix" all
    blockfile $checkpoint write bonds all
    #blockfile $checkpoint write interactions
    close $checkpoint
    
    if { [info exists save_all_checkpoints] && $save_all_checkpoints } {
        #copy current checkpoint with simulation time label to checkpointing subdir
        file copy -force "$output_checkpointing_path$output_checkpointing_file" "$output_sub_checkpointing_path${simulation_identifier}_[format "%.4f" [setmd time]]_checkpointing.txt"
    }
}


# instructions to write last checkpoint in case the simulation is aborted
############################################################################
proc doOnExit {} {
    #write checkpoint
    puts "\nSimulation aborted - write last checkpoint and quit..."
    
    global list_simulation_time list_end_to_end_distance list_end_to_end_distance_error list_radius_of_gyration list_radius_of_gyration_error list_hydrodynamic_radius list_hydrodynamic_radius_error list_energy_total list_energy_kinetic list_energy_coulomb list_boxlength analyze_file_equi vtf_file_equi output_checkpointing_path output_checkpointing_file save_all_checkpoints output_sub_checkpointing_path simulation_identifier j_start npt_history_avg_pressure piston_mass subiteration_counter_start

    #update analyze file
    foreach t $list_simulation_time re $list_end_to_end_distance err_re $list_end_to_end_distance_error rg $list_radius_of_gyration err_rg $list_radius_of_gyration_error rh $list_hydrodynamic_radius err_rh $list_hydrodynamic_radius_error etot $list_energy_total ekin $list_energy_kinetic ecoul $list_energy_coulomb boxl $list_boxlength {
        puts $analyze_file_equi "$t\t$re\t$err_re\t$rg\t$err_rg\t$rh\t$err_rh\t$etot\t$ekin\t$ecoul\t$boxl"
    }
    flush $analyze_file_equi
    close $analyze_file_equi

    #write vcf data from buffer
    flush $vtf_file_equi
    close $vtf_file_equi

    #reset analyze list variables
    set list_simulation_time ""
    set list_end_to_end_distance ""
    set list_end_to_end_distance_error ""
    set list_radius_of_gyration ""
    set list_radius_of_gyration_error ""
    set list_hydrodynamic_radius ""
    set list_hydrodynamic_radius_error ""
    set list_energy_total ""
    set list_energy_kinetic ""
    set list_energy_coulomb ""
    set list_boxlength ""

    
    set checkpoint [open "$output_checkpointing_path$output_checkpointing_file" "w"]
    blockfile $checkpoint write variable all
    #blockfile $checkpoint write tclvariable all
    blockfile $checkpoint write tclvariable { j_start npt_history_avg_pressure piston_mass subiteration_counter_start }
    blockfile $checkpoint write particles "id type pos q v f fix" all
    blockfile $checkpoint write bonds all
    #blockfile $checkpoint write interactions
    #do not write interactions in checkpointing blockfile - new p3m optimization in the next run necessary (rarely occuring problem: Complex value is not zero)
    close $checkpoint
    
    if { [info exists save_all_checkpoints] && $save_all_checkpoints } {
        #copy current checkpoint with simulation time label to checkpointing subdir
        file copy -force "$output_checkpointing_path$output_checkpointing_file" "$output_sub_checkpointing_path${simulation_identifier}_[format "%.4f" [setmd time]]_checkpointing.txt"
    }

    exit
}


signal trap SIGINT doOnExit
signal trap SIGTERM doOnExit
signal trap SIGUSR1 doOnExit
signal trap 18 doOnExit
#signal 18 is sent when time limit is reached




proc set_piston {new_Vkappa} {
    global NPTGammaV
    global NPTGamma0
    global equi_temp
    global NPTPressure
    global NPTMinPistonMass
    global NPTMaxPistonMass
    global piston_mass
    set piston_mass [expr 4.0/($NPTGamma0*$NPTGamma0*$new_Vkappa)]
    set NPTGammaV [expr $piston_mass*$NPTGamma0]
    #puts "pistonmass $piston_mass\ngamma0 $NPTGamma0\ngammaV $NPTGammaV"
    #lower bound, should be irrelevant for a solution that should self-consistently tend to the right values
    if { [expr $piston_mass < $NPTMinPistonMass] } {
        puts "Piston Mass Lower Bound"
        set piston_mass $NPTMinPistonMass
        set NPTGammaV [expr $piston_mass*$NPTGamma0]
    }
    if { [expr $piston_mass > $NPTMaxPistonMass] } {
        puts "Piston Mass Upper Bound"
        set piston_mass $NPTMaxPistonMass
        set NPTGammaV [expr $piston_mass*$NPTGamma0]
    }
    
    puts "set_piston called, new values calculated: piston_mass=$piston_mass, gamma_0=$NPTGamma0, gamma_v=$NPTGammaV"

    thermostat npt_isotropic $equi_temp $NPTGamma0 $NPTGammaV
    integrate set npt_isotropic $NPTPressure $piston_mass
}


set iteration_counter 0
set distribution_counter 0
set rdf_counter 0
set structurefactor_counter 0
set formfactor_counter 0

if { $add_sim_time } {
    #the given time console parameter is in this case not the absolute maximum simulation time but the additional time
    #therefore increase variable for compatibility with below code
    set simulation_max_time [expr $simulation_max_time+[setmd time]]
}


setmd time_step $equi_dt

if { [info exists vmd_writeHeader] && $vmd_writeHeader } {
    writevsf $vtf_file_equi short
    flush $vtf_file_equi
}

if { !$NPTEqui } {
    puts "\n Equilibration integration (NVT)"
    thermostat langevin $equi_temp $equi_gamma
    puts [integrate set]
    integrate set nvt
    puts [integrate set]
} else {
    puts "\n Equilibration integration (NPT)"
    
    #piston_mass comes from checkpoint (otherwise use NPTInitPistonMass as default)
    set NPTGammaV [expr $piston_mass*$NPTGamma0]
    
    thermostat off
    puts [thermostat]
    thermostat npt_isotropic $equi_temp $NPTGamma0 $NPTGammaV
    puts [thermostat]
    integrate set npt_isotropic $NPTPressure $piston_mass
    
    #npt equilibration runs
    set npt_successful_iterations 0
    #set NPTRequiredSuccessfulIterations 500
    
    puts "j_start = $j_start"
    
    #write out boxlength during NPT equilibration
    if { $j_start == 0 } {
        #write header line only for the first time
        set npt_boxl_file [open "${output_npt_path}${output_npt_file}" "w"]
        puts $npt_boxl_file "\#time\tbox_lx\tbox_ly\tbox_lz\tpiston_mass\tgamma_v\tp_inst\tVkappa\tavgRelPressErr"
        flush $npt_boxl_file
    } else {
        set npt_boxl_file [open "${output_npt_path}${output_npt_file}" "a"]
    }
    
    
    ### main integration loop for NPT runs #############################################################
    
    for { set j $j_start } { $npt_successful_iterations < $NPTRequiredSuccessfulIterations && [setmd time] < $simulation_max_time } { incr j } {
        puts "NPT Equilibration iteration $j, Simulation Time: [format "%.4f" [setmd time]], Successful Iterations: $npt_successful_iterations"
        
        for { set subiteration_counter 0 } { $subiteration_counter < $equi_steps_iteration_factor } { incr subiteration_counter } {
            integrate $equi_steps_per_subiteration
        }
        
        if { $vmd_writeCoordinates && [expr ($j % $vmd_writeCoordinatesInterval)] == 0 } {
            writevcf $vtf_file_equi
        }
        
        set Vkappa [expr abs([analyze Vkappa])]
        
        if { $j <= $NPTRequiredSuccessfulIterations } {
            set npt_avg [expr 1.0*$j]
        } else { 
            set npt_avg [expr 1.0*$NPTRequiredSuccessfulIterations]
        }
        
        set tmp_abs_pressure_err [expr abs([setmd npt_p_inst] - $NPTPressure)]
        set tmp_rel_pressure_err [expr $tmp_abs_pressure_err/$NPTPressure]
        
        #adjusted running averaging
        #set npt_history_avg_rel_pressure_err [expr $npt_avg/($npt_avg+1.0)*$npt_history_avg_rel_pressure_err+1.0/($npt_avg+1.0)*$tmp_rel_pressure_err]
        set npt_history_avg_pressure [expr $npt_avg/($npt_avg+1.0)*$npt_history_avg_pressure+1.0/($npt_avg+1.0)*[setmd npt_p_inst]]
        
        set tmp_pressure_abs_error [expr abs($npt_history_avg_pressure - $NPTPressure)]
        set tmp_pressure_rel_error [expr $tmp_pressure_abs_error/$NPTPressure]
        
        #previous if condition: $npt_history_avg_rel_pressure_err <= $NPTPressureRelativeErrorTolerance
        if { $tmp_pressure_rel_error <= $NPTPressureRelativeErrorTolerance } {
            incr npt_successful_iterations
        } else {
            set npt_successful_iterations 0
        }
        
        puts $npt_boxl_file "[setmd time]\t[lindex [setmd box_l] 0]\t[lindex [setmd box_l] 1]\t[lindex [setmd box_l] 2]\t[setmd npt_piston]\t[setmd nptiso_gammav]\t[setmd npt_p_inst]\t$Vkappa\t$npt_history_avg_pressure"
        flush $npt_boxl_file
        
        if { [expr ($j % $NPTRecalculatePistonMassIterationInterval)] == 0 && $Vkappa > 0} {
            #calculate new piston mass from volume fluctuations from time to time...
            set_piston $Vkappa
            analyze Vkappa reset
        }
        
        if { [expr ($j % $checkpoint_iteration_interval)] == 0 } {
            #write checkpoint
            writeCheckpoint
        }
        
        incr j_start
    }
    
}


### main integration loop for NVT runs #############################################################
#probably adjust boxlength (e. g. from previous npt simulations)
set box_l_init [lindex [setmd box_l] 0]
#check if given monomer concentration corresponds to actual boxlength...
#puts "Boxlength for given monomer concentration: $box_l"
#puts "Actual boxlength: $box_l_init"
if { $box_l_init != $box_l && !$NPTEqui} {
    puts "Actual boxlength does not match the corresponding boxlength from given monomer density. Adjust box..."
    
    #box_l_change_factor * box_l_init = max change of boxlength per iteration
    set box_l_change_factor 0.0001
    
    if { [expr $box_l_init - $box_l] > 0.0 } {
        set box_l_change_factor [expr -1.0*$box_l_change_factor]
    }
    
    set change_boxl 1
} else {
    set change_boxl 0
}

while { [setmd time] < $simulation_max_time && !$NPTEqui} {
    #equilibration iteration for NVT ensemble
    puts "Equilibration iteration $iteration_counter, Simulation Time: [format "%.4f" [setmd time]]"
    
    if { $change_boxl } {
        #linearly change boxsize until sign of the differences would change in the next step
        if { [expr ([lindex [setmd box_l] 0] - $box_l) * ([lindex [setmd box_l] 0]+$box_l_change_factor*$box_l_init - $box_l)] > 0 } {
            change_volume [expr [lindex [setmd box_l] 0]+$box_l_change_factor*$box_l_init] xyz
            puts "New Boxlength: [lindex [setmd box_l] 0]"
        } else {
            change_volume $box_l xyz
            puts "Final Boxlength: [lindex [setmd box_l] 0]"
        }
    }

    
    for { set subiteration_counter $subiteration_counter_start } { $subiteration_counter < $equi_steps_iteration_factor } { incr subiteration_counter } {
        integrate $equi_steps_per_subiteration
    }
    set subiteration_counter_start 0
    #integrate $equi_steps_per_iteration

    lappend list_simulation_time [setmd time]
    lappend list_end_to_end_distance [lindex [analyze re 0 $N_P $MPC] 0]
    lappend list_end_to_end_distance_error [lindex [analyze re 0 $N_P $MPC] 1]
    #lappend list_radius_of_gyration [lindex [analyze rg 0 $N_P $MPC] 0]
    #lappend list_radius_of_gyration_error [lindex [analyze rg 0 $N_P $MPC] 1]
    #lappend list_hydrodynamic_radius [lindex [analyze rh 0 $N_P $MPC] 0]
    #lappend list_hydrodynamic_radius_error [lindex [analyze rh 0 $N_P $MPC] 1]
    #calculation of R_G and R_H for entire chain (side chains included)
    lappend list_radius_of_gyration [lindex [analyze rg 0 $N_P [expr $MPC+$N_side_chains_per_polymer*$createSidechains]] 0]
    lappend list_radius_of_gyration_error [lindex [analyze rg 0 $N_P [expr $MPC+$N_side_chains_per_polymer*$createSidechains]] 1]
    lappend list_hydrodynamic_radius [lindex [analyze rh 0 $N_P [expr $MPC+$N_side_chains_per_polymer*$createSidechains]] 0]
    lappend list_hydrodynamic_radius_error [lindex [analyze rh 0 $N_P [expr $MPC+$N_side_chains_per_polymer*$createSidechains]] 1]
    lappend list_energy_total [analyze energy total]
    lappend list_energy_kinetic [analyze energy kinetic]
    lappend list_energy_coulomb [analyze energy coulomb]
    lappend list_boxlength [lindex [setmd box_l] 0]

    #this might cause problems in case the simulation is interrupted and later resumed (interpreter position on interrupt event unclear)
    if { $vmd_writeCoordinates && [expr ($iteration_counter % $vmd_writeCoordinatesInterval)] == 0 } {
        writevcf $vtf_file_equi
    }
    
    
    incr iteration_counter
    if { [expr ($iteration_counter % $distribution_iteration_interval)] == 0 && $calculateIntegratedCIDistribution } {
        #write integrated counterion distribution P(r) if enabled
        puts "Analyze counterion distribution..."
        if { $createSidechains } {
            #sidechains
            set distr [analyze distribution [list $type_CI] [list $type_cM $type_nM] $distr_r_min $distr_r_max $distr_bins $distr_log_flag $distr_int_flag]
        } else {
            #no sidechains
            set distr [analyze distribution [list $type_CI] [list $type_cM $type_nM $type_cSM $type_nSM] $distr_r_min $distr_r_max $distr_bins $distr_log_flag $distr_int_flag]
        }
        
        set distr_r_list ""
        set distr_pr_list ""
        foreach value [lindex $distr 1] {
            lappend distr_r_list [lindex $value 0]
            lappend distr_pr_list [lindex $value 1]
        }
        #set avg_int_CI_distr [vecadd $avg_int_CI_distr $distr_pr_list]
        set avg_int_CI_distr [vecadd [vecscale [expr $distribution_counter/($distribution_counter+1.0)] $avg_int_CI_distr] [vecscale [expr 1.0/($distribution_counter+1.0)] $distr_pr_list]]
        incr distribution_counter
        
        #puts "Writing distribution debug to blockfile"
        #set debug_distr_file [open "debug_distr.txt" "w"]
        #blockfile $debug_distr_file write tclvariable distr
        #close $debug_distr_file
        #puts "Done"
    }
    
    if { $calculateRDF && [expr ($iteration_counter % $rdf_iteration_interval)] == 0 } {
        #order of the if statements important for compatibility with old configuration files
        #write different rdfs if enabled
        puts "Calculate radial distribution functions..."
        
        set rdf_rlist ""
        set rdf_list1_nb_nb ""
        set rdf_list2_cpb_cpb ""
        set rdf_list3_cpb_nb ""
        set rdf_list4_nb_ci ""
        set rdf_list5_cpb_ci ""
        set rdf_list6_ci_ci ""
        
        set rdf1_nb_nb [analyze rdf $types_rdf_backbone $types_rdf_backbone $rdf_r_min $rdf_r_max $rdf_bins]
        set rdf2_cpb_cpb [analyze rdf $types_rdf_charged_polymer_beads $types_rdf_charged_polymer_beads $rdf_r_min $rdf_r_max $rdf_bins]
        set rdf3_cpb_nb [analyze rdf $types_rdf_charged_polymer_beads $types_rdf_backbone $rdf_r_min $rdf_r_max $rdf_bins]
        set rdf4_nb_ci [analyze rdf $types_rdf_backbone $types_rdf_CI $rdf_r_min $rdf_r_max $rdf_bins]
        set rdf5_cpb_ci [analyze rdf $types_rdf_charged_polymer_beads $types_rdf_CI $rdf_r_min $rdf_r_max $rdf_bins]
        set rdf6_ci_ci [analyze rdf $types_rdf_CI $types_rdf_CI $rdf_r_min $rdf_r_max $rdf_bins]
        
        foreach value [lindex $rdf1_nb_nb 1] {
            lappend rdf_rlist [lindex $value 0]
            if { [regexp "nan" [lindex $value 1]] } {
                lappend rdf_list1_nb_nb 0.0
                #puts "nan detected"
            } else {
                lappend rdf_list1_nb_nb [lindex $value 1]
            }
        }
        foreach value [lindex $rdf2_cpb_cpb 1] {
            if { [regexp "nan" [lindex $value 1]] } {
                lappend rdf_list2_cpb_cpb 0.0
            } else {
                lappend rdf_list2_cpb_cpb [lindex $value 1]
            }
            #lappend rdf_list2_cpb_cpb [lindex $value 1]
        }
        foreach value [lindex $rdf3_cpb_nb 1] {
            if { [regexp "nan" [lindex $value 1]] } {
                lappend rdf_list3_cpb_nb 0.0
            } else {
                lappend rdf_list3_cpb_nb [lindex $value 1]
            }
            #lappend rdf_list3_cpb_nb [lindex $value 1]
        }
        foreach value [lindex $rdf4_nb_ci 1] {
            if { [regexp "nan" [lindex $value 1]] } {
                lappend rdf_list4_nb_ci 0.0
            } else {
                lappend rdf_list4_nb_ci [lindex $value 1]
            }
            #lappend rdf_list4_nb_ci [lindex $value 1]
        }
        foreach value [lindex $rdf5_cpb_ci 1] {
            if { [regexp "nan" [lindex $value 1]] } {
                lappend rdf_list5_cpb_ci 0.0
            } else {
                lappend rdf_list5_cpb_ci [lindex $value 1]
            }
            #lappend rdf_list5_cpb_ci [lindex $value 1]
        }
        foreach value [lindex $rdf6_ci_ci 1] {
            if { [regexp "nan" [lindex $value 1]] } {
                lappend rdf_list6_ci_ci 0.0
            } else {
                lappend rdf_list6_ci_ci [lindex $value 1]
            }
            #lappend rdf_list6_ci_ci [lindex $value 1]
        }
        
        #puts $rdf_list1_nb_nb
        #puts $rdf_list2_cpb_cpb
        
        #set avg_rdf_list1_nb_nb [vecadd $avg_rdf_list1_nb_nb $rdf_list1_nb_nb]
        set avg_rdf_list1_nb_nb [vecadd [vecscale [expr $rdf_counter/($rdf_counter+1.0)] $avg_rdf_list1_nb_nb] [vecscale [expr 1.0/($rdf_counter+1.0)] $rdf_list1_nb_nb]]
        #set avg_rdf_list2_cpb_cpb [vecadd $avg_rdf_list2_cpb_cpb $rdf_list2_cpb_cpb]
        set avg_rdf_list2_cpb_cpb [vecadd [vecscale [expr $rdf_counter/($rdf_counter+1.0)] $avg_rdf_list2_cpb_cpb] [vecscale [expr 1.0/($rdf_counter+1.0)] $rdf_list2_cpb_cpb]]
        #set avg_rdf_list3_cpb_nb [vecadd $avg_rdf_list3_cpb_nb $rdf_list3_cpb_nb]
        set avg_rdf_list3_cpb_nb [vecadd [vecscale [expr $rdf_counter/($rdf_counter+1.0)] $avg_rdf_list3_cpb_nb] [vecscale [expr 1.0/($rdf_counter+1.0)] $rdf_list3_cpb_nb]]
        #set avg_rdf_list4_nb_ci [vecadd $avg_rdf_list4_nb_ci $rdf_list4_nb_ci]
        set avg_rdf_list4_nb_ci [vecadd [vecscale [expr $rdf_counter/($rdf_counter+1.0)] $avg_rdf_list4_nb_ci] [vecscale [expr 1.0/($rdf_counter+1.0)] $rdf_list4_nb_ci]]
        #set avg_rdf_list5_cpb_ci [vecadd $avg_rdf_list5_cpb_ci $rdf_list5_cpb_ci]
        set avg_rdf_list5_cpb_ci [vecadd [vecscale [expr $rdf_counter/($rdf_counter+1.0)] $avg_rdf_list5_cpb_ci] [vecscale [expr 1.0/($rdf_counter+1.0)] $rdf_list5_cpb_ci]]
        #set avg_rdf_list6_ci_ci [vecadd $avg_rdf_list6_ci_ci $rdf_list6_ci_ci]
        set avg_rdf_list6_ci_ci [vecadd [vecscale [expr $rdf_counter/($rdf_counter+1.0)] $avg_rdf_list6_ci_ci] [vecscale [expr 1.0/($rdf_counter+1.0)] $rdf_list6_ci_ci]]
        
        incr rdf_counter
    }
    
    if { [expr ($iteration_counter % $structurefactor_iteration_interval)] == 0 && $calculateStructureFactor } {
        puts "Calculate structure factor..."
        set strufac [analyze structurefactor $strufac_type $strufac_order]
        
        set strufac_q_list ""
        set strufac_sq_list ""
        foreach value $strufac {
            lappend strufac_q_list [lindex $value 0]
            lappend strufac_sq_list [lindex $value 1]
        }
        #set avg_strufac [vecadd $avg_strufac $strufac_sq_list]
        set avg_strufac [vecadd [vecscale [expr $structurefactor_counter/($structurefactor_counter+1.0)] $avg_strufac] [vecscale [expr 1.0/($structurefactor_counter+1.0)] $strufac_sq_list]]
        incr structurefactor_counter
    }
    
    if { [expr ($iteration_counter % $formfactor_iteration_interval)] == 0 && $calculateFormFactor } {
        puts "Calculate form factor (one chain)..."
        ########################## This part has to be adjusted in case N_P > 1 (manual averaging over all chains), easy solution: no ascending IDs, e.g. all first side monomers start for example with at least 1000000, all second side monomers with at least 2000000 and so on ###############################################
        #the following line has to be changed for N_P > 1
        #set formfac [analyze formfactor $formfac_q_min $formfac_q_max $formfac_bins 0 $N_P $MPC]
        set formfac [analyze formfactor $formfac_q_min $formfac_q_max $formfac_bins 0 1 $MPC]
        
        #formfactor for one single polymer chain
        set formfac_polymer [analyze formfactor $formfac_q_min $formfac_q_max $formfac_bins 0 1 [expr $MPC+$N_side_chains_per_polymer*$createSidechains]]
        
        set formfac_q_list ""
        set formfac_sq_list ""
        set formfac_polymer_sq_list ""
        
        foreach value $formfac {
            lappend formfac_q_list [lindex $value 0]
            lappend formfac_sq_list [lindex $value 1]
        }
        foreach value $formfac_polymer {
            lappend formfac_polymer_sq_list [lindex $value 1]
        }
        
        #set avg_formfac [vecadd $avg_formfac $formfac_sq_list]
        set avg_formfac [vecadd [vecscale [expr $formfactor_counter/($formfactor_counter+1.0)] $avg_formfac] [vecscale [expr 1.0/($formfactor_counter+1.0)] $formfac_sq_list]]
        #set avg_formfac_polymer [vecadd $avg_formfac_polymer $formfac_polymer_sq_list]
        set avg_formfac_polymer [vecadd [vecscale [expr $formfactor_counter/($formfactor_counter+1.0)] $avg_formfac_polymer] [vecscale [expr 1.0/($formfactor_counter+1.0)] $formfac_polymer_sq_list]]
        
        incr formfactor_counter
        
        #debugging
        #puts "Writing formfactor debug to blockfile"
        #set debug_formfac_file [open "debug_formfac.txt" "w"]
        #blockfile $debug_formfac_file write tclvariable formfac
        #close $debug_formfac_file
        #puts "Done"
    }
            
    if { [expr ($iteration_counter % $checkpoint_iteration_interval)] == 0 } {
        #write checkpoint
        writeCheckpoint
    }
}
puts "Equilibration complete."


flush $vtf_file_equi
close $vtf_file_equi


# Define paths containing simulation_max_time
#############################################################
# Output paths containing (probably adjusted) simulation_max_time

set output_intCIdistr_path ./output_analyze/

set output_intCIdistr_file "${simulation_identifier}_[format "%.4f" $simulation_max_time]_distribution.txt"

set output_rdf_path ./output_analyze/
set output_rdf_file "${simulation_identifier}_[format "%.4f" $simulation_max_time]_rdf.txt"

set output_strufac_path ./output_analyze/
set output_strufac_file "${simulation_identifier}_[format "%.4f" $simulation_max_time]_structurefactor.txt"

set output_formfac_path ./output_analyze/
set output_formfac_file "${simulation_identifier}_[format "%.4f" $simulation_max_time]_formfactor.txt"

set output_info_path ./simulation_info/
set output_info_file "${simulation_identifier}_[format "%.4f" $simulation_max_time]_all_tcl_variables.txt"

makedirexist $output_intCIdistr_path
makedirexist $output_rdf_path
makedirexist $output_strufac_path
makedirexist $output_formfac_path
makedirexist $output_info_path


#puts "Writing analysis to file '$output_analyze_file'..."
#set analyze_file_equi [open "$output_analyze_path$output_analyze_file" "w"]
#puts $analyze_file_equi "\#time\tr_e\tr_g\tE_tot\tE_kin\tE_coul"

#close $analyze_file_equi
#puts "    Done."

puts "Write the last configuration as checkpoint..."
set checkpoint [open "$output_checkpointing_path$output_checkpointing_file" "w"]
blockfile $checkpoint write variable all
#blockfile $checkpoint write tclvariable all
blockfile $checkpoint write tclvariable { j_start npt_history_avg_pressure piston_mass subiteration_counter_start }
blockfile $checkpoint write particles "id type pos q v f fix" all
blockfile $checkpoint write bonds all
#blockfile $checkpoint write interactions
close $checkpoint
puts "    Done."

puts "Writing info block file '$output_info_file' for tcl variables..."
set info_blockfile [open "$output_info_path$output_info_file" "w"]
blockfile $info_blockfile write tclvariable all
close $info_blockfile
puts "    Done."

#update analyze file
foreach t $list_simulation_time re $list_end_to_end_distance err_re $list_end_to_end_distance_error rg $list_radius_of_gyration err_rg $list_radius_of_gyration_error rh $list_hydrodynamic_radius err_rh $list_hydrodynamic_radius_error etot $list_energy_total ekin $list_energy_kinetic ecoul $list_energy_coulomb boxl $list_boxlength {
    puts $analyze_file_equi "$t\t$re\t$err_re\t$rg\t$err_rg\t$rh\t$err_rh\t$etot\t$ekin\t$ecoul\t$boxl"
}

flush $analyze_file_equi
close $analyze_file_equi

#create file containing the (integrated) counterion distribution P(r)
if { $calculateIntegratedCIDistribution && $distribution_counter} {
    puts "Write averaged counterion distributions to file '$output_intCIdistr_file'..."
    #set avg_int_CI_distr [vecscale [expr 1.0/$distribution_counter] $avg_int_CI_distr]
    set distribution_file [open "$output_intCIdistr_path$output_intCIdistr_file" "w"]
    #puts $analyze_file_equi "\#time\tr_e\tr_g\tE_tot\tE_kin\tE_coul"
    puts $distribution_file "\#r\tP(r)"
    #foreach r $distr_r_list pr $distr_pr_list  ######################################### IMPORTANT BUG FROM PREVIOUS SIMULATIONS FIXED ######################################################
    foreach r $distr_r_list pr $avg_int_CI_distr {
        puts $distribution_file "$r\t$pr"
    }
    close $distribution_file
    puts "    Done."
}

#create file containing the radial distribution functions rdf1-6
if { $calculateRDF && $rdf_counter} {
    puts "Write averaged radial distribution functions to file '$output_rdf_file'..."
    #set avg_rdf_list1_nb_nb [vecscale [expr 1.0/$rdf_counter] $avg_rdf_list1_nb_nb]
    #set avg_rdf_list2_cpb_cpb [vecscale [expr 1.0/$rdf_counter] $avg_rdf_list2_cpb_cpb]
    #set avg_rdf_list3_cpb_nb [vecscale [expr 1.0/$rdf_counter] $avg_rdf_list3_cpb_nb]
    #set avg_rdf_list4_nb_ci [vecscale [expr 1.0/$rdf_counter] $avg_rdf_list4_nb_ci]
    #set avg_rdf_list5_cpb_ci [vecscale [expr 1.0/$rdf_counter] $avg_rdf_list5_cpb_ci]
    #set avg_rdf_list6_ci_ci [vecscale [expr 1.0/$rdf_counter] $avg_rdf_list6_ci_ci]
    
    set rdf_file [open "$output_rdf_path$output_rdf_file" "w"]
    puts $rdf_file "\#r\trdf1\trdf2\trdf3\trdf4\trdf5\trdf6"
    foreach r $rdf_rlist rdf1 $avg_rdf_list1_nb_nb rdf2 $avg_rdf_list2_cpb_cpb rdf3 $avg_rdf_list3_cpb_nb rdf4 $avg_rdf_list4_nb_ci rdf5 $avg_rdf_list5_cpb_ci rdf6 $avg_rdf_list6_ci_ci {
        puts $rdf_file "$r\t$rdf1\t$rdf2\t$rdf3\t$rdf4\t$rdf5\t$rdf6"
    }
    close $rdf_file
    puts "    Done."
}

#create file containing the structure factor S(q)
if { $calculateStructureFactor && $structurefactor_counter} {
    puts "Write averaged structure factor to file '$output_strufac_file'..."
    #set avg_strufac [vecscale [expr 1.0/$structurefactor_counter] $avg_strufac]
    set structurefactor_file [open "$output_strufac_path$output_strufac_file" "w"]
    puts $structurefactor_file "\#q\tS(q)"
    #foreach q $strufac_q_list sq $strufac_sq_list  ######################################### IMPORTANT BUG FROM PREVIOUS SIMULATIONS FIXED ######################################################
    foreach q $strufac_q_list sq $avg_strufac {
        puts $structurefactor_file "$q\t$sq"
    }
    close $structurefactor_file
    puts "    Done."
}

#create file containing the spherically averaged form factor of a single chain S_1(q)
if { $calculateFormFactor && $formfactor_counter } {
    puts "Write averaged spherically averaged form factor of a single chain to file '$output_formfac_file'..."
    #set avg_formfac [vecscale [expr 1.0/$formfactor_counter] $avg_formfac]
    #set avg_formfac_polymer [vecscale [expr 1.0/$formfactor_counter] $avg_formfac_polymer]
    set formfactor_file [open "$output_formfac_path$output_formfac_file" "w"]
    puts $formfactor_file "\#q\tS_1(q)\tS_1_polymer(q)"
    #foreach q $formfac_q_list sq $formfac_sq_list  ######################################### IMPORTANT BUG FROM PREVIOUS SIMULATIONS FIXED ######################################################
    foreach q $formfac_q_list sq $avg_formfac sq_polymer $avg_formfac_polymer {
        puts $formfactor_file "$q\t$sq\t$sq_polymer"
    }
    close $formfactor_file
    puts "    Done."
}
