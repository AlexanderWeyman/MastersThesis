    # counter variable for the number of charged monomers
    set N_charged_monomers 0
    
    #set box_l_xy [expr ($bond_l*$monomer_density)**(-0.5)]
    #set box_l_z [expr $MPC*$bond_l]
    
    set N_side_monomers $createSidechains
    
    #set N_side_chains_per_polymer [expr int(ceil(double($MPC - $start_particle_id)/$side_bead_distance))]
    
    set first_monomer_id [setmd n_part]
    
    
    ### set up neutral backbone ###
    
    #set the first particle of the charged rod
    part [setmd n_part] pos 0.0 0.0 [expr ($bond_l - $box_l_z)/2.0] type $type_nM fix
    
    for {set current_z [expr (3.0*$bond_l - $box_l_z)/2.0]} {$current_z < [expr $box_l_z/2.0]} {set current_z [expr $current_z+$bond_l]} {
        #set up neutral fixed rod with no side beads
        part [setmd n_part] pos 0.0 0.0 $current_z type $type_nM bond $type_FENE [expr ([setmd n_part]-1)] fix
    }
    
    
    ### attach side beads ###
    #counts the number of side chains (-1) beeing attachted to the backbone
    set side_chain_counter 0
    
    #pick particle ids from neutral backbone that should be modified
    for {set current_particle_id [expr $first_monomer_id+$start_particle_id]} {$current_particle_id < [expr $first_monomer_id+$MPC]} {incr current_particle_id $side_bead_distance} {
        set current_particle_id_charged $current_particle_id
        
        
        #sub-loop for attaching side monomers
        for {set side_monomer_counter 0} {$side_monomer_counter < $N_side_monomers} {incr side_monomer_counter} {
            set pos_current_bead [part $current_particle_id_charged print pos]

            set pos_current_bead_x [lindex $pos_current_bead 0]
            set pos_current_bead_y [lindex $pos_current_bead 1]
            set pos_current_bead_z [lindex $pos_current_bead 2]
            
            
            ### try to attach neutral side monomer outside given shield radius ###
            
            set SC_shield_counter 0
            #the error variable is set to 0 after successful addition of side monomer
            set SC_shield_error 1
            while {$SC_shield_counter < $max_try} {
                set rand_phi_current_bead [expr ($phi_max-$phi_min)*[t_random]+$phi_min]
                set rand_theta_current_bead [expr ($theta_max-$theta_min)*[t_random]+$theta_min]
                
                set pos_current_bead_side_x [expr $pos_current_bead_x+$bond_l_side*cos($rand_phi_current_bead)*sin($rand_theta_current_bead)]
                set pos_current_bead_side_y [expr $pos_current_bead_y+$bond_l_side*sin($rand_phi_current_bead)*sin($rand_theta_current_bead)]
                set pos_current_bead_side_z [expr $pos_current_bead_z+$bond_l_side*cos($rand_theta_current_bead)]
                
                #check shield radius
                if {[analyze distto $pos_current_bead_side_x $pos_current_bead_side_y $pos_current_bead_side_z] > $shield} {
                    #found suitable position to place side bead
                    set new_side_bead_id [expr $MPC+$side_monomer_counter*$N_side_chains_per_polymer+$side_chain_counter+$first_monomer_id]
                    #part [setmd n_part] pos $pos_current_bead_side_x $pos_current_bead_side_y $pos_current_bead_side_z type $type_nSM bond $type_FENE $current_particle_id_charged
                    part $new_side_bead_id pos $pos_current_bead_side_x $pos_current_bead_side_y $pos_current_bead_side_z type $type_nSM bond $type_FENE $current_particle_id_charged
                    
                    #puts "New Particle with ID '$new_side_bead_id' connected with Polymer Particle with ID '$current_particle_id_charged'."
                    
                    set SC_shield_error 0
                    #set current_particle_id_charged [expr ([setmd n_part]-1)]
                    set current_particle_id_charged $new_side_bead_id
                    break
                }
                incr SC_shield_counter
            }
            #print error and exit in case the side monomer could not be added successfully in a $shield radius within $max_try attempts
            if {$SC_shield_error} {
                puts "Error: Failed to add side monomer with given shield radius and maximum attempts."
                exit
            }            
        }
        
        
        ### change monomer type and charge of last (side) monomer ###
        
        #different types for charged side monomer and charged backbone monomer
        if {$N_side_monomers} {
            part $current_particle_id_charged type $type_cSM q $val_cM
        } else {
            part $current_particle_id_charged type $type_cM q $val_cM
        }
        incr N_charged_monomers
        incr side_chain_counter
    }
    
    
    #create counterions
    counterions $N_charged_monomers start [setmd n_part] mode $counterion_mode $shield $max_try charge $val_CI type $type_CI
