    # counter variable for the number of charged monomers
    set N_charged_monomers 0
    
    set N_side_monomers $createSidechains
    
    #set N_side_chains_per_polymer [expr int(ceil(double($MPC - $start_particle_id)/$side_bead_distance))]
    #puts "Side chains per Polymer: $N_side_chains_per_polymer"

    for {set current_polymer 0} {$current_polymer < $N_P} {incr current_polymer} {
        #create main chain of the polymer (for multiple polymers use for loop and change start 0 to start [setmd n_part])
        set first_monomer_id [setmd n_part]
        polymer 1 $MPC $bond_l start $first_monomer_id mode $poly_mode $shield $max_try types $type_nM $type_cM bond $type_FENE

        #puts "added polymer, current_polymer=$current_polymer, n_part=[setmd n_part]"
        
        #counts the number of side chains (-1) beeing attachted to the backbone
        set side_chain_counter 0
        
        #current_particle_id runs over the backbone monomer id's where a side chain should be attached
        for {set current_particle_id [expr $first_monomer_id+$start_particle_id]} {$current_particle_id < [expr $first_monomer_id+$MPC]} {incr current_particle_id $side_bead_distance} {
            set current_particle_id_charged $current_particle_id
            
            #sub-loop for attaching side monomers
            for {set side_monomer_counter 0} {$side_monomer_counter < $N_side_monomers} {incr side_monomer_counter} {
                set pos_current_bead [part $current_particle_id print pos]
                
                #positions of the backbone monomer where the current side chain should be attached
                set pos_current_bead_x [lindex $pos_current_bead 0]
                set pos_current_bead_y [lindex $pos_current_bead 1]
                set pos_current_bead_z [lindex $pos_current_bead 2]
                
                #try to find a coordinate with distance $bond_l_side from the current backbone bead and a shield radius of $shield
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
                        #part [setmd n_part] pos $pos_current_bead_side_x $pos_current_bead_side_y $pos_current_bead_side_z type $type_nSM bond $type_FENE $current_particle_id_charged
                        set new_side_bead_id [expr $MPC+$side_monomer_counter*$N_side_chains_per_polymer+$side_chain_counter+$first_monomer_id]
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
    }

    #create counterions
    counterions $N_charged_monomers start [setmd n_part] mode $counterion_mode $shield $max_try charge $val_CI type $type_CI
