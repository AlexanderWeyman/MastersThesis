    # counter variable for the number of charged monomers
    set N_charged_monomers 0
    
    set N_side_monomers $createSidechains
    
    #set N_side_chains_per_polymer [expr int(ceil(double($MPC - $start_particle_id)/$side_bead_distance))]
    #puts "Side chains per Polymer: $N_side_chains_per_polymer"
    
    set mindist_backbone_wall [expr $N_side_monomers*$bond_l_side+$shield]
    set min_z_BB $mindist_backbone_wall
    set max_z_BB [expr [lindex [setmd box_l] 2] - $mindist_backbone_wall]
    
    for {set current_polymer 0} {$current_polymer < $N_P} {incr current_polymer} {
        #create main chain of the polymer (for multiple polymers use for loop and change start 0 to start [setmd n_part])
        puts "Inserting Polymer Chain [expr $current_polymer+1] / $N_P..."
        set first_monomer_id [setmd n_part]
        set polymer_id_counter $first_monomer_id
        #rewritten polymer command with additional contraints distance check
        
        
        
        #check if box is large enough to comply with minimum distances
        if {[expr 2*$mindist_backbone_wall] > [lindex [setmd box_l] 2]} {
            puts "Error: Box is too small to place constraint."
            exit
        }
        
        #place first particle
        set BB_shield_counter 0
        set BB_shield_error 1
        while {$BB_shield_counter < $max_try} {
            set testpos_x [expr [lindex [setmd box_l] 0]*[t_random] ]
            set testpos_y [expr [lindex [setmd box_l] 1]*[t_random] ]
            set testpos_z [expr ($max_z_BB-$min_z_BB)*[t_random]+$min_z_BB]
            
            #check shield radius
            if {[analyze distto $testpos_x $testpos_y $testpos_z] > $shield || [setmd n_part] == 0} {
                #insert first bead of current polymer chain
                part $polymer_id_counter pos $testpos_x $testpos_y $testpos_z type $type_nM
                set BB_shield_error 0
                break
            }
            incr BB_shield_counter
        }
        if {$BB_shield_error} {
            puts "Error: Failed to add first backbone bead with given shield radius and maximum attempts."
            exit
        }
    
        incr polymer_id_counter
        
        while {$polymer_id_counter < [expr $first_monomer_id+$MPC]} {
            #try to place and link next polymer bead
            set BB_shield_counter 0
            set BB_shield_error 1
            while {$BB_shield_counter < $max_try} {
                set rand_phi_current_bead [expr ($phi_max-$phi_min)*[t_random]+$phi_min]
                set rand_theta_current_bead [expr ($theta_max-$theta_min)*[t_random]+$theta_min]
                
                set pos_last_bead [part [expr $polymer_id_counter-1] print pos]
                set last_x [lindex $pos_last_bead 0]
                set last_y [lindex $pos_last_bead 1]
                set last_z [lindex $pos_last_bead 2]
                
                set testpos_x [expr $last_x+$bond_l*cos($rand_phi_current_bead)*sin($rand_theta_current_bead)]
                set testpos_y [expr $last_y+$bond_l*sin($rand_phi_current_bead)*sin($rand_theta_current_bead)]
                set testpos_z [expr $last_z+$bond_l*cos($rand_theta_current_bead)]
                
                #check shield radius
                if {[analyze distto $testpos_x $testpos_y $testpos_z] > $shield && $testpos_z > $min_z_BB && $testpos_z < $max_z_BB} {
                    #insert first bead of current polymer chain
                    part $polymer_id_counter pos $testpos_x $testpos_y $testpos_z type $type_nM bond $type_FENE [expr $polymer_id_counter-1]
                    set BB_shield_error 0
                    break
                }
                incr BB_shield_counter
            }
            if {$BB_shield_error} {
                puts "Error: Failed to add backbone bead with given shield radius and maximum attempts. No of last polymer bead added: $polymer_id_counter."
                exit
            }
            
            incr polymer_id_counter
            
        }

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
    #counterions $N_charged_monomers start [setmd n_part] mode $counterion_mode $shield $max_try charge $val_CI type $type_CI
    
    set mindist_ci_wall [expr $bond_l_side+$shield]
    set min_z_CI $mindist_ci_wall
    set max_z_CI [expr [lindex [setmd box_l] 2] - $mindist_ci_wall]
    
    set ci_counter 0
    while {$ci_counter < $N_charged_monomers} {
        #place first particle
        set CI_shield_counter 0
        set CI_shield_error 1
        while {$CI_shield_counter < $max_try} {
            set testpos_x [expr [lindex [setmd box_l] 0]*[t_random] ]
            set testpos_y [expr [lindex [setmd box_l] 1]*[t_random] ]
            set testpos_z [expr ($max_z_CI-$min_z_CI)*[t_random]+$min_z_CI]
            
            #check shield radius
            if {[analyze distto $testpos_x $testpos_y $testpos_z] > $shield} {
                #insert first bead of current polymer chain
                part [setmd n_part] pos $testpos_x $testpos_y $testpos_z type $type_CI q $val_CI
                set CI_shield_error 0
                break
            }
            incr CI_shield_counter
        }
        if {$CI_shield_error} {
            puts "Error: Failed to add counterion with given shield radius and maximum attempts."
            exit
        }
    
        incr ci_counter
    }
    
