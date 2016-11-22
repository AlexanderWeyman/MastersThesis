    # counter variable for the number of charged monomers
    set N_charged_monomers 0

    for {set current_polymer 0} {$current_polymer < $N_P} {incr current_polymer} {
        #create main chain of the polymer (for multiple polymers use for loop and change start 0 to start [setmd n_part])
        set first_monomer_id [setmd n_part]
        polymer 1 $MPC $bond_l start $first_monomer_id mode $poly_mode $shield $max_try types $type_nM $type_cM bond $type_FENE

        #puts "added polymer, current_polymer=$current_polymer, n_part=[setmd n_part]"

        for {set current_particle_id [expr $first_monomer_id+$start_particle_id]} {$current_particle_id < [expr $first_monomer_id+$MPC]} {incr current_particle_id $side_bead_distance} {
            #puts $current_particle_id

            set pos_current_bead [part $current_particle_id print pos]
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
                    part [setmd n_part] pos $pos_current_bead_side_x $pos_current_bead_side_y $pos_current_bead_side_z type $type_cSM bond $type_FENE $current_particle_id q $val_cM
                    incr N_charged_monomers
                    set SC_shield_error 0
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
    }

    #create counterions
    counterions $N_charged_monomers start [setmd n_part] mode $counterion_mode $shield $max_try charge $val_CI type $type_CI