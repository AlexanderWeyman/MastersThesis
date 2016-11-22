    # counter variable for the number of charged monomers
    set N_charged_monomers 0
    
    #set box_l_xy [expr ($bond_l*$monomer_density)**(-0.5)]
    #set box_l_z [expr $MPC*$bond_l]
    
    set first_monomer_id [setmd n_part]
    
    #set the first particle of the charged rod
    part [setmd n_part] pos 0.0 0.0 [expr ($bond_l - $box_l_z)/2.0] type $type_nM fix
    
    for {set current_z [expr (3.0*$bond_l - $box_l_z)/2.0]} {$current_z < [expr $box_l_z/2.0]} {set current_z [expr $current_z+$bond_l]} {
        #set up neutral fixed rod with no side beads
        part [setmd n_part] pos 0.0 0.0 $current_z type $type_nM bond $type_FENE [expr ([setmd n_part]-1)] fix
    }
    
    
    #change neutral to charged particle
    for {set current_particle_id [expr $first_monomer_id+$start_particle_id]} {$current_particle_id < [expr $first_monomer_id+$MPC]} {incr current_particle_id $side_bead_distance} {
        part $current_particle_id type $type_cM q $val_cM
        incr N_charged_monomers
    }
    
    
    #create counterions
    counterions $N_charged_monomers start [setmd n_part] mode $counterion_mode $shield $max_try charge $val_CI type $type_CI