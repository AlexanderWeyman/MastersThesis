    polymer $N_P $MPC $bond_l start [setmd n_part] mode $poly_mode $shield $max_try charge $val_cM distance $side_bead_distance types $type_nM $type_cM bond $type_FENE

    #create counterions
    counterions [part gc number $type_cM] start [setmd n_part] mode $counterion_mode $shield $max_try charge $val_CI type $type_CI