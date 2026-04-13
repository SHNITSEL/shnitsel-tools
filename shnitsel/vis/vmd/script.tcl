proc l {molid} {
    set old_name [molinfo $molid get name]
    set new_name [string map {.xyz ""} $old_name]
    mol rename $molid $new_name
    # animate style Loop
    mol modstyle 0 $molid Licorice 0.300000 12.000000 12.000000
    mol color Name
    mol representation Licorice 0.300000 12.000000 12.000000
    mol selection all
    mol material Opaque
    mol modrep 0 0
}

proc mtransl {molid x y z} {
    mol fix $molid
    translate by $x $y $z
    mol free $molid
    translate by [expr -$x] [expr -$y] [expr -$z]
}

proc spread {{sep 0.5} {ncols 5} {vsep_ratio 1} {molids "sentinel"} {ignore "ignore"}} {
    display update off

    set mat [lindex [molinfo top get global_matrix] 0]
    set old_x [lindex [lindex $mat 0] 3]
    set old_y [lindex [lindex $mat 1] 3]
    set old_z [lindex [lindex $mat 2] 3]
    
    set old_scale [molinfo top get scale_matrix]
    set old_rotate [molinfo top get rotate_matrix]
    
    display resetview
    display projection Orthographic
    
    if {$molids eq "sentinel"} {
        set molids [molinfo list]
    }
    foreach m $molids {
        set x [expr $sep * ($m % $ncols)]
        set y [expr $vsep_ratio * $sep * ($m / $ncols)]
        mtransl $m $x $y 0
        molinfo $m set scale_matrix $old_scale
        molinfo $m set rotate_matrix $old_rotate
    }

    set mat [lindex [molinfo top get global_matrix] 0]
    set new_x [lindex [lindex $mat 0] 3]
    set new_y [lindex [lindex $mat 1] 3]
    set new_z [lindex [lindex $mat 2] 3]

    translate by [expr $old_x - $new_x] [expr $old_y - $new_y] [expr $old_z - $new_z]

    display update on
}

proc spread_gui {} {
    toplevel .spread
    label .spread.label_scale -text "Separation"
    scale .spread.entry_scale -orient horizontal -length 200 -from 0 -to 2 -resolution 0.1 \
      -variable entry_scale -command {spread $entry_scale $entry_ncols $entry_vsep sentinel}
    
    label .spread.label_vsep -text "Vertical sep."
    scale .spread.entry_vsep -orient horizontal -length 200 -from 0.1 -to 3 -resolution 0.1 \
      -variable entry_vsep -command {spread $entry_scale $entry_ncols $entry_vsep sentinel}
    .spread.entry_vsep set 1  ; # Default value
    
    label .spread.label_ncols -text "Columns"
    scale .spread.entry_ncols -orient horizontal -length 200 -from 1 -to [llength [molinfo list]] -resolution 1 \
      -variable entry_ncols -command {spread $entry_scale $entry_ncols 1 sentinel}
    
    button .spread.button_spread -text "Spread" -command {spread $entry_scale $entry_ncols $entry_vsep}

    bind .spread.entry_scale <Return> {spread $entry_scale $entry_ncols $entry_vsep}
    bind .spread.entry_ncols <Return> {spread $entry_scale $entry_ncols $entry_vsep}
    
    grid .spread.label_scale .spread.entry_scale -sticky nsew
    grid .spread.label_vsep  .spread.entry_vsep  -sticky nsew
    grid .spread.label_ncols .spread.entry_ncols -sticky nsew
    grid .spread.button_spread                   -sticky nsew -columnspan 2
}

foreach m [molinfo list] {
    l $m
}

spread_gui

