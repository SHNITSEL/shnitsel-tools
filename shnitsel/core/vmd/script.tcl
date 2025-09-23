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

proc spread {{scale 0.5} {ncols 5} {molids "sentinel"}} {
    display resetview
    display projection Orthographic
    if {$molids eq "sentinel"} {
        set molids [molinfo list]
    }
    foreach m $molids {
        puts $m
        set scale $scale
        set x [expr $scale * ($m % $ncols)]
        set y [expr $scale * ($m / $ncols)]
        mtransl $m $x $y 0
    }
    translate by [expr $scale * ($ncols - 1) / 2.0] 0 0
}

proc spread_gui {} {
    toplevel .spread
    label .spread.label_scale -text "Separation"
    entry .spread.entry_scale -textvariable entry_scale
    label .spread.label_ncols -text "Columns"
    entry .spread.entry_ncols -textvariable entry_ncols
    button .spread.button_spread -text "Spread" -command {spread $entry_scale $entry_ncols}

    bind .spread.entry_scale <Return> {spread $entry_scale $entry_ncols}
    bind .spread.entry_ncols <Return> {spread $entry_scale $entry_ncols}
    
    grid .spread.label_scale .spread.entry_scale -sticky nsew
    grid .spread.label_ncols .spread.entry_ncols -sticky nsew
    grid .spread.button_spread                   -sticky nsew -columnspan 2
}

foreach m [molinfo list] {
    l $m
}

spread_gui

