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

}

foreach m [molinfo list] {
    l $m
}

spread

