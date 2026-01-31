def parse_molcas_log(f):
    geoms = []
    for line in f:
        # fmt: off
        ## Constraint values
            # *************************************************************
            # * Values of the primitive constraints                       *
            # *************************************************************
            # A        : Energy difference =         0.00166157 hartree,         4.36244406 kJ/mol
            #         Average energy    =      -232.87685409 hartree
            # NADC     : H12               =         0.00000000 hartree
            #
            #
        # fmt: on
        if line.startswith(" * Values of the primitive constraints"):
            ...

        # fmt: off
        ## Energies (only parse once?)
            #**********************************************************************************************************************
            #*                                    Energy Statistics for Geometry Optimization                                     *
            #**********************************************************************************************************************
            #                       Energy     Grad      Grad              Step                 Estimated   Geom       Hessian     
            #Iter      Energy       Change     Norm      Max    Element    Max     Element     Final Energy Update Update   Index  
            #  1   -232.77172162  0.00000000 0.207765 -0.082263 nrc001  -0.234823* nrc001     -232.81168205 RS-RFO  None      0    
            #  2   -232.83465949 -0.06293787 0.118016  0.075937 nrc004   0.232694* nrc010     -232.85821976 RS-RFO  BFGS      0    
        ## etc.
            #
            #        +----------------------------------+----------------------------------+
            #        +    Cartesian Displacements       +    Gradient in internals         +
            #        +  Value      Threshold Converged? +  Value      Threshold Converged? +
            #  +-----+----------------------------------+----------------------------------+
            #  + RMS + 5.2941E-02  1.2000E-03     No    + 2.5018E-03  3.0000E-04     No    +
            #  +-----+----------------------------------+----------------------------------+
            #  + Max + 8.9838E-02  1.8000E-03     No    + 8.0092E-03  4.5000E-04     No    +
            #  +-----+----------------------------------+----------------------------------+

            #  Convergence not reached yet!

            # *****************************************************************************************************************
        # fmt: on
        if (
            line.startswith('* ')
            and ' '.join(line.strip().split()[1:-1])
            == "Energy Statistics for Geometry Optimization"
        ):
            ...

        # fmt: off
        ## Geometries
            #++ Geometry section
            #
            #********************************************************************************
            #  Geometrical information of the new structure
            #********************************************************************************
            #
            #
            # *********************************************************
            # * Nuclear coordinates for the next iteration / Bohr     *
            # *********************************************************
            #  ATOM              X               Y               Z     
            #  C1              -3.489709        0.192996        0.016320
            #  H2              -2.387903       -1.393580       -0.661536
            #  C3              -2.434606        1.898692        1.821819
        ## etc.
            # *********************************************************
            # * Nuclear coordinates for the next iteration / Angstrom *
            # *********************************************************
            #  ATOM              X               Y               Z     
            #  C1              -1.846674        0.102129        0.008636
            #  H2              -1.263624       -0.737451       -0.350070
            #  C3              -1.288338        1.004745        0.964065
        ## etc.
            #--
        # fmt: on
        if line.startswith("++ Geometry section"):
            while (
                not (line := next(f)).startswith(" * Nuclear coordinates")
                or "Angstrom" not in line
            ):
                pass

            next(f)  # line of stars
            next(f)  # header line

            geom = []
            while (stripline := next(f).strip()) != "":
                geom.append(list(map(float, stripline[1:])))
        geoms.append(geom)
    return geoms