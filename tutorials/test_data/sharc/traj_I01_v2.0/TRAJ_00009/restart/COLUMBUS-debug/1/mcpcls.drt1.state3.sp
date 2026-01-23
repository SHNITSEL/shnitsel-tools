

     ******************************************
     **    PROGRAM:              MCPC        **
     **    PROGRAM VERSION:      5.5         **
     **    DISTRIBUTION VERSION: 5.9.a       **
     ******************************************


 original author: Daniel Robertson, FSU
 later revisions: Ron Shepard, ANL;
                  Michal Dallos, University Vienna



 This Version of Program mcpc is Maintained by:
     Thomas Mueller
     Juelich Supercomputing Centre (JSC)
     Institute of Advanced Simulation (IAS)
     D-52425 Juelich, Germany 
     Email: th.mueller@fz-juelich.de



   ******  File header section  ******

 Headers form the restart file:
    Hermit Integral Program : SIFS version  polonium12        13:55:54.332 21-Feb-18
     title                                                                          


   ******  DRT info section  ******

 Informations for the DRT no.  1
 Header form the DRT file: 
     title                                                                          
 Molecular symmetry group:   sym1
 Total number of electrons:   16
 Spin multiplicity:            1
 Number of active orbitals:    4
 Number of active electrons:   6
 Total number of CSFs:        10

   ***  Informations from the DRT number:   1

 
 Symmetry orbital summary:
 Symm.blocks:         1
 Symm.labels:         a  

 List of doubly occupied orbitals:
  1 a    2 a    3 a    4 a    5 a  

 List of active orbitals:
  6 a    7 a    8 a    9 a  


   ******  MCSCF convergence information:  ******

 MCSCF convergence criteria were not satisfied.

 mcscf energy=   -93.9741423714    nuclear repulsion=    22.8629618686
 demc=             0.0000067580    wnorm=                 0.0022028292
 knorm=            0.0056993976    apxde=                 0.0000032204


 MCSCF calculation performmed for   1 symmetry.

 State averaging:
 No,  ssym, navst, wavst
  1    a      3   0.3333 0.3333 0.3333

 Input the DRT No of interest: [  1]:
In the DRT No.: 1 there are  3 states.

 Which one to take? [  1]:
 The CSFs for the state No  3 of the symmetry  a   will be printed
 according to the following print options :

 1) print csf info by sorted index number.
 2) print csf info by contribution threshold.
 3) print csf info by csf number.
 4) set additional print options.
 5) print the entire sorted csf vector.
 6) print the entire csf vector.
 7) print the mcscf molecular orbitals.
 8) print the mcscf natural orbitals and occupation numbers.
 9) export wave function files for cioverlap (all states).
 0) end.

 input menu number [  0]: csfs will be printed based on coefficient magnitudes.

 input the coefficient threshold (end with 0.) [ 0.0000]:
 List of active orbitals:
  6 a    7 a    8 a    9 a  

   csf       coeff       coeff**2    step(*)
  -----  ------------  ------------  ------------
      1 -0.9915948736  0.9832603934  3330
     10  0.1018015662  0.0103635589  0333
      6  0.0640207430  0.0040986555  3033
      2 -0.0343746498  0.0011816165  3312
      9 -0.0321409134  0.0010330383  1233
      5 -0.0055355495  0.0000306423  3123
      8 -0.0052793736  0.0000278718  1323
      7 -0.0019441119  0.0000037796  1332

 input the coefficient threshold (end with 0.) [ 0.0000]:
 1) print csf info by sorted index number.
 2) print csf info by contribution threshold.
 3) print csf info by csf number.
 4) set additional print options.
 5) print the entire sorted csf vector.
 6) print the entire csf vector.
 7) print the mcscf molecular orbitals.
 8) print the mcscf natural orbitals and occupation numbers.
 9) export wave function files for cioverlap (all states).
 0) end.

 input menu number [  0]: