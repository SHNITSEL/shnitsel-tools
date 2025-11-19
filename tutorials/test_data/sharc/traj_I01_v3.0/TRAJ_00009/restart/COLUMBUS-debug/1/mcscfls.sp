

     ******************************************
     **    PROGRAM:              MCSCF       **
     **    PROGRAM VERSION:      5.5         **
     **    DISTRIBUTION VERSION: 5.9.a       **
     ******************************************

 This program allows the csf mixing coefficient and orbital expansion coefficient
 optimization using the graphical unitary group approach and the exponential
 operator mcscf method.
 references:  r. shepard and j. simons, ' int. j. quantum chem. symp. 14, 211 (1980).
              r. shepard, i. shavitt, and j. simons, j. chem. phys. 76, 543 (1982).
              r. shepard in "ab initio methods in quantum chemistry ii" advances in chemical
                  physics 69, edited by k. p. lawley (wiley, new york, 1987) pp. 63-200.
 Original autor: Ron Shepard, ANL
 Later revisions: Michal Dallos, University Vienna

 This Version of Program MCSCF is Maintained by:
     Thomas Mueller
     Juelich Supercomputing Centre (JSC)
     Institute of Advanced Simulation (IAS)
     D-52425 Juelich, Germany 
     Email: th.mueller@fz-juelich.de



     ******************************************
     **    PROGRAM:              MCSCF       **
     **    PROGRAM VERSION:      5.4.0.2     **
     **    DISTRIBUTION VERSION: 5.9.a       **
     ******************************************

 Workspace allocation information:
       393216000 of real*8 words ( 3000.00 MB) of work space has been allocated.

 user input information:

 ======== echo of the mcscf input ========
 ------------------------------------------------------------------------
  &input
   niter=100,
   nmiter=50,
   nciitr=300,
   tol(3)=1.e-6,
   tol(2)=1.e-6,
   tol(1)=1.e-10,
   NSTATE=0,
   npath=1,3,9,10,13,17,19,21,-11,12, 2,
   ncoupl=5,
   tol(9)=1.e-3,
   FCIORB=  1,6,20,1,7,20,1,8,20,1,9,20
   NAVST(1) = 3,
   WAVST(1,1)=1 ,
   WAVST(1,2)=1 ,
   WAVST(1,3)=1 ,
  &end
 ------------------------------------------------------------------------


 ***  Integral file informations  ***


 input integral file : /public/polonium12/scratch/tmp/anna.6654/Singlet_2/TRAJ_0
 000

 Integral file header information:
 Hermit Integral Program : SIFS version  polonium12        13:55:54.332 21-Feb-18

 Core type energy values:
 energy( 1)=  2.286296186856E+01, ietype=   -1,    core energy of type: Nuc.Rep.
 total ao core energy =   22.862961869


   ******  Basis set information:  ******

 Number of irreps:                  1
 Total number of basis functions:  82

 irrep no.              1
 irrep label           A  
 no. of bas.fcions.    82
 warning: resetting tol(6) due to possible inconsistency.


 ***  MCSCF optimization procedure parmeters:  ***


 maximum number of mcscf iterations:        niter=   100

 maximum number of psci micro-iterations:   nmiter=   50
 maximum r,s subspace dimension allowed:    nvrsmx=   30

 tol(1)=  1.0000E-10. . . . delta-emc convergence criterion.
 tol(2)=  1.0000E-06. . . . wnorm convergence criterion.
 tol(3)=  1.0000E-06. . . . knorm convergence criterion.
 tol(4)=  1.0000E-08. . . . apxde convergence criterion.
 tol(5)=  1.0000E-04. . . . small diagonal matrix element tolerance.
 tol(6)=  1.2500E-07. . . . minimum ci-psci residual norm.
 tol(7)=  1.0000E-05. . . . maximum ci-psci residual norm.
 tol(8)=  1.0000E+00. . . . maximum abs(k(xy)) allowed.
 tol(9)=  1.0000E-03. . . . wnorm coupling tolerance.
 tol(10)= 0.0000E+00. . . . maximum psci emergency shift parameter.
 tol(11)= 0.0000E+00. . . . minimum psci emergency shift parameter.
 tol(12)= 0.0000E+00. . . . increment of psci emergency shift parameter.


 *** State averaging informations: ***


 MCSCF calculation performed for  1 DRT.

 DRT  first state   no.of aver.states   weights
  1   ground state          3             0.333 0.333 0.333

 The number of hmc(*) eigenvalues and eigenvectors calculated each iteration per DRT:
 DRT.   no.of eigenv.(=ncol)
    1        4

 orbital coefficients are optimized for the ground state (nstate=0).

 Orbitals included in invariant subspaces:
   symmetry   orbital   mask
       1       6(  6)    20
       1       7(  7)    20
       1       8(  8)    20
       1       9(  9)    20

 npath(*) options:
  2:  orbital-state coupling terms will be included beginning on iteration ncoupl=  5
  3:  print intermediate timing information.
  9:  suppress the drt listing.
 10:  suppress the hmc(*) eigenvector listing.
 12:  diagonalize the hmc(*) matrix iteratively.
        nunitv= 1 nciitr=** mxvadd=20 nvcimx=20
       rtolci(*),wnorm=     1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02
                            1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02
                            1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02
                            1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02
                            1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02
                            1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02
                            1.0000E-02 1.0000E-02 0.0000E+00
   noldv =   0
 13:  get initial orbitals from the formatted file, mocoef.
 17:  print the final natural orbitals and occupations.
 19:  transform the virtual orbitals to diagonalize qvv(*).
 21:  write out the one- and two- electron density for further use (files:mcd1fl, mcd2fl).


   ******  DRT info section  ******


 Informations for the DRT no.  1

 DRT file header:
  title                                                                          
 Molecular symmetry group:    a  
 Total number of electrons:   16
 Spin multiplicity:            1
 Number of active orbitals:    4
 Number of active electrons:   6
 Total number of CSFs:        10
 

 faar:   0 active-active rotations allowed out of:   6 possible.


 Number of active-double rotations:        20
 Number of active-active rotations:         0
 Number of double-virtual rotations:      365
 Number of active-virtual rotations:      292
 lenbfsdef=                131071  lenbfs=                  5329
  number of integrals per class 1:11 (cf adda 
 class  1 (pq|rs):         #          55
 class  2 (pq|ri):         #         200
 class  3 (pq|ia):         #        3650
 class  4 (pi|qa):         #        5840
 class  5 (pq|ra):         #        2920
 class  6 (pq|ij)/(pi|qj): #         400
 class  7 (pq|ab):         #       27010
 class  8 (pa|qb):         #       53290
 class  9 p(bp,ai)         #      106580
 class 10p(ai,jp):        #        7300
 class 11p(ai,bj):        #       79935

 Size of orbital-Hessian matrix B:                   253195
 Size of the orbital-state Hessian matrix C:          20310
 Total size of the state Hessian matrix M:                0
 Size of HESSIAN-matrix for quadratic conv.:         273505


 Source of the initial MO coeficients:

 Input MO coefficient file: /public/polonium12/scratch/tmp/anna.6654/Singlet_2/TRAJ_0000
 

               starting mcscf iteration...   1

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368297 transformed 1/r12    array elements were written in      68 records.

 !timer: 2-e transformation              cpu_time=     0.757 walltime=     0.757

 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

 trial vectors are generated internally.

 trial vector  1 is unit matrix column     2
 ciiter=   5 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0792139062     -116.9421757747        0.0000000000        0.0000001250
    2       -93.9333374345     -116.7962993030        0.0000000000        0.0000001250
    3       -93.9041754386     -116.7671373071        0.0000000000        0.0000001250
    4       -93.7609833952     -116.6239452638        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  1.199323187650188E-002
 Total number of micro iterations:   10

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.99824710 pnorm= 0.0000E+00 rznorm= 2.6993E-06 rpnorm= 0.0000E+00 noldr= 10 nnewr= 10 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.713966  -23.228275   -2.325267   -1.855759   -1.562641

 qvv(*) eigenvalues. symmetry block  1
    -0.422058   -0.199050   -0.162104   -0.123176   -0.084965   -0.042602   -0.025253    0.023203    0.041877    0.065654
     0.082321    0.094262    0.140887    0.200268    0.222244    0.275147    0.322558    0.393705    0.420612    0.451377
     0.476478    0.496502    0.523288    0.559267    0.616130    0.644956    0.665237    0.728686    0.804100    0.820099
     0.908212    0.923383    0.955280    0.966156    1.004424    1.193720    1.275763    1.305708    1.324371    1.335047
     1.356992    1.585823    1.792349    1.896172    1.912869    2.033240    2.482323    2.529065    2.605925    2.773638
     2.778912    2.854493    2.864435    3.051454    3.112204    3.303904    3.392160    3.415811    3.490237    3.808272
     4.181979    4.251899    4.256664    4.332558    4.507785    4.748788    4.852146    5.000047    5.032357    5.233748
     5.934556    6.224147    6.373215

 restrt: restart information saved on the restart file (unit= 13).
 !timer: mcscf iteration                 cpu_time=     0.824 walltime=     0.824

 not all mcscf convergence criteria are satisfied.
 iter=    1 emc=    -93.9722422597 demc= 9.3972E+01 wnorm= 9.5946E-02 knorm= 5.9184E-02 apxde= 1.0382E-03    *not conv.*     

               starting mcscf iteration...   2
 !timer:                                 cpu_time=     0.832 walltime=     0.913

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368298 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0806109071     -116.9435727756        0.0000000000        0.0000100000
    2       -93.9345838773     -116.7975457459        0.0000000000        0.0000100000
    3       -93.9047482275     -116.7677100961        0.0000000000        0.0000100000
    4       -93.7618199779     -116.6247818464        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  4.288970866312067E-004
 Total number of micro iterations:   10

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.99795288 pnorm= 0.0000E+00 rznorm= 1.0655E-06 rpnorm= 0.0000E+00 noldr= 10 nnewr= 10 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.713985  -23.224186   -2.322688   -1.855130   -1.564720

 qvv(*) eigenvalues. symmetry block  1
    -0.420543   -0.199313   -0.161883   -0.123265   -0.084915   -0.042106   -0.025316    0.022940    0.041693    0.066279
     0.082501    0.094172    0.141001    0.199781    0.221868    0.275681    0.322618    0.393692    0.420226    0.451223
     0.476546    0.496954    0.524205    0.559756    0.616876    0.644830    0.664839    0.729703    0.804026    0.819595
     0.908158    0.923390    0.955466    0.966395    1.004466    1.193365    1.276081    1.305998    1.324205    1.334718
     1.357382    1.587291    1.793884    1.895500    1.912966    2.032885    2.482492    2.530580    2.606811    2.775317
     2.779737    2.855193    2.866830    3.052490    3.112364    3.301965    3.393038    3.417015    3.491712    3.806888
     4.181608    4.251845    4.258428    4.336105    4.507720    4.748952    4.851241    5.000484    5.032668    5.236424
     5.933149    6.227700    6.374422

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=    2 emc=    -93.9733143373 demc= 1.0721E-03 wnorm= 3.4312E-03 knorm= 6.3954E-02 apxde= 1.5288E-05    *not conv.*     

               starting mcscf iteration...   3
 !timer:                                 cpu_time=     1.154 walltime=     1.235

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368295 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0806324721     -116.9435943407        0.0000000000        0.0000001250
    2       -93.9346849122     -116.7976467807        0.0000000000        0.0000001250
    3       -93.9047002926     -116.7676621612        0.0000000000        0.0000001250
    4       -93.7631241256     -116.6260859942        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  3.485041672745654E-004
 Total number of micro iterations:   10

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.99931015 pnorm= 0.0000E+00 rznorm= 7.4701E-07 rpnorm= 0.0000E+00 noldr= 10 nnewr= 10 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.713459  -23.223996   -2.322500   -1.854966   -1.585501

 qvv(*) eigenvalues. symmetry block  1
    -0.420374   -0.199305   -0.161935   -0.123283   -0.084902   -0.041933   -0.025262    0.023014    0.041762    0.066345
     0.082561    0.094202    0.141231    0.199767    0.221745    0.275709    0.322649    0.393740    0.420083    0.451124
     0.476384    0.497091    0.524349    0.559843    0.617011    0.644914    0.665019    0.729813    0.803908    0.819304
     0.908057    0.923376    0.955564    0.966485    1.004457    1.193728    1.276250    1.306319    1.324355    1.334908
     1.357583    1.587354    1.794204    1.896137    1.912868    2.033145    2.482555    2.530775    2.606945    2.775465
     2.779925    2.855297    2.866887    3.052086    3.112275    3.300314    3.393172    3.417457    3.491952    3.805819
     4.181725    4.251963    4.258585    4.336276    4.507291    4.748884    4.850950    5.001596    5.033070    5.236582
     5.932278    6.227870    6.375950

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=    3 emc=    -93.9733392256 demc= 2.4888E-05 wnorm= 2.7880E-03 knorm= 3.7138E-02 apxde= 5.7822E-06    *not conv.*     

               starting mcscf iteration...   4
 !timer:                                 cpu_time=     1.476 walltime=     1.557

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368298 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0806740461     -116.9436359146        0.0000000000        0.0000001250
    2       -93.9347305617     -116.7976924303        0.0000000000        0.0000001250
    3       -93.9046466859     -116.7676085545        0.0000000000        0.0000001250
    4       -93.7639001746     -116.6268620431        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  2.586113418762733E-004
 Total number of micro iterations:   11

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.99960403 pnorm= 0.0000E+00 rznorm= 1.1486E-07 rpnorm= 0.0000E+00 noldr= 11 nnewr= 11 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.713651  -23.223912   -2.322418   -1.854888   -1.599979

 qvv(*) eigenvalues. symmetry block  1
    -0.420290   -0.199310   -0.161961   -0.123296   -0.084910   -0.041881   -0.025201    0.023027    0.041796    0.066386
     0.082576    0.094196    0.141323    0.199744    0.221674    0.275747    0.322643    0.393760    0.420016    0.451070
     0.476295    0.497163    0.524402    0.559876    0.617064    0.644979    0.665097    0.729841    0.803835    0.819138
     0.908011    0.923331    0.955601    0.966527    1.004426    1.193876    1.276302    1.306466    1.324362    1.334981
     1.357660    1.587373    1.794353    1.896394    1.912764    2.033242    2.482569    2.530862    2.607006    2.775532
     2.780009    2.855367    2.866913    3.051822    3.112171    3.299425    3.393290    3.417748    3.492008    3.805210
     4.181755    4.252032    4.258658    4.336334    4.507042    4.748830    4.850690    5.002133    5.033226    5.236650
     5.931725    6.227950    6.376735

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=    4 emc=    -93.9733504312 demc= 1.1206E-05 wnorm= 2.0689E-03 knorm= 2.8139E-02 apxde= 4.2018E-06    *not conv.*     

               starting mcscf iteration...   5
 !timer:                                 cpu_time=     1.800 walltime=     1.880

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368298 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0807018561     -116.9436637246        0.0000000000        0.0000001250
    2       -93.9347686383     -116.7977305069        0.0000000000        0.0000001250
    3       -93.9046056029     -116.7675674714        0.0000000000        0.0000001250
    4       -93.7645201940     -116.6274820625        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  2.259973115991537E-004
 Total number of micro iterations:   10

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.99973175 pnorm= 0.0000E+00 rznorm= 3.8896E-07 rpnorm= 0.0000E+00 noldr= 10 nnewr= 10 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.713725  -23.223835   -2.322343   -1.854815   -1.611616

 qvv(*) eigenvalues. symmetry block  1
    -0.420213   -0.199313   -0.161985   -0.123303   -0.084920   -0.041833   -0.025147    0.023046    0.041828    0.066423
     0.082583    0.094194    0.141402    0.199733    0.221610    0.275784    0.322635    0.393772    0.419957    0.451020
     0.476219    0.497215    0.524450    0.559908    0.617112    0.645039    0.665171    0.729869    0.803778    0.819001
     0.907976    0.923294    0.955637    0.966567    1.004396    1.194003    1.276358    1.306595    1.324374    1.335048
     1.357730    1.587399    1.794479    1.896614    1.912702    2.033329    2.482575    2.530944    2.607062    2.775593
     2.780086    2.855431    2.866938    3.051594    3.112077    3.298700    3.393365    3.417960    3.492062    3.804709
     4.181781    4.252099    4.258727    4.336388    4.506817    4.748764    4.850508    5.002571    5.033343    5.236714
     5.931312    6.228023    6.377356

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=    5 emc=    -93.9733586991 demc= 8.2679E-06 wnorm= 1.8080E-03 knorm= 2.3161E-02 apxde= 3.4072E-06    *not conv.*     

               starting mcscf iteration...   6
 !timer:                                 cpu_time=     2.123 walltime=     2.203

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368297 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0807251968     -116.9436870654        0.0000000000        0.0000001250
    2       -93.9348008225     -116.7977626910        0.0000000000        0.0000001250
    3       -93.9045703164     -116.7675321849        0.0000000000        0.0000001250
    4       -93.7650494558     -116.6280113244        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  2.061254938170635E-004
 Total number of micro iterations:   11

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.99979658 pnorm= 0.0000E+00 rznorm= 8.6394E-08 rpnorm= 0.0000E+00 noldr= 11 nnewr= 11 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.713768  -23.223765   -2.322274   -1.854747   -1.621594

 qvv(*) eigenvalues. symmetry block  1
    -0.420144   -0.199314   -0.162006   -0.123308   -0.084930   -0.041791   -0.025099    0.023067    0.041858    0.066456
     0.082588    0.094193    0.141469    0.199727    0.221552    0.275818    0.322627    0.393782    0.419905    0.450974
     0.476154    0.497256    0.524493    0.559938    0.617156    0.645093    0.665239    0.729895    0.803729    0.818882
     0.907948    0.923260    0.955670    0.966603    1.004367    1.194113    1.276410    1.306709    1.324386    1.335107
     1.357793    1.587426    1.794588    1.896806    1.912660    2.033406    2.482578    2.531018    2.607111    2.775648
     2.780156    2.855489    2.866963    3.051394    3.111988    3.298078    3.393421    3.418132    3.492112    3.804277
     4.181801    4.252161    4.258789    4.336437    4.506612    4.748697    4.850361    5.002943    5.033436    5.236774
     5.930974    6.228091    6.377875

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=    6 emc=    -93.9733654452 demc= 6.7461E-06 wnorm= 1.6490E-03 knorm= 2.0169E-02 apxde= 2.9243E-06    *not conv.*     

               starting mcscf iteration...   7
 !timer:                                 cpu_time=     2.447 walltime=     2.527

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368298 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0807458504     -116.9437077189        0.0000000000        0.0000001250
    2       -93.9348292414     -116.7977911100        0.0000000000        0.0000001250
    3       -93.9045386835     -116.7675005521        0.0000000000        0.0000001250
    4       -93.7655205482     -116.6284824168        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  1.928718184230798E-004
 Total number of micro iterations:   11

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.99983625 pnorm= 0.0000E+00 rznorm= 8.5836E-08 rpnorm= 0.0000E+00 noldr= 11 nnewr= 11 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.713797  -23.223699   -2.322210   -1.854685   -1.630622

 qvv(*) eigenvalues. symmetry block  1
    -0.420079   -0.199315   -0.162025   -0.123312   -0.084940   -0.041751   -0.025056    0.023088    0.041885    0.066487
     0.082592    0.094193    0.141528    0.199724    0.221498    0.275850    0.322620    0.393789    0.419858    0.450932
     0.476094    0.497291    0.524533    0.559967    0.617198    0.645143    0.665302    0.729919    0.803685    0.818775
     0.907923    0.923230    0.955702    0.966637    1.004339    1.194214    1.276461    1.306812    1.324398    1.335159
     1.357851    1.587454    1.794686    1.896977    1.912630    2.033474    2.482579    2.531088    2.607157    2.775700
     2.780220    2.855543    2.866986    3.051212    3.111903    3.297524    3.393464    3.418277    3.492159    3.803890
     4.181818    4.252219    4.258847    4.336483    4.506422    4.748630    4.850236    5.003272    5.033511    5.236830
     5.930683    6.228153    6.378329

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=    7 emc=    -93.9733712584 demc= 5.8132E-06 wnorm= 1.5430E-03 knorm= 1.8096E-02 apxde= 2.5977E-06    *not conv.*     

               starting mcscf iteration...   8
 !timer:                                 cpu_time=     2.770 walltime=     2.850

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368297 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0807646947     -116.9437265633        0.0000000000        0.0000001250
    2       -93.9348549864     -116.7978168550        0.0000000000        0.0000001250
    3       -93.9045096312     -116.7674714997        0.0000000000        0.0000001250
    4       -93.7659501953     -116.6289120639        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  1.833521041358599E-004
 Total number of micro iterations:   11

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.99986293 pnorm= 0.0000E+00 rznorm= 8.5044E-08 rpnorm= 0.0000E+00 noldr= 11 nnewr= 11 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.713819  -23.223637   -2.322150   -1.854626   -1.638999

 qvv(*) eigenvalues. symmetry block  1
    -0.420018   -0.199316   -0.162043   -0.123314   -0.084951   -0.041715   -0.025015    0.023109    0.041910    0.066515
     0.082595    0.094194    0.141583    0.199723    0.221448    0.275881    0.322612    0.393795    0.419814    0.450891
     0.476039    0.497321    0.524571    0.559995    0.617237    0.645191    0.665361    0.729943    0.803645    0.818677
     0.907901    0.923201    0.955732    0.966669    1.004312    1.194306    1.276510    1.306908    1.324411    1.335207
     1.357906    1.587481    1.794776    1.897134    1.912608    2.033537    2.482577    2.531153    2.607200    2.775748
     2.780282    2.855594    2.867008    3.051043    3.111820    3.297020    3.393500    3.418403    3.492204    3.803535
     4.181833    4.252275    4.258903    4.336527    4.506243    4.748562    4.850127    5.003570    5.033575    5.236884
     5.930426    6.228212    6.378735

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=    8 emc=    -93.9733764374 demc= 5.1790E-06 wnorm= 1.4668E-03 knorm= 1.6557E-02 apxde= 2.3640E-06    *not conv.*     

               starting mcscf iteration...   9
 !timer:                                 cpu_time=     3.092 walltime=     3.172

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368298 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0807822262     -116.9437440948        0.0000000000        0.0000001250
    2       -93.9348787815     -116.7978406500        0.0000000000        0.0000001250
    3       -93.9044824762     -116.7674443447        0.0000000000        0.0000001250
    4       -93.7663485973     -116.6293104658        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  1.762718856808156E-004
 Total number of micro iterations:   11

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.99988194 pnorm= 0.0000E+00 rznorm= 8.5156E-08 rpnorm= 0.0000E+00 noldr= 11 nnewr= 11 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.713838  -23.223578   -2.322092   -1.854569   -1.646892

 qvv(*) eigenvalues. symmetry block  1
    -0.419961   -0.199316   -0.162060   -0.123316   -0.084961   -0.041680   -0.024977    0.023130    0.041934    0.066542
     0.082597    0.094195    0.141633    0.199723    0.221399    0.275911    0.322605    0.393800    0.419774    0.450853
     0.475988    0.497348    0.524607    0.560021    0.617275    0.645237    0.665418    0.729966    0.803609    0.818585
     0.907881    0.923173    0.955762    0.966700    1.004286    1.194392    1.276557    1.306998    1.324423    1.335251
     1.357959    1.587509    1.794859    1.897279    1.912593    2.033595    2.482575    2.531216    2.607241    2.775794
     2.780340    2.855642    2.867030    3.050885    3.111740    3.296552    3.393530    3.418515    3.492247    3.803205
     4.181845    4.252328    4.258956    4.336569    4.506071    4.748493    4.850028    5.003844    5.033629    5.236936
     5.930195    6.228269    6.379105

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=    9 emc=    -93.9733811613 demc= 4.7238E-06 wnorm= 1.4102E-03 knorm= 1.5365E-02 apxde= 2.1909E-06    *not conv.*     

               starting mcscf iteration...  10
 !timer:                                 cpu_time=     3.412 walltime=     3.492

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368296 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0807987549     -116.9437606235        0.0000000000        0.0000001250
    2       -93.9349011085     -116.7978629771        0.0000000000        0.0000001250
    3       -93.9044567786     -116.7674186471        0.0000000000        0.0000001250
    4       -93.7667226136     -116.6296844822        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  1.709130804464135E-004
 Total number of micro iterations:   11

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.99989607 pnorm= 0.0000E+00 rznorm= 8.5959E-08 rpnorm= 0.0000E+00 noldr= 11 nnewr= 11 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.713854  -23.223521   -2.322037   -1.854514   -1.654411

 qvv(*) eigenvalues. symmetry block  1
    -0.419905   -0.199316   -0.162076   -0.123317   -0.084971   -0.041646   -0.024941    0.023152    0.041957    0.066568
     0.082598    0.094197    0.141679    0.199724    0.221353    0.275940    0.322598    0.393804    0.419735    0.450816
     0.475940    0.497373    0.524642    0.560047    0.617311    0.645280    0.665474    0.729988    0.803574    0.818499
     0.907862    0.923147    0.955790    0.966730    1.004260    1.194473    1.276603    1.307084    1.324435    1.335292
     1.358009    1.587536    1.794938    1.897415    1.912582    2.033650    2.482571    2.531276    2.607280    2.775838
     2.780396    2.855689    2.867051    3.050734    3.111660    3.296114    3.393555    3.418615    3.492289    3.802894
     4.181856    4.252380    4.259007    4.336610    4.505906    4.748425    4.849939    5.004099    5.033676    5.236986
     5.929983    6.228324    6.379446

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   10 emc=    -93.9733855473 demc= 4.3860E-06 wnorm= 1.3673E-03 knorm= 1.4417E-02 apxde= 2.0598E-06    *not conv.*     

               starting mcscf iteration...  11
 !timer:                                 cpu_time=     3.734 walltime=     3.814

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368296 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0808144940     -116.9437763626        0.0000000000        0.0000001250
    2       -93.9349223102     -116.7978841787        0.0000000000        0.0000001250
    3       -93.9044322280     -116.7673940966        0.0000000000        0.0000001250
    4       -93.7670771250     -116.6300389936        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  1.668343985751875E-004
 Total number of micro iterations:   11

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.99990688 pnorm= 0.0000E+00 rznorm= 8.7249E-08 rpnorm= 0.0000E+00 noldr= 11 nnewr= 11 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.713869  -23.223465   -2.321983   -1.854461   -1.661631

 qvv(*) eigenvalues. symmetry block  1
    -0.419851   -0.199316   -0.162092   -0.123317   -0.084981   -0.041614   -0.024906    0.023173    0.041979    0.066593
     0.082600    0.094199    0.141723    0.199727    0.221309    0.275969    0.322591    0.393807    0.419698    0.450779
     0.475893    0.497395    0.524676    0.560073    0.617347    0.645323    0.665527    0.730009    0.803542    0.818417
     0.907845    0.923121    0.955818    0.966759    1.004235    1.194550    1.276648    1.307166    1.324446    1.335330
     1.358057    1.587564    1.795013    1.897544    1.912576    2.033702    2.482566    2.531334    2.607318    2.775881
     2.780450    2.855734    2.867071    3.050590    3.111582    3.295700    3.393576    3.418707    3.492330    3.802599
     4.181865    4.252430    4.259057    4.336650    4.505745    4.748356    4.849857    5.004338    5.033717    5.237035
     5.929787    6.228377    6.379765

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   11 emc=    -93.9733896774 demc= 4.1301E-06 wnorm= 1.3347E-03 knorm= 1.3647E-02 apxde= 1.9594E-06    *not conv.*     

               starting mcscf iteration...  12
 !timer:                                 cpu_time=     4.055 walltime=     4.134

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368298 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0808295981     -116.9437914667        0.0000000000        0.0000001250
    2       -93.9349426403     -116.7979045089        0.0000000000        0.0000001250
    3       -93.9044085951     -116.7673704637        0.0000000000        0.0000001250
    4       -93.7674157508     -116.6303776194        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  1.637459303741970E-004
 Total number of micro iterations:   11

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.99991536 pnorm= 0.0000E+00 rznorm= 8.8889E-08 rpnorm= 0.0000E+00 noldr= 11 nnewr= 11 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.713882  -23.223411   -2.321930   -1.854409   -1.668611

 qvv(*) eigenvalues. symmetry block  1
    -0.419799   -0.199316   -0.162107   -0.123317   -0.084991   -0.041582   -0.024873    0.023195    0.042000    0.066617
     0.082600    0.094201    0.141765    0.199730    0.221265    0.275997    0.322584    0.393809    0.419662    0.450744
     0.475849    0.497415    0.524708    0.560098    0.617382    0.645364    0.665579    0.730030    0.803511    0.818339
     0.907828    0.923096    0.955846    0.966787    1.004210    1.194625    1.276692    1.307244    1.324458    1.335367
     1.358104    1.587591    1.795085    1.897668    1.912573    2.033751    2.482561    2.531391    2.607355    2.775923
     2.780503    2.855778    2.867092    3.050451    3.111505    3.295304    3.393595    3.418790    3.492370    3.802315
     4.181874    4.252480    4.259105    4.336689    4.505589    4.748287    4.849780    5.004566    5.033753    5.237082
     5.929604    6.228429    6.380065

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   12 emc=    -93.9733936112 demc= 3.9338E-06 wnorm= 1.3100E-03 knorm= 1.3011E-02 apxde= 1.8819E-06    *not conv.*     

               starting mcscf iteration...  13
 !timer:                                 cpu_time=     4.372 walltime=     4.452

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368298 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0808441839     -116.9438060525        0.0000000000        0.0000001250
    2       -93.9349622935     -116.7979241621        0.0000000000        0.0000001250
    3       -93.9043857041     -116.7673475727        0.0000000000        0.0000001250
    4       -93.7677412624     -116.6307031309        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  1.614479922282428E-004
 Total number of micro iterations:   11

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.99992213 pnorm= 0.0000E+00 rznorm= 9.0788E-08 rpnorm= 0.0000E+00 noldr= 11 nnewr= 11 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.713895  -23.223358   -2.321879   -1.854359   -1.675394

 qvv(*) eigenvalues. symmetry block  1
    -0.419748   -0.199316   -0.162121   -0.123317   -0.085002   -0.041552   -0.024840    0.023217    0.042020    0.066641
     0.082601    0.094203    0.141805    0.199733    0.221223    0.276024    0.322576    0.393811    0.419628    0.450709
     0.475806    0.497434    0.524741    0.560122    0.617416    0.645404    0.665631    0.730051    0.803481    0.818263
     0.907813    0.923071    0.955873    0.966814    1.004185    1.194696    1.276736    1.307321    1.324469    1.335402
     1.358150    1.587619    1.795154    1.897786    1.912572    2.033799    2.482555    2.531446    2.607391    2.775964
     2.780555    2.855820    2.867112    3.050316    3.111427    3.294924    3.393611    3.418868    3.492409    3.802042
     4.181881    4.252529    4.259153    4.336727    4.505435    4.748217    4.849708    5.004783    5.033785    5.237129
     5.929433    6.228479    6.380350

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   13 emc=    -93.9733973939 demc= 3.7827E-06 wnorm= 1.2916E-03 knorm= 1.2479E-02 apxde= 1.8225E-06    *not conv.*     

               starting mcscf iteration...  14
 !timer:                                 cpu_time=     4.690 walltime=     4.770

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368297 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0808583420     -116.9438202106        0.0000000000        0.0000001250
    2       -93.9349814234     -116.7979432920        0.0000000000        0.0000001250
    3       -93.9043634166     -116.7673252852        0.0000000000        0.0000001250
    4       -93.7680558361     -116.6310177047        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  1.597979989893102E-004
 Total number of micro iterations:   11

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.99992763 pnorm= 0.0000E+00 rznorm= 9.2887E-08 rpnorm= 0.0000E+00 noldr= 11 nnewr= 11 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.713907  -23.223306   -2.321828   -1.854309   -1.682016

 qvv(*) eigenvalues. symmetry block  1
    -0.419698   -0.199316   -0.162136   -0.123316   -0.085012   -0.041521   -0.024809    0.023238    0.042041    0.066664
     0.082601    0.094206    0.141844    0.199738    0.221181    0.276051    0.322569    0.393813    0.419594    0.450675
     0.475764    0.497452    0.524772    0.560147    0.617450    0.645444    0.665681    0.730072    0.803452    0.818189
     0.907798    0.923047    0.955900    0.966842    1.004160    1.194766    1.276780    1.307396    1.324481    1.335435
     1.358195    1.587647    1.795221    1.897900    1.912575    2.033844    2.482548    2.531501    2.607426    2.776005
     2.780606    2.855862    2.867132    3.050185    3.111350    3.294558    3.393625    3.418940    3.492448    3.801778
     4.181887    4.252577    4.259200    4.336764    4.505284    4.748146    4.849641    5.004991    5.033813    5.237176
     5.929271    6.228529    6.380621

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   14 emc=    -93.9734010607 demc= 3.6669E-06 wnorm= 1.2784E-03 knorm= 1.2031E-02 apxde= 1.7774E-06    *not conv.*     

               starting mcscf iteration...  15
 !timer:                                 cpu_time=     5.011 walltime=     5.113

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368298 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0808721446     -116.9438340132        0.0000000000        0.0000001250
    2       -93.9350001543     -116.7979620229        0.0000000000        0.0000001250
    3       -93.9043416212     -116.7673034898        0.0000000000        0.0000001250
    4       -93.7683612182     -116.6313230868        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  1.586911929237368E-004
 Total number of micro iterations:   11

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.99993214 pnorm= 0.0000E+00 rznorm= 9.5147E-08 rpnorm= 0.0000E+00 noldr= 11 nnewr= 11 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.713919  -23.223255   -2.321778   -1.854260   -1.688504

 qvv(*) eigenvalues. symmetry block  1
    -0.419648   -0.199316   -0.162150   -0.123315   -0.085022   -0.041492   -0.024778    0.023260    0.042060    0.066686
     0.082601    0.094208    0.141881    0.199743    0.221140    0.276078    0.322562    0.393814    0.419562    0.450642
     0.475723    0.497469    0.524803    0.560171    0.617483    0.645483    0.665731    0.730092    0.803423    0.818118
     0.907784    0.923023    0.955926    0.966869    1.004135    1.194833    1.276823    1.307469    1.324492    1.335467
     1.358240    1.587675    1.795286    1.898011    1.912579    2.033888    2.482541    2.531555    2.607461    2.776044
     2.780656    2.855904    2.867152    3.050057    3.111273    3.294202    3.393637    3.419007    3.492486    3.801520
     4.181893    4.252624    4.259247    4.336802    4.505135    4.748075    4.849576    5.005191    5.033838    5.237222
     5.929116    6.228579    6.380880

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   15 emc=    -93.9734046400 demc= 3.5793E-06 wnorm= 1.2695E-03 knorm= 1.1649E-02 apxde= 1.7441E-06    *not conv.*     

               starting mcscf iteration...  16
 !timer:                                 cpu_time=     5.333 walltime=     5.435

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368298 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0808856502     -116.9438475187        0.0000000000        0.0000001250
    2       -93.9350185891     -116.7979804576        0.0000000000        0.0000001250
    3       -93.9043202262     -116.7672820947        0.0000000000        0.0000001250
    4       -93.7686588343     -116.6316207029        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  1.580488151917627E-004
 Total number of micro iterations:   11

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.99993588 pnorm= 0.0000E+00 rznorm= 9.7539E-08 rpnorm= 0.0000E+00 noldr= 11 nnewr= 11 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.713931  -23.223204   -2.321728   -1.854211   -1.694881

 qvv(*) eigenvalues. symmetry block  1
    -0.419599   -0.199315   -0.162163   -0.123314   -0.085033   -0.041463   -0.024747    0.023283    0.042080    0.066709
     0.082601    0.094211    0.141917    0.199748    0.221099    0.276105    0.322555    0.393814    0.419530    0.450608
     0.475682    0.497485    0.524834    0.560195    0.617516    0.645522    0.665781    0.730113    0.803396    0.818048
     0.907770    0.922999    0.955953    0.966895    1.004110    1.194900    1.276866    1.307540    1.324503    1.335499
     1.358283    1.587703    1.795350    1.898119    1.912585    2.033931    2.482534    2.531608    2.607495    2.776083
     2.780706    2.855945    2.867171    3.049931    3.111196    3.293856    3.393647    3.419069    3.492525    3.801269
     4.181898    4.252671    4.259293    4.336838    4.504987    4.748004    4.849515    5.005386    5.033859    5.237268
     5.928969    6.228628    6.381130

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   16 emc=    -93.9734081551 demc= 3.5151E-06 wnorm= 1.2644E-03 knorm= 1.1324E-02 apxde= 1.7207E-06    *not conv.*     

               starting mcscf iteration...  17
 !timer:                                 cpu_time=     5.661 walltime=     5.763

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368298 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0808989071     -116.9438607757        0.0000000000        0.0000001250
    2       -93.9350368148     -116.7979986833        0.0000000000        0.0000001250
    3       -93.9042991548     -116.7672610233        0.0000000000        0.0000001250
    4       -93.7689498664     -116.6319117349        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  1.578105274777632E-004
 Total number of micro iterations:   11

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.99993901 pnorm= 0.0000E+00 rznorm= 1.0005E-07 rpnorm= 0.0000E+00 noldr= 11 nnewr= 11 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.713942  -23.223153   -2.321679   -1.854163   -1.701169

 qvv(*) eigenvalues. symmetry block  1
    -0.419551   -0.199315   -0.162177   -0.123312   -0.085043   -0.041434   -0.024717    0.023305    0.042099    0.066731
     0.082600    0.094215    0.141953    0.199754    0.221059    0.276132    0.322548    0.393815    0.419498    0.450575
     0.475643    0.497500    0.524865    0.560219    0.617549    0.645561    0.665830    0.730133    0.803369    0.817980
     0.907756    0.922976    0.955979    0.966922    1.004084    1.194964    1.276909    1.307611    1.324514    1.335529
     1.358327    1.587731    1.795412    1.898224    1.912594    2.033973    2.482526    2.531661    2.607529    2.776122
     2.780755    2.855985    2.867191    3.049806    3.111118    3.293518    3.393656    3.419128    3.492563    3.801022
     4.181902    4.252718    4.259338    4.336875    4.504840    4.747931    4.849456    5.005575    5.033878    5.237313
     5.928828    6.228676    6.381371

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   17 emc=    -93.9734116256 demc= 3.4704E-06 wnorm= 1.2625E-03 knorm= 1.1044E-02 apxde= 1.7057E-06    *not conv.*     

               starting mcscf iteration...  18
 !timer:                                 cpu_time=     5.989 walltime=     6.089

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368298 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0809119559     -116.9438738245        0.0000000000        0.0000001250
    2       -93.9350549062     -116.7980167748        0.0000000000        0.0000001250
    3       -93.9042783419     -116.7672402104        0.0000000000        0.0000001250
    4       -93.7692353070     -116.6321971755        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  1.579293814745236E-004
 Total number of micro iterations:   11

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.99994164 pnorm= 0.0000E+00 rznorm= 1.0265E-07 rpnorm= 0.0000E+00 noldr= 11 nnewr= 11 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.713953  -23.223102   -2.321630   -1.854115   -1.707382

 qvv(*) eigenvalues. symmetry block  1
    -0.419502   -0.199314   -0.162191   -0.123311   -0.085054   -0.041405   -0.024687    0.023328    0.042118    0.066753
     0.082599    0.094218    0.141987    0.199760    0.221019    0.276159    0.322540    0.393815    0.419467    0.450542
     0.475604    0.497514    0.524895    0.560244    0.617582    0.645599    0.665879    0.730153    0.803343    0.817913
     0.907743    0.922952    0.956006    0.966948    1.004059    1.195028    1.276952    1.307680    1.324525    1.335559
     1.358370    1.587760    1.795474    1.898327    1.912603    2.034014    2.482517    2.531713    2.607563    2.776161
     2.780805    2.856025    2.867211    3.049684    3.111040    3.293187    3.393663    3.419184    3.492601    3.800780
     4.181906    4.252765    4.259384    4.336912    4.504694    4.747857    4.849400    5.005758    5.033895    5.237359
     5.928692    6.228725    6.381605

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   18 emc=    -93.9734150680 demc= 3.4425E-06 wnorm= 1.2634E-03 knorm= 1.0804E-02 apxde= 1.6981E-06    *not conv.*     

               starting mcscf iteration...  19
 !timer:                                 cpu_time=     6.308 walltime=     6.409

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368297 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0809248307     -116.9438866993        0.0000000000        0.0000001250
    2       -93.9350729292     -116.7980347978        0.0000000000        0.0000001250
    3       -93.9042577315     -116.7672196001        0.0000000000        0.0000001250
    4       -93.7695159995     -116.6324778680        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  1.583683775723225E-004
 Total number of micro iterations:   11

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.99994386 pnorm= 0.0000E+00 rznorm= 1.0535E-07 rpnorm= 0.0000E+00 noldr= 11 nnewr= 11 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.713964  -23.223052   -2.321581   -1.854068   -1.713535

 qvv(*) eigenvalues. symmetry block  1
    -0.419454   -0.199314   -0.162204   -0.123309   -0.085065   -0.041377   -0.024658    0.023351    0.042137    0.066774
     0.082599    0.094222    0.142021    0.199767    0.220979    0.276186    0.322533    0.393814    0.419436    0.450509
     0.475565    0.497528    0.524925    0.560268    0.617614    0.645637    0.665928    0.730173    0.803317    0.817846
     0.907731    0.922928    0.956032    0.966975    1.004034    1.195091    1.276995    1.307749    1.324536    1.335588
     1.358413    1.587789    1.795534    1.898429    1.912615    2.034054    2.482508    2.531766    2.607596    2.776200
     2.780853    2.856066    2.867231    3.049563    3.110961    3.292862    3.393669    3.419236    3.492639    3.800541
     4.181909    4.252812    4.259430    4.336948    4.504548    4.747783    4.849346    5.005938    5.033909    5.237405
     5.928562    6.228773    6.381831

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   19 emc=    -93.9734184972 demc= 3.4291E-06 wnorm= 1.2669E-03 knorm= 1.0596E-02 apxde= 1.6970E-06    *not conv.*     

               starting mcscf iteration...  20
 !timer:                                 cpu_time=     6.629 walltime=     6.730

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368297 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0809375609     -116.9438994295        0.0000000000        0.0000001250
    2       -93.9350909422     -116.7980528107        0.0000000000        0.0000001250
    3       -93.9042372748     -116.7671991434        0.0000000000        0.0000001250
    4       -93.7697926672     -116.6327545358        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  1.590980492123587E-004
 Total number of micro iterations:   11

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.99994573 pnorm= 0.0000E+00 rznorm= 1.0813E-07 rpnorm= 0.0000E+00 noldr= 11 nnewr= 11 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.713974  -23.223002   -2.321532   -1.854020   -1.719641

 qvv(*) eigenvalues. symmetry block  1
    -0.419406   -0.199313   -0.162217   -0.123307   -0.085076   -0.041348   -0.024629    0.023374    0.042155    0.066796
     0.082597    0.094225    0.142054    0.199774    0.220939    0.276213    0.322525    0.393814    0.419406    0.450476
     0.475527    0.497541    0.524956    0.560292    0.617647    0.645676    0.665977    0.730193    0.803291    0.817781
     0.907718    0.922905    0.956059    0.967001    1.004008    1.195153    1.277038    1.307818    1.324547    1.335616
     1.358456    1.587818    1.795594    1.898528    1.912628    2.034093    2.482498    2.531818    2.607630    2.776238
     2.780902    2.856106    2.867251    3.049442    3.110881    3.292541    3.393674    3.419285    3.492677    3.800305
     4.181911    4.252859    4.259475    4.336985    4.504402    4.747707    4.849294    5.006114    5.033921    5.237450
     5.928435    6.228821    6.382052

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   20 emc=    -93.9734219260 demc= 3.4288E-06 wnorm= 1.2728E-03 knorm= 1.0418E-02 apxde= 1.7018E-06    *not conv.*     

               starting mcscf iteration...  21
 !timer:                                 cpu_time=     6.948 walltime=     7.049

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368298 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0809501717     -116.9439120402        0.0000000000        0.0000001250
    2       -93.9351089981     -116.7980708666        0.0000000000        0.0000001250
    3       -93.9042169289     -116.7671787974        0.0000000000        0.0000001250
    4       -93.7700659364     -116.6330278049        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  1.600947281645077E-004
 Total number of micro iterations:   11

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.99994732 pnorm= 0.0000E+00 rznorm= 1.1099E-07 rpnorm= 0.0000E+00 noldr= 11 nnewr= 11 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.713985  -23.222951   -2.321484   -1.853973   -1.725710

 qvv(*) eigenvalues. symmetry block  1
    -0.419358   -0.199313   -0.162230   -0.123304   -0.085087   -0.041320   -0.024599    0.023398    0.042174    0.066818
     0.082596    0.094229    0.142086    0.199781    0.220900    0.276240    0.322518    0.393813    0.419375    0.450443
     0.475488    0.497554    0.524986    0.560317    0.617680    0.645714    0.666026    0.730214    0.803265    0.817716
     0.907706    0.922881    0.956086    0.967028    1.003983    1.195215    1.277081    1.307886    1.324558    1.335644
     1.358498    1.587848    1.795653    1.898627    1.912643    2.034132    2.482488    2.531870    2.607663    2.776277
     2.780951    2.856146    2.867271    3.049323    3.110801    3.292225    3.393677    3.419332    3.492715    3.800072
     4.181913    4.252907    4.259521    4.337021    4.504256    4.747630    4.849243    5.006287    5.033930    5.237496
     5.928313    6.228870    6.382267

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   21 emc=    -93.9734253662 demc= 3.4402E-06 wnorm= 1.2808E-03 knorm= 1.0264E-02 apxde= 1.7120E-06    *not conv.*     

               starting mcscf iteration...  22
 !timer:                                 cpu_time=     7.267 walltime=     7.367

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368297 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0809626850     -116.9439245536        0.0000000000        0.0000001250
    2       -93.9351271454     -116.7980890139        0.0000000000        0.0000001250
    3       -93.9041966554     -116.7671585240        0.0000000000        0.0000001250
    4       -93.7703363528     -116.6332982214        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  1.613392735584239E-004
 Total number of micro iterations:   11

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.99994867 pnorm= 0.0000E+00 rznorm= 1.1392E-07 rpnorm= 0.0000E+00 noldr= 11 nnewr= 11 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.713995  -23.222901   -2.321435   -1.853926   -1.731752

 qvv(*) eigenvalues. symmetry block  1
    -0.419310   -0.199312   -0.162244   -0.123302   -0.085098   -0.041292   -0.024570    0.023422    0.042192    0.066839
     0.082595    0.094233    0.142118    0.199789    0.220860    0.276267    0.322510    0.393812    0.419345    0.450410
     0.475450    0.497566    0.525016    0.560341    0.617713    0.645752    0.666075    0.730234    0.803240    0.817652
     0.907694    0.922858    0.956113    0.967054    1.003957    1.195276    1.277125    1.307954    1.324569    1.335671
     1.358541    1.587878    1.795712    1.898725    1.912658    2.034170    2.482478    2.531923    2.607697    2.776316
     2.781000    2.856186    2.867292    3.049204    3.110719    3.291912    3.393680    3.419376    3.492754    3.799840
     4.181915    4.252954    4.259567    4.337058    4.504110    4.747552    4.849194    5.006457    5.033938    5.237542
     5.928195    6.228918    6.382478

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   22 emc=    -93.9734288286 demc= 3.4624E-06 wnorm= 1.2907E-03 knorm= 1.0132E-02 apxde= 1.7271E-06    *not conv.*     

               starting mcscf iteration...  23
 !timer:                                 cpu_time=     7.586 walltime=     7.685

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368298 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0809751200     -116.9439369886        0.0000000000        0.0000001250
    2       -93.9351454293     -116.7981072979        0.0000000000        0.0000001250
    3       -93.9041764203     -116.7671382888        0.0000000000        0.0000001250
    4       -93.7706043958     -116.6335662643        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  1.628161234408464E-004
 Total number of micro iterations:   11

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.99994980 pnorm= 0.0000E+00 rznorm= 1.1692E-07 rpnorm= 0.0000E+00 noldr= 11 nnewr= 11 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.714006  -23.222850   -2.321385   -1.853879   -1.737775

 qvv(*) eigenvalues. symmetry block  1
    -0.419262   -0.199312   -0.162257   -0.123299   -0.085109   -0.041263   -0.024541    0.023446    0.042211    0.066861
     0.082593    0.094238    0.142150    0.199798    0.220821    0.276295    0.322502    0.393810    0.419315    0.450376
     0.475412    0.497578    0.525046    0.560366    0.617747    0.645791    0.666124    0.730255    0.803215    0.817588
     0.907683    0.922834    0.956140    0.967081    1.003930    1.195337    1.277169    1.308022    1.324580    1.335698
     1.358584    1.587909    1.795770    1.898821    1.912676    2.034208    2.482467    2.531975    2.607730    2.776354
     2.781049    2.856226    2.867312    3.049085    3.110637    3.291603    3.393682    3.419417    3.492793    3.799609
     4.181916    4.253002    4.259613    4.337095    4.503963    4.747472    4.849146    5.006625    5.033944    5.237588
     5.928080    6.228967    6.382684

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   23 emc=    -93.9734323232 demc= 3.4946E-06 wnorm= 1.3025E-03 knorm= 1.0019E-02 apxde= 1.7470E-06    *not conv.*     

               starting mcscf iteration...  24
 !timer:                                 cpu_time=     7.905 walltime=     8.004

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368298 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0809874936     -116.9439493621        0.0000000000        0.0000001250
    2       -93.9351638924     -116.7981257609        0.0000000000        0.0000001250
    3       -93.9041561924     -116.7671180609        0.0000000000        0.0000001250
    4       -93.7708704882     -116.6338323568        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  1.645125753976951E-004
 Total number of micro iterations:   11

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.99995076 pnorm= 0.0000E+00 rznorm= 1.1997E-07 rpnorm= 0.0000E+00 noldr= 11 nnewr= 11 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.714016  -23.222799   -2.321336   -1.853832   -1.743787

 qvv(*) eigenvalues. symmetry block  1
    -0.419214   -0.199311   -0.162270   -0.123296   -0.085121   -0.041235   -0.024512    0.023470    0.042229    0.066882
     0.082591    0.094242    0.142181    0.199807    0.220781    0.276323    0.322494    0.393808    0.419285    0.450343
     0.475374    0.497589    0.525077    0.560391    0.617780    0.645829    0.666174    0.730275    0.803190    0.817524
     0.907671    0.922810    0.956168    0.967108    1.003904    1.195397    1.277213    1.308089    1.324591    1.335724
     1.358627    1.587940    1.795828    1.898917    1.912694    2.034245    2.482456    2.532028    2.607764    2.776393
     2.781099    2.856266    2.867333    3.048967    3.110553    3.291296    3.393682    3.419457    3.492832    3.799380
     4.181916    4.253050    4.259659    4.337132    4.503815    4.747391    4.849100    5.006790    5.033948    5.237635
     5.927968    6.229016    6.382886

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   24 emc=    -93.9734358594 demc= 3.5362E-06 wnorm= 1.3161E-03 knorm= 9.9239E-03 apxde= 1.7713E-06    *not conv.*     

               starting mcscf iteration...  25
 !timer:                                 cpu_time=     8.224 walltime=     8.324

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368298 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0809998204     -116.9439616890        0.0000000000        0.0000001250
    2       -93.9351825749     -116.7981444435        0.0000000000        0.0000001250
    3       -93.9041359435     -116.7670978120        0.0000000000        0.0000001250
    4       -93.7711350053     -116.6340968738        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  1.664182319150224E-004
 Total number of micro iterations:   11

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.99995155 pnorm= 0.0000E+00 rznorm= 1.2308E-07 rpnorm= 0.0000E+00 noldr= 11 nnewr= 11 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.714027  -23.222747   -2.321286   -1.853785   -1.749793

 qvv(*) eigenvalues. symmetry block  1
    -0.419165   -0.199310   -0.162283   -0.123293   -0.085133   -0.041206   -0.024483    0.023496    0.042248    0.066904
     0.082589    0.094247    0.142212    0.199816    0.220741    0.276351    0.322486    0.393806    0.419255    0.450309
     0.475336    0.497600    0.525108    0.560416    0.617814    0.645868    0.666224    0.730296    0.803165    0.817461
     0.907660    0.922786    0.956195    0.967135    1.003877    1.195457    1.277258    1.308157    1.324602    1.335750
     1.358670    1.587971    1.795886    1.899013    1.912714    2.034282    2.482444    2.532081    2.607797    2.776432
     2.781149    2.856307    2.867354    3.048849    3.110468    3.290990    3.393682    3.419494    3.492871    3.799152
     4.181916    4.253099    4.259706    4.337170    4.503666    4.747308    4.849055    5.006953    5.033950    5.237682
     5.927859    6.229065    6.383084

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   25 emc=    -93.9734394463 demc= 3.5868E-06 wnorm= 1.3313E-03 knorm= 9.8436E-03 apxde= 1.8000E-06    *not conv.*     

               starting mcscf iteration...  26
 !timer:                                 cpu_time=     8.545 walltime=     8.644

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368298 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0810121135     -116.9439739821        0.0000000000        0.0000001250
    2       -93.9352015159     -116.7981633844        0.0000000000        0.0000001250
    3       -93.9041156477     -116.7670775163        0.0000000000        0.0000001250
    4       -93.7713982810     -116.6343601495        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  1.685245666296660E-004
 Total number of micro iterations:   11

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.99995221 pnorm= 0.0000E+00 rznorm= 1.2623E-07 rpnorm= 0.0000E+00 noldr= 11 nnewr= 11 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.714037  -23.222696   -2.321236   -1.853738   -1.755801

 qvv(*) eigenvalues. symmetry block  1
    -0.419116   -0.199310   -0.162296   -0.123290   -0.085144   -0.041178   -0.024454    0.023521    0.042266    0.066926
     0.082587    0.094252    0.142243    0.199825    0.220701    0.276379    0.322478    0.393804    0.419225    0.450275
     0.475298    0.497610    0.525139    0.560442    0.617848    0.645907    0.666275    0.730317    0.803140    0.817398
     0.907649    0.922762    0.956223    0.967162    1.003850    1.195518    1.277303    1.308225    1.324613    1.335776
     1.358713    1.588003    1.795944    1.899108    1.912736    2.034318    2.482431    2.532134    2.607831    2.776472
     2.781199    2.856347    2.867375    3.048730    3.110381    3.290687    3.393680    3.419529    3.492911    3.798924
     4.181915    4.253148    4.259753    4.337208    4.503516    4.747224    4.849011    5.007115    5.033950    5.237730
     5.927752    6.229115    6.383279

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   26 emc=    -93.9734430924 demc= 3.6461E-06 wnorm= 1.3482E-03 knorm= 9.7769E-03 apxde= 1.8328E-06    *not conv.*     

               starting mcscf iteration...  27
 !timer:                                 cpu_time=     8.864 walltime=     8.964

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368297 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0810243845     -116.9439862531        0.0000000000        0.0000001250
    2       -93.9352207529     -116.7981826214        0.0000000000        0.0000001250
    3       -93.9040952811     -116.7670571497        0.0000000000        0.0000001250
    4       -93.7716606135     -116.6346224820        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  1.708245801053747E-004
 Total number of micro iterations:   11

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.99995274 pnorm= 0.0000E+00 rznorm= 1.2941E-07 rpnorm= 0.0000E+00 noldr= 11 nnewr= 11 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.714047  -23.222643   -2.321185   -1.853691   -1.761815

 qvv(*) eigenvalues. symmetry block  1
    -0.419067   -0.199309   -0.162309   -0.123286   -0.085157   -0.041149   -0.024425    0.023547    0.042285    0.066948
     0.082584    0.094257    0.142274    0.199835    0.220660    0.276407    0.322469    0.393802    0.419195    0.450241
     0.475260    0.497621    0.525170    0.560468    0.617882    0.645947    0.666326    0.730339    0.803115    0.817335
     0.907638    0.922737    0.956252    0.967190    1.003822    1.195578    1.277349    1.308293    1.324624    1.335802
     1.358757    1.588036    1.796001    1.899202    1.912758    2.034355    2.482418    2.532188    2.607866    2.776511
     2.781249    2.856388    2.867397    3.048612    3.110293    3.290385    3.393678    3.419562    3.492951    3.798697
     4.181914    4.253197    4.259800    4.337246    4.503364    4.747138    4.848968    5.007275    5.033948    5.237778
     5.927649    6.229165    6.383471

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   27 emc=    -93.9734468062 demc= 3.7138E-06 wnorm= 1.3666E-03 knorm= 9.7224E-03 apxde= 1.8697E-06    *not conv.*     

               starting mcscf iteration...  28
 !timer:                                 cpu_time=     9.184 walltime=     9.283

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368293 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0810366434     -116.9439985120        0.0000000000        0.0000001250
    2       -93.9352403228     -116.7982021914        0.0000000000        0.0000001250
    3       -93.9040748215     -116.7670366901        0.0000000000        0.0000001250
    4       -93.7719222697     -116.6348841382        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  1.733125226373544E-004
 Total number of micro iterations:   11

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.99995316 pnorm= 0.0000E+00 rznorm= 1.3261E-07 rpnorm= 0.0000E+00 noldr= 11 nnewr= 11 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.714058  -23.222591   -2.321134   -1.853645   -1.767840

 qvv(*) eigenvalues. symmetry block  1
    -0.419017   -0.199308   -0.162322   -0.123283   -0.085169   -0.041119   -0.024396    0.023574    0.042304    0.066970
     0.082582    0.094263    0.142304    0.199846    0.220619    0.276436    0.322461    0.393799    0.419165    0.450207
     0.475222    0.497630    0.525201    0.560494    0.617917    0.645986    0.666377    0.730360    0.803090    0.817272
     0.907627    0.922712    0.956281    0.967218    1.003794    1.195638    1.277396    1.308362    1.324636    1.335827
     1.358801    1.588070    1.796059    1.899297    1.912782    2.034391    2.482405    2.532243    2.607900    2.776552
     2.781300    2.856430    2.867418    3.048493    3.110204    3.290084    3.393674    3.419592    3.492992    3.798469
     4.181913    4.253247    4.259848    4.337284    4.503211    4.747050    4.848926    5.007433    5.033945    5.237826
     5.927548    6.229216    6.383660

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   28 emc=    -93.9734505959 demc= 3.7898E-06 wnorm= 1.3865E-03 knorm= 9.6790E-03 apxde= 1.9107E-06    *not conv.*     

               starting mcscf iteration...  29
 !timer:                                 cpu_time=     9.503 walltime=     9.601

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368298 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0810488992     -116.9440107678        0.0000000000        0.0000001250
    2       -93.9352602620     -116.7982221306        0.0000000000        0.0000001250
    3       -93.9040542482     -116.7670161168        0.0000000000        0.0000001250
    4       -93.7721834886     -116.6351453572        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  1.759836671355675E-004
 Total number of micro iterations:   11

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.99995348 pnorm= 0.0000E+00 rznorm= 1.3583E-07 rpnorm= 0.0000E+00 noldr= 11 nnewr= 11 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.714068  -23.222537   -2.321082   -1.853599   -1.773880

 qvv(*) eigenvalues. symmetry block  1
    -0.418967   -0.199308   -0.162336   -0.123279   -0.085182   -0.041090   -0.024367    0.023601    0.042322    0.066992
     0.082579    0.094268    0.142334    0.199857    0.220578    0.276466    0.322452    0.393796    0.419135    0.450172
     0.475183    0.497640    0.525233    0.560520    0.617952    0.646027    0.666429    0.730382    0.803065    0.817209
     0.907616    0.922687    0.956310    0.967246    1.003765    1.195698    1.277443    1.308430    1.324647    1.335852
     1.358845    1.588104    1.796117    1.899392    1.912808    2.034427    2.482391    2.532298    2.607935    2.776592
     2.781351    2.856472    2.867441    3.048374    3.110113    3.289783    3.393669    3.419621    3.493033    3.798241
     4.181911    4.253298    4.259897    4.337323    4.503057    4.746961    4.848886    5.007591    5.033940    5.237875
     5.927449    6.229267    6.383846

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   29 emc=    -93.9734544698 demc= 3.8739E-06 wnorm= 1.4079E-03 knorm= 9.6457E-03 apxde= 1.9557E-06    *not conv.*     

               starting mcscf iteration...  30
 !timer:                                 cpu_time=     9.821 walltime=     9.920

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368298 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0810611596     -116.9440230281        0.0000000000        0.0000001250
    2       -93.9352806064     -116.7982424750        0.0000000000        0.0000001250
    3       -93.9040335418     -116.7669954103        0.0000000000        0.0000001250
    4       -93.7724444844     -116.6354063529        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  1.788341213781179E-004
 Total number of micro iterations:   11

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.99995371 pnorm= 0.0000E+00 rznorm= 1.3902E-07 rpnorm= 0.0000E+00 noldr= 11 nnewr= 11 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.714079  -23.222483   -2.321030   -1.853553   -1.779939

 qvv(*) eigenvalues. symmetry block  1
    -0.418916   -0.199307   -0.162349   -0.123274   -0.085195   -0.041061   -0.024337    0.023628    0.042341    0.067015
     0.082576    0.094274    0.142364    0.199868    0.220537    0.276495    0.322443    0.393792    0.419105    0.450137
     0.475144    0.497649    0.525265    0.560547    0.617987    0.646067    0.666482    0.730404    0.803040    0.817146
     0.907605    0.922662    0.956340    0.967275    1.003736    1.195758    1.277490    1.308499    1.324658    1.335877
     1.358890    1.588138    1.796175    1.899486    1.912834    2.034462    2.482377    2.532353    2.607970    2.776633
     2.781403    2.856514    2.867463    3.048254    3.110020    3.289483    3.393664    3.419648    3.493075    3.798013
     4.181908    4.253349    4.259946    4.337363    4.502900    4.746869    4.848846    5.007747    5.033933    5.237925
     5.927352    6.229319    6.384030

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   30 emc=    -93.9734584359 demc= 3.9661E-06 wnorm= 1.4307E-03 knorm= 9.6216E-03 apxde= 2.0047E-06    *not conv.*     

               starting mcscf iteration...  31
 !timer:                                 cpu_time=    10.142 walltime=    10.241

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368297 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0810734312     -116.9440352998        0.0000000000        0.0000001250
    2       -93.9353013919     -116.7982632605        0.0000000000        0.0000001250
    3       -93.9040126840     -116.7669745525        0.0000000000        0.0000001250
    4       -93.7727054483     -116.6356673169        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  1.818606686787117E-004
 Total number of micro iterations:   11

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.99995386 pnorm= 0.0000E+00 rznorm= 1.4219E-07 rpnorm= 0.0000E+00 noldr= 11 nnewr= 11 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.714089  -23.222429   -2.320977   -1.853508   -1.786020

 qvv(*) eigenvalues. symmetry block  1
    -0.418865   -0.199306   -0.162362   -0.123270   -0.085208   -0.041031   -0.024307    0.023656    0.042360    0.067037
     0.082573    0.094280    0.142394    0.199880    0.220495    0.276526    0.322434    0.393788    0.419075    0.450101
     0.475105    0.497658    0.525298    0.560575    0.618023    0.646108    0.666535    0.730427    0.803015    0.817083
     0.907594    0.922636    0.956370    0.967304    1.003707    1.195818    1.277539    1.308569    1.324670    1.335901
     1.358935    1.588174    1.796233    1.899581    1.912862    2.034498    2.482362    2.532409    2.608005    2.776674
     2.781456    2.856556    2.867486    3.048134    3.109925    3.289184    3.393657    3.419673    3.493117    3.797785
     4.181905    4.253401    4.259995    4.337403    4.502742    4.746775    4.848807    5.007902    5.033923    5.237976
     5.927258    6.229371    6.384211

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   31 emc=    -93.9734625024 demc= 4.0665E-06 wnorm= 1.4549E-03 knorm= 9.6057E-03 apxde= 2.0577E-06    *not conv.*     

               starting mcscf iteration...  32
 !timer:                                 cpu_time=    10.461 walltime=    10.561

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368297 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0810857199     -116.9440475885        0.0000000000        0.0000001250
    2       -93.9353226544     -116.7982845229        0.0000000000        0.0000001250
    3       -93.9039916577     -116.7669535262        0.0000000000        0.0000001250
    4       -93.7729665513     -116.6359284198        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  1.850606316265595E-004
 Total number of micro iterations:   11

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.99995394 pnorm= 0.0000E+00 rznorm= 1.4529E-07 rpnorm= 0.0000E+00 noldr= 11 nnewr= 11 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.714100  -23.222373   -2.320924   -1.853465   -1.792126

 qvv(*) eigenvalues. symmetry block  1
    -0.418813   -0.199305   -0.162376   -0.123265   -0.085221   -0.041000   -0.024277    0.023685    0.042379    0.067060
     0.082569    0.094287    0.142423    0.199893    0.220453    0.276556    0.322425    0.393784    0.419044    0.450065
     0.475066    0.497667    0.525331    0.560603    0.618060    0.646150    0.666589    0.730449    0.802990    0.817019
     0.907584    0.922610    0.956400    0.967333    1.003677    1.195879    1.277588    1.308639    1.324681    1.335926
     1.358981    1.588210    1.796291    1.899676    1.912891    2.034533    2.482347    2.532466    2.608041    2.776716
     2.781509    2.856599    2.867509    3.048013    3.109828    3.288884    3.393649    3.419697    3.493160    3.797555
     4.181902    4.253454    4.260046    4.337443    4.502582    4.746680    4.848768    5.008057    5.033913    5.238027
     5.927166    6.229424    6.384390

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   32 emc=    -93.9734666773 demc= 4.1749E-06 wnorm= 1.4805E-03 knorm= 9.5974E-03 apxde= 2.1148E-06    *not conv.*     

               starting mcscf iteration...  33
 !timer:                                 cpu_time=    10.781 walltime=    10.881

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368298 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0810980304     -116.9440598990        0.0000000000        0.0000001250
    2       -93.9353444299     -116.7983062984        0.0000000000        0.0000001250
    3       -93.9039704466     -116.7669323151        0.0000000000        0.0000001250
    4       -93.7732279447     -116.6361898133        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  1.884317523551420E-004
 Total number of micro iterations:   11

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.99995396 pnorm= 0.0000E+00 rznorm= 1.4829E-07 rpnorm= 0.0000E+00 noldr= 11 nnewr= 11 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.714110  -23.222317   -2.320869   -1.853423   -1.798259

 qvv(*) eigenvalues. symmetry block  1
    -0.418760   -0.199305   -0.162390   -0.123261   -0.085235   -0.040970   -0.024247    0.023714    0.042398    0.067083
     0.082566    0.094293    0.142453    0.199905    0.220410    0.276587    0.322415    0.393780    0.419013    0.450028
     0.475026    0.497676    0.525364    0.560631    0.618097    0.646192    0.666643    0.730472    0.802964    0.816956
     0.907573    0.922583    0.956432    0.967363    1.003646    1.195939    1.277638    1.308710    1.324693    1.335950
     1.359027    1.588247    1.796350    1.899771    1.912921    2.034568    2.482331    2.532523    2.608077    2.776758
     2.781563    2.856643    2.867533    3.047891    3.109729    3.288584    3.393640    3.419718    3.493203    3.797325
     4.181898    4.253507    4.260097    4.337484    4.502420    4.746582    4.848731    5.008210    5.033900    5.238079
     5.927077    6.229478    6.384567

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   33 emc=    -93.9734709690 demc= 4.2916E-06 wnorm= 1.5075E-03 knorm= 9.5961E-03 apxde= 2.1759E-06    *not conv.*     

               starting mcscf iteration...  34
 !timer:                                 cpu_time=    11.101 walltime=    11.200

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368297 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0811103668     -116.9440722353        0.0000000000        0.0000001250
    2       -93.9353667547     -116.7983286233        0.0000000000        0.0000001250
    3       -93.9039490354     -116.7669109039        0.0000000000        0.0000001250
    4       -93.7734897623     -116.6364516308        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  1.919720861424024E-004
 Total number of micro iterations:   11

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.99995391 pnorm= 0.0000E+00 rznorm= 1.5117E-07 rpnorm= 0.0000E+00 noldr= 11 nnewr= 11 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.714121  -23.222260   -2.320814   -1.853385   -1.804420

 qvv(*) eigenvalues. symmetry block  1
    -0.418707   -0.199304   -0.162403   -0.123256   -0.085249   -0.040939   -0.024217    0.023744    0.042417    0.067107
     0.082562    0.094300    0.142482    0.199919    0.220367    0.276619    0.322406    0.393776    0.418983    0.449991
     0.474986    0.497684    0.525398    0.560660    0.618135    0.646235    0.666699    0.730496    0.802939    0.816892
     0.907563    0.922556    0.956463    0.967393    1.003615    1.196000    1.277689    1.308781    1.324704    1.335974
     1.359074    1.588284    1.796408    1.899867    1.912953    2.034604    2.482314    2.532582    2.608114    2.776801
     2.781618    2.856687    2.867557    3.047768    3.109628    3.288284    3.393630    3.419737    3.493247    3.797094
     4.181893    4.253562    4.260148    4.337526    4.502255    4.746481    4.848695    5.008363    5.033885    5.238132
     5.926989    6.229533    6.384741

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   34 emc=    -93.9734753856 demc= 4.4166E-06 wnorm= 1.5358E-03 knorm= 9.6011E-03 apxde= 2.2413E-06    *not conv.*     

               starting mcscf iteration...  35
 !timer:                                 cpu_time=    11.420 walltime=    11.519

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368297 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0811227319     -116.9440846005        0.0000000000        0.0000001250
    2       -93.9353896656     -116.7983515341        0.0000000000        0.0000001250
    3       -93.9039274095     -116.7668892780        0.0000000000        0.0000001250
    4       -93.7737521205     -116.6367139890        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  1.956799041479242E-004
 Total number of micro iterations:   11

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.99995381 pnorm= 0.0000E+00 rznorm= 1.5390E-07 rpnorm= 0.0000E+00 noldr= 11 nnewr= 11 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.714131  -23.222203   -2.320759   -1.853351   -1.810610

 qvv(*) eigenvalues. symmetry block  1
    -0.418653   -0.199303   -0.162417   -0.123250   -0.085263   -0.040908   -0.024186    0.023774    0.042437    0.067131
     0.082558    0.094307    0.142512    0.199933    0.220323    0.276651    0.322396    0.393771    0.418951    0.449954
     0.474945    0.497692    0.525432    0.560689    0.618173    0.646278    0.666755    0.730520    0.802913    0.816828
     0.907552    0.922529    0.956496    0.967423    1.003583    1.196061    1.277741    1.308854    1.324716    1.335998
     1.359121    1.588323    1.796467    1.899963    1.912986    2.034639    2.482297    2.532641    2.608151    2.776845
     2.781673    2.856731    2.867582    3.047645    3.109525    3.287983    3.393619    3.419755    3.493292    3.796861
     4.181888    4.253617    4.260201    4.337568    4.502088    4.746379    4.848659    5.008515    5.033868    5.238186
     5.926904    6.229588    6.384913

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   35 emc=    -93.9734799356 demc= 4.5500E-06 wnorm= 1.5654E-03 knorm= 9.6118E-03 apxde= 2.3108E-06    *not conv.*     

               starting mcscf iteration...  36
 !timer:                                 cpu_time=    11.739 walltime=    11.838

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368297 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0811351281     -116.9440969967        0.0000000000        0.0000001250
    2       -93.9354131996     -116.7983750682        0.0000000000        0.0000001250
    3       -93.9039055551     -116.7668674237        0.0000000000        0.0000001250
    4       -93.7740151197     -116.6369769882        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  1.995536034816842E-004
 Total number of micro iterations:   11

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.99995365 pnorm= 0.0000E+00 rznorm= 1.5647E-07 rpnorm= 0.0000E+00 noldr= 11 nnewr= 11 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.714142  -23.222144   -2.320702   -1.853324   -1.816828

 qvv(*) eigenvalues. symmetry block  1
    -0.418598   -0.199302   -0.162431   -0.123245   -0.085278   -0.040876   -0.024155    0.023806    0.042457    0.067154
     0.082554    0.094315    0.142541    0.199947    0.220279    0.276684    0.322385    0.393765    0.418920    0.449916
     0.474904    0.497700    0.525466    0.560719    0.618212    0.646322    0.666812    0.730544    0.802887    0.816763
     0.907542    0.922501    0.956528    0.967454    1.003551    1.196123    1.277793    1.308926    1.324728    1.336021
     1.359169    1.588362    1.796527    1.900059    1.913020    2.034673    2.482279    2.532700    2.608189    2.776889
     2.781729    2.856777    2.867607    3.047520    3.109420    3.287682    3.393606    3.419770    3.493338    3.796628
     4.181883    4.253673    4.260254    4.337611    4.501919    4.746274    4.848624    5.008667    5.033849    5.238240
     5.926820    6.229644    6.385083

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   36 emc=    -93.9734846276 demc= 4.6920E-06 wnorm= 1.5964E-03 knorm= 9.6279E-03 apxde= 2.3846E-06    *not conv.*     

               starting mcscf iteration...  37
 !timer:                                 cpu_time=    12.058 walltime=    12.157

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368297 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0811475567     -116.9441094253        0.0000000000        0.0000001250
    2       -93.9354373945     -116.7983992631        0.0000000000        0.0000001250
    3       -93.9038834592     -116.7668453278        0.0000000000        0.0000001250
    4       -93.7742788443     -116.6372407129        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  2.035916217236785E-004
 Total number of micro iterations:   11

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.99995345 pnorm= 0.0000E+00 rznorm= 1.5896E-07 rpnorm= 0.0000E+00 noldr= 11 nnewr= 11 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.714153  -23.222085   -2.320645   -1.853311   -1.823069

 qvv(*) eigenvalues. symmetry block  1
    -0.418543   -0.199302   -0.162445   -0.123239   -0.085293   -0.040844   -0.024124    0.023837    0.042476    0.067179
     0.082550    0.094322    0.142570    0.199962    0.220234    0.276717    0.322375    0.393760    0.418889    0.449877
     0.474862    0.497707    0.525502    0.560749    0.618251    0.646367    0.666870    0.730568    0.802861    0.816699
     0.907532    0.922473    0.956562    0.967486    1.003518    1.196185    1.277847    1.309000    1.324740    1.336045
     1.359217    1.588402    1.796587    1.900157    1.913056    2.034708    2.482260    2.532761    2.608227    2.776933
     2.781786    2.856822    2.867632    3.047395    3.109312    3.287380    3.393592    3.419784    3.493384    3.796393
     4.181876    4.253730    4.260308    4.337655    4.501747    4.746166    4.848590    5.008818    5.033829    5.238296
     5.926739    6.229701    6.385251

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   37 emc=    -93.9734894702 demc= 4.8426E-06 wnorm= 1.6287E-03 knorm= 9.6487E-03 apxde= 2.4628E-06    *not conv.*     

               starting mcscf iteration...  38
 !timer:                                 cpu_time=    12.379 walltime=    12.478

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368297 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0811600184     -116.9441218870        0.0000000000        0.0000001250
    2       -93.9354622883     -116.7984241569        0.0000000000        0.0000001250
    3       -93.9038611096     -116.7668229781        0.0000000000        0.0000001250
    4       -93.7745433637     -116.6375052322        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  2.077923550701873E-004
 Total number of micro iterations:   11

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.99995321 pnorm= 0.0000E+00 rznorm= 1.6165E-07 rpnorm= 0.0000E+00 noldr= 11 nnewr= 11 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.714164  -23.222025   -2.320586   -1.853320   -1.829325

 qvv(*) eigenvalues. symmetry block  1
    -0.418487   -0.199301   -0.162459   -0.123233   -0.085308   -0.040812   -0.024092    0.023870    0.042496    0.067203
     0.082545    0.094330    0.142599    0.199977    0.220189    0.276751    0.322364    0.393754    0.418857    0.449838
     0.474820    0.497715    0.525537    0.560780    0.618291    0.646412    0.666929    0.730593    0.802835    0.816633
     0.907521    0.922444    0.956596    0.967518    1.003484    1.196247    1.277901    1.309074    1.324751    1.336068
     1.359266    1.588444    1.796647    1.900255    1.913093    2.034743    2.482241    2.532823    2.608266    2.776979
     2.781844    2.856869    2.867658    3.047269    3.109202    3.287078    3.393577    3.419796    3.493431    3.796156
     4.181870    4.253788    4.260363    4.337699    4.501572    4.746056    4.848557    5.008968    5.033806    5.238352
     5.926661    6.229759    6.385416

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   38 emc=    -93.9734944721 demc= 5.0019E-06 wnorm= 1.6623E-03 knorm= 9.6739E-03 apxde= 2.5454E-06    *not conv.*     

               starting mcscf iteration...  39
 !timer:                                 cpu_time=    12.698 walltime=    12.796

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368298 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0811725128     -116.9441343814        0.0000000000        0.0000001250
    2       -93.9354879198     -116.7984497884        0.0000000000        0.0000001250
    3       -93.9038384945     -116.7668003630        0.0000000000        0.0000001250
    4       -93.7748087317     -116.6377706002        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  2.121540774662896E-004
 Total number of micro iterations:   11

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.99995292 pnorm= 0.0000E+00 rznorm= 1.6565E-07 rpnorm= 0.0000E+00 noldr= 11 nnewr= 11 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.714174  -23.221963   -2.320527   -1.853380   -1.835569

 qvv(*) eigenvalues. symmetry block  1
    -0.418430   -0.199300   -0.162474   -0.123227   -0.085324   -0.040779   -0.024060    0.023903    0.042517    0.067228
     0.082540    0.094338    0.142629    0.199993    0.220143    0.276786    0.322353    0.393747    0.418825    0.449798
     0.474778    0.497722    0.525573    0.560812    0.618332    0.646458    0.666989    0.730619    0.802808    0.816568
     0.907511    0.922415    0.956631    0.967551    1.003450    1.196309    1.277957    1.309150    1.324763    1.336091
     1.359316    1.588486    1.796708    1.900353    1.913131    2.034778    2.482222    2.532885    2.608305    2.777025
     2.781903    2.856916    2.867685    3.047141    3.109089    3.286775    3.393561    3.419805    3.493479    3.795918
     4.181862    4.253847    4.260419    4.337744    4.501394    4.745942    4.848525    5.009118    5.033780    5.238410
     5.926584    6.229818    6.385580

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   39 emc=    -93.9734996424 demc= 5.1703E-06 wnorm= 1.6972E-03 knorm= 9.7031E-03 apxde= 2.6324E-06    *not conv.*     

               starting mcscf iteration...  40
 !timer:                                 cpu_time=    13.017 walltime=    13.115

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368298 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0811850389     -116.9441469075        0.0000000000        0.0000001250
    2       -93.9355143281     -116.7984761967        0.0000000000        0.0000001250
    3       -93.9038156031     -116.7667774717        0.0000000000        0.0000001250
    4       -93.7750749873     -116.6380368558        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  2.166748614872068E-004
 Total number of micro iterations:   11

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.99995261 pnorm= 0.0000E+00 rznorm= 1.7494E-07 rpnorm= 0.0000E+00 noldr= 11 nnewr= 11 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.714185  -23.221901   -2.320467   -1.853571   -1.841721

 qvv(*) eigenvalues. symmetry block  1
    -0.418372   -0.199299   -0.162488   -0.123220   -0.085340   -0.040745   -0.024028    0.023937    0.042537    0.067253
     0.082535    0.094347    0.142658    0.200010    0.220097    0.276821    0.322342    0.393740    0.418792    0.449758
     0.474734    0.497729    0.525610    0.560844    0.618373    0.646505    0.667050    0.730645    0.802782    0.816502
     0.907501    0.922386    0.956666    0.967584    1.003415    1.196372    1.278013    1.309226    1.324775    1.336114
     1.359366    1.588529    1.796769    1.900453    1.913171    2.034812    2.482201    2.532948    2.608345    2.777071
     2.781962    2.856964    2.867712    3.047013    3.108973    3.286471    3.393543    3.419813    3.493528    3.795679
     4.181854    4.253907    4.260476    4.337790    4.501214    4.745826    4.848493    5.009267    5.033753    5.238468
     5.926510    6.229878    6.385741

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   40 emc=    -93.9735049901 demc= 5.3477E-06 wnorm= 1.7334E-03 knorm= 9.7357E-03 apxde= 2.7241E-06    *not conv.*     

               starting mcscf iteration...  41
 !timer:                                 cpu_time=    13.336 walltime=    13.434

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368298 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0811975947     -116.9441594633        0.0000000000        0.0000001250
    2       -93.9355415529     -116.7985034214        0.0000000000        0.0000001250
    3       -93.9037924255     -116.7667542940        0.0000000000        0.0000001250
    4       -93.7753421545     -116.6383040231        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  2.213524991793521E-004
 Total number of micro iterations:   11

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.99995226 pnorm= 0.0000E+00 rznorm= 1.9998E-07 rpnorm= 0.0000E+00 noldr= 11 nnewr= 11 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.714197  -23.221838   -2.320406   -1.854263   -1.847412

 qvv(*) eigenvalues. symmetry block  1
    -0.418313   -0.199298   -0.162503   -0.123213   -0.085356   -0.040711   -0.023995    0.023972    0.042558    0.067279
     0.082530    0.094356    0.142687    0.200027    0.220050    0.276856    0.322331    0.393733    0.418760    0.449717
     0.474691    0.497735    0.525647    0.560877    0.618415    0.646552    0.667111    0.730671    0.802755    0.816436
     0.907491    0.922355    0.956702    0.967618    1.003379    1.196436    1.278071    1.309303    1.324788    1.336137
     1.359417    1.588573    1.796830    1.900553    1.913212    2.034847    2.482180    2.533013    2.608386    2.777119
     2.782023    2.857012    2.867739    3.046883    3.108855    3.286166    3.393523    3.419819    3.493578    3.795438
     4.181845    4.253968    4.260534    4.337836    4.501031    4.745708    4.848463    5.009416    5.033723    5.238527
     5.926438    6.229939    6.385900

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   41 emc=    -93.9735105244 demc= 5.5343E-06 wnorm= 1.7708E-03 knorm= 9.7715E-03 apxde= 2.8203E-06    *not conv.*     

               starting mcscf iteration...  42
 !timer:                                 cpu_time=    13.656 walltime=    13.754

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368298 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0812101773     -116.9441720459        0.0000000000        0.0000001250
    2       -93.9355696337     -116.7985315023        0.0000000000        0.0000001250
    3       -93.9037689528     -116.7667308214        0.0000000000        0.0000001250
    4       -93.7756102431     -116.6385721117        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  2.261844195557251E-004
 Total number of micro iterations:   11

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.99995188 pnorm= 0.0000E+00 rznorm= 2.1238E-07 rpnorm= 0.0000E+00 noldr= 11 nnewr= 11 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.714208  -23.221774   -2.320344   -1.857147   -1.850950

 qvv(*) eigenvalues. symmetry block  1
    -0.418253   -0.199298   -0.162517   -0.123206   -0.085373   -0.040677   -0.023962    0.024008    0.042578    0.067305
     0.082524    0.094365    0.142715    0.200045    0.220002    0.276893    0.322319    0.393726    0.418727    0.449675
     0.474646    0.497741    0.525685    0.560910    0.618458    0.646600    0.667174    0.730698    0.802727    0.816369
     0.907481    0.922325    0.956739    0.967652    1.003343    1.196499    1.278129    1.309381    1.324800    1.336160
     1.359469    1.588618    1.796892    1.900654    1.913254    2.034881    2.482158    2.533078    2.608427    2.777167
     2.782084    2.857061    2.867768    3.046752    3.108733    3.285860    3.393502    3.419822    3.493628    3.795196
     4.181836    4.254031    4.260593    4.337883    4.500844    4.745586    4.848433    5.009563    5.033691    5.238588
     5.926369    6.230000    6.386056

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   42 emc=    -93.9735162546 demc= 5.7302E-06 wnorm= 1.8095E-03 knorm= 9.8100E-03 apxde= 2.9213E-06    *not conv.*     

               starting mcscf iteration...  43
 !timer:                                 cpu_time=    13.975 walltime=    14.073

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368298 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0812227849     -116.9441846534        0.0000000000        0.0000001250
    2       -93.9355986128     -116.7985604814        0.0000000000        0.0000001250
    3       -93.9037451730     -116.7667070415        0.0000000000        0.0000001250
    4       -93.7758792410     -116.6388411095        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  2.311675434642656E-004
 Total number of micro iterations:   11

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.99995148 pnorm= 0.0000E+00 rznorm= 2.3868E-07 rpnorm= 0.0000E+00 noldr= 11 nnewr= 11 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.714219  -23.221708   -2.320280   -1.862675   -1.851883

 qvv(*) eigenvalues. symmetry block  1
    -0.418193   -0.199297   -0.162532   -0.123199   -0.085390   -0.040642   -0.023928    0.024044    0.042599    0.067331
     0.082518    0.094375    0.142744    0.200064    0.219953    0.276930    0.322307    0.393718    0.418693    0.449633
     0.474601    0.497748    0.525724    0.560944    0.618502    0.646649    0.667238    0.730725    0.802700    0.816302
     0.907471    0.922293    0.956776    0.967687    1.003306    1.196563    1.278189    1.309460    1.324812    1.336182
     1.359522    1.588664    1.796955    1.900756    1.913298    2.034915    2.482135    2.533145    2.608468    2.777216
     2.782147    2.857111    2.867796    3.046621    3.108609    3.285554    3.393479    3.419823    3.493680    3.794951
     4.181826    4.254094    4.260652    4.337931    4.500655    4.745461    4.848404    5.009711    5.033656    5.238650
     5.926302    6.230063    6.386210

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   43 emc=    -93.9735221902 demc= 5.9356E-06 wnorm= 1.8493E-03 knorm= 9.8509E-03 apxde= 3.0269E-06    *not conv.*     

               starting mcscf iteration...  44
 !timer:                                 cpu_time=    14.294 walltime=    14.392

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368298 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0812354081     -116.9441972766        0.0000000000        0.0000001250
    2       -93.9356285277     -116.7985903963        0.0000000000        0.0000001250
    3       -93.9037210863     -116.7666829549        0.0000000000        0.0000001250
    4       -93.7761491351     -116.6391110036        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  2.362983813500269E-004
 Total number of micro iterations:   11

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.99995105 pnorm= 0.0000E+00 rznorm= 2.2344E-07 rpnorm= 0.0000E+00 noldr= 11 nnewr= 11 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.714230  -23.221642   -2.320216   -1.868937   -1.852121

 qvv(*) eigenvalues. symmetry block  1
    -0.418131   -0.199296   -0.162547   -0.123191   -0.085407   -0.040607   -0.023895    0.024082    0.042620    0.067358
     0.082512    0.094385    0.142773    0.200083    0.219904    0.276968    0.322295    0.393709    0.418660    0.449590
     0.474555    0.497753    0.525763    0.560979    0.618546    0.646699    0.667303    0.730753    0.802672    0.816234
     0.907461    0.922262    0.956814    0.967723    1.003268    1.196628    1.278250    1.309539    1.324824    1.336204
     1.359575    1.588711    1.797018    1.900859    1.913343    2.034949    2.482112    2.533212    2.608511    2.777266
     2.782210    2.857162    2.867826    3.046488    3.108482    3.285246    3.393454    3.419823    3.493732    3.794706
     4.181815    4.254159    4.260713    4.337980    4.500463    4.745333    4.848377    5.009857    5.033619    5.238712
     5.926237    6.230127    6.386361

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   44 emc=    -93.9735283407 demc= 6.1505E-06 wnorm= 1.8904E-03 knorm= 9.8939E-03 apxde= 3.1373E-06    *not conv.*     

               starting mcscf iteration...  45
 !timer:                                 cpu_time=    14.615 walltime=    14.712

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368297 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0812480456     -116.9442099142        0.0000000000        0.0000001250
    2       -93.9356594227     -116.7986212913        0.0000000000        0.0000001250
    3       -93.9036966784     -116.7666585469        0.0000000000        0.0000001250
    4       -93.7764198805     -116.6393817491        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  2.415727104510627E-004
 Total number of micro iterations:   11

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.99995061 pnorm= 0.0000E+00 rznorm= 2.0387E-07 rpnorm= 0.0000E+00 noldr= 11 nnewr= 11 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.714242  -23.221574   -2.320151   -1.875404   -1.852191

 qvv(*) eigenvalues. symmetry block  1
    -0.418068   -0.199295   -0.162562   -0.123183   -0.085425   -0.040571   -0.023860    0.024120    0.042642    0.067385
     0.082506    0.094395    0.142802    0.200103    0.219855    0.277006    0.322282    0.393700    0.418626    0.449547
     0.474509    0.497759    0.525803    0.561014    0.618591    0.646749    0.667369    0.730781    0.802644    0.816166
     0.907451    0.922229    0.956853    0.967759    1.003229    1.196693    1.278312    1.309620    1.324837    1.336226
     1.359629    1.588760    1.797081    1.900963    1.913389    2.034983    2.482088    2.533281    2.608554    2.777316
     2.782274    2.857213    2.867856    3.046353    3.108351    3.284938    3.393427    3.419820    3.493786    3.794458
     4.181803    4.254224    4.260775    4.338030    4.500267    4.745201    4.848350    5.010003    5.033579    5.238776
     5.926176    6.230192    6.386510

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   45 emc=    -93.9735347156 demc= 6.3749E-06 wnorm= 1.9326E-03 knorm= 9.9384E-03 apxde= 3.2524E-06    *not conv.*     

               starting mcscf iteration...  46
 !timer:                                 cpu_time=    14.935 walltime=    15.081

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368298 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0812606900     -116.9442225585        0.0000000000        0.0000001250
    2       -93.9356913367     -116.7986532053        0.0000000000        0.0000001250
    3       -93.9036719465     -116.7666338150        0.0000000000        0.0000001250
    4       -93.7766914233     -116.6396532918        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  2.469856212736199E-004
 Total number of micro iterations:   11

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.99995016 pnorm= 0.0000E+00 rznorm= 2.0123E-07 rpnorm= 0.0000E+00 noldr= 11 nnewr= 11 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.714254  -23.221506   -2.320085   -1.881969   -1.852200

 qvv(*) eigenvalues. symmetry block  1
    -0.418005   -0.199295   -0.162577   -0.123175   -0.085443   -0.040535   -0.023826    0.024159    0.042663    0.067412
     0.082500    0.094406    0.142830    0.200123    0.219804    0.277046    0.322269    0.393691    0.418592    0.449503
     0.474462    0.497764    0.525843    0.561051    0.618637    0.646801    0.667436    0.730809    0.802616    0.816098
     0.907442    0.922196    0.956893    0.967796    1.003190    1.196758    1.278375    1.309702    1.324849    1.336248
     1.359683    1.588809    1.797145    1.901069    1.913436    2.035017    2.482063    2.533350    2.608597    2.777367
     2.782340    2.857265    2.867887    3.046218    3.108218    3.284630    3.393399    3.419814    3.493840    3.794209
     4.181791    4.254291    4.260838    4.338080    4.500068    4.745067    4.848324    5.010147    5.033536    5.238841
     5.926117    6.230258    6.386656

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   46 emc=    -93.9735413244 demc= 6.6088E-06 wnorm= 1.9759E-03 knorm= 9.9842E-03 apxde= 3.3722E-06    *not conv.*     

               starting mcscf iteration...  47
 !timer:                                 cpu_time=    15.254 walltime=    15.400

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368297 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0812733339     -116.9442352025        0.0000000000        0.0000001250
    2       -93.9357243102     -116.7986861788        0.0000000000        0.0000001250
    3       -93.9036468858     -116.7666087543        0.0000000000        0.0000001250
    4       -93.7769636937     -116.6399255623        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  2.525314143788155E-004
 Total number of micro iterations:   11

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.99994969 pnorm= 0.0000E+00 rznorm= 2.0480E-07 rpnorm= 0.0000E+00 noldr= 11 nnewr= 11 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.714266  -23.221436   -2.320017   -1.888597   -1.852179

 qvv(*) eigenvalues. symmetry block  1
    -0.417940   -0.199294   -0.162593   -0.123166   -0.085462   -0.040498   -0.023791    0.024199    0.042685    0.067440
     0.082493    0.094416    0.142859    0.200145    0.219753    0.277086    0.322256    0.393681    0.418558    0.449458
     0.474414    0.497769    0.525884    0.561088    0.618684    0.646853    0.667505    0.730839    0.802588    0.816029
     0.907432    0.922163    0.956934    0.967834    1.003149    1.196823    1.278440    1.309785    1.324862    1.336269
     1.359739    1.588860    1.797209    1.901175    1.913485    2.035051    2.482037    2.533421    2.608642    2.777420
     2.782407    2.857318    2.867918    3.046082    3.108081    3.284321    3.393368    3.419807    3.493896    3.793958
     4.181778    4.254359    4.260902    4.338132    4.499866    4.744929    4.848299    5.010291    5.033490    5.238908
     5.926062    6.230325    6.386799

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   47 emc=    -93.9735481766 demc= 6.8522E-06 wnorm= 2.0203E-03 knorm= 1.0031E-02 apxde= 3.4967E-06    *not conv.*     

               starting mcscf iteration...  48
 !timer:                                 cpu_time=    15.573 walltime=    15.719

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368298 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0812859689     -116.9442478375        0.0000000000        0.0000001250
    2       -93.9357583836     -116.7987202522        0.0000000000        0.0000001250
    3       -93.9036214926     -116.7665833611        0.0000000000        0.0000001250
    4       -93.7772366026     -116.6401984712        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  2.582034906356280E-004
 Total number of micro iterations:   11

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.99994921 pnorm= 0.0000E+00 rznorm= 2.1016E-07 rpnorm= 0.0000E+00 noldr= 11 nnewr= 11 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.714278  -23.221365   -2.319948   -1.895275   -1.852142

 qvv(*) eigenvalues. symmetry block  1
    -0.417875   -0.199293   -0.162608   -0.123157   -0.085481   -0.040461   -0.023756    0.024240    0.042707    0.067468
     0.082486    0.094428    0.142887    0.200167    0.219701    0.277126    0.322243    0.393671    0.418524    0.449412
     0.474366    0.497774    0.525926    0.561125    0.618731    0.646906    0.667574    0.730868    0.802559    0.815959
     0.907422    0.922128    0.956975    0.967872    1.003108    1.196889    1.278506    1.309869    1.324874    1.336290
     1.359795    1.588912    1.797274    1.901283    1.913535    2.035084    2.482010    2.533493    2.608687    2.777472
     2.782474    2.857371    2.867950    3.045945    3.107941    3.284011    3.393335    3.419797    3.493952    3.793706
     4.181764    4.254428    4.260967    4.338184    4.499660    4.744787    4.848276    5.010434    5.033442    5.238975
     5.926009    6.230394    6.386938

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   48 emc=    -93.9735552817 demc= 7.1051E-06 wnorm= 2.0656E-03 knorm= 1.0078E-02 apxde= 3.6259E-06    *not conv.*     

               starting mcscf iteration...  49
 !timer:                                 cpu_time=    15.892 walltime=    16.038

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368298 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0812985857     -116.9442604542        0.0000000000        0.0000001250
    2       -93.9357935964     -116.7987554650        0.0000000000        0.0000001250
    3       -93.9035957644     -116.7665576330        0.0000000000        0.0000001250
    4       -93.7775100435     -116.6404719121        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  2.639942326720670E-004
 Total number of micro iterations:   11

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.99994873 pnorm= 0.0000E+00 rznorm= 2.1606E-07 rpnorm= 0.0000E+00 noldr= 11 nnewr= 11 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.714290  -23.221292   -2.319879   -1.901994   -1.852094

 qvv(*) eigenvalues. symmetry block  1
    -0.417808   -0.199293   -0.162624   -0.123148   -0.085501   -0.040423   -0.023720    0.024281    0.042729    0.067496
     0.082479    0.094440    0.142915    0.200190    0.219649    0.277168    0.322229    0.393660    0.418489    0.449366
     0.474317    0.497779    0.525968    0.561164    0.618780    0.646960    0.667645    0.730899    0.802530    0.815890
     0.907413    0.922094    0.957017    0.967911    1.003067    1.196956    1.278573    1.309954    1.324886    1.336311
     1.359852    1.588965    1.797340    1.901392    1.913586    2.035117    2.481983    2.533566    2.608732    2.777526
     2.782543    2.857426    2.867982    3.045806    3.107797    3.283702    3.393300    3.419785    3.494009    3.793452
     4.181749    4.254499    4.261033    4.338237    4.499451    4.744642    4.848253    5.010575    5.033390    5.239044
     5.925960    6.230463    6.387074

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   49 emc=    -93.9735626488 demc= 7.3671E-06 wnorm= 2.1120E-03 knorm= 1.0126E-02 apxde= 3.7596E-06    *not conv.*     

               starting mcscf iteration...  50
 !timer:                                 cpu_time=    16.210 walltime=    16.356

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368297 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0813111740     -116.9442730425        0.0000000000        0.0000001250
    2       -93.9358299868     -116.7987918554        0.0000000000        0.0000001250
    3       -93.9035697003     -116.7665315689        0.0000000000        0.0000001250
    4       -93.7777838919     -116.6407457605        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  2.698949196148492E-004
 Total number of micro iterations:   11

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.99994825 pnorm= 0.0000E+00 rznorm= 2.2209E-07 rpnorm= 0.0000E+00 noldr= 11 nnewr= 11 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.714303  -23.221219   -2.319808   -1.908747   -1.852040

 qvv(*) eigenvalues. symmetry block  1
    -0.417741   -0.199292   -0.162640   -0.123138   -0.085521   -0.040384   -0.023684    0.024324    0.042752    0.067525
     0.082471    0.094452    0.142943    0.200213    0.219596    0.277210    0.322215    0.393649    0.418454    0.449319
     0.474267    0.497783    0.526011    0.561203    0.618829    0.647015    0.667717    0.730929    0.802501    0.815820
     0.907404    0.922058    0.957060    0.967951    1.003024    1.197022    1.278641    1.310040    1.324899    1.336331
     1.359910    1.589019    1.797405    1.901502    1.913639    2.035150    2.481954    2.533640    2.608779    2.777581
     2.782613    2.857481    2.868016    3.045667    3.107650    3.283392    3.393262    3.419770    3.494068    3.793198
     4.181733    4.254571    4.261100    4.338291    4.499238    4.744493    4.848232    5.010716    5.033336    5.239114
     5.925914    6.230534    6.387207

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   50 emc=    -93.9735702870 demc= 7.6382E-06 wnorm= 2.1592E-03 knorm= 1.0173E-02 apxde= 3.8976E-06    *not conv.*     

               starting mcscf iteration...  51
 !timer:                                 cpu_time=    16.532 walltime=    16.679

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368298 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0813237225     -116.9442855910        0.0000000000        0.0000001250
    2       -93.9358675918     -116.7988294603        0.0000000000        0.0000001250
    3       -93.9035433007     -116.7665051693        0.0000000000        0.0000001250
    4       -93.7780580051     -116.6410198737        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  2.758956315691307E-004
 Total number of micro iterations:   11

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.99994778 pnorm= 0.0000E+00 rznorm= 2.2810E-07 rpnorm= 0.0000E+00 noldr= 11 nnewr= 11 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.714316  -23.221145   -2.319736   -1.915530   -1.851980

 qvv(*) eigenvalues. symmetry block  1
    -0.417672   -0.199291   -0.162656   -0.123128   -0.085541   -0.040345   -0.023647    0.024368    0.042774    0.067554
     0.082463    0.094464    0.142970    0.200237    0.219542    0.277253    0.322200    0.393637    0.418419    0.449272
     0.474216    0.497787    0.526054    0.561243    0.618879    0.647070    0.667790    0.730961    0.802472    0.815749
     0.907394    0.922022    0.957104    0.967991    1.002980    1.197089    1.278711    1.310127    1.324911    1.336352
     1.359968    1.589075    1.797472    1.901613    1.913692    2.035183    2.481925    2.533715    2.608826    2.777636
     2.782684    2.857537    2.868050    3.045527    3.107499    3.283082    3.393222    3.419753    3.494127    3.792941
     4.181716    4.254644    4.261169    4.338346    4.499023    4.744341    4.848212    5.010855    5.033278    5.239185
     5.925872    6.230605    6.387335

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   51 emc=    -93.9735782050 demc= 7.9179E-06 wnorm= 2.2072E-03 knorm= 1.0220E-02 apxde= 4.0398E-06    *not conv.*     

               starting mcscf iteration...  52
 !timer:                                 cpu_time=    16.853 walltime=    17.000

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368296 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0813362189     -116.9442980875        0.0000000000        0.0000001250
    2       -93.9359064464     -116.7988683149        0.0000000000        0.0000001250
    3       -93.9035165676     -116.7664784361        0.0000000000        0.0000001250
    4       -93.7783322212     -116.6412940898        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  2.819851511895624E-004
 Total number of micro iterations:   10

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.99994731 pnorm= 0.0000E+00 rznorm= 6.2018E-07 rpnorm= 0.0000E+00 noldr= 10 nnewr= 10 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.714329  -23.221069   -2.319663   -1.922338   -1.851916

 qvv(*) eigenvalues. symmetry block  1
    -0.417603   -0.199291   -0.162672   -0.123118   -0.085562   -0.040306   -0.023611    0.024412    0.042797    0.067584
     0.082455    0.094477    0.142998    0.200262    0.219487    0.277297    0.322185    0.393624    0.418384    0.449223
     0.474165    0.497790    0.526099    0.561283    0.618930    0.647127    0.667864    0.730993    0.802443    0.815678
     0.907385    0.921985    0.957149    0.968032    1.002936    1.197156    1.278781    1.310215    1.324923    1.336371
     1.360027    1.589131    1.797538    1.901726    1.913747    2.035215    2.481895    2.533791    2.608873    2.777692
     2.782755    2.857593    2.868085    3.045386    3.107345    3.282773    3.393179    3.419733    3.494188    3.792684
     4.181699    4.254718    4.261238    4.338402    4.498803    4.744185    4.848193    5.010992    5.033218    5.239257
     5.925833    6.230678    6.387460

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   52 emc=    -93.9735864109 demc= 8.2060E-06 wnorm= 2.2559E-03 knorm= 1.0265E-02 apxde= 4.1860E-06    *not conv.*     

               starting mcscf iteration...  53
 !timer:                                 cpu_time=    17.174 walltime=    17.320

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368298 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0813486491     -116.9443105176        0.0000000000        0.0000001250
    2       -93.9359465870     -116.7989084556        0.0000000000        0.0000001250
    3       -93.9034895021     -116.7664513706        0.0000000000        0.0000001250
    4       -93.7786063534     -116.6415682219        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  2.881508072936539E-004
 Total number of micro iterations:   10

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.99994685 pnorm= 0.0000E+00 rznorm= 6.2778E-07 rpnorm= 0.0000E+00 noldr= 10 nnewr= 10 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.714343  -23.220992   -2.319589   -1.929168   -1.851848

 qvv(*) eigenvalues. symmetry block  1
    -0.417532   -0.199290   -0.162688   -0.123107   -0.085584   -0.040266   -0.023573    0.024458    0.042820    0.067614
     0.082447    0.094490    0.143025    0.200288    0.219432    0.277342    0.322170    0.393611    0.418348    0.449175
     0.474113    0.497794    0.526143    0.561325    0.618982    0.647184    0.667939    0.731025    0.802413    0.815607
     0.907376    0.921948    0.957194    0.968074    1.002891    1.197223    1.278854    1.310304    1.324936    1.336390
     1.360087    1.589189    1.797605    1.901840    1.913802    2.035247    2.481863    2.533869    2.608922    2.777750
     2.782828    2.857651    2.868120    3.045244    3.107187    3.282464    3.393134    3.419711    3.494249    3.792426
     4.181680    4.254793    4.261309    4.338459    4.498581    4.744025    4.848175    5.011127    5.033153    5.239331
     5.925798    6.230752    6.387580

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   53 emc=    -93.9735949127 demc= 8.5018E-06 wnorm= 2.3052E-03 knorm= 1.0310E-02 apxde= 4.3358E-06    *not conv.*     

               starting mcscf iteration...  54
 !timer:                                 cpu_time=    17.492 walltime=    17.638

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368298 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0813610007     -116.9443228692        0.0000000000        0.0000001250
    2       -93.9359880372     -116.7989499058        0.0000000000        0.0000001250
    3       -93.9034621145     -116.7664239831        0.0000000000        0.0000001250
    4       -93.7788802129     -116.6418420814        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  2.943786407883970E-004
 Total number of micro iterations:   10

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.99994642 pnorm= 0.0000E+00 rznorm= 6.3019E-07 rpnorm= 0.0000E+00 noldr= 10 nnewr= 10 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.714358  -23.220914   -2.319513   -1.936013   -1.851778

 qvv(*) eigenvalues. symmetry block  1
    -0.417461   -0.199290   -0.162704   -0.123096   -0.085605   -0.040225   -0.023536    0.024504    0.042843    0.067644
     0.082438    0.094504    0.143052    0.200315    0.219377    0.277387    0.322155    0.393598    0.418312    0.449125
     0.474060    0.497797    0.526189    0.561367    0.619034    0.647242    0.668016    0.731058    0.802383    0.815536
     0.907368    0.921910    0.957240    0.968116    1.002845    1.197291    1.278927    1.310394    1.324948    1.336409
     1.360147    1.589249    1.797673    1.901956    1.913859    2.035278    2.481831    2.533947    2.608971    2.777807
     2.782902    2.857709    2.868156    3.045101    3.107026    3.282157    3.393085    3.419687    3.494311    3.792166
     4.181661    4.254870    4.261380    4.338516    4.498355    4.743862    4.848159    5.011261    5.033086    5.239406
     5.925768    6.230827    6.387696

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   54 emc=    -93.9736037175 demc= 8.8047E-06 wnorm= 2.3550E-03 knorm= 1.0352E-02 apxde= 4.4890E-06    *not conv.*     

               starting mcscf iteration...  55
 !timer:                                 cpu_time=    17.811 walltime=    17.956

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368297 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0813732576     -116.9443351262        0.0000000000        0.0000001250
    2       -93.9360308272     -116.7989926958        0.0000000000        0.0000001250
    3       -93.9034344100     -116.7663962785        0.0000000000        0.0000001250
    4       -93.7791535736     -116.6421154422        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  3.006529081316736E-004
 Total number of micro iterations:   10

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.99994600 pnorm= 0.0000E+00 rznorm= 6.3199E-07 rpnorm= 0.0000E+00 noldr= 10 nnewr= 10 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.714372  -23.220835   -2.319437   -1.942870   -1.851705

 qvv(*) eigenvalues. symmetry block  1
    -0.417388   -0.199289   -0.162720   -0.123085   -0.085628   -0.040184   -0.023498    0.024552    0.042866    0.067675
     0.082430    0.094518    0.143078    0.200343    0.219320    0.277433    0.322139    0.393584    0.418277    0.449075
     0.474006    0.497799    0.526235    0.561410    0.619087    0.647301    0.668094    0.731091    0.802353    0.815465
     0.907359    0.921872    0.957287    0.968159    1.002798    1.197358    1.279002    1.310485    1.324960    1.336427
     1.360209    1.589309    1.797740    1.902073    1.913916    2.035309    2.481798    2.534027    2.609021    2.777866
     2.782978    2.857767    2.868193    3.044958    3.106860    3.281850    3.393034    3.419659    3.494375    3.791907
     4.181640    4.254948    4.261453    4.338574    4.498126    4.743694    4.848145    5.011392    5.033014    5.239482
     5.925741    6.230903    6.387807

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   55 emc=    -93.9736128316 demc= 9.1141E-06 wnorm= 2.4052E-03 knorm= 1.0392E-02 apxde= 4.6450E-06    *not conv.*     

               starting mcscf iteration...  56
 !timer:                                 cpu_time=    18.129 walltime=    18.275

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368298 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0813854038     -116.9443472724        0.0000000000        0.0000001250
    2       -93.9360749803     -116.7990368489        0.0000000000        0.0000001250
    3       -93.9034063978     -116.7663682664        0.0000000000        0.0000001250
    4       -93.7794261959     -116.6423880645        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  3.069563226128498E-004
 Total number of micro iterations:   10

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.99994561 pnorm= 0.0000E+00 rznorm= 6.3305E-07 rpnorm= 0.0000E+00 noldr= 10 nnewr= 10 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.714388  -23.220754   -2.319359   -1.949734   -1.851631

 qvv(*) eigenvalues. symmetry block  1
    -0.417315   -0.199289   -0.162736   -0.123073   -0.085650   -0.040143   -0.023460    0.024600    0.042889    0.067706
     0.082420    0.094532    0.143104    0.200371    0.219263    0.277480    0.322123    0.393569    0.418241    0.449024
     0.473952    0.497802    0.526282    0.561454    0.619142    0.647361    0.668172    0.731125    0.802323    0.815393
     0.907351    0.921833    0.957335    0.968203    1.002750    1.197425    1.279077    1.310577    1.324972    1.336445
     1.360270    1.589371    1.797808    1.902191    1.913974    2.035339    2.481764    2.534107    2.609071    2.777926
     2.783054    2.857827    2.868231    3.044814    3.106691    3.281545    3.392979    3.419630    3.494439    3.791647
     4.181619    4.255027    4.261527    4.338634    4.497894    4.743523    4.848131    5.011521    5.032940    5.239559
     5.925719    6.230981    6.387913

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   56 emc=    -93.9736222606 demc= 9.4290E-06 wnorm= 2.4557E-03 knorm= 1.0430E-02 apxde= 4.8034E-06    *not conv.*     

               starting mcscf iteration...  57
 !timer:                                 cpu_time=    18.448 walltime=    18.594

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368298 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0813974223     -116.9443592909        0.0000000000        0.0000001250
    2       -93.9361205153     -116.7990823839        0.0000000000        0.0000001250
    3       -93.9033780898     -116.7663399583        0.0000000000        0.0000001250
    4       -93.7796978204     -116.6426596890        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  3.132698899129876E-004
 Total number of micro iterations:   10

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.99994524 pnorm= 0.0000E+00 rznorm= 6.3332E-07 rpnorm= 0.0000E+00 noldr= 10 nnewr= 10 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.714404  -23.220673   -2.319281   -1.956599   -1.851554

 qvv(*) eigenvalues. symmetry block  1
    -0.417241   -0.199288   -0.162753   -0.123060   -0.085674   -0.040101   -0.023422    0.024650    0.042912    0.067738
     0.082411    0.094547    0.143130    0.200400    0.219206    0.277528    0.322106    0.393553    0.418205    0.448973
     0.473897    0.497804    0.526329    0.561498    0.619196    0.647421    0.668252    0.731159    0.802293    0.815321
     0.907343    0.921793    0.957383    0.968248    1.002702    1.197493    1.279155    1.310669    1.324984    1.336462
     1.360333    1.589434    1.797876    1.902311    1.914033    2.035369    2.481729    2.534189    2.609122    2.777986
     2.783131    2.857887    2.868269    3.044671    3.106519    3.281242    3.392920    3.419597    3.494505    3.791386
     4.181596    4.255108    4.261601    4.338693    4.497659    4.743348    4.848120    5.011648    5.032861    5.239637
     5.925702    6.231059    6.388013

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   57 emc=    -93.9736320091 demc= 9.7485E-06 wnorm= 2.5062E-03 knorm= 1.0465E-02 apxde= 4.9637E-06    *not conv.*     

               starting mcscf iteration...  58
 !timer:                                 cpu_time=    18.768 walltime=    18.914

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368296 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0814092956     -116.9443711641        0.0000000000        0.0000001250
    2       -93.9361674461     -116.7991293147        0.0000000000        0.0000001250
    3       -93.9033494996     -116.7663113682        0.0000000000        0.0000001250
    4       -93.7799681685     -116.6429300371        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  3.195728732167807E-004
 Total number of micro iterations:    9

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.99994492 pnorm= 0.0000E+00 rznorm= 7.8555E-07 rpnorm= 0.0000E+00 noldr=  9 nnewr=  9 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.714420  -23.220591   -2.319201   -1.963459   -1.851475

 qvv(*) eigenvalues. symmetry block  1
    -0.417165   -0.199288   -0.162769   -0.123048   -0.085697   -0.040059   -0.023384    0.024700    0.042936    0.067769
     0.082401    0.094562    0.143155    0.200430    0.219148    0.277576    0.322089    0.393537    0.418170    0.448921
     0.473841    0.497805    0.526377    0.561544    0.619252    0.647482    0.668333    0.731194    0.802263    0.815250
     0.907335    0.921753    0.957433    0.968293    1.002653    1.197560    1.279233    1.310763    1.324996    1.336478
     1.360396    1.589498    1.797945    1.902432    1.914092    2.035398    2.481693    2.534271    2.609174    2.778047
     2.783209    2.857947    2.868308    3.044526    3.106342    3.280941    3.392858    3.419562    3.494571    3.791126
     4.181572    4.255189    4.261677    4.338754    4.497420    4.743170    4.848110    5.011772    5.032779    5.239717
     5.925689    6.231138    6.388108

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   58 emc=    -93.9736420805 demc= 1.0071E-05 wnorm= 2.5566E-03 knorm= 1.0495E-02 apxde= 5.1252E-06    *not conv.*     

               starting mcscf iteration...  59
 !timer:                                 cpu_time=    19.087 walltime=    19.232

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368295 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0814209997     -116.9443828683        0.0000000000        0.0000001250
    2       -93.9362157690     -116.7991776375        0.0000000000        0.0000001250
    3       -93.9033206612     -116.7662825298        0.0000000000        0.0000001250
    4       -93.7802369535     -116.6431988220        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  3.258428145019977E-004
 Total number of micro iterations:    9

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.99994464 pnorm= 0.0000E+00 rznorm= 8.1636E-07 rpnorm= 0.0000E+00 noldr=  9 nnewr=  9 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.714437  -23.220507   -2.319121   -1.970308   -1.851395

 qvv(*) eigenvalues. symmetry block  1
    -0.417089   -0.199288   -0.162785   -0.123035   -0.085721   -0.040016   -0.023345    0.024751    0.042959    0.067801
     0.082392    0.094578    0.143180    0.200461    0.219089    0.277625    0.322072    0.393520    0.418134    0.448868
     0.473785    0.497806    0.526425    0.561590    0.619308    0.647544    0.668415    0.731230    0.802232    0.815178
     0.907327    0.921712    0.957483    0.968338    1.002603    1.197627    1.279312    1.310857    1.325007    1.336494
     1.360459    1.589563    1.798013    1.902555    1.914152    2.035427    2.481657    2.534355    2.609226    2.778108
     2.783287    2.858009    2.868348    3.044382    3.106162    3.280643    3.392793    3.419525    3.494638    3.790866
     4.181547    4.255272    4.261754    4.338816    4.497179    4.742987    4.848102    5.011893    5.032693    5.239797
     5.925682    6.231218    6.388197

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   59 emc=    -93.9736524766 demc= 1.0396E-05 wnorm= 2.6067E-03 knorm= 1.0523E-02 apxde= 5.2873E-06    *not conv.*     

               starting mcscf iteration...  60
 !timer:                                 cpu_time=    19.406 walltime=    19.551

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368297 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0814325263     -116.9443943948        0.0000000000        0.0000001250
    2       -93.9362655094     -116.7992273779        0.0000000000        0.0000001250
    3       -93.9032915591     -116.7662534277        0.0000000000        0.0000001250
    4       -93.7805038382     -116.6434657068        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  3.320553485781843E-004
 Total number of micro iterations:    9

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.99994440 pnorm= 0.0000E+00 rznorm= 8.1400E-07 rpnorm= 0.0000E+00 noldr=  9 nnewr=  9 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.714456  -23.220423   -2.319040   -1.977140   -1.851314

 qvv(*) eigenvalues. symmetry block  1
    -0.417013   -0.199287   -0.162802   -0.123021   -0.085746   -0.039973   -0.023306    0.024803    0.042983    0.067833
     0.082382    0.094594    0.143204    0.200493    0.219031    0.277675    0.322055    0.393503    0.418099    0.448816
     0.473727    0.497807    0.526474    0.561636    0.619366    0.647607    0.668497    0.731265    0.802202    0.815106
     0.907320    0.921671    0.957533    0.968385    1.002552    1.197694    1.279393    1.310952    1.325019    1.336509
     1.360523    1.589630    1.798081    1.902679    1.914212    2.035455    2.481619    2.534439    2.609278    2.778171
     2.783367    2.858070    2.868388    3.044238    3.105978    3.280348    3.392723    3.419485    3.494707    3.790607
     4.181522    4.255355    4.261832    4.338878    4.496935    4.742801    4.848095    5.012010    5.032604    5.239879
     5.925680    6.231300    6.388279

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   60 emc=    -93.9736631982 demc= 1.0722E-05 wnorm= 2.6564E-03 knorm= 1.0545E-02 apxde= 5.4490E-06    *not conv.*     

               starting mcscf iteration...  61
 !timer:                                 cpu_time=    19.725 walltime=    19.869

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368297 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0814438511     -116.9444057197        0.0000000000        0.0000001250
    2       -93.9363166501     -116.7992785187        0.0000000000        0.0000001250
    3       -93.9032622315     -116.7662241001        0.0000000000        0.0000001250
    4       -93.7807685014     -116.6437303699        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  3.381845502497226E-004
 Total number of micro iterations:    9

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.99994421 pnorm= 0.0000E+00 rznorm= 8.0924E-07 rpnorm= 0.0000E+00 noldr=  9 nnewr=  9 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.714475  -23.220338   -2.318957   -1.983949   -1.851231

 qvv(*) eigenvalues. symmetry block  1
    -0.416935   -0.199287   -0.162818   -0.123008   -0.085770   -0.039929   -0.023267    0.024856    0.043006    0.067866
     0.082371    0.094610    0.143227    0.200525    0.218971    0.277725    0.322037    0.393485    0.418063    0.448762
     0.473670    0.497808    0.526523    0.561684    0.619423    0.647670    0.668581    0.731301    0.802172    0.815035
     0.907313    0.921629    0.957585    0.968431    1.002501    1.197761    1.279475    1.311047    1.325030    1.336523
     1.360587    1.589698    1.798149    1.902804    1.914272    2.035482    2.481580    2.534524    2.609331    2.778233
     2.783448    2.858133    2.868429    3.044095    3.105791    3.280057    3.392650    3.419442    3.494775    3.790349
     4.181494    4.255440    4.261910    4.338941    4.496689    4.742612    4.848090    5.012124    5.032510    5.239962
     5.925683    6.231382    6.388354

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   61 emc=    -93.9736742442 demc= 1.1046E-05 wnorm= 2.7055E-03 knorm= 1.0563E-02 apxde= 5.6096E-06    *not conv.*     

               starting mcscf iteration...  62
 !timer:                                 cpu_time=    20.043 walltime=    20.188

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368298 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0814549541     -116.9444168226        0.0000000000        0.0000001250
    2       -93.9363691780     -116.7993310466        0.0000000000        0.0000001250
    3       -93.9032327031     -116.7661945716        0.0000000000        0.0000001250
    4       -93.7810305939     -116.6439924625        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  3.442027157899698E-004
 Total number of micro iterations:    9

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.99994408 pnorm= 0.0000E+00 rznorm= 8.0311E-07 rpnorm= 0.0000E+00 noldr=  9 nnewr=  9 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.714495  -23.220252   -2.318874   -1.990728   -1.851147

 qvv(*) eigenvalues. symmetry block  1
    -0.416857   -0.199287   -0.162835   -0.122994   -0.085796   -0.039886   -0.023228    0.024910    0.043029    0.067898
     0.082361    0.094626    0.143250    0.200558    0.218912    0.277776    0.322019    0.393466    0.418028    0.448709
     0.473611    0.497808    0.526573    0.561731    0.619482    0.647734    0.668666    0.731338    0.802142    0.814964
     0.907306    0.921586    0.957637    0.968479    1.002449    1.197827    1.279557    1.311143    1.325041    1.336537
     1.360652    1.589767    1.798218    1.902931    1.914332    2.035508    2.481540    2.534610    2.609385    2.778297
     2.783529    2.858195    2.868470    3.043952    3.105600    3.279770    3.392572    3.419397    3.494845    3.790092
     4.181466    4.255525    4.261989    4.339004    4.496441    4.742419    4.848087    5.012235    5.032413    5.240045
     5.925691    6.231464    6.388423

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   62 emc=    -93.9736856117 demc= 1.1367E-05 wnorm= 2.7536E-03 knorm= 1.0576E-02 apxde= 5.7681E-06    *not conv.*     

               starting mcscf iteration...  63
 !timer:                                 cpu_time=    20.361 walltime=    20.505

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368298 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0814658152     -116.9444276837        0.0000000000        0.0000001250
    2       -93.9364230715     -116.7993849401        0.0000000000        0.0000001250
    3       -93.9032030010     -116.7661648696        0.0000000000        0.0000001250
    4       -93.7812897528     -116.6442516213        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  3.500805851941273E-004
 Total number of micro iterations:    9

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.99994400 pnorm= 0.0000E+00 rznorm= 7.9564E-07 rpnorm= 0.0000E+00 noldr=  9 nnewr=  9 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.714516  -23.220165   -2.318791   -1.997470   -1.851061

 qvv(*) eigenvalues. symmetry block  1
    -0.416779   -0.199287   -0.162851   -0.122979   -0.085821   -0.039842   -0.023189    0.024964    0.043053    0.067931
     0.082350    0.094643    0.143272    0.200592    0.218852    0.277828    0.322001    0.393447    0.417994    0.448655
     0.473552    0.497807    0.526623    0.561780    0.619540    0.647798    0.668751    0.731375    0.802112    0.814893
     0.907299    0.921544    0.957689    0.968526    1.002396    1.197892    1.279641    1.311239    1.325051    1.336549
     1.360717    1.589836    1.798286    1.903059    1.914393    2.035533    2.481500    2.534696    2.609439    2.778361
     2.783610    2.858258    2.868513    3.043809    3.105406    3.279487    3.392491    3.419350    3.494916    3.789836
     4.181437    4.255612    4.262069    4.339068    4.496190    4.742222    4.848087    5.012341    5.032311    5.240130
     5.925706    6.231548    6.388484

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   63 emc=    -93.9736972959 demc= 1.1684E-05 wnorm= 2.8006E-03 knorm= 1.0583E-02 apxde= 5.9234E-06    *not conv.*     

               starting mcscf iteration...  64
 !timer:                                 cpu_time=    20.680 walltime=    20.824

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368298 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0814764144     -116.9444382830        0.0000000000        0.0000001250
    2       -93.9364782997     -116.7994401682        0.0000000000        0.0000001250
    3       -93.9031731553     -116.7661350239        0.0000000000        0.0000001250
    4       -93.7815456021     -116.6445074706        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  3.557874532705307E-004
 Total number of micro iterations:    9

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.99994400 pnorm= 0.0000E+00 rznorm= 7.8683E-07 rpnorm= 0.0000E+00 noldr=  9 nnewr=  9 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.714538  -23.220077   -2.318706   -2.004166   -1.850975

 qvv(*) eigenvalues. symmetry block  1
    -0.416699   -0.199287   -0.162867   -0.122964   -0.085847   -0.039798   -0.023150    0.025020    0.043076    0.067964
     0.082339    0.094660    0.143293    0.200627    0.218792    0.277880    0.321982    0.393427    0.417959    0.448600
     0.473493    0.497807    0.526673    0.561829    0.619600    0.647863    0.668837    0.731412    0.802083    0.814822
     0.907293    0.921500    0.957742    0.968575    1.002343    1.197957    1.279725    1.311336    1.325062    1.336561
     1.360782    1.589907    1.798353    1.903188    1.914453    2.035557    2.481458    2.534783    2.609493    2.778425
     2.783693    2.858321    2.868555    3.043667    3.105208    3.279210    3.392404    3.419300    3.494987    3.789583
     4.181406    4.255699    4.262150    4.339133    4.495937    4.742022    4.848088    5.012443    5.032206    5.240215
     5.925726    6.231632    6.388537

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   64 emc=    -93.9737092898 demc= 1.1994E-05 wnorm= 2.8463E-03 knorm= 1.0583E-02 apxde= 6.0745E-06    *not conv.*     

               starting mcscf iteration...  65
 !timer:                                 cpu_time=    20.998 walltime=    21.142

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368298 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0814867318     -116.9444486004        0.0000000000        0.0000001250
    2       -93.9365348219     -116.7994966905        0.0000000000        0.0000001250
    3       -93.9031431989     -116.7661050674        0.0000000000        0.0000001250
    4       -93.7817977554     -116.6447596239        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  3.612913484065350E-004
 Total number of micro iterations:    9

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.99994406 pnorm= 0.0000E+00 rznorm= 7.7666E-07 rpnorm= 0.0000E+00 noldr=  9 nnewr=  9 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.714562  -23.219989   -2.318622   -2.010811   -1.850888

 qvv(*) eigenvalues. symmetry block  1
    -0.416620   -0.199287   -0.162883   -0.122949   -0.085874   -0.039753   -0.023111    0.025076    0.043099    0.067997
     0.082328    0.094678    0.143314    0.200663    0.218732    0.277933    0.321964    0.393406    0.417925    0.448546
     0.473433    0.497805    0.526723    0.561879    0.619660    0.647928    0.668923    0.731450    0.802053    0.814752
     0.907288    0.921457    0.957796    0.968623    1.002289    1.198021    1.279811    1.311432    1.325072    1.336571
     1.360847    1.589979    1.798420    1.903319    1.914513    2.035580    2.481416    2.534870    2.609548    2.778490
     2.783776    2.858385    2.868599    3.043527    3.105007    3.278938    3.392313    3.419248    3.495058    3.789332
     4.181375    4.255787    4.262231    4.339198    4.495683    4.741819    4.848091    5.012541    5.032097    5.240301
     5.925753    6.231717    6.388582

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   65 emc=    -93.9737215842 demc= 1.2294E-05 wnorm= 2.8903E-03 knorm= 1.0577E-02 apxde= 6.2202E-06    *not conv.*     

               starting mcscf iteration...  66
 !timer:                                 cpu_time=    21.317 walltime=    21.460

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368297 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0814967480     -116.9444586166        0.0000000000        0.0000001250
    2       -93.9365925874     -116.7995544559        0.0000000000        0.0000001250
    3       -93.9031131674     -116.7660750359        0.0000000000        0.0000001250
    4       -93.7820458182     -116.6450076867        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  3.665592549707244E-004
 Total number of micro iterations:    9

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.99994420 pnorm= 0.0000E+00 rznorm= 7.6516E-07 rpnorm= 0.0000E+00 noldr=  9 nnewr=  9 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.714586  -23.219901   -2.318536   -2.017397   -1.850801

 qvv(*) eigenvalues. symmetry block  1
    -0.416540   -0.199287   -0.162899   -0.122933   -0.085900   -0.039709   -0.023072    0.025133    0.043122    0.068030
     0.082317    0.094696    0.143333    0.200699    0.218671    0.277986    0.321945    0.393385    0.417892    0.448491
     0.473373    0.497804    0.526774    0.561929    0.619720    0.647993    0.669010    0.731487    0.802024    0.814683
     0.907282    0.921413    0.957849    0.968672    1.002235    1.198084    1.279897    1.311529    1.325081    1.336581
     1.360912    1.590052    1.798487    1.903450    1.914572    2.035602    2.481372    2.534958    2.609603    2.778555
     2.783859    2.858448    2.868643    3.043388    3.104804    3.278673    3.392218    3.419193    3.495130    3.789085
     4.181342    4.255875    4.262312    4.339263    4.495428    4.741613    4.848096    5.012634    5.031983    5.240387
     5.925785    6.231802    6.388619

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   66 emc=    -93.9737341676 demc= 1.2583E-05 wnorm= 2.9325E-03 knorm= 1.0564E-02 apxde= 6.3593E-06    *not conv.*     

               starting mcscf iteration...  67
 !timer:                                 cpu_time=    21.634 walltime=    21.778

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368297 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0815064439     -116.9444683125        0.0000000000        0.0000001250
    2       -93.9366515343     -116.7996134028        0.0000000000        0.0000001250
    3       -93.9030830994     -116.7660449679        0.0000000000        0.0000001250
    4       -93.7822893905     -116.6452512590        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  3.715573819393393E-004
 Total number of micro iterations:    9

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.99994441 pnorm= 0.0000E+00 rznorm= 7.5233E-07 rpnorm= 0.0000E+00 noldr=  9 nnewr=  9 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.714613  -23.219812   -2.318451   -2.023914   -1.850712

 qvv(*) eigenvalues. symmetry block  1
    -0.416460   -0.199287   -0.162915   -0.122917   -0.085927   -0.039664   -0.023033    0.025190    0.043145    0.068063
     0.082306    0.094714    0.143352    0.200736    0.218611    0.278039    0.321926    0.393363    0.417859    0.448436
     0.473312    0.497802    0.526825    0.561980    0.619780    0.648059    0.669097    0.731525    0.801996    0.814614
     0.907277    0.921369    0.957904    0.968721    1.002181    1.198147    1.279984    1.311626    1.325091    1.336590
     1.360977    1.590125    1.798553    1.903583    1.914631    2.035623    2.481328    2.535046    2.609658    2.778620
     2.783942    2.858512    2.868687    3.043250    3.104597    3.278414    3.392117    3.419137    3.495203    3.788840
     4.181308    4.255964    4.262394    4.339328    4.495171    4.741405    4.848103    5.012721    5.031866    5.240474
     5.925824    6.231888    6.388648

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   67 emc=    -93.9737470259 demc= 1.2858E-05 wnorm= 2.9725E-03 knorm= 1.0544E-02 apxde= 6.4905E-06    *not conv.*     

               starting mcscf iteration...  68
 !timer:                                 cpu_time=    21.953 walltime=    22.096

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368295 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0815158015     -116.9444776700        0.0000000000        0.0000001250
    2       -93.9367115898     -116.7996734584        0.0000000000        0.0000001250
    3       -93.9030530360     -116.7660149046        0.0000000000        0.0000001250
    4       -93.7825280699     -116.6454899385        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  3.762514769457748E-004
 Total number of micro iterations:    9

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.99994471 pnorm= 0.0000E+00 rznorm= 7.3822E-07 rpnorm= 0.0000E+00 noldr=  9 nnewr=  9 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.714640  -23.219723   -2.318365   -2.030357   -1.850624

 qvv(*) eigenvalues. symmetry block  1
    -0.416380   -0.199287   -0.162931   -0.122901   -0.085954   -0.039620   -0.022995    0.025248    0.043168    0.068096
     0.082295    0.094732    0.143369    0.200774    0.218551    0.278093    0.321907    0.393340    0.417827    0.448381
     0.473251    0.497799    0.526876    0.562031    0.619841    0.648125    0.669185    0.731563    0.801968    0.814546
     0.907273    0.921324    0.957958    0.968771    1.002126    1.198208    1.280071    1.311722    1.325099    1.336597
     1.361041    1.590200    1.798618    1.903717    1.914689    2.035642    2.481284    2.535134    2.609713    2.778686
     2.784026    2.858575    2.868732    3.043114    3.104388    3.278163    3.392012    3.419078    3.495276    3.788599
     4.181273    4.256054    4.262477    4.339394    4.494914    4.741193    4.848112    5.012804    5.031746    5.240562
     5.925870    6.231973    6.388667

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   68 emc=    -93.9737601424 demc= 1.3117E-05 wnorm= 3.0100E-03 knorm= 1.0516E-02 apxde= 6.6126E-06    *not conv.*     

               starting mcscf iteration...  69
 !timer:                                 cpu_time=    22.270 walltime=    22.414

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368297 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0815248035     -116.9444866721        0.0000000000        0.0000001250
    2       -93.9367726699     -116.7997345384        0.0000000000        0.0000001250
    3       -93.9030230210     -116.7659848895        0.0000000000        0.0000001250
    4       -93.7827614550     -116.6457233235        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  3.806071847593282E-004
 Total number of micro iterations:    9

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.99994509 pnorm= 0.0000E+00 rznorm= 7.2287E-07 rpnorm= 0.0000E+00 noldr=  9 nnewr=  9 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.714669  -23.219633   -2.318279   -2.036717   -1.850535

 qvv(*) eigenvalues. symmetry block  1
    -0.416300   -0.199287   -0.162946   -0.122885   -0.085982   -0.039576   -0.022956    0.025306    0.043190    0.068129
     0.082283    0.094750    0.143386    0.200812    0.218491    0.278147    0.321887    0.393317    0.417796    0.448327
     0.473190    0.497796    0.526927    0.562082    0.619902    0.648191    0.669273    0.731602    0.801940    0.814480
     0.907269    0.921280    0.958013    0.968820    1.002070    1.198268    1.280158    1.311818    1.325108    1.336603
     1.361106    1.590274    1.798683    1.903851    1.914745    2.035660    2.481238    2.535222    2.609768    2.778751
     2.784110    2.858639    2.868777    3.042981    3.104176    3.277920    3.391902    3.419018    3.495349    3.788363
     4.181237    4.256143    4.262559    4.339460    4.494657    4.740980    4.848124    5.012881    5.031621    5.240650
     5.925921    6.232059    6.388678

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   69 emc=    -93.9737734981 demc= 1.3356E-05 wnorm= 3.0449E-03 knorm= 1.0479E-02 apxde= 6.7241E-06    *not conv.*     

               starting mcscf iteration...  70
 !timer:                                 cpu_time=    22.590 walltime=    22.733

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368298 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0815334340     -116.9444953026        0.0000000000        0.0000001250
    2       -93.9368346790     -116.7997965476        0.0000000000        0.0000001250
    3       -93.9029931003     -116.7659549688        0.0000000000        0.0000001250
    4       -93.7829891487     -116.6459510173        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  3.845904465539758E-004
 Total number of micro iterations:    9

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.99994556 pnorm= 0.0000E+00 rznorm= 7.0636E-07 rpnorm= 0.0000E+00 noldr=  9 nnewr=  9 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.714700  -23.219544   -2.318193   -2.042986   -1.850446

 qvv(*) eigenvalues. symmetry block  1
    -0.416220   -0.199287   -0.162961   -0.122868   -0.086009   -0.039531   -0.022919    0.025365    0.043212    0.068161
     0.082272    0.094768    0.143402    0.200851    0.218431    0.278201    0.321868    0.393293    0.417765    0.448272
     0.473129    0.497792    0.526978    0.562133    0.619963    0.648257    0.669360    0.731640    0.801913    0.814414
     0.907265    0.921235    0.958067    0.968870    1.002015    1.198326    1.280246    1.311914    1.325116    1.336608
     1.361170    1.590350    1.798746    1.903986    1.914801    2.035677    2.481192    2.535309    2.609823    2.778816
     2.784193    2.858702    2.868822    3.042849    3.103963    3.277685    3.391786    3.418956    3.495422    3.788131
     4.181200    4.256233    4.262642    4.339526    4.494400    4.740764    4.848138    5.012952    5.031493    5.240738
     5.925980    6.232145    6.388679

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   70 emc=    -93.9737870711 demc= 1.3573E-05 wnorm= 3.0767E-03 knorm= 1.0434E-02 apxde= 6.8240E-06    *not conv.*     

               starting mcscf iteration...  71
 !timer:                                 cpu_time=    22.909 walltime=    23.051

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368297 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0815416787     -116.9445035472        0.0000000000        0.0000001250
    2       -93.9368975105     -116.7998593791        0.0000000000        0.0000001250
    3       -93.9029633220     -116.7659251905        0.0000000000        0.0000001250
    4       -93.7832107625     -116.6461726310        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  3.881679346575865E-004
 Total number of micro iterations:    9

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.99994612 pnorm= 0.0000E+00 rznorm= 6.8876E-07 rpnorm= 0.0000E+00 noldr=  9 nnewr=  9 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.714733  -23.219455   -2.318107   -2.049157   -1.850356

 qvv(*) eigenvalues. symmetry block  1
    -0.416140   -0.199287   -0.162976   -0.122851   -0.086037   -0.039487   -0.022881    0.025424    0.043234    0.068194
     0.082260    0.094787    0.143416    0.200890    0.218371    0.278255    0.321848    0.393268    0.417736    0.448218
     0.473068    0.497788    0.527028    0.562185    0.620024    0.648322    0.669448    0.731678    0.801887    0.814349
     0.907262    0.921191    0.958122    0.968919    1.001960    1.198384    1.280334    1.312008    1.325124    1.336612
     1.361234    1.590426    1.798809    1.904122    1.914856    2.035692    2.481145    2.535397    2.609878    2.778882
     2.784276    2.858765    2.868867    3.042720    3.103747    3.277460    3.391666    3.418893    3.495495    3.787904
     4.181162    4.256323    4.262724    4.339592    4.494143    4.740546    4.848154    5.013017    5.031361    5.240825
     5.926045    6.232231    6.388671

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   71 emc=    -93.9738008371 demc= 1.3766E-05 wnorm= 3.1053E-03 knorm= 1.0381E-02 apxde= 6.9108E-06    *not conv.*     

               starting mcscf iteration...  72
 !timer:                                 cpu_time=    23.229 walltime=    23.373

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368297 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0815495246     -116.9445113932        0.0000000000        0.0000001250
    2       -93.9369610469     -116.7999229154        0.0000000000        0.0000001250
    3       -93.9029337358     -116.7658956044        0.0000000000        0.0000001250
    4       -93.7834259196     -116.6463877882        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  3.913075160139438E-004
 Total number of micro iterations:    9

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.99994676 pnorm= 0.0000E+00 rznorm= 6.7019E-07 rpnorm= 0.0000E+00 noldr=  9 nnewr=  9 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.714767  -23.219366   -2.318022   -2.055222   -1.850267

 qvv(*) eigenvalues. symmetry block  1
    -0.416060   -0.199288   -0.162991   -0.122833   -0.086065   -0.039444   -0.022844    0.025483    0.043255    0.068226
     0.082249    0.094806    0.143429    0.200930    0.218312    0.278309    0.321829    0.393243    0.417707    0.448164
     0.473007    0.497783    0.527079    0.562236    0.620084    0.648388    0.669535    0.731716    0.801861    0.814285
     0.907260    0.921146    0.958176    0.968969    1.001904    1.198439    1.280422    1.312102    1.325131    1.336615
     1.361297    1.590502    1.798870    1.904257    1.914909    2.035706    2.481098    2.535484    2.609933    2.778947
     2.784359    2.858827    2.868913    3.042594    3.103530    3.277244    3.391541    3.418828    3.495568    3.787683
     4.181123    4.256412    4.262806    4.339657    4.493887    4.740328    4.848172    5.013076    5.031227    5.240913
     5.926117    6.232317    6.388653

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   72 emc=    -93.9738147691 demc= 1.3932E-05 wnorm= 3.1305E-03 knorm= 1.0318E-02 apxde= 6.9834E-06    *not conv.*     

               starting mcscf iteration...  73
 !timer:                                 cpu_time=    23.547 walltime=    23.691

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368298 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0815569609     -116.9445188295        0.0000000000        0.0000001250
    2       -93.9370251602     -116.7999870287        0.0000000000        0.0000001250
    3       -93.9029043931     -116.7658662617        0.0000000000        0.0000001250
    4       -93.7836342600     -116.6465961286        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  3.939787347088447E-004
 Total number of micro iterations:    9

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.99994750 pnorm= 0.0000E+00 rznorm= 6.5078E-07 rpnorm= 0.0000E+00 noldr=  9 nnewr=  9 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.714802  -23.219278   -2.317937   -2.061174   -1.850179

 qvv(*) eigenvalues. symmetry block  1
    -0.415982   -0.199288   -0.163006   -0.122816   -0.086093   -0.039400   -0.022808    0.025543    0.043276    0.068258
     0.082237    0.094824    0.143442    0.200970    0.218254    0.278363    0.321809    0.393218    0.417680    0.448110
     0.472945    0.497778    0.527129    0.562288    0.620145    0.648453    0.669622    0.731754    0.801837    0.814223
     0.907258    0.921101    0.958231    0.969018    1.001849    1.198493    1.280510    1.312195    1.325138    1.336617
     1.361359    1.590578    1.798930    1.904393    1.914960    2.035719    2.481050    2.535570    2.609987    2.779011
     2.784442    2.858889    2.868959    3.042471    3.103312    3.277038    3.391410    3.418761    3.495641    3.787468
     4.181083    4.256502    4.262888    4.339722    4.493633    4.740107    4.848192    5.013129    5.031089    5.241001
     5.926195    6.232402    6.388625

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   73 emc=    -93.9738288381 demc= 1.4069E-05 wnorm= 3.1518E-03 knorm= 1.0247E-02 apxde= 7.0408E-06    *not conv.*     

               starting mcscf iteration...  74
 !timer:                                 cpu_time=    23.866 walltime=    24.009

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368298 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0815639788     -116.9445258474        0.0000000000        0.0000001250
    2       -93.9370897131     -116.8000515816        0.0000000000        0.0000001250
    3       -93.9028753459     -116.7658372145        0.0000000000        0.0000001250
    4       -93.7838354436     -116.6467973121        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  3.961533032060988E-004
 Total number of micro iterations:    9

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.99994832 pnorm= 0.0000E+00 rznorm= 6.3067E-07 rpnorm= 0.0000E+00 noldr=  9 nnewr=  9 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.714840  -23.219190   -2.317852   -2.067005   -1.850090

 qvv(*) eigenvalues. symmetry block  1
    -0.415903   -0.199288   -0.163020   -0.122798   -0.086121   -0.039358   -0.022772    0.025602    0.043297    0.068289
     0.082226    0.094843    0.143452    0.201010    0.218196    0.278417    0.321790    0.393192    0.417653    0.448057
     0.472885    0.497773    0.527178    0.562340    0.620205    0.648517    0.669708    0.731792    0.801813    0.814162
     0.907257    0.921057    0.958285    0.969067    1.001794    1.198546    1.280598    1.312287    1.325144    1.336617
     1.361420    1.590655    1.798989    1.904530    1.915010    2.035729    2.481001    2.535656    2.610041    2.779076
     2.784524    2.858951    2.869005    3.042351    3.103093    3.276843    3.391274    3.418694    3.495714    3.787259
     4.181043    4.256590    4.262969    4.339787    4.493380    4.739887    4.848215    5.013175    5.030948    5.241088
     5.926280    6.232486    6.388588

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   74 emc=    -93.9738430126 demc= 1.4175E-05 wnorm= 3.1692E-03 knorm= 1.0166E-02 apxde= 7.0818E-06    *not conv.*     

               starting mcscf iteration...  75
 !timer:                                 cpu_time=    24.185 walltime=    24.327

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368296 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0815705717     -116.9445324402        0.0000000000        0.0000001250
    2       -93.9371545596     -116.8001164281        0.0000000000        0.0000001250
    3       -93.9028466470     -116.7658085155        0.0000000000        0.0000001250
    4       -93.7840291546     -116.6469910231        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  3.978055896724654E-004
 Total number of micro iterations:    9

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.99994924 pnorm= 0.0000E+00 rznorm= 6.1002E-07 rpnorm= 0.0000E+00 noldr=  9 nnewr=  9 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.714879  -23.219103   -2.317769   -2.072710   -1.850003

 qvv(*) eigenvalues. symmetry block  1
    -0.415826   -0.199289   -0.163033   -0.122780   -0.086148   -0.039315   -0.022737    0.025662    0.043317    0.068321
     0.082214    0.094861    0.143462    0.201051    0.218139    0.278470    0.321771    0.393166    0.417628    0.448005
     0.472824    0.497767    0.527228    0.562391    0.620265    0.648581    0.669794    0.731830    0.801790    0.814102
     0.907256    0.921013    0.958339    0.969116    1.001739    1.198596    1.280685    1.312377    1.325150    1.336616
     1.361481    1.590731    1.799046    1.904665    1.915058    2.035739    2.480953    2.535741    2.610094    2.779139
     2.784605    2.859011    2.869051    3.042235    3.102873    3.276658    3.391134    3.418626    3.495786    3.787057
     4.181001    4.256679    4.263050    4.339852    4.493129    4.739665    4.848240    5.013214    5.030804    5.241175
     5.926371    6.232570    6.388541

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   75 emc=    -93.9738572594 demc= 1.4247E-05 wnorm= 3.1824E-03 knorm= 1.0076E-02 apxde= 7.1057E-06    *not conv.*     

               starting mcscf iteration...  76
 !timer:                                 cpu_time=    24.503 walltime=    24.646

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368298 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0815767351     -116.9445386037        0.0000000000        0.0000001250
    2       -93.9372195464     -116.8001814150        0.0000000000        0.0000001250
    3       -93.9028183490     -116.7657802176        0.0000000000        0.0000001250
    4       -93.7842151054     -116.6471769740        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  3.989130878186282E-004
 Total number of micro iterations:    9

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.99995023 pnorm= 0.0000E+00 rznorm= 5.8901E-07 rpnorm= 0.0000E+00 noldr=  9 nnewr=  9 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.714920  -23.219017   -2.317686   -2.078282   -1.849916

 qvv(*) eigenvalues. symmetry block  1
    -0.415749   -0.199289   -0.163047   -0.122762   -0.086176   -0.039274   -0.022703    0.025721    0.043337    0.068351
     0.082203    0.094880    0.143471    0.201092    0.218082    0.278523    0.321751    0.393139    0.417604    0.447953
     0.472764    0.497760    0.527276    0.562442    0.620324    0.648644    0.669879    0.731867    0.801768    0.814045
     0.907256    0.920969    0.958392    0.969164    1.001684    1.198645    1.280772    1.312466    1.325155    1.336614
     1.361541    1.590807    1.799101    1.904801    1.915105    2.035747    2.480904    2.535825    2.610147    2.779202
     2.784685    2.859071    2.869096    3.042122    3.102654    3.276485    3.390989    3.418557    3.495857    3.786863
     4.180959    4.256766    4.263129    4.339915    4.492881    4.739444    4.848267    5.013247    5.030658    5.241261
     5.926468    6.232653    6.388484

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   76 emc=    -93.9738715435 demc= 1.4284E-05 wnorm= 3.1913E-03 knorm= 9.9765E-03 apxde= 7.1116E-06    *not conv.*     

               starting mcscf iteration...  77
 !timer:                                 cpu_time=    24.821 walltime=    24.963

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368298 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0815824674     -116.9445443359        0.0000000000        0.0000001250
    2       -93.9372845143     -116.8002463829        0.0000000000        0.0000001250
    3       -93.9027905042     -116.7657523728        0.0000000000        0.0000001250
    4       -93.7843930402     -116.6473549088        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  3.994568552756365E-004
 Total number of micro iterations:    9

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.99995131 pnorm= 0.0000E+00 rznorm= 5.6781E-07 rpnorm= 0.0000E+00 noldr=  9 nnewr=  9 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.714963  -23.218931   -2.317604   -2.083716   -1.849830

 qvv(*) eigenvalues. symmetry block  1
    -0.415674   -0.199289   -0.163060   -0.122744   -0.086204   -0.039232   -0.022669    0.025781    0.043356    0.068381
     0.082192    0.094898    0.143478    0.201133    0.218026    0.278576    0.321732    0.393112    0.417581    0.447901
     0.472705    0.497753    0.527324    0.562493    0.620382    0.648707    0.669962    0.731904    0.801747    0.813988
     0.907256    0.920926    0.958444    0.969212    1.001630    1.198691    1.280857    1.312553    1.325160    1.336610
     1.361599    1.590883    1.799155    1.904935    1.915149    2.035753    2.480854    2.535908    2.610199    2.779265
     2.784765    2.859130    2.869142    3.042014    3.102434    3.276324    3.390839    3.418488    3.495928    3.786676
     4.180916    4.256853    4.263208    4.339978    4.492636    4.739223    4.848296    5.013273    5.030509    5.241346
     5.926572    6.232735    6.388417

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   77 emc=    -93.9738858286 demc= 1.4285E-05 wnorm= 3.1957E-03 knorm= 9.8678E-03 apxde= 7.0992E-06    *not conv.*     

               starting mcscf iteration...  78
 !timer:                                 cpu_time=    25.141 walltime=    25.283

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368298 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0815877690     -116.9445496376        0.0000000000        0.0000001250
    2       -93.9373492997     -116.8003111683        0.0000000000        0.0000001250
    3       -93.9027631638     -116.7657250323        0.0000000000        0.0000001250
    4       -93.7845627380     -116.6475246066        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  3.994219063093869E-004
 Total number of micro iterations:    9

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.99995247 pnorm= 0.0000E+00 rznorm= 5.4662E-07 rpnorm= 0.0000E+00 noldr=  9 nnewr=  9 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.715008  -23.218847   -2.317523   -2.089005   -1.849745

 qvv(*) eigenvalues. symmetry block  1
    -0.415599   -0.199290   -0.163072   -0.122726   -0.086232   -0.039192   -0.022636    0.025840    0.043374    0.068411
     0.082181    0.094916    0.143484    0.201174    0.217972    0.278628    0.321713    0.393085    0.417560    0.447851
     0.472646    0.497745    0.527371    0.562543    0.620440    0.648768    0.670045    0.731941    0.801727    0.813934
     0.907257    0.920883    0.958496    0.969259    1.001576    1.198735    1.280943    1.312638    1.325165    1.336606
     1.361656    1.590959    1.799207    1.905069    1.915191    2.035757    2.480805    2.535989    2.610250    2.779326
     2.784843    2.859187    2.869187    3.041909    3.102215    3.276175    3.390685    3.418418    3.495998    3.786496
     4.180873    4.256938    4.263286    4.340040    4.492395    4.739002    4.848327    5.013293    5.030358    5.241430
     5.926681    6.232816    6.388341

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   78 emc=    -93.9739000775 demc= 1.4249E-05 wnorm= 3.1954E-03 knorm= 9.7498E-03 apxde= 7.0680E-06    *not conv.*     

               starting mcscf iteration...  79
 !timer:                                 cpu_time=    25.460 walltime=    25.603

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368298 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0815926431     -116.9445545117        0.0000000000        0.0000001250
    2       -93.9374137365     -116.8003756050        0.0000000000        0.0000001250
    3       -93.9027363773     -116.7656982459        0.0000000000        0.0000001250
    4       -93.7847240159     -116.6476858845        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  3.987975453359164E-004
 Total number of micro iterations:    9

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.99995370 pnorm= 0.0000E+00 rznorm= 5.2562E-07 rpnorm= 0.0000E+00 noldr=  9 nnewr=  9 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.715054  -23.218764   -2.317443   -2.094145   -1.849661

 qvv(*) eigenvalues. symmetry block  1
    -0.415526   -0.199290   -0.163084   -0.122708   -0.086259   -0.039152   -0.022604    0.025898    0.043392    0.068440
     0.082170    0.094934    0.143488    0.201215    0.217918    0.278680    0.321695    0.393057    0.417539    0.447802
     0.472587    0.497737    0.527417    0.562593    0.620497    0.648829    0.670127    0.731977    0.801709    0.813881
     0.907258    0.920840    0.958548    0.969305    1.001523    1.198778    1.281027    1.312722    1.325169    1.336600
     1.361712    1.591034    1.799257    1.905202    1.915231    2.035760    2.480756    2.536069    2.610301    2.779386
     2.784920    2.859244    2.869231    3.041808    3.101998    3.276038    3.390526    3.418348    3.496067    3.786325
     4.180830    4.257023    4.263363    4.340101    4.492157    4.738783    4.848360    5.013305    5.030206    5.241513
     5.926796    6.232896    6.388256

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   79 emc=    -93.9739142523 demc= 1.4175E-05 wnorm= 3.1904E-03 knorm= 9.6229E-03 apxde= 7.0178E-06    *not conv.*     

               starting mcscf iteration...  80
 !timer:                                 cpu_time=    25.780 walltime=    25.922

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368298 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0815970953     -116.9445589639        0.0000000000        0.0000001250
    2       -93.9374776577     -116.8004395262        0.0000000000        0.0000001250
    3       -93.9027101923     -116.7656720609        0.0000000000        0.0000001250
    4       -93.7848767312     -116.6478385997        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  3.975776289879370E-004
 Total number of micro iterations:    9

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.99995500 pnorm= 0.0000E+00 rznorm= 5.0498E-07 rpnorm= 0.0000E+00 noldr=  9 nnewr=  9 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.715101  -23.218683   -2.317365   -2.099132   -1.849579

 qvv(*) eigenvalues. symmetry block  1
    -0.415455   -0.199290   -0.163096   -0.122690   -0.086286   -0.039114   -0.022574    0.025956    0.043410    0.068468
     0.082159    0.094952    0.143492    0.201256    0.217866    0.278731    0.321676    0.393029    0.417521    0.447753
     0.472530    0.497729    0.527463    0.562642    0.620553    0.648888    0.670206    0.732012    0.801692    0.813831
     0.907260    0.920799    0.958598    0.969351    1.001470    1.198818    1.281110    1.312803    1.325173    1.336592
     1.361766    1.591108    1.799305    1.905333    1.915270    2.035762    2.480707    2.536148    2.610350    2.779446
     2.784995    2.859299    2.869276    3.041712    3.101781    3.275913    3.390363    3.418279    3.496135    3.786162
     4.180785    4.257106    4.263438    4.340162    4.491923    4.738565    4.848395    5.013312    5.030053    5.241594
     5.926917    6.232974    6.388162

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   80 emc=    -93.9739283151 demc= 1.4063E-05 wnorm= 3.1806E-03 knorm= 9.4872E-03 apxde= 6.9489E-06    *not conv.*     

               starting mcscf iteration...  81
 !timer:                                 cpu_time=    26.098 walltime=    26.240

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368298 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0816011335     -116.9445630020        0.0000000000        0.0000001250
    2       -93.9375408975     -116.8005027661        0.0000000000        0.0000001250
    3       -93.9026846535     -116.7656465221        0.0000000000        0.0000001250
    4       -93.7850207832     -116.6479826517        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  3.957607468349785E-004
 Total number of micro iterations:    9

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.99995635 pnorm= 0.0000E+00 rznorm= 4.8487E-07 rpnorm= 0.0000E+00 noldr=  9 nnewr=  9 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.715150  -23.218603   -2.317289   -2.103962   -1.849498

 qvv(*) eigenvalues. symmetry block  1
    -0.415384   -0.199290   -0.163107   -0.122672   -0.086313   -0.039076   -0.022544    0.026013    0.043426    0.068495
     0.082149    0.094970    0.143494    0.201297    0.217814    0.278781    0.321658    0.393001    0.417504    0.447706
     0.472474    0.497720    0.527507    0.562690    0.620608    0.648946    0.670285    0.732047    0.801675    0.813782
     0.907263    0.920758    0.958647    0.969396    1.001419    1.198855    1.281191    1.312882    1.325176    1.336584
     1.361819    1.591181    1.799351    1.905463    1.915306    2.035762    2.480657    2.536224    2.610398    2.779504
     2.785069    2.859353    2.869320    3.041619    3.101567    3.275800    3.390197    3.418209    3.496202    3.786008
     4.180741    4.257188    4.263513    4.340221    4.491694    4.738348    4.848432    5.013311    5.029898    5.241675
     5.927042    6.233051    6.388058

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   81 emc=    -93.9739422282 demc= 1.3913E-05 wnorm= 3.1661E-03 knorm= 9.3430E-03 apxde= 6.8614E-06    *not conv.*     

               starting mcscf iteration...  82
 !timer:                                 cpu_time=    26.416 walltime=    26.558

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368297 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0816047678     -116.9445666364        0.0000000000        0.0000001250
    2       -93.9376032935     -116.8005651621        0.0000000000        0.0000001250
    3       -93.9026598026     -116.7656216712        0.0000000000        0.0000001250
    4       -93.7851561147     -116.6481179833        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  3.933503128017584E-004
 Total number of micro iterations:    8

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.99995776 pnorm= 0.0000E+00 rznorm= 1.2284E-06 rpnorm= 0.0000E+00 noldr=  8 nnewr=  8 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.715201  -23.218525   -2.317214   -2.108633   -1.849419

 qvv(*) eigenvalues. symmetry block  1
    -0.415316   -0.199291   -0.163118   -0.122654   -0.086339   -0.039039   -0.022515    0.026070    0.043442    0.068522
     0.082139    0.094987    0.143495    0.201338    0.217764    0.278830    0.321641    0.392974    0.417488    0.447660
     0.472418    0.497710    0.527550    0.562738    0.620662    0.649003    0.670362    0.732081    0.801661    0.813735
     0.907266    0.920717    0.958695    0.969440    1.001368    1.198891    1.281272    1.312958    1.325179    1.336575
     1.361871    1.591254    1.799395    1.905590    1.915340    2.035760    2.480609    2.536299    2.610446    2.779561
     2.785141    2.859406    2.869363    3.041532    3.101354    3.275700    3.390027    3.418140    3.496267    3.785862
     4.180697    4.257268    4.263585    4.340278    4.491470    4.738134    4.848470    5.013304    5.029742    5.241754
     5.927173    6.233126    6.387947

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   82 emc=    -93.9739559546 demc= 1.3726E-05 wnorm= 3.1468E-03 knorm= 9.1908E-03 apxde= 6.7560E-06    *not conv.*     

               starting mcscf iteration...  83
 !timer:                                 cpu_time=    26.734 walltime=    26.876

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368298 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0816080146     -116.9445698832        0.0000000000        0.0000001250
    2       -93.9376646810     -116.8006265496        0.0000000000        0.0000001250
    3       -93.9026356807     -116.7655975492        0.0000000000        0.0000001250
    4       -93.7852827087     -116.6482445773        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  3.903543279075610E-004
 Total number of micro iterations:    8

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.99995922 pnorm= 0.0000E+00 rznorm= 1.2073E-06 rpnorm= 0.0000E+00 noldr=  8 nnewr=  8 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.715253  -23.218449   -2.317140   -2.113142   -1.849342

 qvv(*) eigenvalues. symmetry block  1
    -0.415249   -0.199291   -0.163128   -0.122636   -0.086366   -0.039003   -0.022487    0.026126    0.043458    0.068548
     0.082129    0.095003    0.143495    0.201378    0.217715    0.278879    0.321623    0.392946    0.417473    0.447614
     0.472364    0.497701    0.527592    0.562785    0.620715    0.649059    0.670437    0.732115    0.801647    0.813690
     0.907270    0.920678    0.958743    0.969483    1.001318    1.198924    1.281351    1.313032    1.325182    1.336564
     1.361921    1.591325    1.799437    1.905716    1.915372    2.035758    2.480560    2.536373    2.610492    2.779616
     2.785211    2.859457    2.869405    3.041449    3.101145    3.275611    3.389854    3.418072    3.496331    3.785725
     4.180652    4.257346    4.263656    4.340335    4.491251    4.737923    4.848510    5.013291    5.029586    5.241831
     5.927307    6.233199    6.387828

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   83 emc=    -93.9739694588 demc= 1.3504E-05 wnorm= 3.1228E-03 knorm= 9.0307E-03 apxde= 6.6333E-06    *not conv.*     

               starting mcscf iteration...  84
 !timer:                                 cpu_time=    27.052 walltime=    27.194

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368298 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0816108801     -116.9445727486        0.0000000000        0.0000001250
    2       -93.9377249240     -116.8006867925        0.0000000000        0.0000001250
    3       -93.9026123154     -116.7655741839        0.0000000000        0.0000001250
    4       -93.7854006032     -116.6483624717        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  3.867861150105219E-004
 Total number of micro iterations:    8

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.99996072 pnorm= 0.0000E+00 rznorm= 1.1799E-06 rpnorm= 0.0000E+00 noldr=  8 nnewr=  8 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.715306  -23.218374   -2.317069   -2.117488   -1.849266

 qvv(*) eigenvalues. symmetry block  1
    -0.415184   -0.199291   -0.163138   -0.122618   -0.086391   -0.038968   -0.022460    0.026181    0.043472    0.068573
     0.082119    0.095020    0.143494    0.201418    0.217668    0.278926    0.321607    0.392918    0.417461    0.447571
     0.472311    0.497690    0.527633    0.562831    0.620766    0.649113    0.670510    0.732147    0.801635    0.813647
     0.907274    0.920639    0.958789    0.969524    1.001269    1.198955    1.281428    1.313103    1.325184    1.336553
     1.361969    1.591395    1.799477    1.905839    1.915402    2.035753    2.480512    2.536444    2.610537    2.779671
     2.785280    2.859507    2.869447    3.041370    3.100938    3.275535    3.389679    3.418004    3.496394    3.785597
     4.180608    4.257423    4.263725    4.340390    4.491037    4.737714    4.848552    5.013272    5.029430    5.241906
     5.927446    6.233271    6.387701

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   84 emc=    -93.9739827065 demc= 1.3248E-05 wnorm= 3.0943E-03 knorm= 8.8636E-03 apxde= 6.4943E-06    *not conv.*     

               starting mcscf iteration...  85
 !timer:                                 cpu_time=    27.371 walltime=    27.513

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368298 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0816133842     -116.9445752528        0.0000000000        0.0000001250
    2       -93.9377838731     -116.8007457416        0.0000000000        0.0000001250
    3       -93.9025897394     -116.7655516080        0.0000000000        0.0000001250
    4       -93.7855098684     -116.6484717370        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  3.826631586069754E-004
 Total number of micro iterations:    8

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.99996224 pnorm= 0.0000E+00 rznorm= 1.1517E-06 rpnorm= 0.0000E+00 noldr=  8 nnewr=  8 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.715360  -23.218302   -2.317000   -2.121669   -1.849193

 qvv(*) eigenvalues. symmetry block  1
    -0.415121   -0.199291   -0.163148   -0.122600   -0.086416   -0.038935   -0.022435    0.026235    0.043486    0.068598
     0.082109    0.095036    0.143491    0.201457    0.217622    0.278972    0.321590    0.392890    0.417449    0.447528
     0.472259    0.497680    0.527673    0.562875    0.620816    0.649166    0.670581    0.732179    0.801624    0.813607
     0.907279    0.920601    0.958833    0.969565    1.001221    1.198984    1.281503    1.313172    1.325186    1.336541
     1.362015    1.591464    1.799515    1.905960    1.915430    2.035748    2.480464    2.536513    2.610581    2.779723
     2.785347    2.859555    2.869488    3.041296    3.100734    3.275471    3.389501    3.417937    3.496455    3.785478
     4.180564    4.257497    4.263792    4.340444    4.490830    4.737509    4.848595    5.013247    5.029274    5.241980
     5.927588    6.233341    6.387567

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   85 emc=    -93.9739956656 demc= 1.2959E-05 wnorm= 3.0613E-03 knorm= 8.6897E-03 apxde= 6.3401E-06    *not conv.*     

               starting mcscf iteration...  86
 !timer:                                 cpu_time=    27.691 walltime=    27.831

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368297 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0816155443     -116.9445774128        0.0000000000        0.0000001250
    2       -93.9378413968     -116.8008032654        0.0000000000        0.0000001250
    3       -93.9025679780     -116.7655298465        0.0000000000        0.0000001250
    4       -93.7856106200     -116.6485724886        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  3.780071065418603E-004
 Total number of micro iterations:    8

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.99996379 pnorm= 0.0000E+00 rznorm= 1.1229E-06 rpnorm= 0.0000E+00 noldr=  8 nnewr=  8 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.715414  -23.218232   -2.316932   -2.125685   -1.849121

 qvv(*) eigenvalues. symmetry block  1
    -0.415059   -0.199291   -0.163157   -0.122583   -0.086441   -0.038902   -0.022410    0.026288    0.043500    0.068621
     0.082100    0.095051    0.143488    0.201496    0.217577    0.279017    0.321574    0.392862    0.417439    0.447487
     0.472208    0.497669    0.527712    0.562919    0.620865    0.649217    0.670650    0.732211    0.801614    0.813568
     0.907284    0.920565    0.958877    0.969605    1.001174    1.199010    1.281577    1.313238    1.325188    1.336528
     1.362060    1.591531    1.799550    1.906078    1.915456    2.035741    2.480417    2.536579    2.610623    2.779774
     2.785411    2.859602    2.869528    3.041227    3.100535    3.275419    3.389321    3.417871    3.496515    3.785368
     4.180519    4.257570    4.263858    4.340496    4.490628    4.737307    4.848639    5.013217    5.029119    5.242051
     5.927734    6.233408    6.387427

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   86 emc=    -93.9740083063 demc= 1.2641E-05 wnorm= 3.0241E-03 knorm= 8.5097E-03 apxde= 6.1721E-06    *not conv.*     

               starting mcscf iteration...  87
 !timer:                                 cpu_time=    28.008 walltime=    28.149

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368297 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0816173788     -116.9445792474        0.0000000000        0.0000001250
    2       -93.9378973751     -116.8008592437        0.0000000000        0.0000001250
    3       -93.9025470512     -116.7655089197        0.0000000000        0.0000001250
    4       -93.7857030135     -116.6486648821        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  3.728436232545713E-004
 Total number of micro iterations:    8

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.99996535 pnorm= 0.0000E+00 rznorm= 1.0938E-06 rpnorm= 0.0000E+00 noldr=  8 nnewr=  8 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.715470  -23.218163   -2.316867   -2.129537   -1.849052

 qvv(*) eigenvalues. symmetry block  1
    -0.415000   -0.199291   -0.163165   -0.122566   -0.086465   -0.038871   -0.022387    0.026340    0.043513    0.068644
     0.082091    0.095066    0.143483    0.201534    0.217534    0.279061    0.321558    0.392835    0.417430    0.447447
     0.472159    0.497658    0.527749    0.562962    0.620912    0.649266    0.670717    0.732241    0.801605    0.813531
     0.907290    0.920529    0.958919    0.969643    1.001129    1.199035    1.281649    1.313301    1.325189    1.336514
     1.362103    1.591597    1.799583    1.906193    1.915481    2.035733    2.480371    2.536644    2.610664    2.779824
     2.785474    2.859647    2.869567    3.041162    3.100339    3.275378    3.389140    3.417807    3.496572    3.785266
     4.180476    4.257641    4.263921    4.340546    4.490432    4.737109    4.848685    5.013181    5.028965    5.242121
     5.927882    6.233474    6.387281

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   87 emc=    -93.9740206017 demc= 1.2295E-05 wnorm= 2.9827E-03 knorm= 8.3241E-03 apxde= 5.9915E-06    *not conv.*     

               starting mcscf iteration...  88
 !timer:                                 cpu_time=    28.326 walltime=    28.467

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368298 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0816189074     -116.9445807759        0.0000000000        0.0000001250
    2       -93.9379517008     -116.8009135693        0.0000000000        0.0000001250
    3       -93.9025269742     -116.7654888428        0.0000000000        0.0000001250
    4       -93.7857872416     -116.6487491101        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  3.672018594458051E-004
 Total number of micro iterations:    8

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.99996692 pnorm= 0.0000E+00 rznorm= 1.0645E-06 rpnorm= 0.0000E+00 noldr=  8 nnewr=  8 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.715526  -23.218097   -2.316804   -2.133226   -1.848985

 qvv(*) eigenvalues. symmetry block  1
    -0.414943   -0.199291   -0.163173   -0.122549   -0.086489   -0.038841   -0.022365    0.026391    0.043525    0.068665
     0.082083    0.095081    0.143478    0.201572    0.217493    0.279104    0.321543    0.392808    0.417423    0.447408
     0.472111    0.497646    0.527785    0.563003    0.620958    0.649314    0.670782    0.732270    0.801598    0.813497
     0.907296    0.920494    0.958960    0.969680    1.001084    1.199057    1.281718    1.313361    1.325191    1.336499
     1.362144    1.591662    1.799614    1.906304    1.915504    2.035724    2.480325    2.536706    2.610703    2.779872
     2.785535    2.859690    2.869605    3.041101    3.100147    3.275347    3.388957    3.417744    3.496628    3.785172
     4.180432    4.257709    4.263983    4.340595    4.490243    4.736915    4.848731    5.013141    5.028811    5.242189
     5.928032    6.233537    6.387129

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   88 emc=    -93.9740325275 demc= 1.1926E-05 wnorm= 2.9376E-03 knorm= 8.1337E-03 apxde= 5.8000E-06    *not conv.*     

               starting mcscf iteration...  89
 !timer:                                 cpu_time=    28.644 walltime=    28.785

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368298 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0816201500     -116.9445820186        0.0000000000        0.0000001250
    2       -93.9380042802     -116.8009661488        0.0000000000        0.0000001250
    3       -93.9025077575     -116.7654696261        0.0000000000        0.0000001250
    4       -93.7858635313     -116.6488253999        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  3.611139776481748E-004
 Total number of micro iterations:    8

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.99996849 pnorm= 0.0000E+00 rznorm= 1.0350E-06 rpnorm= 0.0000E+00 noldr=  8 nnewr=  8 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.715582  -23.218034   -2.316743   -2.136753   -1.848920

 qvv(*) eigenvalues. symmetry block  1
    -0.414888   -0.199291   -0.163181   -0.122532   -0.086512   -0.038811   -0.022343    0.026440    0.043536    0.068686
     0.082074    0.095095    0.143472    0.201609    0.217452    0.279146    0.321529    0.392781    0.417417    0.447371
     0.472065    0.497634    0.527819    0.563044    0.621002    0.649360    0.670845    0.732299    0.801592    0.813464
     0.907303    0.920461    0.958999    0.969716    1.001041    1.199077    1.281786    1.313419    1.325192    1.336483
     1.362184    1.591725    1.799643    1.906412    1.915526    2.035714    2.480280    2.536767    2.610741    2.779919
     2.785593    2.859732    2.869643    3.041045    3.099960    3.275327    3.388774    3.417682    3.496683    3.785087
     4.180390    4.257775    4.264042    4.340643    4.490061    4.736726    4.848779    5.013096    5.028659    5.242254
     5.928185    6.233599    6.386973

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   89 emc=    -93.9740440626 demc= 1.1535E-05 wnorm= 2.8889E-03 knorm= 7.9390E-03 apxde= 5.5992E-06    *not conv.*     

               starting mcscf iteration...  90
 !timer:                                 cpu_time=    28.963 walltime=    29.104

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368298 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0816211273     -116.9445829959        0.0000000000        0.0000001250
    2       -93.9380550341     -116.8010169026        0.0000000000        0.0000001250
    3       -93.9024894064     -116.7654512750        0.0000000000        0.0000001250
    4       -93.7859321408     -116.6488940094        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  3.546146405619398E-004
 Total number of micro iterations:    8

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.99997004 pnorm= 0.0000E+00 rznorm= 1.0053E-06 rpnorm= 0.0000E+00 noldr=  8 nnewr=  8 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.715639  -23.217972   -2.316684   -2.140120   -1.848857

 qvv(*) eigenvalues. symmetry block  1
    -0.414835   -0.199291   -0.163188   -0.122516   -0.086535   -0.038784   -0.022323    0.026488    0.043547    0.068706
     0.082066    0.095109    0.143465    0.201645    0.217414    0.279186    0.321515    0.392755    0.417412    0.447335
     0.472020    0.497622    0.527853    0.563083    0.621045    0.649404    0.670905    0.732326    0.801587    0.813433
     0.907310    0.920428    0.959037    0.969751    1.000999    1.199095    1.281852    1.313473    1.325193    1.336468
     1.362221    1.591786    1.799670    1.906517    1.915546    2.035703    2.480236    2.536824    2.610778    2.779963
     2.785649    2.859772    2.869679    3.040993    3.099777    3.275317    3.388591    3.417621    3.496735    3.785010
     4.180347    4.257839    4.264100    4.340688    4.489885    4.736541    4.848827    5.013047    5.028509    5.242317
     5.928338    6.233658    6.386813

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   90 emc=    -93.9740551893 demc= 1.1127E-05 wnorm= 2.8369E-03 knorm= 7.7407E-03 apxde= 5.3906E-06    *not conv.*     

               starting mcscf iteration...  91
 !timer:                                 cpu_time=    29.281 walltime=    29.421

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368297 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0816218600     -116.9445837285        0.0000000000        0.0000001250
    2       -93.9381038975     -116.8010657661        0.0000000000        0.0000001250
    3       -93.9024719215     -116.7654337900        0.0000000000        0.0000001250
    4       -93.7859933554     -116.6489552239        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  3.477404753433874E-004
 Total number of micro iterations:    9

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.99997158 pnorm= 0.0000E+00 rznorm= 3.3247E-07 rpnorm= 0.0000E+00 noldr=  9 nnewr=  9 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.715695  -23.217913   -2.316627   -2.143331   -1.848797

 qvv(*) eigenvalues. symmetry block  1
    -0.414784   -0.199291   -0.163195   -0.122500   -0.086557   -0.038757   -0.022304    0.026535    0.043557    0.068725
     0.082059    0.095122    0.143457    0.201680    0.217376    0.279225    0.321501    0.392729    0.417408    0.447301
     0.471976    0.497610    0.527885    0.563121    0.621086    0.649447    0.670963    0.732353    0.801583    0.813404
     0.907317    0.920397    0.959074    0.969784    1.000959    1.199111    1.281915    1.313525    1.325194    1.336451
     1.362257    1.591845    1.799695    1.906618    1.915565    2.035691    2.480192    2.536880    2.610813    2.780006
     2.785703    2.859810    2.869714    3.040945    3.099600    3.275316    3.388407    3.417562    3.496786    3.784941
     4.180306    4.257901    4.264155    4.340732    4.489715    4.736361    4.848876    5.012994    5.028361    5.242378
     5.928493    6.233715    6.386649

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   91 emc=    -93.9740658930 demc= 1.0704E-05 wnorm= 2.7819E-03 knorm= 7.5393E-03 apxde= 5.1760E-06    *not conv.*     

               starting mcscf iteration...  92
 !timer:                                 cpu_time=    29.602 walltime=    29.743

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368298 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0816223645     -116.9445842331        0.0000000000        0.0000001250
    2       -93.9381508256     -116.8011126942        0.0000000000        0.0000001250
    3       -93.9024552978     -116.7654171663        0.0000000000        0.0000001250
    4       -93.7860474897     -116.6490093583        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  3.405297931276908E-004
 Total number of micro iterations:    9

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.99997309 pnorm= 0.0000E+00 rznorm= 3.2169E-07 rpnorm= 0.0000E+00 noldr=  9 nnewr=  9 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.715752  -23.217856   -2.316573   -2.146387   -1.848739

 qvv(*) eigenvalues. symmetry block  1
    -0.414735   -0.199291   -0.163201   -0.122484   -0.086578   -0.038731   -0.022287    0.026581    0.043567    0.068743
     0.082052    0.095134    0.143449    0.201715    0.217341    0.279262    0.321488    0.392704    0.417406    0.447268
     0.471935    0.497598    0.527915    0.563157    0.621126    0.649488    0.671019    0.732379    0.801580    0.813377
     0.907324    0.920366    0.959109    0.969816    1.000920    1.199125    1.281977    1.313574    1.325195    1.336434
     1.362291    1.591903    1.799718    1.906715    1.915583    2.035679    2.480150    2.536933    2.610847    2.780048
     2.785755    2.859847    2.869748    3.040900    3.099427    3.275325    3.388225    3.417505    3.496835    3.784879
     4.180265    4.257960    4.264208    4.340775    4.489552    4.736186    4.848925    5.012938    5.028215    5.242437
     5.928648    6.233769    6.386483

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   92 emc=    -93.9740761626 demc= 1.0270E-05 wnorm= 2.7242E-03 knorm= 7.3359E-03 apxde= 4.9570E-06    *not conv.*     

               starting mcscf iteration...  93
 !timer:                                 cpu_time=    29.920 walltime=    30.061

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368297 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0816226693     -116.9445845379        0.0000000000        0.0000001250
    2       -93.9381957722     -116.8011576408        0.0000000000        0.0000001250
    3       -93.9024395293     -116.7654013979        0.0000000000        0.0000001250
    4       -93.7860948596     -116.6490567281        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  3.330210299227858E-004
 Total number of micro iterations:    9

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.99997458 pnorm= 0.0000E+00 rznorm= 3.1215E-07 rpnorm= 0.0000E+00 noldr=  9 nnewr=  9 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.715808  -23.217801   -2.316520   -2.149294   -1.848683

 qvv(*) eigenvalues. symmetry block  1
    -0.414688   -0.199291   -0.163207   -0.122469   -0.086598   -0.038707   -0.022270    0.026625    0.043575    0.068760
     0.082044    0.095146    0.143440    0.201748    0.217306    0.279298    0.321476    0.392679    0.417404    0.447236
     0.471894    0.497586    0.527944    0.563193    0.621164    0.649527    0.671073    0.732403    0.801578    0.813352
     0.907332    0.920337    0.959142    0.969847    1.000882    1.199138    1.282036    1.313620    1.325196    1.336417
     1.362324    1.591959    1.799739    1.906808    1.915600    2.035666    2.480109    2.536984    2.610880    2.780087
     2.785805    2.859882    2.869781    3.040860    3.099260    3.275341    3.388043    3.417449    3.496882    3.784824
     4.180225    4.258018    4.264259    4.340815    4.489396    4.736015    4.848974    5.012878    5.028072    5.242494
     5.928803    6.233822    6.386314

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   93 emc=    -93.9740859903 demc= 9.8277E-06 wnorm= 2.6642E-03 knorm= 7.1308E-03 apxde= 4.7352E-06    *not conv.*     

               starting mcscf iteration...  94
 !timer:                                 cpu_time=    30.239 walltime=    30.395

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368298 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0816227901     -116.9445846587        0.0000000000        0.0000001250
    2       -93.9382387210     -116.8012005895        0.0000000000        0.0000001250
    3       -93.9024246028     -116.7653864714        0.0000000000        0.0000001250
    4       -93.7861358143     -116.6490976829        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  3.252536816137372E-004
 Total number of micro iterations:    9

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.99997602 pnorm= 0.0000E+00 rznorm= 3.0332E-07 rpnorm= 0.0000E+00 noldr=  9 nnewr=  9 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.715864  -23.217749   -2.316470   -2.152055   -1.848629

 qvv(*) eigenvalues. symmetry block  1
    -0.414644   -0.199290   -0.163213   -0.122454   -0.086618   -0.038684   -0.022254    0.026668    0.043584    0.068776
     0.082038    0.095158    0.143430    0.201781    0.217274    0.279333    0.321464    0.392654    0.417404    0.447206
     0.471856    0.497573    0.527972    0.563227    0.621200    0.649564    0.671124    0.732427    0.801578    0.813328
     0.907340    0.920309    0.959174    0.969877    1.000845    1.199149    1.282093    1.313664    1.325197    1.336400
     1.362354    1.592013    1.799758    1.906897    1.915616    2.035652    2.480068    2.537032    2.610910    2.780125
     2.785852    2.859916    2.869812    3.040823    3.099097    3.275365    3.387863    3.417395    3.496927    3.784777
     4.180186    4.258073    4.264308    4.340854    4.489246    4.735850    4.849024    5.012816    5.027931    5.242548
     5.928958    6.233872    6.386143

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   94 emc=    -93.9740953713 demc= 9.3810E-06 wnorm= 2.6020E-03 knorm= 6.9246E-03 apxde= 4.5121E-06    *not conv.*     

               starting mcscf iteration...  95
 !timer:                                 cpu_time=    30.558 walltime=    30.713

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368298 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0816227460     -116.9445846146        0.0000000000        0.0000001250
    2       -93.9382796637     -116.8012415322        0.0000000000        0.0000001250
    3       -93.9024105022     -116.7653723708        0.0000000000        0.0000001250
    4       -93.7861707086     -116.6491325771        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  3.172669696165985E-004
 Total number of micro iterations:    9

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.99997743 pnorm= 0.0000E+00 rznorm= 2.9512E-07 rpnorm= 0.0000E+00 noldr=  9 nnewr=  9 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.715919  -23.217699   -2.316422   -2.154673   -1.848578

 qvv(*) eigenvalues. symmetry block  1
    -0.414601   -0.199290   -0.163218   -0.122439   -0.086638   -0.038662   -0.022239    0.026710    0.043592    0.068792
     0.082031    0.095169    0.143420    0.201813    0.217242    0.279366    0.321452    0.392630    0.417404    0.447177
     0.471819    0.497561    0.527999    0.563259    0.621235    0.649600    0.671174    0.732450    0.801578    0.813307
     0.907349    0.920282    0.959205    0.969905    1.000810    1.199158    1.282148    1.313704    1.325197    1.336382
     1.362384    1.592065    1.799776    1.906982    1.915631    2.035638    2.480029    2.537079    2.610940    2.780162
     2.785898    2.859948    2.869843    3.040790    3.098941    3.275397    3.387684    3.417343    3.496970    3.784735
     4.180148    4.258125    4.264355    4.340892    4.489103    4.735690    4.849074    5.012752    5.027793    5.242600
     5.929113    6.233920    6.385971

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   95 emc=    -93.9741043040 demc= 8.9327E-06 wnorm= 2.5381E-03 knorm= 6.7181E-03 apxde= 4.2892E-06    *not conv.*     

               starting mcscf iteration...  96
 !timer:                                 cpu_time=    30.876 walltime=    31.031

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368298 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0816225554     -116.9445844240        0.0000000000        0.0000001250
    2       -93.9383186047     -116.8012804732        0.0000000000        0.0000001250
    3       -93.9023972085     -116.7653590771        0.0000000000        0.0000001250
    4       -93.7861999042     -116.6491617728        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  3.090995070151140E-004
 Total number of micro iterations:    9

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.99997880 pnorm= 0.0000E+00 rznorm= 2.8751E-07 rpnorm= 0.0000E+00 noldr=  9 nnewr=  9 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.715974  -23.217651   -2.316377   -2.157155   -1.848529

 qvv(*) eigenvalues. symmetry block  1
    -0.414560   -0.199289   -0.163223   -0.122425   -0.086656   -0.038641   -0.022225    0.026750    0.043599    0.068806
     0.082025    0.095179    0.143410    0.201844    0.217212    0.279399    0.321442    0.392607    0.417406    0.447149
     0.471783    0.497548    0.528024    0.563291    0.621269    0.649635    0.671221    0.732472    0.801579    0.813287
     0.907357    0.920257    0.959234    0.969932    1.000776    1.199165    1.282200    1.313743    1.325198    1.336365
     1.362411    1.592115    1.799792    1.907064    1.915646    2.035623    2.479991    2.537123    2.610968    2.780197
     2.785941    2.859978    2.869872    3.040760    3.098789    3.275434    3.387508    3.417292    3.497011    3.784700
     4.180110    4.258176    4.264399    4.340928    4.488966    4.735536    4.849123    5.012685    5.027658    5.242650
     5.929266    6.233966    6.385799

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   96 emc=    -93.9741127895 demc= 8.4856E-06 wnorm= 2.4728E-03 knorm= 6.5118E-03 apxde= 4.0679E-06    *not conv.*     

               starting mcscf iteration...  97
 !timer:                                 cpu_time=    31.194 walltime=    31.349

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368296 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0816222356     -116.9445841042        0.0000000000        0.0000001250
    2       -93.9383555601     -116.8013174286        0.0000000000        0.0000001250
    3       -93.9023846997     -116.7653465682        0.0000000000        0.0000001250
    4       -93.7862237670     -116.6491856355        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  3.007889408313569E-004
 Total number of micro iterations:    9

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.99998012 pnorm= 0.0000E+00 rznorm= 2.8045E-07 rpnorm= 0.0000E+00 noldr=  9 nnewr=  9 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.716028  -23.217605   -2.316333   -2.159504   -1.848482

 qvv(*) eigenvalues. symmetry block  1
    -0.414522   -0.199289   -0.163227   -0.122411   -0.086674   -0.038621   -0.022212    0.026789    0.043606    0.068820
     0.082019    0.095189    0.143399    0.201874    0.217184    0.279429    0.321431    0.392584    0.417408    0.447123
     0.471749    0.497535    0.528048    0.563321    0.621301    0.649667    0.671265    0.732493    0.801580    0.813268
     0.907365    0.920232    0.959262    0.969958    1.000744    1.199172    1.282251    1.313778    1.325199    1.336347
     1.362437    1.592164    1.799806    1.907141    1.915660    2.035608    2.479954    2.537165    2.610995    2.780230
     2.785983    2.860007    2.869901    3.040733    3.098643    3.275478    3.387334    3.417244    3.497051    3.784671
     4.180074    4.258224    4.264442    4.340962    4.488835    4.735386    4.849173    5.012617    5.027526    5.242698
     5.929418    6.234010    6.385626

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   97 emc=    -93.9741208318 demc= 8.0423E-06 wnorm= 2.4063E-03 knorm= 6.3062E-03 apxde= 3.8493E-06    *not conv.*     

               starting mcscf iteration...  98
 !timer:                                 cpu_time=    31.512 walltime=    31.667

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368297 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0816218031     -116.9445836717        0.0000000000        0.0000001250
    2       -93.9383905565     -116.8013524250        0.0000000000        0.0000001250
    3       -93.9023729513     -116.7653348198        0.0000000000        0.0000001250
    4       -93.7862426629     -116.6492045314        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  2.923716030266469E-004
 Total number of micro iterations:    9

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.99998138 pnorm= 0.0000E+00 rznorm= 2.7390E-07 rpnorm= 0.0000E+00 noldr=  9 nnewr=  9 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.716081  -23.217561   -2.316292   -2.161725   -1.848438

 qvv(*) eigenvalues. symmetry block  1
    -0.414485   -0.199288   -0.163231   -0.122398   -0.086692   -0.038602   -0.022200    0.026826    0.043612    0.068833
     0.082014    0.095199    0.143388    0.201903    0.217156    0.279459    0.321421    0.392562    0.417410    0.447097
     0.471716    0.497523    0.528071    0.563350    0.621331    0.649698    0.671308    0.732513    0.801583    0.813251
     0.907374    0.920208    0.959289    0.969982    1.000712    1.199177    1.282299    1.313812    1.325200    1.336329
     1.362462    1.592210    1.799818    1.907214    1.915674    2.035593    2.479918    2.537205    2.611020    2.780262
     2.786022    2.860034    2.869928    3.040708    3.098502    3.275527    3.387162    3.417196    3.497089    3.784647
     4.180039    4.258270    4.264483    4.340994    4.488710    4.735242    4.849222    5.012547    5.027398    5.242744
     5.929568    6.234052    6.385454

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   98 emc=    -93.9741284370 demc= 7.6052E-06 wnorm= 2.3390E-03 knorm= 6.1020E-03 apxde= 3.6345E-06    *not conv.*     

               starting mcscf iteration...  99
 !timer:                                 cpu_time=    31.831 walltime=    31.986

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368297 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0816212733     -116.9445831418        0.0000000000        0.0000001250
    2       -93.9384236301     -116.8013854987        0.0000000000        0.0000001250
    3       -93.9023619369     -116.7653238055        0.0000000000        0.0000001250
    4       -93.7862569551     -116.6492188237        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  2.838822152856613E-004
 Total number of micro iterations:    9

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.99998260 pnorm= 0.0000E+00 rznorm= 2.6780E-07 rpnorm= 0.0000E+00 noldr=  9 nnewr=  9 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.716133  -23.217520   -2.316252   -2.163823   -1.848395

 qvv(*) eigenvalues. symmetry block  1
    -0.414450   -0.199287   -0.163235   -0.122385   -0.086708   -0.038585   -0.022189    0.026862    0.043618    0.068845
     0.082009    0.095208    0.143377    0.201931    0.217130    0.279487    0.321412    0.392541    0.417414    0.447074
     0.471685    0.497510    0.528093    0.563378    0.621360    0.649728    0.671349    0.732532    0.801586    0.813235
     0.907383    0.920186    0.959314    0.970005    1.000682    1.199180    1.282345    1.313843    1.325200    1.336311
     1.362485    1.592255    1.799830    1.907283    1.915687    2.035578    2.479883    2.537243    2.611045    2.780292
     2.786059    2.860060    2.869954    3.040687    3.098367    3.275580    3.386994    3.417151    3.497125    3.784627
     4.180005    4.258314    4.264522    4.341025    4.488592    4.735103    4.849271    5.012477    5.027273    5.242788
     5.929716    6.234091    6.385283

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   99 emc=    -93.9741356134 demc= 7.1765E-06 wnorm= 2.2711E-03 knorm= 5.8996E-03 apxde= 3.4246E-06    *not conv.*     

               starting mcscf iteration... 100
 !timer:                                 cpu_time=    32.149 walltime=    32.304

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    82, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 393184042

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 392998378
 address segment size,           sizesg = 392846892
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     368297 transformed 1/r12    array elements were written in      68 records.


 mosort: allocated sort2 space, avc2is=   393045895 available sort2 space, avcisx=   393046147

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.0816206603     -116.9445825289        0.0000000000        0.0000001250
    2       -93.9384548255     -116.8014166940        0.0000000000        0.0000001250
    3       -93.9023516285     -116.7653134971        0.0000000000        0.0000001250
    4       -93.7862670014     -116.6492288700        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  2.753536493179471E-004
 Total number of micro iterations:    9

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.99998376 pnorm= 0.0000E+00 rznorm= 2.6213E-07 rpnorm= 0.0000E+00 noldr=  9 nnewr=  9 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.716184  -23.217481   -2.316214   -2.165804   -1.848355

 qvv(*) eigenvalues. symmetry block  1
    -0.414417   -0.199287   -0.163239   -0.122373   -0.086724   -0.038568   -0.022179    0.026897    0.043623    0.068857
     0.082004    0.095216    0.143365    0.201958    0.217106    0.279514    0.321403    0.392520    0.417418    0.447051
     0.471656    0.497498    0.528113    0.563404    0.621388    0.649756    0.671388    0.732551    0.801589    0.813220
     0.907392    0.920164    0.959338    0.970028    1.000654    1.199183    1.282390    1.313872    1.325201    1.336294
     1.362506    1.592298    1.799839    1.907348    1.915700    2.035562    2.479849    2.537279    2.611067    2.780320
     2.786095    2.860084    2.869979    3.040668    3.098237    3.275638    3.386829    3.417107    3.497159    3.784613
     4.179971    4.258356    4.264559    4.341055    4.488479    4.734970    4.849319    5.012406    5.027151    5.242829
     5.929862    6.234129    6.385113

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=  100 emc=    -93.9741423714 demc= 6.7580E-06 wnorm= 2.2028E-03 knorm= 5.6994E-03 apxde= 3.2204E-06    *not conv.*     

 final mcscf convergence values:
 iter=  100 emc=    -93.9741423714 demc= 6.7580E-06 wnorm= 2.2028E-03 knorm= 5.6994E-03 apxde= 3.2204E-06    *not conv.*     




   ---------Individual total energies for all states:----------
   DRT #1 state # 1 wt 0.333 total energy=      -94.081620660, rel. (eV)=   0.000000
   DRT #1 state # 2 wt 0.333 total energy=      -93.938454825, rel. (eV)=   3.895742
   DRT #1 state # 3 wt 0.333 total energy=      -93.902351629, rel. (eV)=   4.878161
   ------------------------------------------------------------


 MO-coefficient print-out skipped (no flag 32)
 They may be found in the MOCOEF directory.

          natural orbitals of the final iteration,block  1    -  A  
               MO    1        MO    2        MO    3        MO    4        MO    5        MO    6        MO    7        MO    8
  occ(*)=     2.00000000     2.00000000     2.00000000     2.00000000     2.00000000     1.99268301     1.66369250     1.33339293
               MO    9        MO   10        MO   11        MO   12        MO   13        MO   14        MO   15        MO   16
  occ(*)=     1.01023156     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
               MO   17        MO   18        MO   19        MO   20        MO   21        MO   22        MO   23        MO   24
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
               MO   25        MO   26        MO   27        MO   28        MO   29        MO   30        MO   31        MO   32
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
               MO   33        MO   34        MO   35        MO   36        MO   37        MO   38        MO   39        MO   40
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
               MO   41        MO   42        MO   43        MO   44        MO   45        MO   46        MO   47        MO   48
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
               MO   49        MO   50        MO   51        MO   52        MO   53        MO   54        MO   55        MO   56
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
               MO   57        MO   58        MO   59        MO   60        MO   61        MO   62        MO   63        MO   64
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
               MO   65        MO   66        MO   67        MO   68        MO   69        MO   70        MO   71        MO   72
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
               MO   73        MO   74        MO   75        MO   76        MO   77        MO   78        MO   79        MO   80
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
               MO   81        MO   82
  occ(*)=     0.00000000     0.00000000
 d1(*), fmc(*), and qmc(*) written to the 1-particle density matrix file.
        180 d2(*) elements written to the 2-particle density matrix file: mcd2fl                                                      


          Mulliken population analysis


  NOTE: For HERMIT use spherical harmonics basis sets !!!
 

                        A   partial gross atomic populations
   ao class       1A         2A         3A         4A         5A         6A  
     2_ s       0.000003   1.993633   1.337153   0.000267   0.004661  -0.000630
     2_ p       0.000004  -0.000000   0.038323   0.000087   0.804725  -0.000168
     2_ d       0.000000   0.000001  -0.003187   0.000034  -0.000243   0.000195
     4_ s       1.998047  -0.000001   0.000411   0.060983   0.000486   1.683404
     4_ p       0.000061  -0.000006   0.000877   0.911653   0.000892   0.346595
     4_ d       0.000001  -0.000000   0.000027   0.003190  -0.000080   0.003152
     5_ s       0.000001   0.002536   0.373589   0.000082   0.584462  -0.000015
     5_ p      -0.000000   0.000781   0.012885   0.000009  -0.023066   0.000025
     6_ s       0.000417  -0.000002   0.000094   0.199469   0.000414   0.006182
     6_ p       0.000272  -0.000000  -0.000004   0.033059   0.000018   0.003345
     7_ s       0.001031   0.000000  -0.000170   0.759369  -0.000571  -0.055118
     7_ p       0.000176   0.000000  -0.000013   0.032119  -0.000064   0.002750
     8_ s      -0.000012   0.002268   0.245516  -0.000280   0.639050   0.002136
     8_ p       0.000000   0.000789  -0.005501  -0.000041  -0.010686   0.000831
 
   ao class       7A         8A         9A        10A        11A        12A  
     2_ s       0.000243   0.095657   0.000030   0.000000   0.000000   0.000000
     2_ p       0.000449   1.200235   0.000114   0.000000   0.000000   0.000000
     2_ d       0.000129  -0.000583   0.000079   0.000000   0.000000   0.000000
     4_ s       0.048537   0.000416   0.000008   0.000000   0.000000   0.000000
     4_ p       1.128912   0.001321   0.962724   0.000000   0.000000   0.000000
     4_ d       0.014738   0.000059   0.002536   0.000000   0.000000   0.000000
     5_ s       0.000024   0.009893  -0.000031   0.000000   0.000000   0.000000
     5_ p       0.000013   0.007281  -0.000027   0.000000   0.000000   0.000000
     6_ s       0.357948   0.000273   0.000002   0.000000   0.000000   0.000000
     6_ p       0.015999   0.000004   0.021869   0.000000   0.000000   0.000000
     7_ s       0.066279  -0.000063  -0.000002   0.000000   0.000000   0.000000
     7_ p       0.031288   0.000026   0.022741   0.000000   0.000000   0.000000
     8_ s      -0.000731   0.012912  -0.000082   0.000000   0.000000   0.000000
     8_ p      -0.000134   0.005961   0.000270   0.000000   0.000000   0.000000
 
   ao class      13A        14A        15A        16A        17A        18A  
 
   ao class      19A        20A        21A        22A        23A        24A  
 
   ao class      25A        26A        27A        28A        29A        30A  
 
   ao class      31A        32A        33A        34A        35A        36A  
 
   ao class      37A        38A        39A        40A        41A        42A  
 
   ao class      43A        44A        45A        46A        47A        48A  
 
   ao class      49A        50A        51A        52A        53A        54A  
 
   ao class      55A        56A        57A        58A        59A        60A  
 
   ao class      61A        62A        63A        64A        65A        66A  
 
   ao class      67A        68A        69A        70A        71A        72A  
 
   ao class      73A        74A        75A        76A        77A        78A  
 
   ao class      79A        80A        81A        82A  


                        gross atomic populations
     ao            2_         4_         5_         6_         7_         8_
      s         3.431017   3.792290   0.970540   0.564798   0.770756   0.900778
      p         2.043768   3.353028  -0.002100   0.074562   0.089022  -0.008511
      d        -0.003574   0.023625   0.000000   0.000000   0.000000   0.000000
    total       5.471211   7.168943   0.968440   0.639361   0.859778   0.892267
 

 Total number of electrons:   16.00000000

 !timer: mcscf                           cpu_time=    32.492 walltime=    32.646
