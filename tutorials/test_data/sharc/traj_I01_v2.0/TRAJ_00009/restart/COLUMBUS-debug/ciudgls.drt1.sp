1
 program ciudg      
 multireference single and double excitation configuration
 interaction based on the graphical unitary group approach.


 references:  h. lischka, r. shepard, f. b. brown, and i. shavitt,
                  int. j. quantum chem. s 15, 91 (1981).
              r. shepard, r. a. bair, r. a. eades, a. f. wagner,
                  m. j. davis, l. b. harding, and t. h. dunning,
                  int j. quantum chem. s 17, 613 (1983).
              r. ahlrichs, h.-j. boehm, c. ehrhardt, p. scharf,
                  h. schiffer, h. lischka, and m. schindler,
                  j. comp. chem. 6, 200 (1985).
              r. shepard, i. shavitt, r. m. pitzer, d. c. comeau, m. pepper
                  h. lischka, p. g. szalay, r. ahlrichs, f. b. brown, and
                  j.-g. zhao, int. j. quantum chem. symp. 22, 149 (1988).

 This Version of Program CIUDG is Maintained by:
     Thomas Mueller
     Juelich Supercomputing Centre (JSC)
     Institute of Advanced Simulation (IAS)
     D-52425 Juelich, Germany 
     Email: th.mueller@fz-juelich.de



     ******************************************
     **    PROGRAM:              CIUDG       **
     **    PROGRAM VERSION:      2009-03.    **
     **    DISTRIBUTION VERSION: 5.9.a       **
     ******************************************


================ Computing sorting integral file structure ================

                    -----z----- -----y----- -----x----- -----w----- ---total---

                CSFs        28        7592      275940      378140      661700
      internal walks        28         104         105         140         377
valid internal walks        28         104         105         140         377
 getinfoarray: info=                     6 :                     1
                  8192                  6552                  8192
                  5460                     0
 icd(3)=                  6481 ci%nnlev=                  3240  l2rec=
                  8192  n2max=                  5460
 lcore1,lcore2=             393170787             393152761
 lencor,maxblo             393216000                 60000
========================================
 current settings:
 minbl3        8100
 minbl4        8100
 locmaxbl3   394204
 locmaxbuf   197102
 maxbl3       60000
 maxbl3       60000
 maxbl4       60000
 maxbuf       30006
========================================

 sorted 4-external integrals:   128 records of integral w-combinations 
                                128 records of integral x-combinations of length= 30006
                                        have been written (wstat,xstat=   95.31   90.62)

 sorted 3-external integrals:    56 records of integral w-combinations 
                                 49 records of integral x-combinations of length= 30006
                                        have been written (wstat,xstat=   82.14   91.84)
 Orig.  diagonal integrals:  1electron:        80
                             0ext.    :        56
                             2ext.    :      1022
                             4ext.    :      5402


 Orig. off-diag. integrals:  4ext.    :   3840822
                             3ext.    :   1416492
                             2ext.    :    226884
                             1ext.    :     18396
                             0ext.    :       546
                             2ext. SO :         0
                             1ext. SO :         0
                             0ext. SO :         0
                             1electron:      2701


 Sorted integrals            3ext.  w :   1379700 x :   1342908
                             4ext.  w :   3646350 x :   3451878


Cycle #  1 sortfile size=   8355585(     255 records of    32767) #buckets=   2
distributed memory consumption per node=         0 available core 393170787
Cycle #  2 sortfile size=   8355585(     255 records of    32767) #buckets=   2
distributed memory consumption per node=         0 available core 393170787
 minimum size of srtscr:   8290051 WP (   253 records)
 maximum size of srtscr:   8355585 WP (   255 records)
 compressed index vector length=                     1
 echo of the input for program ciudg:
 ------------------------------------------------------------------------
  &input
   NTYPE = 0,
   GSET = 0,
    DAVCOR =10,
   NCOREL = 12
   NROOT = 3
   IVMODE = 3
   NBKITR = 1
   NVBKMN = 3
   RTOLBK = 1e-4,1e-4,1e-4,
   NITER = 90
   NVCIMN = 5
   RTOLCI = 0.000100, 0.000100, 0.000100,
   NVCIMX = 8
   NVRFMX = 8
   NVBKMX = 8
   iden=2
   CSFPRN = 10,
  /&end
 transition
   1  1  1  2
   1  1  1  3
   1  2  1  3
 ------------------------------------------------------------------------
transition densities: resetting nroot to    3
lodens (list->root)=  1  2  3
invlodens (root->list)=  1  2  3
 USING SEGMENTS OF EQUAL SIZE

****************  list of control variables  ****************
 lvlprt =    0      nroot  =    3      noldv  =   0      noldhv =   0
 nunitv =    3      nbkitr =    1      niter  =  90      davcor =  10
 csfprn =   10      ivmode =    3      istrt  =   0      vout   =   0
 iortls =    0      nvbkmx =    8      ibktv  =  -1      ibkthv =  -1
 nvcimx =    8      icitv  =   -1      icithv =  -1      frcsub =   0
 nvbkmn =    3      nvcimn =    5      maxseg =   4      nrfitr =  30
 ncorel =   12      nvrfmx =    8      nvrfmn =   5      iden   =   2
 itran  =    0      froot  =    0      rtmode =   0      ncouple=   1
 skipso =    F      dalton2=    0      molcas =   0      finalv =   0
 finalw =    0      cosmocalc=   0    with_tsklst=   0
 nsegwx =    1     1     1     1
 nseg0x =    1     1     1     1
 nseg1x =    1     1     1     1
 nseg2x =    1     1     1     1
 nseg3x =    1     1     1     1
 nseg4x =    1     1     1     1
 no0ex  =      0    no1ex  =      0    no2ex  =     0    no3ex  =     0
 no4ex  =      0    nodiag =      0
 cdg4ex =    1      c3ex1ex=    1      c2ex0ex=   1
 fileloc=    0     0     0     0     0     0     0     1     1     1
 directhd=   1      noaqccshift_zyxw=      0
 critical_crit=-1.00000    critical_delta= 0.05000

 ctol   = 0.010000    lrtshift=1.000000    smalld =0.001000


 convergence tolerances of bk and full diagonalization steps
 root #       rtolbk        rtol
 ------      --------      ------
    1        1.000E-04    1.000E-04
    2        1.000E-04    1.000E-04
    3        1.000E-04    1.000E-04
 Computing density:                    .drt1.state1
 Computing density:                    .drt1.state2
 Computing density:                    .drt1.state3
 Computing transition density:          drt1.state1-> drt1.state2 (.trd1to2)
 Computing transition density:          drt1.state1-> drt1.state3 (.trd1to3)
 Computing transition density:          drt1.state2-> drt1.state3 (.trd2to3)
 using                      1  nodes and                      1  cores.
 szdg/szodg per processor=                  9192                250482
 Main memory management:
 global                1 DP per process
 vdisk                 0 DP per process
 stack                 0 DP per process
 core          393215999 DP per process

********** Integral sort section *************


 workspace allocation information: lencor= 393215999

 echo of the input for program cisrt:
 ------------------------------------------------------------------------
  &input
  maxbl3=60000
  maxbl4=60000
  &end
 ------------------------------------------------------------------------
 
 ( 6) listing file:                    ciudgls             
 ( 5) input file:                      cisrtin   
 (17) cidrt file:                      cidrtfl             
 (11) transformed integrals file:      moints    
 (12) diagonal integral file:          diagint             
 (13) off-diagonal integral file:      ofdgint             
 (31) 4-external w integrals file:     fil4w               
 (32) 4-external x integrals file:     fil4x               
 (33) 3-external w integrals file:     fil3w               
 (34) 3-external x integrals file:     fil3x               
 (21) scratch da sorting file:         srtscr              
 (12) 2-e integral file [fsplit=2]:    moints2   

 input integral file header information:
 Hermit Integral Program : SIFS version  polonium20        05:27:01.829 17-Feb-18
  cidrt_title                                                                    
 MO-coefficients from mcscf.x                                                    
  with dummy occupation 1.0 for active orbitals                                  
  total ao core energy =   22.862961869                                          
 MCSCF energy =     -93.974142371                                                
 SIFS file created by program tran.      polonium20        05:27:36.619 17-Feb-18

 input energy(*) values:
 energy( 1)=  2.286296186856E+01, ietype=   -1,    core energy of type: Nuc.Rep.

 total core energy =   2.286296186856E+01

 nsym = 1 nmot=  80

 symmetry  =    1
 slabel(*) =  A  
 nmpsy(*)  =   80

 info(*) =          1      8192      6552      8192      5460         0

 orbital labels, i:molab(i)=
   1:tout:001   2:tout:002   3:tout:003   4:tout:004   5:tout:005   6:tout:006   7:tout:007   8:tout:008   9:tout:009  10:tout:010
  11:tout:011  12:tout:012  13:tout:013  14:tout:014  15:tout:015  16:tout:016  17:tout:017  18:tout:018  19:tout:019  20:tout:020
  21:tout:021  22:tout:022  23:tout:023  24:tout:024  25:tout:025  26:tout:026  27:tout:027  28:tout:028  29:tout:029  30:tout:030
  31:tout:031  32:tout:032  33:tout:033  34:tout:034  35:tout:035  36:tout:036  37:tout:037  38:tout:038  39:tout:039  40:tout:040
  41:tout:041  42:tout:042  43:tout:043  44:tout:044  45:tout:045  46:tout:046  47:tout:047  48:tout:048  49:tout:049  50:tout:050
  51:tout:051  52:tout:052  53:tout:053  54:tout:054  55:tout:055  56:tout:056  57:tout:057  58:tout:058  59:tout:059  60:tout:060
  61:tout:061  62:tout:062  63:tout:063  64:tout:064  65:tout:065  66:tout:066  67:tout:067  68:tout:068  69:tout:069  70:tout:070
  71:tout:071  72:tout:072  73:tout:073  74:tout:074  75:tout:075  76:tout:076  77:tout:077  78:tout:078  79:tout:079  80:tout:080

 input parameters:
 prnopt=  0
 ldamin=    4095 ldamax=   32767 ldainc=      64
 maxbuf=   30006 maxbl3=   60000 maxbl4=   60000 intmxo=     766
  Using 32 bit compression 

 drt information:
  cidrt_title                                                                    
 nmotd =  82 nfctd =   2 nfvtc =   0 nmot  =  80
 nlevel =  80 niot  =   7 lowinl=  74
 orbital-to-level map(*)
   -1  -1  74  75  76  77  78  79  80   1   2   3   4   5   6   7   8   9  10  11
   12  13  14  15  16  17  18  19  20  21  22  23  24  25  26  27  28  29  30  31
   32  33  34  35  36  37  38  39  40  41  42  43  44  45  46  47  48  49  50  51
   52  53  54  55  56  57  58  59  60  61  62  63  64  65  66  67  68  69  70  71
   72  73
 compressed map(*)
   74  75  76  77  78  79  80   1   2   3   4   5   6   7   8   9  10  11  12  13
   14  15  16  17  18  19  20  21  22  23  24  25  26  27  28  29  30  31  32  33
   34  35  36  37  38  39  40  41  42  43  44  45  46  47  48  49  50  51  52  53
   54  55  56  57  58  59  60  61  62  63  64  65  66  67  68  69  70  71  72  73
 levsym(*)
    1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1
    1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1
    1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1
    1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1
 repartitioning mu(*)=
   2.  2.  2.  0.  0.  0.  0.
  ... transition density matrix calculation requested 
  ... setting rmu(*) to zero 

 new core energy added to the energy(*) list.
 from the integral file: h1_core= -8.497667908435E+01

 indxdg: diagonal integral statistics.
 total number of integrals contributing to diagonal matrix elements:      6480
 number with all external indices:      5402
 number with half external - half internal indices:      1022
 number with all internal indices:        56

 indxof: off-diagonal integral statistics.
    4-external integrals: num=    3840822 strt=          1
    3-external integrals: num=    1416492 strt=    3840823
    2-external integrals: num=     226884 strt=    5257315
    1-external integrals: num=      18396 strt=    5484199
    0-external integrals: num=        546 strt=    5502595

 total number of off-diagonal integrals:     5503140


 indxof(2nd)  ittp=   3 numx(ittp)=      226884
 indxof(2nd)  ittp=   4 numx(ittp)=       18396
 indxof(2nd)  ittp=   5 numx(ittp)=         546

 intermediate da file sorting parameters:
 nbuk=   2 lendar=   32767 nipbk=   21844 nipsg= 393038133
 pro2e        1    3241    6481    9721   12961   12989   13017   16257   59945  103633
   136400  144592  150052  171891

 pro2e:   5250420 integrals read in   962 records.
 pro1e        1    3241    6481    9721   12961   12989   13017   16257   59945  103633
   136400  144592  150052  171891
 pro1e: eref =    0.000000000000000E+00
 total size of srtscr:                   254  records of                  32767 
 WP =              66582544 Bytes
 putdg        1    3241    6481    9721   10487   43254   65099   16257   59945  103633
   136400  144592  150052  171891

 putf:      12 buffers of length     766 written to file 12
 diagonal integral file completed.

 putd34:   128 records of integral w-combinations and
           128 records of integral x-combinations of length= 30006 have been written.
 wstat,xstat=   95.31   90.62
 prep4e:   128 blocks of linear combinations of 4-external integrals processed.
 number of sorted 4-external integrals    7098228
 number of original 4-external integrals  3840822


 putf34: external integral file complete. nfilw=    31 nfilx=    32 nrecw=   128 nrecx=   128 lbufp= 30006

 putd34:    56 records of integral w-combinations and
            49 records of integral x-combinations of length= 30006 have been written.
 wstat,xstat=   82.14   91.84
 prep3e:    56 blocks of linear combinations of 3-external integrals processed.
 number of sorted 3-external integrals    2722608
 number of original 3-external integrals  1416492


 putf34: external integral file complete. nfilw=    33 nfilx=    34 nrecw=    56 nrecx=    49 lbufp= 30006
ptofdgf: num,ittp,ipos,istrtx,numx,maxrd    226884         3   5257315   5257315    226884   5503140
ptofdgf: num,ittp,ipos,istrtx,numx,maxrd     18396         4   5484199   5484199     18396   5503140
ptofdgf: num,ittp,ipos,istrtx,numx,maxrd       546         5   5502595   5502595       546   5503140

 putf:     327 buffers of length     766 written to file 13
 off-diagonal files sort completed.
 executing brd_struct for cisrtinfo
cisrtinfo:
bufszi   766
 diagfile 4ext:    5402 2ext:    1022 0ext:      56
 fil4w,fil4x  : 3840822 fil3w,fil3x : 1416492
 ofdgint  2ext:  226884 1ext:   18396 0ext:     546so0ext:       0so1ext:       0so2ext:       0
buffer minbl4    8100 minbl3    8100 maxbl2    8103nbas:  73   0   0   0   0   0   0   0 maxbuf 30006
 CIUDG version 5.9.7 ( 5-Oct-2004)

 workspace allocation information: lcore= 393215999

 core energy values from the integral file:
 energy( 1)=  2.286296186856E+01, ietype=   -1,    core energy of type: Nuc.Rep.
 energy( 2)= -8.497667908435E+01, ietype=    6,   fcore energy of type: H1(*)   

 total core repulsion energy = -6.211371721580E+01
 nmot  =    82 niot  =     7 nfct  =     2 nfvt  =     0
 nrow  =    43 nsym  =     1 ssym  =     1 lenbuf=  1600
 nwalk,xbar:        377       28      104      105      140
 nvalwt,nvalw:      377       28      104      105      140
 ncsft:          661700
 total number of valid internal walks:     377
 nvalz,nvaly,nvalx,nvalw =       28     104     105     140

 cisrt info file parameters:
 file number  12 blocksize    766
 mxbld   6128
 nd4ext,nd2ext,nd0ext  5402  1022    56
 n4ext,n3ext,n2ext,n1ext,n0ext,n2int,n1int,n0int  3840822  1416492   226884    18396      546        0        0        0
 minbl4,minbl3,maxbl2  8100  8100  8103
 maxbuf 30006
 number of external orbitals per symmetry block:  73
 nmsym   1 number of internal orbitals   7
 bummer (warning):transition densities: resetting ref occupation number to 0                      0
 executing brd_struct for drt
 executing brd_struct for orbinf
 executing brd_struct for momap
 calcthrxt: niot,maxw1=                     7                   119
 block size     0
 pthz,pthy,pthx,pthw:    28   104   105   140 total internal walks:     377
 maxlp3,n2lp,n1lp,n0lp   119     0     0     0
 orbsym(*)= 1 1 1 1 1 1 1

 setref:       10 references kept,
                0 references were marked as invalid, out of
               10 total.
 nmb.of records onel     4
 nmb.of records 2-ext   297
 nmb.of records 1-ext    25
 nmb.of records 0-ext     1
 nmb.of records 2-int     0
 nmb.of records 1-int     0
 nmb.of records 0-int     0
 ---------memory usage in DP -----------------
 < n-ex core usage >
     routines:
    fourex            68808
    threx             60727
    twoex             57122
    onex              11717
    allin               766
    diagon             6792
               =======
   maximum            68808
 
  __ static summary __ 
   reflst                28
   hrfspc                28
               -------
   static->              28
 
  __ core required  __ 
   totstc                28
   max n-ex           68808
               -------
   totnec->           68836
 
  __ core available __ 
   totspc         393215999
   totnec -           68836
               -------
   totvec->       393147163

 number of external paths / symmetry
 vertex x    2628
 vertex w    2701
segment: free space=   393147163
 reducing frespc by                  1237 to              393145926 
  for index/conft/indsym storage .
 resegmenting ...



                   segmentation summary for type all-internal
 -------------------------------------------------------------------------------
 seg.      no. of|    no. of|  starting|  internal|  starting|  starting|
  no.    internal|        ci|       csf|     walks|      walk|       DRT|
            paths|  elements|    number|     /seg.|    number|    record|
 -------------------------------------------------------------------------------
  Z 1          28|        28|         0|        28|         0|         1|
 -------------------------------------------------------------------------------
  Y 2         104|      7592|        28|       104|        28|         2|
 -------------------------------------------------------------------------------
  X 3         105|    275940|      7620|       105|       132|         3|
 -------------------------------------------------------------------------------
  W 4         140|    378140|    283560|       140|       237|         4|
 -------------------------------------------------------------------------------
max. additional memory requirements:index=           4DP  conft+indsym=         560DP  drtbuffer=         673 DP

dimension of the ci-matrix ->>>    661700

 executing brd_struct for civct
 gentasklist: ntask=                    20
                    TASKLIST
----------------------------------------------------------------------------------------------------
TASK# BRA# KET#  T-TYPE    DESCR.   SEGMENTTYPE    SEGEL              SEGCI          VWALKS   
----------------------------------------------------------------------------------------------------
     1  3   1    24      two-ext xz   2X  3 1     105      28     275940         28     105      28
     2  4   1    25      two-ext wz   2X  4 1     140      28     378140         28     140      28
     3  4   3    26      two-ext wx*  WX  4 3     140     105     378140     275940     140     105
     4  4   3    27      two-ext wx+  WX  4 3     140     105     378140     275940     140     105
     5  2   1    11      one-ext yz   1X  2 1     104      28       7592         28     104      28
     6  3   2    15      1ex3ex yx    3X  3 2     105     104     275940       7592     105     104
     7  4   2    16      1ex3ex yw    3X  4 2     140     104     378140       7592     140     104
     8  1   1     1      allint zz    OX  1 1      28      28         28         28      28      28
     9  2   2     5      0ex2ex yy    OX  2 2     104     104       7592       7592     104     104
    10  3   3     6      0ex2ex xx*   OX  3 3     105     105     275940     275940     105     105
    11  3   3    18      0ex2ex xx+   OX  3 3     105     105     275940     275940     105     105
    12  4   4     7      0ex2ex ww*   OX  4 4     140     140     378140     378140     140     140
    13  4   4    19      0ex2ex ww+   OX  4 4     140     140     378140     378140     140     140
    14  2   2    42      four-ext y   4X  2 2     104     104       7592       7592     104     104
    15  3   3    43      four-ext x   4X  3 3     105     105     275940     275940     105     105
    16  4   4    44      four-ext w   4X  4 4     140     140     378140     378140     140     140
    17  1   1    75      dg-024ext z  OX  1 1      28      28         28         28      28      28
    18  2   2    76      dg-024ext y  OX  2 2     104     104       7592       7592     104     104
    19  3   3    77      dg-024ext x  OX  3 3     105     105     275940     275940     105     105
    20  4   4    78      dg-024ext w  OX  4 4     140     140     378140     378140     140     140
----------------------------------------------------------------------------------------------------
REDTASK #   1 TIME=  19.000 N=  1 (task/type/sgbra)=(   1/24/0) (
REDTASK #   2 TIME=  18.000 N=  1 (task/type/sgbra)=(   2/25/0) (
REDTASK #   3 TIME=  17.000 N=  1 (task/type/sgbra)=(   3/26/1) (
REDTASK #   4 TIME=  16.000 N=  1 (task/type/sgbra)=(   4/27/2) (
REDTASK #   5 TIME=  15.000 N=  1 (task/type/sgbra)=(   5/11/0) (
REDTASK #   6 TIME=  14.000 N=  1 (task/type/sgbra)=(   6/15/0) (
REDTASK #   7 TIME=  13.000 N=  1 (task/type/sgbra)=(   7/16/0) (
REDTASK #   8 TIME=  12.000 N=  1 (task/type/sgbra)=(   8/ 1/0) (
REDTASK #   9 TIME=  11.000 N=  1 (task/type/sgbra)=(   9/ 5/0) (
REDTASK #  10 TIME=  10.000 N=  1 (task/type/sgbra)=(  10/ 6/1) (
REDTASK #  11 TIME=   9.000 N=  1 (task/type/sgbra)=(  11/18/2) (
REDTASK #  12 TIME=   8.000 N=  1 (task/type/sgbra)=(  12/ 7/1) (
REDTASK #  13 TIME=   7.000 N=  1 (task/type/sgbra)=(  13/19/2) (
REDTASK #  14 TIME=   6.000 N=  1 (task/type/sgbra)=(  14/42/1) (
REDTASK #  15 TIME=   5.000 N=  1 (task/type/sgbra)=(  15/43/1) (
REDTASK #  16 TIME=   4.000 N=  1 (task/type/sgbra)=(  16/44/1) (
REDTASK #  17 TIME=   3.000 N=  1 (task/type/sgbra)=(  17/75/1) (
REDTASK #  18 TIME=   2.000 N=  1 (task/type/sgbra)=(  18/76/1) (
REDTASK #  19 TIME=   1.000 N=  1 (task/type/sgbra)=(  19/77/1) (
REDTASK #  20 TIME=   0.000 N=  1 (task/type/sgbra)=(  20/78/1) (
 initializing v-file: 1:                661700

    ---------trial vector generation----------

    trial vectors will be created by: 

    (ivmode= 3) diagonalizing h in the reference space.                     

      3 vectors will be written to unit 11 beginning with logical record   1

            3 vectors will be created
 ========= Executing OUT-OF-CORE method ========


====================================================================================================
Diagonal     counts:  0x:         735 2x:           0 4x:           0
All internal counts: zz :        1113 yy:           0 xx:           0 ww:           0
One-external counts: yz :           0 yx:           0 yw:           0
Two-external counts: yy :           0 ww:           0 xx:           0 xz:           0 wz:           0 wx:           0
Three-ext.   counts: yx :           0 yw:           0

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================


 reference space has dimension      10
 dsyevx: computed roots 1 to    6(converged:   6)

    root           eigenvalues
    ----           ------------
       1         -94.0816206603
       2         -93.9384548255
       3         -93.9023516285
       4         -93.7862670014
       5         -93.6680015243
       6         -93.6116790275

 strefv generated    3 initial ci vector(s).
    ---------end of vector generation---------

 ufvoutnew: ... writing  recamt=                    28

         vector  1 from unit 11 written to unit 49 filename cirefv              
 ufvoutnew: ... writing  recamt=                    28

         vector  2 from unit 11 written to unit 49 filename cirefv              
 ufvoutnew: ... writing  recamt=                    28

         vector  3 from unit 11 written to unit 49 filename cirefv              

 ************************************************************************
 beginning the bk-type iterative procedure (nzcsf=    28)...
 ************************************************************************

               initial diagonalization conditions:

 number of configuration state functions:            661700
 number of initial trial vectors:                         3
 number of initial matrix-vector products:                0
 maximum dimension of the subspace vectors:               8
 number of roots to converge:                             3
 number of iterations:                                    1
 residual norm convergence criteria:               0.000100  0.000100  0.000100

          starting bk iteration   1

 ========= Executing OUT-OF-CORE method ========


====================================================================================================
Diagonal     counts:  0x:        9055 2x:        2233 4x:         349
All internal counts: zz :        1113 yy:           0 xx:           0 ww:           0
One-external counts: yz :        3142 yx:           0 yw:           0
Two-external counts: yy :           0 ww:           0 xx:           0 xz:         403 wz:         604 wx:           0
Three-ext.   counts: yx :           0 yw:           0

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================


 ========= Executing OUT-OF-CORE method ========


====================================================================================================
Diagonal     counts:  0x:        9055 2x:        2233 4x:         349
All internal counts: zz :        1113 yy:           0 xx:           0 ww:           0
One-external counts: yz :        3142 yx:           0 yw:           0
Two-external counts: yy :           0 ww:           0 xx:           0 xz:         403 wz:         604 wx:           0
Three-ext.   counts: yx :           0 yw:           0

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================


 ========= Executing OUT-OF-CORE method ========


====================================================================================================
Diagonal     counts:  0x:        9055 2x:        2233 4x:         349
All internal counts: zz :        1113 yy:           0 xx:           0 ww:           0
One-external counts: yz :        3142 yx:           0 yw:           0
Two-external counts: yy :           0 ww:           0 xx:           0 xz:         403 wz:         604 wx:           0
Three-ext.   counts: yx :           0 yw:           0

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================


 Final subspace hamiltonian 

                ht   1         ht   2         ht   3
   ht   1   -31.96790344
   ht   2    -0.00000000   -31.82473761
   ht   3    -0.00000000     0.00000000   -31.78863441

          calcsovref: tciref block   1

              civs   1       civs   2       civs   3
 refs   1   1.000000        0.00000      -1.490778E-18
 refs   2    0.00000       1.000000       2.077772E-18
 refs   3    0.00000        0.00000       1.000000    

          calcsovref: scrb block   1

                ci   1         ci   2         ci   3
 civs   1  -1.000000      -6.112187E-14  -6.938894E-18
 civs   2  -6.116442E-14   1.000000      -1.517883E-18
 civs   3    0.00000        0.00000       1.000000    

          calcsovref: sovref block   1

              v      1       v      2       v      3
 ref    1  -1.000000      -6.112187E-14  -8.429672E-18
 ref    2  -6.116442E-14   1.000000       5.598888E-19
 ref    3    0.00000        0.00000       1.000000    

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3
 ref:   1    -1.00000000    -0.00000000    -0.00000000
 ref:   2    -0.00000000     1.00000000     0.00000000
 ref:   3     0.00000000     0.00000000     1.00000000

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  1  1    -94.0816206603 -3.5527E-15  3.0928E-01  9.4410E-01  1.0000E-04
 mr-sdci #  1  2    -93.9384548255  1.0658E-14  0.0000E+00  9.4220E-01  1.0000E-04
 mr-sdci #  1  3    -93.9023516285  2.1316E-14  0.0000E+00  9.5573E-01  1.0000E-04
 
 root number  1 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.003999
time for cinew                         0.038994
time for eigenvalue solver             0.000000
time for vector access                 0.000000

 mr-sdci  convergence not reached after  1 iterations.

 final mr-sdci  convergence information:

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  1  1    -94.0816206603 -3.5527E-15  3.0928E-01  9.4410E-01  1.0000E-04
 mr-sdci #  1  2    -93.9384548255  1.0658E-14  0.0000E+00  9.4220E-01  1.0000E-04
 mr-sdci #  1  3    -93.9023516285  2.1316E-14  0.0000E+00  9.5573E-01  1.0000E-04
 
diagon:itrnv=   2
 expansion vectors are not transformed.
 matrix-vector products are not transformed.

    4 expansion eigenvectors written to unit nvfile (= 11)
    3 matrix-vector products written to unit nhvfil (= 10)

 ************************************************************************
 beginning the ci iterative diagonalization procedure... 
 ************************************************************************

               initial diagonalization conditions:

 number of configuration state functions:            661700
 number of initial trial vectors:                         4
 number of initial matrix-vector products:                3
 maximum dimension of the subspace vectors:               8
 number of roots to converge:                             3
 number of iterations:                                   90
 residual norm convergence criteria:               0.000100  0.000100  0.000100

          starting ci iteration   1

 ========= Executing OUT-OF-CORE method ========


====================================================================================================
Diagonal     counts:  0x:        9055 2x:        2233 4x:         349
All internal counts: zz :        1113 yy:        8560 xx:       15526 ww:       21792
One-external counts: yz :        3142 yx:        7351 yw:        8696
Two-external counts: yy :        1938 ww:        4100 xx:        3262 xz:         403 wz:         604 wx:        4304
Three-ext.   counts: yx :        4024 yw:        4720

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================


 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4
   ht   1   -31.96790344
   ht   2    -0.00000000   -31.82473761
   ht   3    -0.00000000     0.00000000   -31.78863441
   ht   4     0.30927858    -0.00011004    -0.00069780    -4.16556153

          calcsovref: tciref block   1

              civs   1       civs   2       civs   3       civs   4
 refs   1   1.000000        0.00000      -1.490778E-18   1.996800E-12
 refs   2    0.00000       1.000000       2.077772E-18  -8.020441E-15
 refs   3    0.00000        0.00000       1.000000      -6.915184E-14

          calcsovref: scrb block   1

                ci   1         ci   2         ci   3         ci   4
 civs   1   0.951317       2.455464E-04   1.453094E-03   0.308209    
 civs   2  -2.204670E-04    1.00000      -2.567109E-06  -1.161826E-04
 civs   3  -1.285115E-03   2.196876E-06   0.999999      -7.480012E-04
 civs   4  -0.822798       1.136640E-04   8.422633E-04    2.53964    

          calcsovref: sovref block   1

              v      1       v      2       v      3       v      4
 ref    1   0.951317       2.455464E-04   1.453094E-03   0.308209    
 ref    2  -2.204670E-04    1.00000      -2.567109E-06  -1.161826E-04
 ref    3  -1.285115E-03   2.196876E-06   0.999999      -7.480012E-04

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4
 ref:   1     0.95131737     0.00024555     0.00145309     0.30820948
 ref:   2    -0.00022047     0.99999997    -0.00000257    -0.00011618
 ref:   3    -0.00128512     0.00000220     0.99999889    -0.00074800

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  1  1    -94.3491167983  2.6750E-01  1.1965E-02  1.8206E-01  1.0000E-04
 mr-sdci #  1  2    -93.9384548380  1.2507E-08  0.0000E+00  9.4220E-01  1.0000E-04
 mr-sdci #  1  3    -93.9023522162  5.8773E-07  0.0000E+00  9.5573E-01  1.0000E-04
 mr-sdci #  1  4    -91.5331731734  2.9419E+01  0.0000E+00  1.4567E+00  1.0000E-04
 
 root number  1 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.004999
time for cinew                         0.052992
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration   2

 ========= Executing OUT-OF-CORE method ========


====================================================================================================
Diagonal     counts:  0x:        9055 2x:        2233 4x:         349
All internal counts: zz :        1113 yy:        8560 xx:       15526 ww:       21792
One-external counts: yz :        3142 yx:        7351 yw:        8696
Two-external counts: yy :        1938 ww:        4100 xx:        3262 xz:         403 wz:         604 wx:        4304
Three-ext.   counts: yx :        4024 yw:        4720

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================


 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5
   ht   1   -31.96790344
   ht   2    -0.00000000   -31.82473761
   ht   3    -0.00000000     0.00000000   -31.78863441
   ht   4     0.30927858    -0.00011004    -0.00069780    -4.16556153
   ht   5    -0.00105410    -0.00012959     0.00081720    -0.12877435    -0.16437773

          calcsovref: tciref block   1

              civs   1       civs   2       civs   3       civs   4       civs   5
 refs   1   1.000000        0.00000      -1.490778E-18   1.996800E-12  -3.874397E-06
 refs   2    0.00000       1.000000       2.077772E-18  -8.020441E-15   3.743258E-06
 refs   3    0.00000        0.00000       1.000000      -6.915184E-14  -2.476237E-05

          calcsovref: scrb block   1

                ci   1         ci   2         ci   3         ci   4         ci   5
 civs   1  -0.948845      -2.207165E-04  -1.546191E-03  -8.037467E-02   0.305337    
 civs   2   2.021152E-04   -1.00000       1.625212E-06  -9.695495E-05  -1.203046E-04
 civs   3   1.336407E-03  -1.322428E-06  -0.999999       7.323115E-04  -7.181208E-04
 civs   4   0.852015      -1.055556E-04  -8.824270E-04  -0.487169        2.51935    
 civs   5  -0.950613      -8.892616E-04   3.624703E-03    13.6068       0.630084    

          calcsovref: sovref block   1

              v      1       v      2       v      3       v      4       v      5
 ref    1  -0.948841      -2.207131E-04  -1.546205E-03  -8.042739E-02   0.305334    
 ref    2   1.985568E-04   -1.00000       1.638780E-06  -4.602106E-05  -1.179461E-04
 ref    3   1.359946E-03  -1.300408E-06  -0.999999       3.953740E-04  -7.337231E-04

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5
 ref:   1    -0.94884144    -0.00022071    -0.00154620    -0.08042739     0.30533412
 ref:   2     0.00019856    -0.99999997     0.00000164    -0.00004602    -0.00011795
 ref:   3     0.00135995    -0.00000130    -0.99999873     0.00039537    -0.00073372

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  2  1    -94.3605184298  1.1402E-02  7.4873E-04  4.4661E-02  1.0000E-04
 mr-sdci #  2  2    -93.9384548464  8.4144E-09  0.0000E+00  9.4220E-01  1.0000E-04
 mr-sdci #  2  3    -93.9023523531  1.3689E-07  0.0000E+00  9.5573E-01  1.0000E-04
 mr-sdci #  2  4    -92.0089520819  4.7578E-01  0.0000E+00  1.2062E+00  1.0000E-04
 mr-sdci #  2  5    -91.5321538049  2.9418E+01  0.0000E+00  1.4631E+00  1.0000E-04
 
 root number  1 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000000
time for cinew                         0.058991
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration   3

 ========= Executing OUT-OF-CORE method ========


====================================================================================================
Diagonal     counts:  0x:        9055 2x:        2233 4x:         349
All internal counts: zz :        1113 yy:        8560 xx:       15526 ww:       21792
One-external counts: yz :        3142 yx:        7351 yw:        8696
Two-external counts: yy :        1938 ww:        4100 xx:        3262 xz:         403 wz:         604 wx:        4304
Three-ext.   counts: yx :        4024 yw:        4720

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================


 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6
   ht   1   -31.96790344
   ht   2    -0.00000000   -31.82473761
   ht   3    -0.00000000     0.00000000   -31.78863441
   ht   4     0.30927858    -0.00011004    -0.00069780    -4.16556153
   ht   5    -0.00105410    -0.00012959     0.00081720    -0.12877435    -0.16437773
   ht   6    -0.00100305    -0.00045284     0.00010649    -0.03410901    -0.00641388    -0.01115798

          calcsovref: tciref block   1

              civs   1       civs   2       civs   3       civs   4       civs   5       civs   6
 refs   1   1.000000        0.00000      -1.490778E-18   1.996800E-12  -3.874397E-06   2.223589E-05
 refs   2    0.00000       1.000000       2.077772E-18  -8.020441E-15   3.743258E-06   1.381677E-05
 refs   3    0.00000        0.00000       1.000000      -6.915184E-14  -2.476237E-05  -6.692768E-07

          calcsovref: scrb block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 civs   1   0.947492       2.908527E-04  -1.099931E-03   0.153209      -2.040467E-02   0.279946    
 civs   2  -2.197730E-04    1.00000      -5.645731E-05  -1.081344E-03   4.107258E-04   3.294668E-04
 civs   3  -1.136563E-03  -5.323126E-05  -0.999995       2.381789E-03  -1.041316E-03  -1.461546E-03
 civs   4  -0.854073       1.549380E-04  -5.119345E-04   0.744734      -0.203465        2.46960    
 civs   5    1.01862      -7.062582E-04  -7.044922E-03   -6.35794       -12.1586      -0.843648    
 civs   6   -1.09903       2.407065E-02   0.159120        45.8893       -14.2728       -22.6601    

          calcsovref: sovref block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.947464       2.913907E-04  -1.096365E-03   0.154254      -2.067493E-02   0.279445    
 ref    2  -2.311451E-04    1.00000      -5.428515E-05  -4.711015E-04   1.680098E-04   1.322018E-05
 ref    3  -1.161051E-03  -5.322988E-05  -0.999995       2.508514E-03  -7.306886E-04  -1.425489E-03

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 ref:   1     0.94746365     0.00029139    -0.00109637     0.15425405    -0.02067493     0.27944514
 ref:   2    -0.00023115     0.99999985    -0.00005429    -0.00047110     0.00016801     0.00001322
 ref:   3    -0.00116105    -0.00005323    -0.99999490     0.00250851    -0.00073069    -0.00142549

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  3  1    -94.3613414855  8.2306E-04  6.5420E-05  1.2914E-02  1.0000E-04
 mr-sdci #  3  2    -93.9384551512  3.0477E-07  0.0000E+00  9.4220E-01  1.0000E-04
 mr-sdci #  3  3    -93.9023653332  1.2980E-05  0.0000E+00  9.5571E-01  1.0000E-04
 mr-sdci #  3  4    -92.6268135141  6.1786E-01  0.0000E+00  9.8557E-01  1.0000E-04
 mr-sdci #  3  5    -91.9324339894  4.0028E-01  0.0000E+00  1.2906E+00  1.0000E-04
 mr-sdci #  3  6    -91.3211142422  2.9207E+01  0.0000E+00  1.3683E+00  1.0000E-04
 
 root number  1 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.005000
time for cinew                         0.061990
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration   4

 ========= Executing OUT-OF-CORE method ========


====================================================================================================
Diagonal     counts:  0x:        9055 2x:        2233 4x:         349
All internal counts: zz :        1113 yy:        8560 xx:       15526 ww:       21792
One-external counts: yz :        3142 yx:        7351 yw:        8696
Two-external counts: yy :        1938 ww:        4100 xx:        3262 xz:         403 wz:         604 wx:        4304
Three-ext.   counts: yx :        4024 yw:        4720

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================


 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7
   ht   1   -31.96790344
   ht   2    -0.00000000   -31.82473761
   ht   3    -0.00000000     0.00000000   -31.78863441
   ht   4     0.30927858    -0.00011004    -0.00069780    -4.16556153
   ht   5    -0.00105410    -0.00012959     0.00081720    -0.12877435    -0.16437773
   ht   6    -0.00100305    -0.00045284     0.00010649    -0.03410901    -0.00641388    -0.01115798
   ht   7    -0.00009887    -0.00018299     0.00019130     0.00322697    -0.00313922    -0.00047006    -0.00109023

          calcsovref: tciref block   1

              civs   1       civs   2       civs   3       civs   4       civs   5       civs   6       civs   7
 refs   1   1.000000        0.00000      -1.490778E-18   1.996800E-12  -3.874397E-06   2.223589E-05   3.959956E-06
 refs   2    0.00000       1.000000       2.077772E-18  -8.020441E-15   3.743258E-06   1.381677E-05   5.726729E-06
 refs   3    0.00000        0.00000       1.000000      -6.915184E-14  -2.476237E-05  -6.692768E-07  -6.800381E-06

          calcsovref: scrb block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 civs   1  -0.947154       2.976679E-04   6.965407E-04   0.126307       8.724259E-02  -6.046625E-02   0.275099    
 civs   2   2.255426E-04    1.00000       6.617745E-05  -1.585794E-04  -1.356904E-03  -5.232384E-04   8.735296E-05
 civs   3   1.053783E-03  -6.528963E-05   0.999986       4.936667E-03   1.080147E-03  -5.393292E-04  -1.633490E-03
 civs   4   0.854148       1.599943E-04   1.658808E-04   0.446457       0.547123      -0.531087        2.44642    
 civs   5   -1.02440      -8.506089E-04   1.610945E-02   -5.81949        1.79534       -12.0205       -4.04976    
 civs   6    1.19289       2.544339E-02  -0.242094        32.4097        36.4702        3.65930       -21.7507    
 civs   7   -1.07535      -1.298149E-02   0.769694       -111.001        82.4572        95.3529        41.7231    

          calcsovref: sovref block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1  -0.947128       2.981856E-04   6.941430E-04   0.126611       8.837311E-02  -5.996072E-02   0.274796    
 ref    2   2.320317E-04    1.00000       6.730063E-05  -3.682400E-04  -3.740733E-04   2.838555E-05   1.060648E-05
 ref    3   1.085664E-03  -6.519732E-05   0.999980       5.813931E-03   4.505408E-04  -8.925588E-04  -1.802384E-03

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 ref:   1    -0.94712789     0.00029819     0.00069414     0.12661075     0.08837311    -0.05996072     0.27479595
 ref:   2     0.00023203     0.99999983     0.00006730    -0.00036824    -0.00037407     0.00002839     0.00001061
 ref:   3     0.00108566    -0.00006520     0.99998038     0.00581393     0.00045054    -0.00089256    -0.00180238

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  4  1    -94.3614118370  7.0351E-05  5.9457E-06  3.8850E-03  1.0000E-04
 mr-sdci #  4  2    -93.9384551587  7.5462E-09  0.0000E+00  9.4220E-01  1.0000E-04
 mr-sdci #  4  3    -93.9023910392  2.5706E-05  0.0000E+00  9.5563E-01  1.0000E-04
 mr-sdci #  4  4    -93.0360017915  4.0919E-01  0.0000E+00  8.0309E-01  1.0000E-04
 mr-sdci #  4  5    -92.2855137242  3.5308E-01  0.0000E+00  1.2117E+00  1.0000E-04
 mr-sdci #  4  6    -91.6046689655  2.8355E-01  0.0000E+00  1.2520E+00  1.0000E-04
 mr-sdci #  4  7    -91.2835161212  2.9170E+01  0.0000E+00  1.3446E+00  1.0000E-04
 
 root number  1 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.004999
time for cinew                         0.063991
time for eigenvalue solver             0.000999
time for vector access                 0.000000

          starting ci iteration   5

 ========= Executing OUT-OF-CORE method ========


====================================================================================================
Diagonal     counts:  0x:        9055 2x:        2233 4x:         349
All internal counts: zz :        1113 yy:        8560 xx:       15526 ww:       21792
One-external counts: yz :        3142 yx:        7351 yw:        8696
Two-external counts: yy :        1938 ww:        4100 xx:        3262 xz:         403 wz:         604 wx:        4304
Three-ext.   counts: yx :        4024 yw:        4720

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================


 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7         ht   8
   ht   1   -31.96790344
   ht   2    -0.00000000   -31.82473761
   ht   3    -0.00000000     0.00000000   -31.78863441
   ht   4     0.30927858    -0.00011004    -0.00069780    -4.16556153
   ht   5    -0.00105410    -0.00012959     0.00081720    -0.12877435    -0.16437773
   ht   6    -0.00100305    -0.00045284     0.00010649    -0.03410901    -0.00641388    -0.01115798
   ht   7    -0.00009887    -0.00018299     0.00019130     0.00322697    -0.00313922    -0.00047006    -0.00109023
   ht   8    -0.00057492    -0.00035793     0.00008340    -0.00015672     0.00008282    -0.00026174    -0.00002724    -0.00010184

          calcsovref: tciref block   1

              civs   1       civs   2       civs   3       civs   4       civs   5       civs   6       civs   7       civs   8
 refs   1   1.000000        0.00000      -1.490778E-18   1.996800E-12  -3.874397E-06   2.223589E-05   3.959956E-06   1.794067E-05
 refs   2    0.00000       1.000000       2.077772E-18  -8.020441E-15   3.743258E-06   1.381677E-05   5.726729E-06   1.124281E-05
 refs   3    0.00000        0.00000       1.000000      -6.915184E-14  -2.476237E-05  -6.692768E-07  -6.800381E-06  -1.734287E-06

          calcsovref: scrb block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 civs   1  -0.947124       3.028352E-04   5.214680E-05   9.529578E-02  -1.156021E-02  -0.170245       8.397195E-02   0.240426    
 civs   2   2.146949E-04   0.999999       2.267379E-04  -3.228966E-03  -3.915474E-03  -3.349669E-03   9.149225E-04  -1.759894E-03
 civs   3   9.843814E-04  -1.294176E-04   0.999784       1.920919E-02   5.589154E-03   5.778431E-03  -5.466151E-05   3.875731E-04
 civs   4   0.854142       1.648641E-04  -5.242117E-04   0.277907      -7.740590E-03   -1.08300       0.713511        2.23776    
 civs   5   -1.02490      -9.862501E-04   3.364272E-02   -3.98909        4.82948        3.68182        11.5207       -3.66108    
 civs   6    1.20112       2.683799E-02  -0.417184        27.2458        15.3090       -34.8366       -2.89323       -28.7262    
 civs   7   -1.16970      -2.163441E-02    1.82006       -100.660        89.4861       -39.5663       -92.1060        43.8698    
 civs   8    1.03817       7.043335E-02   -8.32269        264.273        262.084        384.260       -38.1117        174.393    

          calcsovref: sovref block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8
 ref    1  -0.947079       3.046137E-04  -9.936734E-05   0.100260      -6.182180E-03  -0.164296       8.281450E-02   0.243103    
 ref    2   2.324275E-04    1.00000       1.379522E-04  -4.727285E-04  -2.268475E-04   2.763608E-04  -3.787589E-05   4.139103E-05
 ref    3   1.015110E-03  -1.293862E-04   0.999785       1.951594E-02   4.396250E-03   5.313225E-03   3.544481E-04  -1.033221E-04

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1    -0.94707911     0.00030461    -0.00009937     0.10025970    -0.00618218    -0.16429646     0.08281450     0.24310339
 ref:   2     0.00023243     0.99999979     0.00013795    -0.00047273    -0.00022685     0.00027636    -0.00003788     0.00004139
 ref:   3     0.00101511    -0.00012939     0.99978517     0.01951594     0.00439625     0.00531322     0.00035445    -0.00010332

 trial vector basis is being transformed.  new dimension:   5

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  5  1    -94.3614180097  6.1727E-06  9.1307E-07  1.4375E-03  1.0000E-04
 mr-sdci #  5  2    -93.9384551790  2.0297E-08  0.0000E+00  9.4220E-01  1.0000E-04
 mr-sdci #  5  3    -93.9026663917  2.7535E-04  0.0000E+00  9.5506E-01  1.0000E-04
 mr-sdci #  5  4    -93.2883947322  2.5239E-01  0.0000E+00  6.9884E-01  1.0000E-04
 mr-sdci #  5  5    -92.3927523507  1.0724E-01  0.0000E+00  1.1898E+00  1.0000E-04
 mr-sdci #  5  6    -92.0685668016  4.6390E-01  0.0000E+00  1.2032E+00  1.0000E-04
 mr-sdci #  5  7    -91.6008019112  3.1729E-01  0.0000E+00  1.2718E+00  1.0000E-04
 mr-sdci #  5  8    -91.1784445518  2.9065E+01  0.0000E+00  1.3644E+00  1.0000E-04
 
 root number  1 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.005999
time for cinew                         0.092986
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration   6

 ========= Executing OUT-OF-CORE method ========


====================================================================================================
Diagonal     counts:  0x:        9055 2x:        2233 4x:         349
All internal counts: zz :        1113 yy:        8560 xx:       15526 ww:       21792
One-external counts: yz :        3142 yx:        7351 yw:        8696
Two-external counts: yy :        1938 ww:        4100 xx:        3262 xz:         403 wz:         604 wx:        4304
Three-ext.   counts: yx :        4024 yw:        4720

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================


 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6
   ht   1   -32.24770079
   ht   2     0.00000000   -31.82473796
   ht   3    -0.00000000     0.00000000   -31.78894918
   ht   4    -0.00000000     0.00000000     0.00000000   -31.17467752
   ht   5    -0.00000000    -0.00000000     0.00000000    -0.00000000   -30.27903513
   ht   6    -0.00036590     0.00027138    -0.00006530    -0.00253427     0.00365340    -0.00001798

          calcsovref: tciref block   1

              civs   1       civs   2       civs   3       civs   4       civs   5       civs   6
 refs   1  -0.947079       3.046137E-04  -9.936734E-05   0.100260      -6.182180E-03  -1.124322E-05
 refs   2   2.324275E-04    1.00000       1.379522E-04  -4.727285E-04  -2.268475E-04  -8.547803E-06
 refs   3   1.015110E-03  -1.293862E-04   0.999785       1.951594E-02   4.396250E-03   3.490874E-06

          calcsovref: scrb block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 civs   1   0.999990      -2.419704E-06   4.364544E-04  -5.419703E-03  -4.925132E-04  -1.457485E-02
 civs   2   6.781501E-06  -0.999998      -4.116585E-04   3.877833E-03   3.600916E-04   1.068838E-02
 civs   3  -4.823339E-05  -1.280027E-04   0.999649       2.257822E-02   5.844028E-04   1.435200E-02
 civs   4   1.794763E-04   7.257393E-05  -1.413868E-02   0.910860      -1.741782E-02  -0.422731    
 civs   5   9.672938E-05   2.079628E-05  -3.794081E-03   5.004641E-02  -0.994360       0.184334    
 civs   6   0.828057       0.179545       -32.7864        440.573        41.8023        1242.65    

          calcsovref: sovref block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1  -0.947062      -2.971798E-04  -1.538274E-03   9.119149E-02   4.397513E-03  -4.368855E-02
 ref    2   2.320153E-04   -1.00000       1.414311E-05  -3.281834E-04   2.365417E-04   2.230249E-04
 ref    3   9.736952E-04   3.542928E-06   0.999028       4.210165E-02  -3.981722E-03   1.123108E-02

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 ref:   1    -0.94706193    -0.00029718    -0.00153827     0.09119149     0.00439751    -0.04368855
 ref:   2     0.00023202    -0.99999983     0.00001414    -0.00032818     0.00023654     0.00022302
 ref:   3     0.00097370     0.00000354     0.99902765     0.04210165    -0.00398172     0.01123108

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  6  1    -94.3614187657  7.5607E-07  2.0549E-07  6.6691E-04  1.0000E-04
 mr-sdci #  6  2    -93.9384552043  2.5310E-08  0.0000E+00  9.4220E-01  1.0000E-04
 mr-sdci #  6  3    -93.9034999524  8.3356E-04  0.0000E+00  9.5308E-01  1.0000E-04
 mr-sdci #  6  4    -93.4329412282  1.4455E-01  0.0000E+00  6.3037E-01  1.0000E-04
 mr-sdci #  6  5    -92.3930494808  2.9713E-04  0.0000E+00  1.1995E+00  1.0000E-04
 mr-sdci #  6  6    -92.1382863697  6.9720E-02  0.0000E+00  1.5679E+00  1.0000E-04
 
 root number  1 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.002999
time for cinew                         0.055992
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration   7

 ========= Executing OUT-OF-CORE method ========


====================================================================================================
Diagonal     counts:  0x:        9055 2x:        2233 4x:         349
All internal counts: zz :        1113 yy:        8560 xx:       15526 ww:       21792
One-external counts: yz :        3142 yx:        7351 yw:        8696
Two-external counts: yy :        1938 ww:        4100 xx:        3262 xz:         403 wz:         604 wx:        4304
Three-ext.   counts: yx :        4024 yw:        4720

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================


 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7
   ht   1   -32.24770079
   ht   2     0.00000000   -31.82473796
   ht   3    -0.00000000     0.00000000   -31.78894918
   ht   4    -0.00000000     0.00000000     0.00000000   -31.17467752
   ht   5    -0.00000000    -0.00000000     0.00000000    -0.00000000   -30.27903513
   ht   6    -0.00036590     0.00027138    -0.00006530    -0.00253427     0.00365340    -0.00001798
   ht   7    -0.00117114    -0.00015435    -0.00003441    -0.00191726     0.00071737    -0.00000002    -0.00000397

          calcsovref: tciref block   1

              civs   1       civs   2       civs   3       civs   4       civs   5       civs   6       civs   7
 refs   1  -0.947079       3.046137E-04  -9.936734E-05   0.100260      -6.182180E-03  -1.124322E-05  -4.644504E-05
 refs   2   2.324275E-04    1.00000       1.379522E-04  -4.727285E-04  -2.268475E-04  -8.547803E-06   4.900642E-06
 refs   3   1.015110E-03  -1.293862E-04   0.999785       1.951594E-02   4.396250E-03   3.490874E-06   3.460732E-06

          calcsovref: scrb block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 civs   1  -0.999950       2.044146E-04  -7.663269E-03   4.865267E-02   7.350561E-02   6.500148E-03  -5.624819E-02
 civs   2  -2.425954E-06    1.00000       3.104885E-03  -1.981642E-03   2.318514E-03  -3.433840E-04  -1.667782E-02
 civs   3   9.054533E-05   2.961438E-03  -0.994581      -9.961779E-02  -3.002144E-02  -2.513550E-03  -7.519065E-04
 civs   4  -2.195841E-04  -8.443870E-04   4.327716E-02  -0.673808       0.758002       6.448644E-02   6.772911E-02
 civs   5  -1.377008E-04  -2.620756E-04   1.223457E-02  -8.105070E-02  -2.030704E-02   -1.01111      -3.481016E-02
 civs   6   -1.06218       -1.81904        90.1716       -638.615       -659.496       -114.709       -939.230    
 civs   7   -1.04028       -4.95360        177.919       -1119.13       -1807.92       -141.773        1849.96    

          calcsovref: sovref block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1   0.947071       2.782072E-04   2.343551E-03  -5.396524E-02   9.789430E-02   1.443460E-02  -1.508944E-02
 ref    2  -2.307130E-04   0.999993       3.043812E-03  -1.672864E-03  -1.244977E-03   1.423908E-04   3.792399E-04
 ref    3  -9.367319E-04   2.790499E-03  -0.992547      -0.119155      -2.379575E-02  -6.584023E-03   3.485570E-03

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 ref:   1     0.94707061     0.00027821     0.00234355    -0.05396524     0.09789430     0.01443460    -0.01508944
 ref:   2    -0.00023071     0.99999305     0.00304381    -0.00167286    -0.00124498     0.00014239     0.00037924
 ref:   3    -0.00093673     0.00279050    -0.99254708    -0.11915539    -0.02379575    -0.00658402     0.00348557

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  7  1    -94.3614189795  2.1377E-07  4.0878E-08  3.0072E-04  1.0000E-04
 mr-sdci #  7  2    -93.9384578755  2.6712E-06  0.0000E+00  9.4220E-01  1.0000E-04
 mr-sdci #  7  3    -93.9071150126  3.6151E-03  0.0000E+00  9.4054E-01  1.0000E-04
 mr-sdci #  7  4    -93.6151880196  1.8225E-01  0.0000E+00  5.9809E-01  1.0000E-04
 mr-sdci #  7  5    -92.7906533819  3.9760E-01  0.0000E+00  1.0271E+00  1.0000E-04
 mr-sdci #  7  6    -92.3891388588  2.5085E-01  0.0000E+00  1.1766E+00  1.0000E-04
 mr-sdci #  7  7    -91.5637149417 -3.7087E-02  0.0000E+00  1.7631E+00  1.0000E-04
 
 root number  1 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.004999
time for cinew                         0.065989
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration   8

 ========= Executing OUT-OF-CORE method ========


====================================================================================================
Diagonal     counts:  0x:        9055 2x:        2233 4x:         349
All internal counts: zz :        1113 yy:        8560 xx:       15526 ww:       21792
One-external counts: yz :        3142 yx:        7351 yw:        8696
Two-external counts: yy :        1938 ww:        4100 xx:        3262 xz:         403 wz:         604 wx:        4304
Three-ext.   counts: yx :        4024 yw:        4720

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================


 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7         ht   8
   ht   1   -32.24770079
   ht   2     0.00000000   -31.82473796
   ht   3    -0.00000000     0.00000000   -31.78894918
   ht   4    -0.00000000     0.00000000     0.00000000   -31.17467752
   ht   5    -0.00000000    -0.00000000     0.00000000    -0.00000000   -30.27903513
   ht   6    -0.00036590     0.00027138    -0.00006530    -0.00253427     0.00365340    -0.00001798
   ht   7    -0.00117114    -0.00015435    -0.00003441    -0.00191726     0.00071737    -0.00000002    -0.00000397
   ht   8    -0.00065341    -0.00003729    -0.00000431     0.00050894     0.00064872     0.00000010     0.00000017    -0.00000084

          calcsovref: tciref block   1

              civs   1       civs   2       civs   3       civs   4       civs   5       civs   6       civs   7       civs   8
 refs   1  -0.947079       3.046137E-04  -9.936734E-05   0.100260      -6.182180E-03  -1.124322E-05  -4.644504E-05  -2.553269E-05
 refs   2   2.324275E-04    1.00000       1.379522E-04  -4.727285E-04  -2.268475E-04  -8.547803E-06   4.900642E-06   1.169643E-06
 refs   3   1.015110E-03  -1.293862E-04   0.999785       1.951594E-02   4.396250E-03   3.490874E-06   3.460732E-06  -6.356555E-07

          calcsovref: scrb block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 civs   1  -0.999963       1.091275E-04   3.536899E-03  -1.567028E-02   2.412776E-02  -7.540756E-02  -7.517566E-02   0.129002    
 civs   2  -2.145356E-06   0.999780      -1.989568E-02   7.024211E-03   4.773135E-04  -1.462694E-03   2.853581E-03   1.853331E-02
 civs   3   1.137121E-04   1.722537E-02   0.963332       0.264128       4.371939E-02   5.115150E-03  -5.859279E-04   5.007138E-03
 civs   4  -2.221419E-04  -3.245624E-03  -8.495317E-02   0.457512      -0.830168      -0.302909      -0.218213      -2.595991E-02
 civs   5  -1.488546E-04  -1.291047E-03  -3.108480E-02   0.105621       0.134663      -0.661623       0.739184      -0.185590    
 civs   6   -1.10924       -6.84651       -172.812        581.325        55.6215        550.479        844.411        601.976    
 civs   7   -1.24939       -19.0816       -430.813        1363.90        1122.51        778.846        92.2440       -2030.46    
 civs   8    1.05117        33.3688        711.674       -2037.95       -3235.80        2004.58        3055.99       -3072.28    

          calcsovref: sovref block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8
 ref    1   0.947066      -6.715094E-06  -7.995568E-03   4.218599E-02  -7.706177E-02  -4.840821E-02  -4.705532E-02   4.235472E-02
 ref    2  -2.298216E-04   0.999788      -1.951644E-02   6.148028E-03   2.091697E-03   2.698650E-04  -4.198889E-04  -7.123546E-05
 ref    3  -9.152386E-04   1.691226E-02   0.958790       0.281493       3.426065E-02  -4.397045E-04  -3.470627E-04   8.395560E-04

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1     0.94706621    -0.00000672    -0.00799557     0.04218599    -0.07706177    -0.04840821    -0.04705532     0.04235472
 ref:   2    -0.00022982     0.99978826    -0.01951644     0.00614803     0.00209170     0.00026987    -0.00041989    -0.00007124
 ref:   3    -0.00091524     0.01691226     0.95879050     0.28149286     0.03426065    -0.00043970    -0.00034706     0.00083956

 trial vector basis is being transformed.  new dimension:   5

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  8  1    -94.3614190225  4.2970E-08  8.8880E-09  1.3561E-04  1.0000E-04
 mr-sdci #  8  2    -93.9384739833  1.6108E-05  0.0000E+00  9.4210E-01  1.0000E-04
 mr-sdci #  8  3    -93.9173036804  1.0189E-02  0.0000E+00  8.9273E-01  1.0000E-04
 mr-sdci #  8  4    -93.7381102350  1.2292E-01  0.0000E+00  6.1901E-01  1.0000E-04
 mr-sdci #  8  5    -93.1270772849  3.3642E-01  0.0000E+00  7.6477E-01  1.0000E-04
 mr-sdci #  8  6    -92.5312533380  1.4211E-01  0.0000E+00  1.2002E+00  1.0000E-04
 mr-sdci #  8  7    -92.1411667696  5.7745E-01  0.0000E+00  1.2544E+00  1.0000E-04
 mr-sdci #  8  8    -91.2399003555  6.1456E-02  0.0000E+00  1.7961E+00  1.0000E-04
 
 root number  1 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.004999
time for cinew                         0.092985
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration   9

 ========= Executing OUT-OF-CORE method ========


====================================================================================================
Diagonal     counts:  0x:        9055 2x:        2233 4x:         349
All internal counts: zz :        1113 yy:        8560 xx:       15526 ww:       21792
One-external counts: yz :        3142 yx:        7351 yw:        8696
Two-external counts: yy :        1938 ww:        4100 xx:        3262 xz:         403 wz:         604 wx:        4304
Three-ext.   counts: yx :        4024 yw:        4720

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================


 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6
   ht   1   -32.24770181
   ht   2    -0.00000000   -31.82475677
   ht   3     0.00000000     0.00000000   -31.80358646
   ht   4    -0.00000000    -0.00000000     0.00000000   -31.62439302
   ht   5    -0.00000000    -0.00000000     0.00000000     0.00000000   -31.01336007
   ht   6    -0.00065564    -0.00000979    -0.00012688     0.00044627     0.00036507    -0.00000023

                v:   1         v:   2         v:   3         v:   4         v:   5         v:   6

   eig(s)   6.774393E-09   1.000000       1.000000       1.000000       1.000000       1.000000    
 
   x:   1  -2.033124E-05  -0.526226      -0.351099      -6.596782E-03   0.149672       0.759849    
   x:   2  -2.855707E-07  -8.920400E-03  -3.386808E-03   0.998673      -4.947436E-02   1.067277E-02
   x:   3  -3.640009E-06  -5.810353E-02  -3.826700E-02  -5.099578E-02  -0.986941       0.136040    
   x:   4   1.329527E-05  -0.847643       0.184223      -2.903672E-03  -2.558144E-02  -0.496890    
   x:   5   1.060563E-05   3.363987E-02  -0.917232       3.819582E-04  -2.107156E-02  -0.396370    
   x:   6    1.00000        0.00000      -1.048716E-11    0.00000        0.00000       2.675694E-05
 bummer (warning):overlap matrix: # small eigenvalues=                      1

          calcsovref: tciref block   1

              civs   1       civs   2       civs   3       civs   4       civs   5       civs   6
 refs   1   0.947066      -6.715094E-06  -7.995568E-03   4.218599E-02  -7.706177E-02   2.322380E-05
 refs   2  -2.298216E-04   0.999788      -1.951644E-02   6.148028E-03   2.091697E-03   1.129789E-07
 refs   3  -9.152386E-04   1.691226E-02   0.958790       0.281493       3.426065E-02  -5.299626E-07

          calcsovref: scrb block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 civs   1  -0.999973       8.286079E-02   3.210872E-03  -6.794868E-02  -0.115473      -0.190378    
 civs   2  -1.798954E-06  -8.099443E-02   0.996256      -2.809976E-02  -8.987060E-03  -7.461795E-03
 civs   3  -2.788470E-05  -0.792122      -8.268506E-02  -0.579072      -0.141618      -0.110819    
 civs   4   3.667513E-05   0.392627       1.824257E-02  -0.713621       0.503504       0.329899    
 civs   5   2.436099E-05   0.130722       5.356802E-03  -0.131143      -0.703847       0.697773    
 civs   6   -1.30399       -4070.95       -157.765        3339.22        5677.02        9361.53    

          calcsovref: sovref block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1  -0.947071      -3.244698E-03   3.882060E-04  -2.171015E-03   9.909419E-02  -1.858455E-03
 ref    2   2.286902E-04  -6.330953E-02   0.997764      -2.106118E-02  -3.930025E-03  -7.082535E-04
 ref    3   9.002978E-04  -0.643767      -5.702936E-02  -0.762764      -2.121873E-02   5.604648E-03

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 ref:   1    -0.94707148    -0.00324470     0.00038821    -0.00217102     0.09909419    -0.00185845
 ref:   2     0.00022869    -0.06330953     0.99776351    -0.02106118    -0.00393003    -0.00070825
 ref:   3     0.00090030    -0.64376691    -0.05702936    -0.76276400    -0.02121873     0.00560465

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  9  1    -94.3614190341  1.1590E-08  0.0000E+00  9.2554E-05  1.0000E-04
 mr-sdci #  9  2    -93.9733710356  3.4897E-02  1.7459E-01  6.0990E-01  1.0000E-04
 mr-sdci #  9  3    -93.9383624522  2.1059E-02  0.0000E+00  9.4168E-01  1.0000E-04
 mr-sdci #  9  4    -93.8518430621  1.1373E-01  0.0000E+00  8.8325E-01  1.0000E-04
 mr-sdci #  9  5    -93.3956821868  2.6860E-01  0.0000E+00  6.7561E-01  1.0000E-04
 mr-sdci #  9  6    -92.5616077764  3.0354E-02  0.0000E+00  1.1795E+00  1.0000E-04
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.004999
time for cinew                         0.118982
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration  10

 ========= Executing OUT-OF-CORE method ========


====================================================================================================
Diagonal     counts:  0x:        9055 2x:        2233 4x:         349
All internal counts: zz :        1113 yy:        8560 xx:       15526 ww:       21792
One-external counts: yz :        3142 yx:        7351 yw:        8696
Two-external counts: yy :        1938 ww:        4100 xx:        3262 xz:         403 wz:         604 wx:        4304
Three-ext.   counts: yx :        4024 yw:        4720

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================


 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7
   ht   1   -32.24770181
   ht   2    -0.00000000   -31.82475677
   ht   3     0.00000000     0.00000000   -31.80358646
   ht   4    -0.00000000    -0.00000000     0.00000000   -31.62439302
   ht   5    -0.00000000    -0.00000000     0.00000000     0.00000000   -31.01336007
   ht   6    -0.00065564    -0.00000979    -0.00012688     0.00044627     0.00036507    -0.00000023
   ht   7     9.42484914     0.78796112     0.61949218    -2.10593927    -7.31873020    -0.00013996   -29.61433581

                v:   1         v:   2         v:   3         v:   4         v:   5         v:   6         v:   7

   eig(s)   6.514375E-09   0.586632       1.000000       1.000000       1.000000       1.000000        1.34890    
 
   x:   1  -2.562915E-05   0.520660       0.295317       0.323822      -0.445353       0.131602       0.566729    
   x:   2  -7.227700E-07   4.296640E-02  -4.703910E-02  -0.819862      -0.553728      -0.122335       4.676812E-02
   x:   3  -3.897095E-06   2.526558E-02   0.796857      -0.358168       0.409408       0.260239       2.750111E-02
   x:   4   1.452534E-05  -0.120887       0.438458       0.181056      -0.137825      -0.850911      -0.131583    
   x:   5   1.479518E-05  -0.411734       0.288701       0.248797      -0.555370       0.419453      -0.448164    
   x:   6    1.00000       3.467008E-05   1.482399E-11  -4.098720E-13   1.582276E-12   1.045250E-05   1.094394E-05
   x:   7  -1.812713E-05   0.736404       1.020018E-15    0.00000        0.00000       6.779277E-10  -0.676542    
 bummer (warning):overlap matrix: # small eigenvalues=                      1

          calcsovref: tciref block   1

              civs   1       civs   2       civs   3       civs   4       civs   5       civs   6       civs   7
 refs   1   0.947066      -6.715094E-06  -7.995568E-03   4.218599E-02  -7.706177E-02   2.322380E-05  -0.294778    
 refs   2  -2.298216E-04   0.999788      -1.951644E-02   6.148028E-03   2.091697E-03   1.129789E-07  -2.289343E-02
 refs   3  -9.152386E-04   1.691226E-02   0.958790       0.281493       3.426065E-02  -5.299626E-07   1.360774E-02

          calcsovref: scrb block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 civs   1   0.999953      -0.239491       4.005126E-03   2.634810E-02  -0.282473      -8.634771E-03  -0.266344    
 civs   2   3.405803E-06   7.800969E-02   0.993782      -1.902488E-02  -7.875898E-02  -1.638693E-02  -1.713234E-02
 civs   3   5.230280E-05   0.762709      -0.100610       0.275918      -0.560631      -2.262289E-04  -0.142478    
 civs   4  -3.374685E-05  -0.118420       3.285702E-02   0.819515       0.170005      -0.461389       0.331235    
 civs   5  -2.120692E-05  -5.707865E-02   9.473535E-03   0.160476      -0.172501       0.673284       0.767423    
 civs   6    1.46866        3677.86       -272.058       -4070.58        1725.48       -5630.85        9415.55    
 civs   7  -5.723762E-05  -0.563185      -5.241457E-03  -0.193271      -0.846402      -0.421424      -0.256172    

          calcsovref: sovref block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1   0.947073       1.791875E-02   4.737812E-04   7.390433E-03   4.700062E-02  -8.606778E-02  -2.093348E-03
 ref    2  -2.262021E-04   7.562416E-02   0.995845      -1.507304E-02  -4.747950E-02  -8.793714E-03  -3.716764E-03
 ref    3  -8.767742E-04   0.687914      -7.001376E-02   0.499915      -0.509088      -0.110047      -2.559577E-02

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 ref:   1     0.94707284     0.01791875     0.00047378     0.00739043     0.04700062    -0.08606778    -0.00209335
 ref:   2    -0.00022620     0.07562416     0.99584537    -0.01507304    -0.04747950    -0.00879371    -0.00371676
 ref:   3    -0.00087677     0.68791413    -0.07001376     0.49991515    -0.50908775    -0.11004722    -0.02559577

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 10  1    -94.3614190354  1.3527E-09  0.0000E+00  9.5897E-05  1.0000E-04
 mr-sdci # 10  2    -94.0904193510  1.1705E-01  5.0718E-01  4.5044E-01  1.0000E-04
 mr-sdci # 10  3    -93.9383887590  2.6307E-05  0.0000E+00  9.4067E-01  1.0000E-04
 mr-sdci # 10  4    -93.8692282470  1.7385E-02  0.0000E+00  6.9122E-01  1.0000E-04
 mr-sdci # 10  5    -93.6190728135  2.2339E-01  0.0000E+00  1.0439E+00  1.0000E-04
 mr-sdci # 10  6    -93.3472745337  7.8567E-01  0.0000E+00  8.1616E-01  1.0000E-04
 mr-sdci # 10  7    -92.5040210150  3.6285E-01  0.0000E+00  1.1373E+00  1.0000E-04
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.001001
time for cinew                         0.060991
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration  11

 ========= Executing OUT-OF-CORE method ========


====================================================================================================
Diagonal     counts:  0x:        9055 2x:        2233 4x:         349
All internal counts: zz :        1113 yy:        8560 xx:       15526 ww:       21792
One-external counts: yz :        3142 yx:        7351 yw:        8696
Two-external counts: yy :        1938 ww:        4100 xx:        3262 xz:         403 wz:         604 wx:        4304
Three-ext.   counts: yx :        4024 yw:        4720

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================


 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7         ht   8
   ht   1   -32.24770181
   ht   2    -0.00000000   -31.82475677
   ht   3     0.00000000     0.00000000   -31.80358646
   ht   4    -0.00000000    -0.00000000     0.00000000   -31.62439302
   ht   5    -0.00000000    -0.00000000     0.00000000     0.00000000   -31.01336007
   ht   6    -0.00065564    -0.00000979    -0.00012688     0.00044627     0.00036507    -0.00000023
   ht   7     9.42484914     0.78796112     0.61949218    -2.10593927    -7.31873020    -0.00013996   -29.61433581
   ht   8   203.16802903    -0.90083385    -8.41591072     6.75673813   -17.97835736     0.00448318   -72.48086895 -1431.94060847

                v:   1         v:   2         v:   3         v:   4         v:   5         v:   6         v:   7         v:   8

   eig(s)   6.451479E-09   0.102288       0.670744       1.000000       1.000000       1.000000        1.14869        45.8382    
 
   x:   1  -2.669910E-06   0.984608       1.291496E-02   9.374634E-04   3.653606E-02  -7.406331E-02  -6.451572E-02  -0.139273    
   x:   2  -8.617500E-07  -4.894662E-03  -6.403245E-02   0.828964       0.535322       0.114044       9.549293E-02   5.897085E-04
   x:   3  -4.896053E-06  -4.173323E-02  -6.701430E-02  -0.502776       0.810116      -0.272644       0.101624       5.796522E-03
   x:   4   1.543880E-05   3.590227E-02   0.196880      -0.213835       0.219501       0.882751      -0.294571      -4.763227E-03
   x:   5   1.301260E-05  -8.502468E-02   0.513909       0.119622       8.733085E-02  -0.357666      -0.760685       1.285316E-02
   x:   6    1.00000       2.087432E-06  -2.645099E-05  -5.466891E-12   6.691839E-12  -1.040878E-05   3.844669E-06  -3.111540E-06
   x:   7  -1.958352E-05  -1.843508E-02  -0.828762      -4.219067E-16   5.085022E-16  -6.386597E-10  -0.556982       5.083086E-02
   x:   8   3.711751E-06   0.141151       3.912077E-02   3.666367E-18  -3.241566E-18  -2.689081E-12   2.736076E-02   0.988836    
 bummer (warning):overlap matrix: # small eigenvalues=                      1

          calcsovref: tciref block   1

              civs   1       civs   2       civs   3       civs   4       civs   5       civs   6       civs   7       civs   8
 refs   1   0.947066      -6.715094E-06  -7.995568E-03   4.218599E-02  -7.706177E-02   2.322380E-05  -0.294778       -6.67464    
 refs   2  -2.298216E-04   0.999788      -1.951644E-02   6.148028E-03   2.091697E-03   1.129789E-07  -2.289343E-02   2.428703E-02
 refs   3  -9.152386E-04   1.691226E-02   0.958790       0.281493       3.426065E-02  -5.299626E-07   1.360774E-02   0.221820    

          calcsovref: scrb block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 civs   1   -1.00001       3.708253E-02   6.860364E-04  -2.191809E-02   8.367596E-02   0.240919       0.742428        2.80567    
 civs   2  -4.023393E-06   7.214844E-02   0.993959      -1.725809E-02   8.194994E-02   1.561985E-02  -2.352719E-02  -1.367236E-02
 civs   3  -6.039823E-05   0.721858      -9.869767E-02   0.297387       0.594727      -1.090598E-02  -0.191315      -0.100771    
 civs   4   3.102364E-05  -9.592556E-02   3.270848E-02   0.809740      -0.214674       0.483049       0.332564      -7.846164E-03
 civs   5   2.583776E-05  -7.709695E-02   9.705652E-03   0.163412       0.185151      -0.648155       0.705869      -0.398615    
 civs   6   -1.63414        4231.13       -274.675       -4109.45       -2045.87        5866.12        8538.61       -2938.20    
 civs   7   7.478670E-05  -0.599538      -5.736323E-03  -0.192955       0.832166       0.440461      -0.268477      -4.010405E-02
 civs   8  -9.831138E-06   4.737451E-02  -5.123120E-04  -7.801150E-03  -3.192781E-02   3.674538E-02   0.157857       0.437705    

          calcsovref: sovref block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8
 ref    1  -0.947067      -9.972555E-03   7.955291E-04   1.194285E-02  -2.854384E-02   5.971188E-02  -0.111906      -0.289587    
 ref    2   2.250886E-04   7.263961E-02   0.995984      -1.396959E-02   4.931614E-02   8.859593E-03  -5.493051E-03  -2.012897E-03
 ref    3   8.665887E-04   0.663761      -6.832733E-02   0.516217       0.522769       0.114392      -3.987307E-02  -1.717926E-02

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1    -0.94706726    -0.00997256     0.00079553     0.01194285    -0.02854384     0.05971188    -0.11190621    -0.28958700
 ref:   2     0.00022509     0.07263961     0.99598353    -0.01396959     0.04931614     0.00885959    -0.00549305    -0.00201290
 ref:   3     0.00086659     0.66376139    -0.06832733     0.51621669     0.52276852     0.11439160    -0.03987307    -0.01717926

 trial vector basis is being transformed.  new dimension:   5

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 11  1    -94.3614190365  1.0597E-09  0.0000E+00  8.9351E-05  1.0000E-04
 mr-sdci # 11  2    -94.1145935040  2.4174E-02  7.0385E-02  3.5759E-01  1.0000E-04
 mr-sdci # 11  3    -93.9383917218  2.9628E-06  0.0000E+00  9.4070E-01  1.0000E-04
 mr-sdci # 11  4    -93.8698638164  6.3557E-04  0.0000E+00  7.0294E-01  1.0000E-04
 mr-sdci # 11  5    -93.6281370737  9.0643E-03  0.0000E+00  1.0205E+00  1.0000E-04
 mr-sdci # 11  6    -93.3578653907  1.0591E-02  0.0000E+00  8.4191E-01  1.0000E-04
 mr-sdci # 11  7    -92.6509767994  1.4696E-01  0.0000E+00  1.0315E+00  1.0000E-04
 mr-sdci # 11  8    -91.3997388262  1.5984E-01  0.0000E+00  1.3547E+00  1.0000E-04
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.004999
time for cinew                         0.090986
time for eigenvalue solver             0.000999
time for vector access                 0.000000

          starting ci iteration  12

 ========= Executing OUT-OF-CORE method ========


====================================================================================================
Diagonal     counts:  0x:        9055 2x:        2233 4x:         349
All internal counts: zz :        1113 yy:        8560 xx:       15526 ww:       21792
One-external counts: yz :        3142 yx:        7351 yw:        8696
Two-external counts: yy :        1938 ww:        4100 xx:        3262 xz:         403 wz:         604 wx:        4304
Three-ext.   counts: yx :        4024 yw:        4720

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================


 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6
   ht   1   -32.24770182
   ht   2     0.00000000   -32.00087629
   ht   3    -0.00000000     0.00000000   -31.82467451
   ht   4    -0.00000000     0.00000000    -0.00000000   -31.75614660
   ht   5    -0.00000000     0.00000000    -0.00000000    -0.00000000   -31.51441986
   ht   6     9.13149911     4.55413331     0.04697075     1.87468227    -6.08881894    -7.93357347

          calcsovref: tciref block   1

              civs   1       civs   2       civs   3       civs   4       civs   5       civs   6
 refs   1  -0.947067      -9.972555E-03   7.955291E-04   1.194285E-02  -2.854384E-02   0.307063    
 refs   2   2.250886E-04   7.263961E-02   0.995984      -1.396959E-02   4.931614E-02  -5.471817E-04
 refs   3   8.665887E-04   0.663761      -6.832733E-02   0.516217       0.522769      -8.890512E-03

          calcsovref: scrb block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 civs   1  -0.999978      -0.182117       1.758110E-02  -0.122497      -0.243013       0.786369    
 civs   2   3.369319E-05   -1.05958      -1.587219E-02   7.128037E-02   2.694892E-02   0.245517    
 civs   3   1.068079E-06  -1.567036E-02   0.998808       4.199195E-02   1.802898E-02  -8.515080E-03
 civs   4   1.084652E-05  -0.123389       3.835103E-02  -0.993306       0.160275       6.175464E-02
 civs   5  -4.995999E-06   1.044907E-02   6.738147E-03  -7.274921E-02  -0.758978      -0.863834    
 civs   6   7.841564E-05  -0.643042       6.208272E-02  -0.432567      -0.858155        2.77699    

          calcsovref: sovref block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.947070      -1.619506E-02   3.631401E-03  -2.727648E-02  -1.003432E-02   0.130907    
 ref    2  -2.220132E-04  -9.002520E-02   0.993410       5.749856E-02  -1.933983E-02  -3.545278E-02
 ref    3  -8.419879E-04  -0.754913      -5.599802E-02  -0.502608      -0.289959      -0.280167    

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 ref:   1     0.94707025    -0.01619506     0.00363140    -0.02727648    -0.01003432     0.13090732
 ref:   2    -0.00022201    -0.09002520     0.99341008     0.05749856    -0.01933983    -0.03545278
 ref:   3    -0.00084199    -0.75491273    -0.05599802    -0.50260833    -0.28995853    -0.28016730

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 12  1    -94.3614190372  7.0305E-10  0.0000E+00  8.3170E-05  1.0000E-04
 mr-sdci # 12  2    -94.1614155124  4.6822E-02  5.2577E-02  2.6484E-01  1.0000E-04
 mr-sdci # 12  3    -93.9387069008  3.1518E-04  0.0000E+00  9.3774E-01  1.0000E-04
 mr-sdci # 12  4    -93.8870278975  1.7164E-02  0.0000E+00  6.7576E-01  1.0000E-04
 mr-sdci # 12  5    -93.7144851204  8.6348E-02  0.0000E+00  8.3499E-01  1.0000E-04
 mr-sdci # 12  6    -92.8599799341 -4.9789E-01  0.0000E+00  1.3934E+00  1.0000E-04
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.003998
time for cinew                         0.067991
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration  13

 ========= Executing OUT-OF-CORE method ========


====================================================================================================
Diagonal     counts:  0x:        9055 2x:        2233 4x:         349
All internal counts: zz :        1113 yy:        8560 xx:       15526 ww:       21792
One-external counts: yz :        3142 yx:        7351 yw:        8696
Two-external counts: yy :        1938 ww:        4100 xx:        3262 xz:         403 wz:         604 wx:        4304
Three-ext.   counts: yx :        4024 yw:        4720

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================


 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7
   ht   1   -32.24770182
   ht   2     0.00000000   -32.00087629
   ht   3    -0.00000000     0.00000000   -31.82467451
   ht   4    -0.00000000     0.00000000    -0.00000000   -31.75614660
   ht   5    -0.00000000     0.00000000    -0.00000000    -0.00000000   -31.51441986
   ht   6     9.13149911     4.55413331     0.04697075     1.87468227    -6.08881894    -7.93357347
   ht   7    17.60759166    -0.77712593     0.01071075     0.32009563     0.85619681    -5.50067712   -11.69288888

          calcsovref: tciref block   1

              civs   1       civs   2       civs   3       civs   4       civs   5       civs   6       civs   7
 refs   1  -0.947067      -9.972555E-03   7.955291E-04   1.194285E-02  -2.854384E-02   0.307063       0.595574    
 refs   2   2.250886E-04   7.263961E-02   0.995984      -1.396959E-02   4.931614E-02  -5.471817E-04   1.537101E-03
 refs   3   8.665887E-04   0.663761      -6.832733E-02   0.516217       0.522769      -8.890512E-03   8.825039E-03

          calcsovref: scrb block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 civs   1  -0.999989      -2.928167E-02  -1.855031E-02   5.356326E-02  -8.159421E-02  -0.358080        2.06042    
 civs   2   4.412491E-05   -1.07114      -4.871164E-02   0.132350       3.450207E-02  -0.255735      -9.940027E-02
 civs   3   1.664103E-06  -2.627153E-02   0.980923       0.190006       2.800356E-02   1.305147E-02   8.229028E-03
 civs   4   1.551124E-05  -0.195827       0.168858      -0.937833       0.281085      -2.009848E-02   8.647580E-02
 civs   5  -6.453106E-06   1.896246E-02   2.200341E-02  -9.633256E-02  -0.695715       0.941809       0.167281    
 civs   6   9.857113E-05  -0.854846       0.177927      -0.620076      -0.917403       -2.81531      -0.148692    
 civs   7  -3.093050E-05   0.389776      -0.126255       0.419697       0.326360       0.804270        3.85070    

          calcsovref: sovref block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1   0.947069       5.161416E-03  -3.363401E-04  -7.892036E-04   1.283971E-02  -7.091298E-02   0.293618    
 ref    2  -2.208599E-04  -9.924250E-02   0.971875       0.208204      -6.854038E-03   4.384586E-02   1.448125E-02
 ref    3  -8.339200E-04  -0.789351      -3.399064E-03  -0.450356      -0.186644       0.343151       0.102639    

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 ref:   1     0.94706860     0.00516142    -0.00033634    -0.00078920     0.01283971    -0.07091298     0.29361809
 ref:   2    -0.00022086    -0.09924250     0.97187539     0.20820394    -0.00685404     0.04384586     0.01448125
 ref:   3    -0.00083392    -0.78935128    -0.00339906    -0.45035563    -0.18664361     0.34315128     0.10263914

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 13  1    -94.3614190373  1.3154E-10  0.0000E+00  8.0958E-05  1.0000E-04
 mr-sdci # 13  2    -94.1820932448  2.0678E-02  1.2064E-02  1.4633E-01  1.0000E-04
 mr-sdci # 13  3    -93.9398883899  1.1815E-03  0.0000E+00  9.1830E-01  1.0000E-04
 mr-sdci # 13  4    -93.9120792581  2.5051E-02  0.0000E+00  6.2626E-01  1.0000E-04
 mr-sdci # 13  5    -93.7307332613  1.6248E-02  0.0000E+00  7.7608E-01  1.0000E-04
 mr-sdci # 13  6    -92.9167084380  5.6729E-02  0.0000E+00  1.3833E+00  1.0000E-04
 mr-sdci # 13  7    -91.6057113297 -1.0453E+00  0.0000E+00  1.4770E+00  1.0000E-04
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.003000
time for cinew                         0.064991
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration  14

 ========= Executing OUT-OF-CORE method ========


====================================================================================================
Diagonal     counts:  0x:        9055 2x:        2233 4x:         349
All internal counts: zz :        1113 yy:        8560 xx:       15526 ww:       21792
One-external counts: yz :        3142 yx:        7351 yw:        8696
Two-external counts: yy :        1938 ww:        4100 xx:        3262 xz:         403 wz:         604 wx:        4304
Three-ext.   counts: yx :        4024 yw:        4720

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================


 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7         ht   8
   ht   1   -32.24770182
   ht   2     0.00000000   -32.00087629
   ht   3    -0.00000000     0.00000000   -31.82467451
   ht   4    -0.00000000     0.00000000    -0.00000000   -31.75614660
   ht   5    -0.00000000     0.00000000    -0.00000000    -0.00000000   -31.51441986
   ht   6     9.13149911     4.55413331     0.04697075     1.87468227    -6.08881894    -7.93357347
   ht   7    17.60759166    -0.77712593     0.01071075     0.32009563     0.85619681    -5.50067712   -11.69288888
   ht   8    -3.49234731    -0.09634524     0.07432779     0.94938721    -1.07910618     0.42494672     2.07248117    -0.80159635

          calcsovref: tciref block   1

              civs   1       civs   2       civs   3       civs   4       civs   5       civs   6       civs   7       civs   8
 refs   1  -0.947067      -9.972555E-03   7.955291E-04   1.194285E-02  -2.854384E-02   0.307063       0.595574      -0.118483    
 refs   2   2.250886E-04   7.263961E-02   0.995984      -1.396959E-02   4.931614E-02  -5.471817E-04   1.537101E-03   1.469871E-04
 refs   3   8.665887E-04   0.663761      -6.832733E-02   0.516217       0.522769      -8.890512E-03   8.825039E-03   6.537002E-03

          calcsovref: scrb block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 civs   1  -0.999994      -6.889971E-02  -3.551038E-02  -5.574968E-03   8.600325E-02   0.124391      -0.192031       -2.21018    
 civs   2   3.549831E-05   -1.05504       0.262966       2.650207E-02  -3.342371E-02   0.130529      -0.307288       1.827728E-02
 civs   3   4.064385E-07  -4.580396E-02  -0.249503       0.965483      -3.006203E-02   2.496835E-02   4.885491E-02   4.722836E-03
 civs   4   9.180977E-06  -0.286906      -0.781500      -0.246278      -0.302303       0.325717       0.367975       2.381738E-02
 civs   5  -1.058416E-05  -3.276777E-02  -0.192486      -3.747016E-02   0.698643      -0.805198       0.478384      -0.113153    
 civs   6   8.126753E-05  -0.975449      -0.554423      -0.107519       0.902346        2.02126       -2.23568      -0.360867    
 civs   7  -9.003444E-06   0.597038       0.688293       0.111139      -0.331340       0.118221        2.37918       -3.28259    
 civs   8   1.090320E-04    1.09553        2.34850       0.330686      -0.105192        4.73268        7.92278        2.91481    

          calcsovref: sovref block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8
 ref    1   0.947068      -5.003619E-04  -1.160130E-02  -2.092329E-03  -1.249019E-02   3.810597E-02  -3.252761E-02  -0.314664    
 ref    2  -2.227963E-04  -0.118268      -0.226276       0.965400       5.309085E-03  -1.011043E-02   5.079069E-02  -4.798722E-03
 ref    3  -8.439315E-04  -0.841355      -0.286130      -0.191005       0.177484      -0.153739       0.325233      -4.367067E-02

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1     0.94706804    -0.00050036    -0.01160130    -0.00209233    -0.01249019     0.03810597    -0.03252761    -0.31466369
 ref:   2    -0.00022280    -0.11826844    -0.22627568     0.96540026     0.00530909    -0.01011043     0.05079069    -0.00479872
 ref:   3    -0.00084393    -0.84135541    -0.28612985    -0.19100529     0.17748360    -0.15373944     0.32523274    -0.04367067

 trial vector basis is being transformed.  new dimension:   5

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 14  1    -94.3614190375  1.5792E-10  0.0000E+00  8.1002E-05  1.0000E-04
 mr-sdci # 14  2    -94.1955120106  1.3419E-02  4.5306E-03  9.7996E-02  1.0000E-04
 mr-sdci # 14  3    -93.9987447183  5.8856E-02  0.0000E+00  4.4308E-01  1.0000E-04
 mr-sdci # 14  4    -93.9359346087  2.3855E-02  0.0000E+00  9.3391E-01  1.0000E-04
 mr-sdci # 14  5    -93.7309099801  1.7672E-04  0.0000E+00  7.7255E-01  1.0000E-04
 mr-sdci # 14  6    -93.1081315853  1.9142E-01  0.0000E+00  1.0424E+00  1.0000E-04
 mr-sdci # 14  7    -92.3805959743  7.7488E-01  0.0000E+00  1.5279E+00  1.0000E-04
 mr-sdci # 14  8    -91.5196613206  1.1992E-01  0.0000E+00  1.4196E+00  1.0000E-04
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.005999
time for cinew                         0.092987
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration  15

 ========= Executing OUT-OF-CORE method ========


====================================================================================================
Diagonal     counts:  0x:        9055 2x:        2233 4x:         349
All internal counts: zz :        1113 yy:        8560 xx:       15526 ww:       21792
One-external counts: yz :        3142 yx:        7351 yw:        8696
Two-external counts: yy :        1938 ww:        4100 xx:        3262 xz:         403 wz:         604 wx:        4304
Three-ext.   counts: yx :        4024 yw:        4720

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================


 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6
   ht   1   -32.24770182
   ht   2     0.00000000   -32.08179479
   ht   3     0.00000000     0.00000000   -31.88502750
   ht   4    -0.00000000    -0.00000000     0.00000000   -31.82221739
   ht   5    -0.00000000     0.00000000     0.00000000    -0.00000000   -31.61719276
   ht   6    -0.18051743     0.10819650     0.01827559     0.02374465     0.61308650    -0.11733044

          calcsovref: tciref block   1

              civs   1       civs   2       civs   3       civs   4       civs   5       civs   6
 refs   1   0.947068      -5.003619E-04  -1.160130E-02  -2.092329E-03  -1.249019E-02   6.857908E-03
 refs   2  -2.227963E-04  -0.118268      -0.226276       0.965400       5.309085E-03  -1.113876E-04
 refs   3  -8.439315E-04  -0.841355      -0.286130      -0.191005       0.177484       2.166661E-03

          calcsovref: scrb block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 civs   1  -0.999999      -6.227848E-03   1.764154E-02   1.451540E-04  -1.795580E-02  -9.214136E-02
 civs   2  -8.026460E-06   0.997872       8.368157E-02   3.612818E-04  -2.460910E-02   1.585705E-02
 civs   3  -1.039215E-05   7.876327E-02  -0.965041       5.793841E-03  -0.214007      -0.130232    
 civs   4  -2.989274E-07   1.542884E-03  -7.240163E-03  -0.999999      -1.311559E-03   1.053780E-02
 civs   5   2.456188E-06  -2.122924E-02   0.119398       1.764267E-03  -0.896882       0.533392    
 civs   6  -2.604135E-04    1.11092       -3.14916      -2.591613E-02    3.20632        16.4583    

          calcsovref: sovref block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1  -0.947068       5.692517E-04   4.788808E-03   1.962634E-03   1.868336E-02   2.042396E-02
 ref    2   2.258502E-04  -0.134585       0.202460      -0.966741       4.495420E-02   3.878522E-02
 ref    3   8.535857E-04  -0.863752       0.221457       0.189300      -7.003030E-02   0.152315    

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 ref:   1    -0.94706835     0.00056925     0.00478881     0.00196263     0.01868336     0.02042396
 ref:   2     0.00022585    -0.13458458     0.20245958    -0.96674057     0.04495420     0.03878522
 ref:   3     0.00085359    -0.86375210     0.22145725     0.18930014    -0.07003030     0.15231463

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 15  1    -94.3614190378  3.4522E-10  0.0000E+00  8.0258E-05  1.0000E-04
 mr-sdci # 15  2    -94.2005757080  5.0637E-03  2.8447E-03  6.8780E-02  1.0000E-04
 mr-sdci # 15  3    -94.0449247610  4.6180E-02  0.0000E+00  3.4768E-01  1.0000E-04
 mr-sdci # 15  4    -93.9359389802  4.3715E-06  0.0000E+00  9.3468E-01  1.0000E-04
 mr-sdci # 15  5    -93.7905907756  5.9681E-02  0.0000E+00  6.3930E-01  1.0000E-04
 mr-sdci # 15  6    -92.4189842141 -6.8915E-01  0.0000E+00  1.4304E+00  1.0000E-04
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.004999
time for cinew                         0.058992
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration  16

 ========= Executing OUT-OF-CORE method ========


====================================================================================================
Diagonal     counts:  0x:        9055 2x:        2233 4x:         349
All internal counts: zz :        1113 yy:        8560 xx:       15526 ww:       21792
One-external counts: yz :        3142 yx:        7351 yw:        8696
Two-external counts: yy :        1938 ww:        4100 xx:        3262 xz:         403 wz:         604 wx:        4304
Three-ext.   counts: yx :        4024 yw:        4720

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================


 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7
   ht   1   -32.24770182
   ht   2     0.00000000   -32.08179479
   ht   3     0.00000000     0.00000000   -31.88502750
   ht   4    -0.00000000    -0.00000000     0.00000000   -31.82221739
   ht   5    -0.00000000     0.00000000     0.00000000    -0.00000000   -31.61719276
   ht   6    -0.18051743     0.10819650     0.01827559     0.02374465     0.61308650    -0.11733044
   ht   7     0.67648169    -0.15857548    -0.15272601    -0.08466646    -0.40381920     0.03384715    -0.10611005

          calcsovref: tciref block   1

              civs   1       civs   2       civs   3       civs   4       civs   5       civs   6       civs   7
 refs   1   0.947068      -5.003619E-04  -1.160130E-02  -2.092329E-03  -1.249019E-02   6.857908E-03  -2.317411E-02
 refs   2  -2.227963E-04  -0.118268      -0.226276       0.965400       5.309085E-03  -1.113876E-04   1.034460E-03
 refs   3  -8.439315E-04  -0.841355      -0.286130      -0.191005       0.177484       2.166661E-03  -7.814412E-03

          calcsovref: scrb block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 civs   1   -1.00001       1.633619E-02   5.782423E-02  -9.653349E-03   3.438193E-02   0.316305       0.232204    
 civs   2  -1.390311E-05   0.978189      -0.202336       3.157069E-03  -4.350054E-02  -8.912358E-02  -5.172735E-02
 civs   3  -2.594714E-05   0.153033       0.884649       4.917489E-02  -0.346254      -0.277696       4.874068E-02
 civs   4   3.011587E-06  -8.360014E-03  -4.266059E-02   0.999793       2.322457E-02  -1.262932E-02  -2.182698E-02
 civs   5   1.426034E-05  -5.747900E-02  -0.215036       1.470742E-02  -0.875786       0.212778      -0.520621    
 civs   6  -5.416342E-04    1.88404        4.59664      -5.695046E-02    3.70794        8.85948       -13.9017    
 civs   7  -5.686409E-04    1.28265        3.98486      -0.475433        2.62908        17.4433        7.35922    

          calcsovref: sovref block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1  -0.947067      -2.861662E-03  -1.344545E-02  -1.362869E-03   1.199301E-02  -4.327848E-02  -3.995851E-02
 ref    2   2.327690E-04  -0.157579      -0.214973       0.953295       0.103564       7.930081E-02  -1.963735E-02
 ref    3   8.682865E-04  -0.881352      -0.134134      -0.201482      -3.674071E-02   7.723794E-02  -0.146482    

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 ref:   1    -0.94706687    -0.00286166    -0.01344545    -0.00136287     0.01199301    -0.04327848    -0.03995851
 ref:   2     0.00023277    -0.15757913    -0.21497345     0.95329475     0.10356391     0.07930081    -0.01963735
 ref:   3     0.00086829    -0.88135161    -0.13413395    -0.20148206    -0.03674071     0.07723794    -0.14648195

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 16  1    -94.3614190388  9.6712E-10  0.0000E+00  6.2214E-05  1.0000E-04
 mr-sdci # 16  2    -94.2042482015  3.6725E-03  1.4208E-03  4.8156E-02  1.0000E-04
 mr-sdci # 16  3    -94.0906366488  4.5712E-02  0.0000E+00  2.2503E-01  1.0000E-04
 mr-sdci # 16  4    -93.9366704325  7.3145E-04  0.0000E+00  9.2447E-01  1.0000E-04
 mr-sdci # 16  5    -93.8092848026  1.8694E-02  0.0000E+00  5.9880E-01  1.0000E-04
 mr-sdci # 16  6    -93.0133218126  5.9434E-01  0.0000E+00  1.1997E+00  1.0000E-04
 mr-sdci # 16  7    -92.3163123521 -6.4284E-02  0.0000E+00  1.5121E+00  1.0000E-04
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.005001
time for cinew                         0.065990
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration  17

 ========= Executing OUT-OF-CORE method ========


====================================================================================================
Diagonal     counts:  0x:        9055 2x:        2233 4x:         349
All internal counts: zz :        1113 yy:        8560 xx:       15526 ww:       21792
One-external counts: yz :        3142 yx:        7351 yw:        8696
Two-external counts: yy :        1938 ww:        4100 xx:        3262 xz:         403 wz:         604 wx:        4304
Three-ext.   counts: yx :        4024 yw:        4720

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================


 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7         ht   8
   ht   1   -32.24770182
   ht   2     0.00000000   -32.08179479
   ht   3     0.00000000     0.00000000   -31.88502750
   ht   4    -0.00000000    -0.00000000     0.00000000   -31.82221739
   ht   5    -0.00000000     0.00000000     0.00000000    -0.00000000   -31.61719276
   ht   6    -0.18051743     0.10819650     0.01827559     0.02374465     0.61308650    -0.11733044
   ht   7     0.67648169    -0.15857548    -0.15272601    -0.08466646    -0.40381920     0.03384715    -0.10611005
   ht   8    -1.36185854    -0.18148017    -0.24169056    -0.02823591     0.04484724    -0.00854911     0.02167312    -0.10129822

          calcsovref: tciref block   1

              civs   1       civs   2       civs   3       civs   4       civs   5       civs   6       civs   7       civs   8
 refs   1   0.947068      -5.003619E-04  -1.160130E-02  -2.092329E-03  -1.249019E-02   6.857908E-03  -2.317411E-02   4.664058E-02
 refs   2  -2.227963E-04  -0.118268      -0.226276       0.965400       5.309085E-03  -1.113876E-04   1.034460E-03  -2.336810E-03
 refs   3  -8.439315E-04  -0.841355      -0.286130      -0.191005       0.177484       2.166661E-03  -7.814412E-03  -1.392801E-02

          calcsovref: scrb block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 civs   1   -1.00001      -1.030494E-02  -4.222572E-02   5.117614E-02  -6.479419E-02   0.270854       0.413138       -1.14144    
 civs   2  -1.376319E-05   0.963622      -0.261979       1.564786E-02  -6.422582E-02   0.152608       6.839729E-02  -5.552176E-02
 civs   3  -2.618804E-05   0.171873       0.820114       0.143713      -0.425617       0.359832      -4.879830E-02  -0.109863    
 civs   4   3.506242E-06  -1.904672E-02  -9.167817E-02   0.991919       9.163724E-02  -3.567273E-02  -1.058839E-04   2.875414E-02
 civs   5   1.485421E-05  -7.178031E-02  -0.243587       4.859052E-02  -0.845247      -0.245454       0.453228       0.296827    
 civs   6  -5.508635E-04    2.10740        4.75718      -8.652821E-02    3.61421       -4.21698        15.4767        3.61194    
 civs   7  -5.886301E-04    1.66548        4.86175      -0.981809        3.49134       -14.5458        1.57212       -11.4729    
 civs   8  -5.753503E-05   0.791340        2.78336       -1.68802        2.78916       -13.0797       -11.0536        20.8507    

          calcsovref: sovref block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8
 ref    1  -0.947067       1.465790E-03   3.635620E-03  -1.246106E-02   2.793615E-02  -4.647208E-02  -5.970173E-02   0.179650    
 ref    2   2.334020E-04  -0.171984      -0.246383       0.926415       0.184587      -0.119284       3.089720E-02   2.053593E-05
 ref    3   8.691845E-04  -0.888490      -0.106379      -0.204171      -4.994711E-02   1.836906E-02   0.211731      -6.662775E-02

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1    -0.94706721     0.00146579     0.00363562    -0.01246106     0.02793615    -0.04647208    -0.05970173     0.17964952
 ref:   2     0.00023340    -0.17198442    -0.24638272     0.92641512     0.18458747    -0.11928431     0.03089720     0.00002054
 ref:   3     0.00086918    -0.88849020    -0.10637910    -0.20417147    -0.04994711     0.01836906     0.21173101    -0.06662775

 trial vector basis is being transformed.  new dimension:   5

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 17  1    -94.3614190388  7.2404E-12  0.0000E+00  6.2426E-05  1.0000E-04
 mr-sdci # 17  2    -94.2053739322  1.1257E-03  8.0305E-04  3.4850E-02  1.0000E-04
 mr-sdci # 17  3    -94.1052126593  1.4576E-02  0.0000E+00  2.0000E-01  1.0000E-04
 mr-sdci # 17  4    -93.9415479473  4.8775E-03  0.0000E+00  8.9087E-01  1.0000E-04
 mr-sdci # 17  5    -93.8221295937  1.2845E-02  0.0000E+00  5.9331E-01  1.0000E-04
 mr-sdci # 17  6    -93.3060923713  2.9277E-01  0.0000E+00  8.8990E-01  1.0000E-04
 mr-sdci # 17  7    -92.4736241315  1.5731E-01  0.0000E+00  1.4089E+00  1.0000E-04
 mr-sdci # 17  8    -91.8639710171  3.4431E-01  0.0000E+00  1.7034E+00  1.0000E-04
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.007000
time for cinew                         0.095985
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration  18

 ========= Executing OUT-OF-CORE method ========


====================================================================================================
Diagonal     counts:  0x:        9055 2x:        2233 4x:         349
All internal counts: zz :        1113 yy:        8560 xx:       15526 ww:       21792
One-external counts: yz :        3142 yx:        7351 yw:        8696
Two-external counts: yy :        1938 ww:        4100 xx:        3262 xz:         403 wz:         604 wx:        4304
Three-ext.   counts: yx :        4024 yw:        4720

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================


 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6
   ht   1   -32.24770182
   ht   2     0.00000000   -32.09165672
   ht   3     0.00000000    -0.00000000   -31.99149544
   ht   4    -0.00000000    -0.00000000     0.00000000   -31.82783073
   ht   5    -0.00000000     0.00000000     0.00000000     0.00000000   -31.70841238
   ht   6    -1.32210016     0.05151828    -0.15175640     0.00292876    -0.23334709    -0.08292347

          calcsovref: tciref block   1

              civs   1       civs   2       civs   3       civs   4       civs   5       civs   6
 refs   1  -0.947067       1.465790E-03   3.635620E-03  -1.246106E-02   2.793615E-02  -4.578841E-02
 refs   2   2.334020E-04  -0.171984      -0.246383       0.926415       0.184587      -5.927735E-04
 refs   3   8.691845E-04  -0.888490      -0.106379      -0.204171      -4.994711E-02  -1.082350E-03

          calcsovref: scrb block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 civs   1   -1.00000      -2.350611E-02   0.105163      -0.108602       0.100822       -1.38033    
 civs   2   2.040542E-07    1.00053       1.853035E-02  -4.058397E-03   1.267734E-03   4.064376E-02
 civs   3   2.921278E-07   1.826675E-02  -0.982987      -7.550936E-02   4.432773E-02  -0.224615    
 civs   4  -2.664447E-07  -8.051794E-03   5.495592E-02  -0.992175      -8.246797E-02   7.552742E-02
 civs   5  -1.481762E-08   1.287621E-03  -1.334562E-02   5.540455E-02  -0.975580      -0.324855    
 civs   6   3.011133E-05   0.573348       -2.56508        2.64894       -2.45919        33.6680    

          calcsovref: sovref block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.947067      -2.321562E-03   1.324993E-02  -4.806823E-03  -8.946443E-03  -0.245113    
 ref    2  -2.337768E-04  -0.184143       0.288998      -0.891232      -0.266138       3.807696E-02
 ref    3  -8.693754E-04  -0.889962       8.041911E-02   0.208483       6.247241E-02  -4.905254E-02

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 ref:   1     0.94706701    -0.00232156     0.01324993    -0.00480682    -0.00894644    -0.24511291
 ref:   2    -0.00023378    -0.18414252     0.28899772    -0.89123191    -0.26613770     0.03807696
 ref:   3    -0.00086938    -0.88996233     0.08041911     0.20848343     0.06247241    -0.04905254

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 18  1    -94.3614190388  1.4637E-12  0.0000E+00  6.2782E-05  1.0000E-04
 mr-sdci # 18  2    -94.2058345494  4.6062E-04  4.7085E-04  2.7750E-02  1.0000E-04
 mr-sdci # 18  3    -94.1146791500  9.4665E-03  0.0000E+00  1.9434E-01  1.0000E-04
 mr-sdci # 18  4    -93.9514404655  9.8925E-03  0.0000E+00  8.5029E-01  1.0000E-04
 mr-sdci # 18  5    -93.8311440137  9.0144E-03  0.0000E+00  6.1974E-01  1.0000E-04
 mr-sdci # 18  6    -92.3077047764 -9.9839E-01  0.0000E+00  1.6115E+00  1.0000E-04
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.004002
time for cinew                         0.057991
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration  19

 ========= Executing OUT-OF-CORE method ========


====================================================================================================
Diagonal     counts:  0x:        9055 2x:        2233 4x:         349
All internal counts: zz :        1113 yy:        8560 xx:       15526 ww:       21792
One-external counts: yz :        3142 yx:        7351 yw:        8696
Two-external counts: yy :        1938 ww:        4100 xx:        3262 xz:         403 wz:         604 wx:        4304
Three-ext.   counts: yx :        4024 yw:        4720

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================


 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7
   ht   1   -32.24770182
   ht   2     0.00000000   -32.09165672
   ht   3     0.00000000    -0.00000000   -31.99149544
   ht   4    -0.00000000    -0.00000000     0.00000000   -31.82783073
   ht   5    -0.00000000     0.00000000     0.00000000     0.00000000   -31.70841238
   ht   6    -1.32210016     0.05151828    -0.15175640     0.00292876    -0.23334709    -0.08292347
   ht   7     0.95062896    -0.05962998    -0.02232946    -0.06810408    -0.14241871     0.04020106    -0.04257606

          calcsovref: tciref block   1

              civs   1       civs   2       civs   3       civs   4       civs   5       civs   6       civs   7
 refs   1  -0.947067       1.465790E-03   3.635620E-03  -1.246106E-02   2.793615E-02  -4.578841E-02   3.264826E-02
 refs   2   2.334020E-04  -0.171984      -0.246383       0.926415       0.184587      -5.927735E-04   2.811117E-03
 refs   3   8.691845E-04  -0.888490      -0.106379      -0.204171      -4.994711E-02  -1.082350E-03  -4.886264E-03

          calcsovref: scrb block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 civs   1   0.999993      -1.273016E-02   2.499131E-02  -2.424251E-02   4.349265E-02   0.122217        1.82968    
 civs   2  -1.465314E-06   0.993636       0.103039       3.862585E-02   1.038306E-02  -5.199661E-02  -8.530582E-02
 civs   3  -6.130447E-06   6.552142E-02  -0.875407       0.392234       9.112302E-02  -0.285247       0.117282    
 civs   4   1.067370E-05  -4.515372E-02   0.274196       0.897883      -0.239828       0.255138      -7.202819E-02
 civs   5   3.834643E-06  -4.179890E-03   8.147804E-03  -5.411037E-02  -0.905780      -0.542565       7.949133E-02
 civs   6  -4.499964E-04    1.69635       -7.56717       -7.50523       -4.17274        21.9301       -23.7894    
 civs   7  -8.542952E-04    1.92724       -9.67583       -11.2600       -4.32781        34.6450        28.9818    

          calcsovref: sovref block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1  -0.947068      -5.556119E-04   7.001915E-04  -1.222484E-02  -1.339190E-02  -8.241469E-03   0.306073    
 ref    2   2.436242E-04  -0.225227       0.430779       0.691332      -0.423296       0.299856       2.972024E-02
 ref    3   8.734232E-04  -0.891642       6.766272E-04  -0.193542       0.100989      -0.141364      -4.022137E-02

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 ref:   1    -0.94706817    -0.00055561     0.00070019    -0.01222484    -0.01339190    -0.00824147     0.30607295
 ref:   2     0.00024362    -0.22522675     0.43077926     0.69133174    -0.42329561     0.29985559     0.02972024
 ref:   3     0.00087342    -0.89164220     0.00067663    -0.19354219     0.10098939    -0.14136409    -0.04022137

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 19  1    -94.3614190391  2.9598E-10  0.0000E+00  6.1906E-05  1.0000E-04
 mr-sdci # 19  2    -94.2067449302  9.1038E-04  4.5870E-04  3.0557E-02  1.0000E-04
 mr-sdci # 19  3    -94.1429237787  2.8245E-02  0.0000E+00  2.0744E-01  1.0000E-04
 mr-sdci # 19  4    -94.0038480352  5.2408E-02  0.0000E+00  5.9687E-01  1.0000E-04
 mr-sdci # 19  5    -93.8386624981  7.5185E-03  0.0000E+00  7.1480E-01  1.0000E-04
 mr-sdci # 19  6    -93.4511578546  1.1435E+00  0.0000E+00  1.0386E+00  1.0000E-04
 mr-sdci # 19  7    -91.6016595645 -8.7196E-01  0.0000E+00  1.4917E+00  1.0000E-04
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.003998
time for cinew                         0.062988
time for eigenvalue solver             0.001003
time for vector access                 0.000000

          starting ci iteration  20

 ========= Executing OUT-OF-CORE method ========


====================================================================================================
Diagonal     counts:  0x:        9055 2x:        2233 4x:         349
All internal counts: zz :        1113 yy:        8560 xx:       15526 ww:       21792
One-external counts: yz :        3142 yx:        7351 yw:        8696
Two-external counts: yy :        1938 ww:        4100 xx:        3262 xz:         403 wz:         604 wx:        4304
Three-ext.   counts: yx :        4024 yw:        4720

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================


 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7         ht   8
   ht   1   -32.24770182
   ht   2     0.00000000   -32.09165672
   ht   3     0.00000000    -0.00000000   -31.99149544
   ht   4    -0.00000000    -0.00000000     0.00000000   -31.82783073
   ht   5    -0.00000000     0.00000000     0.00000000     0.00000000   -31.70841238
   ht   6    -1.32210016     0.05151828    -0.15175640     0.00292876    -0.23334709    -0.08292347
   ht   7     0.95062896    -0.05962998    -0.02232946    -0.06810408    -0.14241871     0.04020106    -0.04257606
   ht   8    -0.39323712    -0.13187092     0.10044908    -0.02631554     0.14198863    -0.01583427     0.01110851    -0.01873720

          calcsovref: tciref block   1

              civs   1       civs   2       civs   3       civs   4       civs   5       civs   6       civs   7       civs   8
 refs   1  -0.947067       1.465790E-03   3.635620E-03  -1.246106E-02   2.793615E-02  -4.578841E-02   3.264826E-02  -1.336744E-02
 refs   2   2.334020E-04  -0.171984      -0.246383       0.926415       0.184587      -5.927735E-04   2.811117E-03   3.822641E-04
 refs   3   8.691845E-04  -0.888490      -0.106379      -0.204171      -4.994711E-02  -1.082350E-03  -4.886264E-03  -6.198321E-03

          calcsovref: scrb block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 civs   1  -0.999997       2.232410E-02   5.971570E-02   1.699681E-02   3.383016E-02  -1.603425E-02  -7.796595E-02   -1.91922    
 civs   2  -5.112039E-07  -0.976638       0.208373       7.653120E-02   2.072011E-03   0.113517      -0.185242       3.143703E-02
 civs   3   1.062506E-05  -0.107579      -0.773394       0.565546       9.457770E-02   0.212691       0.266015      -5.835605E-02
 civs   4  -1.898406E-05   8.944641E-02   0.433183       0.787278      -0.181188      -0.398274       1.826938E-02   7.440738E-02
 civs   5  -3.618314E-06   5.284087E-03   1.235385E-02  -1.428292E-02  -0.913152       0.384015       0.456414       1.787480E-02
 civs   6   5.849060E-04   -2.48507       -9.19125       -6.30918       -4.21944       -16.8597       -17.8842        20.9150    
 civs   7   1.147924E-03   -3.28929       -13.5606       -11.6133       -3.94771       -30.3870       -9.18762       -30.6459    
 civs   8   5.703748E-04   -1.42804       -6.77940       -8.25769        1.86833       -15.4612        44.3115        12.9837    

          calcsovref: sovref block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8
 ref    1   0.947068       1.555227E-03   4.633144E-03  -4.022276E-03  -1.560377E-02   1.839073E-02   1.365532E-02  -0.314725    
 ref    2  -2.510880E-04   0.269998       0.523052       0.542148      -0.367945      -0.451351       6.918486E-02  -1.282903E-02
 ref    3  -8.755788E-04   0.888284      -7.364374E-02  -0.173413       8.300607E-02   0.201198      -0.100714       7.152812E-03

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1     0.94706770     0.00155523     0.00463314    -0.00402228    -0.01560377     0.01839073     0.01365532    -0.31472467
 ref:   2    -0.00025109     0.26999781     0.52305217     0.54214788    -0.36794475    -0.45135064     0.06918486    -0.01282903
 ref:   3    -0.00087558     0.88828365    -0.07364374    -0.17341270     0.08300607     0.20119769    -0.10071432     0.00715281

 trial vector basis is being transformed.  new dimension:   5

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 20  1    -94.3614190392  1.6826E-10  0.0000E+00  6.0339E-05  1.0000E-04
 mr-sdci # 20  2    -94.2074018375  6.5691E-04  4.3127E-04  2.7680E-02  1.0000E-04
 mr-sdci # 20  3    -94.1620135897  1.9090E-02  0.0000E+00  1.6397E-01  1.0000E-04
 mr-sdci # 20  4    -94.0434946023  3.9647E-02  0.0000E+00  4.0736E-01  1.0000E-04
 mr-sdci # 20  5    -93.8403833710  1.7209E-03  0.0000E+00  6.7890E-01  1.0000E-04
 mr-sdci # 20  6    -93.5605821036  1.0942E-01  0.0000E+00  1.0822E+00  1.0000E-04
 mr-sdci # 20  7    -92.5932100097  9.9155E-01  0.0000E+00  1.4290E+00  1.0000E-04
 mr-sdci # 20  8    -91.5232692074 -3.4070E-01  0.0000E+00  1.4380E+00  1.0000E-04
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.003998
time for cinew                         0.096985
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration  21

 ========= Executing OUT-OF-CORE method ========


====================================================================================================
Diagonal     counts:  0x:        9055 2x:        2233 4x:         349
All internal counts: zz :        1113 yy:        8560 xx:       15526 ww:       21792
One-external counts: yz :        3142 yx:        7351 yw:        8696
Two-external counts: yy :        1938 ww:        4100 xx:        3262 xz:         403 wz:         604 wx:        4304
Three-ext.   counts: yx :        4024 yw:        4720

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================


 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6
   ht   1   -32.24770182
   ht   2    -0.00000000   -32.09368462
   ht   3    -0.00000000    -0.00000000   -32.04829637
   ht   4    -0.00000000    -0.00000000    -0.00000000   -31.92977739
   ht   5    -0.00000000     0.00000000     0.00000000     0.00000000   -31.72666616
   ht   6     0.09977476     0.12262544    -0.09204135    -0.08920151    -0.01586822    -0.01331631

          calcsovref: tciref block   1

              civs   1       civs   2       civs   3       civs   4       civs   5       civs   6
 refs   1   0.947068       1.555227E-03   4.633144E-03  -4.022276E-03  -1.560377E-02  -3.574339E-03
 refs   2  -2.510880E-04   0.269998       0.523052       0.542148      -0.367945       2.570315E-03
 refs   3  -8.755788E-04   0.888284      -7.364374E-02  -0.173413       8.300607E-02  -4.717794E-03

          calcsovref: scrb block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 civs   1  -0.999999       4.349918E-03  -1.784239E-02   2.180103E-02  -1.635927E-02   0.153024    
 civs   2   1.108839E-06    1.00189       5.674789E-02   4.896637E-03  -1.382106E-02   0.174345    
 civs   3   1.498845E-06   6.777046E-02  -0.960361      -0.199429       5.506772E-02  -0.224394    
 civs   4   1.473291E-06   2.745727E-02  -0.143726       0.951561       0.116104      -0.279300    
 civs   5  -9.113090E-07  -1.052251E-02   4.697712E-02  -8.197590E-02   0.991511       9.343202E-02
 civs   6   1.671345E-04    1.40567       -5.76591        7.04555       -5.28712        49.4573    

          calcsovref: sovref block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1  -0.947068       1.021242E-03  -8.047743E-04  -8.000818E-03  -1.230005E-02  -3.295602E-02
 ref    2   2.537348E-04   0.328326      -0.597018       0.461163      -0.290389      -0.129014    
 ref    3   8.753333E-04   0.872703       0.177174      -0.186040       7.079290E-02  -5.880309E-03

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 ref:   1    -0.94706780     0.00102124    -0.00080477    -0.00800082    -0.01230005    -0.03295602
 ref:   2     0.00025373     0.32832571    -0.59701841     0.46116347    -0.29038931    -0.12901440
 ref:   3     0.00087533     0.87270326     0.17717406    -0.18603951     0.07079290    -0.00588031

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 21  1    -94.3614190393  1.3991E-11  0.0000E+00  6.0878E-05  1.0000E-04
 mr-sdci # 21  2    -94.2080101818  6.0834E-04  5.7546E-04  3.4125E-02  1.0000E-04
 mr-sdci # 21  3    -94.1758679405  1.3854E-02  0.0000E+00  1.6700E-01  1.0000E-04
 mr-sdci # 21  4    -94.0700254730  2.6531E-02  0.0000E+00  3.4486E-01  1.0000E-04
 mr-sdci # 21  5    -93.8539626200  1.3579E-02  0.0000E+00  6.2926E-01  1.0000E-04
 mr-sdci # 21  6    -92.8122920393 -7.4829E-01  0.0000E+00  1.3794E+00  1.0000E-04
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.003998
time for cinew                         0.057991
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration  22

 ========= Executing OUT-OF-CORE method ========


====================================================================================================
Diagonal     counts:  0x:        9055 2x:        2233 4x:         349
All internal counts: zz :        1113 yy:        8560 xx:       15526 ww:       21792
One-external counts: yz :        3142 yx:        7351 yw:        8696
Two-external counts: yy :        1938 ww:        4100 xx:        3262 xz:         403 wz:         604 wx:        4304
Three-ext.   counts: yx :        4024 yw:        4720

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================


 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7
   ht   1   -32.24770182
   ht   2    -0.00000000   -32.09368462
   ht   3    -0.00000000    -0.00000000   -32.04829637
   ht   4    -0.00000000    -0.00000000    -0.00000000   -31.92977739
   ht   5    -0.00000000     0.00000000     0.00000000     0.00000000   -31.72666616
   ht   6     0.09977476     0.12262544    -0.09204135    -0.08920151    -0.01586822    -0.01331631
   ht   7    -0.28823833    -0.09953581     0.13516177     0.10834721     0.02009428     0.00462473    -0.01827327

          calcsovref: tciref block   1

              civs   1       civs   2       civs   3       civs   4       civs   5       civs   6       civs   7
 refs   1   0.947068       1.555227E-03   4.633144E-03  -4.022276E-03  -1.560377E-02  -3.574339E-03   9.932463E-03
 refs   2  -2.510880E-04   0.269998       0.523052       0.542148      -0.367945       2.570315E-03  -5.988132E-03
 refs   3  -8.755788E-04   0.888284      -7.364374E-02  -0.173413       8.300607E-02  -4.717794E-03   6.695042E-03

          calcsovref: scrb block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 civs   1   -1.00000      -5.905187E-03   1.686640E-02   2.586109E-02  -4.309850E-02   0.102647      -0.398985    
 civs   2   1.446078E-06   0.974891       0.233953       2.016051E-02  -2.811545E-03  -6.692198E-02  -0.195234    
 civs   3   3.584580E-06   0.213212      -0.880139       0.388427      -7.694038E-02   0.116510       0.239730    
 civs   4   3.899328E-06   0.103962      -0.293543      -0.870386      -0.285844       0.241574       0.202513    
 civs   5  -2.620193E-06  -4.895237E-02   0.129394       0.189762      -0.929504      -0.275965       7.684483E-02
 civs   6   2.946413E-04    4.07612       -10.0465       -10.0752        8.76932       -43.4680       -22.3265    
 civs   7   1.640520E-04    2.07192       -5.36529       -6.38129        7.85756       -26.5311        36.9095    

          calcsovref: sovref block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1  -0.947068       3.266981E-03  -5.959522E-03  -5.068056E-04   2.117591E-02  -7.165139E-03   6.733191E-02
 ref    2   2.562065E-04   0.447186      -0.597644      -0.320780       0.121532       0.322500      -0.124109    
 ref    3   8.754144E-04   0.822833       0.345740       0.160777      -1.314418E-02  -0.105468       0.132975    

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 ref:   1    -0.94706761     0.00326698    -0.00595952    -0.00050681     0.02117591    -0.00716514     0.06733191
 ref:   2     0.00025621     0.44718577    -0.59764447    -0.32077995     0.12153215     0.32249973    -0.12410908
 ref:   3     0.00087541     0.82283311     0.34573967     0.16077733    -0.01314418    -0.10546799     0.13297469

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 22  1    -94.3614190393  1.5930E-11  0.0000E+00  6.1760E-05  1.0000E-04
 mr-sdci # 22  2    -94.2092223862  1.2122E-03  5.0173E-04  2.7495E-02  1.0000E-04
 mr-sdci # 22  3    -94.1904642328  1.4596E-02  0.0000E+00  8.1682E-02  1.0000E-04
 mr-sdci # 22  4    -94.0959763357  2.5951E-02  0.0000E+00  1.9141E-01  1.0000E-04
 mr-sdci # 22  5    -93.8939104858  3.9948E-02  0.0000E+00  4.8057E-01  1.0000E-04
 mr-sdci # 22  6    -93.1338088765  3.2152E-01  0.0000E+00  1.1589E+00  1.0000E-04
 mr-sdci # 22  7    -92.2205149200 -3.7270E-01  0.0000E+00  1.6445E+00  1.0000E-04
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.002998
time for cinew                         0.067989
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration  23

 ========= Executing OUT-OF-CORE method ========


====================================================================================================
Diagonal     counts:  0x:        9055 2x:        2233 4x:         349
All internal counts: zz :        1113 yy:        8560 xx:       15526 ww:       21792
One-external counts: yz :        3142 yx:        7351 yw:        8696
Two-external counts: yy :        1938 ww:        4100 xx:        3262 xz:         403 wz:         604 wx:        4304
Three-ext.   counts: yx :        4024 yw:        4720

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================


 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7         ht   8
   ht   1   -32.24770182
   ht   2    -0.00000000   -32.09368462
   ht   3    -0.00000000    -0.00000000   -32.04829637
   ht   4    -0.00000000    -0.00000000    -0.00000000   -31.92977739
   ht   5    -0.00000000     0.00000000     0.00000000     0.00000000   -31.72666616
   ht   6     0.09977476     0.12262544    -0.09204135    -0.08920151    -0.01586822    -0.01331631
   ht   7    -0.28823833    -0.09953581     0.13516177     0.10834721     0.02009428     0.00462473    -0.01827327
   ht   8     0.85455009    -0.05180857    -0.17467752    -0.15410353     0.25752396    -0.00306881     0.00967637    -0.04306658

          calcsovref: tciref block   1

              civs   1       civs   2       civs   3       civs   4       civs   5       civs   6       civs   7       civs   8
 refs   1   0.947068       1.555227E-03   4.633144E-03  -4.022276E-03  -1.560377E-02  -3.574339E-03   9.932463E-03  -2.953297E-02
 refs   2  -2.510880E-04   0.269998       0.523052       0.542148      -0.367945       2.570315E-03  -5.988132E-03   1.500798E-02
 refs   3  -8.755788E-04   0.888284      -7.364374E-02  -0.173413       8.300607E-02  -4.717794E-03   6.695042E-03  -2.383111E-03

          calcsovref: scrb block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 civs   1  -0.999993       1.447989E-02  -3.310846E-02  -4.639380E-02   0.135251      -4.602689E-02   0.898767      -0.761791    
 civs   2   1.133413E-06   0.954868       0.303608       2.987663E-02  -2.134911E-02  -4.741457E-02  -0.150075      -0.154730    
 civs   3   3.479392E-06   0.258850      -0.833944       0.449047      -0.113349       0.138586      -0.104512       0.300578    
 civs   4   4.417116E-06   0.128028      -0.301948      -0.818571      -0.391916       0.262347      -8.903053E-02   0.258650    
 civs   5  -3.518361E-06  -6.429836E-02   0.144485       0.216698      -0.830119      -0.363140       0.470838      -8.209995E-02
 civs   6   3.753571E-04    4.99813       -10.4480       -10.1267        7.57915       -41.8060       -15.5051       -20.0913    
 civs   7   2.793612E-04    2.81234       -6.17422       -7.27184        8.76504       -27.7157        10.3762        34.1454    
 civs   8   3.037081E-04   0.911281       -2.11174       -3.02091        7.17519       -6.20404        39.2263       -14.8843    

          calcsovref: sovref block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8
 ref    1  -0.947069       4.172346E-05   2.597923E-03   1.128562E-02  -9.874341E-03   1.895624E-02  -0.156502       0.130462    
 ref    2   2.607515E-04   0.495952      -0.592653      -0.308388       0.102556       0.300945       0.169800      -0.193423    
 ref    3   8.746411E-04   0.794656       0.408476       0.159737      -5.850695E-03  -0.101461      -2.273950E-02   0.148282    

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1    -0.94706883     0.00004172     0.00259792     0.01128562    -0.00987434     0.01895624    -0.15650170     0.13046239
 ref:   2     0.00026075     0.49595158    -0.59265332    -0.30838790     0.10255594     0.30094467     0.16980024    -0.19342325
 ref:   3     0.00087464     0.79465635     0.40847633     0.15973742    -0.00585070    -0.10146114    -0.02273950     0.14828156

 trial vector basis is being transformed.  new dimension:   5

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 23  1    -94.3614190393  6.9583E-11  0.0000E+00  6.0221E-05  1.0000E-04
 mr-sdci # 23  2    -94.2096807600  4.5837E-04  5.1214E-04  2.6838E-02  1.0000E-04
 mr-sdci # 23  3    -94.1933690676  2.9048E-03  0.0000E+00  7.4438E-02  1.0000E-04
 mr-sdci # 23  4    -94.1015902470  5.6139E-03  0.0000E+00  1.9102E-01  1.0000E-04
 mr-sdci # 23  5    -93.9266991212  3.2789E-02  0.0000E+00  4.4857E-01  1.0000E-04
 mr-sdci # 23  6    -93.1410570285  7.2482E-03  0.0000E+00  1.1189E+00  1.0000E-04
 mr-sdci # 23  7    -92.8360686012  6.1555E-01  0.0000E+00  1.3859E+00  1.0000E-04
 mr-sdci # 23  8    -92.1346202815  6.1135E-01  0.0000E+00  1.6852E+00  1.0000E-04
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.007000
time for cinew                         0.092983
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration  24

 ========= Executing OUT-OF-CORE method ========


====================================================================================================
Diagonal     counts:  0x:        9055 2x:        2233 4x:         349
All internal counts: zz :        1113 yy:        8560 xx:       15526 ww:       21792
One-external counts: yz :        3142 yx:        7351 yw:        8696
Two-external counts: yy :        1938 ww:        4100 xx:        3262 xz:         403 wz:         604 wx:        4304
Three-ext.   counts: yx :        4024 yw:        4720

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================


 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6
   ht   1   -32.24770182
   ht   2     0.00000000   -32.09596354
   ht   3    -0.00000000    -0.00000000   -32.07965185
   ht   4     0.00000000     0.00000000    -0.00000000   -31.98787303
   ht   5     0.00000000     0.00000000    -0.00000000     0.00000000   -31.81298191
   ht   6     1.03317422     0.03249642     0.14407906     0.15430879    -0.20332932    -0.05333976

          calcsovref: tciref block   1

              civs   1       civs   2       civs   3       civs   4       civs   5       civs   6
 refs   1  -0.947069       4.172346E-05   2.597923E-03   1.128562E-02  -9.874341E-03   3.556667E-02
 refs   2   2.607515E-04   0.495952      -0.592653      -0.308388       0.102556       1.020876E-02
 refs   3   8.746411E-04   0.794656       0.408476       0.159737      -5.850695E-03  -6.924202E-03

          calcsovref: scrb block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 civs   1   -1.00000       2.184961E-02   4.901139E-02   9.519366E-02  -0.192786        1.31476    
 civs   2  -6.068279E-07   0.998944      -5.314722E-02  -1.203163E-02   6.217792E-03   3.056509E-02
 civs   3   4.518596E-07  -5.013523E-02  -0.989719       5.887541E-02  -6.046105E-02   0.212435    
 civs   4   6.262000E-07  -1.194660E-02  -3.212544E-02  -0.977904      -0.132122       0.253226    
 civs   5  -9.352517E-07   9.366005E-03   2.267499E-02   7.374101E-02  -0.946013      -0.407293    
 civs   6  -1.378037E-04   0.681995        1.52980        2.97126       -6.01735        41.0367    

          calcsovref: sovref block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.947068       3.247329E-03   4.832692E-03   3.910918E-03  -2.374218E-02   0.221804    
 ref    2  -2.630172E-04   0.536754       0.588064       0.298634      -7.883822E-02   0.188673    
 ref    3  -8.738830E-04   0.766672      -0.462325      -0.162641       6.170832E-03  -0.129101    

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 ref:   1     0.94706813     0.00324733     0.00483269     0.00391092    -0.02374218     0.22180407
 ref:   2    -0.00026302     0.53675350     0.58806432     0.29863414    -0.07883822     0.18867324
 ref:   3    -0.00087388     0.76667200    -0.46232451    -0.16264139     0.00617083    -0.12910050

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 24  1    -94.3614190394  1.8936E-11  0.0000E+00  5.9902E-05  1.0000E-04
 mr-sdci # 24  2    -94.2100306496  3.4989E-04  4.2942E-04  2.5760E-02  1.0000E-04
 mr-sdci # 24  3    -94.1953631873  1.9941E-03  0.0000E+00  6.9266E-02  1.0000E-04
 mr-sdci # 24  4    -94.1088221021  7.2319E-03  0.0000E+00  1.9116E-01  1.0000E-04
 mr-sdci # 24  5    -93.9612808079  3.4582E-02  0.0000E+00  3.9442E-01  1.0000E-04
 mr-sdci # 24  6    -92.4039048489 -7.3715E-01  0.0000E+00  1.6294E+00  1.0000E-04
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.005997
time for cinew                         0.057995
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration  25

 ========= Executing OUT-OF-CORE method ========


====================================================================================================
Diagonal     counts:  0x:        9055 2x:        2233 4x:         349
All internal counts: zz :        1113 yy:        8560 xx:       15526 ww:       21792
One-external counts: yz :        3142 yx:        7351 yw:        8696
Two-external counts: yy :        1938 ww:        4100 xx:        3262 xz:         403 wz:         604 wx:        4304
Three-ext.   counts: yx :        4024 yw:        4720

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================


 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7
   ht   1   -32.24770182
   ht   2     0.00000000   -32.09596354
   ht   3    -0.00000000    -0.00000000   -32.07965185
   ht   4     0.00000000     0.00000000    -0.00000000   -31.98787303
   ht   5     0.00000000     0.00000000    -0.00000000     0.00000000   -31.81298191
   ht   6     1.03317422     0.03249642     0.14407906     0.15430879    -0.20332932    -0.05333976
   ht   7    -0.88941378    -0.03006478     0.00080247     0.05387599    -0.10147742     0.02786818    -0.03814033

          calcsovref: tciref block   1

              civs   1       civs   2       civs   3       civs   4       civs   5       civs   6       civs   7
 refs   1  -0.947069       4.172346E-05   2.597923E-03   1.128562E-02  -9.874341E-03   3.556667E-02  -3.081331E-02
 refs   2   2.607515E-04   0.495952      -0.592653      -0.308388       0.102556       1.020876E-02   5.480085E-03
 refs   3   8.746411E-04   0.794656       0.408476       0.159737      -5.850695E-03  -6.924202E-03  -1.340049E-03

          calcsovref: scrb block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 civs   1   -1.00000      -1.681104E-02  -2.067396E-02   2.275946E-02  -2.470516E-02  -0.164931       -1.82717    
 civs   2  -4.922089E-07  -0.950912       0.299141      -7.237104E-02   3.224765E-02  -2.680359E-02  -5.442676E-02
 civs   3   2.793228E-07   0.250396       0.898523       0.305831      -0.151162       0.160794      -0.148814    
 civs   4   3.967992E-07   7.250572E-02   0.121731      -0.786053      -0.542192       0.319641      -0.114890    
 civs   5  -7.552710E-07  -4.510972E-02  -7.051796E-02   0.212521      -0.680126      -0.745703       0.154172    
 civs   6  -9.800815E-05   -3.42216       -5.02684        10.1091       -11.7909        22.9251       -30.5441    
 civs   7   6.486096E-05   -3.36575       -5.08974        10.9178       -12.8010        32.6105        30.7670    

          calcsovref: sovref block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1   0.947068      -2.091505E-04   2.040290E-03  -8.600992E-03  -1.320358E-03  -2.188145E-02  -0.307141    
 ref    2  -2.620073E-04  -0.700375      -0.508141       0.210098       1.250663E-02   0.129064      -3.124440E-02
 ref    3  -8.742630E-04  -0.613330       0.666206      -0.143999      -1.997427E-02  -0.102779       4.537444E-02

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 ref:   1     0.94706802    -0.00020915     0.00204029    -0.00860099    -0.00132036    -0.02188145    -0.30714139
 ref:   2    -0.00026201    -0.70037543    -0.50814136     0.21009817     0.01250663     0.12906386    -0.03124440
 ref:   3    -0.00087426    -0.61333036     0.66620620    -0.14399871    -0.01997427    -0.10277928     0.04537444

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 25  1    -94.3614190394  1.3642E-12  0.0000E+00  5.9971E-05  1.0000E-04
 mr-sdci # 25  2    -94.2115256547  1.4950E-03  9.8127E-04  4.2432E-02  1.0000E-04
 mr-sdci # 25  3    -94.2010809701  5.7178E-03  0.0000E+00  6.7367E-02  1.0000E-04
 mr-sdci # 25  4    -94.1382641496  2.9442E-02  0.0000E+00  1.6307E-01  1.0000E-04
 mr-sdci # 25  5    -94.0230857490  6.1805E-02  0.0000E+00  2.8968E-01  1.0000E-04
 mr-sdci # 25  6    -93.6172220329  1.2133E+00  0.0000E+00  7.2742E-01  1.0000E-04
 mr-sdci # 25  7    -91.5168119258 -1.3193E+00  0.0000E+00  1.4540E+00  1.0000E-04
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.007000
time for cinew                         0.061989
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration  26

 ========= Executing OUT-OF-CORE method ========


====================================================================================================
Diagonal     counts:  0x:        9055 2x:        2233 4x:         349
All internal counts: zz :        1113 yy:        8560 xx:       15526 ww:       21792
One-external counts: yz :        3142 yx:        7351 yw:        8696
Two-external counts: yy :        1938 ww:        4100 xx:        3262 xz:         403 wz:         604 wx:        4304
Three-ext.   counts: yx :        4024 yw:        4720

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================


 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7         ht   8
   ht   1   -32.24770182
   ht   2     0.00000000   -32.09596354
   ht   3    -0.00000000    -0.00000000   -32.07965185
   ht   4     0.00000000     0.00000000    -0.00000000   -31.98787303
   ht   5     0.00000000     0.00000000    -0.00000000     0.00000000   -31.81298191
   ht   6     1.03317422     0.03249642     0.14407906     0.15430879    -0.20332932    -0.05333976
   ht   7    -0.88941378    -0.03006478     0.00080247     0.05387599    -0.10147742     0.02786818    -0.03814033
   ht   8    -0.44547888     0.21099692     0.04996597    -0.04532793     0.08105736     0.01828197    -0.01013936    -0.03910605

          calcsovref: tciref block   1

              civs   1       civs   2       civs   3       civs   4       civs   5       civs   6       civs   7       civs   8
 refs   1  -0.947069       4.172346E-05   2.597923E-03   1.128562E-02  -9.874341E-03   3.556667E-02  -3.081331E-02  -1.536415E-02
 refs   2   2.607515E-04   0.495952      -0.592653      -0.308388       0.102556       1.020876E-02   5.480085E-03  -9.939263E-03
 refs   3   8.746411E-04   0.794656       0.408476       0.159737      -5.850695E-03  -6.924202E-03  -1.340049E-03  -4.381600E-03

          calcsovref: scrb block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 civs   1   -1.00000      -2.840310E-02  -2.702816E-02  -3.674579E-02  -3.284075E-02   7.559479E-02  -0.219159        1.84829    
 civs   2  -1.275126E-06  -0.870474       0.475242       0.118206       5.829300E-02   9.570646E-02   0.193381       3.859668E-02
 civs   3  -2.990198E-07   0.372815       0.808335      -0.415916      -0.153028      -0.108284       0.142956       0.140262    
 civs   4   1.154774E-08   0.116079       0.141098       0.682138      -0.660849      -0.299568       9.869332E-02   0.111640    
 civs   5  -9.281112E-07  -6.526730E-02  -7.437034E-02  -0.206725      -0.563419       0.819521      -0.177981      -0.149071    
 civs   6  -2.558569E-05   -5.37179       -5.77999       -11.2634       -11.6031       -15.7723        18.2271        29.5749    
 civs   7   1.875962E-04   -5.85531       -6.46121       -13.3757       -14.0298       -26.3618        14.3402       -31.4637    
 civs   8  -1.499981E-04    1.28786        1.45129        3.24227        3.47779        10.5799        29.5069       -2.38596    

          calcsovref: sovref block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8
 ref    1   0.947068      -6.355909E-04   1.262426E-03   5.197385E-03  -4.999168E-03   5.430795E-03  -3.612980E-02   0.310680    
 ref    2  -2.590487E-04  -0.794889      -0.403352      -0.146967       3.570603E-02  -0.122547      -6.617692E-02   3.999522E-02
 ref    3  -8.751897E-04  -0.481144       0.773111       0.115890      -3.457563E-02   7.742128E-02  -4.603161E-02  -4.387890E-02

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1     0.94706817    -0.00063559     0.00126243     0.00519739    -0.00499917     0.00543079    -0.03612980     0.31068017
 ref:   2    -0.00025905    -0.79488899    -0.40335199    -0.14696687     0.03570603    -0.12254672    -0.06617692     0.03999522
 ref:   3    -0.00087519    -0.48114376     0.77311150     0.11588973    -0.03457563     0.07742128    -0.04603161    -0.04387890

 trial vector basis is being transformed.  new dimension:   5

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 26  1    -94.3614190394  2.8180E-11  0.0000E+00  6.0851E-05  1.0000E-04
 mr-sdci # 26  2    -94.2128114361  1.2858E-03  8.1425E-04  3.8185E-02  1.0000E-04
 mr-sdci # 26  3    -94.2033219388  2.2410E-03  0.0000E+00  4.7298E-02  1.0000E-04
 mr-sdci # 26  4    -94.1495514011  1.1287E-02  0.0000E+00  1.2023E-01  1.0000E-04
 mr-sdci # 26  5    -94.0360150653  1.2929E-02  0.0000E+00  2.2987E-01  1.0000E-04
 mr-sdci # 26  6    -93.7530216096  1.3580E-01  0.0000E+00  5.8650E-01  1.0000E-04
 mr-sdci # 26  7    -92.6239620592  1.1072E+00  0.0000E+00  1.4544E+00  1.0000E-04
 mr-sdci # 26  8    -91.5100879522 -6.2453E-01  0.0000E+00  1.4469E+00  1.0000E-04
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.009998
time for cinew                         0.096985
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration  27

 ========= Executing OUT-OF-CORE method ========


====================================================================================================
Diagonal     counts:  0x:        9055 2x:        2233 4x:         349
All internal counts: zz :        1113 yy:        8560 xx:       15526 ww:       21792
One-external counts: yz :        3142 yx:        7351 yw:        8696
Two-external counts: yy :        1938 ww:        4100 xx:        3262 xz:         403 wz:         604 wx:        4304
Three-ext.   counts: yx :        4024 yw:        4720

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================


 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6
   ht   1   -32.24770182
   ht   2     0.00000000   -32.09909422
   ht   3    -0.00000000     0.00000000   -32.08960472
   ht   4    -0.00000000    -0.00000000    -0.00000000   -32.03583419
   ht   5    -0.00000000    -0.00000000    -0.00000000    -0.00000000   -31.92229785
   ht   6    -0.11514228    -0.01586278    -0.31801998    -0.02647725     0.04140936    -0.03050035

          calcsovref: tciref block   1

              civs   1       civs   2       civs   3       civs   4       civs   5       civs   6
 refs   1   0.947068      -6.355909E-04   1.262426E-03   5.197385E-03  -4.999168E-03   3.935752E-03
 refs   2  -2.590487E-04  -0.794889      -0.403352      -0.146967       3.570603E-02  -1.102338E-02
 refs   3  -8.751897E-04  -0.481144       0.773111       0.115890      -3.457563E-02   1.363250E-02

          calcsovref: scrb block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 civs   1   -1.00000       4.552240E-03  -4.450577E-03  -1.019935E-02   6.415490E-03  -0.120839    
 civs   2   2.692081E-06  -0.990482      -0.123768      -4.232277E-02   9.231721E-03  -4.478689E-02
 civs   3  -1.954511E-06  -0.106895       0.978025      -8.787323E-02   2.851070E-02  -0.369799    
 civs   4   5.294905E-06  -4.364559E-02   5.003701E-02   0.991205       3.769641E-02  -0.111120    
 civs   5   2.850255E-06  -1.096988E-02   1.122369E-02   3.457196E-02   -1.00026      -1.101553E-02
 civs   6   5.372979E-04   -1.27610        1.24754        2.85839       -1.79752        33.8473    

          calcsovref: sovref block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1  -0.947068      -3.885578E-04   2.212304E-03   6.485235E-03   4.227804E-03   1.781064E-02
 ref    2   2.510984E-04   0.850528      -0.316811      -0.106860      -4.028025E-02  -0.172383    
 ref    3   8.802249E-04   0.371843       0.838094       0.105079       3.204296E-02   0.184685    

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 ref:   1    -0.94706786    -0.00038856     0.00221230     0.00648524     0.00422780     0.01781064
 ref:   2     0.00025110     0.85052788    -0.31681094    -0.10686050    -0.04028025    -0.17238305
 ref:   3     0.00088022     0.37184335     0.83809434     0.10507850     0.03204296     0.18468531

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 27  1    -94.3614190397  2.5993E-10  0.0000E+00  6.0272E-05  1.0000E-04
 mr-sdci # 27  2    -94.2138598610  1.0484E-03  9.0235E-04  3.7019E-02  1.0000E-04
 mr-sdci # 27  3    -94.2045650480  1.2431E-03  0.0000E+00  3.8451E-02  1.0000E-04
 mr-sdci # 27  4    -94.1560203375  6.4689E-03  0.0000E+00  9.0798E-02  1.0000E-04
 mr-sdci # 27  5    -94.0383405361  2.3255E-03  0.0000E+00  2.1069E-01  1.0000E-04
 mr-sdci # 27  6    -93.2597861549 -4.9324E-01  0.0000E+00  1.3105E+00  1.0000E-04
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.003998
time for cinew                         0.054993
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration  28

 ========= Executing OUT-OF-CORE method ========


====================================================================================================
Diagonal     counts:  0x:        9055 2x:        2233 4x:         349
All internal counts: zz :        1113 yy:        8560 xx:       15526 ww:       21792
One-external counts: yz :        3142 yx:        7351 yw:        8696
Two-external counts: yy :        1938 ww:        4100 xx:        3262 xz:         403 wz:         604 wx:        4304
Three-ext.   counts: yx :        4024 yw:        4720

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================


 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7
   ht   1   -32.24770182
   ht   2     0.00000000   -32.09909422
   ht   3    -0.00000000     0.00000000   -32.08960472
   ht   4    -0.00000000    -0.00000000    -0.00000000   -32.03583419
   ht   5    -0.00000000    -0.00000000    -0.00000000    -0.00000000   -31.92229785
   ht   6    -0.11514228    -0.01586278    -0.31801998    -0.02647725     0.04140936    -0.03050035
   ht   7    -0.40486937     0.04319255    -0.26179560    -0.26763357    -0.09677465    -0.00111568    -0.04646973

          calcsovref: tciref block   1

              civs   1       civs   2       civs   3       civs   4       civs   5       civs   6       civs   7
 refs   1   0.947068      -6.355909E-04   1.262426E-03   5.197385E-03  -4.999168E-03   3.935752E-03   1.422720E-02
 refs   2  -2.590487E-04  -0.794889      -0.403352      -0.146967       3.570603E-02  -1.102338E-02   2.019938E-03
 refs   3  -8.751897E-04  -0.481144       0.773111       0.115890      -3.457563E-02   1.363250E-02   6.554241E-03

          calcsovref: scrb block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 civs   1  -0.999996       1.420037E-02   1.082295E-02   2.603855E-02   1.370775E-03  -0.159994      -0.359164    
 civs   2   4.317143E-06   0.939575      -0.308738      -0.108384       1.294919E-02   0.107864       6.440492E-03
 civs   3   1.572714E-06   0.264033       0.930294      -0.183141       3.314484E-02   0.157779      -0.432695    
 civs   4   1.500162E-05   0.128153       0.107917       0.982707       4.425520E-02   0.105506      -0.232655    
 civs   5   5.030895E-06   2.183671E-02   1.622361E-02   4.553146E-02   -1.00087      -1.040543E-02  -6.709183E-02
 civs   6   1.032232E-03    3.48283        2.53948        6.10621       -2.37240       -22.3106        24.6326    
 civs   7  -6.354145E-04   -2.12049       -1.58353       -3.80916       0.565216        19.0859        21.6023    

          calcsovref: sovref block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1  -0.947069      -2.719294E-03  -4.339060E-04  -7.832820E-04   5.269636E-03   3.293556E-02   6.271106E-02
 ref    2   2.402943E-04  -0.914090      -0.176299      -5.778793E-02  -3.861022E-02   0.119274      -2.660037E-02
 ref    3   8.857965E-04  -0.200278       0.903945       8.112564E-02   3.049038E-02  -9.624618E-02   0.115442    

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 ref:   1    -0.94706903    -0.00271929    -0.00043391    -0.00078328     0.00526964     0.03293556     0.06271106
 ref:   2     0.00024029    -0.91409028    -0.17629912    -0.05778793    -0.03861022     0.11927440    -0.02660037
 ref:   3     0.00088580    -0.20027779     0.90394512     0.08112564     0.03049038    -0.09624618     0.11544245

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 28  1    -94.3614190400  3.1896E-10  0.0000E+00  5.9991E-05  1.0000E-04
 mr-sdci # 28  2    -94.2158155143  1.9557E-03  8.5530E-04  3.8694E-02  1.0000E-04
 mr-sdci # 28  3    -94.2060407803  1.4757E-03  0.0000E+00  3.0236E-02  1.0000E-04
 mr-sdci # 28  4    -94.1632853246  7.2650E-03  0.0000E+00  7.8726E-02  1.0000E-04
 mr-sdci # 28  5    -94.0384077978  6.7262E-05  0.0000E+00  2.0764E-01  1.0000E-04
 mr-sdci # 28  6    -93.9595067711  6.9972E-01  0.0000E+00  5.0846E-01  1.0000E-04
 mr-sdci # 28  7    -92.4026165347 -2.2135E-01  0.0000E+00  1.5575E+00  1.0000E-04
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.004997
time for cinew                         0.065990
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration  29

 ========= Executing OUT-OF-CORE method ========


====================================================================================================
Diagonal     counts:  0x:        9055 2x:        2233 4x:         349
All internal counts: zz :        1113 yy:        8560 xx:       15526 ww:       21792
One-external counts: yz :        3142 yx:        7351 yw:        8696
Two-external counts: yy :        1938 ww:        4100 xx:        3262 xz:         403 wz:         604 wx:        4304
Three-ext.   counts: yx :        4024 yw:        4720

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================


 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7         ht   8
   ht   1   -32.24770182
   ht   2     0.00000000   -32.09909422
   ht   3    -0.00000000     0.00000000   -32.08960472
   ht   4    -0.00000000    -0.00000000    -0.00000000   -32.03583419
   ht   5    -0.00000000    -0.00000000    -0.00000000    -0.00000000   -31.92229785
   ht   6    -0.11514228    -0.01586278    -0.31801998    -0.02647725     0.04140936    -0.03050035
   ht   7    -0.40486937     0.04319255    -0.26179560    -0.26763357    -0.09677465    -0.00111568    -0.04646973
   ht   8    -1.15141747    -0.10480191    -0.05545329     0.13721045     0.11165462    -0.01225159    -0.00509939    -0.06920409

          calcsovref: tciref block   1

              civs   1       civs   2       civs   3       civs   4       civs   5       civs   6       civs   7       civs   8
 refs   1   0.947068      -6.355909E-04   1.262426E-03   5.197385E-03  -4.999168E-03   3.935752E-03   1.422720E-02   3.974717E-02
 refs   2  -2.590487E-04  -0.794889      -0.403352      -0.146967       3.570603E-02  -1.102338E-02   2.019938E-03  -8.348786E-03
 refs   3  -8.751897E-04  -0.481144       0.773111       0.115890      -3.457563E-02   1.363250E-02   6.554241E-03   2.429901E-03

          calcsovref: scrb block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 civs   1   -1.00001      -1.258528E-02   8.666361E-03  -2.452225E-02   0.118472      -5.240900E-02  -0.521799       -1.29443    
 civs   2   3.337102E-06   0.908502       0.372488      -0.142786       0.116355      -8.480738E-02  -5.417127E-02  -6.678122E-02
 civs   3   1.899037E-06   0.307319      -0.903690      -0.229405       0.122799      -0.114753       0.347116      -0.264978    
 civs   4   1.922147E-05   0.160997      -0.122774       0.967006       0.135068      -0.125057       0.240025      -6.536716E-02
 civs   5   5.572784E-06   2.300867E-02  -1.510387E-02   3.945204E-02   0.505718       0.863016       0.115098       3.962476E-02
 civs   6   1.164241E-03    4.31641       -2.85357        7.26140       -13.7663        11.8006       -28.7615        3.96323    
 civs   7  -8.203690E-04   -2.95791        2.00817       -5.22870        14.2143       -9.73874       -10.0312        21.5088    
 civs   8   4.255389E-04   0.961443      -0.663819        1.79988       -6.94047        3.71284        21.0173        28.2937    

          calcsovref: sovref block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8
 ref    1  -0.947067       1.733395E-03  -2.777732E-03   7.134863E-03  -1.735923E-02   7.739718E-04   8.643023E-02   0.219461    
 ref    2   2.349628E-04  -0.930536       0.126975      -4.030675E-02   9.456231E-02  -1.784584E-02  -6.663468E-03  -6.513973E-02
 ref    3   8.886237E-04  -0.139862      -0.918939       7.116252E-02  -7.435158E-02   1.386464E-02  -8.805072E-02   8.321476E-02

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1    -0.94706684     0.00173340    -0.00277773     0.00713486    -0.01735923     0.00077397     0.08643023     0.21946105
 ref:   2     0.00023496    -0.93053561     0.12697543    -0.04030675     0.09456231    -0.01784584    -0.00666347    -0.06513973
 ref:   3     0.00088862    -0.13986169    -0.91893933     0.07116252    -0.07435158     0.01386464    -0.08805072     0.08321476

 trial vector basis is being transformed.  new dimension:   5

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 29  1    -94.3614190402  2.2079E-10  0.0000E+00  6.1779E-05  1.0000E-04
 mr-sdci # 29  2    -94.2166405072  8.2499E-04  9.6962E-04  3.9676E-02  1.0000E-04
 mr-sdci # 29  3    -94.2064696123  4.2883E-04  0.0000E+00  2.8711E-02  1.0000E-04
 mr-sdci # 29  4    -94.1660737338  2.7884E-03  0.0000E+00  8.0720E-02  1.0000E-04
 mr-sdci # 29  5    -94.0483504034  9.9426E-03  0.0000E+00  3.7351E-01  1.0000E-04
 mr-sdci # 29  6    -94.0350995382  7.5593E-02  0.0000E+00  2.4015E-01  1.0000E-04
 mr-sdci # 29  7    -92.5920685078  1.8945E-01  0.0000E+00  1.3848E+00  1.0000E-04
 mr-sdci # 29  8    -92.0651105970  5.5502E-01  0.0000E+00  1.5411E+00  1.0000E-04
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.007999
time for cinew                         0.095982
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration  30

 ========= Executing OUT-OF-CORE method ========


====================================================================================================
Diagonal     counts:  0x:        9055 2x:        2233 4x:         349
All internal counts: zz :        1113 yy:        8560 xx:       15526 ww:       21792
One-external counts: yz :        3142 yx:        7351 yw:        8696
Two-external counts: yy :        1938 ww:        4100 xx:        3262 xz:         403 wz:         604 wx:        4304
Three-ext.   counts: yx :        4024 yw:        4720

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================


 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6
   ht   1   -32.24770182
   ht   2    -0.00000000   -32.10292329
   ht   3    -0.00000000     0.00000000   -32.09275240
   ht   4     0.00000000     0.00000000    -0.00000000   -32.05235652
   ht   5    -0.00000000    -0.00000000     0.00000000     0.00000000   -31.93463319
   ht   6    -1.40618287    -0.01870656     0.13429742     0.03133295     0.24799192    -0.09099374

          calcsovref: tciref block   1

              civs   1       civs   2       civs   3       civs   4       civs   5       civs   6
 refs   1  -0.947067       1.733395E-03  -2.777732E-03   7.134863E-03  -1.735923E-02  -4.865875E-02
 refs   2   2.349628E-04  -0.930536       0.126975      -4.030675E-02   9.456231E-02  -3.075378E-03
 refs   3   8.886237E-04  -0.139862      -0.918939       7.116252E-02  -7.435158E-02   6.402232E-03

          calcsovref: scrb block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 civs   1  -0.999984      -2.917869E-02  -1.954251E-02   5.775774E-02  -0.192549       -1.43270    
 civs   2  -2.285023E-06   0.997654      -4.430269E-02   2.740255E-02  -3.551614E-02  -3.235555E-02
 civs   3   9.936746E-08  -3.966195E-02  -0.996572      -2.945411E-02   4.379835E-02   0.146889    
 civs   4  -3.836990E-06   2.438416E-02   2.042980E-02  -0.995348      -9.696731E-02   6.799401E-03
 civs   5   7.508715E-06  -2.931365E-02  -2.115869E-02   8.562832E-02  -0.951922       0.383487    
 civs   6  -3.718643E-04   0.669075       0.448119       -1.32444        4.41544        32.8557    

          calcsovref: sovref block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.947069      -2.399863E-03  -9.234532E-05   1.286253E-03  -1.684312E-02  -0.248923    
 ref    2  -2.308118E-04  -0.939208      -8.952189E-02   2.306416E-02  -6.112157E-02  -1.663163E-02
 ref    3  -8.915931E-04  -9.491431E-02   0.927864      -6.239211E-02   5.669337E-02   5.059104E-02

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 ref:   1     0.94706941    -0.00239986    -0.00009235     0.00128625    -0.01684312    -0.24892313
 ref:   2    -0.00023081    -0.93920837    -0.08952189     0.02306416    -0.06112157    -0.01663163
 ref:   3    -0.00089159    -0.09491431     0.92786378    -0.06239211     0.05669337     0.05059104

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 30  1    -94.3614190405  2.5765E-10  0.0000E+00  5.7433E-05  1.0000E-04
 mr-sdci # 30  2    -94.2172905419  6.5003E-04  5.4549E-04  2.9487E-02  1.0000E-04
 mr-sdci # 30  3    -94.2067777279  3.0812E-04  0.0000E+00  2.0400E-02  1.0000E-04
 mr-sdci # 30  4    -94.1684902983  2.4166E-03  0.0000E+00  6.4204E-02  1.0000E-04
 mr-sdci # 30  5    -94.0872044307  3.8854E-02  0.0000E+00  2.2900E-01  1.0000E-04
 mr-sdci # 30  6    -91.9735810886 -2.0615E+00  0.0000E+00  1.5684E+00  1.0000E-04
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.003998
time for cinew                         0.064991
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration  31

 ========= Executing OUT-OF-CORE method ========


====================================================================================================
Diagonal     counts:  0x:        9055 2x:        2233 4x:         349
All internal counts: zz :        1113 yy:        8560 xx:       15526 ww:       21792
One-external counts: yz :        3142 yx:        7351 yw:        8696
Two-external counts: yy :        1938 ww:        4100 xx:        3262 xz:         403 wz:         604 wx:        4304
Three-ext.   counts: yx :        4024 yw:        4720

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================


 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7
   ht   1   -32.24770182
   ht   2    -0.00000000   -32.10292329
   ht   3    -0.00000000     0.00000000   -32.09275240
   ht   4     0.00000000     0.00000000    -0.00000000   -32.05235652
   ht   5    -0.00000000    -0.00000000     0.00000000     0.00000000   -31.93463319
   ht   6    -1.40618287    -0.01870656     0.13429742     0.03133295     0.24799192    -0.09099374
   ht   7     0.96758692    -0.01067628    -0.05230637     0.00771321     0.16546375     0.04457184    -0.04458557

          calcsovref: tciref block   1

              civs   1       civs   2       civs   3       civs   4       civs   5       civs   6       civs   7
 refs   1  -0.947067       1.733395E-03  -2.777732E-03   7.134863E-03  -1.735923E-02  -4.865875E-02   3.356424E-02
 refs   2   2.349628E-04  -0.930536       0.126975      -4.030675E-02   9.456231E-02  -3.075378E-03  -1.167679E-03
 refs   3   8.886237E-04  -0.139862      -0.918939       7.116252E-02  -7.435158E-02   6.402232E-03  -1.165207E-03

          calcsovref: scrb block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 civs   1  -0.999993      -2.050122E-02   1.280697E-02   4.138480E-02   8.007470E-02   0.253395        1.79804    
 civs   2  -4.361500E-06   0.989037       9.737994E-02   7.053815E-02   7.513555E-02  -4.667351E-02   1.516435E-02
 civs   3   1.616965E-06  -8.454848E-02   0.990139      -6.994210E-02  -7.120522E-02   4.005346E-02  -0.149212    
 civs   4  -7.859723E-06   5.309786E-02  -4.226352E-02  -0.967724       0.245565      -1.015195E-02  -2.812759E-03
 civs   5   1.163998E-05  -6.017152E-02   4.152068E-02   0.187591       0.839524       0.615562      -0.150941    
 civs   6  -7.557300E-04    1.49574      -0.946518       -3.12697       -7.07469        21.2835       -25.1345    
 civs   7  -8.653660E-04    1.49077      -0.948904       -3.16558       -7.61365        39.3768        23.3971    

          calcsovref: sovref block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1   0.947068       4.426597E-05  -1.525753E-03  -3.134338E-03   3.697968E-04   3.509172E-02   0.308498    
 ref    2  -2.259452E-04  -0.945246       4.475954E-02  -4.451574E-03   2.119823E-02  -4.239276E-03   3.182931E-03
 ref    3  -8.947482E-04  -4.456023E-02  -0.934535      -4.470072E-02  -2.637140E-02   1.383548E-02  -4.056299E-02

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 ref:   1     0.94706768     0.00004427    -0.00152575    -0.00313434     0.00036980     0.03509172     0.30849795
 ref:   2    -0.00022595    -0.94524554     0.04475954    -0.00445157     0.02119823    -0.00423928     0.00318293
 ref:   3    -0.00089475    -0.04456023    -0.93453488    -0.04470072    -0.02637140     0.01383548    -0.04056299

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 31  1    -94.3614190408  3.8924E-10  0.0000E+00  5.3824E-05  1.0000E-04
 mr-sdci # 31  2    -94.2181063105  8.1577E-04  2.5163E-04  2.0011E-02  1.0000E-04
 mr-sdci # 31  3    -94.2071250717  3.4734E-04  0.0000E+00  1.2497E-02  1.0000E-04
 mr-sdci # 31  4    -94.1718504075  3.3601E-03  0.0000E+00  5.4301E-02  1.0000E-04
 mr-sdci # 31  5    -94.1166804471  2.9476E-02  0.0000E+00  1.5986E-01  1.0000E-04
 mr-sdci # 31  6    -93.2692157369  1.2956E+00  0.0000E+00  8.8432E-01  1.0000E-04
 mr-sdci # 31  7    -91.5284246645 -1.0636E+00  0.0000E+00  1.4651E+00  1.0000E-04
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.005997
time for cinew                         0.065994
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration  32

 ========= Executing OUT-OF-CORE method ========


====================================================================================================
Diagonal     counts:  0x:        9055 2x:        2233 4x:         349
All internal counts: zz :        1113 yy:        8560 xx:       15526 ww:       21792
One-external counts: yz :        3142 yx:        7351 yw:        8696
Two-external counts: yy :        1938 ww:        4100 xx:        3262 xz:         403 wz:         604 wx:        4304
Three-ext.   counts: yx :        4024 yw:        4720

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================


 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7         ht   8
   ht   1   -32.24770182
   ht   2    -0.00000000   -32.10292329
   ht   3    -0.00000000     0.00000000   -32.09275240
   ht   4     0.00000000     0.00000000    -0.00000000   -32.05235652
   ht   5    -0.00000000    -0.00000000     0.00000000     0.00000000   -31.93463319
   ht   6    -1.40618287    -0.01870656     0.13429742     0.03133295     0.24799192    -0.09099374
   ht   7     0.96758692    -0.01067628    -0.05230637     0.00771321     0.16546375     0.04457184    -0.04458557
   ht   8    -0.19624623    -0.02342904    -0.06530262     0.06622358     0.18007884    -0.01092297     0.00409095    -0.01182197

          calcsovref: tciref block   1

              civs   1       civs   2       civs   3       civs   4       civs   5       civs   6       civs   7       civs   8
 refs   1  -0.947067       1.733395E-03  -2.777732E-03   7.134863E-03  -1.735923E-02  -4.865875E-02   3.356424E-02  -6.759765E-03
 refs   2   2.349628E-04  -0.930536       0.126975      -4.030675E-02   9.456231E-02  -3.075378E-03  -1.167679E-03  -1.492594E-03
 refs   3   8.886237E-04  -0.139862      -0.918939       7.116252E-02  -7.435158E-02   6.402232E-03  -1.165207E-03  -2.096386E-03

          calcsovref: scrb block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 civs   1  -0.999992      -2.344285E-02   1.441697E-02   4.862920E-02  -8.237391E-02  -0.128051      -0.163874       -1.83129    
 civs   2  -4.350701E-06   0.984805       0.111538       8.836518E-02  -8.903151E-02   6.223265E-02  -5.501143E-03  -1.904317E-02
 civs   3   2.175826E-06  -9.646558E-02   0.988676      -8.054997E-02   6.659510E-02   4.468084E-02  -0.142520       0.129707    
 civs   4  -8.617019E-06   6.236458E-02  -4.874743E-02  -0.952933      -0.298352      -4.904480E-02   9.401738E-02   1.824235E-02
 civs   5   1.100513E-05  -6.574467E-02   4.512145E-02   0.216003      -0.764991      -0.744843      -8.924789E-02   0.172741    
 civs   6  -7.841966E-04    1.71551       -1.07081       -3.73693        7.35046       -13.8409       -18.0622        23.7796    
 civs   7  -9.335144E-04    1.89005       -1.18510       -4.25345        9.02424       -30.9768       -22.2198       -24.8923    
 civs   8  -2.229646E-04   0.877079      -0.538392       -2.18313        5.35644       -32.5103        46.7986        7.80067    

          calcsovref: sovref block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8
 ref    1   0.947068      -2.022883E-04  -1.370961E-03  -2.398908E-03  -2.155707E-03  -1.263147E-02  -2.544830E-02  -0.314214    
 ref    2  -2.254138E-04  -0.946173       3.346284E-02  -1.389036E-02  -1.016798E-02   6.538099E-03  -1.360251E-02  -6.349192E-03
 ref    3  -8.949049E-04  -3.284298E-02  -0.935291      -3.656030E-02   1.214395E-02   1.764911E-02  -4.293882E-02   3.519211E-02

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1     0.94706776    -0.00020229    -0.00137096    -0.00239891    -0.00215571    -0.01263147    -0.02544830    -0.31421430
 ref:   2    -0.00022541    -0.94617319     0.03346284    -0.01389036    -0.01016798     0.00653810    -0.01360251    -0.00634919
 ref:   3    -0.00089490    -0.03284298    -0.93529062    -0.03656030     0.01214395     0.01764911    -0.04293882     0.03519211

 trial vector basis is being transformed.  new dimension:   5

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 32  1    -94.3614190409  1.7835E-11  0.0000E+00  5.3476E-05  1.0000E-04
 mr-sdci # 32  2    -94.2183270908  2.2078E-04  6.8646E-05  1.2358E-02  1.0000E-04
 mr-sdci # 32  3    -94.2072078365  8.2765E-05  0.0000E+00  8.0779E-03  1.0000E-04
 mr-sdci # 32  4    -94.1731046293  1.2542E-03  0.0000E+00  4.4765E-02  1.0000E-04
 mr-sdci # 32  5    -94.1252651739  8.5847E-03  0.0000E+00  1.3598E-01  1.0000E-04
 mr-sdci # 32  6    -93.5632266431  2.9401E-01  0.0000E+00  6.4684E-01  1.0000E-04
 mr-sdci # 32  7    -92.6608998597  1.1325E+00  0.0000E+00  1.4141E+00  1.0000E-04
 mr-sdci # 32  8    -91.5037546886 -5.6136E-01  0.0000E+00  1.4489E+00  1.0000E-04
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.005997
time for cinew                         0.101990
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration  33

 ========= Executing OUT-OF-CORE method ========


====================================================================================================
Diagonal     counts:  0x:        9055 2x:        2233 4x:         349
All internal counts: zz :        1113 yy:        8560 xx:       15526 ww:       21792
One-external counts: yz :        3142 yx:        7351 yw:        8696
Two-external counts: yy :        1938 ww:        4100 xx:        3262 xz:         403 wz:         604 wx:        4304
Three-ext.   counts: yx :        4024 yw:        4720

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================


 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6
   ht   1   -32.24770183
   ht   2     0.00000000   -32.10460988
   ht   3     0.00000000    -0.00000000   -32.09349062
   ht   4    -0.00000000    -0.00000000     0.00000000   -32.05938741
   ht   5    -0.00000000     0.00000000     0.00000000    -0.00000000   -32.01154796
   ht   6    -0.01513464     0.00792443     0.05207150    -0.01361657    -0.05100160    -0.00178559

          calcsovref: tciref block   1

              civs   1       civs   2       civs   3       civs   4       civs   5       civs   6
 refs   1   0.947068      -2.022883E-04  -1.370961E-03  -2.398908E-03  -2.155707E-03   5.527051E-04
 refs   2  -2.254138E-04  -0.946173       3.346284E-02  -1.389036E-02  -1.016798E-02   1.380517E-04
 refs   3  -8.949049E-04  -3.284298E-02  -0.935291      -3.656030E-02   1.214395E-02   2.079991E-03

          calcsovref: scrb block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 civs   1  -1.000000      -5.128320E-04  -3.141472E-04  -1.503872E-03   4.597919E-03  -6.482766E-02
 civs   2   3.769939E-07    1.00020      -3.978168E-03  -4.132666E-03   5.307332E-03   2.716200E-02
 civs   3   6.995416E-07  -2.279526E-03  -0.998874       9.165960E-03  -2.125815E-02   0.228355    
 civs   4  -7.301070E-07  -4.928440E-03  -3.912270E-03   -1.00021      -3.927609E-02  -3.932000E-02
 civs   5   5.725800E-07   5.709719E-03   4.117476E-03   3.685019E-02  -0.981081      -0.289106    
 civs   6   5.173480E-04    1.09332       0.669711        3.20571       -9.80032        138.137    

          calcsovref: sovref block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1  -0.947068      -8.109190E-05   1.443365E-03   2.655796E-03   1.175070E-03   1.535186E-02
 ref    2   2.251563E-04  -0.946273      -2.955613E-02   1.817847E-02   3.434168E-03   4.511910E-03
 ref    3   8.953482E-04  -2.819330E-02   0.935954       3.524772E-02  -1.115869E-02   7.083818E-02

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 ref:   1    -0.94706771    -0.00008109     0.00144336     0.00265580     0.00117507     0.01535186
 ref:   2     0.00022516    -0.94627281    -0.02955613     0.01817847     0.00343417     0.00451191
 ref:   3     0.00089535    -0.02819330     0.93595449     0.03524772    -0.01115869     0.07083818

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 33  1    -94.3614190409  1.9689E-11  0.0000E+00  5.2999E-05  1.0000E-04
 mr-sdci # 33  2    -94.2184021488  7.5058E-05  3.6277E-05  7.3478E-03  1.0000E-04
 mr-sdci # 33  3    -94.2072356225  2.7786E-05  0.0000E+00  5.4888E-03  1.0000E-04
 mr-sdci # 33  4    -94.1736973862  5.9276E-04  0.0000E+00  3.7820E-02  1.0000E-04
 mr-sdci # 33  5    -94.1314899912  6.2248E-03  0.0000E+00  1.1089E-01  1.0000E-04
 mr-sdci # 33  6    -92.9067602919 -6.5647E-01  0.0000E+00  1.3450E+00  1.0000E-04
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.003998
time for cinew                         0.055992
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration  34

 ========= Executing OUT-OF-CORE method ========


====================================================================================================
Diagonal     counts:  0x:        9055 2x:        2233 4x:         349
All internal counts: zz :        1113 yy:        8560 xx:       15526 ww:       21792
One-external counts: yz :        3142 yx:        7351 yw:        8696
Two-external counts: yy :        1938 ww:        4100 xx:        3262 xz:         403 wz:         604 wx:        4304
Three-ext.   counts: yx :        4024 yw:        4720

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================


 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7
   ht   1   -32.24770183
   ht   2     0.00000000   -32.10460988
   ht   3     0.00000000    -0.00000000   -32.09349062
   ht   4    -0.00000000    -0.00000000     0.00000000   -32.05938741
   ht   5    -0.00000000     0.00000000     0.00000000    -0.00000000   -32.01154796
   ht   6    -0.01513464     0.00792443     0.05207150    -0.01361657    -0.05100160    -0.00178559
   ht   7     0.05920404    -0.01351398    -0.03791669     0.01008482    -0.03355150     0.00019621    -0.00152565

          calcsovref: tciref block   1

              civs   1       civs   2       civs   3       civs   4       civs   5       civs   6       civs   7
 refs   1   0.947068      -2.022883E-04  -1.370961E-03  -2.398908E-03  -2.155707E-03   5.527051E-04  -2.062706E-03
 refs   2  -2.254138E-04  -0.946173       3.346284E-02  -1.389036E-02  -1.016798E-02   1.380517E-04  -4.802157E-04
 refs   3  -8.949049E-04  -3.284298E-02  -0.935291      -3.656030E-02   1.214395E-02   2.079991E-03  -1.493897E-03

          calcsovref: scrb block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 civs   1  -0.999998       1.242526E-03   7.411773E-04   3.910319E-03  -1.249749E-02   0.176034       0.219773    
 civs   2   2.717405E-07   0.999756      -6.380047E-03  -8.179813E-03   1.315854E-02  -3.742701E-02  -6.043091E-02
 civs   3   8.559309E-09  -4.704821E-03  -0.999024       1.069215E-02  -2.042035E-02   2.176990E-02  -0.273980    
 civs   4  -1.131837E-06  -7.125958E-03  -5.714758E-03  -0.997959      -6.902555E-02   2.684464E-02   6.123180E-02
 civs   5  -7.078179E-08   7.659764E-03   5.650039E-03   5.569836E-02  -0.954872      -0.398832       9.392605E-02
 civs   6   1.097079E-03    1.67711        1.02551        5.16917       -15.0897        97.4254       -98.6496    
 civs   7   1.253661E-03    1.10513       0.665650        3.45042       -10.6625        120.783        94.4897    

          calcsovref: sovref block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1  -0.947068      -3.710705E-04   1.268150E-03   1.704096E-03   4.066745E-03  -2.780400E-02  -4.124980E-02
 ref    2   2.247224E-04  -0.946379      -2.754987E-02   2.044876E-02   5.742674E-04  -4.768616E-03  -1.283931E-02
 ref    3   8.953359E-04  -2.624475E-02   0.936003       3.302421E-02  -5.852212E-03  -2.908086E-03  -8.940753E-02

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 ref:   1    -0.94706805    -0.00037107     0.00126815     0.00170410     0.00406674    -0.02780400    -0.04124980
 ref:   2     0.00022472    -0.94637852    -0.02754987     0.02044876     0.00057427    -0.00476862    -0.01283931
 ref:   3     0.00089534    -0.02624475     0.93600281     0.03302421    -0.00585221    -0.00290809    -0.08940753

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 34  1    -94.3614190409  6.4084E-11  0.0000E+00  5.0813E-05  1.0000E-04
 mr-sdci # 34  2    -94.2184422412  4.0092E-05  9.7395E-06  4.4920E-03  1.0000E-04
 mr-sdci # 34  3    -94.2072498771  1.4255E-05  0.0000E+00  4.8371E-03  1.0000E-04
 mr-sdci # 34  4    -94.1740509733  3.5359E-04  0.0000E+00  3.4850E-02  1.0000E-04
 mr-sdci # 34  5    -94.1349730334  3.4830E-03  0.0000E+00  1.0004E-01  1.0000E-04
 mr-sdci # 34  6    -93.6014149888  6.9465E-01  0.0000E+00  8.3259E-01  1.0000E-04
 mr-sdci # 34  7    -92.4837075798 -1.7719E-01  0.0000E+00  1.5450E+00  1.0000E-04
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.006004
time for cinew                         0.068985
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration  35

 ========= Executing OUT-OF-CORE method ========


====================================================================================================
Diagonal     counts:  0x:        9055 2x:        2233 4x:         349
All internal counts: zz :        1113 yy:        8560 xx:       15526 ww:       21792
One-external counts: yz :        3142 yx:        7351 yw:        8696
Two-external counts: yy :        1938 ww:        4100 xx:        3262 xz:         403 wz:         604 wx:        4304
Three-ext.   counts: yx :        4024 yw:        4720

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================


 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7         ht   8
   ht   1   -32.24770183
   ht   2     0.00000000   -32.10460988
   ht   3     0.00000000    -0.00000000   -32.09349062
   ht   4    -0.00000000    -0.00000000     0.00000000   -32.05938741
   ht   5    -0.00000000     0.00000000     0.00000000    -0.00000000   -32.01154796
   ht   6    -0.01513464     0.00792443     0.05207150    -0.01361657    -0.05100160    -0.00178559
   ht   7     0.05920404    -0.01351398    -0.03791669     0.01008482    -0.03355150     0.00019621    -0.00152565
   ht   8    -0.08964245     0.00281021    -0.00124570    -0.00455431    -0.01104892    -0.00027219     0.00006514    -0.00050000

          calcsovref: tciref block   1

              civs   1       civs   2       civs   3       civs   4       civs   5       civs   6       civs   7       civs   8
 refs   1   0.947068      -2.022883E-04  -1.370961E-03  -2.398908E-03  -2.155707E-03   5.527051E-04  -2.062706E-03   3.083137E-03
 refs   2  -2.254138E-04  -0.946173       3.346284E-02  -1.389036E-02  -1.016798E-02   1.380517E-04  -4.802157E-04   6.771937E-05
 refs   3  -8.949049E-04  -3.284298E-02  -0.935291      -3.656030E-02   1.214395E-02   2.079991E-03  -1.493897E-03  -1.437988E-05

          calcsovref: scrb block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 civs   1  -0.999994      -1.554909E-03  -1.365981E-03  -5.677284E-03   2.350318E-02   0.254175      -0.104869       -1.08823    
 civs   2   1.581843E-07   0.999691      -7.205599E-03  -9.351187E-03   1.539070E-02   2.569699E-02   5.508619E-02   4.964193E-02
 civs   3   3.289266E-07  -5.580337E-03  -0.999090       1.132644E-02  -1.996179E-02   6.248140E-03   0.279508      -3.621022E-02
 civs   4  -8.041959E-07  -7.883140E-03  -6.475155E-03  -0.997034      -8.145784E-02  -2.005704E-02  -5.598677E-02  -5.070308E-02
 civs   5   7.154060E-07   7.918762E-03   6.009563E-03   6.209622E-02  -0.932352       0.432035      -0.110472       0.105642    
 civs   6   9.536732E-04    1.84558        1.15927        5.91638       -17.4293       -67.6925        108.252       -74.0414    
 civs   7   8.903842E-04    1.42415       0.913905        4.72580       -15.2853       -100.518       -84.6179       -83.4627    
 civs   8  -1.850541E-03    1.18855       0.899370        4.16509       -15.6086       -146.395       -36.4364        348.853    

          calcsovref: sovref block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8
 ref    1  -0.947069       8.157857E-05   1.608571E-03   3.231148E-03  -1.739006E-03  -4.160807E-02   2.269479E-02   0.176101    
 ref    2   2.248563E-04  -0.946387      -2.680343E-02   2.127532E-02  -7.468446E-04   7.351231E-04   1.226843E-02   5.176510E-03
 ref    3   8.953046E-04  -2.553367E-02   0.936020       3.211090E-02  -3.394471E-03   1.053357E-02   8.966753E-02   2.010062E-03

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1    -0.94706878     0.00008158     0.00160857     0.00323115    -0.00173901    -0.04160807     0.02269479     0.17610147
 ref:   2     0.00022486    -0.94638691    -0.02680343     0.02127532    -0.00074684     0.00073512     0.01226843     0.00517651
 ref:   3     0.00089530    -0.02553367     0.93602031     0.03211090    -0.00339447     0.01053357     0.08966753     0.00201006

 trial vector basis is being transformed.  new dimension:   5

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 35  1    -94.3614190410  3.4653E-11  0.0000E+00  5.0378E-05  1.0000E-04
 mr-sdci # 35  2    -94.2184538172  1.1576E-05  1.0676E-05  3.9903E-03  1.0000E-04
 mr-sdci # 35  3    -94.2072563633  6.4862E-06  0.0000E+00  4.7514E-03  1.0000E-04
 mr-sdci # 35  4    -94.1741795840  1.2861E-04  0.0000E+00  3.5135E-02  1.0000E-04
 mr-sdci # 35  5    -94.1367263313  1.7533E-03  0.0000E+00  9.7505E-02  1.0000E-04
 mr-sdci # 35  6    -93.8623412396  2.6093E-01  0.0000E+00  5.6177E-01  1.0000E-04
 mr-sdci # 35  7    -92.4880275395  4.3200E-03  0.0000E+00  1.5128E+00  1.0000E-04
 mr-sdci # 35  8    -92.1069616609  6.0321E-01  0.0000E+00  1.3140E+00  1.0000E-04
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.008003
time for cinew                         0.100983
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration  36

 ========= Executing OUT-OF-CORE method ========


====================================================================================================
Diagonal     counts:  0x:        9055 2x:        2233 4x:         349
All internal counts: zz :        1113 yy:        8560 xx:       15526 ww:       21792
One-external counts: yz :        3142 yx:        7351 yw:        8696
Two-external counts: yy :        1938 ww:        4100 xx:        3262 xz:         403 wz:         604 wx:        4304
Three-ext.   counts: yx :        4024 yw:        4720

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================


 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6
   ht   1   -32.24770183
   ht   2     0.00000000   -32.10473660
   ht   3    -0.00000000    -0.00000000   -32.09353915
   ht   4     0.00000000     0.00000000    -0.00000000   -32.06046237
   ht   5     0.00000000     0.00000000     0.00000000     0.00000000   -32.02300912
   ht   6    -0.16838668    -0.00187926    -0.00625196    -0.00039141    -0.00264343    -0.00119359

          calcsovref: tciref block   1

              civs   1       civs   2       civs   3       civs   4       civs   5       civs   6
 refs   1  -0.947069       8.157857E-05   1.608571E-03   3.231148E-03  -1.739006E-03  -5.843971E-03
 refs   2   2.248563E-04  -0.946387      -2.680343E-02   2.127532E-02  -7.468446E-04  -8.078814E-05
 refs   3   8.953046E-04  -2.553367E-02   0.936020       3.211090E-02  -3.394471E-03   2.637683E-04

          calcsovref: scrb block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 civs   1   -1.00000      -2.297640E-03   1.847411E-03   8.269665E-03  -3.087193E-02   -1.61059    
 civs   2   1.412363E-08   0.999973       3.580014E-04   4.745720E-04  -1.124341E-03  -1.936155E-02
 civs   3  -1.193711E-07   2.494897E-04  -0.999930       7.172903E-04  -1.873807E-03  -6.113581E-02
 civs   4   1.613869E-07   3.658825E-04  -3.951054E-04  -0.999952      -6.087958E-03  -8.426769E-03
 civs   5  -6.101095E-07  -7.937684E-04   7.346788E-04   6.062366E-03   -1.00031      -7.616172E-03
 civs   6   8.575305E-04   0.440034      -0.353807       -1.58376        5.91239        308.444    

          calcsovref: sovref block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.947068      -3.109805E-04  -1.292971E-03  -1.816849E-03  -3.597219E-03  -0.277313    
 ref    2  -2.249329E-04  -0.946396       2.648280E-02  -2.161737E-02   1.247247E-03  -5.492178E-03
 ref    3  -8.951872E-04  -2.517100E-02  -0.936071      -3.188098E-02   3.006687E-03   2.294102E-02

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 ref:   1     0.94706801    -0.00031098    -0.00129297    -0.00181685    -0.00359722    -0.27731337
 ref:   2    -0.00022493    -0.94639570     0.02648280    -0.02161737     0.00124725    -0.00549218
 ref:   3    -0.00089519    -0.02517100    -0.93607090    -0.03188098     0.00300669     0.02294102

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 36  1    -94.3614190410  1.9099E-11  0.0000E+00  5.0375E-05  1.0000E-04
 mr-sdci # 36  2    -94.2184585148  4.6976E-06  2.2125E-06  1.9456E-03  1.0000E-04
 mr-sdci # 36  3    -94.2072593815  3.0182E-06  0.0000E+00  3.9160E-03  1.0000E-04
 mr-sdci # 36  4    -94.1742386702  5.9086E-05  0.0000E+00  3.2932E-02  1.0000E-04
 mr-sdci # 36  5    -94.1375561857  8.2985E-04  0.0000E+00  8.2452E-02  1.0000E-04
 mr-sdci # 36  6    -91.8819856796 -1.9804E+00  0.0000E+00  1.5914E+00  1.0000E-04
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.005997
time for cinew                         0.064995
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration  37

 ========= Executing OUT-OF-CORE method ========


====================================================================================================
Diagonal     counts:  0x:        9055 2x:        2233 4x:         349
All internal counts: zz :        1113 yy:        8560 xx:       15526 ww:       21792
One-external counts: yz :        3142 yx:        7351 yw:        8696
Two-external counts: yy :        1938 ww:        4100 xx:        3262 xz:         403 wz:         604 wx:        4304
Three-ext.   counts: yx :        4024 yw:        4720

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================


 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7
   ht   1   -32.24770183
   ht   2     0.00000000   -32.10473660
   ht   3    -0.00000000    -0.00000000   -32.09353915
   ht   4     0.00000000     0.00000000    -0.00000000   -32.06046237
   ht   5     0.00000000     0.00000000     0.00000000     0.00000000   -32.02300912
   ht   6    -0.16838668    -0.00187926    -0.00625196    -0.00039141    -0.00264343    -0.00119359
   ht   7     0.04918815    -0.00116728     0.00370726    -0.00109974    -0.00353775     0.00027389    -0.00013586

          calcsovref: tciref block   1

              civs   1       civs   2       civs   3       civs   4       civs   5       civs   6       civs   7
 refs   1  -0.947069       8.157857E-05   1.608571E-03   3.231148E-03  -1.739006E-03  -5.843971E-03   1.696498E-03
 refs   2   2.248563E-04  -0.946387      -2.680343E-02   2.127532E-02  -7.468446E-04  -8.078814E-05  -4.103276E-05
 refs   3   8.953046E-04  -2.553367E-02   0.936020       3.211090E-02  -3.394471E-03   2.637683E-04  -1.249190E-04

          calcsovref: scrb block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 civs   1   -1.00000      -1.616162E-03   1.334101E-03   6.078865E-03  -2.467461E-02   0.107485        1.82379    
 civs   2  -4.119859E-08   0.999894       6.469330E-04   9.920438E-04  -2.234472E-03  -3.670720E-02   5.120105E-03
 civs   3   2.633309E-08   5.911326E-04   -1.00002       5.962585E-04  -1.144564E-03   4.049433E-02   8.664869E-02
 civs   4   1.580771E-07   5.460495E-04  -5.941904E-04  -0.999701      -9.422690E-03  -3.392596E-02  -3.802406E-03
 civs   5  -9.475131E-07  -1.330092E-03   1.185747E-03   9.292483E-03   -1.00180      -5.828705E-02  -2.615256E-02
 civs   6   1.188509E-03   0.758687      -0.597219       -2.67949        8.95872        171.702       -260.371    
 civs   7   1.671521E-03    1.53761       -1.16979       -5.18728        14.4914        658.255        304.347    

          calcsovref: sovref block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1   0.947068      -2.079842E-04  -1.370421E-03  -2.143735E-03  -2.691686E-03   1.156553E-02   0.310838    
 ref    2  -2.249794E-04  -0.946415       2.627462E-02  -2.180006E-02   1.369141E-03  -7.881795E-03   1.727301E-03
 ref    3  -8.951691E-04  -2.494902E-02  -0.936084      -3.165341E-02   2.614402E-03   1.106406E-03  -2.412258E-02

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 ref:   1     0.94706814    -0.00020798    -0.00137042    -0.00214373    -0.00269169     0.01156553     0.31083773
 ref:   2    -0.00022498    -0.94641502     0.02627462    -0.02180006     0.00136914    -0.00788179     0.00172730
 ref:   3    -0.00089517    -0.02494902    -0.93608449    -0.03165341     0.00261440     0.00110641    -0.02412258

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 37  1    -94.3614190410  4.8530E-12  0.0000E+00  5.0417E-05  1.0000E-04
 mr-sdci # 37  2    -94.2184619168  3.4019E-06  1.3948E-06  1.4262E-03  1.0000E-04
 mr-sdci # 37  3    -94.2072613177  1.9361E-06  0.0000E+00  3.8112E-03  1.0000E-04
 mr-sdci # 37  4    -94.1742747815  3.6111E-05  0.0000E+00  3.2326E-02  1.0000E-04
 mr-sdci # 37  5    -94.1378240910  2.6791E-04  0.0000E+00  7.7575E-02  1.0000E-04
 mr-sdci # 37  6    -93.5592102991  1.6772E+00  0.0000E+00  9.2904E-01  1.0000E-04
 mr-sdci # 37  7    -91.5235904540 -9.6444E-01  0.0000E+00  1.4596E+00  1.0000E-04
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.004997
time for cinew                         0.065994
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration  38

 ========= Executing OUT-OF-CORE method ========


====================================================================================================
Diagonal     counts:  0x:        9055 2x:        2233 4x:         349
All internal counts: zz :        1113 yy:        8560 xx:       15526 ww:       21792
One-external counts: yz :        3142 yx:        7351 yw:        8696
Two-external counts: yy :        1938 ww:        4100 xx:        3262 xz:         403 wz:         604 wx:        4304
Three-ext.   counts: yx :        4024 yw:        4720

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================


 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7         ht   8
   ht   1   -32.24770183
   ht   2     0.00000000   -32.10473660
   ht   3    -0.00000000    -0.00000000   -32.09353915
   ht   4     0.00000000     0.00000000    -0.00000000   -32.06046237
   ht   5     0.00000000     0.00000000     0.00000000     0.00000000   -32.02300912
   ht   6    -0.16838668    -0.00187926    -0.00625196    -0.00039141    -0.00264343    -0.00119359
   ht   7     0.04918815    -0.00116728     0.00370726    -0.00109974    -0.00353775     0.00027389    -0.00013586
   ht   8    -0.01430888     0.00173168     0.00223742     0.00034126     0.00487330    -0.00008495     0.00001353    -0.00005919

          calcsovref: tciref block   1

              civs   1       civs   2       civs   3       civs   4       civs   5       civs   6       civs   7       civs   8
 refs   1  -0.947069       8.157857E-05   1.608571E-03   3.231148E-03  -1.739006E-03  -5.843971E-03   1.696498E-03  -5.006507E-04
 refs   2   2.248563E-04  -0.946387      -2.680343E-02   2.127532E-02  -7.468446E-04  -8.078814E-05  -4.103276E-05   6.498989E-05
 refs   3   8.953046E-04  -2.553367E-02   0.936020       3.211090E-02  -3.394471E-03   2.637683E-04  -1.249190E-04  -9.619495E-05

          calcsovref: scrb block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 civs   1   -1.00000      -1.783381E-03   1.486835E-03  -6.413174E-03   2.486479E-02   8.397104E-03  -7.729187E-02    1.85702    
 civs   2   1.240647E-07   0.999914       7.462738E-04  -1.023069E-03   2.237923E-03   9.610927E-03   5.957933E-02  -3.431684E-03
 civs   3   5.308834E-07   8.229556E-04   -1.00012      -4.219852E-04   1.005555E-03  -6.073277E-02   3.389773E-02   8.028857E-02
 civs   4   1.511565E-07   5.850309E-04  -6.493176E-04   0.999637       9.560620E-03   2.717780E-02   2.567203E-02  -6.855713E-03
 civs   5  -4.695792E-07  -1.244536E-03   1.110062E-03  -9.167148E-03    1.00167      -1.706353E-02   0.140809      -4.796215E-02
 civs   6   1.615647E-03   0.891691      -0.723876        2.99453       -9.16939       -123.471       -151.215       -245.443    
 civs   7   3.858651E-03    2.17939       -1.77851        6.68234       -15.4798       -529.594       -372.767        343.923    
 civs   8   4.407789E-03    1.01800      -0.946379        2.18570       -1.34731       -386.475        672.271       -114.495    

          calcsovref: sovref block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8
 ref    1   0.947068      -2.474113E-04  -1.333997E-03   2.061222E-03   2.740872E-03   8.660429E-03  -1.217492E-02   0.316617    
 ref    2  -2.249878E-04  -0.946409       2.615612E-02   2.187847E-02  -1.395590E-03  -2.862852E-04   1.433279E-02  -3.210029E-04
 ref    3  -8.952883E-04  -2.487472E-02  -0.936056       3.150062E-02  -2.542055E-03   1.461135E-02  -2.750433E-02  -1.984422E-02

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1     0.94706795    -0.00024741    -0.00133400     0.00206122     0.00274087     0.00866043    -0.01217492     0.31661675
 ref:   2    -0.00022499    -0.94640948     0.02615612     0.02187847    -0.00139559    -0.00028629     0.01433279    -0.00032100
 ref:   3    -0.00089529    -0.02487472    -0.93605576     0.03150062    -0.00254206     0.01461135    -0.02750433    -0.01984422

 trial vector basis is being transformed.  new dimension:   5

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 38  1    -94.3614190410  3.4085E-11  0.0000E+00  4.9989E-05  1.0000E-04
 mr-sdci # 38  2    -94.2184633367  1.4200E-06  5.3094E-07  1.0991E-03  1.0000E-04
 mr-sdci # 38  3    -94.2072625157  1.1980E-06  0.0000E+00  3.6794E-03  1.0000E-04
 mr-sdci # 38  4    -94.1742806957  5.9142E-06  0.0000E+00  3.2641E-02  1.0000E-04
 mr-sdci # 38  5    -94.1378261274  2.0364E-06  0.0000E+00  7.7470E-02  1.0000E-04
 mr-sdci # 38  6    -93.8901711857  3.3096E-01  0.0000E+00  5.4094E-01  1.0000E-04
 mr-sdci # 38  7    -92.5434974896  1.0199E+00  0.0000E+00  1.4236E+00  1.0000E-04
 mr-sdci # 38  8    -91.4976957896 -6.0927E-01  0.0000E+00  1.4416E+00  1.0000E-04
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.015999
time for cinew                         0.104988
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration  39

 ========= Executing OUT-OF-CORE method ========


====================================================================================================
Diagonal     counts:  0x:        9055 2x:        2233 4x:         349
All internal counts: zz :        1113 yy:        8560 xx:       15526 ww:       21792
One-external counts: yz :        3142 yx:        7351 yw:        8696
Two-external counts: yy :        1938 ww:        4100 xx:        3262 xz:         403 wz:         604 wx:        4304
Three-ext.   counts: yx :        4024 yw:        4720

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================


 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6
   ht   1   -32.24770183
   ht   2    -0.00000000   -32.10474612
   ht   3     0.00000000     0.00000000   -32.09354530
   ht   4     0.00000000     0.00000000     0.00000000   -32.06056348
   ht   5     0.00000000     0.00000000     0.00000000    -0.00000000   -32.02410891
   ht   6    -0.00159085    -0.00039897     0.00205586     0.00066114    -0.00017461    -0.00001208

          calcsovref: tciref block   1

              civs   1       civs   2       civs   3       civs   4       civs   5       civs   6
 refs   1   0.947068      -2.474113E-04  -1.333997E-03   2.061222E-03   2.740872E-03   5.396734E-05
 refs   2  -2.249878E-04  -0.946409       2.615612E-02   2.187847E-02  -1.395590E-03  -1.743050E-05
 refs   3  -8.952883E-04  -2.487472E-02  -0.936056       3.150062E-02  -2.542055E-03   7.575094E-05

          calcsovref: scrb block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 civs   1   -1.00000      -5.291795E-05  -4.231609E-05  -2.893463E-04  -2.017878E-04  -7.954945E-02
 civs   2   3.941790E-08   0.999986      -5.129054E-05  -1.432620E-04  -7.768182E-05  -2.067644E-02
 civs   3  -2.779151E-07   2.838701E-05  -0.999945       4.504022E-04   2.867006E-04   0.103803    
 civs   4  -1.617950E-07   9.027403E-05   9.067544E-05    1.00011      -2.300020E-04   2.975713E-02
 civs   5   6.349725E-08  -3.117533E-05  -2.818926E-05  -3.377838E-04   -1.00002      -6.350421E-03
 civs   6  -4.532571E-03    1.07231       0.857499        5.86369        4.08947        1612.47    

          calcsovref: sovref block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1  -0.947068      -2.395924E-04   1.340246E-03   2.102385E-03  -2.712171E-03   1.159270E-02
 ref    2   2.250186E-04  -0.946413      -2.611905E-02   2.192667E-02   1.400367E-03  -5.144948E-03
 ref    3   8.951986E-04  -2.481675E-02   0.936073       3.153149E-02   2.578385E-03   2.652032E-02

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 ref:   1    -0.94706798    -0.00023959     0.00134025     0.00210238    -0.00271217     0.01159270
 ref:   2     0.00022502    -0.94641259    -0.02611905     0.02192667     0.00140037    -0.00514495
 ref:   3     0.00089520    -0.02481675     0.93607341     0.03153149     0.00257839     0.02652032

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 39  1    -94.3614190410  1.1255E-11  0.0000E+00  4.9834E-05  1.0000E-04
 mr-sdci # 39  2    -94.2184639060  5.6933E-07  2.4687E-07  5.9376E-04  1.0000E-04
 mr-sdci # 39  3    -94.2072628766  3.6089E-07  0.0000E+00  3.5809E-03  1.0000E-04
 mr-sdci # 39  4    -94.1742971413  1.6446E-05  0.0000E+00  3.2319E-02  1.0000E-04
 mr-sdci # 39  5    -94.1378338971  7.7697E-06  0.0000E+00  7.7875E-02  1.0000E-04
 mr-sdci # 39  6    -92.9304003157 -9.5977E-01  0.0000E+00  1.3247E+00  1.0000E-04
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.006996
time for cinew                         0.062988
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration  40

 ========= Executing OUT-OF-CORE method ========


====================================================================================================
Diagonal     counts:  0x:        9055 2x:        2233 4x:         349
All internal counts: zz :        1113 yy:        8560 xx:       15526 ww:       21792
One-external counts: yz :        3142 yx:        7351 yw:        8696
Two-external counts: yy :        1938 ww:        4100 xx:        3262 xz:         403 wz:         604 wx:        4304
Three-ext.   counts: yx :        4024 yw:        4720

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================


 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7
   ht   1   -32.24770183
   ht   2    -0.00000000   -32.10474612
   ht   3     0.00000000     0.00000000   -32.09354530
   ht   4     0.00000000     0.00000000     0.00000000   -32.06056348
   ht   5     0.00000000     0.00000000     0.00000000    -0.00000000   -32.02410891
   ht   6    -0.00159085    -0.00039897     0.00205586     0.00066114    -0.00017461    -0.00001208
   ht   7     0.00564789     0.00060845     0.00004066    -0.00002843     0.00127986     0.00000115    -0.00001056

          calcsovref: tciref block   1

              civs   1       civs   2       civs   3       civs   4       civs   5       civs   6       civs   7
 refs   1   0.947068      -2.474113E-04  -1.333997E-03   2.061222E-03   2.740872E-03   5.396734E-05  -1.964916E-04
 refs   2  -2.249878E-04  -0.946409       2.615612E-02   2.187847E-02  -1.395590E-03  -1.743050E-05   2.333810E-05
 refs   3  -8.952883E-04  -2.487472E-02  -0.936056       3.150062E-02  -2.542055E-03   7.575094E-05  -3.002939E-06

          calcsovref: scrb block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 civs   1  -1.000000       1.052448E-04   4.521470E-05  -5.609196E-04   6.297393E-04   0.195514       0.253060    
 civs   2   7.143967E-09    1.00000      -5.650650E-05   1.024981E-04  -2.453418E-05   1.224706E-02   3.604853E-02
 civs   3  -3.369468E-07   4.830943E-05  -0.999927      -6.529713E-04   4.777847E-04   7.402341E-02  -7.409094E-02
 civs   4  -1.911068E-07   1.277322E-04   1.169701E-04   -1.00015      -3.719203E-04   1.635564E-02  -2.326074E-02
 civs   5  -6.401732E-09  -6.505095E-06  -1.585285E-05   3.015626E-04  -0.999805       5.582234E-02   4.977636E-02
 civs   6  -5.437462E-03    1.55550        1.12712       -8.54674        6.78623        1117.15       -1173.76    
 civs   7  -2.118164E-03    1.03919       0.575732       -5.61062        5.50747        1431.01        1114.27    

          calcsovref: sovref block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1  -0.947068      -2.678023E-04   1.324633E-03  -1.949902E-03  -2.861266E-03  -3.564195E-02  -4.244594E-02
 ref    2   2.250134E-04  -0.946408      -2.610436E-02  -2.197820E-02   1.433003E-03   4.506018E-03   9.774156E-03
 ref    3   8.951920E-04  -2.480128E-02   0.936076      -3.152766E-02   2.580181E-03   1.093115E-02  -2.488843E-02

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 ref:   1    -0.94706792    -0.00026780     0.00132463    -0.00194990    -0.00286127    -0.03564195    -0.04244594
 ref:   2     0.00022501    -0.94640811    -0.02610436    -0.02197820     0.00143300     0.00450602     0.00977416
 ref:   3     0.00089519    -0.02480128     0.93607598    -0.03152766     0.00258018     0.01093115    -0.02488843

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 40  1    -94.3614190411  1.3003E-12  0.0000E+00  4.9811E-05  1.0000E-04
 mr-sdci # 40  2    -94.2184641626  2.5655E-07  7.1587E-08  3.6609E-04  1.0000E-04
 mr-sdci # 40  3    -94.2072629540  7.7409E-08  0.0000E+00  3.5707E-03  1.0000E-04
 mr-sdci # 40  4    -94.1743041192  6.9779E-06  0.0000E+00  3.2422E-02  1.0000E-04
 mr-sdci # 40  5    -94.1378402204  6.3232E-06  0.0000E+00  7.7945E-02  1.0000E-04
 mr-sdci # 40  6    -93.6305839696  7.0018E-01  0.0000E+00  8.4135E-01  1.0000E-04
 mr-sdci # 40  7    -92.5058791361 -3.7618E-02  0.0000E+00  1.5363E+00  1.0000E-04
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.006004
time for cinew                         0.071983
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration  41

 ========= Executing OUT-OF-CORE method ========


====================================================================================================
Diagonal     counts:  0x:        9055 2x:        2233 4x:         349
All internal counts: zz :        1113 yy:        8560 xx:       15526 ww:       21792
One-external counts: yz :        3142 yx:        7351 yw:        8696
Two-external counts: yy :        1938 ww:        4100 xx:        3262 xz:         403 wz:         604 wx:        4304
Three-ext.   counts: yx :        4024 yw:        4720

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================


 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7         ht   8
   ht   1   -32.24770183
   ht   2    -0.00000000   -32.10474612
   ht   3     0.00000000     0.00000000   -32.09354530
   ht   4     0.00000000     0.00000000     0.00000000   -32.06056348
   ht   5     0.00000000     0.00000000     0.00000000    -0.00000000   -32.02410891
   ht   6    -0.00159085    -0.00039897     0.00205586     0.00066114    -0.00017461    -0.00001208
   ht   7     0.00564789     0.00060845     0.00004066    -0.00002843     0.00127986     0.00000115    -0.00001056
   ht   8    -0.00810229    -0.00015331    -0.00124631    -0.00011016    -0.00046573    -0.00000179     0.00000073    -0.00000395

          calcsovref: tciref block   1

              civs   1       civs   2       civs   3       civs   4       civs   5       civs   6       civs   7       civs   8
 refs   1   0.947068      -2.474113E-04  -1.333997E-03   2.061222E-03   2.740872E-03   5.396734E-05  -1.964916E-04   2.800351E-04
 refs   2  -2.249878E-04  -0.946409       2.615612E-02   2.187847E-02  -1.395590E-03  -1.743050E-05   2.333810E-05  -4.479941E-06
 refs   3  -8.952883E-04  -2.487472E-02  -0.936056       3.150062E-02  -2.542055E-03   7.575094E-05  -3.002939E-06  -4.648886E-05

          calcsovref: scrb block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 civs   1   -1.00000      -1.576535E-04   2.548620E-05   1.689937E-03  -2.647566E-03   0.267449       3.689499E-02   -1.13762    
 civs   2   8.711654E-09   0.999998       5.908860E-05   1.308792E-04  -5.050200E-05  -3.575978E-03  -2.849856E-02  -3.235047E-02
 civs   3  -2.976153E-07   8.508116E-06   0.999937      -3.466678E-04   4.100848E-05   1.665273E-02   0.127555      -0.189992    
 civs   4  -1.911355E-07   1.411792E-04  -1.221248E-04   -1.00013      -6.540232E-04  -9.202671E-05   3.081616E-02  -2.427777E-02
 civs   5   3.884337E-09  -2.175325E-05   2.040341E-05   5.490291E-04  -0.999792      -2.851072E-02  -2.441925E-02  -9.837957E-02
 civs   6  -5.534110E-03    1.72664       -1.17517       -10.2847        9.70793       -767.289        1402.89       -607.145    
 civs   7  -2.377336E-03    1.40726      -0.677582       -9.14714        11.2101       -1202.02       -810.558       -1142.63    
 civs   8  -1.180857E-03    1.26932      -0.342953       -11.0825        16.4453       -1751.69       -987.322        3850.51    

          calcsovref: sovref block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8
 ref    1  -0.947068      -2.243790E-04  -1.336307E-03  -2.320248E-03  -2.322630E-03  -4.256533E-02  -6.731114E-03   0.192571    
 ref    2   2.250138E-04  -0.946407       2.610204E-02  -2.199986E-02   1.449188E-03  -3.033638E-03  -7.939173E-03  -7.824527E-03
 ref    3   8.952033E-04  -2.481043E-02  -0.936073      -3.142263E-02   2.423363E-03   1.125223E-02   3.691348E-02  -4.241459E-02

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1    -0.94706797    -0.00022438    -0.00133631    -0.00232025    -0.00232263    -0.04256533    -0.00673111     0.19257097
 ref:   2     0.00022501    -0.94640733     0.02610204    -0.02199986     0.00144919    -0.00303364    -0.00793917    -0.00782453
 ref:   3     0.00089520    -0.02481043    -0.93607308    -0.03142263     0.00242336     0.01125223     0.03691348    -0.04241459

 trial vector basis is being transformed.  new dimension:   5

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 41  1    -94.3614190411  9.9476E-14  0.0000E+00  4.9784E-05  1.0000E-04
 mr-sdci # 41  2    -94.2184642534  9.0866E-08  1.0545E-07  3.9684E-04  1.0000E-04
 mr-sdci # 41  3    -94.2072629604  6.4807E-09  0.0000E+00  3.5711E-03  1.0000E-04
 mr-sdci # 41  4    -94.1743103986  6.2795E-06  0.0000E+00  3.2586E-02  1.0000E-04
 mr-sdci # 41  5    -94.1378527735  1.2553E-05  0.0000E+00  7.7937E-02  1.0000E-04
 mr-sdci # 41  6    -93.9091758056  2.7859E-01  0.0000E+00  6.0753E-01  1.0000E-04
 mr-sdci # 41  7    -92.5293380550  2.3459E-02  0.0000E+00  1.4236E+00  1.0000E-04
 mr-sdci # 41  8    -92.1669149319  6.6922E-01  0.0000E+00  1.3092E+00  1.0000E-04
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.006996
time for cinew                         0.097984
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration  42

 ========= Executing OUT-OF-CORE method ========


====================================================================================================
Diagonal     counts:  0x:        9055 2x:        2233 4x:         349
All internal counts: zz :        1113 yy:        8560 xx:       15526 ww:       21792
One-external counts: yz :        3142 yx:        7351 yw:        8696
Two-external counts: yy :        1938 ww:        4100 xx:        3262 xz:         403 wz:         604 wx:        4304
Three-ext.   counts: yx :        4024 yw:        4720

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================


 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6
   ht   1   -32.24770183
   ht   2     0.00000000   -32.10474704
   ht   3    -0.00000000     0.00000000   -32.09354574
   ht   4     0.00000000     0.00000000     0.00000000   -32.06059318
   ht   5    -0.00000000     0.00000000     0.00000000    -0.00000000   -32.02413556
   ht   6    -0.01683985    -0.00017652     0.00109047     0.00007751     0.00075746    -0.00001193

          calcsovref: tciref block   1

              civs   1       civs   2       civs   3       civs   4       civs   5       civs   6
 refs   1  -0.947068      -2.243790E-04  -1.336307E-03  -2.320248E-03  -2.322630E-03  -5.838749E-04
 refs   2   2.250138E-04  -0.946407       2.610204E-02  -2.199986E-02   1.449188E-03  -7.598884E-06
 refs   3   8.952033E-04  -2.481043E-02  -0.936073      -3.142263E-02   2.423363E-03   4.396791E-05

          calcsovref: scrb block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 civs   1   -1.00000      -2.234662E-04   6.197661E-05  -2.173519E-03   2.504144E-03   -1.62156    
 civs   2  -4.924885E-09   0.999998       1.769416E-06  -3.281228E-05   3.262350E-05  -1.720094E-02
 civs   3   3.497802E-08   1.342786E-05   0.999996       1.450965E-04  -1.649427E-04   0.105544    
 civs   4  -3.098811E-09  -8.734339E-06   3.342097E-06  -0.999989      -1.438587E-04   8.749855E-03
 civs   5   2.980998E-08   1.619865E-05  -4.763870E-06   2.289429E-04   -1.00011       7.201368E-02
 civs   6   1.035214E-03   0.427931      -0.118683        4.16222       -4.79536        3105.23    

          calcsovref: sovref block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.947068      -2.626344E-04  -1.325698E-03   1.947758E-03   2.751731E-03  -0.277662    
 ref    2  -2.250161E-04  -0.946408       2.610110E-02   2.200267E-02  -1.444364E-03  -5.015244E-03
 ref    3  -8.951907E-04  -2.480401E-02  -0.936075       3.146890E-02  -2.474126E-03   3.660777E-02

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 ref:   1     0.94706787    -0.00026263    -0.00132570     0.00194776     0.00275173    -0.27766244
 ref:   2    -0.00022502    -0.94640783     0.02610110     0.02200267    -0.00144436    -0.00501524
 ref:   3    -0.00089519    -0.02480401    -0.93607463     0.03146890    -0.00247413     0.03660777

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 42  1    -94.3614190411  2.6290E-13  0.0000E+00  4.9830E-05  1.0000E-04
 mr-sdci # 42  2    -94.2184642986  4.5124E-08  1.9676E-08  1.8392E-04  1.0000E-04
 mr-sdci # 42  3    -94.2072629639  3.4543E-09  0.0000E+00  3.5731E-03  1.0000E-04
 mr-sdci # 42  4    -94.1743145882  4.1896E-06  0.0000E+00  3.2382E-02  1.0000E-04
 mr-sdci # 42  5    -94.1378582491  5.4756E-06  0.0000E+00  7.7671E-02  1.0000E-04
 mr-sdci # 42  6    -91.8420837371 -2.0671E+00  0.0000E+00  1.5804E+00  1.0000E-04
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.004005
time for cinew                         0.063988
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration  43

 ========= Executing OUT-OF-CORE method ========


====================================================================================================
Diagonal     counts:  0x:        9055 2x:        2233 4x:         349
All internal counts: zz :        1113 yy:        8560 xx:       15526 ww:       21792
One-external counts: yz :        3142 yx:        7351 yw:        8696
Two-external counts: yy :        1938 ww:        4100 xx:        3262 xz:         403 wz:         604 wx:        4304
Three-ext.   counts: yx :        4024 yw:        4720

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================


 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7
   ht   1   -32.24770183
   ht   2     0.00000000   -32.10474704
   ht   3    -0.00000000     0.00000000   -32.09354574
   ht   4     0.00000000     0.00000000     0.00000000   -32.06059318
   ht   5    -0.00000000     0.00000000     0.00000000    -0.00000000   -32.02413556
   ht   6    -0.01683985    -0.00017652     0.00109047     0.00007751     0.00075746    -0.00001193
   ht   7     0.00454997    -0.00009142    -0.00023019     0.00008853     0.00054039     0.00000253    -0.00000120

          calcsovref: tciref block   1

              civs   1       civs   2       civs   3       civs   4       civs   5       civs   6       civs   7
 refs   1  -0.947068      -2.243790E-04  -1.336307E-03  -2.320248E-03  -2.322630E-03  -5.838749E-04   1.572610E-04
 refs   2   2.250138E-04  -0.946407       2.610204E-02  -2.199986E-02   1.449188E-03  -7.598884E-06  -3.080519E-06
 refs   3   8.952033E-04  -2.481043E-02  -0.936073      -3.142263E-02   2.423363E-03   4.396791E-05  -7.440382E-06

          calcsovref: scrb block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 civs   1  -1.000000      -1.569372E-04  -3.146766E-06  -1.598407E-03   3.312338E-03   0.127118        1.82000    
 civs   2  -3.307201E-08   0.999992      -9.422748E-06  -9.160983E-05  -4.843157E-05  -2.928084E-02   6.193906E-03
 civs   3   2.595222E-08   9.815783E-06  -0.999997       1.420062E-04  -1.734041E-04   5.995988E-03  -0.112381    
 civs   4   1.489776E-08  -9.415605E-06  -5.657939E-06  -0.999943       7.929191E-06   2.567479E-02   9.104654E-04
 civs   5   1.823349E-07   4.432609E-05   2.923027E-05   4.345504E-04  -0.999692       0.159291      -1.075169E-02
 civs   6   2.399719E-03   0.712294       0.372772        6.72375       -1.06205        1649.78       -2670.61    
 civs   7   7.604919E-03    1.52396        1.35735        13.5565        19.5452        7006.92        3015.02    

          calcsovref: sovref block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1   0.947068      -2.520718E-04   1.335036E-03   2.038810E-03   2.878939E-03   1.782884E-02   0.309952    
 ref    2  -2.250236E-04  -0.946409      -2.609990E-02   2.199643E-02  -1.459000E-03  -6.558748E-03   2.584376E-03
 ref    3  -8.951780E-04  -2.479917E-02   0.936077       3.148457E-02  -2.448500E-03   1.521013E-02  -3.323623E-02

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 ref:   1     0.94706793    -0.00025207     0.00133504     0.00203881     0.00287894     0.01782884     0.30995171
 ref:   2    -0.00022502    -0.94640910    -0.02609990     0.02199643    -0.00145900    -0.00655875     0.00258438
 ref:   3    -0.00089518    -0.02479917     0.93607709     0.03148457    -0.00244850     0.01521013    -0.03323623

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 43  1    -94.3614190411  9.0949E-13  0.0000E+00  4.9742E-05  1.0000E-04
 mr-sdci # 43  2    -94.2184643286  2.9985E-08  1.2798E-08  1.3393E-04  1.0000E-04
 mr-sdci # 43  3    -94.2072629873  2.3398E-08  0.0000E+00  3.5693E-03  1.0000E-04
 mr-sdci # 43  4    -94.1743168078  2.2196E-06  0.0000E+00  3.2312E-02  1.0000E-04
 mr-sdci # 43  5    -94.1378626002  4.3510E-06  0.0000E+00  7.7986E-02  1.0000E-04
 mr-sdci # 43  6    -93.5556230418  1.7135E+00  0.0000E+00  9.3874E-01  1.0000E-04
 mr-sdci # 43  7    -91.5248226542 -1.0045E+00  0.0000E+00  1.4587E+00  1.0000E-04
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.005997
time for cinew                         0.066994
time for eigenvalue solver             0.000999
time for vector access                 0.000000

          starting ci iteration  44

 ========= Executing OUT-OF-CORE method ========


====================================================================================================
Diagonal     counts:  0x:        9055 2x:        2233 4x:         349
All internal counts: zz :        1113 yy:        8560 xx:       15526 ww:       21792
One-external counts: yz :        3142 yx:        7351 yw:        8696
Two-external counts: yy :        1938 ww:        4100 xx:        3262 xz:         403 wz:         604 wx:        4304
Three-ext.   counts: yx :        4024 yw:        4720

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================


 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7         ht   8
   ht   1   -32.24770183
   ht   2     0.00000000   -32.10474704
   ht   3    -0.00000000     0.00000000   -32.09354574
   ht   4     0.00000000     0.00000000     0.00000000   -32.06059318
   ht   5    -0.00000000     0.00000000     0.00000000    -0.00000000   -32.02413556
   ht   6    -0.01683985    -0.00017652     0.00109047     0.00007751     0.00075746    -0.00001193
   ht   7     0.00454997    -0.00009142    -0.00023019     0.00008853     0.00054039     0.00000253    -0.00000120
   ht   8    -0.00141603     0.00018806    -0.00054959    -0.00006270    -0.00024478    -0.00000082     0.00000011    -0.00000056

          calcsovref: tciref block   1

              civs   1       civs   2       civs   3       civs   4       civs   5       civs   6       civs   7       civs   8
 refs   1  -0.947068      -2.243790E-04  -1.336307E-03  -2.320248E-03  -2.322630E-03  -5.838749E-04   1.572610E-04  -4.953981E-05
 refs   2   2.250138E-04  -0.946407       2.610204E-02  -2.199986E-02   1.449188E-03  -7.598884E-06  -3.080519E-06   7.657390E-06
 refs   3   8.952033E-04  -2.481043E-02  -0.936073      -3.142263E-02   2.423363E-03   4.396791E-05  -7.440382E-06  -2.262274E-05

          calcsovref: scrb block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 civs   1  -1.000000      -1.729822E-04  -2.461409E-05  -1.707911E-03   2.543490E-03   2.439224E-03  -9.689968E-02    1.85728    
 civs   2  -5.990916E-08   0.999995      -6.442142E-06  -7.015249E-05   1.098803E-04  -1.162369E-03   6.074422E-02  -3.061898E-03
 civs   3   1.462329E-07  -1.000035E-05   -1.00002       1.041188E-05  -1.189191E-03   6.924682E-02  -0.140534      -8.967862E-02
 civs   4   1.879133E-08  -1.101628E-05  -8.702870E-06  -0.999938       1.915663E-04  -1.347016E-02  -2.952563E-02   4.989826E-03
 civs   5   1.779323E-07   4.344777E-05   2.680471E-05   3.699312E-04  -0.999118      -0.101938      -0.156615       9.703904E-03
 civs   6   1.794000E-03   0.837790       0.548428        7.75614        7.50377       -1162.35       -1510.93       -2516.34    
 civs   7   4.161795E-03    2.19652        2.29461        18.9968        64.0374       -5545.24       -4121.17        3473.19    
 civs   8  -6.747124E-03    1.03401        1.41151        7.69714        58.6032       -4050.51        6932.99       -1211.44    

          calcsovref: sovref block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8
 ref    1   0.947068      -2.555758E-04   1.330320E-03   2.014283E-03   2.698972E-03   5.141144E-03  -1.698923E-02   0.316551    
 ref    2  -2.250316E-04  -0.946408      -2.609672E-02   2.200678E-02  -1.392124E-03  -2.044819E-03   1.650910E-02   2.482480E-05
 ref    3  -8.951383E-04  -2.480355E-02   0.936069       3.143758E-02  -2.786827E-03   1.717348E-02  -6.210782E-02  -2.352262E-02

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1     0.94706796    -0.00025558     0.00133032     0.00201428     0.00269897     0.00514114    -0.01698923     0.31655143
 ref:   2    -0.00022503    -0.94640792    -0.02609672     0.02200678    -0.00139212    -0.00204482     0.01650910     0.00002482
 ref:   3    -0.00089514    -0.02480355     0.93606918     0.03143758    -0.00278683     0.01717348    -0.06210782    -0.02352262

 trial vector basis is being transformed.  new dimension:   5

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 44  1    -94.3614190411  7.1054E-13  0.0000E+00  4.9709E-05  1.0000E-04
 mr-sdci # 44  2    -94.2184643418  1.3233E-08  4.9317E-09  1.0658E-04  1.0000E-04
 mr-sdci # 44  3    -94.2072630113  2.4035E-08  0.0000E+00  3.5661E-03  1.0000E-04
 mr-sdci # 44  4    -94.1743174660  6.5822E-07  0.0000E+00  3.2355E-02  1.0000E-04
 mr-sdci # 44  5    -94.1378969853  3.4385E-05  0.0000E+00  7.8523E-02  1.0000E-04
 mr-sdci # 44  6    -93.9018074751  3.4618E-01  0.0000E+00  5.3770E-01  1.0000E-04
 mr-sdci # 44  7    -92.5260426714  1.0012E+00  0.0000E+00  1.4062E+00  1.0000E-04
 mr-sdci # 44  8    -91.4980959691 -6.6882E-01  0.0000E+00  1.4423E+00  1.0000E-04
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.005005
time for cinew                         0.099983
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration  45

 ========= Executing OUT-OF-CORE method ========


====================================================================================================
Diagonal     counts:  0x:        9055 2x:        2233 4x:         349
All internal counts: zz :        1113 yy:        8560 xx:       15526 ww:       21792
One-external counts: yz :        3142 yx:        7351 yw:        8696
Two-external counts: yy :        1938 ww:        4100 xx:        3262 xz:         403 wz:         604 wx:        4304
Three-ext.   counts: yx :        4024 yw:        4720

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================


 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6
   ht   1   -32.24770183
   ht   2    -0.00000000   -32.10474713
   ht   3     0.00000000    -0.00000000   -32.09354580
   ht   4    -0.00000000     0.00000000     0.00000000   -32.06060025
   ht   5     0.00000000    -0.00000000    -0.00000000    -0.00000000   -32.02417977
   ht   6    -0.00015797    -0.00006362    -0.00038835     0.00003602    -0.00003360    -0.00000011

                v:   1         v:   2         v:   3         v:   4         v:   5         v:   6

   eig(s)   3.372183E-09    1.00000       1.000000       1.000000       1.000000       1.000000    
 
   x:   1   4.898764E-06  -0.123299       0.167200       0.829276       0.365136       0.368536    
   x:   2   1.981477E-06   0.112310      -0.165108      -0.403835       0.880240       0.149067    
   x:   3   1.210024E-05   6.402326E-02   7.892264E-02  -0.287839      -0.279578       0.910306    
   x:   4  -1.125239E-06   0.833429       0.543948       3.947546E-02   2.813787E-02  -8.465216E-02
   x:   5   1.041358E-06   0.522954      -0.801670       0.254577      -0.113567       7.834176E-02
   x:   6   -1.00000        0.00000        0.00000        0.00000        0.00000       1.329250E-05
 bummer (warning):overlap matrix: # small eigenvalues=                      1

          calcsovref: tciref block   1

              civs   1       civs   2       civs   3       civs   4       civs   5       civs   6
 refs   1   0.947068      -2.555758E-04   1.330320E-03   2.014283E-03   2.698972E-03   5.277118E-06
 refs   2  -2.250316E-04  -0.946408      -2.609672E-02   2.200678E-02  -1.392124E-03  -2.705259E-06
 refs   3  -8.951383E-04  -2.480355E-02   0.936069       3.143758E-02  -2.786827E-03   1.480771E-05

          calcsovref: scrb block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 civs   1   -1.00000      -5.441685E-06   7.785796E-06   6.684358E-05   2.946236E-04  -8.435979E-02
 civs   2   4.782791E-08   0.999998       3.848378E-06   2.855780E-05   1.228394E-04  -3.418619E-02
 civs   3   2.961298E-07  -1.274529E-05  -0.999981       1.679830E-04   7.337236E-04  -0.208462    
 civs   4  -3.532749E-08   2.722441E-06  -4.611650E-06   -1.00002       2.902982E-05   1.858505E-02
 civs   5  -1.976801E-09   2.300116E-06  -4.089629E-06  -7.972333E-05  -0.999931      -2.142472E-02
 civs   6  -2.456519E-02    1.11064       -1.58908       -13.6432       -60.1358        17220.3    

          calcsovref: sovref block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1  -0.947068      -2.548732E-04  -1.331328E-03  -2.023004E-03  -2.736098E-03   1.069032E-02
 ref    2   2.250443E-04  -0.946408       2.609678E-02  -2.200152E-02   1.419880E-03  -8.333391E-03
 ref    3   8.950494E-04  -2.479890E-02  -0.936075      -3.148338E-02   2.580581E-03   6.142589E-02

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 ref:   1    -0.94706797    -0.00025487    -0.00133133    -0.00202300    -0.00273610     0.01069032
 ref:   2     0.00022504    -0.94640844     0.02609678    -0.02200152     0.00141988    -0.00833339
 ref:   3     0.00089505    -0.02479890    -0.93607494    -0.03148338     0.00258058     0.06142589

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 45  1    -94.3614190411  2.9914E-12  0.0000E+00  4.9636E-05  1.0000E-04
 mr-sdci # 45  2    -94.2184643473  5.4774E-09  0.0000E+00  5.6384E-05  1.0000E-04
 mr-sdci # 45  3    -94.2072630225  1.1117E-08  6.2193E-06  3.5637E-03  1.0000E-04
 mr-sdci # 45  4    -94.1743182647  7.9866E-07  0.0000E+00  3.2330E-02  1.0000E-04
 mr-sdci # 45  5    -94.1379120647  1.5079E-05  0.0000E+00  7.7759E-02  1.0000E-04
 mr-sdci # 45  6    -92.9014068257 -1.0004E+00  0.0000E+00  1.2882E+00  1.0000E-04
 
 root number  3 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.002998
time for cinew                         0.126984
time for eigenvalue solver             0.000999
time for vector access                 0.000000

          starting ci iteration  46

 ========= Executing OUT-OF-CORE method ========


====================================================================================================
Diagonal     counts:  0x:        9055 2x:        2233 4x:         349
All internal counts: zz :        1113 yy:        8560 xx:       15526 ww:       21792
One-external counts: yz :        3142 yx:        7351 yw:        8696
Two-external counts: yy :        1938 ww:        4100 xx:        3262 xz:         403 wz:         604 wx:        4304
Three-ext.   counts: yx :        4024 yw:        4720

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================


 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7
   ht   1   -32.24770183
   ht   2    -0.00000000   -32.10474713
   ht   3     0.00000000    -0.00000000   -32.09354580
   ht   4    -0.00000000     0.00000000     0.00000000   -32.06060025
   ht   5     0.00000000    -0.00000000    -0.00000000    -0.00000000   -32.02417977
   ht   6    -0.00015797    -0.00006362    -0.00038835     0.00003602    -0.00003360    -0.00000011
   ht   7    -0.00104337    -0.01602477     0.00196438     0.00579142    -0.03898310     0.00000017    -0.00020494

                v:   1         v:   2         v:   3         v:   4         v:   5         v:   6         v:   7

   eig(s)   3.362060E-09   4.869720E-06   1.000000       1.000000       1.000000       1.000000        1.00000    
 
   x:   1  -4.945424E-06  -3.234810E-05   0.303906       0.844918       0.240235       0.368033       2.435909E-02
   x:   2  -2.701380E-06  -4.991387E-04   0.258419      -0.378465       0.797980       0.109724       0.375728    
   x:   3  -1.201222E-05   6.103200E-05  -8.190160E-02  -0.306506      -0.223784       0.920412      -4.591865E-02
   x:   4   1.385042E-06   1.801325E-04   0.913161      -0.203209      -0.318541      -7.062702E-02  -0.135596    
   x:   5  -2.795066E-06  -1.215923E-03   1.700239E-02   8.739373E-02  -0.392386      -1.912410E-02   0.915281    
   x:   6   0.999999      -1.442293E-03  -3.010566E-11   1.163747E-11   1.464973E-11   1.321708E-05   1.413901E-06
   x:   7   1.442290E-03   0.999998       1.777085E-14    0.00000      -3.827028E-14  -1.406720E-08   1.328471E-03
 bummer (warning):overlap matrix: # small eigenvalues=                      1

          calcsovref: tciref block   1

              civs   1       civs   2       civs   3       civs   4       civs   5       civs   6       civs   7
 refs   1   0.947068      -2.555758E-04   1.330320E-03   2.014283E-03   2.698972E-03   5.277118E-06   6.044811E-05
 refs   2  -2.250316E-04  -0.946408      -2.609672E-02   2.200678E-02  -1.392124E-03  -2.705259E-06  -6.447652E-04
 refs   3  -8.951383E-04  -2.480355E-02   0.936069       3.143758E-02  -2.786827E-03   1.480771E-05  -6.389277E-05

          calcsovref: scrb block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 civs   1   -1.00000      -5.429751E-06  -1.484885E-05   6.591018E-06   4.559808E-04   7.690481E-02  -3.968574E-02
 civs   2  -5.670003E-07   0.999998      -3.847528E-04  -1.011137E-03   2.978113E-03  -2.370703E-02  -0.229702    
 civs   3   3.276418E-07  -1.254565E-05  -0.999925       7.081030E-04  -1.928196E-04   0.206503      -3.448276E-02
 civs   4   7.921584E-08   2.800025E-06  -2.500321E-04  -0.999626      -3.582878E-03   4.072284E-03   8.916967E-02
 civs   5  -1.257931E-06   2.574445E-06  -4.570804E-04  -9.323719E-05  -0.992876      -0.120873      -0.552731    
 civs   6  -2.506068E-02    1.11091       -2.11644       -15.1174       -55.2165       -16444.0        5198.67    
 civs   7   1.233707E-03  -4.105626E-04   0.781109        2.08916       -5.74050        112.901        439.527    

          calcsovref: sovref block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1  -0.947068      -2.548842E-04  -1.309874E-03  -1.959830E-03  -2.894521E-03  -7.155707E-03   1.511811E-02
 ref    2   2.248355E-04  -0.946408       2.595613E-02  -2.236609E-02   2.340431E-03  -1.102102E-02  -7.642329E-02
 ref    3   8.950151E-04  -2.479868E-02  -0.936077      -3.109501E-02   1.948716E-03  -5.642716E-02   2.669606E-02

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 ref:   1    -0.94706794    -0.00025488    -0.00130987    -0.00195983    -0.00289452    -0.00715571     0.01511811
 ref:   2     0.00022484    -0.94640838     0.02595613    -0.02236609     0.00234043    -0.01102102    -0.07642329
 ref:   3     0.00089502    -0.02479868    -0.93607728    -0.03109501     0.00194872    -0.05642716     0.02669606

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 46  1    -94.3614190411  1.3273E-11  0.0000E+00  4.9810E-05  1.0000E-04
 mr-sdci # 46  2    -94.2184643473  1.3429E-12  0.0000E+00  5.6358E-05  1.0000E-04
 mr-sdci # 46  3    -94.2072678804  4.8579E-06  1.3009E-06  1.4929E-03  1.0000E-04
 mr-sdci # 46  4    -94.1743522456  3.3981E-05  0.0000E+00  3.2047E-02  1.0000E-04
 mr-sdci # 46  5    -94.1381648112  2.5275E-04  0.0000E+00  7.8672E-02  1.0000E-04
 mr-sdci # 46  6    -92.9260507955  2.4644E-02  0.0000E+00  1.2645E+00  1.0000E-04
 mr-sdci # 46  7    -92.5279307875  1.8881E-03  0.0000E+00  1.2856E+00  1.0000E-04
 
 root number  3 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.004997
time for cinew                         0.064987
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration  47

 ========= Executing OUT-OF-CORE method ========


====================================================================================================
Diagonal     counts:  0x:        9055 2x:        2233 4x:         349
All internal counts: zz :        1113 yy:        8560 xx:       15526 ww:       21792
One-external counts: yz :        3142 yx:        7351 yw:        8696
Two-external counts: yy :        1938 ww:        4100 xx:        3262 xz:         403 wz:         604 wx:        4304
Three-ext.   counts: yx :        4024 yw:        4720

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================


 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7         ht   8
   ht   1   -32.24770183
   ht   2    -0.00000000   -32.10474713
   ht   3     0.00000000    -0.00000000   -32.09354580
   ht   4    -0.00000000     0.00000000     0.00000000   -32.06060025
   ht   5     0.00000000    -0.00000000    -0.00000000    -0.00000000   -32.02417977
   ht   6    -0.00015797    -0.00006362    -0.00038835     0.00003602    -0.00003360    -0.00000011
   ht   7    -0.00104337    -0.01602477     0.00196438     0.00579142    -0.03898310     0.00000017    -0.00020494
   ht   8    -0.00054225     0.02504850     0.00268176    -0.00027039     0.01340312    -0.00000004     0.00004665    -0.00005786

                v:   1         v:   2         v:   3         v:   4         v:   5         v:   6         v:   7         v:   8

   eig(s)   3.350842E-09   9.643056E-07   4.965586E-06   1.000000       1.000000       1.000000       1.000000        1.00000    
 
   x:   1   4.874738E-06  -2.163618E-05  -2.935498E-05  -0.910527       0.101931      -0.395860       6.069583E-02   1.269581E-02
   x:   2   5.128927E-06   6.935355E-04  -6.138896E-04  -0.111876      -8.917907E-02   0.126839      -0.815540       0.546202    
   x:   3   1.232085E-05   9.195900E-05   4.736357E-05   0.383123       4.828962E-02  -0.899993      -0.202128      -6.445592E-03
   x:   4  -1.336989E-06   1.953870E-05   1.792693E-04   4.821528E-02   0.973014       0.110799      -0.166145      -0.105061    
   x:   5   3.692912E-06   2.249835E-04  -1.266007E-03   9.652199E-02   0.180467      -7.030075E-02   0.512592       0.830916    
   x:   6  -0.999994       3.185419E-03  -1.528331E-03   2.702105E-12   3.757259E-12  -1.277591E-05  -3.043687E-06   2.049974E-06
   x:   7  -1.016837E-03   0.154792       0.987946      -7.360521E-16  -8.002875E-14  -2.408409E-08   2.604277E-04   1.302695E-03
   x:   8   3.383591E-03   0.987942      -0.154786       1.584630E-14  -6.895911E-14  -7.439060E-08   4.384125E-04  -7.738158E-04
 bummer (warning):overlap matrix: # small eigenvalues=                      1

          calcsovref: tciref block   1

              civs   1       civs   2       civs   3       civs   4       civs   5       civs   6       civs   7       civs   8
 refs   1   0.947068      -2.555758E-04   1.330320E-03   2.014283E-03   2.698972E-03   5.277118E-06   6.044811E-05   3.071898E-05
 refs   2  -2.250316E-04  -0.946408      -2.609672E-02   2.200678E-02  -1.392124E-03  -2.705259E-06  -6.447652E-04   9.281798E-04
 refs   3  -8.951383E-04  -2.480355E-02   0.936069       3.143758E-02  -2.786827E-03   1.480771E-05  -6.389277E-05  -7.175875E-05

          calcsovref: scrb block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 civs   1    1.00000       5.441390E-06   4.099588E-05  -7.517857E-05   4.611685E-04  -3.680534E-02  -7.530790E-02   2.691957E-02
 civs   2   1.748816E-07  -0.999998      -3.453104E-04   1.253693E-03   2.835727E-03   0.517528       6.886854E-03   0.560937    
 civs   3  -3.780581E-07   1.255577E-05   0.999822       1.159333E-03  -2.167921E-04   6.590531E-02  -0.208728       8.348067E-02
 civs   4  -8.364019E-08  -2.771549E-06   3.464325E-04  -0.999518      -3.612332E-03   5.856886E-02  -6.882944E-03  -6.981122E-02
 civs   5   1.158882E-06  -2.654107E-06   2.787087E-04   4.641718E-04  -0.992909      -3.504550E-02   0.126624       0.613354    
 civs   6   2.504278E-02   -1.11096        1.98056       -14.6114       -55.2541        2376.73        16327.0       -5119.21    
 civs   7  -1.350207E-03   3.099667E-04   -1.00717        2.79598       -5.78536        322.608       -128.117       -318.094    
 civs   8  -5.771554E-04  -4.821678E-04   -1.08077        3.35638      -0.211289        875.755       -31.5908        502.435    

          calcsovref: sovref block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8
 ref    1   0.947068       2.548739E-04   1.286817E-03  -1.887027E-03  -2.899153E-03   2.406729E-02   6.171312E-03  -3.831401E-03
 ref    2  -2.249236E-04   0.946408      -2.611715E-02  -2.186147E-02   2.308115E-03   0.108256       7.733097E-03   0.149846    
 ref    3  -8.950038E-04   2.479874E-02   0.936092      -3.100540E-02   1.946440E-03   2.565848E-03   5.616173E-02  -3.123160E-02

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1     0.94706793     0.00025487     0.00128682    -0.00188703    -0.00289915     0.02406729     0.00617131    -0.00383140
 ref:   2    -0.00022492     0.94640830    -0.02611715    -0.02186147     0.00230811     0.10825572     0.00773310     0.14984564
 ref:   3    -0.00089500     0.02479874     0.93609205    -0.03100540     0.00194644     0.00256585     0.05616173    -0.03123160

 trial vector basis is being transformed.  new dimension:   5

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 47  1    -94.3614190411  4.7606E-13  0.0000E+00  4.9662E-05  1.0000E-04
 mr-sdci # 47  2    -94.2184643473  3.1974E-13  0.0000E+00  5.6342E-05  1.0000E-04
 mr-sdci # 47  3    -94.2072692863  1.4060E-06  2.9767E-07  8.0248E-04  1.0000E-04
 mr-sdci # 47  4    -94.1743654336  1.3188E-05  0.0000E+00  3.0758E-02  1.0000E-04
 mr-sdci # 47  5    -94.1381648619  5.0631E-08  0.0000E+00  7.8646E-02  1.0000E-04
 mr-sdci # 47  6    -93.1030421730  1.7699E-01  0.0000E+00  9.6100E-01  1.0000E-04
 mr-sdci # 47  7    -92.9257951167  3.9786E-01  0.0000E+00  1.2625E+00  1.0000E-04
 mr-sdci # 47  8    -92.3389887513  8.4089E-01  0.0000E+00  1.3428E+00  1.0000E-04
 
 root number  3 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.006004
time for cinew                         0.099983
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration  48

 ========= Executing OUT-OF-CORE method ========


====================================================================================================
Diagonal     counts:  0x:        9055 2x:        2233 4x:         349
All internal counts: zz :        1113 yy:        8560 xx:       15526 ww:       21792
One-external counts: yz :        3142 yx:        7351 yw:        8696
Two-external counts: yy :        1938 ww:        4100 xx:        3262 xz:         403 wz:         604 wx:        4304
Three-ext.   counts: yx :        4024 yw:        4720

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================


 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6
   ht   1   -32.24770183
   ht   2     0.00000000   -32.10474713
   ht   3     0.00000000    -0.00000000   -32.09355207
   ht   4     0.00000000    -0.00000000    -0.00000000   -32.06064822
   ht   5     0.00000000    -0.00000000    -0.00000000    -0.00000000   -32.02444765
   ht   6     0.00186236    -0.01390800     0.00147693    -0.00117314     0.00285600    -0.00001213

          calcsovref: tciref block   1

              civs   1       civs   2       civs   3       civs   4       civs   5       civs   6
 refs   1   0.947068       2.548739E-04   1.286817E-03  -1.887027E-03  -2.899153E-03  -6.857880E-05
 refs   2  -2.249236E-04   0.946408      -2.611715E-02  -2.186147E-02   2.308115E-03   4.946513E-04
 refs   3  -8.950038E-04   2.479874E-02   0.936092      -3.100540E-02   1.946440E-03  -3.300131E-05

          calcsovref: scrb block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 civs   1  -0.999999       8.111922E-08  -4.531688E-05   7.631258E-04   5.174003E-04   0.133411    
 civs   2  -6.899875E-06   0.999999       3.397202E-04  -5.721072E-03  -3.879225E-03   -1.00069    
 civs   3   7.638839E-07   1.019219E-07   -1.00004       4.881461E-04   3.735775E-04   0.105986    
 civs   4  -1.597060E-07   1.060872E-07  -8.915369E-05   0.999502      -1.549264E-03  -8.988094E-02
 civs   5   1.183550E-06   6.738614E-08  -3.243061E-05  -2.443544E-05  -0.999194       0.209652    
 civs   6   1.592757E-02   1.403515E-03  -0.784112        13.2060        8.95453        2309.96    

          calcsovref: sovref block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1  -0.947068       2.548541E-04  -1.275659E-03  -2.069763E-03   2.775154E-03  -3.262205E-02
 ref    2   2.262581E-04   0.946408       2.605363E-02  -2.074568E-02  -1.524219E-03   0.195213    
 ref    3   8.950286E-04   2.479877E-02  -0.936089      -3.111142E-02  -1.939307E-03   1.240580E-03

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 ref:   1    -0.94706815     0.00025485    -0.00127566    -0.00206976     0.00277515    -0.03262205
 ref:   2     0.00022626     0.94640842     0.02605363    -0.02074568    -0.00152422     0.19521305
 ref:   3     0.00089503     0.02479877    -0.93608873    -0.03111142    -0.00193931     0.00124058

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 48  1    -94.3614190412  1.0388E-10  0.0000E+00  4.6089E-05  1.0000E-04
 mr-sdci # 48  2    -94.2184643473  7.8160E-13  0.0000E+00  5.6341E-05  1.0000E-04
 mr-sdci # 48  3    -94.2072695198  2.3341E-07  8.1684E-08  3.7454E-04  1.0000E-04
 mr-sdci # 48  4    -94.1744306740  6.5240E-05  0.0000E+00  2.7897E-02  1.0000E-04
 mr-sdci # 48  5    -94.1381943908  2.9529E-05  0.0000E+00  7.7988E-02  1.0000E-04
 mr-sdci # 48  6    -92.1766031891 -9.2644E-01  0.0000E+00  1.5863E+00  1.0000E-04
 
 root number  3 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.005005
time for cinew                         0.062988
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration  49

 ========= Executing OUT-OF-CORE method ========


====================================================================================================
Diagonal     counts:  0x:        9055 2x:        2233 4x:         349
All internal counts: zz :        1113 yy:        8560 xx:       15526 ww:       21792
One-external counts: yz :        3142 yx:        7351 yw:        8696
Two-external counts: yy :        1938 ww:        4100 xx:        3262 xz:         403 wz:         604 wx:        4304
Three-ext.   counts: yx :        4024 yw:        4720

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================


 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7
   ht   1   -32.24770183
   ht   2     0.00000000   -32.10474713
   ht   3     0.00000000    -0.00000000   -32.09355207
   ht   4     0.00000000    -0.00000000    -0.00000000   -32.06064822
   ht   5     0.00000000    -0.00000000    -0.00000000    -0.00000000   -32.02444765
   ht   6     0.00186236    -0.01390800     0.00147693    -0.00117314     0.00285600    -0.00001213
   ht   7     0.00133090    -0.00664441     0.00028303    -0.00009135     0.00087536    -0.00000289    -0.00000342

          calcsovref: tciref block   1

              civs   1       civs   2       civs   3       civs   4       civs   5       civs   6       civs   7
 refs   1   0.947068       2.548739E-04   1.286817E-03  -1.887027E-03  -2.899153E-03  -6.857880E-05  -4.180178E-05
 refs   2  -2.249236E-04   0.946408      -2.611715E-02  -2.186147E-02   2.308115E-03   4.946513E-04   2.467603E-04
 refs   3  -8.950038E-04   2.479874E-02   0.936092      -3.100540E-02   1.946440E-03  -3.300131E-05  -4.989340E-06

          calcsovref: scrb block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 civs   1   0.999999      -1.609686E-07  -1.681466E-05  -8.195787E-04   1.529696E-03   4.408603E-02  -0.207514    
 civs   2   8.206219E-06  -0.999999       2.409713E-04   5.914327E-03  -7.286296E-03  -1.637836E-02    1.30167    
 civs   3  -7.180129E-07  -6.749748E-08   -1.00004      -4.848087E-04   2.485132E-04  -3.808112E-02  -0.105398    
 civs   4   5.215265E-08  -1.416773E-07  -7.498286E-05  -0.999519      -1.387996E-03   4.635560E-02   7.777764E-02
 civs   5  -1.072205E-06   3.230019E-08  -7.451562E-05   1.904898E-04  -0.999072      -3.474952E-02  -0.235053    
 civs   6  -1.248570E-02  -5.494583E-04   -1.09167       -12.5807       -2.66443       -1425.49       -1818.08    
 civs   7  -1.351663E-02  -3.129313E-03    1.12081       -2.24253        40.7825        3062.95       -2483.84    

          calcsovref: sovref block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1   0.947068      -2.548576E-04  -1.274361E-03   2.066765E-03   2.824211E-03   1.143407E-02   3.271177E-02
 ref    2  -2.266531E-04  -0.946409       2.608431E-02   2.068520E-02  -4.327290E-04   3.508505E-02  -0.279769    
 ref    3  -8.949955E-04  -2.479875E-02  -0.936091       3.111081E-02  -1.966575E-03  -5.837037E-03   3.326459E-03

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 ref:   1     0.94706814    -0.00025486    -0.00127436     0.00206676     0.00282421     0.01143407     0.03271177
 ref:   2    -0.00022665    -0.94640851     0.02608431     0.02068520    -0.00043273     0.03508505    -0.27976886
 ref:   3    -0.00089500    -0.02479875    -0.93609108     0.03111081    -0.00196658    -0.00583704     0.00332646

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 49  1    -94.3614190412  1.5504E-11  0.0000E+00  4.5212E-05  1.0000E-04
 mr-sdci # 49  2    -94.2184643473  6.9633E-13  0.0000E+00  5.6383E-05  1.0000E-04
 mr-sdci # 49  3    -94.2072696113  9.1552E-08  2.7963E-08  2.4466E-04  1.0000E-04
 mr-sdci # 49  4    -94.1744310268  3.5285E-07  0.0000E+00  2.7805E-02  1.0000E-04
 mr-sdci # 49  5    -94.1383066299  1.1224E-04  0.0000E+00  7.6058E-02  1.0000E-04
 mr-sdci # 49  6    -93.3926890447  1.2161E+00  0.0000E+00  9.9843E-01  1.0000E-04
 mr-sdci # 49  7    -91.3769812722 -1.5488E+00  0.0000E+00  1.4535E+00  1.0000E-04
 
 root number  3 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.004005
time for cinew                         0.064987
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration  50

 ========= Executing OUT-OF-CORE method ========


====================================================================================================
Diagonal     counts:  0x:        9055 2x:        2233 4x:         349
All internal counts: zz :        1113 yy:        8560 xx:       15526 ww:       21792
One-external counts: yz :        3142 yx:        7351 yw:        8696
Two-external counts: yy :        1938 ww:        4100 xx:        3262 xz:         403 wz:         604 wx:        4304
Three-ext.   counts: yx :        4024 yw:        4720

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================


 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7         ht   8
   ht   1   -32.24770183
   ht   2     0.00000000   -32.10474713
   ht   3     0.00000000    -0.00000000   -32.09355207
   ht   4     0.00000000    -0.00000000    -0.00000000   -32.06064822
   ht   5     0.00000000    -0.00000000    -0.00000000    -0.00000000   -32.02444765
   ht   6     0.00186236    -0.01390800     0.00147693    -0.00117314     0.00285600    -0.00001213
   ht   7     0.00133090    -0.00664441     0.00028303    -0.00009135     0.00087536    -0.00000289    -0.00000342
   ht   8    -0.00318738     0.00101072     0.00027686    -0.00045328     0.00047311     0.00000099     0.00000019    -0.00000096

          calcsovref: tciref block   1

              civs   1       civs   2       civs   3       civs   4       civs   5       civs   6       civs   7       civs   8
 refs   1   0.947068       2.548739E-04   1.286817E-03  -1.887027E-03  -2.899153E-03  -6.857880E-05  -4.180178E-05   1.116886E-04
 refs   2  -2.249236E-04   0.946408      -2.611715E-02  -2.186147E-02   2.308115E-03   4.946513E-04   2.467603E-04  -4.340107E-05
 refs   3  -8.950038E-04   2.479874E-02   0.936092      -3.100540E-02   1.946440E-03  -3.300131E-05  -4.989340E-06  -8.582012E-06

          calcsovref: scrb block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 civs   1   0.999997       1.980126E-06   1.089129E-04  -1.572497E-03   1.425651E-03   0.244309      -0.580592       0.394788    
 civs   2   8.246124E-06  -0.999999      -2.422681E-04   5.923512E-03  -7.285283E-03   9.658919E-03  -0.295683       -1.27019    
 civs   3  -5.607247E-07  -1.849843E-07    1.00003      -4.173862E-04   2.570783E-04  -1.903177E-03   0.119186       7.445938E-02
 civs   4  -2.138501E-07   2.476598E-07   9.340388E-05  -0.999611      -1.393499E-03   7.476562E-03  -0.146942      -3.718894E-02
 civs   5  -7.273828E-07  -3.613315E-07   5.768272E-05   3.238962E-04  -0.999053      -3.257225E-02   0.208561       0.182490    
 civs   6  -1.380011E-02   1.623484E-03    1.18765       -13.4280       -2.79203        1048.33        1547.84        1458.96    
 civs   7  -7.939839E-03  -1.112990E-02   -1.47057       0.760178        41.2213       -2694.99       -833.042        2782.31    
 civs   8   1.983938E-02  -2.373206E-02   -1.02171        8.37592        1.16129       -2984.47        6430.58       -1979.91    

          calcsovref: sovref block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8
 ref    1   0.947069      -2.552948E-04   1.255513E-03   2.221660E-03   2.845750E-03  -6.111288E-02   9.678757E-02  -6.428875E-02
 ref    2  -2.267476E-04  -0.946408      -2.608029E-02   2.065291E-02  -4.370465E-04  -8.030454E-03   1.859254E-03   0.291252    
 ref    3  -8.949915E-04  -2.479870E-02   0.936091       3.111903E-02  -1.966176E-03   2.406782E-03   7.606210E-03  -5.681378E-03

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1     0.94706850    -0.00025529     0.00125551     0.00222166     0.00284575    -0.06111288     0.09678757    -0.06428875
 ref:   2    -0.00022675    -0.94640841    -0.02608029     0.02065291    -0.00043705    -0.00803045     0.00185925     0.29125173
 ref:   3    -0.00089499    -0.02479870     0.93609062     0.03111903    -0.00196618     0.00240678     0.00760621    -0.00568138

 trial vector basis is being transformed.  new dimension:   5

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 50  1    -94.3614190412  1.2669E-11  0.0000E+00  4.4967E-05  1.0000E-04
 mr-sdci # 50  2    -94.2184643473  1.5646E-11  0.0000E+00  5.6318E-05  1.0000E-04
 mr-sdci # 50  3    -94.2072696399  2.8570E-08  9.2609E-09  1.2142E-04  1.0000E-04
 mr-sdci # 50  4    -94.1744328675  1.8406E-06  0.0000E+00  2.7753E-02  1.0000E-04
 mr-sdci # 50  5    -94.1383066635  3.3628E-08  0.0000E+00  7.6058E-02  1.0000E-04
 mr-sdci # 50  6    -93.7228705139  3.3018E-01  0.0000E+00  4.5206E-01  1.0000E-04
 mr-sdci # 50  7    -91.7427859985  3.6580E-01  0.0000E+00  1.3539E+00  1.0000E-04
 mr-sdci # 50  8    -91.3434313276 -9.9556E-01  0.0000E+00  1.4600E+00  1.0000E-04
 
 root number  3 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.005997
time for cinew                         0.101982
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration  51

 ========= Executing OUT-OF-CORE method ========


====================================================================================================
Diagonal     counts:  0x:        9055 2x:        2233 4x:         349
All internal counts: zz :        1113 yy:        8560 xx:       15526 ww:       21792
One-external counts: yz :        3142 yx:        7351 yw:        8696
Two-external counts: yy :        1938 ww:        4100 xx:        3262 xz:         403 wz:         604 wx:        4304
Three-ext.   counts: yx :        4024 yw:        4720

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================


 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6
   ht   1   -32.24770183
   ht   2     0.00000000   -32.10474713
   ht   3    -0.00000000     0.00000000   -32.09355242
   ht   4    -0.00000000     0.00000000    -0.00000000   -32.06071565
   ht   5     0.00000000    -0.00000000     0.00000000     0.00000000   -32.02458945
   ht   6    -0.00528646     0.00003619     0.00005211     0.00006709    -0.00050824    -0.00000115

                v:   1         v:   2         v:   3         v:   4         v:   5         v:   6

   eig(s)   9.224126E-09    1.00000       1.000000       1.000000       1.000000       1.000000    
 
   x:   1   1.639329E-04   5.953391E-03  -5.621445E-02  -3.787033E-03  -7.975313E-02   0.995203    
   x:   2  -1.127151E-06   6.539626E-03  -2.175945E-02  -0.999474      -2.210206E-02  -6.842701E-03
   x:   3  -1.623984E-06  -0.764219      -0.598022       2.750232E-03   0.241317      -9.858878E-03
   x:   4  -2.081511E-06  -0.606581       0.539932      -2.735881E-03  -0.583408      -1.263643E-02
   x:   5   1.585686E-05  -0.218976       0.589245      -3.197145E-02   0.771072       9.626376E-02
   x:   6   -1.00000        0.00000        0.00000        0.00000        0.00000       1.647231E-04
 bummer (warning):overlap matrix: # small eigenvalues=                      1

          calcsovref: tciref block   1

              civs   1       civs   2       civs   3       civs   4       civs   5       civs   6
 refs   1   0.947069      -2.552948E-04   1.255513E-03   2.221660E-03   2.845750E-03   1.824507E-04
 refs   2  -2.267476E-04  -0.946408      -2.608029E-02   2.065291E-02  -4.370465E-04   2.418424E-06
 refs   3  -8.949915E-04  -2.479870E-02   0.936091       3.111903E-02  -1.966176E-03  -2.011108E-06

          calcsovref: scrb block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 civs   1    1.00000      -3.584338E-06   6.752882E-05  -2.646664E-03   3.199394E-03   -1.70687    
 civs   2  -2.601995E-08   -1.00000      -4.824898E-07   1.837905E-05  -2.211860E-05   1.173804E-02
 civs   3  -3.899550E-08   5.359752E-08   -1.00000       2.166519E-05  -2.907375E-05   1.686945E-02
 civs   4  -3.853217E-09  -1.323120E-07   3.634902E-06  -0.999965      -2.341275E-04   2.322304E-02
 civs   5   3.229069E-07  -2.303609E-07   3.984148E-06  -6.541372E-05  -0.999689      -0.166977    
 civs   6  -2.315565E-02   2.186518E-02  -0.411939        16.1451       -19.5168        10412.0    

          calcsovref: sovref block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.947068       2.558886E-04  -1.266699E-03  -1.782636E-03  -3.376220E-03   0.282752    
 ref    2  -2.267790E-04   0.946408       2.607982E-02  -2.063048E-02   4.058412E-04   1.457143E-02
 ref    3  -8.949849E-04   2.479871E-02  -0.936090      -3.112810E-02   1.967998E-03  -2.860880E-03

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 ref:   1     0.94706787     0.00025589    -0.00126670    -0.00178264    -0.00337622     0.28275207
 ref:   2    -0.00022678     0.94640844     0.02607982    -0.02063048     0.00040584     0.01457143
 ref:   3    -0.00089498     0.02479871    -0.93609036    -0.03112810     0.00196800    -0.00286088

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 51  1    -94.3614190412  1.2882E-11  0.0000E+00  4.5000E-05  1.0000E-04
 mr-sdci # 51  2    -94.2184643473  1.0793E-11  0.0000E+00  5.5897E-05  1.0000E-04
 mr-sdci # 51  3    -94.2072696437  3.8149E-09  1.7757E-09  6.0108E-05  1.0000E-04
 mr-sdci # 51  4    -94.1744386490  5.7815E-06  0.0000E+00  2.7510E-02  1.0000E-04
 mr-sdci # 51  5    -94.1383149883  8.3248E-06  0.0000E+00  7.5950E-02  1.0000E-04
 mr-sdci # 51  6    -91.7693170107 -1.9536E+00  0.0000E+00  1.5526E+00  1.0000E-04
 
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.002998
time for cinew                         0.066994
time for eigenvalue solver             0.000000
time for vector access                 0.000000

 mr-sdci  convergence criteria satisfied after 51 iterations.

 final mr-sdci  convergence information:

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 51  1    -94.3614190412  1.2882E-11  0.0000E+00  4.5000E-05  1.0000E-04
 mr-sdci # 51  2    -94.2184643473  1.0793E-11  0.0000E+00  5.5897E-05  1.0000E-04
 mr-sdci # 51  3    -94.2072696437  3.8149E-09  1.7757E-09  6.0108E-05  1.0000E-04
 mr-sdci # 51  4    -94.1744386490  5.7815E-06  0.0000E+00  2.7510E-02  1.0000E-04
 mr-sdci # 51  5    -94.1383149883  8.3248E-06  0.0000E+00  7.5950E-02  1.0000E-04
 mr-sdci # 51  6    -91.7693170107 -1.9536E+00  0.0000E+00  1.5526E+00  1.0000E-04

####################CIUDGINFO####################

   ci vector at position   1 energy=  -94.361419041214
   ci vector at position   2 energy=  -94.218464347297
   ci vector at position   3 energy=  -94.207269643688

################END OF CIUDGINFO################

 
diagon:itrnv=   0
    3 of the   7 expansion vectors are transformed.
    3 of the   6 matrix-vector products are transformed.

    3 expansion eigenvectors written to unit nvfile (= 11)
    3 matrix-vector products written to unit nhvfil (= 10)
maximum overlap with reference    1(overlap= 0.94707)

 information on vector: 1 from unit 11 written to unit 48 filename civout              
maximum overlap with reference    2(overlap= 0.94641)

 information on vector: 2 from unit 11 written to unit 48 filename civout              
maximum overlap with reference    3(overlap= 0.93609)

 information on vector: 3 from unit 11 written to unit 48 filename civout              


 --- list of ci coefficients ( ctol =   1.00E-02 )  total energy( 1) =       -94.3614190412

                                                       internal orbitals

                                          level       1    2    3    4    5    6    7

                                          orbital     3    4    5    6    7    8    9

                                         symmetry   a    a    a    a    a    a    a  

 path  s ms    csf#    c(i)    ext. orb.(sym)
 z*  1  1       1 -0.031739                        +-   +-   +-   +-   +-   +-      
 z*  1  1       2  0.946293                        +-   +-   +-   +-   +-   +     - 
 z*  1  1       3  0.018921                        +-   +-   +-   +-   +-        +- 
 y   1  1      37  0.012234              9( a  )   +-   +-   +-   +-   +-    -      
 y   1  1      39 -0.015089             11( a  )   +-   +-   +-   +-   +-    -      
 y   1  1      47  0.013179             19( a  )   +-   +-   +-   +-   +-    -      
 y   1  1      74 -0.021597             46( a  )   +-   +-   +-   +-   +-    -      
 y   1  1     107  0.033577              6( a  )   +-   +-   +-   +-   +-         - 
 y   1  1     108  0.021276              7( a  )   +-   +-   +-   +-   +-         - 
 y   1  1     123  0.010239             22( a  )   +-   +-   +-   +-   +-         - 
 y   1  1     124 -0.010287             23( a  )   +-   +-   +-   +-   +-         - 
 y   1  1     126 -0.021718             25( a  )   +-   +-   +-   +-   +-         - 
 y   1  1     128 -0.011367             27( a  )   +-   +-   +-   +-   +-         - 
 y   1  1     129 -0.014053             28( a  )   +-   +-   +-   +-   +-         - 
 y   1  1     137 -0.015323             36( a  )   +-   +-   +-   +-   +-         - 
 y   1  1     139 -0.017760             38( a  )   +-   +-   +-   +-   +-         - 
 y   1  1     142 -0.028031             41( a  )   +-   +-   +-   +-   +-         - 
 y   1  1     250  0.012099              3( a  )   +-   +-   +-   +-    -   +     - 
 y   1  1     251  0.020327              4( a  )   +-   +-   +-   +-    -   +     - 
 y   1  1     252 -0.020645              5( a  )   +-   +-   +-   +-    -   +     - 
 y   1  1     254  0.011603              7( a  )   +-   +-   +-   +-    -   +     - 
 y   1  1     258  0.015272             11( a  )   +-   +-   +-   +-    -   +     - 
 y   1  1     259 -0.018663             12( a  )   +-   +-   +-   +-    -   +     - 
 y   1  1     261  0.014263             14( a  )   +-   +-   +-   +-    -   +     - 
 y   1  1     262 -0.042082             15( a  )   +-   +-   +-   +-    -   +     - 
 y   1  1     279 -0.019852             32( a  )   +-   +-   +-   +-    -   +     - 
 y   1  1     292 -0.013596             45( a  )   +-   +-   +-   +-    -   +     - 
 y   1  1     408 -0.018597             15( a  )   +-   +-   +-   +-   +     -    - 
 y   1  1     687 -0.012278              2( a  )   +-   +-   +-    -   +-   +     - 
 y   1  1     693 -0.010360              8( a  )   +-   +-   +-    -   +-   +     - 
 y   1  1     697  0.013849             12( a  )   +-   +-   +-    -   +-   +     - 
 y   1  1     698  0.018609             13( a  )   +-   +-   +-    -   +-   +     - 
 y   1  1     699 -0.012882             14( a  )   +-   +-   +-    -   +-   +     - 
 y   1  1     729  0.010965             44( a  )   +-   +-   +-    -   +-   +     - 
 y   1  1    1058 -0.010547              8( a  )   +-   +-   +-   +    +-    -    - 
 y   1  1    1068  0.015017             18( a  )   +-   +-   +-   +    +-    -    - 
 y   1  1    1086 -0.029177             36( a  )   +-   +-   +-   +    +-    -    - 
 y   1  1    1111  0.016318             61( a  )   +-   +-   +-   +    +-    -    - 
 y   1  1    1115  0.010513             65( a  )   +-   +-   +-   +    +-    -    - 
 y   1  1    1235  0.013764             39( a  )   +-   +-   +-   +     -    -   +- 
 y   1  1    1236  0.010832             40( a  )   +-   +-   +-   +     -    -   +- 
 y   1  1    1263  0.016608             67( a  )   +-   +-   +-   +     -    -   +- 
 y   1  1    1361 -0.011121             19( a  )   +-   +-   +-        +-    -   +- 
 y   1  1    1388  0.012649             46( a  )   +-   +-   +-        +-    -   +- 
 y   1  1    1410 -0.010745             68( a  )   +-   +-   +-        +-    -   +- 
 y   1  1    1603  0.011873             42( a  )   +-   +-    -   +-   +-   +     - 
 y   1  1    3024 -0.011099              3( a  )   +-    -   +-   +-   +-   +     - 
 y   1  1    3032 -0.010194             11( a  )   +-    -   +-   +-   +-   +     - 
 y   1  1    3033  0.015749             12( a  )   +-    -   +-   +-   +-   +     - 
 y   1  1    3497  0.010046             38( a  )   +-    -   +-   +    +-    -   +- 
 y   1  1    3527 -0.010717             68( a  )   +-    -   +-   +    +-    -   +- 
 y   1  1    6314 -0.018487              8( a  )   +    +-   +-   +-   +-    -    - 
 y   1  1    6320  0.012130             14( a  )   +    +-   +-   +-   +-    -    - 
 y   1  1    6325  0.010303             19( a  )   +    +-   +-   +-   +-    -    - 
 y   1  1    6328 -0.015184             22( a  )   +    +-   +-   +-   +-    -    - 
 y   1  1    6357 -0.028596             51( a  )   +    +-   +-   +-   +-    -    - 
 y   1  1    6778 -0.013202             34( a  )   +    +-    -   +-   +-   +-    - 
 y   1  1    6792  0.016305             48( a  )   +    +-    -   +-   +-   +-    - 
 x   1  1    7788  0.010242   15( a  )  19( a  )   +-   +-   +-   +-    -    -      
 x   1  1    8625 -0.014744   15( a  )  46( a  )   +-   +-   +-   +-    -    -      
 x   1  1   23523 -0.010554   15( a  )  17( a  )   +-   +-   +-    -    -   +     - 
 x   1  1  100634  0.010125   44( a  )  46( a  )   +-    -   +-   +-   +-    -      
 x   1  1  107619 -0.012106   15( a  )  17( a  )   +-    -   +-   +-    -   +     - 
 x   1  1  182245 -0.011670    1( a  )  50( a  )    -   +-   +-   +-   +-         - 
 x   1  1  234757  0.010869    1( a  )  49( a  )    -   +-    -   +-   +-   +     - 
 w   1  1  287311  0.013419   15( a  )  46( a  )   +-   +-   +-   +-   +     -      
 w   1  1  294474 -0.011257    5( a  )  15( a  )   +-   +-   +-   +-        +     - 
 w   1  1  294476  0.010479    7( a  )  15( a  )   +-   +-   +-   +-        +     - 
 w   1  1  294477 -0.012407    8( a  )  15( a  )   +-   +-   +-   +-        +     - 
 w   1  1  294480  0.011779   11( a  )  15( a  )   +-   +-   +-   +-        +     - 
 w   1  1  294481 -0.016475   12( a  )  15( a  )   +-   +-   +-   +-        +     - 
 w   1  1  294484 -0.030262   15( a  )  15( a  )   +-   +-   +-   +-        +     - 
 w   1  1  294785  0.010395   15( a  )  29( a  )   +-   +-   +-   +-        +     - 
 w   1  1  294814 -0.014043   15( a  )  30( a  )   +-   +-   +-   +-        +     - 
 w   1  1  294974 -0.024212   15( a  )  35( a  )   +-   +-   +-   +-        +     - 
 w   1  1  294994 -0.015384   35( a  )  35( a  )   +-   +-   +-   +-        +     - 
 w   1  1  295325 -0.010016   15( a  )  44( a  )   +-   +-   +-   +-        +     - 
 w   1  1  295369 -0.015679   15( a  )  45( a  )   +-   +-   +-   +-        +     - 
 w   1  1  296169 -0.010534   35( a  )  60( a  )   +-   +-   +-   +-        +     - 
 w   1  1  300893 -0.011429   46( a  )  47( a  )   +-   +-   +-   +    +-    -      
 w   1  1  394302 -0.015077    1( a  )   1( a  )   +-   +-        +-   +-   +     - 
 w   1  1  394337 -0.011033    8( a  )   8( a  )   +-   +-        +-   +-   +     - 
 w   1  1  394400  0.011593    8( a  )  14( a  )   +-   +-        +-   +-   +     - 
 w   1  1  395196 -0.010361   34( a  )  42( a  )   +-   +-        +-   +-   +     - 
 w   1  1  395204 -0.018005   42( a  )  42( a  )   +-   +-        +-   +-   +     - 
 w   1  1  395526 -0.011805   49( a  )  49( a  )   +-   +-        +-   +-   +     - 
 w   1  1  426881 -0.013253   15( a  )  18( a  )   +-   +    +-   +-    -   +     - 
 w   1  1  427134 -0.011328   15( a  )  29( a  )   +-   +    +-   +-    -   +     - 
 w   1  1  427674  0.011080   15( a  )  44( a  )   +-   +    +-   +-    -   +     - 
 w   1  1  502951 -0.011726   15( a  )  35( a  )   +-        +-   +-   +-   +     - 
 w   1  1  502971 -0.011210   35( a  )  35( a  )   +-        +-   +-   +-   +     - 
 w   1  1  529352  0.014998    1( a  )   1( a  )   +    +-   +-   +-   +-         - 
 w   1  1  529605 -0.010896    1( a  )  23( a  )   +    +-   +-   +-   +-         - 
 w   1  1  530577  0.020896    1( a  )  50( a  )   +    +-   +-   +-   +-         - 
 w   1  1  583492 -0.013954    1( a  )  16( a  )   +    +-    -   +-   +-   +     - 
 w   1  1  583900 -0.015300    1( a  )  33( a  )   +    +-    -   +-   +-   +     - 
 w   1  1  584240  0.013086    8( a  )  42( a  )   +    +-    -   +-   +-   +     - 
 w   1  1  584261 -0.010873   29( a  )  42( a  )   +    +-    -   +-   +-   +     - 
 w   1  1  584316 -0.011088   42( a  )  43( a  )   +    +-    -   +-   +-   +     - 
 w   1  1  584548 -0.027531    1( a  )  49( a  )   +    +-    -   +-   +-   +     - 
 w   1  1  584739 -0.011134   42( a  )  52( a  )   +    +-    -   +-   +-   +     - 
 w   1  1  637392 -0.046703    1( a  )   1( a  )        +-   +-   +-   +-   +     - 
 w   1  1  637407  0.013714    1( a  )   6( a  )        +-   +-   +-   +-   +     - 
 w   1  1  637413 -0.020531    1( a  )   7( a  )        +-   +-   +-   +-   +     - 
 w   1  1  637717 -0.012077    1( a  )  26( a  )        +-   +-   +-   +-   +     - 
 w   1  1  638058 -0.013026    1( a  )  37( a  )        +-   +-   +-   +-   +     - 

 ci coefficient statistics:
           rq > 0.1                1
      0.1> rq > 0.01             105
     0.01> rq > 0.001          10094
    0.001> rq > 0.0001         42994
   0.0001> rq > 0.00001       139752
  0.00001> rq > 0.000001      229950
 0.000001> rq                 238804
           all                661700
 ========= Executing OUT-OF-CORE method ========


====================================================================================================
Diagonal     counts:  0x:         735 2x:           0 4x:           0
All internal counts: zz :        1113 yy:           0 xx:           0 ww:           0
One-external counts: yz :           0 yx:           0 yw:           0
Two-external counts: yy :           0 ww:           0 xx:           0 xz:           0 wz:           0 wx:           0
Three-ext.   counts: yx :           0 yw:           0

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================


  iref  icsf         v(icsf)             hv(icsf)
     1     1     -0.031738769911      2.986219791897
     2     2      0.946292736709    -89.028760296632
     3     3      0.018920907233     -1.779480434946
     4     4     -0.009155881546      0.861471559757
     5     5      0.002594361026     -0.244117660783
     6     6      0.000635165725     -0.059704465625
     7     7      0.000589998989     -0.055699555190
     8     8     -0.003514609068      0.330217340417
     9     9     -0.000179901869      0.016934910754
    10    10      0.001231392511     -0.115684860213

 number of reference csfs (nref) is    10.  root number (iroot) is  1.
 c0**2 =   0.89694051  c**2 (all zwalks) =   0.89695348

 pople ci energy extrapolation is computed with 12 correlated electrons.

 eref      =    -94.081619227920   "relaxed" cnot**2         =   0.896940507373
 eci       =    -94.361419041214   deltae = eci - eref       =  -0.279799813294
 eci+dv1   =    -94.390255068009   dv1 = (1-cnot**2)*deltae  =  -0.028836026795
 eci+dv2   =    -94.393568360548   dv2 = dv1 / cnot**2       =  -0.032149319334
 eci+dv3   =    -94.397741898493   dv3 = dv1 / (2*cnot**2-1) =  -0.036322857279
 eci+pople =    -94.390370475136   ( 12e- scaled deltae )    =  -0.308751247216


 --- list of ci coefficients ( ctol =   1.00E-02 )  total energy( 2) =       -94.2184643473

                                                       internal orbitals

                                          level       1    2    3    4    5    6    7

                                          orbital     3    4    5    6    7    8    9

                                         symmetry   a    a    a    a    a    a    a  

 path  s ms    csf#    c(i)    ext. orb.(sym)
 z*  1  1       1 -0.030512                        +-   +-   +-   +-   +-   +-      
 z*  1  1       3  0.019911                        +-   +-   +-   +-   +-        +- 
 z*  1  1       4 -0.025784                        +-   +-   +-   +-   +    +-    - 
 z*  1  1       5  0.945755                        +-   +-   +-   +-   +     -   +- 
 z*  1  1       6 -0.011939                        +-   +-   +-   +-        +-   +- 
 z   1  1      17 -0.010976                        +-   +    +-   +-   +-    -   +- 
 y   1  1      43 -0.015058             15( a  )   +-   +-   +-   +-   +-    -      
 y   1  1     258  0.011083             11( a  )   +-   +-   +-   +-    -   +     - 
 y   1  1     326 -0.034025              6( a  )   +-   +-   +-   +-    -        +- 
 y   1  1     327 -0.021178              7( a  )   +-   +-   +-   +-    -        +- 
 y   1  1     342 -0.010109             22( a  )   +-   +-   +-   +-    -        +- 
 y   1  1     343  0.010027             23( a  )   +-   +-   +-   +-    -        +- 
 y   1  1     345  0.021875             25( a  )   +-   +-   +-   +-    -        +- 
 y   1  1     347  0.011618             27( a  )   +-   +-   +-   +-    -        +- 
 y   1  1     348  0.014227             28( a  )   +-   +-   +-   +-    -        +- 
 y   1  1     356  0.016011             36( a  )   +-   +-   +-   +-    -        +- 
 y   1  1     358  0.017967             38( a  )   +-   +-   +-   +-    -        +- 
 y   1  1     361  0.028458             41( a  )   +-   +-   +-   +-    -        +- 
 y   1  1     541 -0.010206              2( a  )   +-   +-   +-   +-         -   +- 
 y   1  1     542  0.011580              3( a  )   +-   +-   +-   +-         -   +- 
 y   1  1     543  0.010030              4( a  )   +-   +-   +-   +-         -   +- 
 y   1  1     544 -0.011483              5( a  )   +-   +-   +-   +-         -   +- 
 y   1  1     547 -0.010937              8( a  )   +-   +-   +-   +-         -   +- 
 y   1  1     550  0.013727             11( a  )   +-   +-   +-   +-         -   +- 
 y   1  1     551 -0.015875             12( a  )   +-   +-   +-   +-         -   +- 
 y   1  1     554 -0.035965             15( a  )   +-   +-   +-   +-         -   +- 
 y   1  1     569 -0.012340             30( a  )   +-   +-   +-   +-         -   +- 
 y   1  1     574 -0.011040             35( a  )   +-   +-   +-   +-         -   +- 
 y   1  1     584 -0.014691             45( a  )   +-   +-   +-   +-         -   +- 
 y   1  1     906 -0.011569              2( a  )   +-   +-   +-    -   +     -   +- 
 y   1  1     909 -0.012093              5( a  )   +-   +-   +-    -   +     -   +- 
 y   1  1     916  0.015601             12( a  )   +-   +-   +-    -   +     -   +- 
 y   1  1     917  0.020862             13( a  )   +-   +-   +-    -   +     -   +- 
 y   1  1     918 -0.015854             14( a  )   +-   +-   +-    -   +     -   +- 
 y   1  1     921 -0.017958             17( a  )   +-   +-   +-    -   +     -   +- 
 y   1  1     948  0.013521             44( a  )   +-   +-   +-    -   +     -   +- 
 y   1  1    1069 -0.012747             19( a  )   +-   +-   +-   +    +-    -    - 
 y   1  1    1071  0.010941             21( a  )   +-   +-   +-   +    +-    -    - 
 y   1  1    1087 -0.012453             37( a  )   +-   +-   +-   +    +-    -    - 
 y   1  1    1089 -0.016429             39( a  )   +-   +-   +-   +    +-    -    - 
 y   1  1    1090 -0.011825             40( a  )   +-   +-   +-   +    +-    -    - 
 y   1  1    1117 -0.019357             67( a  )   +-   +-   +-   +    +-    -    - 
 y   1  1    1208 -0.011997             12( a  )   +-   +-   +-   +     -    -   +- 
 y   1  1    1211 -0.024362             15( a  )   +-   +-   +-   +     -    -   +- 
 y   1  1    1231 -0.026151             35( a  )   +-   +-   +-   +     -    -   +- 
 y   1  1    1256 -0.014220             60( a  )   +-   +-   +-   +     -    -   +- 
 y   1  1    1354  0.011180             12( a  )   +-   +-   +-        +-    -   +- 
 y   1  1    1357  0.023082             15( a  )   +-   +-   +-        +-    -   +- 
 y   1  1    1377  0.011058             35( a  )   +-   +-   +-        +-    -   +- 
 y   1  1    1822  0.011810             42( a  )   +-   +-    -   +-   +     -   +- 
 y   1  1    3243 -0.014103              3( a  )   +-    -   +-   +-   +     -   +- 
 y   1  1    3245 -0.010866              5( a  )   +-    -   +-   +-   +     -   +- 
 y   1  1    3251 -0.014267             11( a  )   +-    -   +-   +-   +     -   +- 
 y   1  1    3252  0.018288             12( a  )   +-    -   +-   +-   +     -   +- 
 y   1  1    3257 -0.016154             17( a  )   +-    -   +-   +-   +     -   +- 
 y   1  1    3258 -0.012529             18( a  )   +-    -   +-   +-   +     -   +- 
 y   1  1    3269 -0.010099             29( a  )   +-    -   +-   +-   +     -   +- 
 y   1  1    3284  0.012559             44( a  )   +-    -   +-   +-   +     -   +- 
 y   1  1    4131  0.018698             15( a  )   +-   +    +-   +-    -    -   +- 
 y   1  1    4151  0.012364             35( a  )   +-   +    +-   +-    -    -   +- 
 y   1  1    4161  0.010437             45( a  )   +-   +    +-   +-    -    -   +- 
 y   1  1    4808 -0.012390             35( a  )   +-        +-   +-   +-    -   +- 
 y   1  1    6460  0.018695              8( a  )   +    +-   +-   +-    -    -   +- 
 y   1  1    6463  0.010039             11( a  )   +    +-   +-   +-    -    -   +- 
 y   1  1    6466 -0.012173             14( a  )   +    +-   +-   +-    -    -   +- 
 y   1  1    6471 -0.010329             19( a  )   +    +-   +-   +-    -    -   +- 
 y   1  1    6474  0.015243             22( a  )   +    +-   +-   +-    -    -   +- 
 y   1  1    6503  0.028645             51( a  )   +    +-   +-   +-    -    -   +- 
 y   1  1    6924 -0.013225             34( a  )   +    +-    -   +-    -   +-   +- 
 y   1  1    6938  0.016313             48( a  )   +    +-    -   +-    -   +-   +- 
 x   1  1   13044  0.011755   15( a  )  19( a  )   +-   +-   +-   +-         -    - 
 x   1  1   13881 -0.014498   15( a  )  46( a  )   +-   +-   +-   +-         -    - 
 x   1  1  113774  0.011159   44( a  )  46( a  )   +-    -   +-   +-   +     -    - 
 x   1  1  190129  0.011665    1( a  )  50( a  )    -   +-   +-   +-    -        +- 
 x   1  1  242641  0.010859    1( a  )  49( a  )    -   +-    -   +-   +     -   +- 
 w   1  1  286927 -0.012870   36( a  )  36( a  )   +-   +-   +-   +-   +     -      
 w   1  1  287305  0.012197    9( a  )  46( a  )   +-   +-   +-   +-   +     -      
 w   1  1  287307 -0.015723   11( a  )  46( a  )   +-   +-   +-   +-   +     -      
 w   1  1  287315  0.014578   19( a  )  46( a  )   +-   +-   +-   +-   +     -      
 w   1  1  287342 -0.031801   46( a  )  46( a  )   +-   +-   +-   +-   +     -      
 w   1  1  288127  0.011018   36( a  )  61( a  )   +-   +-   +-   +-   +     -      
 w   1  1  295414  0.012161   15( a  )  46( a  )   +-   +-   +-   +-        +     - 
 w   1  1  314319 -0.010163   13( a  )  46( a  )   +-   +-   +-   +    +     -    - 
 w   1  1  314398 -0.013921   46( a  )  47( a  )   +-   +-   +-   +    +     -    - 
 w   1  1  332298 -0.011497   15( a  )  15( a  )   +-   +-   +-        +     -   +- 
 w   1  1  332788 -0.011256   15( a  )  35( a  )   +-   +-   +-        +     -   +- 
 w   1  1  402405 -0.015056    1( a  )   1( a  )   +-   +-        +-   +     -   +- 
 w   1  1  402440 -0.011035    8( a  )   8( a  )   +-   +-        +-   +     -   +- 
 w   1  1  402503  0.011595    8( a  )  14( a  )   +-   +-        +-   +     -   +- 
 w   1  1  403299 -0.010355   34( a  )  42( a  )   +-   +-        +-   +     -   +- 
 w   1  1  403307 -0.017992   42( a  )  42( a  )   +-   +-        +-   +     -   +- 
 w   1  1  403629 -0.011797   49( a  )  49( a  )   +-   +-        +-   +     -   +- 
 w   1  1  433194 -0.010065   44( a  )  46( a  )   +-   +    +-   +-   +     -    - 
 w   1  1  511054 -0.012575   15( a  )  35( a  )   +-        +-   +-   +     -   +- 
 w   1  1  511074 -0.011557   35( a  )  35( a  )   +-        +-   +-   +     -   +- 
 w   1  1  537455 -0.014972    1( a  )   1( a  )   +    +-   +-   +-    -        +- 
 w   1  1  537708  0.010881    1( a  )  23( a  )   +    +-   +-   +-    -        +- 
 w   1  1  538680 -0.020871    1( a  )  50( a  )   +    +-   +-   +-    -        +- 
 w   1  1  591595 -0.013939    1( a  )  16( a  )   +    +-    -   +-   +     -   +- 
 w   1  1  592003 -0.015282    1( a  )  33( a  )   +    +-    -   +-   +     -   +- 
 w   1  1  592343  0.013090    8( a  )  42( a  )   +    +-    -   +-   +     -   +- 
 w   1  1  592364 -0.010868   29( a  )  42( a  )   +    +-    -   +-   +     -   +- 
 w   1  1  592419 -0.011088   42( a  )  43( a  )   +    +-    -   +-   +     -   +- 
 w   1  1  592651 -0.027506    1( a  )  49( a  )   +    +-    -   +-   +     -   +- 
 w   1  1  592842 -0.011127   42( a  )  52( a  )   +    +-    -   +-   +     -   +- 
 w   1  1  645495 -0.046650    1( a  )   1( a  )        +-   +-   +-   +     -   +- 
 w   1  1  645510  0.013702    1( a  )   6( a  )        +-   +-   +-   +     -   +- 
 w   1  1  645516 -0.020508    1( a  )   7( a  )        +-   +-   +-   +     -   +- 
 w   1  1  645820 -0.012063    1( a  )  26( a  )        +-   +-   +-   +     -   +- 
 w   1  1  646161 -0.013010    1( a  )  37( a  )        +-   +-   +-   +     -   +- 

 ci coefficient statistics:
           rq > 0.1                1
      0.1> rq > 0.01             109
     0.01> rq > 0.001          10413
    0.001> rq > 0.0001         52017
   0.0001> rq > 0.00001       173836
  0.00001> rq > 0.000001      241748
 0.000001> rq                 183576
           all                661700
 ========= Executing OUT-OF-CORE method ========


====================================================================================================
Diagonal     counts:  0x:         735 2x:           0 4x:           0
All internal counts: zz :        1113 yy:           0 xx:           0 ww:           0
One-external counts: yz :           0 yx:           0 yw:           0
Two-external counts: yy :           0 ww:           0 xx:           0 xz:           0 wz:           0 wx:           0
Three-ext.   counts: yx :           0 yw:           0

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================


  iref  icsf         v(icsf)             hv(icsf)
     1     1     -0.030511856946      2.865454323940
     2     2     -0.004255160799      0.399646202668
     3     3      0.019910631113     -1.870507682821
     4     4     -0.025783658549      2.421932703493
     5     5      0.945754951600    -88.842695739679
     6     6     -0.011938914218      1.121561477727
     7     7      0.000089945652     -0.008110397652
     8     8     -0.005869460921      0.546314908141
     9     9     -0.001949500940      0.182644968190
    10    10      0.003188613121     -0.298815738229

 number of reference csfs (nref) is    10.  root number (iroot) is  2.
 c0**2 =   0.89665370  c**2 (all zwalks) =   0.89678412

 pople ci energy extrapolation is computed with 12 correlated electrons.

 eref      =    -93.938323931345   "relaxed" cnot**2         =   0.896653702706
 eci       =    -94.218464347297   deltae = eci - eref       =  -0.280140415952
 eci+dv1   =    -94.247415822008   dv1 = (1-cnot**2)*deltae  =  -0.028951474711
 eci+dv2   =    -94.250752704133   dv2 = dv1 / cnot**2       =  -0.032288356836
 eci+dv3   =    -94.254958995542   dv3 = dv1 / (2*cnot**2-1) =  -0.036494648246
 eci+pople =    -94.247548058270   ( 12e- scaled deltae )    =  -0.309224126925


 --- list of ci coefficients ( ctol =   1.00E-02 )  total energy( 3) =       -94.2072696437

                                                       internal orbitals

                                          level       1    2    3    4    5    6    7

                                          orbital     3    4    5    6    7    8    9

                                         symmetry   a    a    a    a    a    a    a  

 path  s ms    csf#    c(i)    ext. orb.(sym)
 z*  1  1       1  0.925983                        +-   +-   +-   +-   +-   +-      
 z*  1  1       2  0.030844                        +-   +-   +-   +-   +-   +     - 
 z*  1  1       5  0.030894                        +-   +-   +-   +-   +     -   +- 
 z*  1  1       6 -0.062933                        +-   +-   +-   +-        +-   +- 
 z*  1  1       8  0.029083                        +-   +-   +-   +    +-    -   +- 
 z*  1  1       9  0.033912                        +-   +-   +-   +     -   +-   +- 
 z*  1  1      10 -0.111830                        +-   +-   +-        +-   +-   +- 
 z   1  1      17  0.018970                        +-   +    +-   +-   +-    -   +- 
 z   1  1      19 -0.023633                        +-   +    +-    -   +-   +-   +- 
 z   1  1      21 -0.033580                        +-        +-   +-   +-   +-   +- 
 y   1  1      30 -0.013178              2( a  )   +-   +-   +-   +-   +-    -      
 y   1  1      32  0.016378              4( a  )   +-   +-   +-   +-   +-    -      
 y   1  1      34 -0.062615              6( a  )   +-   +-   +-   +-   +-    -      
 y   1  1      35 -0.039208              7( a  )   +-   +-   +-   +-   +-    -      
 y   1  1      37  0.018830              9( a  )   +-   +-   +-   +-   +-    -      
 y   1  1      38  0.010193             10( a  )   +-   +-   +-   +-   +-    -      
 y   1  1      40  0.010664             12( a  )   +-   +-   +-   +-   +-    -      
 y   1  1      42 -0.013475             14( a  )   +-   +-   +-   +-   +-    -      
 y   1  1      45  0.014300             17( a  )   +-   +-   +-   +-   +-    -      
 y   1  1      47  0.014115             19( a  )   +-   +-   +-   +-   +-    -      
 y   1  1      49  0.016356             21( a  )   +-   +-   +-   +-   +-    -      
 y   1  1      50 -0.020168             22( a  )   +-   +-   +-   +-   +-    -      
 y   1  1      51  0.019528             23( a  )   +-   +-   +-   +-   +-    -      
 y   1  1      53  0.027780             25( a  )   +-   +-   +-   +-   +-    -      
 y   1  1      55  0.015104             27( a  )   +-   +-   +-   +-   +-    -      
 y   1  1      56  0.015354             28( a  )   +-   +-   +-   +-   +-    -      
 y   1  1      64  0.014611             36( a  )   +-   +-   +-   +-   +-    -      
 y   1  1      66  0.018308             38( a  )   +-   +-   +-   +-   +-    -      
 y   1  1      69  0.027088             41( a  )   +-   +-   +-   +-   +-    -      
 y   1  1      79 -0.015755             51( a  )   +-   +-   +-   +-   +-    -      
 y   1  1     178 -0.015151              4( a  )   +-   +-   +-   +-    -   +-      
 y   1  1     179  0.014697              5( a  )   +-   +-   +-   +-    -   +-      
 y   1  1     185 -0.012834             11( a  )   +-   +-   +-   +-    -   +-      
 y   1  1     186  0.013726             12( a  )   +-   +-   +-   +-    -   +-      
 y   1  1     188 -0.012744             14( a  )   +-   +-   +-   +-    -   +-      
 y   1  1     189  0.030776             15( a  )   +-   +-   +-   +-    -   +-      
 y   1  1     204  0.010046             30( a  )   +-   +-   +-   +-    -   +-      
 y   1  1     206  0.021892             32( a  )   +-   +-   +-   +-    -   +-      
 y   1  1     218  0.013696             44( a  )   +-   +-   +-   +-    -   +-      
 y   1  1     219  0.017440             45( a  )   +-   +-   +-   +-    -   +-      
 y   1  1     229 -0.015916             55( a  )   +-   +-   +-   +-    -   +-      
 y   1  1     614  0.018787              2( a  )   +-   +-   +-    -   +-   +-      
 y   1  1     617  0.010361              5( a  )   +-   +-   +-    -   +-   +-      
 y   1  1     620  0.014338              8( a  )   +-   +-   +-    -   +-   +-      
 y   1  1     621 -0.013815              9( a  )   +-   +-   +-    -   +-   +-      
 y   1  1     624 -0.012514             12( a  )   +-   +-   +-    -   +-   +-      
 y   1  1     625 -0.026173             13( a  )   +-   +-   +-    -   +-   +-      
 y   1  1     626  0.014833             14( a  )   +-   +-   +-    -   +-   +-      
 y   1  1     629  0.010402             17( a  )   +-   +-   +-    -   +-   +-      
 y   1  1     635 -0.010723             23( a  )   +-   +-   +-    -   +-   +-      
 y   1  1     641 -0.010173             29( a  )   +-   +-   +-    -   +-   +-      
 y   1  1     648  0.011264             36( a  )   +-   +-   +-    -   +-   +-      
 y   1  1     656 -0.020373             44( a  )   +-   +-   +-    -   +-   +-      
 y   1  1     657  0.011346             45( a  )   +-   +-   +-    -   +-   +-      
 y   1  1     659 -0.018627             47( a  )   +-   +-   +-    -   +-   +-      
 y   1  1    1142  0.014508             19( a  )   +-   +-   +-   +     -   +-    - 
 y   1  1    1144 -0.012814             21( a  )   +-   +-   +-   +     -   +-    - 
 y   1  1    1160  0.013275             37( a  )   +-   +-   +-   +     -   +-    - 
 y   1  1    1162  0.020095             39( a  )   +-   +-   +-   +     -   +-    - 
 y   1  1    1163  0.015727             40( a  )   +-   +-   +-   +     -   +-    - 
 y   1  1    1169 -0.010475             46( a  )   +-   +-   +-   +     -   +-    - 
 y   1  1    1179 -0.012724             56( a  )   +-   +-   +-   +     -   +-    - 
 y   1  1    1190  0.028398             67( a  )   +-   +-   +-   +     -   +-    - 
 y   1  1    1278 -0.012846              9( a  )   +-   +-   +-        +-   +-    - 
 y   1  1    1279  0.010697             10( a  )   +-   +-   +-        +-   +-    - 
 y   1  1    1280  0.016930             11( a  )   +-   +-   +-        +-   +-    - 
 y   1  1    1285  0.010196             16( a  )   +-   +-   +-        +-   +-    - 
 y   1  1    1288 -0.020619             19( a  )   +-   +-   +-        +-   +-    - 
 y   1  1    1291 -0.013380             22( a  )   +-   +-   +-        +-   +-    - 
 y   1  1    1315  0.030144             46( a  )   +-   +-   +-        +-   +-    - 
 y   1  1    1337 -0.019775             68( a  )   +-   +-   +-        +-   +-    - 
 y   1  1    1490 -0.013031              2( a  )   +-   +-    -   +-   +-   +-      
 y   1  1    1492 -0.025152              4( a  )   +-   +-    -   +-   +-   +-      
 y   1  1    1493 -0.011660              5( a  )   +-   +-    -   +-   +-   +-      
 y   1  1    1496 -0.022122              8( a  )   +-   +-    -   +-   +-   +-      
 y   1  1    1497  0.013644              9( a  )   +-   +-    -   +-   +-   +-      
 y   1  1    1500  0.025216             12( a  )   +-   +-    -   +-   +-   +-      
 y   1  1    1502  0.020562             14( a  )   +-   +-    -   +-   +-   +-      
 y   1  1    1504 -0.020655             16( a  )   +-   +-    -   +-   +-   +-      
 y   1  1    1505 -0.010759             17( a  )   +-   +-    -   +-   +-   +-      
 y   1  1    1506  0.010616             18( a  )   +-   +-    -   +-   +-   +-      
 y   1  1    1511  0.017239             23( a  )   +-   +-    -   +-   +-   +-      
 y   1  1    1512  0.024988             24( a  )   +-   +-    -   +-   +-   +-      
 y   1  1    1522 -0.012203             34( a  )   +-   +-    -   +-   +-   +-      
 y   1  1    1530 -0.024755             42( a  )   +-   +-    -   +-   +-   +-      
 y   1  1    2950 -0.010445              2( a  )   +-    -   +-   +-   +-   +-      
 y   1  1    2951  0.017683              3( a  )   +-    -   +-   +-   +-   +-      
 y   1  1    2953  0.011535              5( a  )   +-    -   +-   +-   +-   +-      
 y   1  1    2959  0.018137             11( a  )   +-    -   +-   +-   +-   +-      
 y   1  1    2960 -0.025057             12( a  )   +-    -   +-   +-   +-   +-      
 y   1  1    2965  0.014299             17( a  )   +-    -   +-   +-   +-   +-      
 y   1  1    2966  0.014697             18( a  )   +-    -   +-   +-   +-   +-      
 y   1  1    2977  0.016486             29( a  )   +-    -   +-   +-   +-   +-      
 y   1  1    2991  0.011150             43( a  )   +-    -   +-   +-   +-   +-      
 y   1  1    2992 -0.026123             44( a  )   +-    -   +-   +-   +-   +-      
 y   1  1    2993  0.010276             45( a  )   +-    -   +-   +-   +-   +-      
 y   1  1    3424  0.013025             38( a  )   +-    -   +-   +    +-   +-    - 
 y   1  1    3454 -0.016984             68( a  )   +-    -   +-   +    +-   +-    - 
 y   1  1    4257 -0.010996             68( a  )   +-   +    +-    -   +-   +-    - 
 y   1  1    4757  0.010756             57( a  )   +-        +-   +-   +-   +-    - 
 y   1  1    4994  0.019712              2( a  )    -   +-   +-   +-   +-   +-      
 y   1  1    4995  0.017963              3( a  )    -   +-   +-   +-   +-   +-      
 y   1  1    5000  0.012552              8( a  )    -   +-   +-   +-   +-   +-      
 y   1  1    5001 -0.010593              9( a  )    -   +-   +-   +-   +-   +-      
 y   1  1    5002 -0.021168             10( a  )    -   +-   +-   +-   +-   +-      
 y   1  1    5003  0.013056             11( a  )    -   +-   +-   +-   +-   +-      
 y   1  1    5017  0.010657             25( a  )    -   +-   +-   +-   +-   +-      
 y   1  1    5019 -0.011501             27( a  )    -   +-   +-   +-   +-   +-      
 y   1  1    5020 -0.011987             28( a  )    -   +-   +-   +-   +-   +-      
 y   1  1    5044 -0.013272             52( a  )    -   +-   +-   +-   +-   +-      
 x   1  1  179617 -0.014071    1( a  )  50( a  )    -   +-   +-   +-   +-    -      
 w   1  1  283561 -0.066723    1( a  )   1( a  )   +-   +-   +-   +-   +-           
 w   1  1  283576  0.012879    1( a  )   6( a  )   +-   +-   +-   +-   +-           
 w   1  1  283581 -0.011206    6( a  )   6( a  )   +-   +-   +-   +-   +-           
 w   1  1  283582 -0.019169    1( a  )   7( a  )   +-   +-   +-   +-   +-           
 w   1  1  283814  0.013679    1( a  )  23( a  )   +-   +-   +-   +-   +-           
 w   1  1  284386  0.013019    6( a  )  41( a  )   +-   +-   +-   +-   +-           
 w   1  1  284416 -0.010892   36( a  )  41( a  )   +-   +-   +-   +-   +-           
 w   1  1  284418 -0.012878   38( a  )  41( a  )   +-   +-   +-   +-   +-           
 w   1  1  284421 -0.014804   41( a  )  41( a  )   +-   +-   +-   +-   +-           
 w   1  1  284786 -0.014857    1( a  )  50( a  )   +-   +-   +-   +-   +-           
 w   1  1  284835 -0.011734   50( a  )  50( a  )   +-   +-   +-   +-   +-           
 w   1  1  284886 -0.014591   51( a  )  51( a  )   +-   +-   +-   +-   +-           
 w   1  1  291776 -0.010689    8( a  )  15( a  )   +-   +-   +-   +-        +-      
 w   1  1  291780 -0.013327   12( a  )  15( a  )   +-   +-   +-   +-        +-      
 w   1  1  291783 -0.024855   15( a  )  15( a  )   +-   +-   +-   +-        +-      
 w   1  1  292084  0.010450   15( a  )  29( a  )   +-   +-   +-   +-        +-      
 w   1  1  292113 -0.013308   15( a  )  30( a  )   +-   +-   +-   +-        +-      
 w   1  1  292273 -0.022763   15( a  )  35( a  )   +-   +-   +-   +-        +-      
 w   1  1  292293 -0.014452   35( a  )  35( a  )   +-   +-   +-   +-        +-      
 w   1  1  292624 -0.010130   15( a  )  44( a  )   +-   +-   +-   +-        +-      
 w   1  1  292668 -0.015759   15( a  )  45( a  )   +-   +-   +-   +-        +-      
 w   1  1  293468 -0.011891   35( a  )  60( a  )   +-   +-   +-   +-        +-      
 w   1  1  391601 -0.012528    1( a  )   1( a  )   +-   +-        +-   +-   +-      
 w   1  1  391636 -0.012858    8( a  )   8( a  )   +-   +-        +-   +-   +-      
 w   1  1  391663 -0.010294    8( a  )  11( a  )   +-   +-        +-   +-   +-      
 w   1  1  391699  0.012936    8( a  )  14( a  )   +-   +-        +-   +-   +-      
 w   1  1  392503 -0.014992   42( a  )  42( a  )   +-   +-        +-   +-   +-      
 w   1  1  392825 -0.010981   49( a  )  49( a  )   +-   +-        +-   +-   +-      
 w   1  1  424180 -0.011028   15( a  )  18( a  )   +-   +    +-   +-    -   +-      
 w   1  1  424433 -0.010127   15( a  )  29( a  )   +-   +    +-   +-    -   +-      
 w   1  1  424973  0.011316   15( a  )  44( a  )   +-   +    +-   +-    -   +-      
 w   1  1  500250 -0.011014   15( a  )  35( a  )   +-        +-   +-   +-   +-      
 w   1  1  500270 -0.010450   35( a  )  35( a  )   +-        +-   +-   +-   +-      
 w   1  1  526651  0.016986    1( a  )   1( a  )   +    +-   +-   +-   +-    -      
 w   1  1  526804  0.011246    1( a  )  18( a  )   +    +-   +-   +-   +-    -      
 w   1  1  526841 -0.012791    1( a  )  20( a  )   +    +-   +-   +-   +-    -      
 w   1  1  526882 -0.011516    1( a  )  22( a  )   +    +-   +-   +-   +-    -      
 w   1  1  526904 -0.016508    1( a  )  23( a  )   +    +-   +-   +-   +-    -      
 w   1  1  527876  0.027458    1( a  )  50( a  )   +    +-   +-   +-   +-    -      
 w   1  1  580791 -0.014870    1( a  )  16( a  )   +    +-    -   +-   +-   +-      
 w   1  1  581199 -0.014856    1( a  )  33( a  )   +    +-    -   +-   +-   +-      
 w   1  1  581539  0.013092    8( a  )  42( a  )   +    +-    -   +-   +-   +-      
 w   1  1  581847 -0.025050    1( a  )  49( a  )   +    +-    -   +-   +-   +-      
 w   1  1  582038 -0.010070   42( a  )  52( a  )   +    +-    -   +-   +-   +-      
 w   1  1  634691 -0.039647    1( a  )   1( a  )        +-   +-   +-   +-   +-      
 w   1  1  634706  0.011264    1( a  )   6( a  )        +-   +-   +-   +-   +-      
 w   1  1  634712 -0.016785    1( a  )   7( a  )        +-   +-   +-   +-   +-      

 ci coefficient statistics:
           rq > 0.1                2
      0.1> rq > 0.01             156
     0.01> rq > 0.001           9097
    0.001> rq > 0.0001         58015
   0.0001> rq > 0.00001       195636
  0.00001> rq > 0.000001      242268
 0.000001> rq                 156526
           all                661700
 ========= Executing OUT-OF-CORE method ========


====================================================================================================
Diagonal     counts:  0x:         735 2x:           0 4x:           0
All internal counts: zz :        1113 yy:           0 xx:           0 ww:           0
One-external counts: yz :           0 yx:           0 yw:           0
Two-external counts: yy :           0 ww:           0 xx:           0 xz:           0 wz:           0 wx:           0
Three-ext.   counts: yx :           0 yw:           0

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================


  iref  icsf         v(icsf)             hv(icsf)
     1     1      0.925982583828    -86.953867768047
     2     2      0.030843550277     -2.896078131630
     3     3      0.000319517309     -0.030123607513
     4     4     -0.000545712515      0.051330284519
     5     5      0.030893599200     -2.902038708624
     6     6     -0.062933107801      5.907156865919
     7     7      0.001095785921     -0.103087371240
     8     8      0.029082557393     -2.725307690272
     9     9      0.033912148293     -3.181800384929
    10    10     -0.111829973374     10.484985182263

 number of reference csfs (nref) is    10.  root number (iroot) is  3.
 c0**2 =   0.87781343  c**2 (all zwalks) =   0.87987556

 pople ci energy extrapolation is computed with 12 correlated electrons.

 eref      =    -93.901900266531   "relaxed" cnot**2         =   0.877813433206
 eci       =    -94.207269643688   deltae = eci - eref       =  -0.305369377156
 eci+dv1   =    -94.244581679486   dv1 = (1-cnot**2)*deltae  =  -0.037312035799
 eci+dv2   =    -94.249775299659   dv2 = dv1 / cnot**2       =  -0.042505655971
 eci+dv3   =    -94.256648559998   dv3 = dv1 / (2*cnot**2-1) =  -0.049378916311
 eci+pople =    -94.246187375507   ( 12e- scaled deltae )    =  -0.344287108976
 passed aftci ... 
 readint2: molcas,dalton2=                     0                     0
 files%faoints=aoints              
                       Size (real*8) of d2temp for two-external contributions     206444
 bummer (warning):inf%bufszi too small ... resetting to                    8030
 
                       Size (real*8) of d2temp for all-internal contributions        546
                       Size (real*8) of d2temp for one-external contributions      18396
                       Size (real*8) of d2temp for two-external contributions     206444
size_thrext:  lsym   l1    ksym   k1strt   k1       cnt3 
                1    1    1    1   74   202356
                1    2    1    1   74   202356
                1    3    1    1   74   202356
                1    4    1    1   74   202356
                1    5    1    1   74   202356
                1    6    1    1   74   202356
                1    7    1    1   74   202356
                       Size (real*8) of d2temp for three-external contributions    1416492
                       Size (real*8) of d2temp for four-external contributions    3848925
 serial operation: forcing vdisk for temporary dd012 matrices
location of d2temp files... fileloc(dd012)=       1
location of d2temp files... fileloc(d3)=       1
location of d2temp files... fileloc(d4)=       1
 files%dd012ext =  unit=  22  vdsk=   1  filestart=       1
 files%d3ext =     unit=  23  vdsk=   1  filestart=  248931
 files%d4ext =     unit=  24  vdsk=   1  filestart= 1688931
            0xdiag    0ext      1ext      2ext      3ext      4ext
d2off                  8031     16061     40151         1         1
d2rec                     1         3        26        24        65
recsize                8030      8030      8030     60000     60000
d2bufferlen=          60000
maxbl3=               60000
maxbl4=               60000
  allocated                5588931  DP for d2temp 
sifcfg setup: record length 4096 DP
# d1 elements per record  3272
# d2 elements per record  2730
  The MR-CISD density will be calculated.
 item #                     1 suffix=:.drt1.state1:
 read_civout: repnuc=  -62.1137172157958     
================================================================================
  Reading record                      1  of civout
 INFO:ref#  1vector#  1 method:  0 last record  0max overlap with ref# 95% root-following 0
 MR-CISD energy:   -94.36141904   -32.24770183
 residuum:     0.00004500
 deltae:     0.00000000

          sovref  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1     0.94706787     0.00025589    -0.00126670    -0.00178264    -0.00337622     0.28275207     0.09678757    -0.06428875
 ref:   2    -0.00022678     0.94640844     0.02607982    -0.02063048     0.00040584     0.01457143     0.00185925     0.29125173
 ref:   3    -0.00089498     0.02479871    -0.93609036    -0.03112810     0.00196800    -0.00286088     0.00760621    -0.00568138

                ci   9         ci  10         ci  11         ci  12         ci  13         ci  14         ci  15         ci  16

                ci  17         ci  18         ci  19         ci  20         ci  21         ci  22         ci  23         ci  24

                ci  25         ci  26         ci  27         ci  28         ci  29         ci  30         ci  31         ci  32

                ci  33         ci  34         ci  35         ci  36         ci  37         ci  38         ci  39         ci  40

          tciref  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1     0.94706787     0.00025589    -0.00126670    -0.00178264    -0.00337622     0.00000000     0.00000000     0.00000000
 ref:   2    -0.00022678     0.94640844     0.02607982    -0.02063048     0.00040584     0.00000000     0.00000000     0.00000000
 ref:   3    -0.00089498     0.02479871    -0.93609036    -0.03112810     0.00196800     0.00000000     0.00000000     0.00000000

                ci   9         ci  10         ci  11         ci  12         ci  13         ci  14         ci  15         ci  16

                ci  17         ci  18         ci  19         ci  20         ci  21         ci  22         ci  23         ci  24

                ci  25         ci  26         ci  27         ci  28         ci  29         ci  30         ci  31         ci  32

                ci  33         ci  34         ci  35         ci  36         ci  37         ci  38         ci  39         ci  40
 computing final density
 =========== Executing IN-CORE method ==========
--------------------------------------------------------------------------------
  1e-density for root #    1
--------------------------------------------------------------------------------
================================================================================
   DYZ=     242  DYX=     503  DYW=     590
   D0Z=     147  D0Y=     928  D0X=     875  D0W=    1184
  DDZI=     189 DDYI=     692 DDXI=     675 DDWI=     866
  DDZE=       0 DDYE=     104 DDXE=     105 DDWE=     140
================================================================================
--------------------------------------------------------------------------------
  2e-density  for root #    1
--------------------------------------------------------------------------------
 call to prpd23 iwx=1
 call to prpd23 iwx=2
 starting prpd24 .... 
d2temp:fourex offset=     0 N= 59550 init,final=   T   Fmode=   1
d2temp:fourex offset= 59550 N= 59292 init,final=   F   Fmode=   1
d2temp:fourex offset=118842 N= 59913 init,final=   F   Fmode=   1
d2temp:fourex offset=178755 N= 59058 init,final=   F   Fmode=   1
d2temp:fourex offset=237813 N= 59667 init,final=   F   Fmode=   1
d2temp:fourex offset=297480 N= 59730 init,final=   F   Fmode=   1
d2temp:fourex offset=357210 N= 59871 init,final=   F   Fmode=   1
d2temp:fourex offset=417081 N= 59544 init,final=   F   Fmode=   1
d2temp:fourex offset=476625 N= 59835 init,final=   F   Fmode=   1
d2temp:fourex offset=536460 N= 59460 init,final=   F   Fmode=   1
d2temp:fourex offset=595920 N= 59343 init,final=   F   Fmode=   1
d2temp:fourex offset=655263 N= 58947 init,final=   F   Fmode=   1
d2temp:fourex offset=714210 N= 58800 init,final=   F   Fmode=   1
d2temp:fourex offset=773010 N= 58800 init,final=   F   Fmode=   1
d2temp:fourex offset=831810 N= 58845 init,final=   F   Fmode=   1
d2temp:fourex offset=890655 N= 59778 init,final=   F   Fmode=   1
d2temp:fourex offset=950433 N= 57120 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58182 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59445 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59955 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58695 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58620 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58140 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56676 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59124 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56793 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59799 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 55743 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59250 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58128 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57915 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 55545 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58545 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 55800 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58632 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57285 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59235 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59400 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58692 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57828 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57120 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59007 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56460 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56628 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59040 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57795 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59907 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58548 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 54660 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59535 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57525 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58692 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59958 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56793 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 54714 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 60000 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56403 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56658 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59787 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59556 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56082 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56250 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59892 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58140 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 55935 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58860 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N=     0 init,final=   F   Tmode=   1
 starting prpd24 .... 
d2temp:fourex offset=     0 N= 59550 init,final=   T   Fmode=   1
d2temp:fourex offset= 59550 N= 59292 init,final=   F   Fmode=   1
d2temp:fourex offset=118842 N= 59913 init,final=   F   Fmode=   1
d2temp:fourex offset=178755 N= 59058 init,final=   F   Fmode=   1
d2temp:fourex offset=237813 N= 59667 init,final=   F   Fmode=   1
d2temp:fourex offset=297480 N= 59730 init,final=   F   Fmode=   1
d2temp:fourex offset=357210 N= 59871 init,final=   F   Fmode=   1
d2temp:fourex offset=417081 N= 59544 init,final=   F   Fmode=   1
d2temp:fourex offset=476625 N= 59835 init,final=   F   Fmode=   1
d2temp:fourex offset=536460 N= 59460 init,final=   F   Fmode=   1
d2temp:fourex offset=595920 N= 59343 init,final=   F   Fmode=   1
d2temp:fourex offset=655263 N= 58947 init,final=   F   Fmode=   1
d2temp:fourex offset=714210 N= 58800 init,final=   F   Fmode=   1
d2temp:fourex offset=773010 N= 58800 init,final=   F   Fmode=   1
d2temp:fourex offset=831810 N= 58845 init,final=   F   Fmode=   1
d2temp:fourex offset=890655 N= 59778 init,final=   F   Fmode=   1
d2temp:fourex offset=950433 N= 57120 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58182 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59445 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59955 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58695 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58620 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58140 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56676 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59124 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56793 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59799 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 55743 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59250 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58128 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57915 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 55545 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58545 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 55800 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58632 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57285 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59235 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59400 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58692 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57828 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57120 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59007 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56460 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56628 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59040 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57795 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59907 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58548 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 54660 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59535 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57525 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58692 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59958 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56793 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 54714 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 60000 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56403 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56658 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59787 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59556 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56082 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56250 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59892 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58140 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 55935 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58860 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N=     0 init,final=   F   Tmode=   1
================================================================================
   DYZ=     242  DYX=    7351  DYW=    8696
   D0Z=     147  D0Y=    8560  D0X=    7763  D0W=   10896
  DDZI=     735 DDYI=    2660 DDXI=    2520 DDWI=    3140
  DDZE=       0 DDYE=       0 DDXE=       0 DDWE=       0
================================================================================
  entering dump_fourext .... 
d2temp:fourex offset=     1 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset= 60001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=120001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=180001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=240001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=300001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=360001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=420001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=480001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=540001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=600001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=660001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=720001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=780001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=840001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=900001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=960001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
Trace of MO density:    12.000000
   12  correlated and     4  frozen core electrons

           D1 (irrep)  block   1

                MO   1         MO   2         MO   3         MO   4         MO   5         MO   6         MO   7         MO   8
   MO   1     2.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
   MO   2     0.00000000     2.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
   MO   3     0.00000000     0.00000000     1.97131025     0.00026457     0.00084040     0.00052108     0.00033513     0.00185792
   MO   4     0.00000000     0.00000000     0.00026457     1.97023055    -0.00018771    -0.00678785    -0.00295539     0.00069909
   MO   5     0.00000000     0.00000000     0.00084040    -0.00018771     1.97280835    -0.00044915    -0.00031298     0.00018916
   MO   6     0.00000000     0.00000000     0.00052108    -0.00678785    -0.00044915     1.97224467     0.00103260     0.00032271
   MO   7     0.00000000     0.00000000     0.00033513    -0.00295539    -0.00031298     0.00103260     1.95616082    -0.00892835
   MO   8     0.00000000     0.00000000     0.00185792     0.00069909     0.00018916     0.00032271    -0.00892835     0.98969953
   MO   9     0.00000000     0.00000000     0.00025945    -0.00179311     0.00021232    -0.00362991    -0.00306874    -0.01846815
   MO  10     0.00000000     0.00000000     0.00163180    -0.00522755    -0.00095553    -0.00737688    -0.00135370     0.00252766
   MO  11     0.00000000     0.00000000    -0.01087764    -0.00398611     0.00583174     0.01558757    -0.00141654    -0.00406095
   MO  12     0.00000000     0.00000000    -0.01021136     0.01347856    -0.00160056     0.00217651     0.01533856     0.00284465
   MO  13     0.00000000     0.00000000    -0.00113543     0.00475340     0.01193671     0.00349134     0.02645155     0.00775877
   MO  14     0.00000000     0.00000000    -0.00105078     0.01140070     0.00551773     0.01185740    -0.02667644     0.00324318
   MO  15     0.00000000     0.00000000    -0.00114359     0.00819739     0.00259355     0.01005624     0.00327899    -0.03276492
   MO  16     0.00000000     0.00000000     0.00011057     0.00789135     0.00022410     0.00924134     0.01468841    -0.02046069
   MO  17     0.00000000     0.00000000    -0.00946033     0.00166270     0.01017051     0.01309114     0.00657490     0.00223805
   MO  18     0.00000000     0.00000000     0.00520253     0.00143423    -0.00629647    -0.01050976     0.00123135     0.00801135
   MO  19     0.00000000     0.00000000     0.01183524    -0.00324536     0.00291126    -0.00417117    -0.00133602     0.00588413
   MO  20     0.00000000     0.00000000    -0.00952743     0.01205643     0.00183935     0.00267231     0.01933413     0.00525158
   MO  21     0.00000000     0.00000000     0.00119784    -0.01959008    -0.01167678    -0.01793089    -0.02357197     0.00443486
   MO  22     0.00000000     0.00000000    -0.00058629    -0.01069362     0.00402804    -0.02356123     0.00955323    -0.00213196
   MO  23     0.00000000     0.00000000     0.00587215     0.00796641    -0.01162999     0.01605783     0.01914719    -0.00715217
   MO  24     0.00000000     0.00000000    -0.00044050     0.01114643     0.00066299     0.00690299    -0.05276551    -0.00284883
   MO  25     0.00000000     0.00000000     0.00136709    -0.00018331     0.01000111     0.00088498     0.00551973     0.00245683
   MO  26     0.00000000     0.00000000    -0.00064425     0.00764307     0.00583315     0.01158642    -0.00359539     0.00762929
   MO  27     0.00000000     0.00000000    -0.00196176     0.00767425    -0.00499395     0.00064082    -0.00424527     0.00211506
   MO  28     0.00000000     0.00000000    -0.00388728    -0.00152916     0.00259312    -0.00089443    -0.00029653     0.00692062
   MO  29     0.00000000     0.00000000    -0.00001650     0.00833668     0.00391486     0.00473134     0.00256787    -0.00122348
   MO  30     0.00000000     0.00000000    -0.00363031    -0.00222882     0.00015599    -0.00202322    -0.00326695     0.00799170
   MO  31     0.00000000     0.00000000     0.00405392     0.00056170    -0.00137007     0.00007466     0.00227885    -0.00937169
   MO  32     0.00000000     0.00000000    -0.00453788    -0.00534678    -0.00755001    -0.00824449    -0.00416804     0.00951437
   MO  33     0.00000000     0.00000000     0.00173540    -0.00107420    -0.01052183     0.00140018    -0.00184069    -0.00160595
   MO  34     0.00000000     0.00000000    -0.00223421    -0.00300248     0.00101553    -0.00336958    -0.00382214     0.02080621
   MO  35     0.00000000     0.00000000     0.00047431     0.00302665     0.00046227     0.00220872     0.00721618    -0.00262871
   MO  36     0.00000000     0.00000000     0.00690568    -0.00492005    -0.00573663    -0.00368882    -0.00400646     0.01116463
   MO  37     0.00000000     0.00000000     0.00767683    -0.00294108     0.00520740    -0.00082468     0.00151348     0.01358280
   MO  38     0.00000000     0.00000000     0.00864251     0.00576237    -0.00262235    -0.00409621    -0.00613703    -0.00253969
   MO  39     0.00000000     0.00000000     0.00580198    -0.00181323    -0.00241357     0.00219053    -0.00602762     0.00071694
   MO  40     0.00000000     0.00000000    -0.00352450     0.00159925     0.00108103    -0.00290425    -0.00793361     0.00772675
   MO  41     0.00000000     0.00000000    -0.00021155     0.00324678    -0.00218189     0.00471924    -0.02523488     0.00147296
   MO  42     0.00000000     0.00000000    -0.00013272     0.00161141    -0.00048702     0.00159305    -0.00241259     0.00426660
   MO  43     0.00000000     0.00000000    -0.00017376     0.00390021     0.00777343     0.00265757    -0.00785414     0.00300187
   MO  44     0.00000000     0.00000000    -0.00101670     0.00424833    -0.00073852    -0.00376129    -0.00260337     0.00501781
   MO  45     0.00000000     0.00000000     0.00232117     0.00843678    -0.00401754     0.00695390     0.00878074     0.01482260
   MO  46     0.00000000     0.00000000     0.00064055    -0.00082068     0.00016832    -0.00020393     0.00190640     0.00409839
   MO  47     0.00000000     0.00000000     0.00031244    -0.00199360    -0.00122526    -0.00116071     0.00051125     0.01713801
   MO  48     0.00000000     0.00000000     0.00273685     0.00283991    -0.00231670     0.00157565     0.00292477     0.00726347
   MO  49     0.00000000     0.00000000    -0.00133355    -0.00349598     0.00105568    -0.00128546    -0.00074938    -0.00213356
   MO  50     0.00000000     0.00000000     0.00163771    -0.00165576    -0.00066729     0.00187119     0.00445344     0.02652017
   MO  51     0.00000000     0.00000000    -0.00109406    -0.00137697     0.01212055    -0.00441896    -0.00686149     0.00694130
   MO  52     0.00000000     0.00000000     0.00541067     0.00390706     0.00051692     0.00574270    -0.00778809    -0.00529396
   MO  53     0.00000000     0.00000000     0.00075449    -0.00802948    -0.00219190    -0.01328355    -0.01045869    -0.00055408
   MO  54     0.00000000     0.00000000    -0.00259478     0.00553502    -0.00236117     0.00619440    -0.01634538     0.00448992
   MO  55     0.00000000     0.00000000    -0.00036746     0.00012516    -0.00024554     0.00033201    -0.00057805    -0.00050751
   MO  56     0.00000000     0.00000000    -0.00000550     0.00030826     0.00058549    -0.00743258    -0.00750125     0.00068068
   MO  57     0.00000000     0.00000000    -0.00046191     0.00229284     0.00128762     0.00420614     0.00326392    -0.00067771
   MO  58     0.00000000     0.00000000    -0.00009121     0.00125157     0.00007329    -0.00077822    -0.00176106     0.00006199
   MO  59     0.00000000     0.00000000    -0.00018329    -0.00076835     0.00018466    -0.00007397     0.00016602     0.00153539
   MO  60     0.00000000     0.00000000     0.00430751     0.00013739    -0.00029991     0.00018985    -0.00009292    -0.00826585
   MO  61     0.00000000     0.00000000     0.00294934    -0.00091126    -0.00052871    -0.00019817     0.00034792     0.00163963
   MO  62     0.00000000     0.00000000     0.00038746    -0.00072647     0.00360309     0.00064959    -0.00002554    -0.00072825
   MO  63     0.00000000     0.00000000    -0.00005428     0.00233495     0.00124194     0.00172079     0.00183376    -0.00107919
   MO  64     0.00000000     0.00000000    -0.00005773    -0.00140853    -0.00022940     0.00080227     0.00402295     0.00045236
   MO  65     0.00000000     0.00000000     0.00067212    -0.00006572    -0.00032781     0.00001007    -0.00036941     0.00098652
   MO  66     0.00000000     0.00000000    -0.00066887     0.00004712     0.00021195    -0.00004889     0.00032150    -0.00112930
   MO  67     0.00000000     0.00000000    -0.00098172     0.00000938     0.00003927    -0.00005056     0.00013030    -0.00113879
   MO  68     0.00000000     0.00000000    -0.00401741     0.00014991     0.00066979     0.00039779     0.00006638    -0.00452375
   MO  69     0.00000000     0.00000000    -0.00001515     0.00033233    -0.00013210    -0.00250179     0.01103641    -0.00017130
   MO  70     0.00000000     0.00000000    -0.00130536    -0.00511920    -0.00040439     0.00114580    -0.00072048    -0.00065492
   MO  71     0.00000000     0.00000000     0.00540679    -0.00084306    -0.00350184    -0.00019623    -0.00046976     0.00143248
   MO  72     0.00000000     0.00000000    -0.00082422    -0.00013484     0.00007168    -0.00018157    -0.00004678    -0.00023362
   MO  73     0.00000000     0.00000000    -0.00278960     0.00232008    -0.00433636    -0.00011023     0.00038499    -0.00135565
   MO  74     0.00000000     0.00000000    -0.00064280    -0.00519291    -0.00090042    -0.00069284    -0.00121955    -0.00005924
   MO  75     0.00000000     0.00000000     0.00014720     0.00668975    -0.00014701     0.00231066    -0.00003803     0.00006729
   MO  76     0.00000000     0.00000000     0.00019191     0.00016877    -0.00008685     0.00002197    -0.00016426    -0.00002576
   MO  77     0.00000000     0.00000000     0.00004101    -0.00001614     0.00012024    -0.00006342    -0.00067038    -0.00004312
   MO  78     0.00000000     0.00000000    -0.00020393    -0.00043526    -0.00012047    -0.00041980    -0.00358214    -0.00006277
   MO  79     0.00000000     0.00000000     0.00043885     0.00011372     0.00271515    -0.00012264    -0.00034280     0.00018945
   MO  80     0.00000000     0.00000000    -0.00009409    -0.00036091    -0.00007693    -0.00017181    -0.00284421     0.00018763
   MO  81     0.00000000     0.00000000     0.00024665     0.00006902     0.00045533     0.00011368    -0.00029439    -0.00213061
   MO  82     0.00000000     0.00000000     0.00003723     0.00122112    -0.00010028     0.00149375    -0.00184186     0.00014611

                MO   9         MO  10         MO  11         MO  12         MO  13         MO  14         MO  15         MO  16
   MO   3     0.00025945     0.00163180    -0.01087764    -0.01021136    -0.00113543    -0.00105078    -0.00114359     0.00011057
   MO   4    -0.00179311    -0.00522755    -0.00398611     0.01347856     0.00475340     0.01140070     0.00819739     0.00789135
   MO   5     0.00021232    -0.00095553     0.00583174    -0.00160056     0.01193671     0.00551773     0.00259355     0.00022410
   MO   6    -0.00362991    -0.00737688     0.01558757     0.00217651     0.00349134     0.01185740     0.01005624     0.00924134
   MO   7    -0.00306874    -0.00135370    -0.00141654     0.01533856     0.02645155    -0.02667644     0.00327899     0.01468841
   MO   8    -0.01846815     0.00252766    -0.00406095     0.00284465     0.00775877     0.00324318    -0.03276492    -0.02046069
   MO   9     0.98968117     0.00188993    -0.00057440    -0.00018496    -0.00232746    -0.00252140     0.00761379    -0.00277516
   MO  10     0.00188993     0.01019634    -0.00009510    -0.00028158    -0.00059456     0.00065769    -0.00228015     0.00293117
   MO  11    -0.00057440    -0.00009510     0.00241625    -0.00094855    -0.00041020     0.00056807     0.00032304    -0.00033025
   MO  12    -0.00018496    -0.00028158    -0.00094855     0.00214950     0.00101435    -0.00044093     0.00026746     0.00089152
   MO  13    -0.00232746    -0.00059456    -0.00041020     0.00101435     0.00161496    -0.00067170    -0.00000063     0.00054862
   MO  14    -0.00252140     0.00065769     0.00056807    -0.00044093    -0.00067170     0.00206750     0.00002036    -0.00037082
   MO  15     0.00761379    -0.00228015     0.00032304     0.00026746    -0.00000063     0.00002036     0.00248504     0.00051955
   MO  16    -0.00277516     0.00293117    -0.00033025     0.00089152     0.00054862    -0.00037082     0.00051955     0.00264398
   MO  17    -0.00676593     0.00007805     0.00266451    -0.00033683     0.00001779     0.00081264     0.00023780    -0.00036073
   MO  18    -0.01372877    -0.00057884    -0.00163613     0.00120406     0.00047963    -0.00039160    -0.00055213     0.00018287
   MO  19     0.01118102     0.00062883    -0.00054825    -0.00080037    -0.00007796     0.00018973    -0.00016336    -0.00041819
   MO  20     0.01749990     0.00033122    -0.00091662     0.00257508     0.00123236    -0.00074069     0.00067130     0.00108545
   MO  21     0.00991886     0.00065580     0.00093204    -0.00204942    -0.00205671     0.00030245    -0.00073965    -0.00194919
   MO  22    -0.00659839     0.00036310    -0.00100100     0.00010123     0.00043105    -0.00083624    -0.00057916    -0.00020781
   MO  23    -0.00226411    -0.00010476    -0.00010167    -0.00012807     0.00001824     0.00026168     0.00054454     0.00053808
   MO  24    -0.00081211     0.00020138     0.00183378    -0.00220156    -0.00276072     0.00440468     0.00048114    -0.00254060
   MO  25     0.00825622     0.00012642     0.00024596    -0.00054213     0.00046348     0.00017271     0.00011746    -0.00019628
   MO  26    -0.00041921    -0.00033442    -0.00045628     0.00100474     0.00080643     0.00138658     0.00046254     0.00032460
   MO  27    -0.00146564     0.00018397    -0.00129401     0.00162881     0.00008404     0.00049869     0.00018877     0.00042473
   MO  28    -0.01514341     0.00001126     0.00029004    -0.00031143     0.00013730     0.00017526    -0.00081601    -0.00026933
   MO  29     0.00038361    -0.00018828    -0.00045718     0.00021333     0.00039065     0.00036382     0.00024774     0.00013038
   MO  30    -0.00077312     0.00016047     0.00007320    -0.00011755    -0.00004175     0.00006179    -0.00058324    -0.00033141
   MO  31    -0.00678735    -0.00036201    -0.00012199     0.00035554     0.00004146    -0.00010551     0.00047164     0.00042346
   MO  32     0.00510303    -0.00039050    -0.00041543     0.00034049    -0.00042410    -0.00022725    -0.00047675    -0.00083054
   MO  33    -0.00061942     0.00032248     0.00029644    -0.00060175    -0.00083251    -0.00011815    -0.00021372    -0.00009218
   MO  34    -0.00306652     0.00061236     0.00017494    -0.00004923     0.00011659     0.00016793    -0.00141138    -0.00069691
   MO  35    -0.00093173     0.00135349    -0.00026753     0.00065719     0.00034763    -0.00014739     0.00002623     0.00096877
   MO  36    -0.00196098     0.00002746     0.00028781    -0.00058950    -0.00033005     0.00006241    -0.00061970    -0.00068068
   MO  37     0.00143468    -0.00021899     0.00017572    -0.00051753     0.00027718     0.00004099    -0.00054885    -0.00072078
   MO  38    -0.00091656     0.00007798    -0.00267712     0.00091965     0.00034598    -0.00023156    -0.00002739     0.00082239
   MO  39    -0.00052584     0.00006486     0.00080722    -0.00164006    -0.00103757     0.00106168    -0.00012383    -0.00099296
   MO  40     0.00106365    -0.00003453    -0.00024821     0.00090170     0.00006585     0.00010753    -0.00025735    -0.00016897
   MO  41     0.00137688    -0.00004401     0.00081607    -0.00083941    -0.00132768     0.00125901    -0.00005212    -0.00068220
   MO  42     0.00290023    -0.00005095     0.00007640     0.00009812    -0.00007639    -0.00012850    -0.00018049    -0.00004581
   MO  43    -0.00116881     0.00006949     0.00004149     0.00007918     0.00032286     0.00085910     0.00007355     0.00008376
   MO  44     0.00101715    -0.00006113     0.00042626    -0.00025346    -0.00075645     0.00182332     0.00012820    -0.00158339
   MO  45    -0.00166724    -0.00058921     0.00076508    -0.00036542     0.00014859    -0.00047796    -0.00080121    -0.00034964
   MO  46     0.00663977     0.00211360    -0.00006156    -0.00002340     0.00005544    -0.00004424    -0.00047850     0.00058604
   MO  47    -0.00627989     0.00156882    -0.00052633     0.00042427     0.00025021     0.00015900    -0.00147384     0.00007124
   MO  48     0.00102388    -0.00032327    -0.00013116    -0.00026584    -0.00009047    -0.00013532    -0.00042817    -0.00045598
   MO  49    -0.00297010    -0.00107729    -0.00026136     0.00046117     0.00033910    -0.00006810     0.00034281     0.00000189
   MO  50     0.00445833    -0.00088394    -0.00095021     0.00090948     0.00085941    -0.00025958    -0.00100530    -0.00066008
   MO  51     0.00041695     0.00006403     0.00089420    -0.00120700     0.00045251     0.00084016    -0.00046823    -0.00101799
   MO  52     0.00091282    -0.00009713    -0.00052488    -0.00042691    -0.00032310     0.00141880     0.00063281     0.00002427
   MO  53    -0.00006854     0.00022034     0.00052425    -0.00195406    -0.00155124    -0.00053499    -0.00098343    -0.00151020
   MO  54    -0.00162728    -0.00013764     0.00087116    -0.00016842    -0.00111373     0.00209615     0.00017839    -0.00084228
   MO  55     0.02612201     0.00012874     0.00002064    -0.00000980    -0.00029762    -0.00020701     0.00089531    -0.00029331
   MO  56    -0.00007456     0.00017642    -0.00120641     0.00020046    -0.00039545     0.00000673    -0.00041419    -0.00066752
   MO  57    -0.00007661    -0.00008195     0.00039711     0.00001090     0.00029574     0.00013199     0.00026406     0.00027483
   MO  58     0.00028952    -0.00008892    -0.00021477     0.00005852    -0.00013803     0.00007753    -0.00005027    -0.00015898
   MO  59     0.00001366    -0.00055288    -0.00005300    -0.00004244     0.00007425    -0.00002937     0.00004334    -0.00017839
   MO  60     0.00019811    -0.00009302     0.00026017     0.00013263    -0.00003934    -0.00006654     0.00051021     0.00027436
   MO  61    -0.00030998     0.00007255    -0.00056110    -0.00043223     0.00004396    -0.00012735    -0.00009766    -0.00000287
   MO  62     0.00019289     0.00003853     0.00000695     0.00015964    -0.00003267     0.00003789     0.00008650     0.00014514
   MO  63     0.00025560    -0.00029213     0.00040323    -0.00033768    -0.00000035    -0.00012920     0.00015201    -0.00001152
   MO  64     0.00000998     0.00006017    -0.00023242     0.00047426     0.00063893    -0.00083469    -0.00006175     0.00027657
   MO  65     0.00500213    -0.00020993    -0.00006491    -0.00001914    -0.00002092    -0.00001061     0.00000120    -0.00012984
   MO  66    -0.00291056     0.00053187     0.00007052    -0.00001579    -0.00002049     0.00004739    -0.00006897     0.00022529
   MO  67     0.00317805     0.00092377     0.00002895    -0.00004487    -0.00008408     0.00005207    -0.00010274     0.00027658
   MO  68     0.00040767    -0.00037383    -0.00010601     0.00026994     0.00010134    -0.00009126     0.00035110     0.00016158
   MO  69     0.00007434    -0.00005844     0.00061110    -0.00034895    -0.00001697     0.00021311     0.00013651    -0.00038716
   MO  70     0.00007919    -0.00003728     0.00026266    -0.00001188     0.00007124     0.00003757     0.00016198     0.00016244
   MO  71    -0.00010077     0.00003267    -0.00031098    -0.00017734    -0.00007868    -0.00002312    -0.00017114    -0.00007731
   MO  72    -0.00020654     0.00007809    -0.00002613     0.00005888    -0.00000286    -0.00000243    -0.00000194     0.00001628
   MO  73     0.00001298    -0.00001175     0.00020288     0.00000160    -0.00009661    -0.00005208     0.00005661    -0.00000346
   MO  74     0.00007397     0.00008912    -0.00045173     0.00022272    -0.00004856    -0.00004371    -0.00007585    -0.00009023
   MO  75    -0.00007818    -0.00007574     0.00039249    -0.00019715    -0.00000262    -0.00013134     0.00000281     0.00019374
   MO  76     0.00209412    -0.00000168    -0.00000531    -0.00001826    -0.00001555    -0.00001712     0.00001970    -0.00001299
   MO  77     0.00127234    -0.00000667    -0.00000359     0.00000283    -0.00000254    -0.00000887     0.00003722     0.00001931
   MO  78    -0.00027848    -0.00001873    -0.00003716     0.00006535     0.00000886    -0.00000424     0.00000362     0.00008506
   MO  79     0.00003664     0.00001303     0.00007303    -0.00008479     0.00020350     0.00013462     0.00002469    -0.00002997
   MO  80    -0.00000043     0.00003993     0.00002224    -0.00019977    -0.00026079     0.00011306    -0.00009709    -0.00020255
   MO  81    -0.00003612    -0.00000062    -0.00027401    -0.00026277    -0.00003443     0.00004400    -0.00001036     0.00000659
   MO  82    -0.00000649    -0.00005051    -0.00002440     0.00014798    -0.00005796     0.00034336     0.00012192     0.00005945

                MO  17         MO  18         MO  19         MO  20         MO  21         MO  22         MO  23         MO  24
   MO   3    -0.00946033     0.00520253     0.01183524    -0.00952743     0.00119784    -0.00058629     0.00587215    -0.00044050
   MO   4     0.00166270     0.00143423    -0.00324536     0.01205643    -0.01959008    -0.01069362     0.00796641     0.01114643
   MO   5     0.01017051    -0.00629647     0.00291126     0.00183935    -0.01167678     0.00402804    -0.01162999     0.00066299
   MO   6     0.01309114    -0.01050976    -0.00417117     0.00267231    -0.01793089    -0.02356123     0.01605783     0.00690299
   MO   7     0.00657490     0.00123135    -0.00133602     0.01933413    -0.02357197     0.00955323     0.01914719    -0.05276551
   MO   8     0.00223805     0.00801135     0.00588413     0.00525158     0.00443486    -0.00213196    -0.00715217    -0.00284883
   MO   9    -0.00676593    -0.01372877     0.01118102     0.01749990     0.00991886    -0.00659839    -0.00226411    -0.00081211
   MO  10     0.00007805    -0.00057884     0.00062883     0.00033122     0.00065580     0.00036310    -0.00010476     0.00020138
   MO  11     0.00266451    -0.00163613    -0.00054825    -0.00091662     0.00093204    -0.00100100    -0.00010167     0.00183378
   MO  12    -0.00033683     0.00120406    -0.00080037     0.00257508    -0.00204942     0.00010123    -0.00012807    -0.00220156
   MO  13     0.00001779     0.00047963    -0.00007796     0.00123236    -0.00205671     0.00043105     0.00001824    -0.00276072
   MO  14     0.00081264    -0.00039160     0.00018973    -0.00074069     0.00030245    -0.00083624     0.00026168     0.00440468
   MO  15     0.00023780    -0.00055213    -0.00016336     0.00067130    -0.00073965    -0.00057916     0.00054454     0.00048114
   MO  16    -0.00036073     0.00018287    -0.00041819     0.00108545    -0.00194919    -0.00020781     0.00053808    -0.00254060
   MO  17     0.00503796    -0.00080698    -0.00101570     0.00020378     0.00053165    -0.00092235    -0.00105600     0.00277884
   MO  18    -0.00080698     0.00237702    -0.00054674     0.00071340    -0.00119548     0.00079817    -0.00038333    -0.00144946
   MO  19    -0.00101570    -0.00054674     0.00149882    -0.00030375     0.00075943     0.00015749     0.00050275     0.00084867
   MO  20     0.00020378     0.00071340    -0.00030375     0.00469561    -0.00189284    -0.00046361    -0.00089574    -0.00291070
   MO  21     0.00053165    -0.00119548     0.00075943    -0.00189284     0.00441628     0.00028801    -0.00072339     0.00376022
   MO  22    -0.00092235     0.00079817     0.00015749    -0.00046361     0.00028801     0.00242443    -0.00078599    -0.00105046
   MO  23    -0.00105600    -0.00038333     0.00050275    -0.00089574    -0.00072339    -0.00078599     0.00421183     0.00106398
   MO  24     0.00277884    -0.00144946     0.00084867    -0.00291070     0.00376022    -0.00105046     0.00106398     0.01336934
   MO  25    -0.00016364    -0.00095321     0.00080230    -0.00027761     0.00012627     0.00017161    -0.00011537     0.00025437
   MO  26    -0.00020331     0.00070914     0.00031986     0.00111331    -0.00174217    -0.00049903     0.00086099     0.00154101
   MO  27    -0.00046920     0.00161965    -0.00039228     0.00192062    -0.00114942    -0.00019654     0.00040518     0.00067753
   MO  28     0.00060167     0.00049026    -0.00056284    -0.00142176    -0.00025031     0.00052243    -0.00012751     0.00023047
   MO  29    -0.00059198    -0.00007455     0.00034573     0.00002100    -0.00078613     0.00027454     0.00042467     0.00041984
   MO  30    -0.00000549     0.00003111    -0.00001849    -0.00022048     0.00030504     0.00004094    -0.00019069     0.00002841
   MO  31     0.00043888     0.00074992    -0.00057553     0.00006235    -0.00058592     0.00008355     0.00011551    -0.00010533
   MO  32    -0.00047906     0.00044061    -0.00000358     0.00069375     0.00105578     0.00000421    -0.00006449     0.00032033
   MO  33    -0.00046057    -0.00035287    -0.00005573    -0.00119182     0.00126459    -0.00016188     0.00083441     0.00052223
   MO  34    -0.00014103    -0.00009562    -0.00010373    -0.00045139     0.00036698     0.00039246    -0.00048425    -0.00005410
   MO  35     0.00049660     0.00050246    -0.00018165     0.00115949    -0.00082987    -0.00018228    -0.00001165    -0.00082632
   MO  36     0.00097807     0.00027049     0.00015807    -0.00048529     0.00111582    -0.00004208    -0.00075559     0.00077413
   MO  37     0.00071903    -0.00008140     0.00059892    -0.00016696     0.00036556    -0.00002174    -0.00029966     0.00031694
   MO  38    -0.00388098     0.00145673     0.00057879     0.00042347    -0.00169493     0.00068944     0.00073396    -0.00205197
   MO  39     0.00027202    -0.00106139     0.00080326    -0.00230082     0.00175889    -0.00048946     0.00162111     0.00392886
   MO  40     0.00079171     0.00077577    -0.00051858     0.00159046    -0.00018829    -0.00002328    -0.00185167    -0.00022197
   MO  41     0.00068871    -0.00059247    -0.00015529    -0.00087940     0.00121224    -0.00147959    -0.00051868     0.00296066
   MO  42     0.00016278     0.00000756    -0.00025438     0.00035712    -0.00002990    -0.00059459    -0.00051508    -0.00064211
   MO  43     0.00023904     0.00007125     0.00005812    -0.00001231    -0.00056960     0.00004461    -0.00109813     0.00085725
   MO  44     0.00197787    -0.00008882     0.00055014    -0.00020519     0.00174975     0.00025329     0.00114519     0.00708351
   MO  45     0.00049430    -0.00071898    -0.00024853    -0.00064986     0.00012489    -0.00082214     0.00045002    -0.00129935
   MO  46    -0.00028528    -0.00037120     0.00038176     0.00036945     0.00002569    -0.00001677     0.00014592    -0.00055805
   MO  47    -0.00036660     0.00076800     0.00004293     0.00031408    -0.00028983     0.00001892     0.00028328    -0.00044044
   MO  48    -0.00048570    -0.00020786     0.00017900    -0.00064817     0.00029143     0.00016663     0.00025962    -0.00008699
   MO  49    -0.00013684     0.00066451    -0.00023524     0.00042927    -0.00079090    -0.00002829     0.00022561    -0.00046339
   MO  50    -0.00125191     0.00079030     0.00036631     0.00122504    -0.00090100    -0.00020249     0.00051987    -0.00174857
   MO  51     0.00119339    -0.00107233     0.00069986    -0.00126718     0.00059754     0.00093022    -0.00214969     0.00173126
   MO  52    -0.00153436    -0.00027375     0.00096031    -0.00109013    -0.00032113    -0.00038804     0.00214408     0.00312348
   MO  53    -0.00038384    -0.00079899     0.00030322    -0.00258668     0.00338000     0.00082162    -0.00079773     0.00110272
   MO  54     0.00167053    -0.00012807    -0.00024533    -0.00025656     0.00091401    -0.00141343     0.00046109     0.00554905
   MO  55    -0.00063338    -0.00158752     0.00125992     0.00210425     0.00123910    -0.00089791    -0.00032092     0.00009262
   MO  56    -0.00092540     0.00080437     0.00029621     0.00018133     0.00084518     0.00146419    -0.00117737     0.00099277
   MO  57     0.00048000    -0.00056965    -0.00000220     0.00019709    -0.00096185    -0.00082975     0.00097570    -0.00012491
   MO  58    -0.00008120     0.00020902    -0.00017797    -0.00003572     0.00003485    -0.00013119    -0.00011125     0.00038668
   MO  59    -0.00016498     0.00002829     0.00011925    -0.00007930    -0.00004857     0.00005840     0.00010584    -0.00011540
   MO  60     0.00104328     0.00016764    -0.00026336     0.00052855    -0.00012201    -0.00019175    -0.00062203     0.00006523
   MO  61    -0.00018155     0.00048370     0.00067725    -0.00020122     0.00000411    -0.00013979     0.00022463    -0.00017344
   MO  62     0.00003125     0.00022166    -0.00012583     0.00021134     0.00001640    -0.00017543    -0.00004967    -0.00004662
   MO  63     0.00040681    -0.00039886     0.00006676    -0.00045258     0.00024596    -0.00006518     0.00016830    -0.00005656
   MO  64    -0.00007442     0.00027412    -0.00003160     0.00068607    -0.00045517     0.00069315     0.00042396    -0.00179017
   MO  65    -0.00011379    -0.00003366     0.00014804     0.00016691     0.00015385    -0.00015905    -0.00000877     0.00001676
   MO  66     0.00013107     0.00002388    -0.00010366    -0.00014648    -0.00004598     0.00009865     0.00003146     0.00003815
   MO  67    -0.00002118    -0.00018583     0.00012907     0.00012139     0.00013471    -0.00004719     0.00004045     0.00005968
   MO  68    -0.00012710     0.00007053    -0.00017413     0.00028658    -0.00034588     0.00007500     0.00000252    -0.00026927
   MO  69     0.00104903    -0.00046359     0.00022396    -0.00038402     0.00058944     0.00005800     0.00095617     0.00174774
   MO  70     0.00013304     0.00002099    -0.00001902     0.00012945    -0.00022561    -0.00038912     0.00031806    -0.00032283
   MO  71    -0.00078625    -0.00000185     0.00028410    -0.00040946     0.00005796     0.00006176     0.00060866    -0.00004784
   MO  72     0.00006052     0.00003034    -0.00000600     0.00010150     0.00003478     0.00012616    -0.00014965     0.00003842
   MO  73     0.00041737    -0.00010600    -0.00015905     0.00010722     0.00009515    -0.00017760    -0.00008732     0.00009000
   MO  74    -0.00030607     0.00039027     0.00002976     0.00027618     0.00007417     0.00042366    -0.00031711     0.00007847
   MO  75     0.00012637    -0.00041715    -0.00017719    -0.00028243    -0.00014794    -0.00047136    -0.00010971    -0.00082481
   MO  76    -0.00003108    -0.00002171    -0.00001739    -0.00008583    -0.00006281     0.00006048     0.00003295    -0.00001766
   MO  77    -0.00002030     0.00002342    -0.00004054    -0.00004200    -0.00005493     0.00003098     0.00001867    -0.00005354
   MO  78    -0.00012947     0.00012003    -0.00003191     0.00014187    -0.00011750    -0.00016863     0.00003185    -0.00025099
   MO  79     0.00016347    -0.00012515     0.00012384    -0.00005284    -0.00016187     0.00020073    -0.00039453     0.00010162
   MO  80    -0.00008935     0.00002351     0.00003004    -0.00025823     0.00043642    -0.00000239    -0.00027063     0.00040691
   MO  81    -0.00068941    -0.00004442     0.00035072    -0.00051273    -0.00001988     0.00009122     0.00050462     0.00007956
   MO  82     0.00004478     0.00002662    -0.00006635     0.00014959    -0.00023111    -0.00027960     0.00005873     0.00064850

                MO  25         MO  26         MO  27         MO  28         MO  29         MO  30         MO  31         MO  32
   MO   3     0.00136709    -0.00064425    -0.00196176    -0.00388728    -0.00001650    -0.00363031     0.00405392    -0.00453788
   MO   4    -0.00018331     0.00764307     0.00767425    -0.00152916     0.00833668    -0.00222882     0.00056170    -0.00534678
   MO   5     0.01000111     0.00583315    -0.00499395     0.00259312     0.00391486     0.00015599    -0.00137007    -0.00755001
   MO   6     0.00088498     0.01158642     0.00064082    -0.00089443     0.00473134    -0.00202322     0.00007466    -0.00824449
   MO   7     0.00551973    -0.00359539    -0.00424527    -0.00029653     0.00256787    -0.00326695     0.00227885    -0.00416804
   MO   8     0.00245683     0.00762929     0.00211506     0.00692062    -0.00122348     0.00799170    -0.00937169     0.00951437
   MO   9     0.00825622    -0.00041921    -0.00146564    -0.01514341     0.00038361    -0.00077312    -0.00678735     0.00510303
   MO  10     0.00012642    -0.00033442     0.00018397     0.00001126    -0.00018828     0.00016047    -0.00036201    -0.00039050
   MO  11     0.00024596    -0.00045628    -0.00129401     0.00029004    -0.00045718     0.00007320    -0.00012199    -0.00041543
   MO  12    -0.00054213     0.00100474     0.00162881    -0.00031143     0.00021333    -0.00011755     0.00035554     0.00034049
   MO  13     0.00046348     0.00080643     0.00008404     0.00013730     0.00039065    -0.00004175     0.00004146    -0.00042410
   MO  14     0.00017271     0.00138658     0.00049869     0.00017526     0.00036382     0.00006179    -0.00010551    -0.00022725
   MO  15     0.00011746     0.00046254     0.00018877    -0.00081601     0.00024774    -0.00058324     0.00047164    -0.00047675
   MO  16    -0.00019628     0.00032460     0.00042473    -0.00026933     0.00013038    -0.00033141     0.00042346    -0.00083054
   MO  17    -0.00016364    -0.00020331    -0.00046920     0.00060167    -0.00059198    -0.00000549     0.00043888    -0.00047906
   MO  18    -0.00095321     0.00070914     0.00161965     0.00049026    -0.00007455     0.00003111     0.00074992     0.00044061
   MO  19     0.00080230     0.00031986    -0.00039228    -0.00056284     0.00034573    -0.00001849    -0.00057553    -0.00000358
   MO  20    -0.00027761     0.00111331     0.00192062    -0.00142176     0.00002100    -0.00022048     0.00006235     0.00069375
   MO  21     0.00012627    -0.00174217    -0.00114942    -0.00025031    -0.00078613     0.00030504    -0.00058592     0.00105578
   MO  22     0.00017161    -0.00049903    -0.00019654     0.00052243     0.00027454     0.00004094     0.00008355     0.00000421
   MO  23    -0.00011537     0.00086099     0.00040518    -0.00012751     0.00042467    -0.00019069     0.00011551    -0.00006449
   MO  24     0.00025437     0.00154101     0.00067753     0.00023047     0.00041984     0.00002841    -0.00010533     0.00032033
   MO  25     0.00222145     0.00029194    -0.00081158    -0.00031622     0.00034245    -0.00002199    -0.00051604    -0.00030759
   MO  26     0.00029194     0.00400010     0.00144552    -0.00013771     0.00047643    -0.00008658     0.00013742    -0.00003868
   MO  27    -0.00081158     0.00144552     0.00334483    -0.00033275     0.00041063    -0.00020100     0.00054221     0.00061566
   MO  28    -0.00031622    -0.00013771    -0.00033275     0.00221296    -0.00011442     0.00014619     0.00037387    -0.00042404
   MO  29     0.00034245     0.00047643     0.00041063    -0.00011442     0.00168040    -0.00010142    -0.00008441    -0.00035355
   MO  30    -0.00002199    -0.00008658    -0.00020100     0.00014619    -0.00010142     0.00115547    -0.00000208     0.00014283
   MO  31    -0.00051604     0.00013742     0.00054221     0.00037387    -0.00008441    -0.00000208     0.00139972    -0.00018470
   MO  32    -0.00030759    -0.00003868     0.00061566    -0.00042404    -0.00035355     0.00014283    -0.00018470     0.00169208
   MO  33    -0.00018764    -0.00071431    -0.00048286     0.00004753    -0.00009482     0.00003474    -0.00011924     0.00016601
   MO  34     0.00012923    -0.00007772    -0.00029930     0.00060629     0.00002653     0.00037789    -0.00038649     0.00023379
   MO  35    -0.00029567     0.00025869     0.00060027    -0.00015715    -0.00004521    -0.00016072     0.00041505     0.00006899
   MO  36    -0.00000120    -0.00015978    -0.00027691     0.00023520    -0.00037351     0.00008937     0.00004116     0.00012487
   MO  37     0.00021188     0.00008217    -0.00030901     0.00011357    -0.00004728     0.00010320    -0.00019854    -0.00006901
   MO  38    -0.00018408     0.00072105     0.00115903    -0.00036462     0.00071341    -0.00008957     0.00006704     0.00001457
   MO  39     0.00031254    -0.00006776    -0.00070569     0.00014152     0.00021701     0.00004162    -0.00027925     0.00001353
   MO  40    -0.00092610     0.00021073     0.00092356    -0.00013702    -0.00007630     0.00005361     0.00009944     0.00033925
   MO  41    -0.00058369    -0.00035878    -0.00010276     0.00000413    -0.00043010     0.00019082    -0.00012128     0.00017728
   MO  42     0.00037167    -0.00064493     0.00027027     0.00013335    -0.00031467     0.00000858     0.00009230     0.00017398
   MO  43     0.00098201     0.00136124    -0.00015231     0.00024546     0.00033290    -0.00012118    -0.00012002    -0.00054023
   MO  44     0.00020159     0.00111956     0.00128711     0.00007357     0.00031626    -0.00006839     0.00020001     0.00075260
   MO  45     0.00009270    -0.00167825    -0.00210187     0.00029507    -0.00052780     0.00035006    -0.00063990    -0.00020679
   MO  46     0.00002373     0.00017353    -0.00038030     0.00056070    -0.00045256    -0.00013448    -0.00009260    -0.00020238
   MO  47    -0.00012917     0.00049590     0.00089319    -0.00035053     0.00002666     0.00017058    -0.00056800     0.00075049
   MO  48    -0.00001406    -0.00071814    -0.00000125     0.00049205     0.00065845    -0.00080926    -0.00060610    -0.00016323
   MO  49    -0.00032131     0.00101166     0.00039121     0.00014922    -0.00082718    -0.00063259    -0.00003315     0.00041268
   MO  50    -0.00003938     0.00144736     0.00099112     0.00023012     0.00023342     0.00057876    -0.00003804     0.00068867
   MO  51     0.00198935     0.00037368    -0.00191848     0.00069296     0.00038471     0.00024976    -0.00071747    -0.00093914
   MO  52     0.00061125     0.00200152     0.00052568    -0.00047311     0.00079647    -0.00023436    -0.00001333    -0.00037865
   MO  53    -0.00022057    -0.00290271    -0.00186189     0.00043747    -0.00119314     0.00046567    -0.00026601     0.00073574
   MO  54    -0.00078837     0.00143809     0.00128021     0.00019479    -0.00003001     0.00010795     0.00007826     0.00044032
   MO  55     0.00110740    -0.00005778    -0.00016671    -0.00232788    -0.00002092    -0.00025331    -0.00111763     0.00080957
   MO  56    -0.00063951    -0.00064177     0.00092368    -0.00013079     0.00062616     0.00002951    -0.00006016     0.00047707
   MO  57    -0.00041103    -0.00040646     0.00008002    -0.00003319     0.00035501     0.00000774    -0.00002174    -0.00044324
   MO  58     0.00101624    -0.00068190     0.00098333     0.00024663    -0.00004808    -0.00006762     0.00027666     0.00018798
   MO  59    -0.00008238     0.00019316     0.00017995    -0.00001783    -0.00050950     0.00034407    -0.00047167    -0.00046478
   MO  60    -0.00032523    -0.00022941     0.00021406    -0.00059554    -0.00043153    -0.00046323     0.00075403    -0.00051116
   MO  61     0.00001257     0.00007156     0.00000832     0.00004032    -0.00011334    -0.00011493     0.00024250    -0.00011235
   MO  62    -0.00001270     0.00048256     0.00016348    -0.00004118    -0.00047470    -0.00011587     0.00025051     0.00027660
   MO  63     0.00016191    -0.00039423    -0.00041999    -0.00007743     0.00047989    -0.00021323    -0.00018993    -0.00049606
   MO  64     0.00006735    -0.00001536     0.00012397    -0.00002704    -0.00011433    -0.00000530     0.00012959     0.00022988
   MO  65     0.00020471    -0.00002838    -0.00003750    -0.00076171     0.00010817     0.00046627    -0.00012407     0.00030605
   MO  66    -0.00009937     0.00000106    -0.00000498     0.00035369     0.00011747     0.00057750     0.00057768    -0.00032635
   MO  67     0.00015555    -0.00008400     0.00001122    -0.00038570    -0.00008843    -0.00026551    -0.00038075     0.00009620
   MO  68    -0.00000663     0.00000114     0.00005620    -0.00011285     0.00000758    -0.00008032     0.00013740    -0.00007651
   MO  69     0.00029335     0.00014814    -0.00014831     0.00006758     0.00007013    -0.00007113     0.00003354     0.00004771
   MO  70    -0.00019169     0.00096154     0.00077042    -0.00016082     0.00013895    -0.00008998     0.00014403     0.00009035
   MO  71    -0.00001782    -0.00001153    -0.00000682    -0.00005128     0.00014986     0.00004685    -0.00010591     0.00006565
   MO  72    -0.00021051     0.00000386    -0.00007300    -0.00004710     0.00004652     0.00001520    -0.00003243    -0.00002315
   MO  73    -0.00022499    -0.00046762    -0.00016870     0.00000434    -0.00023982     0.00004719     0.00005592     0.00001939
   MO  74    -0.00015421    -0.00003181     0.00075015    -0.00010884     0.00039734    -0.00008294     0.00007051     0.00016860
   MO  75     0.00004317    -0.00055448    -0.00050197     0.00002621     0.00014049    -0.00000723    -0.00017070    -0.00036775
   MO  76    -0.00020785    -0.00000198     0.00002238     0.00059647    -0.00015274    -0.00061521    -0.00009904    -0.00013093
   MO  77    -0.00017546     0.00007554    -0.00005423     0.00050209     0.00000142     0.00041716     0.00051409    -0.00022581
   MO  78    -0.00005979     0.00025500    -0.00022764    -0.00008311    -0.00059417     0.00006028    -0.00001532     0.00021114
   MO  79     0.00043658     0.00028876    -0.00028646     0.00007580     0.00016713     0.00001777    -0.00012804    -0.00029313
   MO  80    -0.00011103     0.00006757     0.00006239     0.00002532    -0.00019460     0.00007863    -0.00001017     0.00021407
   MO  81     0.00020324     0.00008025    -0.00017694    -0.00010809     0.00015727    -0.00011595     0.00004354    -0.00023519
   MO  82    -0.00005708     0.00026564     0.00025486    -0.00003617     0.00016083    -0.00003021     0.00000174    -0.00006334

                MO  33         MO  34         MO  35         MO  36         MO  37         MO  38         MO  39         MO  40
   MO   3     0.00173540    -0.00223421     0.00047431     0.00690568     0.00767683     0.00864251     0.00580198    -0.00352450
   MO   4    -0.00107420    -0.00300248     0.00302665    -0.00492005    -0.00294108     0.00576237    -0.00181323     0.00159925
   MO   5    -0.01052183     0.00101553     0.00046227    -0.00573663     0.00520740    -0.00262235    -0.00241357     0.00108103
   MO   6     0.00140018    -0.00336958     0.00220872    -0.00368882    -0.00082468    -0.00409621     0.00219053    -0.00290425
   MO   7    -0.00184069    -0.00382214     0.00721618    -0.00400646     0.00151348    -0.00613703    -0.00602762    -0.00793361
   MO   8    -0.00160595     0.02080621    -0.00262871     0.01116463     0.01358280    -0.00253969     0.00071694     0.00772675
   MO   9    -0.00061942    -0.00306652    -0.00093173    -0.00196098     0.00143468    -0.00091656    -0.00052584     0.00106365
   MO  10     0.00032248     0.00061236     0.00135349     0.00002746    -0.00021899     0.00007798     0.00006486    -0.00003453
   MO  11     0.00029644     0.00017494    -0.00026753     0.00028781     0.00017572    -0.00267712     0.00080722    -0.00024821
   MO  12    -0.00060175    -0.00004923     0.00065719    -0.00058950    -0.00051753     0.00091965    -0.00164006     0.00090170
   MO  13    -0.00083251     0.00011659     0.00034763    -0.00033005     0.00027718     0.00034598    -0.00103757     0.00006585
   MO  14    -0.00011815     0.00016793    -0.00014739     0.00006241     0.00004099    -0.00023156     0.00106168     0.00010753
   MO  15    -0.00021372    -0.00141138     0.00002623    -0.00061970    -0.00054885    -0.00002739    -0.00012383    -0.00025735
   MO  16    -0.00009218    -0.00069691     0.00096877    -0.00068068    -0.00072078     0.00082239    -0.00099296    -0.00016897
   MO  17    -0.00046057    -0.00014103     0.00049660     0.00097807     0.00071903    -0.00388098     0.00027202     0.00079171
   MO  18    -0.00035287    -0.00009562     0.00050246     0.00027049    -0.00008140     0.00145673    -0.00106139     0.00077577
   MO  19    -0.00005573    -0.00010373    -0.00018165     0.00015807     0.00059892     0.00057879     0.00080326    -0.00051858
   MO  20    -0.00119182    -0.00045139     0.00115949    -0.00048529    -0.00016696     0.00042347    -0.00230082     0.00159046
   MO  21     0.00126459     0.00036698    -0.00082987     0.00111582     0.00036556    -0.00169493     0.00175889    -0.00018829
   MO  22    -0.00016188     0.00039246    -0.00018228    -0.00004208    -0.00002174     0.00068944    -0.00048946    -0.00002328
   MO  23     0.00083441    -0.00048425    -0.00001165    -0.00075559    -0.00029966     0.00073396     0.00162111    -0.00185167
   MO  24     0.00052223    -0.00005410    -0.00082632     0.00077413     0.00031694    -0.00205197     0.00392886    -0.00022197
   MO  25    -0.00018764     0.00012923    -0.00029567    -0.00000120     0.00021188    -0.00018408     0.00031254    -0.00092610
   MO  26    -0.00071431    -0.00007772     0.00025869    -0.00015978     0.00008217     0.00072105    -0.00006776     0.00021073
   MO  27    -0.00048286    -0.00029930     0.00060027    -0.00027691    -0.00030901     0.00115903    -0.00070569     0.00092356
   MO  28     0.00004753     0.00060629    -0.00015715     0.00023520     0.00011357    -0.00036462     0.00014152    -0.00013702
   MO  29    -0.00009482     0.00002653    -0.00004521    -0.00037351    -0.00004728     0.00071341     0.00021701    -0.00007630
   MO  30     0.00003474     0.00037789    -0.00016072     0.00008937     0.00010320    -0.00008957     0.00004162     0.00005361
   MO  31    -0.00011924    -0.00038649     0.00041505     0.00004116    -0.00019854     0.00006704    -0.00027925     0.00009944
   MO  32     0.00016601     0.00023379     0.00006899     0.00012487    -0.00006901     0.00001457     0.00001353     0.00033925
   MO  33     0.00171552     0.00001980    -0.00040458     0.00030263    -0.00030508    -0.00023803     0.00087158    -0.00057875
   MO  34     0.00001980     0.00147561    -0.00041180     0.00010333     0.00024181    -0.00019121     0.00006177     0.00018644
   MO  35    -0.00040458    -0.00041180     0.00092420    -0.00002646    -0.00003052     0.00007618    -0.00061453     0.00036049
   MO  36     0.00030263     0.00010333    -0.00002646     0.00121897     0.00061879    -0.00096702     0.00031394     0.00032677
   MO  37    -0.00030508     0.00024181    -0.00003052     0.00061879     0.00104676    -0.00066586     0.00031482     0.00018134
   MO  38    -0.00023803    -0.00019121     0.00007618    -0.00096702    -0.00066586     0.00398786    -0.00075125    -0.00018915
   MO  39     0.00087158     0.00006177    -0.00061453     0.00031394     0.00031482    -0.00075125     0.00240123    -0.00111476
   MO  40    -0.00057875     0.00018644     0.00036049     0.00032677     0.00018134    -0.00018915    -0.00111476     0.00176480
   MO  41     0.00022353     0.00009247    -0.00030193     0.00028272     0.00000358    -0.00054443     0.00088732     0.00045642
   MO  42    -0.00006700     0.00009918     0.00009063     0.00010817    -0.00000640    -0.00011599    -0.00034417     0.00012222
   MO  43    -0.00033391     0.00016223    -0.00009660     0.00026399    -0.00002358     0.00015414    -0.00028044     0.00024976
   MO  44    -0.00004014    -0.00010885     0.00001414     0.00047323     0.00043060    -0.00134534     0.00188228     0.00000071
   MO  45     0.00090978     0.00049025    -0.00025635     0.00026289     0.00018836    -0.00057653     0.00026850    -0.00043182
   MO  46    -0.00003148     0.00028934     0.00024473    -0.00007251     0.00006259     0.00013796    -0.00015340    -0.00010557
   MO  47    -0.00012186     0.00068529     0.00031068     0.00013464     0.00026915     0.00043631    -0.00004924     0.00024705
   MO  48     0.00059640     0.00046426    -0.00042942     0.00010717     0.00016953     0.00008090     0.00044735    -0.00007679
   MO  49    -0.00072698    -0.00026405     0.00018730    -0.00023601    -0.00007129     0.00037823    -0.00046785     0.00002895
   MO  50    -0.00032502     0.00094395    -0.00004280     0.00001492     0.00051526     0.00095073    -0.00048567     0.00034352
   MO  51    -0.00120543     0.00089054    -0.00064722     0.00055583     0.00079615    -0.00117443     0.00023671    -0.00005918
   MO  52     0.00020064    -0.00029886    -0.00038207    -0.00054359    -0.00027648     0.00152163     0.00143364    -0.00112639
   MO  53     0.00091160     0.00060647    -0.00090792     0.00048780     0.00020101    -0.00086402     0.00123297    -0.00054367
   MO  54     0.00018913    -0.00001303    -0.00005385     0.00039147     0.00009879    -0.00105594     0.00126067     0.00068833
   MO  55    -0.00008470    -0.00047043    -0.00006266    -0.00019495     0.00013511    -0.00020078    -0.00002140     0.00014587
   MO  56    -0.00023604     0.00021536    -0.00012550    -0.00003704     0.00003538     0.00105547    -0.00023920     0.00091225
   MO  57    -0.00062273    -0.00007541     0.00019616    -0.00068058     0.00023447    -0.00011174     0.00023484    -0.00013451
   MO  58    -0.00010725     0.00000320    -0.00000605    -0.00004880    -0.00018539     0.00019416    -0.00006096    -0.00038576
   MO  59    -0.00000898    -0.00002709    -0.00038130     0.00001119     0.00015377     0.00011917    -0.00007080    -0.00002190
   MO  60    -0.00037137    -0.00039301     0.00033500     0.00035115     0.00019695    -0.00049891    -0.00038338     0.00052740
   MO  61    -0.00020257    -0.00062978     0.00036999     0.00055242     0.00058970     0.00038953     0.00020444    -0.00011239
   MO  62     0.00047758    -0.00018958     0.00014514     0.00000832     0.00002925    -0.00003352    -0.00005870     0.00009835
   MO  63     0.00082016    -0.00017096    -0.00020698     0.00031981     0.00019645    -0.00061035     0.00033597    -0.00016149
   MO  64    -0.00019730     0.00000695     0.00030580    -0.00018215     0.00013968    -0.00006009    -0.00069036    -0.00023381
   MO  65     0.00000487    -0.00010155     0.00003313     0.00002436     0.00006744     0.00006270     0.00006857    -0.00000041
   MO  66     0.00012901     0.00002109    -0.00000886     0.00001608    -0.00006960    -0.00006057    -0.00001152    -0.00003526
   MO  67    -0.00000014    -0.00013921     0.00000327    -0.00006437    -0.00005158    -0.00002975     0.00003133    -0.00006281
   MO  68    -0.00019114    -0.00031279     0.00012551    -0.00036739    -0.00035938     0.00026671    -0.00030966    -0.00001581
   MO  69     0.00018284    -0.00011270    -0.00003812     0.00020536     0.00025045    -0.00122329     0.00110950    -0.00057131
   MO  70    -0.00014533    -0.00004311     0.00004544     0.00005555     0.00017683    -0.00040924     0.00007388     0.00002191
   MO  71     0.00021111     0.00010915    -0.00016937    -0.00020243    -0.00000817     0.00060406     0.00036019    -0.00036961
   MO  72    -0.00000735    -0.00000876     0.00002550     0.00001835     0.00000511    -0.00000551    -0.00006050     0.00021837
   MO  73    -0.00006426    -0.00009191     0.00010430    -0.00003544     0.00004275    -0.00030358    -0.00005291     0.00013202
   MO  74    -0.00008227     0.00003117     0.00003813    -0.00001674    -0.00003205     0.00038385    -0.00013308     0.00043104
   MO  75     0.00029801     0.00009336    -0.00013103     0.00000195    -0.00001128    -0.00030607     0.00000105    -0.00005082
   MO  76    -0.00003610     0.00007060    -0.00005305    -0.00004648     0.00002519    -0.00001214    -0.00003382     0.00004599
   MO  77     0.00002126     0.00006829    -0.00001745    -0.00001278    -0.00000847     0.00005058    -0.00006336     0.00001687
   MO  78    -0.00016159    -0.00009220     0.00009370    -0.00004316    -0.00007808     0.00024387    -0.00020408    -0.00006038
   MO  79    -0.00046250     0.00012143    -0.00006903     0.00009186     0.00006885    -0.00008106    -0.00013544    -0.00001181
   MO  80     0.00004469     0.00010111    -0.00012556     0.00013470     0.00009185    -0.00009511     0.00006956     0.00015147
   MO  81     0.00013226    -0.00003633    -0.00014187    -0.00006324     0.00001875     0.00056671     0.00032469    -0.00040496
   MO  82    -0.00003195    -0.00006559     0.00006500    -0.00010943    -0.00009458     0.00025122    -0.00005983     0.00015549

                MO  41         MO  42         MO  43         MO  44         MO  45         MO  46         MO  47         MO  48
   MO   3    -0.00021155    -0.00013272    -0.00017376    -0.00101670     0.00232117     0.00064055     0.00031244     0.00273685
   MO   4     0.00324678     0.00161141     0.00390021     0.00424833     0.00843678    -0.00082068    -0.00199360     0.00283991
   MO   5    -0.00218189    -0.00048702     0.00777343    -0.00073852    -0.00401754     0.00016832    -0.00122526    -0.00231670
   MO   6     0.00471924     0.00159305     0.00265757    -0.00376129     0.00695390    -0.00020393    -0.00116071     0.00157565
   MO   7    -0.02523488    -0.00241259    -0.00785414    -0.00260337     0.00878074     0.00190640     0.00051125     0.00292477
   MO   8     0.00147296     0.00426660     0.00300187     0.00501781     0.01482260     0.00409839     0.01713801     0.00726347
   MO   9     0.00137688     0.00290023    -0.00116881     0.00101715    -0.00166724     0.00663977    -0.00627989     0.00102388
   MO  10    -0.00004401    -0.00005095     0.00006949    -0.00006113    -0.00058921     0.00211360     0.00156882    -0.00032327
   MO  11     0.00081607     0.00007640     0.00004149     0.00042626     0.00076508    -0.00006156    -0.00052633    -0.00013116
   MO  12    -0.00083941     0.00009812     0.00007918    -0.00025346    -0.00036542    -0.00002340     0.00042427    -0.00026584
   MO  13    -0.00132768    -0.00007639     0.00032286    -0.00075645     0.00014859     0.00005544     0.00025021    -0.00009047
   MO  14     0.00125901    -0.00012850     0.00085910     0.00182332    -0.00047796    -0.00004424     0.00015900    -0.00013532
   MO  15    -0.00005212    -0.00018049     0.00007355     0.00012820    -0.00080121    -0.00047850    -0.00147384    -0.00042817
   MO  16    -0.00068220    -0.00004581     0.00008376    -0.00158339    -0.00034964     0.00058604     0.00007124    -0.00045598
   MO  17     0.00068871     0.00016278     0.00023904     0.00197787     0.00049430    -0.00028528    -0.00036660    -0.00048570
   MO  18    -0.00059247     0.00000756     0.00007125    -0.00008882    -0.00071898    -0.00037120     0.00076800    -0.00020786
   MO  19    -0.00015529    -0.00025438     0.00005812     0.00055014    -0.00024853     0.00038176     0.00004293     0.00017900
   MO  20    -0.00087940     0.00035712    -0.00001231    -0.00020519    -0.00064986     0.00036945     0.00031408    -0.00064817
   MO  21     0.00121224    -0.00002990    -0.00056960     0.00174975     0.00012489     0.00002569    -0.00028983     0.00029143
   MO  22    -0.00147959    -0.00059459     0.00004461     0.00025329    -0.00082214    -0.00001677     0.00001892     0.00016663
   MO  23    -0.00051868    -0.00051508    -0.00109813     0.00114519     0.00045002     0.00014592     0.00028328     0.00025962
   MO  24     0.00296066    -0.00064211     0.00085725     0.00708351    -0.00129935    -0.00055805    -0.00044044    -0.00008699
   MO  25    -0.00058369     0.00037167     0.00098201     0.00020159     0.00009270     0.00002373    -0.00012917    -0.00001406
   MO  26    -0.00035878    -0.00064493     0.00136124     0.00111956    -0.00167825     0.00017353     0.00049590    -0.00071814
   MO  27    -0.00010276     0.00027027    -0.00015231     0.00128711    -0.00210187    -0.00038030     0.00089319    -0.00000125
   MO  28     0.00000413     0.00013335     0.00024546     0.00007357     0.00029507     0.00056070    -0.00035053     0.00049205
   MO  29    -0.00043010    -0.00031467     0.00033290     0.00031626    -0.00052780    -0.00045256     0.00002666     0.00065845
   MO  30     0.00019082     0.00000858    -0.00012118    -0.00006839     0.00035006    -0.00013448     0.00017058    -0.00080926
   MO  31    -0.00012128     0.00009230    -0.00012002     0.00020001    -0.00063990    -0.00009260    -0.00056800    -0.00060610
   MO  32     0.00017728     0.00017398    -0.00054023     0.00075260    -0.00020679    -0.00020238     0.00075049    -0.00016323
   MO  33     0.00022353    -0.00006700    -0.00033391    -0.00004014     0.00090978    -0.00003148    -0.00012186     0.00059640
   MO  34     0.00009247     0.00009918     0.00016223    -0.00010885     0.00049025     0.00028934     0.00068529     0.00046426
   MO  35    -0.00030193     0.00009063    -0.00009660     0.00001414    -0.00025635     0.00024473     0.00031068    -0.00042942
   MO  36     0.00028272     0.00010817     0.00026399     0.00047323     0.00026289    -0.00007251     0.00013464     0.00010717
   MO  37     0.00000358    -0.00000640    -0.00002358     0.00043060     0.00018836     0.00006259     0.00026915     0.00016953
   MO  38    -0.00054443    -0.00011599     0.00015414    -0.00134534    -0.00057653     0.00013796     0.00043631     0.00008090
   MO  39     0.00088732    -0.00034417    -0.00028044     0.00188228     0.00026850    -0.00015340    -0.00004924     0.00044735
   MO  40     0.00045642     0.00012222     0.00024976     0.00000071    -0.00043182    -0.00010557     0.00024705    -0.00007679
   MO  41     0.00236568     0.00044013     0.00007871     0.00032848     0.00039012    -0.00011043    -0.00014752    -0.00006227
   MO  42     0.00044013     0.00098788     0.00005046    -0.00051717     0.00039491     0.00007598     0.00001430     0.00014910
   MO  43     0.00007871     0.00005046     0.00174040     0.00003060    -0.00030948     0.00004567    -0.00015523    -0.00010016
   MO  44     0.00032848    -0.00051717     0.00003060     0.00592289    -0.00094152    -0.00022347     0.00009693    -0.00013184
   MO  45     0.00039012     0.00039491    -0.00030948    -0.00094152     0.00421099     0.00025671     0.00006061     0.00046791
   MO  46    -0.00011043     0.00007598     0.00004567    -0.00022347     0.00025671     0.00169287    -0.00041229     0.00010480
   MO  47    -0.00014752     0.00001430    -0.00015523     0.00009693     0.00006061    -0.00041229     0.00232160     0.00011726
   MO  48    -0.00006227     0.00014910    -0.00010016    -0.00013184     0.00046791     0.00010480     0.00011726     0.00227948
   MO  49    -0.00016815    -0.00012515     0.00006850    -0.00001570    -0.00050824     0.00024863     0.00035747    -0.00031263
   MO  50    -0.00052027     0.00022914    -0.00006198    -0.00028203     0.00000805     0.00051748     0.00080598     0.00005747
   MO  51     0.00001965    -0.00011543     0.00183366     0.00035552    -0.00025660     0.00009356    -0.00047121    -0.00009397
   MO  52     0.00027780    -0.00061602     0.00058929     0.00139871    -0.00096249    -0.00006875     0.00002810    -0.00009641
   MO  53     0.00109350     0.00017317    -0.00116153    -0.00020137     0.00064472     0.00005377    -0.00036595     0.00020382
   MO  54     0.00196996    -0.00006617     0.00033758     0.00293834    -0.00059816    -0.00025562     0.00014835    -0.00013892
   MO  55     0.00023615     0.00028322    -0.00012128     0.00018191    -0.00010419     0.00032702    -0.00032655    -0.00011246
   MO  56    -0.00005349    -0.00036756    -0.00022978     0.00108569    -0.00117584    -0.00030437     0.00018900     0.00036041
   MO  57     0.00039082     0.00009749    -0.00084243    -0.00014392     0.00043760     0.00000208     0.00002968     0.00011606
   MO  58     0.00011247     0.00117823     0.00014335     0.00052340    -0.00035766    -0.00026154     0.00010880     0.00018945
   MO  59    -0.00004226    -0.00013296    -0.00006402    -0.00012862     0.00008384     0.00001339    -0.00002183    -0.00011729
   MO  60     0.00018837     0.00009118    -0.00001035     0.00009450    -0.00029620    -0.00015771    -0.00028983    -0.00028931
   MO  61    -0.00003022     0.00007385    -0.00009010     0.00002549    -0.00006987     0.00015067     0.00005036    -0.00009566
   MO  62     0.00004658     0.00010993     0.00020611     0.00004948     0.00003351     0.00011194     0.00004823    -0.00022480
   MO  63    -0.00021924    -0.00012649     0.00006331    -0.00013346     0.00059736    -0.00027581    -0.00022498     0.00069862
   MO  64    -0.00141034    -0.00019136    -0.00057256    -0.00020001    -0.00018384     0.00016358     0.00021427    -0.00018313
   MO  65     0.00008200    -0.00001247    -0.00008663    -0.00000157     0.00015637    -0.00054849     0.00041006    -0.00051562
   MO  66    -0.00004357    -0.00001455    -0.00000586     0.00001219    -0.00012377     0.00019791    -0.00044106    -0.00046051
   MO  67    -0.00000148    -0.00001325    -0.00001151    -0.00000406    -0.00000929     0.00000782     0.00044851     0.00011843
   MO  68    -0.00013229    -0.00002637    -0.00000115    -0.00012133    -0.00011487    -0.00008584    -0.00020548    -0.00020087
   MO  69    -0.00030966    -0.00032668    -0.00026367     0.00205483     0.00000741    -0.00005623    -0.00007286     0.00003229
   MO  70    -0.00004381    -0.00011256    -0.00014579    -0.00043342    -0.00161622    -0.00020682     0.00039805    -0.00009442
   MO  71    -0.00003272    -0.00012259    -0.00028839    -0.00011765     0.00002968     0.00001267     0.00019798     0.00016357
   MO  72    -0.00001554    -0.00023132    -0.00000827     0.00003764     0.00005056     0.00001006    -0.00003830    -0.00000505
   MO  73     0.00028575     0.00011037    -0.00032768     0.00010103     0.00046114     0.00003504    -0.00012001    -0.00016791
   MO  74    -0.00021019     0.00001179    -0.00005924     0.00029676    -0.00090410    -0.00032194     0.00023181     0.00040150
   MO  75     0.00029450     0.00028356    -0.00002313    -0.00127471     0.00078840    -0.00016958    -0.00014511     0.00054470
   MO  76     0.00001912     0.00008381     0.00003365     0.00000121    -0.00011730     0.00079822    -0.00040645     0.00100322
   MO  77     0.00004085     0.00005789    -0.00005649    -0.00002224    -0.00011476     0.00055340    -0.00087172    -0.00046433
   MO  78     0.00012206    -0.00000026     0.00004967    -0.00022351     0.00035124     0.00033052     0.00013897    -0.00081735
   MO  79    -0.00015213    -0.00007950     0.00046750    -0.00006712    -0.00027325     0.00001982    -0.00011823    -0.00012497
   MO  80     0.00037216     0.00003193    -0.00002068    -0.00018213    -0.00058963    -0.00007901     0.00006179    -0.00009096
   MO  81    -0.00007292    -0.00010618     0.00002957    -0.00011149     0.00002914     0.00004903     0.00001256     0.00015327
   MO  82     0.00032135     0.00006963     0.00022368     0.00025619     0.00018824    -0.00002891    -0.00004745    -0.00000833

                MO  49         MO  50         MO  51         MO  52         MO  53         MO  54         MO  55         MO  56
   MO   3    -0.00133355     0.00163771    -0.00109406     0.00541067     0.00075449    -0.00259478    -0.00036746    -0.00000550
   MO   4    -0.00349598    -0.00165576    -0.00137697     0.00390706    -0.00802948     0.00553502     0.00012516     0.00030826
   MO   5     0.00105568    -0.00066729     0.01212055     0.00051692    -0.00219190    -0.00236117    -0.00024554     0.00058549
   MO   6    -0.00128546     0.00187119    -0.00441896     0.00574270    -0.01328355     0.00619440     0.00033201    -0.00743258
   MO   7    -0.00074938     0.00445344    -0.00686149    -0.00778809    -0.01045869    -0.01634538    -0.00057805    -0.00750125
   MO   8    -0.00213356     0.02652017     0.00694130    -0.00529396    -0.00055408     0.00448992    -0.00050751     0.00068068
   MO   9    -0.00297010     0.00445833     0.00041695     0.00091282    -0.00006854    -0.00162728     0.02612201    -0.00007456
   MO  10    -0.00107729    -0.00088394     0.00006403    -0.00009713     0.00022034    -0.00013764     0.00012874     0.00017642
   MO  11    -0.00026136    -0.00095021     0.00089420    -0.00052488     0.00052425     0.00087116     0.00002064    -0.00120641
   MO  12     0.00046117     0.00090948    -0.00120700    -0.00042691    -0.00195406    -0.00016842    -0.00000980     0.00020046
   MO  13     0.00033910     0.00085941     0.00045251    -0.00032310    -0.00155124    -0.00111373    -0.00029762    -0.00039545
   MO  14    -0.00006810    -0.00025958     0.00084016     0.00141880    -0.00053499     0.00209615    -0.00020701     0.00000673
   MO  15     0.00034281    -0.00100530    -0.00046823     0.00063281    -0.00098343     0.00017839     0.00089531    -0.00041419
   MO  16     0.00000189    -0.00066008    -0.00101799     0.00002427    -0.00151020    -0.00084228    -0.00029331    -0.00066752
   MO  17    -0.00013684    -0.00125191     0.00119339    -0.00153436    -0.00038384     0.00167053    -0.00063338    -0.00092540
   MO  18     0.00066451     0.00079030    -0.00107233    -0.00027375    -0.00079899    -0.00012807    -0.00158752     0.00080437
   MO  19    -0.00023524     0.00036631     0.00069986     0.00096031     0.00030322    -0.00024533     0.00125992     0.00029621
   MO  20     0.00042927     0.00122504    -0.00126718    -0.00109013    -0.00258668    -0.00025656     0.00210425     0.00018133
   MO  21    -0.00079090    -0.00090100     0.00059754    -0.00032113     0.00338000     0.00091401     0.00123910     0.00084518
   MO  22    -0.00002829    -0.00020249     0.00093022    -0.00038804     0.00082162    -0.00141343    -0.00089791     0.00146419
   MO  23     0.00022561     0.00051987    -0.00214969     0.00214408    -0.00079773     0.00046109    -0.00032092    -0.00117737
   MO  24    -0.00046339    -0.00174857     0.00173126     0.00312348     0.00110272     0.00554905     0.00009262     0.00099277
   MO  25    -0.00032131    -0.00003938     0.00198935     0.00061125    -0.00022057    -0.00078837     0.00110740    -0.00063951
   MO  26     0.00101166     0.00144736     0.00037368     0.00200152    -0.00290271     0.00143809    -0.00005778    -0.00064177
   MO  27     0.00039121     0.00099112    -0.00191848     0.00052568    -0.00186189     0.00128021    -0.00016671     0.00092368
   MO  28     0.00014922     0.00023012     0.00069296    -0.00047311     0.00043747     0.00019479    -0.00232788    -0.00013079
   MO  29    -0.00082718     0.00023342     0.00038471     0.00079647    -0.00119314    -0.00003001    -0.00002092     0.00062616
   MO  30    -0.00063259     0.00057876     0.00024976    -0.00023436     0.00046567     0.00010795    -0.00025331     0.00002951
   MO  31    -0.00003315    -0.00003804    -0.00071747    -0.00001333    -0.00026601     0.00007826    -0.00111763    -0.00006016
   MO  32     0.00041268     0.00068867    -0.00093914    -0.00037865     0.00073574     0.00044032     0.00080957     0.00047707
   MO  33    -0.00072698    -0.00032502    -0.00120543     0.00020064     0.00091160     0.00018913    -0.00008470    -0.00023604
   MO  34    -0.00026405     0.00094395     0.00089054    -0.00029886     0.00060647    -0.00001303    -0.00047043     0.00021536
   MO  35     0.00018730    -0.00004280    -0.00064722    -0.00038207    -0.00090792    -0.00005385    -0.00006266    -0.00012550
   MO  36    -0.00023601     0.00001492     0.00055583    -0.00054359     0.00048780     0.00039147    -0.00019495    -0.00003704
   MO  37    -0.00007129     0.00051526     0.00079615    -0.00027648     0.00020101     0.00009879     0.00013511     0.00003538
   MO  38     0.00037823     0.00095073    -0.00117443     0.00152163    -0.00086402    -0.00105594    -0.00020078     0.00105547
   MO  39    -0.00046785    -0.00048567     0.00023671     0.00143364     0.00123297     0.00126067    -0.00002140    -0.00023920
   MO  40     0.00002895     0.00034352    -0.00005918    -0.00112639    -0.00054367     0.00068833     0.00014587     0.00091225
   MO  41    -0.00016815    -0.00052027     0.00001965     0.00027780     0.00109350     0.00196996     0.00023615    -0.00005349
   MO  42    -0.00012515     0.00022914    -0.00011543    -0.00061602     0.00017317    -0.00006617     0.00028322    -0.00036756
   MO  43     0.00006850    -0.00006198     0.00183366     0.00058929    -0.00116153     0.00033758    -0.00012128    -0.00022978
   MO  44    -0.00001570    -0.00028203     0.00035552     0.00139871    -0.00020137     0.00293834     0.00018191     0.00108569
   MO  45    -0.00050824     0.00000805    -0.00025660    -0.00096249     0.00064472    -0.00059816    -0.00010419    -0.00117584
   MO  46     0.00024863     0.00051748     0.00009356    -0.00006875     0.00005377    -0.00025562     0.00032702    -0.00030437
   MO  47     0.00035747     0.00080598    -0.00047121     0.00002810    -0.00036595     0.00014835    -0.00032655     0.00018900
   MO  48    -0.00031263     0.00005747    -0.00009397    -0.00009641     0.00020382    -0.00013892    -0.00011246     0.00036041
   MO  49     0.00217903     0.00009221    -0.00026651     0.00022174    -0.00045585     0.00008444    -0.00033276    -0.00034380
   MO  50     0.00009221     0.00318731    -0.00060926     0.00020288    -0.00114949    -0.00011661     0.00023776    -0.00017412
   MO  51    -0.00026651    -0.00060926     0.00528160     0.00005105     0.00064701    -0.00033870     0.00003626     0.00010888
   MO  52     0.00022174     0.00020288     0.00005105     0.00336804    -0.00103206     0.00134122     0.00018773    -0.00019024
   MO  53    -0.00045585    -0.00114949     0.00064701    -0.00103206     0.00527738    -0.00051097    -0.00003305     0.00097041
   MO  54     0.00008444    -0.00011661    -0.00033870     0.00134122    -0.00051097     0.00406023    -0.00016594     0.00006808
   MO  55    -0.00033276     0.00023776     0.00003626     0.00018773    -0.00003305    -0.00016594     0.00474338    -0.00001749
   MO  56    -0.00034380    -0.00017412     0.00010888    -0.00019024     0.00097041     0.00006808    -0.00001749     0.00302095
   MO  57    -0.00002704     0.00006952    -0.00053103     0.00006511    -0.00053692     0.00022410    -0.00004106    -0.00047760
   MO  58    -0.00016885    -0.00001851     0.00003956    -0.00005763     0.00016582    -0.00002843     0.00000708     0.00006966
   MO  59     0.00020255     0.00015770     0.00002714     0.00009419     0.00001808     0.00006374     0.00000345    -0.00007500
   MO  60     0.00010191    -0.00058190    -0.00007304    -0.00023322    -0.00015313     0.00002420     0.00006575     0.00001174
   MO  61     0.00016201     0.00018597    -0.00009765     0.00017289     0.00001933    -0.00019361    -0.00007241    -0.00005004
   MO  62     0.00026487     0.00012426    -0.00080785     0.00001673    -0.00005440     0.00025318     0.00004550    -0.00030503
   MO  63    -0.00070193    -0.00036789    -0.00014411    -0.00011668    -0.00018755    -0.00018891     0.00004411    -0.00012286
   MO  64     0.00031306     0.00046869    -0.00041618    -0.00055482    -0.00032747    -0.00115316    -0.00004270    -0.00011076
   MO  65    -0.00028742     0.00003640    -0.00004977     0.00007877     0.00003448    -0.00003239     0.00058390     0.00002403
   MO  66    -0.00066914     0.00023627    -0.00003689     0.00001008     0.00002790     0.00001809    -0.00057357    -0.00002052
   MO  67     0.00017293    -0.00041893    -0.00000411     0.00008698     0.00002427    -0.00008268     0.00048516    -0.00001316
   MO  68     0.00019855    -0.00018517    -0.00011682     0.00014408    -0.00015197    -0.00023713     0.00003097     0.00001239
   MO  69    -0.00009186    -0.00025463     0.00012476     0.00029717     0.00000271     0.00046371     0.00002617    -0.00020832
   MO  70     0.00031651     0.00061287    -0.00037525     0.00032660    -0.00060079     0.00033835     0.00001224    -0.00041043
   MO  71    -0.00003595     0.00026903    -0.00035602     0.00052019     0.00028034    -0.00015832    -0.00002582     0.00010651
   MO  72    -0.00005915    -0.00008668     0.00000920    -0.00009365    -0.00002720     0.00003492    -0.00000599     0.00020884
   MO  73     0.00002924    -0.00030039    -0.00013807    -0.00040146     0.00033162     0.00010419     0.00001728     0.00000721
   MO  74    -0.00030927     0.00013189    -0.00023680    -0.00009310    -0.00018345    -0.00001855    -0.00000045     0.00079094
   MO  75    -0.00060456    -0.00025808     0.00008618    -0.00042783     0.00021734    -0.00029309     0.00000980    -0.00051215
   MO  76     0.00071711    -0.00003519    -0.00000335    -0.00010284     0.00003330     0.00007977    -0.00024585     0.00001129
   MO  77    -0.00051718     0.00061915    -0.00005605    -0.00002190     0.00007102     0.00005605    -0.00038503    -0.00003273
   MO  78     0.00092052     0.00007836    -0.00016414     0.00009409     0.00008424    -0.00006220     0.00007524    -0.00029526
   MO  79     0.00005668    -0.00011393     0.00139910     0.00014322    -0.00028861    -0.00021109     0.00001544    -0.00008341
   MO  80     0.00008116    -0.00004784     0.00025037     0.00011224     0.00087031     0.00052850     0.00002618     0.00003317
   MO  81    -0.00005127     0.00010244     0.00006848     0.00075923     0.00014932    -0.00025618    -0.00004081     0.00000216
   MO  82    -0.00002168    -0.00005340    -0.00008076     0.00053120    -0.00064354     0.00086391     0.00002663     0.00017833

                MO  57         MO  58         MO  59         MO  60         MO  61         MO  62         MO  63         MO  64
   MO   3    -0.00046191    -0.00009121    -0.00018329     0.00430751     0.00294934     0.00038746    -0.00005428    -0.00005773
   MO   4     0.00229284     0.00125157    -0.00076835     0.00013739    -0.00091126    -0.00072647     0.00233495    -0.00140853
   MO   5     0.00128762     0.00007329     0.00018466    -0.00029991    -0.00052871     0.00360309     0.00124194    -0.00022940
   MO   6     0.00420614    -0.00077822    -0.00007397     0.00018985    -0.00019817     0.00064959     0.00172079     0.00080227
   MO   7     0.00326392    -0.00176106     0.00016602    -0.00009292     0.00034792    -0.00002554     0.00183376     0.00402295
   MO   8    -0.00067771     0.00006199     0.00153539    -0.00826585     0.00163963    -0.00072825    -0.00107919     0.00045236
   MO   9    -0.00007661     0.00028952     0.00001366     0.00019811    -0.00030998     0.00019289     0.00025560     0.00000998
   MO  10    -0.00008195    -0.00008892    -0.00055288    -0.00009302     0.00007255     0.00003853    -0.00029213     0.00006017
   MO  11     0.00039711    -0.00021477    -0.00005300     0.00026017    -0.00056110     0.00000695     0.00040323    -0.00023242
   MO  12     0.00001090     0.00005852    -0.00004244     0.00013263    -0.00043223     0.00015964    -0.00033768     0.00047426
   MO  13     0.00029574    -0.00013803     0.00007425    -0.00003934     0.00004396    -0.00003267    -0.00000035     0.00063893
   MO  14     0.00013199     0.00007753    -0.00002937    -0.00006654    -0.00012735     0.00003789    -0.00012920    -0.00083469
   MO  15     0.00026406    -0.00005027     0.00004334     0.00051021    -0.00009766     0.00008650     0.00015201    -0.00006175
   MO  16     0.00027483    -0.00015898    -0.00017839     0.00027436    -0.00000287     0.00014514    -0.00001152     0.00027657
   MO  17     0.00048000    -0.00008120    -0.00016498     0.00104328    -0.00018155     0.00003125     0.00040681    -0.00007442
   MO  18    -0.00056965     0.00020902     0.00002829     0.00016764     0.00048370     0.00022166    -0.00039886     0.00027412
   MO  19    -0.00000220    -0.00017797     0.00011925    -0.00026336     0.00067725    -0.00012583     0.00006676    -0.00003160
   MO  20     0.00019709    -0.00003572    -0.00007930     0.00052855    -0.00020122     0.00021134    -0.00045258     0.00068607
   MO  21    -0.00096185     0.00003485    -0.00004857    -0.00012201     0.00000411     0.00001640     0.00024596    -0.00045517
   MO  22    -0.00082975    -0.00013119     0.00005840    -0.00019175    -0.00013979    -0.00017543    -0.00006518     0.00069315
   MO  23     0.00097570    -0.00011125     0.00010584    -0.00062203     0.00022463    -0.00004967     0.00016830     0.00042396
   MO  24    -0.00012491     0.00038668    -0.00011540     0.00006523    -0.00017344    -0.00004662    -0.00005656    -0.00179017
   MO  25    -0.00041103     0.00101624    -0.00008238    -0.00032523     0.00001257    -0.00001270     0.00016191     0.00006735
   MO  26    -0.00040646    -0.00068190     0.00019316    -0.00022941     0.00007156     0.00048256    -0.00039423    -0.00001536
   MO  27     0.00008002     0.00098333     0.00017995     0.00021406     0.00000832     0.00016348    -0.00041999     0.00012397
   MO  28    -0.00003319     0.00024663    -0.00001783    -0.00059554     0.00004032    -0.00004118    -0.00007743    -0.00002704
   MO  29     0.00035501    -0.00004808    -0.00050950    -0.00043153    -0.00011334    -0.00047470     0.00047989    -0.00011433
   MO  30     0.00000774    -0.00006762     0.00034407    -0.00046323    -0.00011493    -0.00011587    -0.00021323    -0.00000530
   MO  31    -0.00002174     0.00027666    -0.00047167     0.00075403     0.00024250     0.00025051    -0.00018993     0.00012959
   MO  32    -0.00044324     0.00018798    -0.00046478    -0.00051116    -0.00011235     0.00027660    -0.00049606     0.00022988
   MO  33    -0.00062273    -0.00010725    -0.00000898    -0.00037137    -0.00020257     0.00047758     0.00082016    -0.00019730
   MO  34    -0.00007541     0.00000320    -0.00002709    -0.00039301    -0.00062978    -0.00018958    -0.00017096     0.00000695
   MO  35     0.00019616    -0.00000605    -0.00038130     0.00033500     0.00036999     0.00014514    -0.00020698     0.00030580
   MO  36    -0.00068058    -0.00004880     0.00001119     0.00035115     0.00055242     0.00000832     0.00031981    -0.00018215
   MO  37     0.00023447    -0.00018539     0.00015377     0.00019695     0.00058970     0.00002925     0.00019645     0.00013968
   MO  38    -0.00011174     0.00019416     0.00011917    -0.00049891     0.00038953    -0.00003352    -0.00061035    -0.00006009
   MO  39     0.00023484    -0.00006096    -0.00007080    -0.00038338     0.00020444    -0.00005870     0.00033597    -0.00069036
   MO  40    -0.00013451    -0.00038576    -0.00002190     0.00052740    -0.00011239     0.00009835    -0.00016149    -0.00023381
   MO  41     0.00039082     0.00011247    -0.00004226     0.00018837    -0.00003022     0.00004658    -0.00021924    -0.00141034
   MO  42     0.00009749     0.00117823    -0.00013296     0.00009118     0.00007385     0.00010993    -0.00012649    -0.00019136
   MO  43    -0.00084243     0.00014335    -0.00006402    -0.00001035    -0.00009010     0.00020611     0.00006331    -0.00057256
   MO  44    -0.00014392     0.00052340    -0.00012862     0.00009450     0.00002549     0.00004948    -0.00013346    -0.00020001
   MO  45     0.00043760    -0.00035766     0.00008384    -0.00029620    -0.00006987     0.00003351     0.00059736    -0.00018384
   MO  46     0.00000208    -0.00026154     0.00001339    -0.00015771     0.00015067     0.00011194    -0.00027581     0.00016358
   MO  47     0.00002968     0.00010880    -0.00002183    -0.00028983     0.00005036     0.00004823    -0.00022498     0.00021427
   MO  48     0.00011606     0.00018945    -0.00011729    -0.00028931    -0.00009566    -0.00022480     0.00069862    -0.00018313
   MO  49    -0.00002704    -0.00016885     0.00020255     0.00010191     0.00016201     0.00026487    -0.00070193     0.00031306
   MO  50     0.00006952    -0.00001851     0.00015770    -0.00058190     0.00018597     0.00012426    -0.00036789     0.00046869
   MO  51    -0.00053103     0.00003956     0.00002714    -0.00007304    -0.00009765    -0.00080785    -0.00014411    -0.00041618
   MO  52     0.00006511    -0.00005763     0.00009419    -0.00023322     0.00017289     0.00001673    -0.00011668    -0.00055482
   MO  53    -0.00053692     0.00016582     0.00001808    -0.00015313     0.00001933    -0.00005440    -0.00018755    -0.00032747
   MO  54     0.00022410    -0.00002843     0.00006374     0.00002420    -0.00019361     0.00025318    -0.00018891    -0.00115316
   MO  55    -0.00004106     0.00000708     0.00000345     0.00006575    -0.00007241     0.00004550     0.00004411    -0.00004270
   MO  56    -0.00047760     0.00006966    -0.00007500     0.00001174    -0.00005004    -0.00030503    -0.00012286    -0.00011076
   MO  57     0.00207464    -0.00005517     0.00000826     0.00008354     0.00000095    -0.00053743    -0.00009716     0.00000117
   MO  58    -0.00005517     0.00245122    -0.00012817    -0.00001027     0.00001146     0.00000788    -0.00014854    -0.00009700
   MO  59     0.00000826    -0.00012817     0.00140437    -0.00002523    -0.00002271    -0.00003116    -0.00000547     0.00003148
   MO  60     0.00008354    -0.00001027    -0.00002523     0.00177532     0.00035671     0.00005726     0.00009481    -0.00007756
   MO  61     0.00000095     0.00001146    -0.00002271     0.00035671     0.00144279     0.00006983    -0.00000395     0.00004111
   MO  62    -0.00053743     0.00000788    -0.00003116     0.00005726     0.00006983     0.00130789     0.00010440     0.00000007
   MO  63    -0.00009716    -0.00014854    -0.00000547     0.00009481    -0.00000395     0.00010440     0.00136463    -0.00005603
   MO  64     0.00000117    -0.00009700     0.00003148    -0.00007756     0.00004111     0.00000007    -0.00005603     0.00166411
   MO  65    -0.00000835    -0.00001777     0.00000770    -0.00003108     0.00000343    -0.00001349    -0.00002140    -0.00001199
   MO  66     0.00000864     0.00001730    -0.00001787     0.00002150     0.00001401    -0.00000579     0.00001754    -0.00000123
   MO  67    -0.00002272     0.00002540    -0.00001414     0.00004783    -0.00000616     0.00000972     0.00003550    -0.00001087
   MO  68     0.00005893     0.00004037    -0.00002189     0.00014884    -0.00001097    -0.00002062    -0.00009129     0.00006784
   MO  69     0.00014450     0.00001783    -0.00005313    -0.00001337    -0.00000853    -0.00005101     0.00026923     0.00003820
   MO  70     0.00016542    -0.00013305     0.00003096     0.00007378     0.00006977     0.00003796     0.00003327     0.00019587
   MO  71     0.00019684    -0.00004421     0.00004079    -0.00028299    -0.00006064    -0.00002486    -0.00004081     0.00002762
   MO  72    -0.00002308    -0.00035999     0.00001447     0.00004744    -0.00001513    -0.00002786     0.00003798    -0.00001669
   MO  73     0.00039809    -0.00000480     0.00000522     0.00018686     0.00007498    -0.00003595    -0.00008868    -0.00002795
   MO  74    -0.00007672     0.00023651    -0.00007904     0.00004430    -0.00001839    -0.00007401     0.00001974    -0.00001537
   MO  75     0.00029728    -0.00001945    -0.00004331    -0.00001514    -0.00013932    -0.00010594     0.00044084    -0.00022838
   MO  76     0.00007095    -0.00001242     0.00000579    -0.00001677     0.00003401     0.00002857    -0.00001374    -0.00000662
   MO  77     0.00005281    -0.00003756     0.00001231    -0.00004335     0.00005477     0.00000386    -0.00011226     0.00001265
   MO  78    -0.00013647    -0.00016520     0.00011438     0.00001504     0.00010430     0.00021079    -0.00036902     0.00006093
   MO  79    -0.00008729    -0.00004833     0.00003181     0.00003938    -0.00002146    -0.00050143    -0.00017279    -0.00001884
   MO  80    -0.00022159    -0.00005387     0.00004588    -0.00001253     0.00003193     0.00010884    -0.00013420    -0.00010113
   MO  81    -0.00002903    -0.00000650     0.00000261     0.00017979     0.00031640    -0.00004511     0.00010718    -0.00008798
   MO  82     0.00012168     0.00012774    -0.00004315     0.00000245    -0.00006107     0.00001802    -0.00003962    -0.00022605

                MO  65         MO  66         MO  67         MO  68         MO  69         MO  70         MO  71         MO  72
   MO   3     0.00067212    -0.00066887    -0.00098172    -0.00401741    -0.00001515    -0.00130536     0.00540679    -0.00082422
   MO   4    -0.00006572     0.00004712     0.00000938     0.00014991     0.00033233    -0.00511920    -0.00084306    -0.00013484
   MO   5    -0.00032781     0.00021195     0.00003927     0.00066979    -0.00013210    -0.00040439    -0.00350184     0.00007168
   MO   6     0.00001007    -0.00004889    -0.00005056     0.00039779    -0.00250179     0.00114580    -0.00019623    -0.00018157
   MO   7    -0.00036941     0.00032150     0.00013030     0.00006638     0.01103641    -0.00072048    -0.00046976    -0.00004678
   MO   8     0.00098652    -0.00112930    -0.00113879    -0.00452375    -0.00017130    -0.00065492     0.00143248    -0.00023362
   MO   9     0.00500213    -0.00291056     0.00317805     0.00040767     0.00007434     0.00007919    -0.00010077    -0.00020654
   MO  10    -0.00020993     0.00053187     0.00092377    -0.00037383    -0.00005844    -0.00003728     0.00003267     0.00007809
   MO  11    -0.00006491     0.00007052     0.00002895    -0.00010601     0.00061110     0.00026266    -0.00031098    -0.00002613
   MO  12    -0.00001914    -0.00001579    -0.00004487     0.00026994    -0.00034895    -0.00001188    -0.00017734     0.00005888
   MO  13    -0.00002092    -0.00002049    -0.00008408     0.00010134    -0.00001697     0.00007124    -0.00007868    -0.00000286
   MO  14    -0.00001061     0.00004739     0.00005207    -0.00009126     0.00021311     0.00003757    -0.00002312    -0.00000243
   MO  15     0.00000120    -0.00006897    -0.00010274     0.00035110     0.00013651     0.00016198    -0.00017114    -0.00000194
   MO  16    -0.00012984     0.00022529     0.00027658     0.00016158    -0.00038716     0.00016244    -0.00007731     0.00001628
   MO  17    -0.00011379     0.00013107    -0.00002118    -0.00012710     0.00104903     0.00013304    -0.00078625     0.00006052
   MO  18    -0.00003366     0.00002388    -0.00018583     0.00007053    -0.00046359     0.00002099    -0.00000185     0.00003034
   MO  19     0.00014804    -0.00010366     0.00012907    -0.00017413     0.00022396    -0.00001902     0.00028410    -0.00000600
   MO  20     0.00016691    -0.00014648     0.00012139     0.00028658    -0.00038402     0.00012945    -0.00040946     0.00010150
   MO  21     0.00015385    -0.00004598     0.00013471    -0.00034588     0.00058944    -0.00022561     0.00005796     0.00003478
   MO  22    -0.00015905     0.00009865    -0.00004719     0.00007500     0.00005800    -0.00038912     0.00006176     0.00012616
   MO  23    -0.00000877     0.00003146     0.00004045     0.00000252     0.00095617     0.00031806     0.00060866    -0.00014965
   MO  24     0.00001676     0.00003815     0.00005968    -0.00026927     0.00174774    -0.00032283    -0.00004784     0.00003842
   MO  25     0.00020471    -0.00009937     0.00015555    -0.00000663     0.00029335    -0.00019169    -0.00001782    -0.00021051
   MO  26    -0.00002838     0.00000106    -0.00008400     0.00000114     0.00014814     0.00096154    -0.00001153     0.00000386
   MO  27    -0.00003750    -0.00000498     0.00001122     0.00005620    -0.00014831     0.00077042    -0.00000682    -0.00007300
   MO  28    -0.00076171     0.00035369    -0.00038570    -0.00011285     0.00006758    -0.00016082    -0.00005128    -0.00004710
   MO  29     0.00010817     0.00011747    -0.00008843     0.00000758     0.00007013     0.00013895     0.00014986     0.00004652
   MO  30     0.00046627     0.00057750    -0.00026551    -0.00008032    -0.00007113    -0.00008998     0.00004685     0.00001520
   MO  31    -0.00012407     0.00057768    -0.00038075     0.00013740     0.00003354     0.00014403    -0.00010591    -0.00003243
   MO  32     0.00030605    -0.00032635     0.00009620    -0.00007651     0.00004771     0.00009035     0.00006565    -0.00002315
   MO  33     0.00000487     0.00012901    -0.00000014    -0.00019114     0.00018284    -0.00014533     0.00021111    -0.00000735
   MO  34    -0.00010155     0.00002109    -0.00013921    -0.00031279    -0.00011270    -0.00004311     0.00010915    -0.00000876
   MO  35     0.00003313    -0.00000886     0.00000327     0.00012551    -0.00003812     0.00004544    -0.00016937     0.00002550
   MO  36     0.00002436     0.00001608    -0.00006437    -0.00036739     0.00020536     0.00005555    -0.00020243     0.00001835
   MO  37     0.00006744    -0.00006960    -0.00005158    -0.00035938     0.00025045     0.00017683    -0.00000817     0.00000511
   MO  38     0.00006270    -0.00006057    -0.00002975     0.00026671    -0.00122329    -0.00040924     0.00060406    -0.00000551
   MO  39     0.00006857    -0.00001152     0.00003133    -0.00030966     0.00110950     0.00007388     0.00036019    -0.00006050
   MO  40    -0.00000041    -0.00003526    -0.00006281    -0.00001581    -0.00057131     0.00002191    -0.00036961     0.00021837
   MO  41     0.00008200    -0.00004357    -0.00000148    -0.00013229    -0.00030966    -0.00004381    -0.00003272    -0.00001554
   MO  42    -0.00001247    -0.00001455    -0.00001325    -0.00002637    -0.00032668    -0.00011256    -0.00012259    -0.00023132
   MO  43    -0.00008663    -0.00000586    -0.00001151    -0.00000115    -0.00026367    -0.00014579    -0.00028839    -0.00000827
   MO  44    -0.00000157     0.00001219    -0.00000406    -0.00012133     0.00205483    -0.00043342    -0.00011765     0.00003764
   MO  45     0.00015637    -0.00012377    -0.00000929    -0.00011487     0.00000741    -0.00161622     0.00002968     0.00005056
   MO  46    -0.00054849     0.00019791     0.00000782    -0.00008584    -0.00005623    -0.00020682     0.00001267     0.00001006
   MO  47     0.00041006    -0.00044106     0.00044851    -0.00020548    -0.00007286     0.00039805     0.00019798    -0.00003830
   MO  48    -0.00051562    -0.00046051     0.00011843    -0.00020087     0.00003229    -0.00009442     0.00016357    -0.00000505
   MO  49    -0.00028742    -0.00066914     0.00017293     0.00019855    -0.00009186     0.00031651    -0.00003595    -0.00005915
   MO  50     0.00003640     0.00023627    -0.00041893    -0.00018517    -0.00025463     0.00061287     0.00026903    -0.00008668
   MO  51    -0.00004977    -0.00003689    -0.00000411    -0.00011682     0.00012476    -0.00037525    -0.00035602     0.00000920
   MO  52     0.00007877     0.00001008     0.00008698     0.00014408     0.00029717     0.00032660     0.00052019    -0.00009365
   MO  53     0.00003448     0.00002790     0.00002427    -0.00015197     0.00000271    -0.00060079     0.00028034    -0.00002720
   MO  54    -0.00003239     0.00001809    -0.00008268    -0.00023713     0.00046371     0.00033835    -0.00015832     0.00003492
   MO  55     0.00058390    -0.00057357     0.00048516     0.00003097     0.00002617     0.00001224    -0.00002582    -0.00000599
   MO  56     0.00002403    -0.00002052    -0.00001316     0.00001239    -0.00020832    -0.00041043     0.00010651     0.00020884
   MO  57    -0.00000835     0.00000864    -0.00002272     0.00005893     0.00014450     0.00016542     0.00019684    -0.00002308
   MO  58    -0.00001777     0.00001730     0.00002540     0.00004037     0.00001783    -0.00013305    -0.00004421    -0.00035999
   MO  59     0.00000770    -0.00001787    -0.00001414    -0.00002189    -0.00005313     0.00003096     0.00004079     0.00001447
   MO  60    -0.00003108     0.00002150     0.00004783     0.00014884    -0.00001337     0.00007378    -0.00028299     0.00004744
   MO  61     0.00000343     0.00001401    -0.00000616    -0.00001097    -0.00000853     0.00006977    -0.00006064    -0.00001513
   MO  62    -0.00001349    -0.00000579     0.00000972    -0.00002062    -0.00005101     0.00003796    -0.00002486    -0.00002786
   MO  63    -0.00002140     0.00001754     0.00003550    -0.00009129     0.00026923     0.00003327    -0.00004081     0.00003798
   MO  64    -0.00001199    -0.00000123    -0.00001087     0.00006784     0.00003820     0.00019587     0.00002762    -0.00001669
   MO  65     0.00087745     0.00001638     0.00003626    -0.00002259    -0.00002210    -0.00000018     0.00003786     0.00000260
   MO  66     0.00001638     0.00107964    -0.00032949     0.00001398     0.00003536     0.00000670    -0.00000359     0.00001317
   MO  67     0.00003626    -0.00032949     0.00070422     0.00003221     0.00003537     0.00001035     0.00001762     0.00000056
   MO  68    -0.00002259     0.00001398     0.00003221     0.00067158    -0.00013400    -0.00015128     0.00001668     0.00000133
   MO  69    -0.00002210     0.00003536     0.00003537    -0.00013400     0.00195547     0.00002175    -0.00002079    -0.00001865
   MO  70    -0.00000018     0.00000670     0.00001035    -0.00015128     0.00002175     0.00176907     0.00006488    -0.00009438
   MO  71     0.00003786    -0.00000359     0.00001762     0.00001668    -0.00002079     0.00006488     0.00056068    -0.00003745
   MO  72     0.00000260     0.00001317     0.00000056     0.00000133    -0.00001865    -0.00009438    -0.00003745     0.00023965
   MO  73     0.00002739    -0.00000630     0.00000051     0.00003521     0.00003645    -0.00026900    -0.00004542     0.00002868
   MO  74     0.00000940    -0.00000732    -0.00000195    -0.00000790    -0.00008136     0.00025672     0.00003233     0.00004803
   MO  75     0.00001213     0.00000351     0.00000111    -0.00004712    -0.00047973    -0.00010757    -0.00002168    -0.00001715
   MO  76    -0.00062959    -0.00025765     0.00003946    -0.00000853    -0.00002270    -0.00000176    -0.00003591    -0.00001789
   MO  77    -0.00015935     0.00065802    -0.00042149    -0.00000551    -0.00005919    -0.00001294    -0.00000284    -0.00000101
   MO  78     0.00005658    -0.00012857     0.00007078     0.00004891    -0.00023385    -0.00008792     0.00000539    -0.00002960
   MO  79    -0.00002242    -0.00001702     0.00000630     0.00003157    -0.00005403     0.00002377    -0.00010800     0.00000393
   MO  80     0.00002905    -0.00001672    -0.00000076    -0.00010257    -0.00046961     0.00043577     0.00003497    -0.00002746
   MO  81    -0.00001807     0.00002989     0.00003979     0.00010090    -0.00004994    -0.00001701     0.00014338    -0.00002364
   MO  82     0.00001480    -0.00001882    -0.00001958     0.00007402    -0.00035624    -0.00025066    -0.00003113     0.00001196

                MO  73         MO  74         MO  75         MO  76         MO  77         MO  78         MO  79         MO  80
   MO   3    -0.00278960    -0.00064280     0.00014720     0.00019191     0.00004101    -0.00020393     0.00043885    -0.00009409
   MO   4     0.00232008    -0.00519291     0.00668975     0.00016877    -0.00001614    -0.00043526     0.00011372    -0.00036091
   MO   5    -0.00433636    -0.00090042    -0.00014701    -0.00008685     0.00012024    -0.00012047     0.00271515    -0.00007693
   MO   6    -0.00011023    -0.00069284     0.00231066     0.00002197    -0.00006342    -0.00041980    -0.00012264    -0.00017181
   MO   7     0.00038499    -0.00121955    -0.00003803    -0.00016426    -0.00067038    -0.00358214    -0.00034280    -0.00284421
   MO   8    -0.00135565    -0.00005924     0.00006729    -0.00002576    -0.00004312    -0.00006277     0.00018945     0.00018763
   MO   9     0.00001298     0.00007397    -0.00007818     0.00209412     0.00127234    -0.00027848     0.00003664    -0.00000043
   MO  10    -0.00001175     0.00008912    -0.00007574    -0.00000168    -0.00000667    -0.00001873     0.00001303     0.00003993
   MO  11     0.00020288    -0.00045173     0.00039249    -0.00000531    -0.00000359    -0.00003716     0.00007303     0.00002224
   MO  12     0.00000160     0.00022272    -0.00019715    -0.00001826     0.00000283     0.00006535    -0.00008479    -0.00019977
   MO  13    -0.00009661    -0.00004856    -0.00000262    -0.00001555    -0.00000254     0.00000886     0.00020350    -0.00026079
   MO  14    -0.00005208    -0.00004371    -0.00013134    -0.00001712    -0.00000887    -0.00000424     0.00013462     0.00011306
   MO  15     0.00005661    -0.00007585     0.00000281     0.00001970     0.00003722     0.00000362     0.00002469    -0.00009709
   MO  16    -0.00000346    -0.00009023     0.00019374    -0.00001299     0.00001931     0.00008506    -0.00002997    -0.00020255
   MO  17     0.00041737    -0.00030607     0.00012637    -0.00003108    -0.00002030    -0.00012947     0.00016347    -0.00008935
   MO  18    -0.00010600     0.00039027    -0.00041715    -0.00002171     0.00002342     0.00012003    -0.00012515     0.00002351
   MO  19    -0.00015905     0.00002976    -0.00017719    -0.00001739    -0.00004054    -0.00003191     0.00012384     0.00003004
   MO  20     0.00010722     0.00027618    -0.00028243    -0.00008583    -0.00004200     0.00014187    -0.00005284    -0.00025823
   MO  21     0.00009515     0.00007417    -0.00014794    -0.00006281    -0.00005493    -0.00011750    -0.00016187     0.00043642
   MO  22    -0.00017760     0.00042366    -0.00047136     0.00006048     0.00003098    -0.00016863     0.00020073    -0.00000239
   MO  23    -0.00008732    -0.00031711    -0.00010971     0.00003295     0.00001867     0.00003185    -0.00039453    -0.00027063
   MO  24     0.00009000     0.00007847    -0.00082481    -0.00001766    -0.00005354    -0.00025099     0.00010162     0.00040691
   MO  25    -0.00022499    -0.00015421     0.00004317    -0.00020785    -0.00017546    -0.00005979     0.00043658    -0.00011103
   MO  26    -0.00046762    -0.00003181    -0.00055448    -0.00000198     0.00007554     0.00025500     0.00028876     0.00006757
   MO  27    -0.00016870     0.00075015    -0.00050197     0.00002238    -0.00005423    -0.00022764    -0.00028646     0.00006239
   MO  28     0.00000434    -0.00010884     0.00002621     0.00059647     0.00050209    -0.00008311     0.00007580     0.00002532
   MO  29    -0.00023982     0.00039734     0.00014049    -0.00015274     0.00000142    -0.00059417     0.00016713    -0.00019460
   MO  30     0.00004719    -0.00008294    -0.00000723    -0.00061521     0.00041716     0.00006028     0.00001777     0.00007863
   MO  31     0.00005592     0.00007051    -0.00017070    -0.00009904     0.00051409    -0.00001532    -0.00012804    -0.00001017
   MO  32     0.00001939     0.00016860    -0.00036775    -0.00013093    -0.00022581     0.00021114    -0.00029313     0.00021407
   MO  33    -0.00006426    -0.00008227     0.00029801    -0.00003610     0.00002126    -0.00016159    -0.00046250     0.00004469
   MO  34    -0.00009191     0.00003117     0.00009336     0.00007060     0.00006829    -0.00009220     0.00012143     0.00010111
   MO  35     0.00010430     0.00003813    -0.00013103    -0.00005305    -0.00001745     0.00009370    -0.00006903    -0.00012556
   MO  36    -0.00003544    -0.00001674     0.00000195    -0.00004648    -0.00001278    -0.00004316     0.00009186     0.00013470
   MO  37     0.00004275    -0.00003205    -0.00001128     0.00002519    -0.00000847    -0.00007808     0.00006885     0.00009185
   MO  38    -0.00030358     0.00038385    -0.00030607    -0.00001214     0.00005058     0.00024387    -0.00008106    -0.00009511
   MO  39    -0.00005291    -0.00013308     0.00000105    -0.00003382    -0.00006336    -0.00020408    -0.00013544     0.00006956
   MO  40     0.00013202     0.00043104    -0.00005082     0.00004599     0.00001687    -0.00006038    -0.00001181     0.00015147
   MO  41     0.00028575    -0.00021019     0.00029450     0.00001912     0.00004085     0.00012206    -0.00015213     0.00037216
   MO  42     0.00011037     0.00001179     0.00028356     0.00008381     0.00005789    -0.00000026    -0.00007950     0.00003193
   MO  43    -0.00032768    -0.00005924    -0.00002313     0.00003365    -0.00005649     0.00004967     0.00046750    -0.00002068
   MO  44     0.00010103     0.00029676    -0.00127471     0.00000121    -0.00002224    -0.00022351    -0.00006712    -0.00018213
   MO  45     0.00046114    -0.00090410     0.00078840    -0.00011730    -0.00011476     0.00035124    -0.00027325    -0.00058963
   MO  46     0.00003504    -0.00032194    -0.00016958     0.00079822     0.00055340     0.00033052     0.00001982    -0.00007901
   MO  47    -0.00012001     0.00023181    -0.00014511    -0.00040645    -0.00087172     0.00013897    -0.00011823     0.00006179
   MO  48    -0.00016791     0.00040150     0.00054470     0.00100322    -0.00046433    -0.00081735    -0.00012497    -0.00009096
   MO  49     0.00002924    -0.00030927    -0.00060456     0.00071711    -0.00051718     0.00092052     0.00005668     0.00008116
   MO  50    -0.00030039     0.00013189    -0.00025808    -0.00003519     0.00061915     0.00007836    -0.00011393    -0.00004784
   MO  51    -0.00013807    -0.00023680     0.00008618    -0.00000335    -0.00005605    -0.00016414     0.00139910     0.00025037
   MO  52    -0.00040146    -0.00009310    -0.00042783    -0.00010284    -0.00002190     0.00009409     0.00014322     0.00011224
   MO  53     0.00033162    -0.00018345     0.00021734     0.00003330     0.00007102     0.00008424    -0.00028861     0.00087031
   MO  54     0.00010419    -0.00001855    -0.00029309     0.00007977     0.00005605    -0.00006220    -0.00021109     0.00052850
   MO  55     0.00001728    -0.00000045     0.00000980    -0.00024585    -0.00038503     0.00007524     0.00001544     0.00002618
   MO  56     0.00000721     0.00079094    -0.00051215     0.00001129    -0.00003273    -0.00029526    -0.00008341     0.00003317
   MO  57     0.00039809    -0.00007672     0.00029728     0.00007095     0.00005281    -0.00013647    -0.00008729    -0.00022159
   MO  58    -0.00000480     0.00023651    -0.00001945    -0.00001242    -0.00003756    -0.00016520    -0.00004833    -0.00005387
   MO  59     0.00000522    -0.00007904    -0.00004331     0.00000579     0.00001231     0.00011438     0.00003181     0.00004588
   MO  60     0.00018686     0.00004430    -0.00001514    -0.00001677    -0.00004335     0.00001504     0.00003938    -0.00001253
   MO  61     0.00007498    -0.00001839    -0.00013932     0.00003401     0.00005477     0.00010430    -0.00002146     0.00003193
   MO  62    -0.00003595    -0.00007401    -0.00010594     0.00002857     0.00000386     0.00021079    -0.00050143     0.00010884
   MO  63    -0.00008868     0.00001974     0.00044084    -0.00001374    -0.00011226    -0.00036902    -0.00017279    -0.00013420
   MO  64    -0.00002795    -0.00001537    -0.00022838    -0.00000662     0.00001265     0.00006093    -0.00001884    -0.00010113
   MO  65     0.00002739     0.00000940     0.00001213    -0.00062959    -0.00015935     0.00005658    -0.00002242     0.00002905
   MO  66    -0.00000630    -0.00000732     0.00000351    -0.00025765     0.00065802    -0.00012857    -0.00001702    -0.00001672
   MO  67     0.00000051    -0.00000195     0.00000111     0.00003946    -0.00042149     0.00007078     0.00000630    -0.00000076
   MO  68     0.00003521    -0.00000790    -0.00004712    -0.00000853    -0.00000551     0.00004891     0.00003157    -0.00010257
   MO  69     0.00003645    -0.00008136    -0.00047973    -0.00002270    -0.00005919    -0.00023385    -0.00005403    -0.00046961
   MO  70    -0.00026900     0.00025672    -0.00010757    -0.00000176    -0.00001294    -0.00008792     0.00002377     0.00043577
   MO  71    -0.00004542     0.00003233    -0.00002168    -0.00003591    -0.00000284     0.00000539    -0.00010800     0.00003497
   MO  72     0.00002868     0.00004803    -0.00001715    -0.00001789    -0.00000101    -0.00002960     0.00000393    -0.00002746
   MO  73     0.00054310    -0.00015564     0.00007572    -0.00000826     0.00002802     0.00011521    -0.00008011    -0.00004942
   MO  74    -0.00015564     0.00090248    -0.00010887    -0.00001415    -0.00009472    -0.00051817    -0.00005730     0.00003777
   MO  75     0.00007572    -0.00010887     0.00119642    -0.00002329    -0.00008639    -0.00042596    -0.00002056     0.00006493
   MO  76    -0.00000826    -0.00001415    -0.00002329     0.00159725    -0.00001368    -0.00001149    -0.00001760     0.00000431
   MO  77     0.00002802    -0.00009472    -0.00008639    -0.00001368     0.00143374    -0.00003986    -0.00002784     0.00002196
   MO  78     0.00011521    -0.00051817    -0.00042596    -0.00001149    -0.00003986     0.00123348    -0.00000482     0.00010759
   MO  79    -0.00008011    -0.00005730    -0.00002056    -0.00001760    -0.00002784    -0.00000482     0.00075030     0.00003446
   MO  80    -0.00004942     0.00003777     0.00006493     0.00000431     0.00002196     0.00010759     0.00003446     0.00095383
   MO  81    -0.00007224    -0.00002837     0.00000607     0.00000163     0.00000499     0.00000878     0.00004319     0.00001891
   MO  82     0.00005232    -0.00007242     0.00014273     0.00000489     0.00000786     0.00002943    -0.00000159     0.00003710

                MO  81         MO  82
   MO   3     0.00024665     0.00003723
   MO   4     0.00006902     0.00122112
   MO   5     0.00045533    -0.00010028
   MO   6     0.00011368     0.00149375
   MO   7    -0.00029439    -0.00184186
   MO   8    -0.00213061     0.00014611
   MO   9    -0.00003612    -0.00000649
   MO  10    -0.00000062    -0.00005051
   MO  11    -0.00027401    -0.00002440
   MO  12    -0.00026277     0.00014798
   MO  13    -0.00003443    -0.00005796
   MO  14     0.00004400     0.00034336
   MO  15    -0.00001036     0.00012192
   MO  16     0.00000659     0.00005945
   MO  17    -0.00068941     0.00004478
   MO  18    -0.00004442     0.00002662
   MO  19     0.00035072    -0.00006635
   MO  20    -0.00051273     0.00014959
   MO  21    -0.00001988    -0.00023111
   MO  22     0.00009122    -0.00027960
   MO  23     0.00050462     0.00005873
   MO  24     0.00007956     0.00064850
   MO  25     0.00020324    -0.00005708
   MO  26     0.00008025     0.00026564
   MO  27    -0.00017694     0.00025486
   MO  28    -0.00010809    -0.00003617
   MO  29     0.00015727     0.00016083
   MO  30    -0.00011595    -0.00003021
   MO  31     0.00004354     0.00000174
   MO  32    -0.00023519    -0.00006334
   MO  33     0.00013226    -0.00003195
   MO  34    -0.00003633    -0.00006559
   MO  35    -0.00014187     0.00006500
   MO  36    -0.00006324    -0.00010943
   MO  37     0.00001875    -0.00009458
   MO  38     0.00056671     0.00025122
   MO  39     0.00032469    -0.00005983
   MO  40    -0.00040496     0.00015549
   MO  41    -0.00007292     0.00032135
   MO  42    -0.00010618     0.00006963
   MO  43     0.00002957     0.00022368
   MO  44    -0.00011149     0.00025619
   MO  45     0.00002914     0.00018824
   MO  46     0.00004903    -0.00002891
   MO  47     0.00001256    -0.00004745
   MO  48     0.00015327    -0.00000833
   MO  49    -0.00005127    -0.00002168
   MO  50     0.00010244    -0.00005340
   MO  51     0.00006848    -0.00008076
   MO  52     0.00075923     0.00053120
   MO  53     0.00014932    -0.00064354
   MO  54    -0.00025618     0.00086391
   MO  55    -0.00004081     0.00002663
   MO  56     0.00000216     0.00017833
   MO  57    -0.00002903     0.00012168
   MO  58    -0.00000650     0.00012774
   MO  59     0.00000261    -0.00004315
   MO  60     0.00017979     0.00000245
   MO  61     0.00031640    -0.00006107
   MO  62    -0.00004511     0.00001802
   MO  63     0.00010718    -0.00003962
   MO  64    -0.00008798    -0.00022605
   MO  65    -0.00001807     0.00001480
   MO  66     0.00002989    -0.00001882
   MO  67     0.00003979    -0.00001958
   MO  68     0.00010090     0.00007402
   MO  69    -0.00004994    -0.00035624
   MO  70    -0.00001701    -0.00025066
   MO  71     0.00014338    -0.00003113
   MO  72    -0.00002364     0.00001196
   MO  73    -0.00007224     0.00005232
   MO  74    -0.00002837    -0.00007242
   MO  75     0.00000607     0.00014273
   MO  76     0.00000163     0.00000489
   MO  77     0.00000499     0.00000786
   MO  78     0.00000878     0.00002943
   MO  79     0.00004319    -0.00000159
   MO  80     0.00001891     0.00003710
   MO  81     0.00051174     0.00003316
   MO  82     0.00003316     0.00091305

Natural orbital populations,block 1
              MO     1       MO     2       MO     3       MO     4       MO     5       MO     6       MO     7       MO     8
  occ(*)=     2.00000000     2.00000000     1.97898219     1.97375579     1.97163433     1.96669788     1.95967991     1.01162155
              MO     9       MO    10       MO    11       MO    12       MO    13       MO    14       MO    15       MO    16
  occ(*)=     0.97442999     0.02952688     0.01938770     0.01350580     0.01308825     0.01078292     0.01019407     0.00816017
              MO    17       MO    18       MO    19       MO    20       MO    21       MO    22       MO    23       MO    24
  occ(*)=     0.00639341     0.00594981     0.00488406     0.00483653     0.00449545     0.00425517     0.00422558     0.00408081
              MO    25       MO    26       MO    27       MO    28       MO    29       MO    30       MO    31       MO    32
  occ(*)=     0.00317742     0.00210219     0.00174500     0.00154565     0.00130483     0.00117759     0.00087553     0.00076968
              MO    33       MO    34       MO    35       MO    36       MO    37       MO    38       MO    39       MO    40
  occ(*)=     0.00067469     0.00061990     0.00060904     0.00052316     0.00042103     0.00038189     0.00036028     0.00035255
              MO    41       MO    42       MO    43       MO    44       MO    45       MO    46       MO    47       MO    48
  occ(*)=     0.00033311     0.00029480     0.00027372     0.00025513     0.00023684     0.00020399     0.00016157     0.00014371
              MO    49       MO    50       MO    51       MO    52       MO    53       MO    54       MO    55       MO    56
  occ(*)=     0.00013378     0.00012522     0.00008743     0.00007718     0.00005988     0.00005443     0.00003989     0.00003706
              MO    57       MO    58       MO    59       MO    60       MO    61       MO    62       MO    63       MO    64
  occ(*)=     0.00003366     0.00003099     0.00002989     0.00002422     0.00002075     0.00001935     0.00001789     0.00001651
              MO    65       MO    66       MO    67       MO    68       MO    69       MO    70       MO    71       MO    72
  occ(*)=     0.00001325     0.00001153     0.00001102     0.00000926     0.00000743     0.00000665     0.00000505     0.00000485
              MO    73       MO    74       MO    75       MO    76       MO    77       MO    78       MO    79       MO    80
  occ(*)=     0.00000419     0.00000306     0.00000215     0.00000072     0.00000042     0.00000027     0.00000021     0.00000013
              MO    81       MO    82       MO
  occ(*)=     0.00000007     0.00000002


 total number of electrons =   16.0000000000

 test slabel:                     1
  a  


          Mulliken population analysis


  NOTE: For HERMIT use spherical harmonics basis sets !!!
 

                        a   partial gross atomic populations
   ao class       1a         2a         3a         4a         5a         6a  
     2_ s       1.993633   0.000003  -0.000522   0.099897   1.256339   0.007233
     2_ p      -0.000000   0.000004  -0.000390   0.736718   0.113736   0.000708
     2_ d       0.000001   0.000000   0.000055  -0.001837  -0.003445   0.000382
     4_ s      -0.000001   1.998050   1.477155   0.000180   0.000774   0.242967
     4_ p      -0.000006   0.000060   0.013694   0.000798   0.009606   1.207609
     4_ d      -0.000000   0.000001  -0.003011  -0.000044   0.000045   0.006800
     5_ s       0.002536   0.000001  -0.000281   0.900124   0.019647   0.000139
     5_ p       0.000781  -0.000000  -0.000135   0.013939  -0.026507  -0.000017
     6_ s      -0.000002   0.000417   0.080102   0.000089   0.001050   0.002525
     6_ p      -0.000000   0.000272   0.008405   0.000004   0.000181   0.032089
     7_ s       0.000000   0.001029   0.365968   0.000374   0.003870   0.460837
     7_ p       0.000000   0.000175   0.032027  -0.000007   0.000002   0.003443
     8_ s       0.002268  -0.000012   0.005410   0.266384   0.569425   0.001363
     8_ p       0.000789   0.000000   0.000507  -0.042863   0.026911   0.000619
 
   ao class       7a         8a         9a        10a        11a        12a  
     2_ s       0.001124   0.035915   0.034122  -0.000000   0.000013   0.003639
     2_ p       0.000129   0.477118   0.432652  -0.000001  -0.000017   0.001579
     2_ d       0.000181  -0.000296  -0.000044  -0.000001   0.000002   0.000418
     4_ s       0.015321   0.000185   0.000063   0.000043   0.000077   0.000171
     4_ p       1.248135   0.475873   0.485697   0.013116   0.009557   0.000056
     4_ d       0.010319   0.001186   0.001469   0.001176   0.001086   0.000022
     5_ s      -0.000028   0.003473   0.003652  -0.000000  -0.000003   0.003031
     5_ p      -0.000050  -0.000124  -0.001501   0.000000   0.000002  -0.000245
     6_ s       0.630814   0.000136   0.000064   0.014666   0.000721   0.000008
     6_ p       0.009991   0.006463   0.006918  -0.000752   0.000503   0.000003
     7_ s       0.003802  -0.000082  -0.000019   0.001261   0.007253   0.000028
     7_ p       0.040027   0.007398   0.007951   0.000016   0.000165   0.000002
     8_ s      -0.000603   0.005093   0.004068   0.000003   0.000028   0.005190
     8_ p       0.000518  -0.000717  -0.000660  -0.000000   0.000002  -0.000396
 
   ao class      13a        14a        15a        16a        17a        18a  
     2_ s       0.000151   0.000153   0.000047  -0.000000   0.000002   0.000005
     2_ p       0.011565   0.006994   0.000677   0.000004  -0.000007   0.000029
     2_ d       0.000027   0.000347   0.000155   0.000000   0.000000   0.000010
     4_ s       0.000023   0.000656   0.005274   0.000000   0.000078  -0.000023
     4_ p       0.000008   0.000273   0.002181   0.005257   0.003917   0.000275
     4_ d       0.000004   0.000052   0.000387   0.001637   0.002001   0.004139
     5_ s       0.000134   0.002180   0.000006  -0.000000  -0.000000   0.000002
     5_ p       0.000521  -0.000054   0.000016   0.000001   0.000001   0.000002
     6_ s       0.000000   0.000002  -0.000068   0.000000  -0.000024   0.000040
     6_ p       0.000000   0.000031   0.000181   0.000668   0.000455   0.000560
     7_ s       0.000006   0.000064   0.000507   0.000000  -0.000333   0.000086
     7_ p      -0.000000   0.000036   0.000259   0.000589   0.000291   0.000833
     8_ s       0.000201  -0.000115   0.000512  -0.000000   0.000005  -0.000017
     8_ p       0.000448   0.000162   0.000061   0.000004   0.000008   0.000009
 
   ao class      19a        20a        21a        22a        23a        24a  
     2_ s       0.000001  -0.000000   0.000000   0.000026   0.000011   0.001096
     2_ p       0.000007   0.000000   0.000008   0.000594   0.001044   0.001788
     2_ d       0.000022   0.000000   0.003187   0.000915   0.001719   0.001349
     4_ s       0.000742   0.000000  -0.000000   0.000001  -0.000004   0.000001
     4_ p       0.000202   0.000016   0.000004   0.000766   0.000392   0.000088
     4_ d       0.003555   0.004193   0.000014   0.001712   0.000833   0.000168
     5_ s      -0.000000   0.000000  -0.000001  -0.000266  -0.000173  -0.000394
     5_ p       0.000004  -0.000000   0.000669   0.000151   0.000477   0.000100
     6_ s       0.000233  -0.000000   0.000000   0.000001   0.000000  -0.000000
     6_ p       0.000063   0.000227   0.000001   0.000106   0.000050   0.000010
     7_ s       0.000013   0.000000   0.000000   0.000000   0.000002   0.000000
     7_ p       0.000029   0.000400   0.000001   0.000087   0.000040   0.000009
     8_ s       0.000002   0.000000  -0.000003  -0.000047  -0.000410  -0.000310
     8_ p       0.000011  -0.000000   0.000615   0.000210   0.000244   0.000175
 
   ao class      25a        26a        27a        28a        29a        30a  
     2_ s       0.001907   0.000000   0.000469   0.000000  -0.000000  -0.000000
     2_ p       0.001458   0.000024   0.000238   0.000003   0.002525  -0.000000
     2_ d       0.000418   0.001928   0.000897   0.000001   0.000135  -0.000000
     4_ s       0.000013   0.000000  -0.000000   0.000004   0.000000   0.000115
     4_ p       0.000036   0.000000  -0.000001   0.000548  -0.000006   0.000010
     4_ d       0.000018   0.000000   0.000001   0.000761   0.000005   0.000353
     5_ s      -0.000333  -0.000000  -0.000006  -0.000001  -0.001175   0.000001
     5_ p      -0.000118   0.000075   0.000084   0.000000   0.000471  -0.000000
     6_ s       0.000010   0.000000  -0.000000   0.000107  -0.000001   0.000181
     6_ p      -0.000000  -0.000000   0.000000  -0.000062  -0.000000   0.000258
     7_ s      -0.000006  -0.000000   0.000001   0.000157  -0.000000   0.000022
     7_ p       0.000003   0.000000   0.000000   0.000030  -0.000001   0.000240
     8_ s      -0.000117   0.000000  -0.000028  -0.000004  -0.001068  -0.000000
     8_ p      -0.000112   0.000074   0.000092   0.000001   0.000419  -0.000001
 
   ao class      31a        32a        33a        34a        35a        36a  
     2_ s      -0.000000   0.000002   0.000211   0.000000   0.000000   0.000000
     2_ p       0.000000  -0.000001   0.000035   0.000000  -0.000001   0.000000
     2_ d       0.000000   0.000003   0.000332   0.000000   0.000001   0.000000
     4_ s      -0.000031   0.000022   0.000001   0.000000   0.000169  -0.000000
     4_ p       0.000161   0.000057   0.000003   0.000003   0.000294   0.000024
     4_ d       0.000117   0.000223   0.000000   0.000152   0.000060   0.000145
     5_ s      -0.000000  -0.000003  -0.000125   0.000000   0.000002  -0.000000
     5_ p       0.000000   0.000001   0.000165   0.000000   0.000000  -0.000000
     6_ s       0.000598  -0.000069   0.000001   0.000000   0.000054  -0.000000
     6_ p       0.000061   0.000085   0.000000   0.000294   0.000006   0.000101
     7_ s       0.000012   0.000233  -0.000001  -0.000000  -0.000029   0.000000
     7_ p      -0.000041   0.000218  -0.000000   0.000171   0.000058   0.000253
     8_ s      -0.000000   0.000001  -0.000097  -0.000000  -0.000004   0.000000
     8_ p      -0.000002  -0.000001   0.000149  -0.000000  -0.000000  -0.000000
 
   ao class      37a        38a        39a        40a        41a        42a  
     2_ s       0.000000   0.000012   0.000163   0.000089   0.000000   0.000007
     2_ p       0.000034  -0.000004   0.000027  -0.000018  -0.000000  -0.000003
     2_ d       0.000034   0.000003   0.000087   0.000072   0.000000   0.000003
     4_ s       0.000001  -0.000001   0.000008   0.000010   0.000000   0.000039
     4_ p       0.000000  -0.000045   0.000002   0.000005   0.000097   0.000011
     4_ d       0.000001   0.000141   0.000005   0.000016   0.000139   0.000180
     5_ s       0.000000   0.000003  -0.000005   0.000012   0.000000  -0.000007
     5_ p       0.000178   0.000001   0.000081   0.000033   0.000000   0.000008
     6_ s       0.000000   0.000019  -0.000003   0.000004   0.000000   0.000010
     6_ p       0.000000   0.000056   0.000001   0.000003   0.000042   0.000041
     7_ s      -0.000000   0.000089   0.000004  -0.000000   0.000000  -0.000001
     7_ p       0.000000   0.000095   0.000000   0.000010   0.000055   0.000007
     8_ s       0.000000   0.000012  -0.000041   0.000024  -0.000000  -0.000005
     8_ p       0.000173   0.000001   0.000030   0.000094   0.000000   0.000004
 
   ao class      43a        44a        45a        46a        47a        48a  
     2_ s       0.000004   0.000013   0.000000  -0.000000   0.000000   0.000001
     2_ p       0.000148   0.000006   0.000000   0.000001   0.000000  -0.000002
     2_ d       0.000029   0.000010   0.000000  -0.000000   0.000032   0.000000
     4_ s       0.000001   0.000051   0.000000  -0.000002  -0.000000  -0.000023
     4_ p      -0.000001   0.000028   0.000005  -0.000009   0.000000  -0.000001
     4_ d       0.000008   0.000107   0.000210   0.000075   0.000000   0.000122
     5_ s       0.000054  -0.000043   0.000000   0.000001   0.000000  -0.000001
     5_ p      -0.000010   0.000021  -0.000000  -0.000000   0.000062   0.000000
     6_ s      -0.000000  -0.000008  -0.000000   0.000009   0.000000   0.000022
     6_ p       0.000002   0.000007   0.000022   0.000058   0.000000   0.000006
     7_ s       0.000001   0.000010  -0.000000   0.000032   0.000000   0.000012
     7_ p      -0.000000   0.000009   0.000001   0.000038   0.000000   0.000008
     8_ s       0.000053   0.000050   0.000000   0.000000  -0.000000  -0.000001
     8_ p      -0.000014  -0.000005  -0.000000   0.000002   0.000067   0.000000
 
   ao class      49a        50a        51a        52a        53a        54a  
     2_ s       0.000000   0.000004   0.000000  -0.000000  -0.000005  -0.000006
     2_ p       0.000001   0.000004  -0.000000  -0.000008   0.000006  -0.000004
     2_ d       0.000003   0.000086   0.000000   0.000100   0.000006   0.000002
     4_ s       0.000000   0.000000  -0.000001  -0.000000  -0.000016  -0.000027
     4_ p       0.000031   0.000000   0.000003   0.000001  -0.000001  -0.000008
     4_ d       0.000069   0.000000   0.000013   0.000001  -0.000000   0.000000
     5_ s       0.000001  -0.000020  -0.000000  -0.000002   0.000000  -0.000001
     5_ p      -0.000000   0.000034  -0.000000  -0.000007   0.000013   0.000004
     6_ s       0.000000  -0.000000   0.000037   0.000000   0.000004   0.000023
     6_ p       0.000015   0.000000   0.000015   0.000001   0.000004   0.000003
     7_ s       0.000000   0.000000   0.000011  -0.000000   0.000026   0.000052
     7_ p       0.000012   0.000001   0.000009   0.000001   0.000001  -0.000001
     8_ s      -0.000002  -0.000022  -0.000000  -0.000013   0.000011   0.000015
     8_ p       0.000004   0.000038  -0.000000   0.000004   0.000011   0.000001
 
   ao class      55a        56a        57a        58a        59a        60a  
     2_ s      -0.000023  -0.000000  -0.000000  -0.000000  -0.000000  -0.000003
     2_ p       0.000001   0.000033  -0.000039  -0.000019   0.000003  -0.000001
     2_ d       0.000039   0.000007   0.000001   0.000004   0.000044   0.000000
     4_ s       0.000000  -0.000000  -0.000000  -0.000000  -0.000000  -0.000000
     4_ p       0.000000   0.000000  -0.000003  -0.000004  -0.000000   0.000001
     4_ d       0.000000   0.000000  -0.000000   0.000000  -0.000000   0.000000
     5_ s       0.000020   0.000000   0.000035   0.000028   0.000001   0.000010
     5_ p       0.000001   0.000004  -0.000012  -0.000010  -0.000018  -0.000003
     6_ s      -0.000000  -0.000000  -0.000011  -0.000015  -0.000001   0.000007
     6_ p       0.000000   0.000000   0.000025   0.000026   0.000001  -0.000006
     7_ s       0.000000   0.000000  -0.000007  -0.000005   0.000000  -0.000004
     7_ p       0.000000   0.000000   0.000010   0.000017   0.000002   0.000024
     8_ s       0.000001   0.000000   0.000046   0.000013   0.000004  -0.000000
     8_ p       0.000001  -0.000008  -0.000012  -0.000004  -0.000005   0.000000
 
   ao class      61a        62a        63a        64a        65a        66a  
     2_ s      -0.000024  -0.000001  -0.000000  -0.000000   0.000001   0.000000
     2_ p      -0.000001   0.000000  -0.000000  -0.000000   0.000001  -0.000002
     2_ d       0.000001   0.000017   0.000003   0.000000   0.000007   0.000002
     4_ s      -0.000000  -0.000000   0.000000   0.000000  -0.000000  -0.000000
     4_ p      -0.000000  -0.000000   0.000001   0.000000   0.000001   0.000005
     4_ d       0.000000  -0.000000   0.000000  -0.000000  -0.000000  -0.000000
     5_ s       0.000030   0.000001   0.000000  -0.000000   0.000003   0.000001
     5_ p      -0.000005  -0.000000  -0.000000  -0.000000  -0.000002   0.000000
     6_ s      -0.000000   0.000000   0.000001   0.000000  -0.000000   0.000000
     6_ p       0.000002   0.000001   0.000011   0.000011   0.000000   0.000001
     7_ s      -0.000000  -0.000000   0.000005   0.000000  -0.000000  -0.000000
     7_ p       0.000002   0.000001  -0.000003   0.000006   0.000002   0.000002
     8_ s       0.000016   0.000002   0.000001   0.000000   0.000001  -0.000000
     8_ p       0.000000  -0.000000  -0.000000   0.000000  -0.000001   0.000002
 
   ao class      67a        68a        69a        70a        71a        72a  
     2_ s      -0.000000  -0.000000  -0.000000   0.000000   0.000000  -0.000000
     2_ p      -0.000000  -0.000001  -0.000000  -0.000000  -0.000000  -0.000001
     2_ d      -0.000000  -0.000003  -0.000004  -0.000000  -0.000000  -0.000000
     4_ s      -0.000000  -0.000000   0.000000   0.000000   0.000000  -0.000000
     4_ p       0.000005   0.000001  -0.000000   0.000007   0.000001   0.000000
     4_ d      -0.000001  -0.000000  -0.000000  -0.000000  -0.000000  -0.000000
     5_ s      -0.000000  -0.000000   0.000000   0.000000   0.000000   0.000000
     5_ p       0.000000   0.000003   0.000003   0.000001   0.000004   0.000001
     6_ s      -0.000000   0.000000   0.000000  -0.000003  -0.000001  -0.000000
     6_ p       0.000002   0.000000  -0.000000   0.000001   0.000000   0.000000
     7_ s      -0.000000  -0.000000   0.000000   0.000001   0.000000   0.000000
     7_ p       0.000004   0.000001  -0.000000  -0.000001  -0.000000  -0.000000
     8_ s       0.000001   0.000004   0.000000  -0.000002  -0.000000  -0.000003
     8_ p       0.000000   0.000004   0.000009   0.000002   0.000001   0.000008
 
   ao class      73a        74a        75a        76a        77a        78a  
     2_ s       0.000000  -0.000001   0.000000   0.000000   0.000000   0.000001
     2_ p       0.000001  -0.000001   0.000000   0.000001   0.000000   0.000000
     2_ d      -0.000000  -0.000001  -0.000000  -0.000000  -0.000000  -0.000000
     4_ s       0.000000   0.000001   0.000004  -0.000000  -0.000000   0.000000
     4_ p       0.000000   0.000001   0.000001  -0.000000  -0.000000  -0.000000
     4_ d      -0.000000  -0.000000  -0.000000  -0.000000  -0.000000  -0.000000
     5_ s      -0.000000  -0.000001  -0.000000   0.000000  -0.000000  -0.000001
     5_ p       0.000003   0.000005   0.000001   0.000000   0.000000   0.000000
     6_ s      -0.000000  -0.000000  -0.000001  -0.000000   0.000000   0.000000
     6_ p       0.000000  -0.000000   0.000000  -0.000000  -0.000000  -0.000000
     7_ s      -0.000000  -0.000001  -0.000003   0.000000   0.000000   0.000000
     7_ p       0.000000   0.000000   0.000001  -0.000000  -0.000000  -0.000000
     8_ s       0.000000   0.000000  -0.000000  -0.000001  -0.000000   0.000000
     8_ p       0.000000   0.000000   0.000000   0.000000  -0.000000   0.000000
 
   ao class      79a        80a        81a        82a  
     2_ s       0.000000  -0.000000  -0.000000  -0.000000
     2_ p       0.000000   0.000000   0.000000  -0.000000
     2_ d      -0.000000  -0.000000  -0.000000  -0.000000
     4_ s      -0.000000   0.000000  -0.000000  -0.000000
     4_ p      -0.000000  -0.000000  -0.000000  -0.000000
     4_ d      -0.000000  -0.000000  -0.000000  -0.000000
     5_ s      -0.000000  -0.000000  -0.000000   0.000000
     5_ p      -0.000000   0.000000  -0.000000  -0.000000
     6_ s       0.000000   0.000000   0.000000  -0.000000
     6_ p      -0.000000  -0.000000  -0.000000  -0.000000
     7_ s      -0.000000   0.000000   0.000000   0.000000
     7_ p       0.000000  -0.000000  -0.000000  -0.000000
     8_ s       0.000000   0.000000  -0.000000  -0.000000
     8_ p      -0.000000  -0.000000   0.000000  -0.000000


                        gross atomic populations
     ao            2_         4_         5_         6_         7_         8_
      s         3.435707   3.742074   0.932256   0.731744   0.845278   0.857276
      p         1.789382   3.478838  -0.010897   0.067522   0.095021  -0.012316
      d         0.007450   0.040666   0.000000   0.000000   0.000000   0.000000
    total       5.232539   7.261577   0.921359   0.799266   0.940299   0.844960
 

 Total number of electrons:   16.00000000

 item #                     2 suffix=:.drt1.state2:
 read_civout: repnuc=  -62.1137172157958     
================================================================================
  Reading record                      1  of civout
 INFO:ref#  1vector#  1 method:  0 last record  0max overlap with ref# 95% root-following 0
 MR-CISD energy:   -94.36141904   -32.24770183
 residuum:     0.00004500
 deltae:     0.00000000
================================================================================
  Reading record                      2  of civout
 INFO:ref#  2vector#  2 method:  0 last record  0max overlap with ref# 95% root-following 0
 MR-CISD energy:   -94.21846435   -32.10474713
 residuum:     0.00005590
 deltae:     0.00000000

          sovref  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1     0.94706787     0.00025589    -0.00126670    -0.00178264    -0.00337622     0.28275207     0.09678757    -0.06428875
 ref:   2    -0.00022678     0.94640844     0.02607982    -0.02063048     0.00040584     0.01457143     0.00185925     0.29125173
 ref:   3    -0.00089498     0.02479871    -0.93609036    -0.03112810     0.00196800    -0.00286088     0.00760621    -0.00568138

                ci   9         ci  10         ci  11         ci  12         ci  13         ci  14         ci  15         ci  16

                ci  17         ci  18         ci  19         ci  20         ci  21         ci  22         ci  23         ci  24

                ci  25         ci  26         ci  27         ci  28         ci  29         ci  30         ci  31         ci  32

                ci  33         ci  34         ci  35         ci  36         ci  37         ci  38         ci  39         ci  40

          tciref  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1     0.94706787     0.00025589    -0.00126670    -0.00178264    -0.00337622     0.00000000     0.00000000     0.00000000
 ref:   2    -0.00022678     0.94640844     0.02607982    -0.02063048     0.00040584     0.00000000     0.00000000     0.00000000
 ref:   3    -0.00089498     0.02479871    -0.93609036    -0.03112810     0.00196800     0.00000000     0.00000000     0.00000000

                ci   9         ci  10         ci  11         ci  12         ci  13         ci  14         ci  15         ci  16

                ci  17         ci  18         ci  19         ci  20         ci  21         ci  22         ci  23         ci  24

                ci  25         ci  26         ci  27         ci  28         ci  29         ci  30         ci  31         ci  32

                ci  33         ci  34         ci  35         ci  36         ci  37         ci  38         ci  39         ci  40
 computing final density
 =========== Executing IN-CORE method ==========
--------------------------------------------------------------------------------
  1e-density for root #    2
--------------------------------------------------------------------------------
================================================================================
   DYZ=     242  DYX=     503  DYW=     590
   D0Z=     147  D0Y=     928  D0X=     875  D0W=    1184
  DDZI=     189 DDYI=     692 DDXI=     675 DDWI=     866
  DDZE=       0 DDYE=     104 DDXE=     105 DDWE=     140
================================================================================
--------------------------------------------------------------------------------
  2e-density  for root #    2
--------------------------------------------------------------------------------
 call to prpd23 iwx=1
 call to prpd23 iwx=2
 starting prpd24 .... 
d2temp:fourex offset=     0 N= 59550 init,final=   T   Fmode=   1
d2temp:fourex offset= 59550 N= 59292 init,final=   F   Fmode=   1
d2temp:fourex offset=118842 N= 59913 init,final=   F   Fmode=   1
d2temp:fourex offset=178755 N= 59058 init,final=   F   Fmode=   1
d2temp:fourex offset=237813 N= 59667 init,final=   F   Fmode=   1
d2temp:fourex offset=297480 N= 59730 init,final=   F   Fmode=   1
d2temp:fourex offset=357210 N= 59871 init,final=   F   Fmode=   1
d2temp:fourex offset=417081 N= 59544 init,final=   F   Fmode=   1
d2temp:fourex offset=476625 N= 59835 init,final=   F   Fmode=   1
d2temp:fourex offset=536460 N= 59460 init,final=   F   Fmode=   1
d2temp:fourex offset=595920 N= 59343 init,final=   F   Fmode=   1
d2temp:fourex offset=655263 N= 58947 init,final=   F   Fmode=   1
d2temp:fourex offset=714210 N= 58800 init,final=   F   Fmode=   1
d2temp:fourex offset=773010 N= 58800 init,final=   F   Fmode=   1
d2temp:fourex offset=831810 N= 58845 init,final=   F   Fmode=   1
d2temp:fourex offset=890655 N= 59778 init,final=   F   Fmode=   1
d2temp:fourex offset=950433 N= 57120 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58182 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59445 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59955 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58695 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58620 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58140 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56676 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59124 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56793 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59799 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 55743 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59250 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58128 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57915 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 55545 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58545 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 55800 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58632 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57285 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59235 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59400 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58692 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57828 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57120 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59007 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56460 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56628 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59040 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57795 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59907 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58548 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 54660 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59535 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57525 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58692 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59958 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56793 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 54714 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 60000 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56403 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56658 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59787 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59556 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56082 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56250 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59892 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58140 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 55935 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58860 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N=     0 init,final=   F   Tmode=   1
 starting prpd24 .... 
d2temp:fourex offset=     0 N= 59550 init,final=   T   Fmode=   1
d2temp:fourex offset= 59550 N= 59292 init,final=   F   Fmode=   1
d2temp:fourex offset=118842 N= 59913 init,final=   F   Fmode=   1
d2temp:fourex offset=178755 N= 59058 init,final=   F   Fmode=   1
d2temp:fourex offset=237813 N= 59667 init,final=   F   Fmode=   1
d2temp:fourex offset=297480 N= 59730 init,final=   F   Fmode=   1
d2temp:fourex offset=357210 N= 59871 init,final=   F   Fmode=   1
d2temp:fourex offset=417081 N= 59544 init,final=   F   Fmode=   1
d2temp:fourex offset=476625 N= 59835 init,final=   F   Fmode=   1
d2temp:fourex offset=536460 N= 59460 init,final=   F   Fmode=   1
d2temp:fourex offset=595920 N= 59343 init,final=   F   Fmode=   1
d2temp:fourex offset=655263 N= 58947 init,final=   F   Fmode=   1
d2temp:fourex offset=714210 N= 58800 init,final=   F   Fmode=   1
d2temp:fourex offset=773010 N= 58800 init,final=   F   Fmode=   1
d2temp:fourex offset=831810 N= 58845 init,final=   F   Fmode=   1
d2temp:fourex offset=890655 N= 59778 init,final=   F   Fmode=   1
d2temp:fourex offset=950433 N= 57120 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58182 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59445 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59955 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58695 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58620 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58140 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56676 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59124 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56793 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59799 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 55743 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59250 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58128 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57915 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 55545 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58545 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 55800 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58632 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57285 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59235 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59400 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58692 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57828 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57120 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59007 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56460 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56628 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59040 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57795 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59907 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58548 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 54660 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59535 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57525 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58692 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59958 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56793 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 54714 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 60000 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56403 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56658 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59787 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59556 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56082 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56250 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59892 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58140 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 55935 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58860 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N=     0 init,final=   F   Tmode=   1
================================================================================
   DYZ=     242  DYX=    7351  DYW=    8696
   D0Z=     147  D0Y=    8560  D0X=    7763  D0W=   10896
  DDZI=     735 DDYI=    2660 DDXI=    2520 DDWI=    3140
  DDZE=       0 DDYE=       0 DDXE=       0 DDWE=       0
================================================================================
  entering dump_fourext .... 
d2temp:fourex offset=     1 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset= 60001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=120001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=180001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=240001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=300001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=360001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=420001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=480001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=540001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=600001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=660001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=720001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=780001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=840001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=900001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=960001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
Trace of MO density:    12.000000
   12  correlated and     4  frozen core electrons

           D1 (irrep)  block   1

                MO   1         MO   2         MO   3         MO   4         MO   5         MO   6         MO   7         MO   8
   MO   1     2.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
   MO   2     0.00000000     2.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
   MO   3     0.00000000     0.00000000     1.97126167     0.00024484     0.00084047     0.00043977    -0.00023507     0.00151514
   MO   4     0.00000000     0.00000000     0.00024484     1.96759580    -0.00014289    -0.00736609    -0.01390207     0.00100388
   MO   5     0.00000000     0.00000000     0.00084047    -0.00014289     1.97279100    -0.00038058     0.00020287     0.00032201
   MO   6     0.00000000     0.00000000     0.00043977    -0.00736609    -0.00038058     1.96777788    -0.00387871     0.00283907
   MO   7     0.00000000     0.00000000    -0.00023507    -0.01390207     0.00020287    -0.00387871     0.99302503     0.01148132
   MO   8     0.00000000     0.00000000     0.00151514     0.00100388     0.00032201     0.00283907     0.01148132     0.99026914
   MO   9     0.00000000     0.00000000     0.00006322     0.00002176    -0.00000652     0.00005341     0.00267930    -0.02646326
   MO  10     0.00000000     0.00000000     0.00199017    -0.00482820    -0.00125030    -0.00733165     0.00027734     0.00261463
   MO  11     0.00000000     0.00000000    -0.01079006    -0.00870963     0.00577413     0.01453823    -0.00970389    -0.00409979
   MO  12     0.00000000     0.00000000    -0.01018284     0.01712549    -0.00166261     0.00393076     0.01033759     0.00308596
   MO  13     0.00000000     0.00000000    -0.00111466     0.00524061     0.01177540     0.00360444     0.00876626     0.00837196
   MO  14     0.00000000     0.00000000    -0.00103878     0.01368422     0.00547219     0.01519827    -0.01002970     0.00266241
   MO  15     0.00000000     0.00000000    -0.00114333     0.00893856     0.00253033     0.01163966    -0.00050508    -0.03297583
   MO  16     0.00000000     0.00000000     0.00014471     0.00862212     0.00018191     0.00984750     0.00778262    -0.02044240
   MO  17     0.00000000     0.00000000    -0.00921086    -0.00288012     0.00986307     0.01226539    -0.01020707     0.00249114
   MO  18     0.00000000     0.00000000     0.00521118     0.00724698    -0.00637221    -0.00815864     0.00796975     0.00780808
   MO  19     0.00000000     0.00000000     0.01182929    -0.00234708     0.00283670    -0.00322494    -0.00182469     0.00626134
   MO  20     0.00000000     0.00000000    -0.00939861     0.01698955     0.00161194     0.00439963     0.01215084     0.00613397
   MO  21     0.00000000     0.00000000     0.00121854    -0.02269490    -0.01160998    -0.01971020    -0.01350976     0.00422034
   MO  22     0.00000000     0.00000000    -0.00057193    -0.01044783     0.00388931    -0.02592923     0.00452987    -0.00230057
   MO  23     0.00000000     0.00000000     0.00574133     0.00907490    -0.01139092     0.01925014    -0.00271233    -0.00695177
   MO  24     0.00000000     0.00000000    -0.00042515     0.01193084     0.00068694     0.01194690    -0.03154540    -0.00404315
   MO  25     0.00000000     0.00000000     0.00125915    -0.00249757     0.01011803    -0.00072990    -0.00212719     0.00290048
   MO  26     0.00000000     0.00000000    -0.00058696     0.01849774     0.00548408     0.02154804     0.00124284     0.00745927
   MO  27     0.00000000     0.00000000    -0.00201805     0.01402388    -0.00489987     0.00865701     0.00621738     0.00176331
   MO  28     0.00000000     0.00000000    -0.00392501    -0.00279088     0.00264009    -0.00232732    -0.00169544     0.00630016
   MO  29     0.00000000     0.00000000     0.00004430     0.00432605     0.00380567     0.00643019     0.00260581    -0.00144235
   MO  30     0.00000000     0.00000000    -0.00366451    -0.00158686     0.00019258    -0.00267368    -0.00112563     0.00796786
   MO  31     0.00000000     0.00000000     0.00407036     0.00274763    -0.00138264     0.00121169     0.00120679    -0.00956595
   MO  32     0.00000000     0.00000000    -0.00455910    -0.00223550    -0.00745743    -0.00748789    -0.00021026     0.00943252
   MO  33     0.00000000     0.00000000     0.00162940    -0.00615438    -0.01024210    -0.00033065    -0.00351337    -0.00172516
   MO  34     0.00000000     0.00000000    -0.00224054    -0.00434707     0.00108249    -0.00398753    -0.00015685     0.02080748
   MO  35     0.00000000     0.00000000     0.00043203     0.00486764     0.00057094     0.00207851     0.00212754    -0.00251925
   MO  36     0.00000000     0.00000000     0.00701228    -0.00669237    -0.00590259    -0.00348761    -0.00254464     0.01126513
   MO  37     0.00000000     0.00000000     0.00766921    -0.00519586     0.00526776    -0.00131916    -0.00189298     0.01383626
   MO  38     0.00000000     0.00000000     0.00856300     0.01133207    -0.00259477    -0.00616522     0.00936610    -0.00281027
   MO  39     0.00000000     0.00000000     0.00577741    -0.00187394    -0.00236893     0.00792387    -0.01091175     0.00058047
   MO  40     0.00000000     0.00000000    -0.00347584     0.00280243     0.00107072    -0.00339873    -0.00057494     0.00756969
   MO  41     0.00000000     0.00000000    -0.00023166     0.00178234    -0.00205850     0.00496408    -0.00450173     0.00115105
   MO  42     0.00000000     0.00000000    -0.00017618     0.00074926    -0.00036697     0.00063957    -0.00112343     0.00430573
   MO  43     0.00000000     0.00000000    -0.00014746     0.00510752     0.00763490     0.00320853     0.00042261     0.00299790
   MO  44     0.00000000     0.00000000    -0.00101885     0.00604691    -0.00072723     0.00199324    -0.00941138     0.00518732
   MO  45     0.00000000     0.00000000     0.00232567     0.00258102    -0.00391767     0.00306315    -0.00163363     0.01576667
   MO  46     0.00000000     0.00000000     0.00062007     0.00129160     0.00012107     0.00033066     0.00267359     0.00449421
   MO  47     0.00000000     0.00000000     0.00034273     0.00016050    -0.00125926    -0.00033124     0.00260646     0.01719093
   MO  48     0.00000000     0.00000000     0.00277864    -0.00298428    -0.00226640    -0.00095427    -0.00327217     0.00742434
   MO  49     0.00000000     0.00000000    -0.00134858     0.00323784     0.00100744     0.00143434     0.00523492    -0.00208108
   MO  50     0.00000000     0.00000000     0.00170203     0.00314216    -0.00074603     0.00328151     0.00680294     0.02694736
   MO  51     0.00000000     0.00000000    -0.00103821    -0.00528618     0.01203439    -0.00544848    -0.00684658     0.00683552
   MO  52     0.00000000     0.00000000     0.00537218     0.00620210     0.00053511     0.00840546    -0.00484379    -0.00552608
   MO  53     0.00000000     0.00000000     0.00076219    -0.01386642    -0.00209496    -0.01592733    -0.00864551    -0.00089343
   MO  54     0.00000000     0.00000000    -0.00253808     0.00405559    -0.00238881     0.00833546    -0.01306518     0.00402252
   MO  55     0.00000000     0.00000000    -0.00036316    -0.00002885    -0.00025245     0.00031460    -0.00053455    -0.00023574
   MO  56     0.00000000     0.00000000    -0.00001666     0.00472008     0.00055434    -0.00675316     0.00655379     0.00049499
   MO  57     0.00000000     0.00000000    -0.00046536    -0.00008908     0.00132106     0.00377204    -0.00094945    -0.00059110
   MO  58     0.00000000     0.00000000    -0.00012866     0.00109655     0.00013829    -0.00101975     0.00148876     0.00003323
   MO  59     0.00000000     0.00000000    -0.00020246     0.00022580     0.00019694     0.00039354    -0.00044869     0.00160348
   MO  60     0.00000000     0.00000000     0.00430169     0.00025858    -0.00033022     0.00026392    -0.00023163    -0.00823525
   MO  61     0.00000000     0.00000000     0.00292584    -0.00027159    -0.00052796    -0.00013437    -0.00085665     0.00165420
   MO  62     0.00000000     0.00000000     0.00033795     0.00067332     0.00366410     0.00109946    -0.00130697    -0.00071298
   MO  63     0.00000000     0.00000000    -0.00006139    -0.00306319     0.00129031     0.00013831     0.00028649    -0.00112732
   MO  64     0.00000000     0.00000000    -0.00005440    -0.00302686    -0.00020169    -0.00192836    -0.00042793     0.00037817
   MO  65     0.00000000     0.00000000     0.00066518    -0.00001382    -0.00031449    -0.00011264    -0.00015959     0.00133711
   MO  66     0.00000000     0.00000000    -0.00065267    -0.00007341     0.00018416     0.00010501     0.00008570    -0.00124205
   MO  67     0.00000000     0.00000000    -0.00094841    -0.00006387    -0.00000964     0.00002484     0.00009038    -0.00098609
   MO  68     0.00000000     0.00000000    -0.00401774     0.00100762     0.00065246     0.00022787    -0.00023032    -0.00449302
   MO  69     0.00000000     0.00000000    -0.00001702     0.00075992    -0.00014837     0.00392641    -0.00072080     0.00015672
   MO  70     0.00000000     0.00000000    -0.00130146    -0.00122690    -0.00044489     0.00034723     0.00186607    -0.00093572
   MO  71     0.00000000     0.00000000     0.00542374    -0.00035372    -0.00352914    -0.00034440     0.00021294     0.00137878
   MO  72     0.00000000     0.00000000    -0.00083959    -0.00013384     0.00009881    -0.00023140    -0.00036066    -0.00023499
   MO  73     0.00000000     0.00000000    -0.00281353     0.00057196    -0.00428857     0.00027600     0.00002951    -0.00127758
   MO  74     0.00000000     0.00000000    -0.00064415    -0.00233385    -0.00091856    -0.00217216    -0.00344984    -0.00029571
   MO  75     0.00000000     0.00000000     0.00015895     0.00447261    -0.00013235     0.00118427    -0.00522824    -0.00000622
   MO  76     0.00000000     0.00000000     0.00019127     0.00011557    -0.00008663     0.00006576    -0.00009230     0.00018996
   MO  77     0.00000000     0.00000000     0.00003730     0.00040518     0.00011641     0.00059121     0.00043925    -0.00000542
   MO  78     0.00000000     0.00000000    -0.00021616     0.00211625    -0.00014632     0.00315006     0.00283056     0.00002249
   MO  79     0.00000000     0.00000000     0.00043446    -0.00031968     0.00273048    -0.00024555    -0.00079361     0.00016085
   MO  80     0.00000000     0.00000000    -0.00006715    -0.00532371    -0.00003751    -0.00429836    -0.00628945    -0.00021617
   MO  81     0.00000000     0.00000000     0.00023353    -0.00008607     0.00046545    -0.00003942    -0.00040734    -0.00215112
   MO  82     0.00000000     0.00000000     0.00003512     0.00153006    -0.00009948     0.00133183    -0.00176742     0.00016353

                MO   9         MO  10         MO  11         MO  12         MO  13         MO  14         MO  15         MO  16
   MO   3     0.00006322     0.00199017    -0.01079006    -0.01018284    -0.00111466    -0.00103878    -0.00114333     0.00014471
   MO   4     0.00002176    -0.00482820    -0.00870963     0.01712549     0.00524061     0.01368422     0.00893856     0.00862212
   MO   5    -0.00000652    -0.00125030     0.00577413    -0.00166261     0.01177540     0.00547219     0.00253033     0.00018191
   MO   6     0.00005341    -0.00733165     0.01453823     0.00393076     0.00360444     0.01519827     0.01163966     0.00984750
   MO   7     0.00267930     0.00027734    -0.00970389     0.01033759     0.00876626    -0.01002970    -0.00050508     0.00778262
   MO   8    -0.02646326     0.00261463    -0.00409979     0.00308596     0.00837196     0.00266241    -0.03297583    -0.02044240
   MO   9     1.96236391    -0.00373454     0.00013639     0.00028759     0.00153063     0.00085736    -0.00358124    -0.00041287
   MO  10    -0.00373454     0.01021379    -0.00008112    -0.00028141    -0.00058233     0.00060147    -0.00222958     0.00292407
   MO  11     0.00013639    -0.00008112     0.00224965    -0.00082640    -0.00032442     0.00037126     0.00028101    -0.00016283
   MO  12     0.00028759    -0.00028141    -0.00082640     0.00191664     0.00064988     0.00007726     0.00027886     0.00060575
   MO  13     0.00153063    -0.00058233    -0.00032442     0.00064988     0.00104542    -0.00000287    -0.00012393     0.00016682
   MO  14     0.00085736     0.00060147     0.00037126     0.00007726    -0.00000287     0.00150655     0.00006969     0.00011820
   MO  15    -0.00358124    -0.00222958     0.00028101     0.00027886    -0.00012393     0.00006969     0.00275021     0.00050070
   MO  16    -0.00041287     0.00292407    -0.00016283     0.00060575     0.00016682     0.00011820     0.00050070     0.00240004
   MO  17     0.00260656     0.00007005     0.00235048    -0.00018432     0.00003768     0.00077993    -0.00000863    -0.00022422
   MO  18     0.00495117    -0.00069573    -0.00147212     0.00113882     0.00052061    -0.00009080    -0.00084543     0.00015132
   MO  19    -0.00366157     0.00069558    -0.00062163    -0.00073099    -0.00008297     0.00004515     0.00009867    -0.00042318
   MO  20    -0.00567929     0.00044791    -0.00078554     0.00225043     0.00062537    -0.00021153     0.00109911     0.00060474
   MO  21    -0.00341278     0.00070683     0.00069656    -0.00162929    -0.00150002    -0.00058630    -0.00052599    -0.00156267
   MO  22     0.00240901     0.00034065    -0.00095978    -0.00007265     0.00025889    -0.00064262    -0.00081584    -0.00037663
   MO  23     0.00050889    -0.00009285    -0.00029218    -0.00016459    -0.00032231     0.00072365     0.00048395     0.00044681
   MO  24    -0.00010302     0.00014856     0.00102481    -0.00085089    -0.00121395     0.00280831     0.00059430    -0.00133647
   MO  25    -0.00400449     0.00020192     0.00014494    -0.00058406     0.00029683     0.00017213     0.00025271    -0.00030779
   MO  26     0.00085483    -0.00041282    -0.00043022     0.00122903     0.00098382     0.00154501     0.00062403     0.00053182
   MO  27     0.00100360     0.00010487    -0.00097760     0.00160648     0.00021186     0.00062301     0.00030482     0.00051508
   MO  28     0.01095128    -0.00006898     0.00022197    -0.00026338     0.00025095     0.00020432    -0.00117914    -0.00017466
   MO  29    -0.00027853    -0.00019341    -0.00023027     0.00002773     0.00031245     0.00043600     0.00026526     0.00008491
   MO  30    -0.00032116     0.00014532     0.00002643    -0.00004987     0.00005097    -0.00000306    -0.00064100    -0.00027818
   MO  31     0.00465142    -0.00040867    -0.00011162     0.00035435     0.00003340     0.00003556     0.00034705     0.00044307
   MO  32    -0.00323993    -0.00038092    -0.00045222     0.00044244    -0.00032822    -0.00034364    -0.00035761    -0.00078565
   MO  33     0.00023954     0.00033235     0.00028319    -0.00060236    -0.00078806    -0.00024372    -0.00025934    -0.00003935
   MO  34     0.00278558     0.00059107     0.00022218    -0.00003639     0.00023294     0.00006761    -0.00150672    -0.00064615
   MO  35    -0.00012757     0.00135663    -0.00026277     0.00055114     0.00014292     0.00011513    -0.00000044     0.00082041
   MO  36     0.00112708     0.00000974     0.00025186    -0.00049468    -0.00016932    -0.00010800    -0.00068586    -0.00055552
   MO  37     0.00044410    -0.00020805     0.00014704    -0.00051720     0.00026298     0.00001411    -0.00056089    -0.00071382
   MO  38     0.00011363     0.00006430    -0.00251662     0.00083970     0.00036233    -0.00014540     0.00000449     0.00068754
   MO  39    -0.00000407     0.00004116     0.00063333    -0.00133494    -0.00071131     0.00073408    -0.00009843    -0.00065428
   MO  40     0.00032968    -0.00004817    -0.00011185     0.00089221     0.00018954     0.00000174    -0.00021679    -0.00016815
   MO  41    -0.00018676    -0.00007301     0.00068141    -0.00035179    -0.00059276     0.00045742     0.00006856    -0.00016403
   MO  42     0.00020485    -0.00003318     0.00012607     0.00005492    -0.00008981    -0.00016619    -0.00012661    -0.00009381
   MO  43    -0.00008920     0.00004743     0.00001478     0.00022543     0.00055654     0.00066615     0.00007443     0.00022647
   MO  44     0.00031482    -0.00005653    -0.00008092     0.00027964    -0.00036928     0.00141042     0.00015044    -0.00119112
   MO  45    -0.00036339    -0.00051948     0.00046765    -0.00041440    -0.00010592    -0.00032553    -0.00101094    -0.00054287
   MO  46     0.00299352     0.00215547    -0.00017962     0.00002043    -0.00002886    -0.00005347    -0.00031195     0.00051845
   MO  47    -0.00415098     0.00151001    -0.00040603     0.00036968     0.00028812     0.00030844    -0.00167831     0.00007745
   MO  48     0.00004515    -0.00031634     0.00003558    -0.00047503    -0.00022573    -0.00006676    -0.00047856    -0.00056920
   MO  49    -0.00240776    -0.00111868    -0.00035487     0.00060084     0.00043481    -0.00000525     0.00027333     0.00008139
   MO  50     0.00414548    -0.00087178    -0.00080395     0.00075716     0.00066138    -0.00001516    -0.00087059    -0.00084764
   MO  51    -0.00041311     0.00006498     0.00073784    -0.00099181     0.00073166     0.00047177    -0.00049499    -0.00080737
   MO  52     0.00008633    -0.00011024    -0.00064882    -0.00015995    -0.00001890     0.00113852     0.00071639     0.00028144
   MO  53     0.00028285     0.00021689     0.00031550    -0.00163557    -0.00113089    -0.00110173    -0.00099908    -0.00119792
   MO  54    -0.00033959    -0.00017051     0.00060763     0.00035803    -0.00044453     0.00141481     0.00019662    -0.00029525
   MO  55     0.01820525     0.00025910    -0.00002522    -0.00001233    -0.00047985    -0.00045899     0.00164376    -0.00047975
   MO  56     0.00016806     0.00015663    -0.00117879     0.00035317    -0.00007544    -0.00032376    -0.00037044    -0.00051138
   MO  57     0.00085674    -0.00006604     0.00042036    -0.00008825     0.00015120     0.00025711     0.00025449     0.00019085
   MO  58     0.00012987    -0.00008742    -0.00020910     0.00008592    -0.00006563    -0.00001776    -0.00004145    -0.00012252
   MO  59    -0.00014420    -0.00055644    -0.00005735    -0.00003181     0.00007514    -0.00001070     0.00004769    -0.00017656
   MO  60    -0.00058408    -0.00009098     0.00026548     0.00013796    -0.00003927    -0.00006986     0.00052331     0.00028272
   MO  61     0.00072422     0.00007027    -0.00056933    -0.00044651     0.00000302    -0.00007559    -0.00009821    -0.00003693
   MO  62     0.00010047     0.00003539    -0.00000658     0.00017381    -0.00004583     0.00007129     0.00010318     0.00013968
   MO  63    -0.00036927    -0.00027432     0.00041495    -0.00041080    -0.00007103    -0.00011963     0.00011266    -0.00003907
   MO  64     0.00002196     0.00008907    -0.00016765     0.00019384     0.00022581    -0.00041372    -0.00013149    -0.00005727
   MO  65    -0.00606633    -0.00021201    -0.00006811    -0.00000921    -0.00000191    -0.00001788    -0.00001518    -0.00011425
   MO  66     0.00160220     0.00052456     0.00007336    -0.00001931    -0.00002057     0.00005877    -0.00008994     0.00023295
   MO  67    -0.00250423     0.00092806     0.00002804    -0.00004338    -0.00008486     0.00004327    -0.00008731     0.00027762
   MO  68    -0.00060533    -0.00037191    -0.00008767     0.00022573     0.00003815    -0.00001886     0.00035479     0.00010366
   MO  69    -0.00004366    -0.00006345     0.00066650    -0.00043121    -0.00014469     0.00035746     0.00013803    -0.00039853
   MO  70    -0.00001430    -0.00005312     0.00023511     0.00005421     0.00011911     0.00003554     0.00019043     0.00020690
   MO  71     0.00015765     0.00003263    -0.00031375    -0.00016917    -0.00006813    -0.00002780    -0.00017186    -0.00007344
   MO  72     0.00008953     0.00007748    -0.00002404     0.00005509    -0.00000546     0.00000019    -0.00000611     0.00001034
   MO  73    -0.00007021    -0.00000673     0.00020645    -0.00000567    -0.00009603    -0.00007119     0.00005544     0.00000628
   MO  74    -0.00003498     0.00008473    -0.00046899     0.00022637    -0.00006865    -0.00000581    -0.00008349    -0.00014153
   MO  75     0.00002178    -0.00006444     0.00038525    -0.00025551    -0.00010380    -0.00002783    -0.00001885     0.00011242
   MO  76    -0.00108236    -0.00000099    -0.00000934    -0.00001949    -0.00002644    -0.00003299     0.00005471    -0.00002161
   MO  77     0.00083315     0.00001229    -0.00000118     0.00000918    -0.00001143    -0.00005647     0.00014424     0.00001486
   MO  78    -0.00013778    -0.00003827     0.00001270     0.00011830     0.00012941    -0.00009343     0.00002172     0.00019727
   MO  79    -0.00006065     0.00001250     0.00006308    -0.00007294     0.00021983     0.00011234     0.00001758    -0.00001868
   MO  80    -0.00002388     0.00005444    -0.00011060    -0.00010031    -0.00015110    -0.00008357    -0.00014681    -0.00014672
   MO  81    -0.00004308    -0.00000002    -0.00027991    -0.00025589    -0.00002527     0.00002780    -0.00001022     0.00001482
   MO  82     0.00001117    -0.00004641    -0.00005549     0.00019869     0.00002558     0.00022235     0.00011713     0.00013050

                MO  17         MO  18         MO  19         MO  20         MO  21         MO  22         MO  23         MO  24
   MO   3    -0.00921086     0.00521118     0.01182929    -0.00939861     0.00121854    -0.00057193     0.00574133    -0.00042515
   MO   4    -0.00288012     0.00724698    -0.00234708     0.01698955    -0.02269490    -0.01044783     0.00907490     0.01193084
   MO   5     0.00986307    -0.00637221     0.00283670     0.00161194    -0.01160998     0.00388931    -0.01139092     0.00068694
   MO   6     0.01226539    -0.00815864    -0.00322494     0.00439963    -0.01971020    -0.02592923     0.01925014     0.01194690
   MO   7    -0.01020707     0.00796975    -0.00182469     0.01215084    -0.01350976     0.00452987    -0.00271233    -0.03154540
   MO   8     0.00249114     0.00780808     0.00626134     0.00613397     0.00422034    -0.00230057    -0.00695177    -0.00404315
   MO   9     0.00260656     0.00495117    -0.00366157    -0.00567929    -0.00341278     0.00240901     0.00050889    -0.00010302
   MO  10     0.00007005    -0.00069573     0.00069558     0.00044791     0.00070683     0.00034065    -0.00009285     0.00014856
   MO  11     0.00235048    -0.00147212    -0.00062163    -0.00078554     0.00069656    -0.00095978    -0.00029218     0.00102481
   MO  12    -0.00018432     0.00113882    -0.00073099     0.00225043    -0.00162929    -0.00007265    -0.00016459    -0.00085089
   MO  13     0.00003768     0.00052061    -0.00008297     0.00062537    -0.00150002     0.00025889    -0.00032231    -0.00121395
   MO  14     0.00077993    -0.00009080     0.00004515    -0.00021153    -0.00058630    -0.00064262     0.00072365     0.00280831
   MO  15    -0.00000863    -0.00084543     0.00009867     0.00109911    -0.00052599    -0.00081584     0.00048395     0.00059430
   MO  16    -0.00022422     0.00015132    -0.00042318     0.00060474    -0.00156267    -0.00037663     0.00044681    -0.00133647
   MO  17     0.00468847    -0.00025681    -0.00129432     0.00002533     0.00016829    -0.00072687    -0.00139070     0.00216424
   MO  18    -0.00025681     0.00283402    -0.00096088    -0.00011131    -0.00147125     0.00099113    -0.00010918    -0.00074145
   MO  19    -0.00129432    -0.00096088     0.00185067     0.00037203     0.00096691    -0.00005047     0.00043543     0.00058747
   MO  20     0.00002533    -0.00011131     0.00037203     0.00519237    -0.00080413    -0.00104472    -0.00110703    -0.00123284
   MO  21     0.00016829    -0.00147125     0.00096691    -0.00080413     0.00400027     0.00042028    -0.00065884     0.00168265
   MO  22    -0.00072687     0.00099113    -0.00005047    -0.00104472     0.00042028     0.00263319    -0.00094527    -0.00045288
   MO  23    -0.00139070    -0.00010918     0.00043543    -0.00110703    -0.00065884    -0.00094527     0.00386364     0.00167601
   MO  24     0.00216424    -0.00074145     0.00058747    -0.00123284     0.00168265    -0.00045288     0.00167601     0.00879287
   MO  25    -0.00047977    -0.00123955     0.00108996     0.00012761     0.00048995     0.00002094    -0.00038623     0.00031183
   MO  26    -0.00013101     0.00079622     0.00028773     0.00131896    -0.00218058    -0.00077122     0.00124954     0.00152486
   MO  27    -0.00009152     0.00152281    -0.00041434     0.00184076    -0.00134055    -0.00047013     0.00092494     0.00086566
   MO  28     0.00081417     0.00114816    -0.00108695    -0.00219552    -0.00074957     0.00094475    -0.00008576     0.00010515
   MO  29    -0.00047135    -0.00019074     0.00038064    -0.00015488    -0.00064122     0.00009981     0.00054835     0.00055070
   MO  30     0.00003173     0.00014995    -0.00010216    -0.00025798     0.00017487     0.00015039    -0.00016048    -0.00012892
   MO  31     0.00058528     0.00103192    -0.00081018    -0.00033088    -0.00077189     0.00019529     0.00016607     0.00012101
   MO  32    -0.00053175     0.00026408     0.00013794     0.00106734     0.00105169    -0.00006128    -0.00001056     0.00011923
   MO  33    -0.00053002    -0.00030596    -0.00009552    -0.00120690     0.00120477    -0.00011508     0.00078116     0.00015054
   MO  34    -0.00001786     0.00000272    -0.00021302    -0.00058574     0.00021552     0.00048916    -0.00039339    -0.00029992
   MO  35     0.00051177     0.00052326    -0.00017197     0.00096942    -0.00061920    -0.00022793    -0.00012371    -0.00016406
   MO  36     0.00096547     0.00039093     0.00006168    -0.00045518     0.00087210     0.00003136    -0.00072105     0.00023500
   MO  37     0.00061555    -0.00008817     0.00062082    -0.00010605     0.00036949    -0.00006332    -0.00038529     0.00015327
   MO  38    -0.00348566     0.00132424     0.00062956     0.00027979    -0.00150861     0.00076842     0.00102016    -0.00134226
   MO  39     0.00003257    -0.00087526     0.00070388    -0.00192277     0.00118181    -0.00046565     0.00169776     0.00266538
   MO  40     0.00102154     0.00064365    -0.00047383     0.00164184    -0.00018710    -0.00000470    -0.00160663    -0.00025593
   MO  41     0.00064451    -0.00051888    -0.00019970    -0.00016420     0.00040582    -0.00126808    -0.00005249     0.00094795
   MO  42     0.00018857    -0.00013829    -0.00015649     0.00044608     0.00010133    -0.00062046    -0.00048629    -0.00059104
   MO  43     0.00031135     0.00015714     0.00000314     0.00010587    -0.00083905     0.00014292    -0.00089947     0.00039366
   MO  44     0.00144770     0.00028034     0.00046372     0.00049705     0.00103168     0.00042235     0.00106231     0.00553388
   MO  45     0.00020923    -0.00045093    -0.00029366    -0.00084923     0.00037855    -0.00059004    -0.00005790    -0.00089065
   MO  46    -0.00049086    -0.00058834     0.00058109     0.00071094     0.00019811    -0.00006979     0.00003066    -0.00043770
   MO  47    -0.00008632     0.00102384    -0.00017851    -0.00013862    -0.00044092     0.00005681     0.00040457    -0.00019449
   MO  48    -0.00039921    -0.00028778     0.00018513    -0.00090178     0.00048072     0.00011725     0.00024409     0.00010878
   MO  49    -0.00005935     0.00090196    -0.00039820     0.00030552    -0.00106887     0.00007850     0.00029513    -0.00038214
   MO  50    -0.00117099     0.00048387     0.00057011     0.00126008    -0.00054147    -0.00044720     0.00055644    -0.00099533
   MO  51     0.00104183    -0.00095275     0.00064613    -0.00096987     0.00024796     0.00109520    -0.00212246     0.00071741
   MO  52    -0.00166633    -0.00018925     0.00096012    -0.00068834    -0.00072632    -0.00037424     0.00233972     0.00222445
   MO  53    -0.00055952    -0.00067637     0.00021795    -0.00214085     0.00286293     0.00108942    -0.00074274    -0.00035855
   MO  54     0.00144510     0.00021636    -0.00040213     0.00034840    -0.00003038    -0.00120153     0.00080378     0.00343577
   MO  55    -0.00127614    -0.00289225     0.00232499     0.00381989     0.00216598    -0.00155386    -0.00053130    -0.00006477
   MO  56    -0.00075907     0.00077014     0.00028942     0.00041091     0.00058895     0.00160404    -0.00085946     0.00036940
   MO  57     0.00042607    -0.00059394     0.00002651     0.00009379    -0.00080376    -0.00090467     0.00086871     0.00014170
   MO  58    -0.00006638     0.00019682    -0.00016790     0.00003083    -0.00001006    -0.00010500    -0.00004843     0.00018885
   MO  59    -0.00015587     0.00002985     0.00011525    -0.00007723    -0.00006306     0.00005545     0.00011319    -0.00007842
   MO  60     0.00104850     0.00016059    -0.00026438     0.00054536    -0.00012338    -0.00019808    -0.00062640     0.00005633
   MO  61    -0.00018575     0.00047711     0.00068576    -0.00022300     0.00004872    -0.00014752     0.00019205    -0.00003112
   MO  62     0.00002705     0.00021493    -0.00011926     0.00023041     0.00001949    -0.00018553    -0.00005718     0.00003119
   MO  63     0.00032091    -0.00037914     0.00007135    -0.00052531     0.00032613    -0.00010196     0.00006747    -0.00010157
   MO  64    -0.00004006     0.00020054     0.00002823     0.00029406     0.00008487     0.00061585     0.00012708    -0.00054683
   MO  65    -0.00010256    -0.00000580     0.00013275     0.00015855     0.00013090    -0.00014961    -0.00000211    -0.00001466
   MO  66     0.00015264     0.00007856    -0.00015615    -0.00023291    -0.00008921     0.00013217     0.00004324     0.00004425
   MO  67    -0.00004158    -0.00022000     0.00016061     0.00017683     0.00015743    -0.00007438     0.00003077     0.00003816
   MO  68    -0.00009165     0.00003762    -0.00015901     0.00022567    -0.00025151     0.00007061    -0.00001607    -0.00001556
   MO  69     0.00098604    -0.00049082     0.00020329    -0.00046968     0.00060840    -0.00010629     0.00085794     0.00179526
   MO  70     0.00010745     0.00005775    -0.00001796     0.00019387    -0.00030451    -0.00039744     0.00035710    -0.00037956
   MO  71    -0.00077954     0.00000154     0.00028314    -0.00040398     0.00004909     0.00006777     0.00061767    -0.00005125
   MO  72     0.00006981     0.00002934    -0.00000882     0.00009289     0.00004024     0.00013283    -0.00015047     0.00005543
   MO  73     0.00040114    -0.00010826    -0.00016001     0.00011024     0.00009142    -0.00018695    -0.00009378     0.00002385
   MO  74    -0.00026375     0.00038201     0.00003696     0.00027496     0.00012321     0.00048403    -0.00034102     0.00026141
   MO  75     0.00009945    -0.00039456    -0.00015660    -0.00035758    -0.00002551    -0.00046142    -0.00019052    -0.00055595
   MO  76    -0.00005200    -0.00006518     0.00000499    -0.00005608    -0.00006086     0.00006333     0.00003974    -0.00002370
   MO  77    -0.00008992    -0.00014025     0.00007481     0.00015651     0.00002465    -0.00004795     0.00002190    -0.00011510
   MO  78    -0.00006884     0.00010651    -0.00007812     0.00016268    -0.00028772    -0.00023157     0.00015464    -0.00054367
   MO  79     0.00015490    -0.00011214     0.00011741    -0.00004069    -0.00018227     0.00021454    -0.00039451     0.00003811
   MO  80    -0.00021423     0.00013003     0.00001821    -0.00011850     0.00034323     0.00017546    -0.00032757    -0.00005509
   MO  81    -0.00069847    -0.00004048     0.00034934    -0.00050241    -0.00002919     0.00009885     0.00050469     0.00003726
   MO  82    -0.00000691     0.00006055    -0.00007742     0.00022232    -0.00032592    -0.00025044     0.00008666     0.00029572

                MO  25         MO  26         MO  27         MO  28         MO  29         MO  30         MO  31         MO  32
   MO   3     0.00125915    -0.00058696    -0.00201805    -0.00392501     0.00004430    -0.00366451     0.00407036    -0.00455910
   MO   4    -0.00249757     0.01849774     0.01402388    -0.00279088     0.00432605    -0.00158686     0.00274763    -0.00223550
   MO   5     0.01011803     0.00548408    -0.00489987     0.00264009     0.00380567     0.00019258    -0.00138264    -0.00745743
   MO   6    -0.00072990     0.02154804     0.00865701    -0.00232732     0.00643019    -0.00267368     0.00121169    -0.00748789
   MO   7    -0.00212719     0.00124284     0.00621738    -0.00169544     0.00260581    -0.00112563     0.00120679    -0.00021026
   MO   8     0.00290048     0.00745927     0.00176331     0.00630016    -0.00144235     0.00796786    -0.00956595     0.00943252
   MO   9    -0.00400449     0.00085483     0.00100360     0.01095128    -0.00027853    -0.00032116     0.00465142    -0.00323993
   MO  10     0.00020192    -0.00041282     0.00010487    -0.00006898    -0.00019341     0.00014532    -0.00040867    -0.00038092
   MO  11     0.00014494    -0.00043022    -0.00097760     0.00022197    -0.00023027     0.00002643    -0.00011162    -0.00045222
   MO  12    -0.00058406     0.00122903     0.00160648    -0.00026338     0.00002773    -0.00004987     0.00035435     0.00044244
   MO  13     0.00029683     0.00098382     0.00021186     0.00025095     0.00031245     0.00005097     0.00003340    -0.00032822
   MO  14     0.00017213     0.00154501     0.00062301     0.00020432     0.00043600    -0.00000306     0.00003556    -0.00034364
   MO  15     0.00025271     0.00062403     0.00030482    -0.00117914     0.00026526    -0.00064100     0.00034705    -0.00035761
   MO  16    -0.00030779     0.00053182     0.00051508    -0.00017466     0.00008491    -0.00027818     0.00044307    -0.00078565
   MO  17    -0.00047977    -0.00013101    -0.00009152     0.00081417    -0.00047135     0.00003173     0.00058528    -0.00053175
   MO  18    -0.00123955     0.00079622     0.00152281     0.00114816    -0.00019074     0.00014995     0.00103192     0.00026408
   MO  19     0.00108996     0.00028773    -0.00041434    -0.00108695     0.00038064    -0.00010216    -0.00081018     0.00013794
   MO  20     0.00012761     0.00131896     0.00184076    -0.00219552    -0.00015488    -0.00025798    -0.00033088     0.00106734
   MO  21     0.00048995    -0.00218058    -0.00134055    -0.00074957    -0.00064122     0.00017487    -0.00077189     0.00105169
   MO  22     0.00002094    -0.00077122    -0.00047013     0.00094475     0.00009981     0.00015039     0.00019529    -0.00006128
   MO  23    -0.00038623     0.00124954     0.00092494    -0.00008576     0.00054835    -0.00016048     0.00016607    -0.00001056
   MO  24     0.00031183     0.00152486     0.00086566     0.00010515     0.00055070    -0.00012892     0.00012101     0.00011923
   MO  25     0.00245399     0.00022333    -0.00082523    -0.00084694     0.00033274    -0.00004228    -0.00077609    -0.00010216
   MO  26     0.00022333     0.00457647     0.00189472    -0.00018639     0.00069693    -0.00016070     0.00022150    -0.00011097
   MO  27    -0.00082523     0.00189472     0.00339100    -0.00024990     0.00034835    -0.00022850     0.00062697     0.00054996
   MO  28    -0.00084694    -0.00018639    -0.00024990     0.00328627    -0.00011047     0.00023282     0.00086728    -0.00076447
   MO  29     0.00033274     0.00069693     0.00034835    -0.00011047     0.00153429    -0.00001152    -0.00003952    -0.00032736
   MO  30    -0.00004228    -0.00016070    -0.00022850     0.00023282    -0.00001152     0.00140771     0.00018871     0.00006133
   MO  31    -0.00077609     0.00022150     0.00062697     0.00086728    -0.00003952     0.00018871     0.00173685    -0.00036388
   MO  32    -0.00010216    -0.00011097     0.00054996    -0.00076447    -0.00032736     0.00006133    -0.00036388     0.00177493
   MO  33    -0.00022258    -0.00075076    -0.00047481     0.00008089    -0.00008268     0.00005984    -0.00008718     0.00014419
   MO  34     0.00006653    -0.00012878    -0.00035001     0.00080851     0.00003010     0.00040061    -0.00029271     0.00014573
   MO  35    -0.00034097     0.00032149     0.00063929    -0.00013419    -0.00008400    -0.00012174     0.00040965     0.00011928
   MO  36    -0.00004003    -0.00019763    -0.00024715     0.00029474    -0.00031266     0.00009094     0.00008988     0.00005994
   MO  37     0.00020187     0.00010799    -0.00022556     0.00005723    -0.00000278     0.00008515    -0.00022560    -0.00005573
   MO  38    -0.00006026     0.00061671     0.00077578    -0.00024326     0.00048866    -0.00003504     0.00008687     0.00002965
   MO  39     0.00026795     0.00002203    -0.00046641     0.00006247     0.00040117    -0.00001660    -0.00023488    -0.00007454
   MO  40    -0.00082141     0.00015644     0.00070720    -0.00010641    -0.00017636     0.00004062     0.00010277     0.00032156
   MO  41    -0.00044387    -0.00044153    -0.00016172    -0.00009360    -0.00033562     0.00008034    -0.00005889     0.00005322
   MO  42     0.00045052    -0.00068392     0.00019748     0.00006818    -0.00034113    -0.00001788     0.00004600     0.00020121
   MO  43     0.00100565     0.00136042    -0.00021172     0.00028731     0.00030509    -0.00015192    -0.00008339    -0.00058045
   MO  44     0.00014777     0.00113504     0.00140046    -0.00000937     0.00025320    -0.00008334     0.00026698     0.00074396
   MO  45    -0.00000502    -0.00190160    -0.00211805     0.00033139    -0.00063762     0.00044298    -0.00067689    -0.00009025
   MO  46     0.00006425     0.00011512    -0.00034997     0.00056244    -0.00043143    -0.00025952    -0.00013255    -0.00019752
   MO  47    -0.00019814     0.00058131     0.00092036    -0.00030841     0.00002470     0.00016301    -0.00061275     0.00075639
   MO  48    -0.00009260    -0.00067059    -0.00005510     0.00066353     0.00049065    -0.00101330    -0.00068898    -0.00014378
   MO  49    -0.00043185     0.00102186     0.00051205     0.00031633    -0.00077001    -0.00086759    -0.00009450     0.00036518
   MO  50     0.00002762     0.00159862     0.00104899     0.00017809     0.00030431     0.00064575    -0.00002233     0.00069865
   MO  51     0.00202756     0.00023213    -0.00196413     0.00064401     0.00039879     0.00021847    -0.00071251    -0.00096716
   MO  52     0.00065327     0.00213573     0.00061885    -0.00059660     0.00083060    -0.00027159     0.00000894    -0.00038647
   MO  53    -0.00015790    -0.00321194    -0.00195146     0.00041901    -0.00109893     0.00042002    -0.00022861     0.00064162
   MO  54    -0.00083266     0.00155985     0.00140602     0.00022908     0.00000536     0.00004693     0.00023676     0.00033523
   MO  55     0.00197381    -0.00012176    -0.00031171    -0.00386252    -0.00002907    -0.00054180    -0.00185625     0.00129270
   MO  56    -0.00051874    -0.00077504     0.00059116    -0.00008174     0.00043986     0.00003216    -0.00001561     0.00043951
   MO  57    -0.00046024    -0.00032925     0.00015865    -0.00003016     0.00035283     0.00002325    -0.00002737    -0.00041564
   MO  58     0.00104610    -0.00069650     0.00089762     0.00024193    -0.00012157    -0.00006303     0.00028054     0.00019495
   MO  59    -0.00008093     0.00018767     0.00020127    -0.00001796    -0.00047439     0.00033672    -0.00047488    -0.00047597
   MO  60    -0.00032210    -0.00023015     0.00021480    -0.00061783    -0.00043523    -0.00047565     0.00074863    -0.00050752
   MO  61     0.00000681     0.00006578     0.00002859     0.00004885    -0.00009292    -0.00011075     0.00024136    -0.00011135
   MO  62    -0.00000789     0.00048263     0.00019783    -0.00005589    -0.00043584    -0.00013402     0.00023959     0.00027451
   MO  63     0.00011580    -0.00031376    -0.00035718    -0.00009861     0.00044023    -0.00020190    -0.00020084    -0.00045443
   MO  64     0.00001169    -0.00001409     0.00016875     0.00000937    -0.00013938     0.00005129     0.00007265     0.00031833
   MO  65     0.00022339    -0.00003441    -0.00004438    -0.00081044     0.00010485     0.00045984    -0.00015063     0.00032114
   MO  66    -0.00018242     0.00001479     0.00001275     0.00054203     0.00014277     0.00067666     0.00072835    -0.00040611
   MO  67     0.00020669    -0.00008793     0.00000632    -0.00050985    -0.00010259    -0.00033509    -0.00047943     0.00014592
   MO  68     0.00000517    -0.00003243     0.00001154    -0.00010250    -0.00001476    -0.00007622     0.00012078    -0.00006685
   MO  69     0.00023609     0.00022980    -0.00004225     0.00003981     0.00013576    -0.00007499     0.00001284     0.00001545
   MO  70    -0.00021563     0.00107769     0.00097353    -0.00019509     0.00027821    -0.00012989     0.00016027     0.00009229
   MO  71    -0.00001373    -0.00001463    -0.00000664    -0.00005414     0.00015864     0.00005346    -0.00010346     0.00006433
   MO  72    -0.00020639    -0.00001476    -0.00009605    -0.00004055     0.00003563     0.00002355    -0.00003033    -0.00002462
   MO  73    -0.00022799    -0.00046320    -0.00018358    -0.00000083    -0.00026116     0.00005285     0.00005651     0.00002009
   MO  74    -0.00012674    -0.00013374     0.00068534    -0.00008911     0.00038877    -0.00007564     0.00006358     0.00017302
   MO  75     0.00001369    -0.00049762    -0.00037447     0.00000948     0.00019528    -0.00000370    -0.00017142    -0.00032837
   MO  76    -0.00025968    -0.00000337     0.00003908     0.00073730    -0.00019912    -0.00082765    -0.00014635    -0.00015432
   MO  77    -0.00011514     0.00007492    -0.00007083     0.00044196     0.00003147     0.00046505     0.00054167    -0.00023814
   MO  78    -0.00004942     0.00026650    -0.00025617    -0.00008107    -0.00054958     0.00003104    -0.00002546     0.00014281
   MO  79     0.00043852     0.00027672    -0.00029364     0.00007262     0.00016324     0.00001915    -0.00012801    -0.00029151
   MO  80    -0.00010235    -0.00004394     0.00001789     0.00002699    -0.00024695     0.00008718     0.00001062     0.00025652
   MO  81     0.00020424     0.00007546    -0.00018061    -0.00010944     0.00015387    -0.00011544     0.00004492    -0.00023348
   MO  82    -0.00005607     0.00026392     0.00023146    -0.00004050     0.00013073    -0.00003164     0.00001228    -0.00006293

                MO  33         MO  34         MO  35         MO  36         MO  37         MO  38         MO  39         MO  40
   MO   3     0.00162940    -0.00224054     0.00043203     0.00701228     0.00766921     0.00856300     0.00577741    -0.00347584
   MO   4    -0.00615438    -0.00434707     0.00486764    -0.00669237    -0.00519586     0.01133207    -0.00187394     0.00280243
   MO   5    -0.01024210     0.00108249     0.00057094    -0.00590259     0.00526776    -0.00259477    -0.00236893     0.00107072
   MO   6    -0.00033065    -0.00398753     0.00207851    -0.00348761    -0.00131916    -0.00616522     0.00792387    -0.00339873
   MO   7    -0.00351337    -0.00015685     0.00212754    -0.00254464    -0.00189298     0.00936610    -0.01091175    -0.00057494
   MO   8    -0.00172516     0.02080748    -0.00251925     0.01126513     0.01383626    -0.00281027     0.00058047     0.00756969
   MO   9     0.00023954     0.00278558    -0.00012757     0.00112708     0.00044410     0.00011363    -0.00000407     0.00032968
   MO  10     0.00033235     0.00059107     0.00135663     0.00000974    -0.00020805     0.00006430     0.00004116    -0.00004817
   MO  11     0.00028319     0.00022218    -0.00026277     0.00025186     0.00014704    -0.00251662     0.00063333    -0.00011185
   MO  12    -0.00060236    -0.00003639     0.00055114    -0.00049468    -0.00051720     0.00083970    -0.00133494     0.00089221
   MO  13    -0.00078806     0.00023294     0.00014292    -0.00016932     0.00026298     0.00036233    -0.00071131     0.00018954
   MO  14    -0.00024372     0.00006761     0.00011513    -0.00010800     0.00001411    -0.00014540     0.00073408     0.00000174
   MO  15    -0.00025934    -0.00150672    -0.00000044    -0.00068586    -0.00056089     0.00000449    -0.00009843    -0.00021679
   MO  16    -0.00003935    -0.00064615     0.00082041    -0.00055552    -0.00071382     0.00068754    -0.00065428    -0.00016815
   MO  17    -0.00053002    -0.00001786     0.00051177     0.00096547     0.00061555    -0.00348566     0.00003257     0.00102154
   MO  18    -0.00030596     0.00000272     0.00052326     0.00039093    -0.00008817     0.00132424    -0.00087526     0.00064365
   MO  19    -0.00009552    -0.00021302    -0.00017197     0.00006168     0.00062082     0.00062956     0.00070388    -0.00047383
   MO  20    -0.00120690    -0.00058574     0.00096942    -0.00045518    -0.00010605     0.00027979    -0.00192277     0.00164184
   MO  21     0.00120477     0.00021552    -0.00061920     0.00087210     0.00036949    -0.00150861     0.00118181    -0.00018710
   MO  22    -0.00011508     0.00048916    -0.00022793     0.00003136    -0.00006332     0.00076842    -0.00046565    -0.00000470
   MO  23     0.00078116    -0.00039339    -0.00012371    -0.00072105    -0.00038529     0.00102016     0.00169776    -0.00160663
   MO  24     0.00015054    -0.00029992    -0.00016406     0.00023500     0.00015327    -0.00134226     0.00266538    -0.00025593
   MO  25    -0.00022258     0.00006653    -0.00034097    -0.00004003     0.00020187    -0.00006026     0.00026795    -0.00082141
   MO  26    -0.00075076    -0.00012878     0.00032149    -0.00019763     0.00010799     0.00061671     0.00002203     0.00015644
   MO  27    -0.00047481    -0.00035001     0.00063929    -0.00024715    -0.00022556     0.00077578    -0.00046641     0.00070720
   MO  28     0.00008089     0.00080851    -0.00013419     0.00029474     0.00005723    -0.00024326     0.00006247    -0.00010641
   MO  29    -0.00008268     0.00003010    -0.00008400    -0.00031266    -0.00000278     0.00048866     0.00040117    -0.00017636
   MO  30     0.00005984     0.00040061    -0.00012174     0.00009094     0.00008515    -0.00003504    -0.00001660     0.00004062
   MO  31    -0.00008718    -0.00029271     0.00040965     0.00008988    -0.00022560     0.00008687    -0.00023488     0.00010277
   MO  32     0.00014419     0.00014573     0.00011928     0.00005994    -0.00005573     0.00002965    -0.00007454     0.00032156
   MO  33     0.00171828     0.00004337    -0.00039951     0.00030022    -0.00030911    -0.00020920     0.00080120    -0.00056143
   MO  34     0.00004337     0.00151880    -0.00038458     0.00012419     0.00026509    -0.00025025     0.00003186     0.00014436
   MO  35    -0.00039951    -0.00038458     0.00085695     0.00002386    -0.00004579     0.00009504    -0.00047831     0.00040896
   MO  36     0.00030022     0.00012419     0.00002386     0.00120206     0.00061928    -0.00092443     0.00018219     0.00032711
   MO  37    -0.00030911     0.00026509    -0.00004579     0.00061928     0.00105633    -0.00062429     0.00026966     0.00023774
   MO  38    -0.00020920    -0.00025025     0.00009504    -0.00092443    -0.00062429     0.00374248    -0.00054248    -0.00038596
   MO  39     0.00080120     0.00003186    -0.00047831     0.00018219     0.00026966    -0.00054248     0.00206476    -0.00106541
   MO  40    -0.00056143     0.00014436     0.00040896     0.00032711     0.00023774    -0.00038596    -0.00106541     0.00163675
   MO  41     0.00014608    -0.00004235    -0.00002857     0.00009245     0.00002630    -0.00053864     0.00043374     0.00026593
   MO  42    -0.00005486     0.00007942     0.00007955     0.00011151     0.00002625    -0.00020239    -0.00030181     0.00009922
   MO  43    -0.00036171     0.00012304    -0.00000526     0.00022316    -0.00001890     0.00013042    -0.00037721     0.00016025
   MO  44    -0.00026184    -0.00018524     0.00024627     0.00022020     0.00027015    -0.00070287     0.00122988     0.00015510
   MO  45     0.00085755     0.00057106    -0.00032678     0.00029686     0.00007636    -0.00018241     0.00019960    -0.00021312
   MO  46    -0.00004650     0.00027004     0.00020249    -0.00013175     0.00006690     0.00022498    -0.00020692    -0.00001292
   MO  47    -0.00009446     0.00072565     0.00031652     0.00022677     0.00028302     0.00030981     0.00009128     0.00016783
   MO  48     0.00059911     0.00050232    -0.00049388     0.00015222     0.00021123    -0.00006986     0.00056977    -0.00008227
   MO  49    -0.00073777    -0.00025679     0.00021050    -0.00024635    -0.00009005     0.00042906    -0.00049987     0.00003779
   MO  50    -0.00025971     0.00096271    -0.00015061     0.00008406     0.00059530     0.00072622    -0.00021577     0.00033324
   MO  51    -0.00126216     0.00085589    -0.00051967     0.00045290     0.00076562    -0.00098071    -0.00009977    -0.00005241
   MO  52     0.00012854    -0.00036970    -0.00025156    -0.00064446    -0.00027981     0.00159916     0.00119906    -0.00117548
   MO  53     0.00085772     0.00055369    -0.00071737     0.00032289     0.00016192    -0.00059079     0.00071713    -0.00054358
   MO  54     0.00004900    -0.00009219     0.00023311     0.00019543     0.00007725    -0.00082501     0.00067126     0.00061223
   MO  55    -0.00016391    -0.00078210    -0.00011897    -0.00038936     0.00024569    -0.00025989    -0.00008642     0.00021896
   MO  56    -0.00027801     0.00012548     0.00002021    -0.00011347     0.00005179     0.00097982    -0.00039586     0.00072210
   MO  57    -0.00061360    -0.00004084     0.00012939    -0.00064196     0.00023471    -0.00012102     0.00033392    -0.00007937
   MO  58    -0.00012611    -0.00001894     0.00002619    -0.00006553    -0.00018222     0.00017996    -0.00009987    -0.00042815
   MO  59    -0.00000126    -0.00002695    -0.00038077     0.00000924     0.00015551     0.00012006    -0.00007186    -0.00002057
   MO  60    -0.00037510    -0.00040035     0.00033835     0.00034888     0.00019558    -0.00050378    -0.00038969     0.00052833
   MO  61    -0.00019540    -0.00062003     0.00035313     0.00055942     0.00058896     0.00039848     0.00022845    -0.00008859
   MO  62     0.00048508    -0.00019192     0.00014216     0.00000433     0.00003100    -0.00002798    -0.00005771     0.00011187
   MO  63     0.00080620    -0.00014161    -0.00024755     0.00034740     0.00018990    -0.00059543     0.00037721    -0.00011803
   MO  64    -0.00014885     0.00009314     0.00012305    -0.00005965     0.00013017    -0.00010476    -0.00033585    -0.00008850
   MO  65    -0.00000036    -0.00010690     0.00004520     0.00002696     0.00006284     0.00006913     0.00006212    -0.00001145
   MO  66     0.00014611     0.00005841    -0.00001319     0.00002686    -0.00007536    -0.00005336    -0.00001034    -0.00003125
   MO  67    -0.00001185    -0.00016320     0.00000561    -0.00007249    -0.00005002    -0.00003110     0.00002715    -0.00006475
   MO  68    -0.00017780    -0.00031339     0.00010328    -0.00035451    -0.00036375     0.00024920    -0.00024732    -0.00001382
   MO  69     0.00020421    -0.00008261    -0.00007875     0.00020818     0.00023364    -0.00116494     0.00108067    -0.00050374
   MO  70    -0.00014899    -0.00003495     0.00004280     0.00006090     0.00020433    -0.00045648     0.00009905     0.00002649
   MO  71     0.00021177     0.00010807    -0.00016669    -0.00020204    -0.00000715     0.00060004     0.00036248    -0.00037406
   MO  72    -0.00000650    -0.00000955     0.00002766     0.00001705     0.00000154     0.00000034    -0.00006266     0.00021655
   MO  73    -0.00006694    -0.00009391     0.00010651    -0.00003831     0.00004012    -0.00029397    -0.00007081     0.00013045
   MO  74    -0.00008841     0.00002716     0.00004511    -0.00002712    -0.00004894     0.00042256    -0.00013128     0.00045031
   MO  75     0.00029177     0.00011644    -0.00016689     0.00002909    -0.00001344    -0.00029299     0.00008175     0.00000346
   MO  76    -0.00004926     0.00007127    -0.00007053    -0.00007155     0.00004047    -0.00003015    -0.00005036     0.00007146
   MO  77     0.00003354     0.00005120    -0.00002764    -0.00003355     0.00000871     0.00002615    -0.00007238     0.00001791
   MO  78    -0.00012002    -0.00009957     0.00010608    -0.00003959    -0.00005081     0.00012897    -0.00020789    -0.00013408
   MO  79    -0.00046588     0.00012017    -0.00005991     0.00008656     0.00006668    -0.00006566    -0.00016085    -0.00001232
   MO  80    -0.00000840     0.00009601    -0.00005791     0.00008915     0.00008006     0.00006218    -0.00015442     0.00017134
   MO  81     0.00013053    -0.00003765    -0.00013741    -0.00006655     0.00001856     0.00057162     0.00031012    -0.00040711
   MO  82    -0.00004114    -0.00007156     0.00009388    -0.00012087    -0.00008394     0.00026144    -0.00015263     0.00012685

                MO  41         MO  42         MO  43         MO  44         MO  45         MO  46         MO  47         MO  48
   MO   3    -0.00023166    -0.00017618    -0.00014746    -0.00101885     0.00232567     0.00062007     0.00034273     0.00277864
   MO   4     0.00178234     0.00074926     0.00510752     0.00604691     0.00258102     0.00129160     0.00016050    -0.00298428
   MO   5    -0.00205850    -0.00036697     0.00763490    -0.00072723    -0.00391767     0.00012107    -0.00125926    -0.00226640
   MO   6     0.00496408     0.00063957     0.00320853     0.00199324     0.00306315     0.00033066    -0.00033124    -0.00095427
   MO   7    -0.00450173    -0.00112343     0.00042261    -0.00941138    -0.00163363     0.00267359     0.00260646    -0.00327217
   MO   8     0.00115105     0.00430573     0.00299790     0.00518732     0.01576667     0.00449421     0.01719093     0.00742434
   MO   9    -0.00018676     0.00020485    -0.00008920     0.00031482    -0.00036339     0.00299352    -0.00415098     0.00004515
   MO  10    -0.00007301    -0.00003318     0.00004743    -0.00005653    -0.00051948     0.00215547     0.00151001    -0.00031634
   MO  11     0.00068141     0.00012607     0.00001478    -0.00008092     0.00046765    -0.00017962    -0.00040603     0.00003558
   MO  12    -0.00035179     0.00005492     0.00022543     0.00027964    -0.00041440     0.00002043     0.00036968    -0.00047503
   MO  13    -0.00059276    -0.00008981     0.00055654    -0.00036928    -0.00010592    -0.00002886     0.00028812    -0.00022573
   MO  14     0.00045742    -0.00016619     0.00066615     0.00141042    -0.00032553    -0.00005347     0.00030844    -0.00006676
   MO  15     0.00006856    -0.00012661     0.00007443     0.00015044    -0.00101094    -0.00031195    -0.00167831    -0.00047856
   MO  16    -0.00016403    -0.00009381     0.00022647    -0.00119112    -0.00054287     0.00051845     0.00007745    -0.00056920
   MO  17     0.00064451     0.00018857     0.00031135     0.00144770     0.00020923    -0.00049086    -0.00008632    -0.00039921
   MO  18    -0.00051888    -0.00013829     0.00015714     0.00028034    -0.00045093    -0.00058834     0.00102384    -0.00028778
   MO  19    -0.00019970    -0.00015649     0.00000314     0.00046372    -0.00029366     0.00058109    -0.00017851     0.00018513
   MO  20    -0.00016420     0.00044608     0.00010587     0.00049705    -0.00084923     0.00071094    -0.00013862    -0.00090178
   MO  21     0.00040582     0.00010133    -0.00083905     0.00103168     0.00037855     0.00019811    -0.00044092     0.00048072
   MO  22    -0.00126808    -0.00062046     0.00014292     0.00042235    -0.00059004    -0.00006979     0.00005681     0.00011725
   MO  23    -0.00005249    -0.00048629    -0.00089947     0.00106231    -0.00005790     0.00003066     0.00040457     0.00024409
   MO  24     0.00094795    -0.00059104     0.00039366     0.00553388    -0.00089065    -0.00043770    -0.00019449     0.00010878
   MO  25    -0.00044387     0.00045052     0.00100565     0.00014777    -0.00000502     0.00006425    -0.00019814    -0.00009260
   MO  26    -0.00044153    -0.00068392     0.00136042     0.00113504    -0.00190160     0.00011512     0.00058131    -0.00067059
   MO  27    -0.00016172     0.00019748    -0.00021172     0.00140046    -0.00211805    -0.00034997     0.00092036    -0.00005510
   MO  28    -0.00009360     0.00006818     0.00028731    -0.00000937     0.00033139     0.00056244    -0.00030841     0.00066353
   MO  29    -0.00033562    -0.00034113     0.00030509     0.00025320    -0.00063762    -0.00043143     0.00002470     0.00049065
   MO  30     0.00008034    -0.00001788    -0.00015192    -0.00008334     0.00044298    -0.00025952     0.00016301    -0.00101330
   MO  31    -0.00005889     0.00004600    -0.00008339     0.00026698    -0.00067689    -0.00013255    -0.00061275    -0.00068898
   MO  32     0.00005322     0.00020121    -0.00058045     0.00074396    -0.00009025    -0.00019752     0.00075639    -0.00014378
   MO  33     0.00014608    -0.00005486    -0.00036171    -0.00026184     0.00085755    -0.00004650    -0.00009446     0.00059911
   MO  34    -0.00004235     0.00007942     0.00012304    -0.00018524     0.00057106     0.00027004     0.00072565     0.00050232
   MO  35    -0.00002857     0.00007955    -0.00000526     0.00024627    -0.00032678     0.00020249     0.00031652    -0.00049388
   MO  36     0.00009245     0.00011151     0.00022316     0.00022020     0.00029686    -0.00013175     0.00022677     0.00015222
   MO  37     0.00002630     0.00002625    -0.00001890     0.00027015     0.00007636     0.00006690     0.00028302     0.00021123
   MO  38    -0.00053864    -0.00020239     0.00013042    -0.00070287    -0.00018241     0.00022498     0.00030981    -0.00006986
   MO  39     0.00043374    -0.00030181    -0.00037721     0.00122988     0.00019960    -0.00020692     0.00009128     0.00056977
   MO  40     0.00026593     0.00009922     0.00016025     0.00015510    -0.00021312    -0.00001292     0.00016783    -0.00008227
   MO  41     0.00141268     0.00041819    -0.00020977    -0.00020000     0.00068229    -0.00000580    -0.00012087     0.00004154
   MO  42     0.00041819     0.00100041     0.00001529    -0.00043014     0.00041334     0.00016948    -0.00008770     0.00017159
   MO  43    -0.00020977     0.00001529     0.00166912    -0.00001957    -0.00014383     0.00005812    -0.00010396    -0.00008765
   MO  44    -0.00020000    -0.00043014    -0.00001957     0.00515532    -0.00081360    -0.00009000     0.00016071    -0.00015967
   MO  45     0.00068229     0.00041334    -0.00014383    -0.00081360     0.00422454     0.00017049     0.00015902     0.00034725
   MO  46    -0.00000580     0.00016948     0.00005812    -0.00009000     0.00017049     0.00209039    -0.00081937     0.00030665
   MO  47    -0.00012087    -0.00008770    -0.00010396     0.00016071     0.00015902    -0.00081937     0.00279611     0.00009146
   MO  48     0.00004154     0.00017159    -0.00008765    -0.00015967     0.00034725     0.00030665     0.00009146     0.00254619
   MO  49    -0.00018180    -0.00015577     0.00013414     0.00010688    -0.00043239     0.00020035     0.00053203     0.00000329
   MO  50    -0.00021626     0.00025297    -0.00003969     0.00007621    -0.00014825     0.00066118     0.00051610     0.00001645
   MO  51    -0.00040117    -0.00009153     0.00173949    -0.00005823    -0.00012157     0.00011148    -0.00039279    -0.00006293
   MO  52    -0.00012476    -0.00061721     0.00047459     0.00105689    -0.00093716    -0.00003904     0.00008079    -0.00010488
   MO  53     0.00046427     0.00021345    -0.00131243    -0.00075242     0.00087236     0.00010931    -0.00032926     0.00028995
   MO  54     0.00103990    -0.00007076     0.00009412     0.00208895    -0.00048917    -0.00017926     0.00030182    -0.00009239
   MO  55     0.00032617     0.00054685    -0.00025925     0.00023669    -0.00028646     0.00096720    -0.00105547    -0.00013710
   MO  56    -0.00049289    -0.00040366    -0.00037627     0.00105005    -0.00078257    -0.00015553     0.00013209     0.00026935
   MO  57     0.00058518     0.00011159    -0.00079593    -0.00011727     0.00028371     0.00001050     0.00000380     0.00010934
   MO  58     0.00001028     0.00117130     0.00010606     0.00048098    -0.00028179    -0.00021045     0.00009599     0.00014362
   MO  59    -0.00004416    -0.00013342    -0.00006174    -0.00010584     0.00009307    -0.00000225    -0.00001862    -0.00009124
   MO  60     0.00018680     0.00009134    -0.00001209     0.00008945    -0.00029946    -0.00016098    -0.00028743    -0.00028653
   MO  61     0.00002908     0.00008094    -0.00007115     0.00007823    -0.00009402     0.00014557     0.00003069    -0.00007515
   MO  62     0.00006227     0.00011351     0.00021378     0.00008874     0.00002420     0.00010282     0.00003774    -0.00018399
   MO  63    -0.00011399    -0.00011563     0.00009180    -0.00025521     0.00042563    -0.00028539    -0.00017581     0.00065056
   MO  64    -0.00080428    -0.00017186    -0.00039583     0.00024245    -0.00042621     0.00007490     0.00017055    -0.00018362
   MO  65     0.00005760    -0.00002822    -0.00008722    -0.00001395     0.00018322    -0.00061630     0.00048014    -0.00054441
   MO  66    -0.00003885    -0.00001843    -0.00001101     0.00000496    -0.00015283     0.00023951    -0.00053929    -0.00050352
   MO  67    -0.00000417    -0.00000995    -0.00000775    -0.00001912     0.00000033    -0.00001954     0.00051344     0.00014875
   MO  68    -0.00004937    -0.00002877     0.00002248     0.00000921    -0.00007851    -0.00008438    -0.00023635    -0.00020026
   MO  69    -0.00017935    -0.00027614    -0.00021878     0.00174443    -0.00009504    -0.00009416    -0.00003125     0.00008356
   MO  70    -0.00008466    -0.00012978    -0.00017000    -0.00041995    -0.00183650    -0.00026916     0.00044155    -0.00001560
   MO  71    -0.00004218    -0.00012625    -0.00029190    -0.00010515     0.00003533     0.00000016     0.00020350     0.00015654
   MO  72    -0.00001665    -0.00023050    -0.00000703     0.00004775     0.00007982     0.00001040    -0.00004109    -0.00001626
   MO  73     0.00028530     0.00011384    -0.00032871     0.00004691     0.00047023     0.00004563    -0.00011946    -0.00018759
   MO  74    -0.00022837     0.00001547    -0.00005232     0.00043543    -0.00082092    -0.00031178     0.00020245     0.00040387
   MO  75     0.00036899     0.00027255     0.00000508    -0.00114720     0.00061319    -0.00021655    -0.00011085     0.00054020
   MO  76     0.00002705     0.00012637     0.00003938     0.00001919    -0.00017517     0.00108049    -0.00055940     0.00133746
   MO  77     0.00004758     0.00009419    -0.00009193    -0.00003019    -0.00014797     0.00069795    -0.00110107    -0.00050704
   MO  78     0.00006678    -0.00000203     0.00001861    -0.00035700     0.00046899     0.00028440     0.00016654    -0.00073624
   MO  79    -0.00018170    -0.00007914     0.00046044    -0.00009690    -0.00025815     0.00001429    -0.00010470    -0.00013050
   MO  80     0.00013215     0.00003152    -0.00008220    -0.00033292    -0.00055922    -0.00003770     0.00008986    -0.00013079
   MO  81    -0.00008978    -0.00010667     0.00002327    -0.00012873     0.00003082     0.00005161     0.00001319     0.00015035
   MO  82     0.00019553     0.00006343     0.00017682     0.00010961     0.00022006     0.00000051    -0.00003576    -0.00002698

                MO  49         MO  50         MO  51         MO  52         MO  53         MO  54         MO  55         MO  56
   MO   3    -0.00134858     0.00170203    -0.00103821     0.00537218     0.00076219    -0.00253808    -0.00036316    -0.00001666
   MO   4     0.00323784     0.00314216    -0.00528618     0.00620210    -0.01386642     0.00405559    -0.00002885     0.00472008
   MO   5     0.00100744    -0.00074603     0.01203439     0.00053511    -0.00209496    -0.00238881    -0.00025245     0.00055434
   MO   6     0.00143434     0.00328151    -0.00544848     0.00840546    -0.01592733     0.00833546     0.00031460    -0.00675316
   MO   7     0.00523492     0.00680294    -0.00684658    -0.00484379    -0.00864551    -0.01306518    -0.00053455     0.00655379
   MO   8    -0.00208108     0.02694736     0.00683552    -0.00552608    -0.00089343     0.00402252    -0.00023574     0.00049499
   MO   9    -0.00240776     0.00414548    -0.00041311     0.00008633     0.00028285    -0.00033959     0.01820525     0.00016806
   MO  10    -0.00111868    -0.00087178     0.00006498    -0.00011024     0.00021689    -0.00017051     0.00025910     0.00015663
   MO  11    -0.00035487    -0.00080395     0.00073784    -0.00064882     0.00031550     0.00060763    -0.00002522    -0.00117879
   MO  12     0.00060084     0.00075716    -0.00099181    -0.00015995    -0.00163557     0.00035803    -0.00001233     0.00035317
   MO  13     0.00043481     0.00066138     0.00073166    -0.00001890    -0.00113089    -0.00044453    -0.00047985    -0.00007544
   MO  14    -0.00000525    -0.00001516     0.00047177     0.00113852    -0.00110173     0.00141481    -0.00045899    -0.00032376
   MO  15     0.00027333    -0.00087059    -0.00049499     0.00071639    -0.00099908     0.00019662     0.00164376    -0.00037044
   MO  16     0.00008139    -0.00084764    -0.00080737     0.00028144    -0.00119792    -0.00029525    -0.00047975    -0.00051138
   MO  17    -0.00005935    -0.00117099     0.00104183    -0.00166633    -0.00055952     0.00144510    -0.00127614    -0.00075907
   MO  18     0.00090196     0.00048387    -0.00095275    -0.00018925    -0.00067637     0.00021636    -0.00289225     0.00077014
   MO  19    -0.00039820     0.00057011     0.00064613     0.00096012     0.00021795    -0.00040213     0.00232499     0.00028942
   MO  20     0.00030552     0.00126008    -0.00096987    -0.00068834    -0.00214085     0.00034840     0.00381989     0.00041091
   MO  21    -0.00106887    -0.00054147     0.00024796    -0.00072632     0.00286293    -0.00003038     0.00216598     0.00058895
   MO  22     0.00007850    -0.00044720     0.00109520    -0.00037424     0.00108942    -0.00120153    -0.00155386     0.00160404
   MO  23     0.00029513     0.00055644    -0.00212246     0.00233972    -0.00074274     0.00080378    -0.00053130    -0.00085946
   MO  24    -0.00038214    -0.00099533     0.00071741     0.00222445    -0.00035855     0.00343577    -0.00006477     0.00036940
   MO  25    -0.00043185     0.00002762     0.00202756     0.00065327    -0.00015790    -0.00083266     0.00197381    -0.00051874
   MO  26     0.00102186     0.00159862     0.00023213     0.00213573    -0.00321194     0.00155985    -0.00012176    -0.00077504
   MO  27     0.00051205     0.00104899    -0.00196413     0.00061885    -0.00195146     0.00140602    -0.00031171     0.00059116
   MO  28     0.00031633     0.00017809     0.00064401    -0.00059660     0.00041901     0.00022908    -0.00386252    -0.00008174
   MO  29    -0.00077001     0.00030431     0.00039879     0.00083060    -0.00109893     0.00000536    -0.00002907     0.00043986
   MO  30    -0.00086759     0.00064575     0.00021847    -0.00027159     0.00042002     0.00004693    -0.00054180     0.00003216
   MO  31    -0.00009450    -0.00002233    -0.00071251     0.00000894    -0.00022861     0.00023676    -0.00185625    -0.00001561
   MO  32     0.00036518     0.00069865    -0.00096716    -0.00038647     0.00064162     0.00033523     0.00129270     0.00043951
   MO  33    -0.00073777    -0.00025971    -0.00126216     0.00012854     0.00085772     0.00004900    -0.00016391    -0.00027801
   MO  34    -0.00025679     0.00096271     0.00085589    -0.00036970     0.00055369    -0.00009219    -0.00078210     0.00012548
   MO  35     0.00021050    -0.00015061    -0.00051967    -0.00025156    -0.00071737     0.00023311    -0.00011897     0.00002021
   MO  36    -0.00024635     0.00008406     0.00045290    -0.00064446     0.00032289     0.00019543    -0.00038936    -0.00011347
   MO  37    -0.00009005     0.00059530     0.00076562    -0.00027981     0.00016192     0.00007725     0.00024569     0.00005179
   MO  38     0.00042906     0.00072622    -0.00098071     0.00159916    -0.00059079    -0.00082501    -0.00025989     0.00097982
   MO  39    -0.00049987    -0.00021577    -0.00009977     0.00119906     0.00071713     0.00067126    -0.00008642    -0.00039586
   MO  40     0.00003779     0.00033324    -0.00005241    -0.00117548    -0.00054358     0.00061223     0.00021896     0.00072210
   MO  41    -0.00018180    -0.00021626    -0.00040117    -0.00012476     0.00046427     0.00103990     0.00032617    -0.00049289
   MO  42    -0.00015577     0.00025297    -0.00009153    -0.00061721     0.00021345    -0.00007076     0.00054685    -0.00040366
   MO  43     0.00013414    -0.00003969     0.00173949     0.00047459    -0.00131243     0.00009412    -0.00025925    -0.00037627
   MO  44     0.00010688     0.00007621    -0.00005823     0.00105689    -0.00075242     0.00208895     0.00023669     0.00105005
   MO  45    -0.00043239    -0.00014825    -0.00012157    -0.00093716     0.00087236    -0.00048917    -0.00028646    -0.00078257
   MO  46     0.00020035     0.00066118     0.00011148    -0.00003904     0.00010931    -0.00017926     0.00096720    -0.00015553
   MO  47     0.00053203     0.00051610    -0.00039279     0.00008079    -0.00032926     0.00030182    -0.00105547     0.00013209
   MO  48     0.00000329     0.00001645    -0.00006293    -0.00010488     0.00028995    -0.00009239    -0.00013710     0.00026935
   MO  49     0.00244270    -0.00014496    -0.00023303     0.00024097    -0.00048078     0.00022123    -0.00076197    -0.00024649
   MO  50    -0.00014496     0.00327870    -0.00045386     0.00039686    -0.00091911     0.00029219     0.00073054    -0.00011520
   MO  51    -0.00023303    -0.00045386     0.00505484    -0.00016297     0.00030778    -0.00084508     0.00002567    -0.00000703
   MO  52     0.00024097     0.00039686    -0.00016297     0.00319768    -0.00131987     0.00089603     0.00031274    -0.00038386
   MO  53    -0.00048078    -0.00091911     0.00030778    -0.00131987     0.00473466    -0.00119715    -0.00004012     0.00083826
   MO  54     0.00022123     0.00029219    -0.00084508     0.00089603    -0.00119715     0.00303019    -0.00044011    -0.00034257
   MO  55    -0.00076197     0.00073054     0.00002567     0.00031274    -0.00004012    -0.00044011     0.00885074    -0.00002127
   MO  56    -0.00024649    -0.00011520    -0.00000703    -0.00038386     0.00083826    -0.00034257    -0.00002127     0.00274403
   MO  57    -0.00002924     0.00007206    -0.00047638     0.00013311    -0.00042761     0.00037434     0.00001316    -0.00041385
   MO  58    -0.00012239     0.00001518     0.00000557    -0.00011065     0.00014124    -0.00015115     0.00004395    -0.00001280
   MO  59     0.00017825     0.00014578     0.00002850     0.00009909     0.00000203     0.00007836    -0.00000866    -0.00006021
   MO  60     0.00010837    -0.00059325    -0.00007633    -0.00023806    -0.00015665     0.00001645     0.00008207     0.00000728
   MO  61     0.00013031     0.00016898    -0.00006510     0.00020463     0.00005975    -0.00011530    -0.00004818     0.00000313
   MO  62     0.00023491     0.00010565    -0.00079494     0.00003031    -0.00005281     0.00028807     0.00008094    -0.00027119
   MO  63    -0.00063142    -0.00032354    -0.00014563    -0.00008376    -0.00014829    -0.00015085     0.00004416    -0.00010655
   MO  64     0.00024277     0.00024351    -0.00009107    -0.00025922     0.00015256    -0.00044405    -0.00001855     0.00024239
   MO  65    -0.00027486     0.00000203    -0.00005370     0.00006944     0.00002373    -0.00004625     0.00045959     0.00001725
   MO  66    -0.00074883     0.00030900    -0.00005347    -0.00000099     0.00003508     0.00004541    -0.00080381    -0.00001485
   MO  67     0.00022880    -0.00046173    -0.00000072     0.00009096     0.00001132    -0.00010503     0.00060811    -0.00001915
   MO  68     0.00018703    -0.00025446    -0.00005431     0.00018339    -0.00006542    -0.00013341     0.00003734     0.00005059
   MO  69    -0.00011679    -0.00017994     0.00005307     0.00027986    -0.00011046     0.00040781     0.00002736    -0.00023197
   MO  70     0.00026578     0.00064649    -0.00037328     0.00035341    -0.00057862     0.00039522     0.00001405    -0.00032446
   MO  71    -0.00004972     0.00026559    -0.00035137     0.00052116     0.00028559    -0.00015491    -0.00002971     0.00011131
   MO  72    -0.00006035    -0.00009438     0.00001098    -0.00009722    -0.00002574     0.00002898    -0.00001492     0.00020634
   MO  73     0.00004479    -0.00028220    -0.00015855    -0.00040947     0.00030136     0.00007414     0.00002203    -0.00002785
   MO  74    -0.00032910     0.00008064    -0.00019868    -0.00010256    -0.00013424    -0.00001339     0.00000213     0.00088130
   MO  75    -0.00061015    -0.00026305     0.00011408    -0.00035719     0.00024298    -0.00018433     0.00000789    -0.00037351
   MO  76     0.00093561    -0.00003156    -0.00001466    -0.00012984     0.00003427     0.00010125    -0.00019122     0.00001786
   MO  77    -0.00069168     0.00077095    -0.00007399    -0.00002643     0.00007589     0.00004893    -0.00016002    -0.00005085
   MO  78     0.00087190     0.00004109    -0.00017252     0.00004876     0.00004555    -0.00009330     0.00003407    -0.00043952
   MO  79     0.00006220    -0.00010758     0.00138377     0.00012647    -0.00030890    -0.00024928     0.00000985    -0.00008740
   MO  80     0.00014152     0.00002873     0.00016676    -0.00001965     0.00079565     0.00023913     0.00001321     0.00009638
   MO  81    -0.00004886     0.00010714     0.00006251     0.00074713     0.00014687    -0.00028018    -0.00004011    -0.00000302
   MO  82     0.00001296    -0.00000106    -0.00013421     0.00043094    -0.00067538     0.00066745     0.00002003     0.00007250

                MO  57         MO  58         MO  59         MO  60         MO  61         MO  62         MO  63         MO  64
   MO   3    -0.00046536    -0.00012866    -0.00020246     0.00430169     0.00292584     0.00033795    -0.00006139    -0.00005440
   MO   4    -0.00008908     0.00109655     0.00022580     0.00025858    -0.00027159     0.00067332    -0.00306319    -0.00302686
   MO   5     0.00132106     0.00013829     0.00019694    -0.00033022    -0.00052796     0.00366410     0.00129031    -0.00020169
   MO   6     0.00377204    -0.00101975     0.00039354     0.00026392    -0.00013437     0.00109946     0.00013831    -0.00192836
   MO   7    -0.00094945     0.00148876    -0.00044869    -0.00023163    -0.00085665    -0.00130697     0.00028649    -0.00042793
   MO   8    -0.00059110     0.00003323     0.00160348    -0.00823525     0.00165420    -0.00071298    -0.00112732     0.00037817
   MO   9     0.00085674     0.00012987    -0.00014420    -0.00058408     0.00072422     0.00010047    -0.00036927     0.00002196
   MO  10    -0.00006604    -0.00008742    -0.00055644    -0.00009098     0.00007027     0.00003539    -0.00027432     0.00008907
   MO  11     0.00042036    -0.00020910    -0.00005735     0.00026548    -0.00056933    -0.00000658     0.00041495    -0.00016765
   MO  12    -0.00008825     0.00008592    -0.00003181     0.00013796    -0.00044651     0.00017381    -0.00041080     0.00019384
   MO  13     0.00015120    -0.00006563     0.00007514    -0.00003927     0.00000302    -0.00004583    -0.00007103     0.00022581
   MO  14     0.00025711    -0.00001776    -0.00001070    -0.00006986    -0.00007559     0.00007129    -0.00011963    -0.00041372
   MO  15     0.00025449    -0.00004145     0.00004769     0.00052331    -0.00009821     0.00010318     0.00011266    -0.00013149
   MO  16     0.00019085    -0.00012252    -0.00017656     0.00028272    -0.00003693     0.00013968    -0.00003907    -0.00005727
   MO  17     0.00042607    -0.00006638    -0.00015587     0.00104850    -0.00018575     0.00002705     0.00032091    -0.00004006
   MO  18    -0.00059394     0.00019682     0.00002985     0.00016059     0.00047711     0.00021493    -0.00037914     0.00020054
   MO  19     0.00002651    -0.00016790     0.00011525    -0.00026438     0.00068576    -0.00011926     0.00007135     0.00002823
   MO  20     0.00009379     0.00003083    -0.00007723     0.00054536    -0.00022300     0.00023041    -0.00052531     0.00029406
   MO  21    -0.00080376    -0.00001006    -0.00006306    -0.00012338     0.00004872     0.00001949     0.00032613     0.00008487
   MO  22    -0.00090467    -0.00010500     0.00005545    -0.00019808    -0.00014752    -0.00018553    -0.00010196     0.00061585
   MO  23     0.00086871    -0.00004843     0.00011319    -0.00062640     0.00019205    -0.00005718     0.00006747     0.00012708
   MO  24     0.00014170     0.00018885    -0.00007842     0.00005633    -0.00003112     0.00003119    -0.00010157    -0.00054683
   MO  25    -0.00046024     0.00104610    -0.00008093    -0.00032210     0.00000681    -0.00000789     0.00011580     0.00001169
   MO  26    -0.00032925    -0.00069650     0.00018767    -0.00023015     0.00006578     0.00048263    -0.00031376    -0.00001409
   MO  27     0.00015865     0.00089762     0.00020127     0.00021480     0.00002859     0.00019783    -0.00035718     0.00016875
   MO  28    -0.00003016     0.00024193    -0.00001796    -0.00061783     0.00004885    -0.00005589    -0.00009861     0.00000937
   MO  29     0.00035283    -0.00012157    -0.00047439    -0.00043523    -0.00009292    -0.00043584     0.00044023    -0.00013938
   MO  30     0.00002325    -0.00006303     0.00033672    -0.00047565    -0.00011075    -0.00013402    -0.00020190     0.00005129
   MO  31    -0.00002737     0.00028054    -0.00047488     0.00074863     0.00024136     0.00023959    -0.00020084     0.00007265
   MO  32    -0.00041564     0.00019495    -0.00047597    -0.00050752    -0.00011135     0.00027451    -0.00045443     0.00031833
   MO  33    -0.00061360    -0.00012611    -0.00000126    -0.00037510    -0.00019540     0.00048508     0.00080620    -0.00014885
   MO  34    -0.00004084    -0.00001894    -0.00002695    -0.00040035    -0.00062003    -0.00019192    -0.00014161     0.00009314
   MO  35     0.00012939     0.00002619    -0.00038077     0.00033835     0.00035313     0.00014216    -0.00024755     0.00012305
   MO  36    -0.00064196    -0.00006553     0.00000924     0.00034888     0.00055942     0.00000433     0.00034740    -0.00005965
   MO  37     0.00023471    -0.00018222     0.00015551     0.00019558     0.00058896     0.00003100     0.00018990     0.00013017
   MO  38    -0.00012102     0.00017996     0.00012006    -0.00050378     0.00039848    -0.00002798    -0.00059543    -0.00010476
   MO  39     0.00033392    -0.00009987    -0.00007186    -0.00038969     0.00022845    -0.00005771     0.00037721    -0.00033585
   MO  40    -0.00007937    -0.00042815    -0.00002057     0.00052833    -0.00008859     0.00011187    -0.00011803    -0.00008850
   MO  41     0.00058518     0.00001028    -0.00004416     0.00018680     0.00002908     0.00006227    -0.00011399    -0.00080428
   MO  42     0.00011159     0.00117130    -0.00013342     0.00009134     0.00008094     0.00011351    -0.00011563    -0.00017186
   MO  43    -0.00079593     0.00010606    -0.00006174    -0.00001209    -0.00007115     0.00021378     0.00009180    -0.00039583
   MO  44    -0.00011727     0.00048098    -0.00010584     0.00008945     0.00007823     0.00008874    -0.00025521     0.00024245
   MO  45     0.00028371    -0.00028179     0.00009307    -0.00029946    -0.00009402     0.00002420     0.00042563    -0.00042621
   MO  46     0.00001050    -0.00021045    -0.00000225    -0.00016098     0.00014557     0.00010282    -0.00028539     0.00007490
   MO  47     0.00000380     0.00009599    -0.00001862    -0.00028743     0.00003069     0.00003774    -0.00017581     0.00017055
   MO  48     0.00010934     0.00014362    -0.00009124    -0.00028653    -0.00007515    -0.00018399     0.00065056    -0.00018362
   MO  49    -0.00002924    -0.00012239     0.00017825     0.00010837     0.00013031     0.00023491    -0.00063142     0.00024277
   MO  50     0.00007206     0.00001518     0.00014578    -0.00059325     0.00016898     0.00010565    -0.00032354     0.00024351
   MO  51    -0.00047638     0.00000557     0.00002850    -0.00007633    -0.00006510    -0.00079494    -0.00014563    -0.00009107
   MO  52     0.00013311    -0.00011065     0.00009909    -0.00023806     0.00020463     0.00003031    -0.00008376    -0.00025922
   MO  53    -0.00042761     0.00014124     0.00000203    -0.00015665     0.00005975    -0.00005281    -0.00014829     0.00015256
   MO  54     0.00037434    -0.00015115     0.00007836     0.00001645    -0.00011530     0.00028807    -0.00015085    -0.00044405
   MO  55     0.00001316     0.00004395    -0.00000866     0.00008207    -0.00004818     0.00008094     0.00004416    -0.00001855
   MO  56    -0.00041385    -0.00001280    -0.00006021     0.00000728     0.00000313    -0.00027119    -0.00010655     0.00024239
   MO  57     0.00204536    -0.00004690     0.00001709     0.00008376    -0.00000965    -0.00053480    -0.00013831    -0.00014833
   MO  58    -0.00004690     0.00241863    -0.00011579    -0.00001165     0.00002893     0.00002498    -0.00016707    -0.00001821
   MO  59     0.00001709    -0.00011579     0.00139710    -0.00002671    -0.00002742    -0.00004184     0.00002032     0.00004171
   MO  60     0.00008376    -0.00001165    -0.00002671     0.00178696     0.00035595     0.00005920     0.00009596    -0.00007284
   MO  61    -0.00000965     0.00002893    -0.00002742     0.00035595     0.00143482     0.00006246    -0.00000280    -0.00001098
   MO  62    -0.00053480     0.00002498    -0.00004184     0.00005920     0.00006246     0.00129785     0.00013077    -0.00000504
   MO  63    -0.00013831    -0.00016707     0.00002032     0.00009596    -0.00000280     0.00013077     0.00127999    -0.00017703
   MO  64    -0.00014833    -0.00001821     0.00004171    -0.00007284    -0.00001098    -0.00000504    -0.00017703     0.00107987
   MO  65    -0.00000963    -0.00002098     0.00000798    -0.00003122     0.00000754    -0.00001267    -0.00001956     0.00000286
   MO  66     0.00002088     0.00001528    -0.00001843     0.00001295     0.00002097    -0.00000908     0.00000932    -0.00000596
   MO  67    -0.00002844     0.00002570    -0.00001322     0.00005326    -0.00000954     0.00001274     0.00003990    -0.00000655
   MO  68     0.00004529     0.00005281    -0.00002417     0.00014919    -0.00001779    -0.00002447    -0.00009536     0.00000092
   MO  69     0.00017560     0.00002245    -0.00005897    -0.00001592    -0.00001917    -0.00005888     0.00030934     0.00009218
   MO  70     0.00013961    -0.00011151     0.00001704     0.00007496     0.00006437     0.00001881     0.00007469     0.00019080
   MO  71     0.00019607    -0.00004207     0.00003943    -0.00028481    -0.00006149    -0.00002824    -0.00003564     0.00002793
   MO  72    -0.00002202    -0.00035903     0.00001394     0.00004733    -0.00001450    -0.00002747     0.00003692    -0.00000955
   MO  73     0.00040891    -0.00001525     0.00001094     0.00018829     0.00007720    -0.00002809    -0.00010171    -0.00002512
   MO  74    -0.00010290     0.00025983    -0.00009064     0.00004387    -0.00001675    -0.00008674     0.00002399     0.00002609
   MO  75     0.00025131     0.00000975    -0.00004739    -0.00001375    -0.00014452    -0.00011351     0.00040855    -0.00029630
   MO  76     0.00009243    -0.00001039     0.00000431    -0.00001467     0.00004071     0.00003924    -0.00001348    -0.00000781
   MO  77     0.00008774    -0.00003579     0.00000725    -0.00004976     0.00005863    -0.00000113    -0.00009974     0.00000759
   MO  78    -0.00005948    -0.00017438     0.00010354     0.00001439     0.00008079     0.00018786    -0.00025691     0.00002999
   MO  79    -0.00008569    -0.00004925     0.00003099     0.00003890    -0.00001969    -0.00050161    -0.00017082     0.00000718
   MO  80    -0.00024972    -0.00005276     0.00004041    -0.00001465     0.00005704     0.00010826    -0.00016397     0.00010748
   MO  81    -0.00002864    -0.00000853     0.00000260     0.00017988     0.00031761    -0.00004494     0.00010809    -0.00007379
   MO  82     0.00014438     0.00009519    -0.00003658    -0.00000030    -0.00004775     0.00002337    -0.00002850    -0.00011762

                MO  65         MO  66         MO  67         MO  68         MO  69         MO  70         MO  71         MO  72
   MO   3     0.00066518    -0.00065267    -0.00094841    -0.00401774    -0.00001702    -0.00130146     0.00542374    -0.00083959
   MO   4    -0.00001382    -0.00007341    -0.00006387     0.00100762     0.00075992    -0.00122690    -0.00035372    -0.00013384
   MO   5    -0.00031449     0.00018416    -0.00000964     0.00065246    -0.00014837    -0.00044489    -0.00352914     0.00009881
   MO   6    -0.00011264     0.00010501     0.00002484     0.00022787     0.00392641     0.00034723    -0.00034440    -0.00023140
   MO   7    -0.00015959     0.00008570     0.00009038    -0.00023032    -0.00072080     0.00186607     0.00021294    -0.00036066
   MO   8     0.00133711    -0.00124205    -0.00098609    -0.00449302     0.00015672    -0.00093572     0.00137878    -0.00023499
   MO   9    -0.00606633     0.00160220    -0.00250423    -0.00060533    -0.00004366    -0.00001430     0.00015765     0.00008953
   MO  10    -0.00021201     0.00052456     0.00092806    -0.00037191    -0.00006345    -0.00005312     0.00003263     0.00007748
   MO  11    -0.00006811     0.00007336     0.00002804    -0.00008767     0.00066650     0.00023511    -0.00031375    -0.00002404
   MO  12    -0.00000921    -0.00001931    -0.00004338     0.00022573    -0.00043121     0.00005421    -0.00016917     0.00005509
   MO  13    -0.00000191    -0.00002057    -0.00008486     0.00003815    -0.00014469     0.00011911    -0.00006813    -0.00000546
   MO  14    -0.00001788     0.00005877     0.00004327    -0.00001886     0.00035746     0.00003554    -0.00002780     0.00000019
   MO  15    -0.00001518    -0.00008994    -0.00008731     0.00035479     0.00013803     0.00019043    -0.00017186    -0.00000611
   MO  16    -0.00011425     0.00023295     0.00027762     0.00010366    -0.00039853     0.00020690    -0.00007344     0.00001034
   MO  17    -0.00010256     0.00015264    -0.00004158    -0.00009165     0.00098604     0.00010745    -0.00077954     0.00006981
   MO  18    -0.00000580     0.00007856    -0.00022000     0.00003762    -0.00049082     0.00005775     0.00000154     0.00002934
   MO  19     0.00013275    -0.00015615     0.00016061    -0.00015901     0.00020329    -0.00001796     0.00028314    -0.00000882
   MO  20     0.00015855    -0.00023291     0.00017683     0.00022567    -0.00046968     0.00019387    -0.00040398     0.00009289
   MO  21     0.00013090    -0.00008921     0.00015743    -0.00025151     0.00060840    -0.00030451     0.00004909     0.00004024
   MO  22    -0.00014961     0.00013217    -0.00007438     0.00007061    -0.00010629    -0.00039744     0.00006777     0.00013283
   MO  23    -0.00000211     0.00004324     0.00003077    -0.00001607     0.00085794     0.00035710     0.00061767    -0.00015047
   MO  24    -0.00001466     0.00004425     0.00003816    -0.00001556     0.00179526    -0.00037956    -0.00005125     0.00005543
   MO  25     0.00022339    -0.00018242     0.00020669     0.00000517     0.00023609    -0.00021563    -0.00001373    -0.00020639
   MO  26    -0.00003441     0.00001479    -0.00008793    -0.00003243     0.00022980     0.00107769    -0.00001463    -0.00001476
   MO  27    -0.00004438     0.00001275     0.00000632     0.00001154    -0.00004225     0.00097353    -0.00000664    -0.00009605
   MO  28    -0.00081044     0.00054203    -0.00050985    -0.00010250     0.00003981    -0.00019509    -0.00005414    -0.00004055
   MO  29     0.00010485     0.00014277    -0.00010259    -0.00001476     0.00013576     0.00027821     0.00015864     0.00003563
   MO  30     0.00045984     0.00067666    -0.00033509    -0.00007622    -0.00007499    -0.00012989     0.00005346     0.00002355
   MO  31    -0.00015063     0.00072835    -0.00047943     0.00012078     0.00001284     0.00016027    -0.00010346    -0.00003033
   MO  32     0.00032114    -0.00040611     0.00014592    -0.00006685     0.00001545     0.00009229     0.00006433    -0.00002462
   MO  33    -0.00000036     0.00014611    -0.00001185    -0.00017780     0.00020421    -0.00014899     0.00021177    -0.00000650
   MO  34    -0.00010690     0.00005841    -0.00016320    -0.00031339    -0.00008261    -0.00003495     0.00010807    -0.00000955
   MO  35     0.00004520    -0.00001319     0.00000561     0.00010328    -0.00007875     0.00004280    -0.00016669     0.00002766
   MO  36     0.00002696     0.00002686    -0.00007249    -0.00035451     0.00020818     0.00006090    -0.00020204     0.00001705
   MO  37     0.00006284    -0.00007536    -0.00005002    -0.00036375     0.00023364     0.00020433    -0.00000715     0.00000154
   MO  38     0.00006913    -0.00005336    -0.00003110     0.00024920    -0.00116494    -0.00045648     0.00060004     0.00000034
   MO  39     0.00006212    -0.00001034     0.00002715    -0.00024732     0.00108067     0.00009905     0.00036248    -0.00006266
   MO  40    -0.00001145    -0.00003125    -0.00006475    -0.00001382    -0.00050374     0.00002649    -0.00037406     0.00021655
   MO  41     0.00005760    -0.00003885    -0.00000417    -0.00004937    -0.00017935    -0.00008466    -0.00004218    -0.00001665
   MO  42    -0.00002822    -0.00001843    -0.00000995    -0.00002877    -0.00027614    -0.00012978    -0.00012625    -0.00023050
   MO  43    -0.00008722    -0.00001101    -0.00000775     0.00002248    -0.00021878    -0.00017000    -0.00029190    -0.00000703
   MO  44    -0.00001395     0.00000496    -0.00001912     0.00000921     0.00174443    -0.00041995    -0.00010515     0.00004775
   MO  45     0.00018322    -0.00015283     0.00000033    -0.00007851    -0.00009504    -0.00183650     0.00003533     0.00007982
   MO  46    -0.00061630     0.00023951    -0.00001954    -0.00008438    -0.00009416    -0.00026916     0.00000016     0.00001040
   MO  47     0.00048014    -0.00053929     0.00051344    -0.00023635    -0.00003125     0.00044155     0.00020350    -0.00004109
   MO  48    -0.00054441    -0.00050352     0.00014875    -0.00020026     0.00008356    -0.00001560     0.00015654    -0.00001626
   MO  49    -0.00027486    -0.00074883     0.00022880     0.00018703    -0.00011679     0.00026578    -0.00004972    -0.00006035
   MO  50     0.00000203     0.00030900    -0.00046173    -0.00025446    -0.00017994     0.00064649     0.00026559    -0.00009438
   MO  51    -0.00005370    -0.00005347    -0.00000072    -0.00005431     0.00005307    -0.00037328    -0.00035137     0.00001098
   MO  52     0.00006944    -0.00000099     0.00009096     0.00018339     0.00027986     0.00035341     0.00052116    -0.00009722
   MO  53     0.00002373     0.00003508     0.00001132    -0.00006542    -0.00011046    -0.00057862     0.00028559    -0.00002574
   MO  54    -0.00004625     0.00004541    -0.00010503    -0.00013341     0.00040781     0.00039522    -0.00015491     0.00002898
   MO  55     0.00045959    -0.00080381     0.00060811     0.00003734     0.00002736     0.00001405    -0.00002971    -0.00001492
   MO  56     0.00001725    -0.00001485    -0.00001915     0.00005059    -0.00023197    -0.00032446     0.00011131     0.00020634
   MO  57    -0.00000963     0.00002088    -0.00002844     0.00004529     0.00017560     0.00013961     0.00019607    -0.00002202
   MO  58    -0.00002098     0.00001528     0.00002570     0.00005281     0.00002245    -0.00011151    -0.00004207    -0.00035903
   MO  59     0.00000798    -0.00001843    -0.00001322    -0.00002417    -0.00005897     0.00001704     0.00003943     0.00001394
   MO  60    -0.00003122     0.00001295     0.00005326     0.00014919    -0.00001592     0.00007496    -0.00028481     0.00004733
   MO  61     0.00000754     0.00002097    -0.00000954    -0.00001779    -0.00001917     0.00006437    -0.00006149    -0.00001450
   MO  62    -0.00001267    -0.00000908     0.00001274    -0.00002447    -0.00005888     0.00001881    -0.00002824    -0.00002747
   MO  63    -0.00001956     0.00000932     0.00003990    -0.00009536     0.00030934     0.00007469    -0.00003564     0.00003692
   MO  64     0.00000286    -0.00000596    -0.00000655     0.00000092     0.00009218     0.00019080     0.00002793    -0.00000955
   MO  65     0.00077807    -0.00005208     0.00004771    -0.00002845    -0.00002119    -0.00000188     0.00004181     0.00000575
   MO  66    -0.00005208     0.00116739    -0.00039847     0.00000779     0.00002954     0.00001037    -0.00000151     0.00001490
   MO  67     0.00004771    -0.00039847     0.00074787     0.00003501     0.00003153     0.00001406     0.00001749    -0.00000059
   MO  68    -0.00002845     0.00000779     0.00003501     0.00066683    -0.00012361    -0.00019172     0.00001449     0.00000443
   MO  69    -0.00002119     0.00002954     0.00003153    -0.00012361     0.00164311     0.00011599    -0.00000933    -0.00002813
   MO  70    -0.00000188     0.00001037     0.00001406    -0.00019172     0.00011599     0.00189535     0.00006845    -0.00010261
   MO  71     0.00004181    -0.00000151     0.00001749     0.00001449    -0.00000933     0.00006845     0.00056146    -0.00003740
   MO  72     0.00000575     0.00001490    -0.00000059     0.00000443    -0.00002813    -0.00010261    -0.00003740     0.00024013
   MO  73     0.00002768    -0.00000600    -0.00000081     0.00004496     0.00001794    -0.00028783    -0.00004556     0.00002943
   MO  74     0.00000953    -0.00000784    -0.00000198    -0.00000924    -0.00009743     0.00027646     0.00003499     0.00005021
   MO  75     0.00001507     0.00000336     0.00000461    -0.00005761    -0.00037761    -0.00011404    -0.00001839    -0.00001484
   MO  76    -0.00070712    -0.00027321     0.00004508    -0.00000418    -0.00001927     0.00000068    -0.00005121    -0.00002451
   MO  77    -0.00019175     0.00075768    -0.00048249    -0.00000959    -0.00004928    -0.00002202    -0.00000632    -0.00000247
   MO  78     0.00006055    -0.00014188     0.00008340     0.00003148    -0.00017894    -0.00014814    -0.00000871    -0.00003323
   MO  79    -0.00002144    -0.00002168     0.00000851     0.00003550    -0.00006459     0.00002642    -0.00010670     0.00000392
   MO  80     0.00002704    -0.00001864    -0.00000315    -0.00006919    -0.00050995     0.00046326     0.00004416    -0.00002180
   MO  81    -0.00001842     0.00003008     0.00003960     0.00010329    -0.00004876    -0.00001910     0.00014315    -0.00002323
   MO  82     0.00001039    -0.00001650    -0.00001984     0.00008810    -0.00032169    -0.00029852    -0.00003451     0.00001219

                MO  73         MO  74         MO  75         MO  76         MO  77         MO  78         MO  79         MO  80
   MO   3    -0.00281353    -0.00064415     0.00015895     0.00019127     0.00003730    -0.00021616     0.00043446    -0.00006715
   MO   4     0.00057196    -0.00233385     0.00447261     0.00011557     0.00040518     0.00211625    -0.00031968    -0.00532371
   MO   5    -0.00428857    -0.00091856    -0.00013235    -0.00008663     0.00011641    -0.00014632     0.00273048    -0.00003751
   MO   6     0.00027600    -0.00217216     0.00118427     0.00006576     0.00059121     0.00315006    -0.00024555    -0.00429836
   MO   7     0.00002951    -0.00344984    -0.00522824    -0.00009230     0.00043925     0.00283056    -0.00079361    -0.00628945
   MO   8    -0.00127758    -0.00029571    -0.00000622     0.00018996    -0.00000542     0.00002249     0.00016085    -0.00021617
   MO   9    -0.00007021    -0.00003498     0.00002178    -0.00108236     0.00083315    -0.00013778    -0.00006065    -0.00002388
   MO  10    -0.00000673     0.00008473    -0.00006444    -0.00000099     0.00001229    -0.00003827     0.00001250     0.00005444
   MO  11     0.00020645    -0.00046899     0.00038525    -0.00000934    -0.00000118     0.00001270     0.00006308    -0.00011060
   MO  12    -0.00000567     0.00022637    -0.00025551    -0.00001949     0.00000918     0.00011830    -0.00007294    -0.00010031
   MO  13    -0.00009603    -0.00006865    -0.00010380    -0.00002644    -0.00001143     0.00012941     0.00021983    -0.00015110
   MO  14    -0.00007119    -0.00000581    -0.00002783    -0.00003299    -0.00005647    -0.00009343     0.00011234    -0.00008357
   MO  15     0.00005544    -0.00008349    -0.00001885     0.00005471     0.00014424     0.00002172     0.00001758    -0.00014681
   MO  16     0.00000628    -0.00014153     0.00011242    -0.00002161     0.00001486     0.00019727    -0.00001868    -0.00014672
   MO  17     0.00040114    -0.00026375     0.00009945    -0.00005200    -0.00008992    -0.00006884     0.00015490    -0.00021423
   MO  18    -0.00010826     0.00038201    -0.00039456    -0.00006518    -0.00014025     0.00010651    -0.00011214     0.00013003
   MO  19    -0.00016001     0.00003696    -0.00015660     0.00000499     0.00007481    -0.00007812     0.00011741     0.00001821
   MO  20     0.00011024     0.00027496    -0.00035758    -0.00005608     0.00015651     0.00016268    -0.00004069    -0.00011850
   MO  21     0.00009142     0.00012321    -0.00002551    -0.00006086     0.00002465    -0.00028772    -0.00018227     0.00034323
   MO  22    -0.00018695     0.00048403    -0.00046142     0.00006333    -0.00004795    -0.00023157     0.00021454     0.00017546
   MO  23    -0.00009378    -0.00034102    -0.00019052     0.00003974     0.00002190     0.00015464    -0.00039451    -0.00032757
   MO  24     0.00002385     0.00026141    -0.00055595    -0.00002370    -0.00011510    -0.00054367     0.00003811    -0.00005509
   MO  25    -0.00022799    -0.00012674     0.00001369    -0.00025968    -0.00011514    -0.00004942     0.00043852    -0.00010235
   MO  26    -0.00046320    -0.00013374    -0.00049762    -0.00000337     0.00007492     0.00026650     0.00027672    -0.00004394
   MO  27    -0.00018358     0.00068534    -0.00037447     0.00003908    -0.00007083    -0.00025617    -0.00029364     0.00001789
   MO  28    -0.00000083    -0.00008911     0.00000948     0.00073730     0.00044196    -0.00008107     0.00007262     0.00002699
   MO  29    -0.00026116     0.00038877     0.00019528    -0.00019912     0.00003147    -0.00054958     0.00016324    -0.00024695
   MO  30     0.00005285    -0.00007564    -0.00000370    -0.00082765     0.00046505     0.00003104     0.00001915     0.00008718
   MO  31     0.00005651     0.00006358    -0.00017142    -0.00014635     0.00054167    -0.00002546    -0.00012801     0.00001062
   MO  32     0.00002009     0.00017302    -0.00032837    -0.00015432    -0.00023814     0.00014281    -0.00029151     0.00025652
   MO  33    -0.00006694    -0.00008841     0.00029177    -0.00004926     0.00003354    -0.00012002    -0.00046588    -0.00000840
   MO  34    -0.00009391     0.00002716     0.00011644     0.00007127     0.00005120    -0.00009957     0.00012017     0.00009601
   MO  35     0.00010651     0.00004511    -0.00016689    -0.00007053    -0.00002764     0.00010608    -0.00005991    -0.00005791
   MO  36    -0.00003831    -0.00002712     0.00002909    -0.00007155    -0.00003355    -0.00003959     0.00008656     0.00008915
   MO  37     0.00004012    -0.00004894    -0.00001344     0.00004047     0.00000871    -0.00005081     0.00006668     0.00008006
   MO  38    -0.00029397     0.00042256    -0.00029299    -0.00003015     0.00002615     0.00012897    -0.00006566     0.00006218
   MO  39    -0.00007081    -0.00013128     0.00008175    -0.00005036    -0.00007238    -0.00020789    -0.00016085    -0.00015442
   MO  40     0.00013045     0.00045031     0.00000346     0.00007146     0.00001791    -0.00013408    -0.00001232     0.00017134
   MO  41     0.00028530    -0.00022837     0.00036899     0.00002705     0.00004758     0.00006678    -0.00018170     0.00013215
   MO  42     0.00011384     0.00001547     0.00027255     0.00012637     0.00009419    -0.00000203    -0.00007914     0.00003152
   MO  43    -0.00032871    -0.00005232     0.00000508     0.00003938    -0.00009193     0.00001861     0.00046044    -0.00008220
   MO  44     0.00004691     0.00043543    -0.00114720     0.00001919    -0.00003019    -0.00035700    -0.00009690    -0.00033292
   MO  45     0.00047023    -0.00082092     0.00061319    -0.00017517    -0.00014797     0.00046899    -0.00025815    -0.00055922
   MO  46     0.00004563    -0.00031178    -0.00021655     0.00108049     0.00069795     0.00028440     0.00001429    -0.00003770
   MO  47    -0.00011946     0.00020245    -0.00011085    -0.00055940    -0.00110107     0.00016654    -0.00010470     0.00008986
   MO  48    -0.00018759     0.00040387     0.00054020     0.00133746    -0.00050704    -0.00073624    -0.00013050    -0.00013079
   MO  49     0.00004479    -0.00032910    -0.00061015     0.00093561    -0.00069168     0.00087190     0.00006220     0.00014152
   MO  50    -0.00028220     0.00008064    -0.00026305    -0.00003156     0.00077095     0.00004109    -0.00010758     0.00002873
   MO  51    -0.00015855    -0.00019868     0.00011408    -0.00001466    -0.00007399    -0.00017252     0.00138377     0.00016676
   MO  52    -0.00040947    -0.00010256    -0.00035719    -0.00012984    -0.00002643     0.00004876     0.00012647    -0.00001965
   MO  53     0.00030136    -0.00013424     0.00024298     0.00003427     0.00007589     0.00004555    -0.00030890     0.00079565
   MO  54     0.00007414    -0.00001339    -0.00018433     0.00010125     0.00004893    -0.00009330    -0.00024928     0.00023913
   MO  55     0.00002203     0.00000213     0.00000789    -0.00019122    -0.00016002     0.00003407     0.00000985     0.00001321
   MO  56    -0.00002785     0.00088130    -0.00037351     0.00001786    -0.00005085    -0.00043952    -0.00008740     0.00009638
   MO  57     0.00040891    -0.00010290     0.00025131     0.00009243     0.00008774    -0.00005948    -0.00008569    -0.00024972
   MO  58    -0.00001525     0.00025983     0.00000975    -0.00001039    -0.00003579    -0.00017438    -0.00004925    -0.00005276
   MO  59     0.00001094    -0.00009064    -0.00004739     0.00000431     0.00000725     0.00010354     0.00003099     0.00004041
   MO  60     0.00018829     0.00004387    -0.00001375    -0.00001467    -0.00004976     0.00001439     0.00003890    -0.00001465
   MO  61     0.00007720    -0.00001675    -0.00014452     0.00004071     0.00005863     0.00008079    -0.00001969     0.00005704
   MO  62    -0.00002809    -0.00008674    -0.00011351     0.00003924    -0.00000113     0.00018786    -0.00050161     0.00010826
   MO  63    -0.00010171     0.00002399     0.00040855    -0.00001348    -0.00009974    -0.00025691    -0.00017082    -0.00016397
   MO  64    -0.00002512     0.00002609    -0.00029630    -0.00000781     0.00000759     0.00002999     0.00000718     0.00010748
   MO  65     0.00002768     0.00000953     0.00001507    -0.00070712    -0.00019175     0.00006055    -0.00002144     0.00002704
   MO  66    -0.00000600    -0.00000784     0.00000336    -0.00027321     0.00075768    -0.00014188    -0.00002168    -0.00001864
   MO  67    -0.00000081    -0.00000198     0.00000461     0.00004508    -0.00048249     0.00008340     0.00000851    -0.00000315
   MO  68     0.00004496    -0.00000924    -0.00005761    -0.00000418    -0.00000959     0.00003148     0.00003550    -0.00006919
   MO  69     0.00001794    -0.00009743    -0.00037761    -0.00001927    -0.00004928    -0.00017894    -0.00006459    -0.00050995
   MO  70    -0.00028783     0.00027646    -0.00011404     0.00000068    -0.00002202    -0.00014814     0.00002642     0.00046326
   MO  71    -0.00004556     0.00003499    -0.00001839    -0.00005121    -0.00000632    -0.00000871    -0.00010670     0.00004416
   MO  72     0.00002943     0.00005021    -0.00001484    -0.00002451    -0.00000247    -0.00003323     0.00000392    -0.00002180
   MO  73     0.00054681    -0.00016727     0.00007330    -0.00001062     0.00003789     0.00015039    -0.00008229    -0.00008103
   MO  74    -0.00016727     0.00095303    -0.00011343    -0.00001318    -0.00010654    -0.00059686    -0.00005225     0.00012401
   MO  75     0.00007330    -0.00011343     0.00103492    -0.00003558    -0.00007141    -0.00030428    -0.00002103    -0.00002445
   MO  76    -0.00001062    -0.00001318    -0.00003558     0.00201633     0.00000968    -0.00002617    -0.00002759     0.00000060
   MO  77     0.00003789    -0.00010654    -0.00007141     0.00000968     0.00173028    -0.00012083    -0.00004096     0.00001187
   MO  78     0.00015039    -0.00059686    -0.00030428    -0.00002617    -0.00012083     0.00110246    -0.00000964     0.00006210
   MO  79    -0.00008229    -0.00005225    -0.00002103    -0.00002759    -0.00004096    -0.00000964     0.00074855     0.00002452
   MO  80    -0.00008103     0.00012401    -0.00002445     0.00000060     0.00001187     0.00006210     0.00002452     0.00088081
   MO  81    -0.00007215    -0.00002828     0.00000178     0.00000049     0.00000503     0.00000621     0.00004203     0.00000799
   MO  82     0.00006732    -0.00011759     0.00011589     0.00000644     0.00001392     0.00006102    -0.00001006    -0.00007391

                MO  81         MO  82
   MO   3     0.00023353     0.00003512
   MO   4    -0.00008607     0.00153006
   MO   5     0.00046545    -0.00009948
   MO   6    -0.00003942     0.00133183
   MO   7    -0.00040734    -0.00176742
   MO   8    -0.00215112     0.00016353
   MO   9    -0.00004308     0.00001117
   MO  10    -0.00000002    -0.00004641
   MO  11    -0.00027991    -0.00005549
   MO  12    -0.00025589     0.00019869
   MO  13    -0.00002527     0.00002558
   MO  14     0.00002780     0.00022235
   MO  15    -0.00001022     0.00011713
   MO  16     0.00001482     0.00013050
   MO  17    -0.00069847    -0.00000691
   MO  18    -0.00004048     0.00006055
   MO  19     0.00034934    -0.00007742
   MO  20    -0.00050241     0.00022232
   MO  21    -0.00002919    -0.00032592
   MO  22     0.00009885    -0.00025044
   MO  23     0.00050469     0.00008666
   MO  24     0.00003726     0.00029572
   MO  25     0.00020424    -0.00005607
   MO  26     0.00007546     0.00026392
   MO  27    -0.00018061     0.00023146
   MO  28    -0.00010944    -0.00004050
   MO  29     0.00015387     0.00013073
   MO  30    -0.00011544    -0.00003164
   MO  31     0.00004492     0.00001228
   MO  32    -0.00023348    -0.00006293
   MO  33     0.00013053    -0.00004114
   MO  34    -0.00003765    -0.00007156
   MO  35    -0.00013741     0.00009388
   MO  36    -0.00006655    -0.00012087
   MO  37     0.00001856    -0.00008394
   MO  38     0.00057162     0.00026144
   MO  39     0.00031012    -0.00015263
   MO  40    -0.00040711     0.00012685
   MO  41    -0.00008978     0.00019553
   MO  42    -0.00010667     0.00006343
   MO  43     0.00002327     0.00017682
   MO  44    -0.00012873     0.00010961
   MO  45     0.00003082     0.00022006
   MO  46     0.00005161     0.00000051
   MO  47     0.00001319    -0.00003576
   MO  48     0.00015035    -0.00002698
   MO  49    -0.00004886     0.00001296
   MO  50     0.00010714    -0.00000106
   MO  51     0.00006251    -0.00013421
   MO  52     0.00074713     0.00043094
   MO  53     0.00014687    -0.00067538
   MO  54    -0.00028018     0.00066745
   MO  55    -0.00004011     0.00002003
   MO  56    -0.00000302     0.00007250
   MO  57    -0.00002864     0.00014438
   MO  58    -0.00000853     0.00009519
   MO  59     0.00000260    -0.00003658
   MO  60     0.00017988    -0.00000030
   MO  61     0.00031761    -0.00004775
   MO  62    -0.00004494     0.00002337
   MO  63     0.00010809    -0.00002850
   MO  64    -0.00007379    -0.00011762
   MO  65    -0.00001842     0.00001039
   MO  66     0.00003008    -0.00001650
   MO  67     0.00003960    -0.00001984
   MO  68     0.00010329     0.00008810
   MO  69    -0.00004876    -0.00032169
   MO  70    -0.00001910    -0.00029852
   MO  71     0.00014315    -0.00003451
   MO  72    -0.00002323     0.00001219
   MO  73    -0.00007215     0.00006732
   MO  74    -0.00002828    -0.00011759
   MO  75     0.00000178     0.00011589
   MO  76     0.00000049     0.00000644
   MO  77     0.00000503     0.00001392
   MO  78     0.00000621     0.00006102
   MO  79     0.00004203    -0.00001006
   MO  80     0.00000799    -0.00007391
   MO  81     0.00051041     0.00002245
   MO  82     0.00002245     0.00082715

Natural orbital populations,block 1
              MO     1       MO     2       MO     3       MO     4       MO     5       MO     6       MO     7       MO     8
  occ(*)=     2.00000000     2.00000000     1.97575616     1.97372268     1.97156734     1.96357639     1.96345475     1.00668458
              MO     9       MO    10       MO    11       MO    12       MO    13       MO    14       MO    15       MO    16
  occ(*)=     0.98338501     0.01948083     0.01793768     0.01779722     0.01351263     0.01307520     0.01075946     0.00982155
              MO    17       MO    18       MO    19       MO    20       MO    21       MO    22       MO    23       MO    24
  occ(*)=     0.00611207     0.00602324     0.00583660     0.00536593     0.00449151     0.00430530     0.00417789     0.00407733
              MO    25       MO    26       MO    27       MO    28       MO    29       MO    30       MO    31       MO    32
  occ(*)=     0.00318970     0.00210491     0.00175265     0.00130789     0.00120101     0.00094466     0.00075676     0.00075257
              MO    33       MO    34       MO    35       MO    36       MO    37       MO    38       MO    39       MO    40
  occ(*)=     0.00070925     0.00067542     0.00065314     0.00054428     0.00053202     0.00042115     0.00036248     0.00035267
              MO    41       MO    42       MO    43       MO    44       MO    45       MO    46       MO    47       MO    48
  occ(*)=     0.00032964     0.00029348     0.00027424     0.00025513     0.00023871     0.00020947     0.00018534     0.00016157
              MO    49       MO    50       MO    51       MO    52       MO    53       MO    54       MO    55       MO    56
  occ(*)=     0.00013557     0.00012569     0.00007843     0.00006766     0.00006029     0.00004262     0.00004174     0.00003839
              MO    57       MO    58       MO    59       MO    60       MO    61       MO    62       MO    63       MO    64
  occ(*)=     0.00003417     0.00003035     0.00002673     0.00002286     0.00002218     0.00002024     0.00001963     0.00001765
              MO    65       MO    66       MO    67       MO    68       MO    69       MO    70       MO    71       MO    72
  occ(*)=     0.00001633     0.00001376     0.00001219     0.00000975     0.00000776     0.00000542     0.00000511     0.00000458
              MO    73       MO    74       MO    75       MO    76       MO    77       MO    78       MO    79       MO    80
  occ(*)=     0.00000341     0.00000305     0.00000276     0.00000080     0.00000049     0.00000031     0.00000026     0.00000016
              MO    81       MO    82       MO
  occ(*)=     0.00000010     0.00000003


 total number of electrons =   16.0000000000

 test slabel:                     1
  a  


          Mulliken population analysis


  NOTE: For HERMIT use spherical harmonics basis sets !!!
 

                        a   partial gross atomic populations
   ao class       1a         2a         3a         4a         5a         6a  
     2_ s       1.993633   0.000003  -0.000583   0.099549   1.260833   0.003619
     2_ p      -0.000000   0.000004   0.000203   0.736462   0.114886  -0.000701
     2_ d       0.000001   0.000000   0.000033  -0.001805  -0.003441   0.000584
     4_ s      -0.000001   1.998050   1.298748   0.000596   0.000721   0.357669
     4_ p      -0.000006   0.000060   0.050044   0.000738   0.004659   1.156373
     4_ d      -0.000000   0.000001  -0.003326  -0.000036   0.000035   0.006634
     5_ s       0.002536   0.000001   0.000095   0.899787   0.019813   0.000049
     5_ p       0.000781  -0.000000  -0.000161   0.013929  -0.026635   0.000090
     6_ s      -0.000002   0.000417   0.202216   0.000104   0.000622   0.027618
     6_ p      -0.000000   0.000272   0.018225  -0.000007   0.000091   0.037650
     7_ s       0.000000   0.001029   0.370227   0.000769   0.001592   0.361122
     7_ p       0.000000   0.000175   0.034014   0.000039  -0.000013   0.012416
     8_ s       0.002268  -0.000012   0.005726   0.266168   0.571732  -0.000095
     8_ p       0.000789   0.000000   0.000294  -0.042571   0.026673   0.000548
 
   ao class       7a         8a         9a        10a        11a        12a  
     2_ s      -0.000004   0.034440   0.036072   0.000018   0.000000  -0.000001
     2_ p       0.000258   0.441751   0.468892  -0.000027  -0.000001   0.000047
     2_ d       0.000002  -0.000060  -0.000270   0.000003  -0.000001   0.000006
     4_ s       0.015175   0.017726   0.017485  -0.000066   0.000146  -0.000000
     4_ p       1.816860   0.323553   0.287456   0.009297   0.007494   0.014872
     4_ d       0.004791   0.002937   0.003003   0.001535   0.001678   0.000611
     5_ s      -0.000050   0.003189   0.003969  -0.000005  -0.000000  -0.000003
     5_ p       0.000029  -0.000455  -0.001366   0.000003   0.000001  -0.000000
     6_ s       0.001284   0.138894   0.126658   0.000032   0.009761   0.000019
     6_ p       0.053158   0.008213   0.006960   0.000565  -0.000867   0.001148
     7_ s       0.016028   0.019848   0.016448   0.007973  -0.000433  -0.000000
     7_ p       0.055210   0.014014   0.013331   0.000116   0.000158   0.001074
     8_ s      -0.000382   0.003850   0.004896   0.000039   0.000003  -0.000007
     8_ p       0.001095  -0.001216  -0.000150  -0.000000  -0.000002   0.000031
 
   ao class      13a        14a        15a        16a        17a        18a  
     2_ s       0.003741   0.000089   0.000117   0.000038  -0.000000   0.000000
     2_ p       0.001393   0.011711   0.007545   0.000136   0.000005   0.000002
     2_ d       0.000459   0.000019   0.000412   0.000078   0.000001   0.000000
     4_ s       0.000085   0.000020   0.000146   0.005876   0.000000   0.000000
     4_ p       0.000030   0.000034   0.000045   0.002106   0.000441   0.000059
     4_ d       0.000013   0.000006   0.000011   0.000326   0.004900   0.005536
     5_ s       0.003155   0.000078   0.002158  -0.000017  -0.000000  -0.000000
     5_ p      -0.000268   0.000533  -0.000034   0.000004   0.000000   0.000000
     6_ s       0.000002   0.000001  -0.000005  -0.000324  -0.000000  -0.000000
     6_ p       0.000002   0.000004   0.000006   0.000059   0.000079   0.000303
     7_ s       0.000033   0.000006   0.000032   0.001013   0.000000   0.000000
     7_ p      -0.000000   0.000003   0.000009   0.000275   0.000686   0.000123
     8_ s       0.005287   0.000117   0.000162   0.000198  -0.000000   0.000000
     8_ p      -0.000417   0.000455   0.000156   0.000053   0.000001   0.000000
 
   ao class      19a        20a        21a        22a        23a        24a  
     2_ s      -0.000001   0.000004   0.000000   0.000005   0.000194   0.000927
     2_ p       0.000010   0.000033   0.000005   0.001020   0.001243   0.001227
     2_ d       0.000005   0.000014   0.003194   0.001619   0.001215   0.001153
     4_ s       0.000961   0.000014   0.000000   0.000034   0.000053   0.000041
     4_ p       0.000717   0.000322   0.000003   0.000790   0.000687   0.000589
     4_ d       0.003850   0.003665   0.000003   0.000554   0.000569   0.000451
     5_ s      -0.000000   0.000001  -0.000001  -0.000367   0.000094  -0.000559
     5_ p       0.000001   0.000003   0.000672   0.000304   0.000404   0.000019
     6_ s       0.000193   0.000054   0.000000  -0.000007  -0.000003  -0.000009
     6_ p       0.000089   0.000285   0.000000   0.000147   0.000151   0.000115
     7_ s      -0.000091   0.000175  -0.000000  -0.000082  -0.000057  -0.000036
     7_ p       0.000102   0.000808   0.000001   0.000081   0.000064   0.000040
     8_ s       0.000001  -0.000015  -0.000003  -0.000243  -0.000422  -0.000101
     8_ p       0.000000   0.000002   0.000617   0.000450  -0.000014   0.000222
 
   ao class      25a        26a        27a        28a        29a        30a  
     2_ s       0.001925   0.000000   0.000468  -0.000000   0.000000   0.000001
     2_ p       0.001404   0.000024   0.000240   0.002531  -0.000001   0.000002
     2_ d       0.000428   0.001931   0.000899   0.000135  -0.000000   0.000001
     4_ s       0.000018   0.000000  -0.000001   0.000002   0.000077  -0.000004
     4_ p       0.000039   0.000000  -0.000001  -0.000007   0.000234   0.000212
     4_ d       0.000031   0.000000   0.000000   0.000005   0.000593   0.000177
     5_ s      -0.000329  -0.000000  -0.000003  -0.001176  -0.000000  -0.000003
     5_ p      -0.000112   0.000076   0.000084   0.000471   0.000000   0.000001
     6_ s       0.000013  -0.000000  -0.000000   0.000001   0.000089  -0.000087
     6_ p       0.000001  -0.000000  -0.000000  -0.000001   0.000121   0.000015
     7_ s      -0.000003  -0.000000   0.000001  -0.000001   0.000131   0.000353
     7_ p       0.000002   0.000000   0.000000  -0.000000  -0.000043   0.000278
     8_ s      -0.000124  -0.000000  -0.000027  -0.001073  -0.000000   0.000001
     8_ p      -0.000104   0.000074   0.000092   0.000422  -0.000000  -0.000001
 
   ao class      31a        32a        33a        34a        35a        36a  
     2_ s      -0.000000  -0.000001   0.000006   0.000208   0.000000   0.000000
     2_ p       0.000000   0.000001  -0.000002   0.000034   0.000000   0.000001
     2_ d       0.000000   0.000002   0.000011   0.000323   0.000000   0.000001
     4_ s       0.000000   0.000048  -0.000013  -0.000000   0.000000  -0.000032
     4_ p       0.000115   0.000275   0.000041  -0.000000   0.000018   0.000165
     4_ d       0.000362   0.000260   0.000173   0.000002   0.000174   0.000098
     5_ s       0.000000   0.000002  -0.000009  -0.000116  -0.000000  -0.000000
     5_ p       0.000000  -0.000000   0.000003   0.000164   0.000000  -0.000000
     6_ s       0.000000   0.000025   0.000141   0.000007   0.000000   0.000296
     6_ p       0.000058  -0.000028   0.000121   0.000002   0.000263  -0.000003
     7_ s      -0.000000   0.000010   0.000007  -0.000002  -0.000000  -0.000004
     7_ p       0.000222   0.000165   0.000225   0.000005   0.000197   0.000027
     8_ s       0.000000  -0.000002   0.000002  -0.000097   0.000000  -0.000001
     8_ p      -0.000000  -0.000003   0.000003   0.000145  -0.000000  -0.000003
 
   ao class      37a        38a        39a        40a        41a        42a  
     2_ s       0.000001   0.000000   0.000185   0.000073   0.000005   0.000011
     2_ p       0.000000   0.000034   0.000033  -0.000022   0.000002  -0.000000
     2_ d       0.000002   0.000034   0.000085   0.000078   0.000002   0.000015
     4_ s      -0.000000   0.000001   0.000008   0.000010   0.000000   0.000055
     4_ p       0.000236   0.000000   0.000004   0.000006   0.000017   0.000031
     4_ d       0.000032   0.000001   0.000009   0.000011   0.000248   0.000124
     5_ s      -0.000000   0.000000  -0.000005   0.000016  -0.000001   0.000005
     5_ p      -0.000001   0.000178   0.000068   0.000043   0.000001   0.000021
     6_ s       0.000002   0.000000   0.000002   0.000003   0.000000   0.000042
     6_ p       0.000126   0.000000   0.000000   0.000000   0.000006   0.000004
     7_ s       0.000000  -0.000000  -0.000001  -0.000003  -0.000000  -0.000023
     7_ p       0.000134   0.000000   0.000001   0.000002   0.000048   0.000006
     8_ s       0.000001   0.000000  -0.000049   0.000028   0.000001  -0.000012
     8_ p      -0.000000   0.000173   0.000021   0.000108   0.000000   0.000014
 
   ao class      43a        44a        45a        46a        47a        48a  
     2_ s       0.000005   0.000008   0.000000  -0.000000   0.000002   0.000000
     2_ p       0.000148   0.000003   0.000004  -0.000001  -0.000003   0.000000
     2_ d       0.000023   0.000002   0.000000   0.000000   0.000003   0.000033
     4_ s       0.000004   0.000053   0.000000  -0.000014   0.000006   0.000000
     4_ p       0.000003   0.000030   0.000007   0.000005  -0.000039  -0.000000
     4_ d       0.000028   0.000094   0.000156   0.000153   0.000105   0.000000
     5_ s       0.000027  -0.000018   0.000002   0.000000  -0.000006   0.000000
     5_ p      -0.000010   0.000005   0.000000   0.000000   0.000004   0.000062
     6_ s      -0.000001   0.000015   0.000000   0.000045  -0.000006   0.000000
     6_ p       0.000005   0.000029   0.000062   0.000025   0.000011   0.000000
     7_ s      -0.000000   0.000010   0.000000  -0.000001   0.000070   0.000000
     7_ p      -0.000000   0.000010   0.000002  -0.000004   0.000022   0.000000
     8_ s       0.000071   0.000024  -0.000001  -0.000002   0.000013   0.000000
     8_ p      -0.000027  -0.000008   0.000005   0.000002   0.000003   0.000067
 
   ao class      49a        50a        51a        52a        53a        54a  
     2_ s       0.000001   0.000004  -0.000000   0.000000  -0.000008  -0.000003
     2_ p      -0.000000   0.000004  -0.000008   0.000002   0.000004   0.000000
     2_ d       0.000004   0.000083   0.000100   0.000003   0.000006   0.000004
     4_ s      -0.000003   0.000000  -0.000000  -0.000017   0.000000  -0.000014
     4_ p      -0.000000   0.000000   0.000000  -0.000000  -0.000002  -0.000008
     4_ d       0.000045   0.000002   0.000001   0.000004   0.000000   0.000000
     5_ s      -0.000003  -0.000019  -0.000003  -0.000000  -0.000001   0.000002
     5_ p       0.000002   0.000033  -0.000008   0.000002   0.000015   0.000000
     6_ s       0.000010   0.000000  -0.000000  -0.000005  -0.000002   0.000047
     6_ p       0.000054   0.000002   0.000001  -0.000001  -0.000001   0.000013
     7_ s      -0.000006  -0.000000   0.000001   0.000085   0.000014  -0.000004
     7_ p       0.000032   0.000001   0.000002  -0.000005  -0.000000   0.000005
     8_ s       0.000001  -0.000022  -0.000013  -0.000000   0.000024   0.000001
     8_ p       0.000000   0.000037   0.000006   0.000001   0.000010  -0.000001
 
   ao class      55a        56a        57a        58a        59a        60a  
     2_ s      -0.000025   0.000000  -0.000000  -0.000001  -0.000000  -0.000000
     2_ p       0.000001   0.000033  -0.000066   0.000006  -0.000000  -0.000000
     2_ d       0.000038   0.000007   0.000000   0.000047   0.000001   0.000000
     4_ s      -0.000001  -0.000000  -0.000000  -0.000000  -0.000000  -0.000000
     4_ p      -0.000000   0.000000   0.000001  -0.000000  -0.000004   0.000001
     4_ d       0.000000   0.000000  -0.000000  -0.000000   0.000001  -0.000000
     5_ s       0.000017   0.000000   0.000062   0.000002   0.000001   0.000000
     5_ p       0.000001   0.000004  -0.000020  -0.000020   0.000000  -0.000000
     6_ s       0.000005   0.000000  -0.000000  -0.000000  -0.000035  -0.000000
     6_ p       0.000001   0.000001   0.000002   0.000000   0.000035   0.000008
     7_ s      -0.000004  -0.000000  -0.000000  -0.000000  -0.000005   0.000000
     7_ p       0.000004   0.000001   0.000003   0.000001   0.000033   0.000014
     8_ s       0.000002   0.000000   0.000072   0.000000  -0.000000  -0.000000
     8_ p       0.000002  -0.000008  -0.000020  -0.000004   0.000000   0.000000
 
   ao class      61a        62a        63a        64a        65a        66a  
     2_ s      -0.000025  -0.000002   0.000000  -0.000002  -0.000000   0.000001
     2_ p      -0.000002   0.000000   0.000000  -0.000001  -0.000000  -0.000000
     2_ d       0.000002   0.000001   0.000017   0.000001   0.000001   0.000004
     4_ s      -0.000000   0.000000  -0.000000  -0.000000   0.000000  -0.000000
     4_ p       0.000001   0.000008   0.000000   0.000004   0.000002   0.000001
     4_ d      -0.000000  -0.000000  -0.000000  -0.000000  -0.000000   0.000000
     5_ s       0.000043   0.000001   0.000000   0.000001  -0.000000   0.000001
     5_ p      -0.000010  -0.000001  -0.000000   0.000000  -0.000000  -0.000001
     6_ s      -0.000000   0.000000  -0.000000   0.000003   0.000004  -0.000000
     6_ p       0.000002   0.000008   0.000001  -0.000004   0.000003   0.000008
     7_ s       0.000000  -0.000000  -0.000000   0.000000  -0.000002   0.000003
     7_ p       0.000002   0.000003   0.000000   0.000012   0.000008  -0.000004
     8_ s       0.000010   0.000002   0.000000   0.000003   0.000000  -0.000000
     8_ p      -0.000000  -0.000000   0.000000   0.000000   0.000001   0.000001
 
   ao class      67a        68a        69a        70a        71a        72a  
     2_ s       0.000000  -0.000000  -0.000000   0.000000  -0.000000   0.000000
     2_ p      -0.000000  -0.000001   0.000000  -0.000001  -0.000001   0.000001
     2_ d       0.000004  -0.000003  -0.000004  -0.000000  -0.000001  -0.000000
     4_ s      -0.000000  -0.000000   0.000000  -0.000000  -0.000000   0.000000
     4_ p       0.000001   0.000000  -0.000000   0.000000   0.000000   0.000001
     4_ d       0.000000  -0.000000  -0.000000  -0.000000  -0.000000  -0.000000
     5_ s       0.000003  -0.000001   0.000000   0.000001   0.000000  -0.000000
     5_ p      -0.000001   0.000003   0.000003   0.000004   0.000002   0.000003
     6_ s      -0.000000   0.000000   0.000000   0.000000  -0.000000  -0.000001
     6_ p       0.000005   0.000000   0.000000   0.000000   0.000000   0.000000
     7_ s       0.000002   0.000000  -0.000000  -0.000000   0.000000   0.000000
     7_ p      -0.000003  -0.000000  -0.000000  -0.000000   0.000000  -0.000000
     8_ s       0.000002   0.000004  -0.000000  -0.000002  -0.000003   0.000000
     8_ p      -0.000001   0.000007   0.000008   0.000003   0.000008   0.000001
 
   ao class      73a        74a        75a        76a        77a        78a  
     2_ s      -0.000000  -0.000000  -0.000000   0.000000   0.000000   0.000001
     2_ p      -0.000000  -0.000000  -0.000000   0.000002   0.000000   0.000000
     2_ d      -0.000000  -0.000000  -0.000000  -0.000000  -0.000000  -0.000000
     4_ s       0.000006   0.000001   0.000000  -0.000000   0.000000  -0.000000
     4_ p       0.000003   0.000002   0.000001  -0.000000  -0.000000  -0.000000
     4_ d      -0.000000  -0.000000  -0.000000  -0.000000  -0.000000  -0.000000
     5_ s      -0.000000  -0.000001  -0.000000   0.000000  -0.000000  -0.000000
     5_ p       0.000000   0.000003   0.000003   0.000000   0.000000   0.000000
     6_ s      -0.000003   0.000000   0.000000  -0.000000   0.000000   0.000000
     6_ p       0.000001  -0.000000   0.000000  -0.000000  -0.000000  -0.000000
     7_ s      -0.000003  -0.000002  -0.000001   0.000000   0.000000  -0.000000
     7_ p       0.000001   0.000000   0.000000  -0.000000  -0.000000  -0.000000
     8_ s       0.000000   0.000000  -0.000000  -0.000001   0.000000  -0.000000
     8_ p       0.000000  -0.000000   0.000000   0.000000  -0.000000   0.000000
 
   ao class      79a        80a        81a        82a  
     2_ s      -0.000000  -0.000000  -0.000000  -0.000000
     2_ p       0.000000   0.000000  -0.000000  -0.000000
     2_ d      -0.000000  -0.000000  -0.000000  -0.000000
     4_ s      -0.000000  -0.000000  -0.000000  -0.000000
     4_ p      -0.000000  -0.000000  -0.000000  -0.000000
     4_ d      -0.000000  -0.000000  -0.000000  -0.000000
     5_ s      -0.000000  -0.000000  -0.000000   0.000000
     5_ p       0.000000   0.000000  -0.000000  -0.000000
     6_ s       0.000000  -0.000000   0.000000   0.000000
     6_ p      -0.000000  -0.000000  -0.000000  -0.000000
     7_ s      -0.000000   0.000000   0.000000  -0.000000
     7_ p      -0.000000  -0.000000  -0.000000  -0.000000
     8_ s       0.000000   0.000000   0.000000  -0.000000
     8_ p      -0.000000  -0.000000  -0.000000  -0.000000


                        gross atomic populations
     ao            2_         4_         5_         6_         7_         8_
      s         3.435532   3.713670   0.932412   0.508134   0.796216   0.857998
      p         1.790511   3.678624  -0.011068   0.127606   0.134132  -0.011948
      d         0.007543   0.040636   0.000000   0.000000   0.000000   0.000000
    total       5.233586   7.432931   0.921344   0.635740   0.930349   0.846050
 

 Total number of electrons:   16.00000000

 item #                     3 suffix=:.drt1.state3:
 read_civout: repnuc=  -62.1137172157958     
================================================================================
  Reading record                      1  of civout
 INFO:ref#  1vector#  1 method:  0 last record  0max overlap with ref# 95% root-following 0
 MR-CISD energy:   -94.36141904   -32.24770183
 residuum:     0.00004500
 deltae:     0.00000000
================================================================================
  Reading record                      2  of civout
 INFO:ref#  2vector#  2 method:  0 last record  0max overlap with ref# 95% root-following 0
 MR-CISD energy:   -94.21846435   -32.10474713
 residuum:     0.00005590
 deltae:     0.00000000
================================================================================
  Reading record                      3  of civout
 INFO:ref#  3vector#  3 method:  0 last record  1max overlap with ref# 94% root-following 0
 MR-CISD energy:   -94.20726964   -32.09355243
 residuum:     0.00006011
 deltae:     0.00000000
 apxde:     0.00000000

          sovref  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1     0.94706787     0.00025589    -0.00126670    -0.00178264    -0.00337622     0.28275207     0.09678757    -0.06428875
 ref:   2    -0.00022678     0.94640844     0.02607982    -0.02063048     0.00040584     0.01457143     0.00185925     0.29125173
 ref:   3    -0.00089498     0.02479871    -0.93609036    -0.03112810     0.00196800    -0.00286088     0.00760621    -0.00568138

                ci   9         ci  10         ci  11         ci  12         ci  13         ci  14         ci  15         ci  16

                ci  17         ci  18         ci  19         ci  20         ci  21         ci  22         ci  23         ci  24

                ci  25         ci  26         ci  27         ci  28         ci  29         ci  30         ci  31         ci  32

                ci  33         ci  34         ci  35         ci  36         ci  37         ci  38         ci  39         ci  40

          tciref  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1     0.94706787     0.00025589    -0.00126670    -0.00178264    -0.00337622     0.00000000     0.00000000     0.00000000
 ref:   2    -0.00022678     0.94640844     0.02607982    -0.02063048     0.00040584     0.00000000     0.00000000     0.00000000
 ref:   3    -0.00089498     0.02479871    -0.93609036    -0.03112810     0.00196800     0.00000000     0.00000000     0.00000000

                ci   9         ci  10         ci  11         ci  12         ci  13         ci  14         ci  15         ci  16

                ci  17         ci  18         ci  19         ci  20         ci  21         ci  22         ci  23         ci  24

                ci  25         ci  26         ci  27         ci  28         ci  29         ci  30         ci  31         ci  32

                ci  33         ci  34         ci  35         ci  36         ci  37         ci  38         ci  39         ci  40
 computing final density
 =========== Executing IN-CORE method ==========
--------------------------------------------------------------------------------
  1e-density for root #    3
--------------------------------------------------------------------------------
================================================================================
   DYZ=     242  DYX=     503  DYW=     590
   D0Z=     147  D0Y=     928  D0X=     875  D0W=    1184
  DDZI=     189 DDYI=     692 DDXI=     675 DDWI=     866
  DDZE=       0 DDYE=     104 DDXE=     105 DDWE=     140
================================================================================
--------------------------------------------------------------------------------
  2e-density  for root #    3
--------------------------------------------------------------------------------
 call to prpd23 iwx=1
 call to prpd23 iwx=2
 starting prpd24 .... 
d2temp:fourex offset=     0 N= 59550 init,final=   T   Fmode=   1
d2temp:fourex offset= 59550 N= 59292 init,final=   F   Fmode=   1
d2temp:fourex offset=118842 N= 59913 init,final=   F   Fmode=   1
d2temp:fourex offset=178755 N= 59058 init,final=   F   Fmode=   1
d2temp:fourex offset=237813 N= 59667 init,final=   F   Fmode=   1
d2temp:fourex offset=297480 N= 59730 init,final=   F   Fmode=   1
d2temp:fourex offset=357210 N= 59871 init,final=   F   Fmode=   1
d2temp:fourex offset=417081 N= 59544 init,final=   F   Fmode=   1
d2temp:fourex offset=476625 N= 59835 init,final=   F   Fmode=   1
d2temp:fourex offset=536460 N= 59460 init,final=   F   Fmode=   1
d2temp:fourex offset=595920 N= 59343 init,final=   F   Fmode=   1
d2temp:fourex offset=655263 N= 58947 init,final=   F   Fmode=   1
d2temp:fourex offset=714210 N= 58800 init,final=   F   Fmode=   1
d2temp:fourex offset=773010 N= 58800 init,final=   F   Fmode=   1
d2temp:fourex offset=831810 N= 58845 init,final=   F   Fmode=   1
d2temp:fourex offset=890655 N= 59778 init,final=   F   Fmode=   1
d2temp:fourex offset=950433 N= 57120 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58182 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59445 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59955 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58695 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58620 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58140 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56676 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59124 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56793 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59799 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 55743 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59250 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58128 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57915 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 55545 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58545 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 55800 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58632 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57285 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59235 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59400 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58692 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57828 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57120 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59007 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56460 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56628 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59040 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57795 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59907 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58548 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 54660 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59535 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57525 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58692 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59958 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56793 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 54714 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 60000 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56403 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56658 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59787 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59556 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56082 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56250 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59892 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58140 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 55935 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58860 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N=     0 init,final=   F   Tmode=   1
 starting prpd24 .... 
d2temp:fourex offset=     0 N= 59550 init,final=   T   Fmode=   1
d2temp:fourex offset= 59550 N= 59292 init,final=   F   Fmode=   1
d2temp:fourex offset=118842 N= 59913 init,final=   F   Fmode=   1
d2temp:fourex offset=178755 N= 59058 init,final=   F   Fmode=   1
d2temp:fourex offset=237813 N= 59667 init,final=   F   Fmode=   1
d2temp:fourex offset=297480 N= 59730 init,final=   F   Fmode=   1
d2temp:fourex offset=357210 N= 59871 init,final=   F   Fmode=   1
d2temp:fourex offset=417081 N= 59544 init,final=   F   Fmode=   1
d2temp:fourex offset=476625 N= 59835 init,final=   F   Fmode=   1
d2temp:fourex offset=536460 N= 59460 init,final=   F   Fmode=   1
d2temp:fourex offset=595920 N= 59343 init,final=   F   Fmode=   1
d2temp:fourex offset=655263 N= 58947 init,final=   F   Fmode=   1
d2temp:fourex offset=714210 N= 58800 init,final=   F   Fmode=   1
d2temp:fourex offset=773010 N= 58800 init,final=   F   Fmode=   1
d2temp:fourex offset=831810 N= 58845 init,final=   F   Fmode=   1
d2temp:fourex offset=890655 N= 59778 init,final=   F   Fmode=   1
d2temp:fourex offset=950433 N= 57120 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58182 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59445 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59955 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58695 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58620 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58140 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56676 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59124 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56793 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59799 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 55743 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59250 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58128 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57915 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 55545 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58545 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 55800 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58632 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57285 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59235 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59400 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58692 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57828 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57120 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59007 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56460 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56628 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59040 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57795 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59907 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58548 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 54660 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59535 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57525 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58692 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59958 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56793 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 54714 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 60000 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56403 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56658 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59787 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59556 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56082 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56250 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59892 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58140 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 55935 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58860 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N=     0 init,final=   F   Tmode=   1
================================================================================
   DYZ=     242  DYX=    7351  DYW=    8696
   D0Z=     147  D0Y=    8560  D0X=    7763  D0W=   10896
  DDZI=     735 DDYI=    2660 DDXI=    2520 DDWI=    3140
  DDZE=       0 DDYE=       0 DDXE=       0 DDWE=       0
================================================================================
  entering dump_fourext .... 
d2temp:fourex offset=     1 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset= 60001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=120001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=180001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=240001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=300001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=360001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=420001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=480001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=540001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=600001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=660001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=720001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=780001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=840001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=900001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=960001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
Trace of MO density:    12.000000
   12  correlated and     4  frozen core electrons

           D1 (irrep)  block   1

                MO   1         MO   2         MO   3         MO   4         MO   5         MO   6         MO   7         MO   8
   MO   1     2.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
   MO   2     0.00000000     2.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
   MO   3     0.00000000     0.00000000     1.97269361     0.00037110     0.00110886     0.00059561     0.00028612    -0.00330671
   MO   4     0.00000000     0.00000000     0.00037110     1.96800517    -0.00061570    -0.01268510    -0.00415357     0.00238627
   MO   5     0.00000000     0.00000000     0.00110886    -0.00061570     1.96675250    -0.00113185    -0.00059931     0.00013620
   MO   6     0.00000000     0.00000000     0.00059561    -0.01268510    -0.00113185     1.94125946    -0.00818012     0.00614934
   MO   7     0.00000000     0.00000000     0.00028612    -0.00415357    -0.00059931    -0.00818012     1.95255150    -0.00151352
   MO   8     0.00000000     0.00000000    -0.00330671     0.00238627     0.00013620     0.00614934    -0.00151352     1.94256220
   MO   9     0.00000000     0.00000000    -0.00045299     0.00098698    -0.00032861     0.00236975    -0.00039378     0.04346387
   MO  10     0.00000000     0.00000000    -0.00422124     0.00364631     0.00294951     0.00541299     0.00013597    -0.00716805
   MO  11     0.00000000     0.00000000     0.02407380     0.01216974    -0.01676336    -0.02426014     0.00243655     0.01760393
   MO  12     0.00000000     0.00000000     0.02179842    -0.02202340     0.00508357     0.00526497    -0.01123363    -0.01064701
   MO  13     0.00000000     0.00000000     0.00169001    -0.00800556    -0.03149243    -0.00090198    -0.01901389    -0.02104150
   MO  14     0.00000000     0.00000000     0.00259241    -0.01483022    -0.01454827    -0.01322347     0.01823180    -0.00626510
   MO  15     0.00000000     0.00000000     0.00410624    -0.01012620    -0.00748713    -0.00913291    -0.00167043     0.07988425
   MO  16     0.00000000     0.00000000     0.00156626    -0.01224299    -0.00147751    -0.00904628    -0.00867321     0.05017181
   MO  17     0.00000000     0.00000000     0.01400702     0.00345399    -0.02779315    -0.01823372    -0.00613127    -0.00004560
   MO  18     0.00000000     0.00000000    -0.01337974    -0.00796239     0.01727648     0.01791507    -0.00208711    -0.02458262
   MO  19     0.00000000     0.00000000    -0.02522136     0.00385682    -0.00770857     0.00475261    -0.00042721    -0.01337778
   MO  20     0.00000000     0.00000000     0.01494959    -0.02220471    -0.00436219     0.00525097    -0.01532147    -0.01201303
   MO  21     0.00000000     0.00000000    -0.00412660     0.03177885     0.03146502     0.01593737     0.01607583    -0.01386461
   MO  22     0.00000000     0.00000000     0.00334823     0.01306606    -0.00838693     0.03341488    -0.01153218    -0.00222711
   MO  23     0.00000000     0.00000000    -0.00633400    -0.01133579     0.02464726    -0.01821794    -0.01690706     0.01699402
   MO  24     0.00000000     0.00000000     0.00204805    -0.01219911    -0.00221096    -0.00720924     0.03534817     0.00787064
   MO  25     0.00000000     0.00000000    -0.00157233     0.00166033    -0.02488280    -0.00126472    -0.00539096    -0.00690569
   MO  26     0.00000000     0.00000000     0.00009117    -0.01675496    -0.01243305    -0.01211669     0.00365284    -0.01805403
   MO  27     0.00000000     0.00000000     0.00245701    -0.01697727     0.01270027     0.00260111     0.00036369    -0.00584693
   MO  28     0.00000000     0.00000000    -0.00083109     0.00355989    -0.00508210     0.00071104     0.00021081    -0.01760243
   MO  29     0.00000000     0.00000000    -0.00228882    -0.01239636    -0.01122877    -0.00765180    -0.00581253     0.00107213
   MO  30     0.00000000     0.00000000     0.00134785     0.00335171    -0.00026216     0.00252409     0.00369411    -0.02055556
   MO  31     0.00000000     0.00000000     0.00028208    -0.00181662     0.00368648     0.00208107    -0.00208731     0.02552716
   MO  32     0.00000000     0.00000000     0.00397217     0.00682687     0.02085795     0.01381361     0.00226131    -0.02461976
   MO  33     0.00000000     0.00000000    -0.00181041     0.00438033     0.03052040    -0.00571027     0.00053048     0.00227709
   MO  34     0.00000000     0.00000000     0.01327847     0.00574170    -0.00323464     0.00267047     0.00316693    -0.03334697
   MO  35     0.00000000     0.00000000    -0.00346106    -0.00654474    -0.00145832     0.00050006    -0.00727066     0.00586909
   MO  36     0.00000000     0.00000000    -0.01414454     0.00883044     0.00937930     0.00167138     0.00324498    -0.01862654
   MO  37     0.00000000     0.00000000    -0.01412872     0.00716323    -0.00747020     0.00051467    -0.00286812    -0.01826338
   MO  38     0.00000000     0.00000000    -0.01024705    -0.01928734     0.00409495     0.01228863     0.00193013     0.00061608
   MO  39     0.00000000     0.00000000    -0.00766617     0.00562753     0.00539209    -0.01043675     0.01118459    -0.00246259
   MO  40     0.00000000     0.00000000     0.00206496    -0.00645337     0.00017238     0.00552710     0.00713975    -0.00749565
   MO  41     0.00000000     0.00000000    -0.00119017    -0.00185402     0.00473388    -0.01148810     0.02663101     0.00007249
   MO  42     0.00000000     0.00000000    -0.00096404    -0.00170653     0.00200561    -0.00306004     0.00291223    -0.00345007
   MO  43     0.00000000     0.00000000    -0.00035858    -0.00880302    -0.01370155    -0.00488346     0.00793385    -0.00364671
   MO  44     0.00000000     0.00000000     0.00097207    -0.01095277     0.00251209     0.00703428     0.00120668    -0.00574462
   MO  45     0.00000000     0.00000000    -0.00376981    -0.00747790     0.00853947    -0.01499804    -0.00861372    -0.01639253
   MO  46     0.00000000     0.00000000    -0.00056906     0.00044361    -0.00077738    -0.00014620    -0.00286624    -0.00392619
   MO  47     0.00000000     0.00000000    -0.00191344    -0.00026206     0.00304261     0.00250317    -0.00120172    -0.02071043
   MO  48     0.00000000     0.00000000    -0.00311218    -0.00108228     0.00382685    -0.00309560    -0.00148889    -0.00809938
   MO  49     0.00000000     0.00000000     0.00150736     0.00043711    -0.00225745     0.00312850    -0.00129391     0.00106164
   MO  50     0.00000000     0.00000000    -0.00370688    -0.00438871     0.00124627    -0.00208101    -0.00615341    -0.03039737
   MO  51     0.00000000     0.00000000     0.00028730     0.01074001    -0.02706081     0.00717756     0.00947244    -0.00749022
   MO  52     0.00000000     0.00000000    -0.00346276    -0.01271998    -0.00271021    -0.01034915     0.00752318     0.00806058
   MO  53     0.00000000     0.00000000    -0.00013988     0.03065356     0.00545204     0.02457756     0.01569596    -0.00070368
   MO  54     0.00000000     0.00000000     0.00236284    -0.01181372     0.00435248    -0.01306882     0.01918020    -0.00503968
   MO  55     0.00000000     0.00000000     0.00041263    -0.00007234     0.00038647    -0.00074541     0.00075385     0.00051326
   MO  56     0.00000000     0.00000000    -0.00098411     0.00224169     0.00068951     0.02150288     0.00560097    -0.00170304
   MO  57     0.00000000     0.00000000     0.00029679    -0.00550519    -0.00240686    -0.01113078    -0.00368218     0.00055158
   MO  58     0.00000000     0.00000000     0.00007677    -0.00111659     0.00036612     0.00288658     0.00059943    -0.00027709
   MO  59     0.00000000     0.00000000    -0.00114449    -0.00003002    -0.00046327    -0.00064941     0.00113618    -0.00209344
   MO  60     0.00000000     0.00000000     0.00203072    -0.00027744    -0.00102987    -0.00023711     0.00040962     0.02008209
   MO  61     0.00000000     0.00000000    -0.01476948     0.00031008     0.00046020     0.00127659    -0.00150803    -0.00119708
   MO  62     0.00000000     0.00000000     0.00036211    -0.00132700     0.00365306    -0.00118911     0.00046191     0.00288377
   MO  63     0.00000000     0.00000000    -0.00254100     0.00282959     0.00111247    -0.00473547    -0.00379477     0.00145175
   MO  64     0.00000000     0.00000000     0.00191863     0.00344912    -0.00007652     0.00708719    -0.01784423    -0.00019457
   MO  65     0.00000000     0.00000000    -0.00066438    -0.00014923     0.00028990     0.00019179     0.00031275    -0.00101424
   MO  66     0.00000000     0.00000000     0.00015286     0.00028668     0.00001282    -0.00029086     0.00003567     0.00092337
   MO  67     0.00000000     0.00000000    -0.00036911     0.00036012     0.00007842     0.00002141     0.00018016     0.00039753
   MO  68     0.00000000     0.00000000     0.00108475    -0.00286239    -0.00113199     0.00066852    -0.00199743     0.00176012
   MO  69     0.00000000     0.00000000     0.00026332     0.00300120     0.00001223    -0.00571194     0.00601470    -0.00005475
   MO  70     0.00000000     0.00000000    -0.00017248     0.00138240    -0.00010005    -0.00348182     0.00024813     0.00065019
   MO  71     0.00000000     0.00000000     0.00213811    -0.00038915    -0.00031367     0.00035169    -0.00023894    -0.00477995
   MO  72     0.00000000     0.00000000    -0.00021647    -0.00070179    -0.00022199     0.00094624     0.00013937     0.00060444
   MO  73     0.00000000     0.00000000    -0.00158573     0.00204165     0.00033429    -0.00058308     0.00020438     0.00304395
   MO  74     0.00000000     0.00000000    -0.00059117    -0.00613993     0.00022670     0.00693192     0.00087138     0.00048569
   MO  75     0.00000000     0.00000000     0.00003261     0.00682067    -0.00006552    -0.00309175     0.00062259     0.00006515
   MO  76     0.00000000     0.00000000     0.00005302     0.00028705     0.00000474    -0.00005970    -0.00002714    -0.00007980
   MO  77     0.00000000     0.00000000     0.00003069     0.00013886     0.00004042    -0.00068024    -0.00054775     0.00012581
   MO  78     0.00000000     0.00000000    -0.00001090    -0.00023442    -0.00011150    -0.00350517    -0.00313048     0.00048915
   MO  79     0.00000000     0.00000000     0.00021153     0.00057011     0.00179102     0.00098956     0.00059341    -0.00040402
   MO  80     0.00000000     0.00000000     0.00014471     0.00327767    -0.00020062     0.00722462     0.00470132    -0.00049368
   MO  81     0.00000000     0.00000000    -0.00203597    -0.00018475    -0.00013124    -0.00011553     0.00018338     0.00271851
   MO  82     0.00000000     0.00000000     0.00015437    -0.00464156     0.00008660    -0.00512600     0.00006750     0.00023947

                MO   9         MO  10         MO  11         MO  12         MO  13         MO  14         MO  15         MO  16
   MO   3    -0.00045299    -0.00422124     0.02407380     0.02179842     0.00169001     0.00259241     0.00410624     0.00156626
   MO   4     0.00098698     0.00364631     0.01216974    -0.02202340    -0.00800556    -0.01483022    -0.01012620    -0.01224299
   MO   5    -0.00032861     0.00294951    -0.01676336     0.00508357    -0.03149243    -0.01454827    -0.00748713    -0.00147751
   MO   6     0.00236975     0.00541299    -0.02426014     0.00526497    -0.00090198    -0.01322347    -0.00913291    -0.00904628
   MO   7    -0.00039378     0.00013597     0.00243655    -0.01123363    -0.01901389     0.01823180    -0.00167043    -0.00867321
   MO   8     0.04346387    -0.00716805     0.01760393    -0.01064701    -0.02104150    -0.00626510     0.07988425     0.05017181
   MO   9     0.05749122     0.00003575     0.00033911    -0.00031386    -0.00107063    -0.00076676     0.00405985     0.00055228
   MO  10     0.00003575     0.02025596    -0.00016957    -0.00035962    -0.00081635     0.00115488    -0.00376794     0.00466463
   MO  11     0.00033911    -0.00016957     0.00287728    -0.00033103    -0.00015747     0.00053071     0.00129883     0.00035921
   MO  12    -0.00031386    -0.00035962    -0.00033103     0.00231176     0.00071078    -0.00014171    -0.00041411     0.00014650
   MO  13    -0.00107063    -0.00081635    -0.00015747     0.00071078     0.00215660     0.00024422    -0.00102140    -0.00059860
   MO  14    -0.00076676     0.00115488     0.00053071    -0.00014171     0.00024422     0.00146674    -0.00033939    -0.00017755
   MO  15     0.00405985    -0.00376794     0.00129883    -0.00041411    -0.00102140    -0.00033939     0.00689167     0.00310014
   MO  16     0.00055228     0.00466463     0.00035921     0.00014650    -0.00059860    -0.00017755     0.00310014     0.00434646
   MO  17    -0.00180358     0.00035465     0.00286461     0.00010255     0.00068204     0.00083518     0.00031327    -0.00022303
   MO  18    -0.00436226    -0.00077278    -0.00206371     0.00095626     0.00030850    -0.00048144    -0.00169902    -0.00083141
   MO  19     0.00281466     0.00113432    -0.00120602    -0.00125113     0.00037513     0.00026523    -0.00107174    -0.00063948
   MO  20     0.00465721     0.00054189    -0.00032597     0.00255312     0.00115889    -0.00023322    -0.00039910     0.00033007
   MO  21     0.00244781     0.00086525     0.00022230    -0.00127785    -0.00194774    -0.00032361    -0.00154867    -0.00169551
   MO  22    -0.00210260     0.00077340    -0.00056173     0.00012403     0.00058841    -0.00035374    -0.00033466     0.00004511
   MO  23    -0.00028400     0.00015501    -0.00074246    -0.00051201    -0.00123664    -0.00012701     0.00101794     0.00081322
   MO  24    -0.00001800     0.00029606     0.00145875    -0.00138318    -0.00160495     0.00261230     0.00088452    -0.00115014
   MO  25     0.00270648     0.00011886     0.00030922    -0.00057948     0.00116295     0.00062170    -0.00015945    -0.00025592
   MO  26    -0.00055732    -0.00039096    -0.00047520     0.00094728     0.00139689     0.00131501    -0.00059359    -0.00036847
   MO  27    -0.00058029     0.00153276    -0.00106388     0.00151815    -0.00016809     0.00022085    -0.00045720     0.00045912
   MO  28    -0.00589115    -0.00001501     0.00011720    -0.00007383     0.00050756     0.00028583    -0.00154407    -0.00090185
   MO  29    -0.00001798    -0.00165782    -0.00022529     0.00001412     0.00065639     0.00040484     0.00057047    -0.00012828
   MO  30    -0.00139148     0.00117931    -0.00009065     0.00017118     0.00027051     0.00018426    -0.00182053    -0.00083235
   MO  31    -0.00217158    -0.00185470     0.00008014     0.00006053    -0.00042844    -0.00036270     0.00205677     0.00092996
   MO  32     0.00122981    -0.00224555    -0.00068936     0.00084601    -0.00037717    -0.00045373    -0.00177278    -0.00206321
   MO  33    -0.00019370     0.00075872    -0.00015996    -0.00041169    -0.00150674    -0.00042576    -0.00017745     0.00017419
   MO  34    -0.00177409     0.00061155     0.00018271     0.00060820     0.00086099     0.00039890    -0.00340593    -0.00205035
   MO  35     0.00007397     0.00064930    -0.00027828     0.00041839     0.00016943    -0.00016346     0.00044459     0.00080961
   MO  36    -0.00078097     0.00026150    -0.00039879    -0.00057685    -0.00014731    -0.00006606    -0.00174665    -0.00123616
   MO  37    -0.00020622     0.00010181    -0.00034708    -0.00069172     0.00079688     0.00018361    -0.00189700    -0.00145903
   MO  38    -0.00014259     0.00009491    -0.00270042     0.00051796     0.00026236    -0.00029806     0.00016402     0.00087172
   MO  39    -0.00006128    -0.00023344     0.00033050    -0.00162315    -0.00098727     0.00071604    -0.00024106    -0.00097258
   MO  40    -0.00014484    -0.00002217    -0.00010475     0.00107216     0.00041261     0.00005675    -0.00088079    -0.00050324
   MO  41     0.00019600    -0.00029141     0.00067517    -0.00063460    -0.00099615     0.00068608    -0.00015169    -0.00065576
   MO  42     0.00023362    -0.00011896    -0.00001522     0.00013921     0.00000580    -0.00012635    -0.00058010    -0.00028429
   MO  43    -0.00015667    -0.00007065     0.00014865     0.00012737     0.00086032     0.00087663     0.00006829     0.00007798
   MO  44    -0.00013116     0.00000044     0.00019288     0.00004413    -0.00046744     0.00141598    -0.00044074    -0.00148341
   MO  45    -0.00026756    -0.00033133    -0.00006568    -0.00015421     0.00012201    -0.00036116    -0.00220073    -0.00133902
   MO  46    -0.00100494     0.00219148    -0.00019485     0.00003212     0.00018598     0.00007139    -0.00095088     0.00033445
   MO  47     0.00059642     0.00169834    -0.00084675     0.00055489     0.00051696     0.00017761    -0.00285585    -0.00110205
   MO  48    -0.00074997    -0.00030981    -0.00029478    -0.00028629    -0.00004042    -0.00005256    -0.00109489    -0.00083835
   MO  49    -0.00009786    -0.00116522    -0.00022947     0.00047603     0.00032216    -0.00009582     0.00050459     0.00003087
   MO  50    -0.00119505    -0.00086936    -0.00145612     0.00106379     0.00136950     0.00004377    -0.00361685    -0.00227303
   MO  51     0.00006536    -0.00007919     0.00113119    -0.00113991     0.00134863     0.00104101    -0.00083964    -0.00123919
   MO  52     0.00069902    -0.00008565    -0.00058955    -0.00055527    -0.00029376     0.00112628     0.00120266     0.00051206
   MO  53    -0.00016675     0.00002696     0.00068163    -0.00189823    -0.00149086    -0.00066680    -0.00087259    -0.00138328
   MO  54    -0.00080457     0.00006021     0.00072288     0.00011360    -0.00081031     0.00160986    -0.00029041    -0.00092024
   MO  55     0.00800581     0.00006045     0.00007345     0.00001237    -0.00015685    -0.00006450     0.00039108    -0.00011190
   MO  56    -0.00019101     0.00021582    -0.00089221     0.00013270    -0.00018732    -0.00008563    -0.00043047    -0.00054538
   MO  57    -0.00019333    -0.00009488     0.00031097    -0.00003597     0.00031153     0.00017856     0.00020117     0.00014434
   MO  58     0.00003479     0.00038653    -0.00015363     0.00004714    -0.00011612     0.00007034    -0.00013456     0.00000640
   MO  59    -0.00006631     0.00208193    -0.00010509    -0.00009230     0.00000875     0.00012920    -0.00045402     0.00035837
   MO  60     0.00056330     0.00003063     0.00057883     0.00006378    -0.00025784    -0.00014157     0.00144658     0.00092950
   MO  61    -0.00008223    -0.00004444    -0.00105770    -0.00080917     0.00003716    -0.00020100    -0.00023306    -0.00011645
   MO  62     0.00008843    -0.00005728    -0.00002230     0.00018624    -0.00007135     0.00002466     0.00020585     0.00021015
   MO  63     0.00013746     0.00001594     0.00029003    -0.00047385    -0.00012555    -0.00013605     0.00018863     0.00009961
   MO  64    -0.00001360     0.00019984    -0.00033396     0.00062638     0.00064480    -0.00084068    -0.00014771     0.00035297
   MO  65     0.00044097    -0.00037934    -0.00009584    -0.00001473     0.00003013    -0.00001350    -0.00010347    -0.00021564
   MO  66    -0.00094462     0.00091554     0.00008513    -0.00004128    -0.00008088     0.00006060    -0.00000765     0.00035864
   MO  67     0.00055940     0.00150114     0.00002943    -0.00008955    -0.00013098     0.00008859    -0.00013005     0.00046182
   MO  68     0.00000727    -0.00061204    -0.00010615     0.00027051     0.00005941    -0.00012591     0.00071517     0.00035492
   MO  69    -0.00000439    -0.00009566     0.00071112    -0.00058257    -0.00035812     0.00054657     0.00009713    -0.00057326
   MO  70     0.00003298    -0.00011442     0.00032577    -0.00002173     0.00006126     0.00008769     0.00024218     0.00018461
   MO  71    -0.00002664     0.00003875    -0.00028119    -0.00008829     0.00008528     0.00002932    -0.00038023    -0.00021656
   MO  72     0.00005524     0.00004680    -0.00002466     0.00006517     0.00000293    -0.00000574     0.00005257     0.00004189
   MO  73     0.00002377     0.00002229     0.00022949    -0.00011468    -0.00012793    -0.00006361     0.00014665     0.00002760
   MO  74    -0.00001316     0.00005522    -0.00046386     0.00037808     0.00008717     0.00000001     0.00002796     0.00008727
   MO  75     0.00013768    -0.00003153     0.00045851    -0.00038976    -0.00007446    -0.00028660    -0.00010456     0.00006159
   MO  76    -0.00257049    -0.00003012     0.00000096    -0.00001569     0.00001662     0.00000661    -0.00009543     0.00000762
   MO  77    -0.00396100    -0.00003594    -0.00001301     0.00002122     0.00005016     0.00001414    -0.00011627     0.00007214
   MO  78     0.00078421    -0.00001616    -0.00012431     0.00016043     0.00009922    -0.00011390     0.00006134     0.00014585
   MO  79     0.00018984     0.00002194     0.00001356    -0.00006949     0.00004018     0.00003979    -0.00002627    -0.00005338
   MO  80     0.00004885     0.00001421     0.00016205    -0.00027513    -0.00022897     0.00006696    -0.00010267    -0.00022412
   MO  81     0.00006522    -0.00000263    -0.00025186    -0.00026530    -0.00006034     0.00002512     0.00008495     0.00006614
   MO  82    -0.00003445    -0.00002315    -0.00012824     0.00024244    -0.00000355     0.00032619     0.00014628     0.00009770

                MO  17         MO  18         MO  19         MO  20         MO  21         MO  22         MO  23         MO  24
   MO   3     0.01400702    -0.01337974    -0.02522136     0.01494959    -0.00412660     0.00334823    -0.00633400     0.00204805
   MO   4     0.00345399    -0.00796239     0.00385682    -0.02220471     0.03177885     0.01306606    -0.01133579    -0.01219911
   MO   5    -0.02779315     0.01727648    -0.00770857    -0.00436219     0.03146502    -0.00838693     0.02464726    -0.00221096
   MO   6    -0.01823372     0.01791507     0.00475261     0.00525097     0.01593737     0.03341488    -0.01821794    -0.00720924
   MO   7    -0.00613127    -0.00208711    -0.00042721    -0.01532147     0.01607583    -0.01153218    -0.01690706     0.03534817
   MO   8    -0.00004560    -0.02458262    -0.01337778    -0.01201303    -0.01386461    -0.00222711     0.01699402     0.00787064
   MO   9    -0.00180358    -0.00436226     0.00281466     0.00465721     0.00244781    -0.00210260    -0.00028400    -0.00001800
   MO  10     0.00035465    -0.00077278     0.00113432     0.00054189     0.00086525     0.00077340     0.00015501     0.00029606
   MO  11     0.00286461    -0.00206371    -0.00120602    -0.00032597     0.00022230    -0.00056173    -0.00074246     0.00145875
   MO  12     0.00010255     0.00095626    -0.00125113     0.00255312    -0.00127785     0.00012403    -0.00051201    -0.00138318
   MO  13     0.00068204     0.00030850     0.00037513     0.00115889    -0.00194774     0.00058841    -0.00123664    -0.00160495
   MO  14     0.00083518    -0.00048144     0.00026523    -0.00023322    -0.00032361    -0.00035374    -0.00012701     0.00261230
   MO  15     0.00031327    -0.00169902    -0.00107174    -0.00039910    -0.00154867    -0.00033466     0.00101794     0.00088452
   MO  16    -0.00022303    -0.00083141    -0.00063948     0.00033007    -0.00169551     0.00004511     0.00081322    -0.00115014
   MO  17     0.00598400    -0.00077926    -0.00111796     0.00137820    -0.00011960    -0.00080856    -0.00274491     0.00204251
   MO  18    -0.00077926     0.00274593     0.00011623     0.00101222    -0.00014433     0.00027947    -0.00049003    -0.00132891
   MO  19    -0.00111796     0.00011623     0.00210432    -0.00084513     0.00039487     0.00023123     0.00045870     0.00046651
   MO  20     0.00137820     0.00101222    -0.00084513     0.00428663    -0.00157573    -0.00035355    -0.00160455    -0.00194591
   MO  21    -0.00011960    -0.00014433     0.00039487    -0.00157573     0.00416338     0.00006853     0.00003580     0.00219398
   MO  22    -0.00080856     0.00027947     0.00023123    -0.00035355     0.00006853     0.00215544    -0.00078647    -0.00069744
   MO  23    -0.00274491    -0.00049003     0.00045870    -0.00160455     0.00003580    -0.00078647     0.00484877     0.00111675
   MO  24     0.00204251    -0.00132891     0.00046651    -0.00194591     0.00219398    -0.00069744     0.00111675     0.00886206
   MO  25     0.00000490    -0.00080055     0.00084507    -0.00061130    -0.00043790     0.00065352    -0.00093886     0.00024530
   MO  26    -0.00008428     0.00076380     0.00047118     0.00099082    -0.00143985     0.00017757    -0.00037495     0.00099270
   MO  27    -0.00016672     0.00142748    -0.00045332     0.00193380    -0.00069237    -0.00043771     0.00068147     0.00036434
   MO  28     0.00032180     0.00035971    -0.00017964    -0.00077851    -0.00000603     0.00043048    -0.00045548     0.00005523
   MO  29    -0.00063571    -0.00033751     0.00035269    -0.00024130    -0.00108259     0.00012837     0.00045851     0.00035721
   MO  30    -0.00006163     0.00015450     0.00014113     0.00009974     0.00042368     0.00016609    -0.00023593    -0.00012031
   MO  31     0.00055565     0.00048418    -0.00072385    -0.00001835    -0.00061534    -0.00026651     0.00014281     0.00002912
   MO  32    -0.00083995     0.00104455    -0.00037134     0.00080676     0.00143295    -0.00033297     0.00013166    -0.00009931
   MO  33    -0.00144451     0.00005536    -0.00022790    -0.00138626     0.00200006    -0.00004964     0.00130825     0.00036380
   MO  34     0.00009987     0.00015125    -0.00005948     0.00025942     0.00049845     0.00037023    -0.00087900    -0.00028364
   MO  35     0.00065937     0.00055681    -0.00014586     0.00107646    -0.00076018    -0.00035888    -0.00023134    -0.00060211
   MO  36     0.00105922     0.00117683     0.00063756    -0.00017374     0.00155529    -0.00017810    -0.00126112     0.00035466
   MO  37     0.00105383     0.00053431     0.00122390     0.00001318     0.00030810    -0.00014292    -0.00085774     0.00009432
   MO  38    -0.00402171     0.00127120     0.00090598     0.00003750    -0.00161069     0.00062149     0.00152230    -0.00186010
   MO  39    -0.00038539    -0.00084886     0.00093069    -0.00232932     0.00163255    -0.00056494     0.00174134     0.00315994
   MO  40     0.00145290     0.00097352    -0.00049984     0.00195934    -0.00025743    -0.00017945    -0.00192521    -0.00045940
   MO  41     0.00083378    -0.00042297    -0.00015300    -0.00058707     0.00095195    -0.00140208    -0.00004483     0.00213127
   MO  42     0.00034464     0.00030749    -0.00026355     0.00035786    -0.00000967    -0.00063633    -0.00048573    -0.00062742
   MO  43     0.00051375     0.00012626     0.00016808     0.00008977    -0.00068766     0.00043454    -0.00169429     0.00056293
   MO  44     0.00162280     0.00020714     0.00037597     0.00009313     0.00135494     0.00002536     0.00078337     0.00557136
   MO  45    -0.00013500     0.00014817     0.00020141    -0.00030454     0.00079587    -0.00063191     0.00009050    -0.00123366
   MO  46    -0.00019378     0.00010050     0.00023460     0.00005021    -0.00008873     0.00023685     0.00000749    -0.00050283
   MO  47    -0.00046670     0.00084784     0.00056887     0.00092882     0.00019940    -0.00025730     0.00009009    -0.00056061
   MO  48    -0.00057615     0.00007216     0.00031849    -0.00070934     0.00050703     0.00002840     0.00031227    -0.00005011
   MO  49    -0.00015118     0.00040648    -0.00017941     0.00063651    -0.00076631    -0.00002506     0.00012764    -0.00044524
   MO  50    -0.00108777     0.00166211     0.00066721     0.00123629    -0.00053579    -0.00020172     0.00006389    -0.00168975
   MO  51     0.00176465    -0.00119715     0.00090472    -0.00100537    -0.00001033     0.00102943    -0.00266489     0.00137849
   MO  52    -0.00182937    -0.00043815     0.00094317    -0.00132344    -0.00057552    -0.00014124     0.00232530     0.00269759
   MO  53    -0.00031348    -0.00099880     0.00023807    -0.00262273     0.00321166     0.00056022    -0.00027870     0.00088784
   MO  54     0.00162142     0.00000191    -0.00026100     0.00012018     0.00090299    -0.00126561     0.00041236     0.00464131
   MO  55    -0.00017759    -0.00068156     0.00049104     0.00092413     0.00057146    -0.00043736    -0.00016036     0.00012520
   MO  56    -0.00059371     0.00056718     0.00031779     0.00017370     0.00061061     0.00114018    -0.00072302     0.00068899
   MO  57     0.00044864    -0.00060827     0.00008149     0.00026294    -0.00114320    -0.00088572     0.00107749    -0.00001893
   MO  58    -0.00000951     0.00017930    -0.00018271    -0.00004233    -0.00000408    -0.00017261    -0.00000450     0.00029656
   MO  59    -0.00008738    -0.00007919     0.00035531    -0.00001882     0.00007350     0.00032154     0.00016545    -0.00001858
   MO  60     0.00148992    -0.00003589    -0.00048555     0.00058939    -0.00031902    -0.00025361    -0.00067586     0.00016906
   MO  61    -0.00031600     0.00085587     0.00114106    -0.00038503     0.00005019    -0.00026851     0.00021257    -0.00027132
   MO  62    -0.00004678     0.00026231    -0.00017708     0.00016609     0.00003900    -0.00009105    -0.00013010    -0.00004203
   MO  63     0.00022514    -0.00032692     0.00016275    -0.00062007     0.00042626     0.00003709     0.00009877    -0.00003580
   MO  64    -0.00025652     0.00035447    -0.00013057     0.00085689    -0.00060554     0.00072610     0.00020968    -0.00216057
   MO  65    -0.00005642     0.00007979     0.00011482     0.00010192     0.00008542    -0.00014613    -0.00002054    -0.00000944
   MO  66     0.00007509    -0.00007411    -0.00005924    -0.00013012     0.00000757     0.00011964     0.00006221     0.00007182
   MO  67    -0.00005105    -0.00020148     0.00012889    -0.00000555     0.00011699     0.00003739     0.00013443     0.00010001
   MO  68    -0.00036089    -0.00003174    -0.00023660     0.00018205    -0.00049625     0.00013447     0.00023681    -0.00031876
   MO  69     0.00099990    -0.00056512     0.00021139    -0.00072385     0.00091614    -0.00016124     0.00067111     0.00247853
   MO  70     0.00026590    -0.00012872    -0.00006552     0.00008761    -0.00029308    -0.00047470     0.00040973    -0.00009692
   MO  71    -0.00061705     0.00000681     0.00029075    -0.00025071    -0.00003513     0.00005922     0.00045254    -0.00008259
   MO  72     0.00004622     0.00002514    -0.00000946     0.00011399    -0.00000094     0.00010360    -0.00012360     0.00001836
   MO  73     0.00039519    -0.00017660    -0.00010063     0.00000276     0.00010958    -0.00013256    -0.00010617     0.00010162
   MO  74    -0.00026621     0.00039306    -0.00001325     0.00046579    -0.00027985     0.00025725    -0.00013741    -0.00003653
   MO  75     0.00013857    -0.00048445    -0.00013240    -0.00047598     0.00017032    -0.00030395    -0.00027178    -0.00097882
   MO  76     0.00005814     0.00015945    -0.00016089    -0.00030851    -0.00017798     0.00015191     0.00004919    -0.00002275
   MO  77     0.00009294     0.00030512    -0.00025650    -0.00036410    -0.00025101     0.00018481     0.00006197    -0.00008788
   MO  78    -0.00022787     0.00013510     0.00001084     0.00032377    -0.00020435    -0.00008303     0.00002220    -0.00053114
   MO  79     0.00004587    -0.00006626     0.00008538    -0.00005245     0.00001108     0.00009565    -0.00018626     0.00007194
   MO  80     0.00007284    -0.00012161     0.00002502    -0.00034241     0.00048349    -0.00002867    -0.00023665     0.00041237
   MO  81    -0.00062456    -0.00004926     0.00031450    -0.00048242    -0.00003066     0.00007806     0.00047954     0.00008353
   MO  82    -0.00002773     0.00013385    -0.00004855     0.00026192    -0.00031913    -0.00022115     0.00011623     0.00061661

                MO  25         MO  26         MO  27         MO  28         MO  29         MO  30         MO  31         MO  32
   MO   3    -0.00157233     0.00009117     0.00245701    -0.00083109    -0.00228882     0.00134785     0.00028208     0.00397217
   MO   4     0.00166033    -0.01675496    -0.01697727     0.00355989    -0.01239636     0.00335171    -0.00181662     0.00682687
   MO   5    -0.02488280    -0.01243305     0.01270027    -0.00508210    -0.01122877    -0.00026216     0.00368648     0.02085795
   MO   6    -0.00126472    -0.01211669     0.00260111     0.00071104    -0.00765180     0.00252409     0.00208107     0.01381361
   MO   7    -0.00539096     0.00365284     0.00036369     0.00021081    -0.00581253     0.00369411    -0.00208731     0.00226131
   MO   8    -0.00690569    -0.01805403    -0.00584693    -0.01760243     0.00107213    -0.02055556     0.02552716    -0.02461976
   MO   9     0.00270648    -0.00055732    -0.00058029    -0.00589115    -0.00001798    -0.00139148    -0.00217158     0.00122981
   MO  10     0.00011886    -0.00039096     0.00153276    -0.00001501    -0.00165782     0.00117931    -0.00185470    -0.00224555
   MO  11     0.00030922    -0.00047520    -0.00106388     0.00011720    -0.00022529    -0.00009065     0.00008014    -0.00068936
   MO  12    -0.00057948     0.00094728     0.00151815    -0.00007383     0.00001412     0.00017118     0.00006053     0.00084601
   MO  13     0.00116295     0.00139689    -0.00016809     0.00050756     0.00065639     0.00027051    -0.00042844    -0.00037717
   MO  14     0.00062170     0.00131501     0.00022085     0.00028583     0.00040484     0.00018426    -0.00036270    -0.00045373
   MO  15    -0.00015945    -0.00059359    -0.00045720    -0.00154407     0.00057047    -0.00182053     0.00205677    -0.00177278
   MO  16    -0.00025592    -0.00036847     0.00045912    -0.00090185    -0.00012828    -0.00083235     0.00092996    -0.00206321
   MO  17     0.00000490    -0.00008428    -0.00016672     0.00032180    -0.00063571    -0.00006163     0.00055565    -0.00083995
   MO  18    -0.00080055     0.00076380     0.00142748     0.00035971    -0.00033751     0.00015450     0.00048418     0.00104455
   MO  19     0.00084507     0.00047118    -0.00045332    -0.00017964     0.00035269     0.00014113    -0.00072385    -0.00037134
   MO  20    -0.00061130     0.00099082     0.00193380    -0.00077851    -0.00024130     0.00009974    -0.00001835     0.00080676
   MO  21    -0.00043790    -0.00143985    -0.00069237    -0.00000603    -0.00108259     0.00042368    -0.00061534     0.00143295
   MO  22     0.00065352     0.00017757    -0.00043771     0.00043048     0.00012837     0.00016609    -0.00026651    -0.00033297
   MO  23    -0.00093886    -0.00037495     0.00068147    -0.00045548     0.00045851    -0.00023593     0.00014281     0.00013166
   MO  24     0.00024530     0.00099270     0.00036434     0.00005523     0.00035721    -0.00012031     0.00002912    -0.00009931
   MO  25     0.00312416     0.00112062    -0.00093824     0.00029038     0.00055755     0.00001623    -0.00060491    -0.00065232
   MO  26     0.00112062     0.00400025     0.00063288     0.00031826     0.00071156     0.00012556    -0.00036636     0.00003603
   MO  27    -0.00093824     0.00063288     0.00310678    -0.00029606    -0.00055244     0.00024361     0.00008320     0.00031643
   MO  28     0.00029038     0.00031826    -0.00029606     0.00198466     0.00013939     0.00054881    -0.00035026     0.00021989
   MO  29     0.00055755     0.00071156    -0.00055244     0.00013939     0.00189570    -0.00028857     0.00012143     0.00009981
   MO  30     0.00001623     0.00012556     0.00024361     0.00054881    -0.00028857     0.00145213    -0.00094356     0.00029110
   MO  31    -0.00060491    -0.00036636     0.00008320    -0.00035026     0.00012143    -0.00094356     0.00233429    -0.00015314
   MO  32    -0.00065232     0.00003603     0.00031643     0.00021989     0.00009981     0.00029110    -0.00015314     0.00321856
   MO  33    -0.00040467    -0.00053318    -0.00019637    -0.00003439    -0.00043999     0.00001928    -0.00017978     0.00048645
   MO  34     0.00042254     0.00052768    -0.00018927     0.00093885    -0.00006160     0.00095935    -0.00109717     0.00100603
   MO  35    -0.00039962     0.00004239     0.00047352    -0.00027173     0.00013214    -0.00041618     0.00079300     0.00020017
   MO  36     0.00010962     0.00024683    -0.00014552     0.00024092    -0.00067113     0.00015315    -0.00003768     0.00035607
   MO  37     0.00025651     0.00024058    -0.00022390     0.00026097    -0.00020820     0.00034364    -0.00041015    -0.00007645
   MO  38    -0.00010903     0.00065139     0.00097457    -0.00041912     0.00075373    -0.00016933     0.00009652    -0.00008885
   MO  39     0.00016870    -0.00024454    -0.00077973     0.00012177     0.00033345    -0.00000501    -0.00023698     0.00018299
   MO  40    -0.00099156     0.00034518     0.00082850    -0.00005952    -0.00031310     0.00018396     0.00005478     0.00045151
   MO  41    -0.00080543    -0.00064992    -0.00001063    -0.00002602    -0.00032406     0.00010591     0.00006152     0.00032734
   MO  42     0.00036941    -0.00058570     0.00043965     0.00025026    -0.00037413     0.00005167     0.00017183     0.00034050
   MO  43     0.00175103     0.00200884    -0.00043031     0.00036389     0.00046578    -0.00020613    -0.00011963    -0.00063552
   MO  44     0.00002086     0.00106807     0.00134304     0.00012773     0.00019007     0.00006445     0.00007569     0.00080252
   MO  45    -0.00016693    -0.00086829    -0.00088785     0.00031045    -0.00024843     0.00054122    -0.00086450     0.00058774
   MO  46    -0.00005764     0.00029933    -0.00000840     0.00069903    -0.00034833     0.00007415    -0.00018078    -0.00026040
   MO  47    -0.00000477     0.00057115     0.00088635    -0.00020755    -0.00018487     0.00060228    -0.00110252     0.00120998
   MO  48     0.00002221    -0.00042940    -0.00021100     0.00054497     0.00035575    -0.00043537    -0.00063750     0.00013108
   MO  49    -0.00023012     0.00073476     0.00043510    -0.00003164    -0.00039907    -0.00056284     0.00003683     0.00040773
   MO  50     0.00003292     0.00173195     0.00097333     0.00083413     0.00016888     0.00106363    -0.00063409     0.00144100
   MO  51     0.00254614     0.00085280    -0.00211794     0.00092396     0.00060415     0.00035173    -0.00090051    -0.00109944
   MO  52     0.00068236     0.00172029     0.00032639    -0.00057948     0.00089233    -0.00045666     0.00019298    -0.00074433
   MO  53    -0.00030961    -0.00308251    -0.00196558     0.00040019    -0.00112201     0.00041461    -0.00025370     0.00067900
   MO  54    -0.00077802     0.00129574     0.00126277     0.00022065    -0.00010981     0.00025342    -0.00013425     0.00058527
   MO  55     0.00051237    -0.00003319    -0.00005788    -0.00125670    -0.00004042    -0.00012967    -0.00064877     0.00046031
   MO  56    -0.00062083    -0.00065635     0.00052399    -0.00007935     0.00033222     0.00011024    -0.00012483     0.00033160
   MO  57    -0.00071788    -0.00064984     0.00030995    -0.00009482     0.00041686     0.00007105     0.00000251    -0.00038991
   MO  58     0.00111245    -0.00067604     0.00105046     0.00027404    -0.00029887     0.00002508     0.00018940     0.00002170
   MO  59    -0.00005547     0.00007713     0.00089238    -0.00006679    -0.00118676     0.00084079    -0.00123025    -0.00144872
   MO  60    -0.00050299    -0.00053082     0.00027753    -0.00099257    -0.00064468    -0.00080257     0.00126363    -0.00105143
   MO  61    -0.00001752    -0.00000366    -0.00002522    -0.00005967    -0.00004418    -0.00025710     0.00042171    -0.00018349
   MO  62     0.00005508     0.00051774     0.00018239    -0.00006086    -0.00038889    -0.00020610     0.00034552     0.00028552
   MO  63     0.00020428    -0.00030915    -0.00047001    -0.00008587     0.00020758    -0.00017331    -0.00020827    -0.00060906
   MO  64     0.00000575    -0.00009944     0.00018929    -0.00003113    -0.00016254     0.00005792     0.00006146     0.00019801
   MO  65     0.00012646    -0.00001346    -0.00003326    -0.00056052     0.00010861     0.00037321    -0.00007098     0.00032146
   MO  66    -0.00006964    -0.00002215     0.00003069     0.00026225     0.00003119     0.00050437     0.00041998    -0.00040968
   MO  67     0.00009357    -0.00011009     0.00007597    -0.00029080    -0.00016131    -0.00016270    -0.00041641    -0.00011513
   MO  68    -0.00000604    -0.00000841     0.00007087    -0.00011322     0.00022187    -0.00016679     0.00019192    -0.00012957
   MO  69     0.00022408     0.00016564    -0.00014036     0.00009463     0.00007814    -0.00003093    -0.00002270     0.00007642
   MO  70    -0.00010719     0.00063861     0.00028446    -0.00008776    -0.00002695    -0.00005473     0.00011352    -0.00002764
   MO  71     0.00005968     0.00006678    -0.00002784     0.00001678     0.00015790     0.00012254    -0.00019560     0.00008037
   MO  72    -0.00014055     0.00000360    -0.00000312    -0.00005128     0.00002908     0.00000572    -0.00000898    -0.00002749
   MO  73    -0.00024203    -0.00049688    -0.00011964    -0.00001672    -0.00016038     0.00002474     0.00007473    -0.00004810
   MO  74    -0.00013174     0.00006639     0.00061436    -0.00011211     0.00029710    -0.00007977     0.00009422     0.00006572
   MO  75     0.00008137    -0.00080353    -0.00075905     0.00006959    -0.00000693     0.00003605    -0.00020135    -0.00030928
   MO  76    -0.00031562    -0.00000319     0.00004072     0.00080531    -0.00016119    -0.00060710    -0.00001433    -0.00018814
   MO  77    -0.00035253     0.00008173     0.00001727     0.00081223     0.00003370     0.00047024     0.00068754    -0.00032976
   MO  78    -0.00003239     0.00029259     0.00000751    -0.00018909    -0.00038669    -0.00000855    -0.00002317     0.00022227
   MO  79     0.00023519     0.00011134    -0.00019479     0.00002731     0.00008697     0.00002879    -0.00011411    -0.00014010
   MO  80    -0.00002871    -0.00010973    -0.00025686     0.00007471    -0.00030690     0.00010670    -0.00005017     0.00015111
   MO  81     0.00016061     0.00004430    -0.00015509    -0.00016212     0.00011865    -0.00016523     0.00012655    -0.00027106
   MO  82    -0.00008729     0.00043818     0.00048884    -0.00007268     0.00025034    -0.00005728     0.00003443    -0.00003472

                MO  33         MO  34         MO  35         MO  36         MO  37         MO  38         MO  39         MO  40
   MO   3    -0.00181041     0.01327847    -0.00346106    -0.01414454    -0.01412872    -0.01024705    -0.00766617     0.00206496
   MO   4     0.00438033     0.00574170    -0.00654474     0.00883044     0.00716323    -0.01928734     0.00562753    -0.00645337
   MO   5     0.03052040    -0.00323464    -0.00145832     0.00937930    -0.00747020     0.00409495     0.00539209     0.00017238
   MO   6    -0.00571027     0.00267047     0.00050006     0.00167138     0.00051467     0.01228863    -0.01043675     0.00552710
   MO   7     0.00053048     0.00316693    -0.00727066     0.00324498    -0.00286812     0.00193013     0.01118459     0.00713975
   MO   8     0.00227709    -0.03334697     0.00586909    -0.01862654    -0.01826338     0.00061608    -0.00246259    -0.00749565
   MO   9    -0.00019370    -0.00177409     0.00007397    -0.00078097    -0.00020622    -0.00014259    -0.00006128    -0.00014484
   MO  10     0.00075872     0.00061155     0.00064930     0.00026150     0.00010181     0.00009491    -0.00023344    -0.00002217
   MO  11    -0.00015996     0.00018271    -0.00027828    -0.00039879    -0.00034708    -0.00270042     0.00033050    -0.00010475
   MO  12    -0.00041169     0.00060820     0.00041839    -0.00057685    -0.00069172     0.00051796    -0.00162315     0.00107216
   MO  13    -0.00150674     0.00086099     0.00016943    -0.00014731     0.00079688     0.00026236    -0.00098727     0.00041261
   MO  14    -0.00042576     0.00039890    -0.00016346    -0.00006606     0.00018361    -0.00029806     0.00071604     0.00005675
   MO  15    -0.00017745    -0.00340593     0.00044459    -0.00174665    -0.00189700     0.00016402    -0.00024106    -0.00088079
   MO  16     0.00017419    -0.00205035     0.00080961    -0.00123616    -0.00145903     0.00087172    -0.00097258    -0.00050324
   MO  17    -0.00144451     0.00009987     0.00065937     0.00105922     0.00105383    -0.00402171    -0.00038539     0.00145290
   MO  18     0.00005536     0.00015125     0.00055681     0.00117683     0.00053431     0.00127120    -0.00084886     0.00097352
   MO  19    -0.00022790    -0.00005948    -0.00014586     0.00063756     0.00122390     0.00090598     0.00093069    -0.00049984
   MO  20    -0.00138626     0.00025942     0.00107646    -0.00017374     0.00001318     0.00003750    -0.00232932     0.00195934
   MO  21     0.00200006     0.00049845    -0.00076018     0.00155529     0.00030810    -0.00161069     0.00163255    -0.00025743
   MO  22    -0.00004964     0.00037023    -0.00035888    -0.00017810    -0.00014292     0.00062149    -0.00056494    -0.00017945
   MO  23     0.00130825    -0.00087900    -0.00023134    -0.00126112    -0.00085774     0.00152230     0.00174134    -0.00192521
   MO  24     0.00036380    -0.00028364    -0.00060211     0.00035466     0.00009432    -0.00186010     0.00315994    -0.00045940
   MO  25    -0.00040467     0.00042254    -0.00039962     0.00010962     0.00025651    -0.00010903     0.00016870    -0.00099156
   MO  26    -0.00053318     0.00052768     0.00004239     0.00024683     0.00024058     0.00065139    -0.00024454     0.00034518
   MO  27    -0.00019637    -0.00018927     0.00047352    -0.00014552    -0.00022390     0.00097457    -0.00077973     0.00082850
   MO  28    -0.00003439     0.00093885    -0.00027173     0.00024092     0.00026097    -0.00041912     0.00012177    -0.00005952
   MO  29    -0.00043999    -0.00006160     0.00013214    -0.00067113    -0.00020820     0.00075373     0.00033345    -0.00031310
   MO  30     0.00001928     0.00095935    -0.00041618     0.00015315     0.00034364    -0.00016933    -0.00000501     0.00018396
   MO  31    -0.00017978    -0.00109717     0.00079300    -0.00003768    -0.00041015     0.00009652    -0.00023698     0.00005478
   MO  32     0.00048645     0.00100603     0.00020017     0.00035607    -0.00007645    -0.00008885     0.00018299     0.00045151
   MO  33     0.00274683    -0.00018452    -0.00059563     0.00053682    -0.00069121    -0.00009750     0.00092762    -0.00072863
   MO  34    -0.00018452     0.00284890    -0.00065103     0.00039502     0.00076490    -0.00049887     0.00002396     0.00051067
   MO  35    -0.00059563    -0.00065103     0.00106717     0.00005602     0.00004594     0.00012618    -0.00059468     0.00048002
   MO  36     0.00053682     0.00039502     0.00005602     0.00233299     0.00128469    -0.00104633     0.00030250     0.00065271
   MO  37    -0.00069121     0.00076490     0.00004594     0.00128469     0.00198875    -0.00061253     0.00035897     0.00058905
   MO  38    -0.00009750    -0.00049887     0.00012618    -0.00104633    -0.00061253     0.00403066    -0.00052283    -0.00049139
   MO  39     0.00092762     0.00002396    -0.00059468     0.00030250     0.00035897    -0.00052283     0.00242739    -0.00114859
   MO  40    -0.00072863     0.00051067     0.00048002     0.00065271     0.00058905    -0.00049139    -0.00114859     0.00192532
   MO  41     0.00006824     0.00008999    -0.00011152     0.00027090     0.00021298    -0.00079147     0.00100376     0.00039202
   MO  42    -0.00012348     0.00029075     0.00017713     0.00031195     0.00019877    -0.00020492    -0.00030004     0.00016183
   MO  43    -0.00023240     0.00022245    -0.00010695     0.00059423     0.00002567    -0.00004308    -0.00038304     0.00019162
   MO  44    -0.00007034     0.00021465     0.00000663     0.00057303     0.00054887    -0.00118342     0.00160083     0.00021821
   MO  45     0.00092508     0.00127377    -0.00040593     0.00087652     0.00080290    -0.00049398     0.00035154     0.00017721
   MO  46    -0.00001482     0.00050331     0.00005230     0.00012148     0.00020827     0.00017265    -0.00017484     0.00003124
   MO  47    -0.00022764     0.00143001     0.00008486     0.00049224     0.00090447     0.00036716    -0.00004424     0.00058134
   MO  48     0.00064789     0.00074038    -0.00043421     0.00038804     0.00036378     0.00006416     0.00053683    -0.00007686
   MO  49    -0.00076618    -0.00028456     0.00027950    -0.00040683    -0.00013108     0.00042908    -0.00049804     0.00009495
   MO  50    -0.00042807     0.00229443    -0.00008767     0.00076935     0.00134213     0.00093656    -0.00040698     0.00084036
   MO  51    -0.00168276     0.00125012    -0.00063956     0.00063261     0.00112365    -0.00141621     0.00024378    -0.00003333
   MO  52     0.00025022    -0.00066405    -0.00036726    -0.00080491    -0.00054421     0.00168394     0.00142078    -0.00131589
   MO  53     0.00100368     0.00057322    -0.00092276     0.00038198     0.00022033    -0.00121884     0.00145482    -0.00077402
   MO  54     0.00020139     0.00020671    -0.00007997     0.00051222     0.00019403    -0.00120773     0.00127177     0.00069909
   MO  55    -0.00002493    -0.00023318     0.00000063    -0.00003767     0.00002220    -0.00015308     0.00002725     0.00005858
   MO  56    -0.00023407     0.00016827    -0.00008192    -0.00007399     0.00015676     0.00080939    -0.00015444     0.00071580
   MO  57    -0.00106303    -0.00003207     0.00025677    -0.00099234     0.00035922     0.00007462     0.00026851    -0.00003051
   MO  58    -0.00007038     0.00000778     0.00000952    -0.00004115    -0.00016934     0.00013943    -0.00006747    -0.00043453
   MO  59     0.00019648     0.00002471    -0.00069169     0.00008804     0.00030731     0.00014263    -0.00021719     0.00000100
   MO  60    -0.00048196    -0.00075325     0.00047169     0.00032065     0.00012016    -0.00062966    -0.00052396     0.00060722
   MO  61    -0.00024361    -0.00084171     0.00051231     0.00092923     0.00094738     0.00057942     0.00029710    -0.00002744
   MO  62     0.00067567    -0.00022118     0.00013895     0.00000343     0.00002782     0.00000181    -0.00007208     0.00013254
   MO  63     0.00093865    -0.00029558    -0.00025473     0.00038703     0.00015938    -0.00056590     0.00035341    -0.00024023
   MO  64    -0.00026044     0.00008451     0.00030760    -0.00027555     0.00004792     0.00023108    -0.00099382    -0.00011575
   MO  65    -0.00003563    -0.00000805     0.00007338     0.00008405     0.00013093     0.00005501     0.00007798     0.00003439
   MO  66     0.00016681    -0.00007399    -0.00008469    -0.00003772    -0.00011685    -0.00004269    -0.00002058    -0.00007147
   MO  67     0.00006011    -0.00021404    -0.00010195    -0.00010654    -0.00010249     0.00001430     0.00003364    -0.00012351
   MO  68    -0.00017055    -0.00058135     0.00015665    -0.00067250    -0.00065900     0.00046172    -0.00031170    -0.00020245
   MO  69     0.00025425    -0.00002694    -0.00019095     0.00028886     0.00025019    -0.00128296     0.00132081    -0.00047276
   MO  70    -0.00013502    -0.00014160     0.00009230    -0.00005483     0.00006014    -0.00034165     0.00009313    -0.00011092
   MO  71     0.00004803     0.00024116    -0.00016373    -0.00021496     0.00006825     0.00053988     0.00028747    -0.00029918
   MO  72    -0.00002757    -0.00003716     0.00002653     0.00000393    -0.00000989     0.00000822    -0.00006625     0.00019518
   MO  73    -0.00014852    -0.00009497     0.00007662    -0.00007566     0.00011737    -0.00032597    -0.00000009     0.00015449
   MO  74    -0.00015448    -0.00005665     0.00014342    -0.00012220    -0.00008325     0.00053753    -0.00027216     0.00036337
   MO  75     0.00034195     0.00015575    -0.00020496     0.00009605     0.00004555    -0.00047449     0.00002679    -0.00012193
   MO  76    -0.00003002     0.00012337    -0.00004436    -0.00001591     0.00001550    -0.00001082    -0.00003126     0.00003858
   MO  77     0.00004211     0.00013585    -0.00000012     0.00002313    -0.00002874     0.00006958    -0.00006781     0.00001108
   MO  78    -0.00018805    -0.00009918     0.00011913    -0.00008068    -0.00005079     0.00028347    -0.00029256    -0.00000077
   MO  79    -0.00028335     0.00008994    -0.00006030     0.00010825     0.00001809    -0.00004115    -0.00005530    -0.00002973
   MO  80     0.00006372     0.00008916    -0.00014048     0.00013641     0.00007956    -0.00024888     0.00016342     0.00000777
   MO  81     0.00012555    -0.00007967    -0.00011370    -0.00006590     0.00000783     0.00048847     0.00027573    -0.00033414
   MO  82    -0.00007452    -0.00008292     0.00008676    -0.00012274    -0.00010348     0.00034125    -0.00006182     0.00020965

                MO  41         MO  42         MO  43         MO  44         MO  45         MO  46         MO  47         MO  48
   MO   3    -0.00119017    -0.00096404    -0.00035858     0.00097207    -0.00376981    -0.00056906    -0.00191344    -0.00311218
   MO   4    -0.00185402    -0.00170653    -0.00880302    -0.01095277    -0.00747790     0.00044361    -0.00026206    -0.00108228
   MO   5     0.00473388     0.00200561    -0.01370155     0.00251209     0.00853947    -0.00077738     0.00304261     0.00382685
   MO   6    -0.01148810    -0.00306004    -0.00488346     0.00703428    -0.01499804    -0.00014620     0.00250317    -0.00309560
   MO   7     0.02663101     0.00291223     0.00793385     0.00120668    -0.00861372    -0.00286624    -0.00120172    -0.00148889
   MO   8     0.00007249    -0.00345007    -0.00364671    -0.00574462    -0.01639253    -0.00392619    -0.02071043    -0.00809938
   MO   9     0.00019600     0.00023362    -0.00015667    -0.00013116    -0.00026756    -0.00100494     0.00059642    -0.00074997
   MO  10    -0.00029141    -0.00011896    -0.00007065     0.00000044    -0.00033133     0.00219148     0.00169834    -0.00030981
   MO  11     0.00067517    -0.00001522     0.00014865     0.00019288    -0.00006568    -0.00019485    -0.00084675    -0.00029478
   MO  12    -0.00063460     0.00013921     0.00012737     0.00004413    -0.00015421     0.00003212     0.00055489    -0.00028629
   MO  13    -0.00099615     0.00000580     0.00086032    -0.00046744     0.00012201     0.00018598     0.00051696    -0.00004042
   MO  14     0.00068608    -0.00012635     0.00087663     0.00141598    -0.00036116     0.00007139     0.00017761    -0.00005256
   MO  15    -0.00015169    -0.00058010     0.00006829    -0.00044074    -0.00220073    -0.00095088    -0.00285585    -0.00109489
   MO  16    -0.00065576    -0.00028429     0.00007798    -0.00148341    -0.00133902     0.00033445    -0.00110205    -0.00083835
   MO  17     0.00083378     0.00034464     0.00051375     0.00162280    -0.00013500    -0.00019378    -0.00046670    -0.00057615
   MO  18    -0.00042297     0.00030749     0.00012626     0.00020714     0.00014817     0.00010050     0.00084784     0.00007216
   MO  19    -0.00015300    -0.00026355     0.00016808     0.00037597     0.00020141     0.00023460     0.00056887     0.00031849
   MO  20    -0.00058707     0.00035786     0.00008977     0.00009313    -0.00030454     0.00005021     0.00092882    -0.00070934
   MO  21     0.00095195    -0.00000967    -0.00068766     0.00135494     0.00079587    -0.00008873     0.00019940     0.00050703
   MO  22    -0.00140208    -0.00063633     0.00043454     0.00002536    -0.00063191     0.00023685    -0.00025730     0.00002840
   MO  23    -0.00004483    -0.00048573    -0.00169429     0.00078337     0.00009050     0.00000749     0.00009009     0.00031227
   MO  24     0.00213127    -0.00062742     0.00056293     0.00557136    -0.00123366    -0.00050283    -0.00056061    -0.00005011
   MO  25    -0.00080543     0.00036941     0.00175103     0.00002086    -0.00016693    -0.00005764    -0.00000477     0.00002221
   MO  26    -0.00064992    -0.00058570     0.00200884     0.00106807    -0.00086829     0.00029933     0.00057115    -0.00042940
   MO  27    -0.00001063     0.00043965    -0.00043031     0.00134304    -0.00088785    -0.00000840     0.00088635    -0.00021100
   MO  28    -0.00002602     0.00025026     0.00036389     0.00012773     0.00031045     0.00069903    -0.00020755     0.00054497
   MO  29    -0.00032406    -0.00037413     0.00046578     0.00019007    -0.00024843    -0.00034833    -0.00018487     0.00035575
   MO  30     0.00010591     0.00005167    -0.00020613     0.00006445     0.00054122     0.00007415     0.00060228    -0.00043537
   MO  31     0.00006152     0.00017183    -0.00011963     0.00007569    -0.00086450    -0.00018078    -0.00110252    -0.00063750
   MO  32     0.00032734     0.00034050    -0.00063552     0.00080252     0.00058774    -0.00026040     0.00120998     0.00013108
   MO  33     0.00006824    -0.00012348    -0.00023240    -0.00007034     0.00092508    -0.00001482    -0.00022764     0.00064789
   MO  34     0.00008999     0.00029075     0.00022245     0.00021465     0.00127377     0.00050331     0.00143001     0.00074038
   MO  35    -0.00011152     0.00017713    -0.00010695     0.00000663    -0.00040593     0.00005230     0.00008486    -0.00043421
   MO  36     0.00027090     0.00031195     0.00059423     0.00057303     0.00087652     0.00012148     0.00049224     0.00038804
   MO  37     0.00021298     0.00019877     0.00002567     0.00054887     0.00080290     0.00020827     0.00090447     0.00036378
   MO  38    -0.00079147    -0.00020492    -0.00004308    -0.00118342    -0.00049398     0.00017265     0.00036716     0.00006416
   MO  39     0.00100376    -0.00030004    -0.00038304     0.00160083     0.00035154    -0.00017484    -0.00004424     0.00053683
   MO  40     0.00039202     0.00016183     0.00019162     0.00021821     0.00017721     0.00003124     0.00058134    -0.00007686
   MO  41     0.00203590     0.00040662    -0.00037024     0.00054059     0.00037821    -0.00020091     0.00003552     0.00005929
   MO  42     0.00040662     0.00109701     0.00001240    -0.00033330     0.00046308     0.00001691     0.00032116     0.00018375
   MO  43    -0.00037024     0.00001240     0.00234471     0.00006336    -0.00036903     0.00009159    -0.00025333    -0.00002185
   MO  44     0.00054059    -0.00033330     0.00006336     0.00519718    -0.00079750    -0.00022685     0.00042095    -0.00001023
   MO  45     0.00037821     0.00046308    -0.00036903    -0.00079750     0.00301027     0.00015145     0.00099788     0.00103380
   MO  46    -0.00020091     0.00001691     0.00009159    -0.00022685     0.00015145     0.00113762     0.00002167     0.00026771
   MO  47     0.00003552     0.00032116    -0.00025333     0.00042095     0.00099788     0.00002167     0.00257972     0.00039556
   MO  48     0.00005929     0.00018375    -0.00002185    -0.00001023     0.00103380     0.00026771     0.00039556     0.00194930
   MO  49    -0.00018974    -0.00008318    -0.00002621    -0.00003356    -0.00073193     0.00022500     0.00024345    -0.00020493
   MO  50    -0.00036674     0.00048886     0.00010350     0.00014629     0.00137906     0.00069233     0.00214644     0.00050322
   MO  51    -0.00004199    -0.00008934     0.00195665     0.00032443    -0.00017020     0.00014414    -0.00024396     0.00003741
   MO  52     0.00011175    -0.00072133     0.00059314     0.00133079    -0.00111877    -0.00013653    -0.00020815    -0.00014230
   MO  53     0.00114192     0.00017118    -0.00136669    -0.00028744     0.00087884     0.00001465    -0.00042139     0.00032540
   MO  54     0.00179982    -0.00005980     0.00028934     0.00310125    -0.00037746    -0.00016573     0.00029125    -0.00005134
   MO  55     0.00011733     0.00002995    -0.00001095     0.00009290     0.00008427    -0.00039615     0.00045321    -0.00018047
   MO  56    -0.00006215    -0.00034889    -0.00035942     0.00101721    -0.00063811    -0.00016852     0.00016393     0.00025946
   MO  57     0.00058689     0.00010269    -0.00123698    -0.00009004     0.00016003    -0.00005907     0.00019608     0.00003621
   MO  58     0.00007238     0.00114138     0.00013784     0.00049979    -0.00029305    -0.00016887     0.00013743     0.00012687
   MO  59    -0.00017601    -0.00017055    -0.00014985    -0.00006168     0.00007233     0.00025816     0.00018022    -0.00011486
   MO  60     0.00022963     0.00009306    -0.00005155     0.00007180    -0.00052115    -0.00020013    -0.00054728    -0.00042356
   MO  61     0.00003833     0.00014464    -0.00007050     0.00000268     0.00006654     0.00013666     0.00014143     0.00000057
   MO  62    -0.00001725     0.00011739     0.00033868     0.00007378    -0.00000171     0.00006280     0.00001625    -0.00018426
   MO  63    -0.00021034    -0.00015091     0.00017033    -0.00021412     0.00057509    -0.00020302    -0.00029361     0.00059393
   MO  64    -0.00143692    -0.00013130    -0.00054036    -0.00058211    -0.00005675     0.00022892     0.00025989    -0.00022242
   MO  65     0.00008575    -0.00000347    -0.00007497     0.00001789     0.00020614    -0.00055489     0.00047490    -0.00046589
   MO  66    -0.00005496    -0.00003361    -0.00001716    -0.00000806    -0.00017051     0.00024454    -0.00050744    -0.00044685
   MO  67    -0.00003613    -0.00007634    -0.00003035    -0.00003219    -0.00003131    -0.00001380     0.00044259     0.00009314
   MO  68    -0.00023832    -0.00012964    -0.00005906    -0.00026744    -0.00044954    -0.00016108    -0.00033834    -0.00024461
   MO  69     0.00021583    -0.00028269    -0.00012497     0.00211768    -0.00018102    -0.00014823    -0.00010498     0.00007744
   MO  70     0.00007371    -0.00007160    -0.00010227    -0.00021912    -0.00068501    -0.00001467     0.00015707    -0.00029148
   MO  71    -0.00003167    -0.00010257    -0.00026565    -0.00008209     0.00010161     0.00004515     0.00024932     0.00012445
   MO  72    -0.00002059    -0.00021235    -0.00000083     0.00002222     0.00000347    -0.00000146    -0.00005048    -0.00001024
   MO  73     0.00031892     0.00009246    -0.00045055     0.00004079     0.00014372    -0.00003230    -0.00010142    -0.00012471
   MO  74    -0.00021970     0.00002304     0.00000592     0.00032864    -0.00041562    -0.00016951     0.00015869     0.00021074
   MO  75     0.00026883     0.00027339    -0.00009845    -0.00140336     0.00085781    -0.00012614    -0.00018812     0.00051152
   MO  76     0.00000287     0.00005583     0.00004076    -0.00001351    -0.00010590     0.00076937    -0.00035661     0.00105373
   MO  77    -0.00001151     0.00001091    -0.00005672    -0.00005561    -0.00016044     0.00047222    -0.00082145    -0.00045255
   MO  78    -0.00005172    -0.00001370     0.00000123    -0.00033979    -0.00003707     0.00021500     0.00022676    -0.00062407
   MO  79    -0.00005351    -0.00005131     0.00020186    -0.00004155    -0.00015439     0.00000937    -0.00006231    -0.00007876
   MO  80     0.00042081     0.00003349    -0.00000545    -0.00010637    -0.00013567     0.00002241    -0.00005748    -0.00017778
   MO  81    -0.00004932    -0.00009078     0.00001818    -0.00008362     0.00000687     0.00003343    -0.00000387     0.00012503
   MO  82     0.00025194     0.00002588     0.00023532     0.00035987    -0.00014709    -0.00009248     0.00003988     0.00003518

                MO  49         MO  50         MO  51         MO  52         MO  53         MO  54         MO  55         MO  56
   MO   3     0.00150736    -0.00370688     0.00028730    -0.00346276    -0.00013988     0.00236284     0.00041263    -0.00098411
   MO   4     0.00043711    -0.00438871     0.01074001    -0.01271998     0.03065356    -0.01181372    -0.00007234     0.00224169
   MO   5    -0.00225745     0.00124627    -0.02706081    -0.00271021     0.00545204     0.00435248     0.00038647     0.00068951
   MO   6     0.00312850    -0.00208101     0.00717756    -0.01034915     0.02457756    -0.01306882    -0.00074541     0.02150288
   MO   7    -0.00129391    -0.00615341     0.00947244     0.00752318     0.01569596     0.01918020     0.00075385     0.00560097
   MO   8     0.00106164    -0.03039737    -0.00749022     0.00806058    -0.00070368    -0.00503968     0.00051326    -0.00170304
   MO   9    -0.00009786    -0.00119505     0.00006536     0.00069902    -0.00016675    -0.00080457     0.00800581    -0.00019101
   MO  10    -0.00116522    -0.00086936    -0.00007919    -0.00008565     0.00002696     0.00006021     0.00006045     0.00021582
   MO  11    -0.00022947    -0.00145612     0.00113119    -0.00058955     0.00068163     0.00072288     0.00007345    -0.00089221
   MO  12     0.00047603     0.00106379    -0.00113991    -0.00055527    -0.00189823     0.00011360     0.00001237     0.00013270
   MO  13     0.00032216     0.00136950     0.00134863    -0.00029376    -0.00149086    -0.00081031    -0.00015685    -0.00018732
   MO  14    -0.00009582     0.00004377     0.00104101     0.00112628    -0.00066680     0.00160986    -0.00006450    -0.00008563
   MO  15     0.00050459    -0.00361685    -0.00083964     0.00120266    -0.00087259    -0.00029041     0.00039108    -0.00043047
   MO  16     0.00003087    -0.00227303    -0.00123919     0.00051206    -0.00138328    -0.00092024    -0.00011190    -0.00054538
   MO  17    -0.00015118    -0.00108777     0.00176465    -0.00182937    -0.00031348     0.00162142    -0.00017759    -0.00059371
   MO  18     0.00040648     0.00166211    -0.00119715    -0.00043815    -0.00099880     0.00000191    -0.00068156     0.00056718
   MO  19    -0.00017941     0.00066721     0.00090472     0.00094317     0.00023807    -0.00026100     0.00049104     0.00031779
   MO  20     0.00063651     0.00123629    -0.00100537    -0.00132344    -0.00262273     0.00012018     0.00092413     0.00017370
   MO  21    -0.00076631    -0.00053579    -0.00001033    -0.00057552     0.00321166     0.00090299     0.00057146     0.00061061
   MO  22    -0.00002506    -0.00020172     0.00102943    -0.00014124     0.00056022    -0.00126561    -0.00043736     0.00114018
   MO  23     0.00012764     0.00006389    -0.00266489     0.00232530    -0.00027870     0.00041236    -0.00016036    -0.00072302
   MO  24    -0.00044524    -0.00168975     0.00137849     0.00269759     0.00088784     0.00464131     0.00012520     0.00068899
   MO  25    -0.00023012     0.00003292     0.00254614     0.00068236    -0.00030961    -0.00077802     0.00051237    -0.00062083
   MO  26     0.00073476     0.00173195     0.00085280     0.00172029    -0.00308251     0.00129574    -0.00003319    -0.00065635
   MO  27     0.00043510     0.00097333    -0.00211794     0.00032639    -0.00196558     0.00126277    -0.00005788     0.00052399
   MO  28    -0.00003164     0.00083413     0.00092396    -0.00057948     0.00040019     0.00022065    -0.00125670    -0.00007935
   MO  29    -0.00039907     0.00016888     0.00060415     0.00089233    -0.00112201    -0.00010981    -0.00004042     0.00033222
   MO  30    -0.00056284     0.00106363     0.00035173    -0.00045666     0.00041461     0.00025342    -0.00012967     0.00011024
   MO  31     0.00003683    -0.00063409    -0.00090051     0.00019298    -0.00025370    -0.00013425    -0.00064877    -0.00012483
   MO  32     0.00040773     0.00144100    -0.00109944    -0.00074433     0.00067900     0.00058527     0.00046031     0.00033160
   MO  33    -0.00076618    -0.00042807    -0.00168276     0.00025022     0.00100368     0.00020139    -0.00002493    -0.00023407
   MO  34    -0.00028456     0.00229443     0.00125012    -0.00066405     0.00057322     0.00020671    -0.00023318     0.00016827
   MO  35     0.00027950    -0.00008767    -0.00063956    -0.00036726    -0.00092276    -0.00007997     0.00000063    -0.00008192
   MO  36    -0.00040683     0.00076935     0.00063261    -0.00080491     0.00038198     0.00051222    -0.00003767    -0.00007399
   MO  37    -0.00013108     0.00134213     0.00112365    -0.00054421     0.00022033     0.00019403     0.00002220     0.00015676
   MO  38     0.00042908     0.00093656    -0.00141621     0.00168394    -0.00121884    -0.00120773    -0.00015308     0.00080939
   MO  39    -0.00049804    -0.00040698     0.00024378     0.00142078     0.00145482     0.00127177     0.00002725    -0.00015444
   MO  40     0.00009495     0.00084036    -0.00003333    -0.00131589    -0.00077402     0.00069909     0.00005858     0.00071580
   MO  41    -0.00018974    -0.00036674    -0.00004199     0.00011175     0.00114192     0.00179982     0.00011733    -0.00006215
   MO  42    -0.00008318     0.00048886    -0.00008934    -0.00072133     0.00017118    -0.00005980     0.00002995    -0.00034889
   MO  43    -0.00002621     0.00010350     0.00195665     0.00059314    -0.00136669     0.00028934    -0.00001095    -0.00035942
   MO  44    -0.00003356     0.00014629     0.00032443     0.00133079    -0.00028744     0.00310125     0.00009290     0.00101721
   MO  45    -0.00073193     0.00137906    -0.00017020    -0.00111877     0.00087884    -0.00037746     0.00008427    -0.00063811
   MO  46     0.00022500     0.00069233     0.00014414    -0.00013653     0.00001465    -0.00016573    -0.00039615    -0.00016852
   MO  47     0.00024345     0.00214644    -0.00024396    -0.00020815    -0.00042139     0.00029125     0.00045321     0.00016393
   MO  48    -0.00020493     0.00050322     0.00003741    -0.00014230     0.00032540    -0.00005134    -0.00018047     0.00025946
   MO  49     0.00175685     0.00008174    -0.00030998     0.00022832    -0.00059267     0.00001813     0.00003057    -0.00026301
   MO  50     0.00008174     0.00472797    -0.00021107    -0.00017772    -0.00130180     0.00017508    -0.00027672    -0.00015550
   MO  51    -0.00030998    -0.00021107     0.00554558    -0.00010289     0.00080826    -0.00030125     0.00003972     0.00010591
   MO  52     0.00022832    -0.00017772    -0.00010289     0.00356026    -0.00116259     0.00122203     0.00008906    -0.00023391
   MO  53    -0.00059267    -0.00130180     0.00080826    -0.00116259     0.00566475    -0.00063997    -0.00004028     0.00096748
   MO  54     0.00001813     0.00017508    -0.00030125     0.00122203    -0.00063997     0.00409049    -0.00000304     0.00002336
   MO  55     0.00003057    -0.00027672     0.00003972     0.00008906    -0.00004028    -0.00000304     0.00173661    -0.00002960
   MO  56    -0.00026301    -0.00015550     0.00010591    -0.00023391     0.00096748     0.00002336    -0.00002960     0.00256886
   MO  57     0.00010434     0.00014695    -0.00056601     0.00010761    -0.00042197     0.00024988    -0.00008761    -0.00027122
   MO  58    -0.00014825    -0.00004283     0.00002837    -0.00005480     0.00015685    -0.00001813    -0.00002421     0.00001425
   MO  59     0.00001195     0.00005866    -0.00002060     0.00009938    -0.00003981     0.00016302     0.00002069     0.00003101
   MO  60     0.00010303    -0.00101156    -0.00018295    -0.00017704    -0.00017773    -0.00002992     0.00004959    -0.00000315
   MO  61     0.00011904     0.00035721    -0.00011164     0.00017138    -0.00003790    -0.00025840    -0.00008332    -0.00000521
   MO  62     0.00021330     0.00008647    -0.00087306     0.00004769    -0.00014905     0.00022315     0.00001259    -0.00024473
   MO  63    -0.00065799    -0.00048733    -0.00010293    -0.00009193     0.00001938    -0.00025448     0.00003682    -0.00013327
   MO  64     0.00036526     0.00057977    -0.00054245    -0.00075803    -0.00050009    -0.00151065    -0.00005825    -0.00009675
   MO  65    -0.00023565     0.00009261    -0.00002061     0.00002912     0.00001264    -0.00001088     0.00040333     0.00002059
   MO  66    -0.00067174     0.00014839    -0.00006128     0.00006381     0.00005313     0.00000401    -0.00038385    -0.00000918
   MO  67     0.00016090    -0.00053481    -0.00003070     0.00014514     0.00005106    -0.00007446     0.00029508     0.00000125
   MO  68     0.00024852    -0.00039379    -0.00026483     0.00034501    -0.00018181    -0.00037335    -0.00000802     0.00003388
   MO  69    -0.00018334    -0.00040681     0.00037373     0.00051026     0.00035157     0.00099650     0.00004177    -0.00010002
   MO  70     0.00040145     0.00033833    -0.00033136     0.00036234    -0.00069370     0.00034058     0.00000470    -0.00063812
   MO  71     0.00000257     0.00036542    -0.00010879     0.00041649     0.00016432    -0.00014203    -0.00000796     0.00008709
   MO  72    -0.00004386    -0.00007315     0.00000320    -0.00007107    -0.00004343     0.00001604     0.00001143     0.00016818
   MO  73    -0.00000203    -0.00034223     0.00004503    -0.00037261     0.00044142     0.00003067     0.00000423     0.00008727
   MO  74    -0.00015076     0.00013924    -0.00032074    -0.00000275    -0.00056081    -0.00003902    -0.00001210     0.00066777
   MO  75    -0.00060184    -0.00033722     0.00023436    -0.00066838     0.00061930    -0.00055618     0.00002238    -0.00052782
   MO  76     0.00078676    -0.00007586    -0.00000618    -0.00014108     0.00004807     0.00010096    -0.00061366     0.00001689
   MO  77    -0.00053193     0.00061471    -0.00009261    -0.00006310     0.00007900     0.00007743    -0.00098248    -0.00001351
   MO  78     0.00079115     0.00019374    -0.00022272     0.00005724     0.00000577    -0.00018911     0.00018651    -0.00024166
   MO  79     0.00003916    -0.00009995     0.00100175     0.00009515    -0.00014512    -0.00015158     0.00004247    -0.00006249
   MO  80     0.00010660    -0.00019023     0.00035268     0.00004813     0.00095699     0.00046309     0.00003433    -0.00002896
   MO  81    -0.00004443     0.00006459     0.00002169     0.00071492     0.00013853    -0.00023862    -0.00003796     0.00000190
   MO  82    -0.00001680     0.00006445    -0.00017923     0.00062908    -0.00082119     0.00094270     0.00002147     0.00018600

                MO  57         MO  58         MO  59         MO  60         MO  61         MO  62         MO  63         MO  64
   MO   3     0.00029679     0.00007677    -0.00114449     0.00203072    -0.01476948     0.00036211    -0.00254100     0.00191863
   MO   4    -0.00550519    -0.00111659    -0.00003002    -0.00027744     0.00031008    -0.00132700     0.00282959     0.00344912
   MO   5    -0.00240686     0.00036612    -0.00046327    -0.00102987     0.00046020     0.00365306     0.00111247    -0.00007652
   MO   6    -0.01113078     0.00288658    -0.00064941    -0.00023711     0.00127659    -0.00118911    -0.00473547     0.00708719
   MO   7    -0.00368218     0.00059943     0.00113618     0.00040962    -0.00150803     0.00046191    -0.00379477    -0.01784423
   MO   8     0.00055158    -0.00027709    -0.00209344     0.02008209    -0.00119708     0.00288377     0.00145175    -0.00019457
   MO   9    -0.00019333     0.00003479    -0.00006631     0.00056330    -0.00008223     0.00008843     0.00013746    -0.00001360
   MO  10    -0.00009488     0.00038653     0.00208193     0.00003063    -0.00004444    -0.00005728     0.00001594     0.00019984
   MO  11     0.00031097    -0.00015363    -0.00010509     0.00057883    -0.00105770    -0.00002230     0.00029003    -0.00033396
   MO  12    -0.00003597     0.00004714    -0.00009230     0.00006378    -0.00080917     0.00018624    -0.00047385     0.00062638
   MO  13     0.00031153    -0.00011612     0.00000875    -0.00025784     0.00003716    -0.00007135    -0.00012555     0.00064480
   MO  14     0.00017856     0.00007034     0.00012920    -0.00014157    -0.00020100     0.00002466    -0.00013605    -0.00084068
   MO  15     0.00020117    -0.00013456    -0.00045402     0.00144658    -0.00023306     0.00020585     0.00018863    -0.00014771
   MO  16     0.00014434     0.00000640     0.00035837     0.00092950    -0.00011645     0.00021015     0.00009961     0.00035297
   MO  17     0.00044864    -0.00000951    -0.00008738     0.00148992    -0.00031600    -0.00004678     0.00022514    -0.00025652
   MO  18    -0.00060827     0.00017930    -0.00007919    -0.00003589     0.00085587     0.00026231    -0.00032692     0.00035447
   MO  19     0.00008149    -0.00018271     0.00035531    -0.00048555     0.00114106    -0.00017708     0.00016275    -0.00013057
   MO  20     0.00026294    -0.00004233    -0.00001882     0.00058939    -0.00038503     0.00016609    -0.00062007     0.00085689
   MO  21    -0.00114320    -0.00000408     0.00007350    -0.00031902     0.00005019     0.00003900     0.00042626    -0.00060554
   MO  22    -0.00088572    -0.00017261     0.00032154    -0.00025361    -0.00026851    -0.00009105     0.00003709     0.00072610
   MO  23     0.00107749    -0.00000450     0.00016545    -0.00067586     0.00021257    -0.00013010     0.00009877     0.00020968
   MO  24    -0.00001893     0.00029656    -0.00001858     0.00016906    -0.00027132    -0.00004203    -0.00003580    -0.00216057
   MO  25    -0.00071788     0.00111245    -0.00005547    -0.00050299    -0.00001752     0.00005508     0.00020428     0.00000575
   MO  26    -0.00064984    -0.00067604     0.00007713    -0.00053082    -0.00000366     0.00051774    -0.00030915    -0.00009944
   MO  27     0.00030995     0.00105046     0.00089238     0.00027753    -0.00002522     0.00018239    -0.00047001     0.00018929
   MO  28    -0.00009482     0.00027404    -0.00006679    -0.00099257    -0.00005967    -0.00006086    -0.00008587    -0.00003113
   MO  29     0.00041686    -0.00029887    -0.00118676    -0.00064468    -0.00004418    -0.00038889     0.00020758    -0.00016254
   MO  30     0.00007105     0.00002508     0.00084079    -0.00080257    -0.00025710    -0.00020610    -0.00017331     0.00005792
   MO  31     0.00000251     0.00018940    -0.00123025     0.00126363     0.00042171     0.00034552    -0.00020827     0.00006146
   MO  32    -0.00038991     0.00002170    -0.00144872    -0.00105143    -0.00018349     0.00028552    -0.00060906     0.00019801
   MO  33    -0.00106303    -0.00007038     0.00019648    -0.00048196    -0.00024361     0.00067567     0.00093865    -0.00026044
   MO  34    -0.00003207     0.00000778     0.00002471    -0.00075325    -0.00084171    -0.00022118    -0.00029558     0.00008451
   MO  35     0.00025677     0.00000952    -0.00069169     0.00047169     0.00051231     0.00013895    -0.00025473     0.00030760
   MO  36    -0.00099234    -0.00004115     0.00008804     0.00032065     0.00092923     0.00000343     0.00038703    -0.00027555
   MO  37     0.00035922    -0.00016934     0.00030731     0.00012016     0.00094738     0.00002782     0.00015938     0.00004792
   MO  38     0.00007462     0.00013943     0.00014263    -0.00062966     0.00057942     0.00000181    -0.00056590     0.00023108
   MO  39     0.00026851    -0.00006747    -0.00021719    -0.00052396     0.00029710    -0.00007208     0.00035341    -0.00099382
   MO  40    -0.00003051    -0.00043453     0.00000100     0.00060722    -0.00002744     0.00013254    -0.00024023    -0.00011575
   MO  41     0.00058689     0.00007238    -0.00017601     0.00022963     0.00003833    -0.00001725    -0.00021034    -0.00143692
   MO  42     0.00010269     0.00114138    -0.00017055     0.00009306     0.00014464     0.00011739    -0.00015091    -0.00013130
   MO  43    -0.00123698     0.00013784    -0.00014985    -0.00005155    -0.00007050     0.00033868     0.00017033    -0.00054036
   MO  44    -0.00009004     0.00049979    -0.00006168     0.00007180     0.00000268     0.00007378    -0.00021412    -0.00058211
   MO  45     0.00016003    -0.00029305     0.00007233    -0.00052115     0.00006654    -0.00000171     0.00057509    -0.00005675
   MO  46    -0.00005907    -0.00016887     0.00025816    -0.00020013     0.00013666     0.00006280    -0.00020302     0.00022892
   MO  47     0.00019608     0.00013743     0.00018022    -0.00054728     0.00014143     0.00001625    -0.00029361     0.00025989
   MO  48     0.00003621     0.00012687    -0.00011486    -0.00042356     0.00000057    -0.00018426     0.00059393    -0.00022242
   MO  49     0.00010434    -0.00014825     0.00001195     0.00010303     0.00011904     0.00021330    -0.00065799     0.00036526
   MO  50     0.00014695    -0.00004283     0.00005866    -0.00101156     0.00035721     0.00008647    -0.00048733     0.00057977
   MO  51    -0.00056601     0.00002837    -0.00002060    -0.00018295    -0.00011164    -0.00087306    -0.00010293    -0.00054245
   MO  52     0.00010761    -0.00005480     0.00009938    -0.00017704     0.00017138     0.00004769    -0.00009193    -0.00075803
   MO  53    -0.00042197     0.00015685    -0.00003981    -0.00017773    -0.00003790    -0.00014905     0.00001938    -0.00050009
   MO  54     0.00024988    -0.00001813     0.00016302    -0.00002992    -0.00025840     0.00022315    -0.00025448    -0.00151065
   MO  55    -0.00008761    -0.00002421     0.00002069     0.00004959    -0.00008332     0.00001259     0.00003682    -0.00005825
   MO  56    -0.00027122     0.00001425     0.00003101    -0.00000315    -0.00000521    -0.00024473    -0.00013327    -0.00009675
   MO  57     0.00240250    -0.00005089     0.00001033     0.00010320     0.00002471    -0.00056429    -0.00026002     0.00002023
   MO  58    -0.00005089     0.00228784     0.00012589     0.00001735    -0.00000285     0.00003255    -0.00011935    -0.00008995
   MO  59     0.00001033     0.00012589     0.00272684     0.00005281    -0.00011669    -0.00010170     0.00018333     0.00010307
   MO  60     0.00010320     0.00001735     0.00005281     0.00239896     0.00037460     0.00008524     0.00010065    -0.00007413
   MO  61     0.00002471    -0.00000285    -0.00011669     0.00037460     0.00174461     0.00004970     0.00004482     0.00000683
   MO  62    -0.00056429     0.00003255    -0.00010170     0.00008524     0.00004970     0.00148335     0.00019706    -0.00000805
   MO  63    -0.00026002    -0.00011935     0.00018333     0.00010065     0.00004482     0.00019706     0.00132356    -0.00010456
   MO  64     0.00002023    -0.00008995     0.00010307    -0.00007413     0.00000683    -0.00000805    -0.00010456     0.00183179
   MO  65     0.00001459    -0.00001806    -0.00005182    -0.00004372     0.00003099    -0.00001205    -0.00004183    -0.00001431
   MO  66    -0.00000634     0.00002118     0.00011234     0.00002996    -0.00000459    -0.00000572     0.00004323    -0.00000395
   MO  67    -0.00002371     0.00002478     0.00017750     0.00002711    -0.00001651    -0.00000773     0.00006396    -0.00001405
   MO  68     0.00008950     0.00003376    -0.00010774     0.00005275    -0.00004186    -0.00006410    -0.00003642     0.00011443
   MO  69     0.00009509     0.00004931    -0.00006963    -0.00001923    -0.00005237    -0.00006703     0.00028247    -0.00031393
   MO  70     0.00032020    -0.00015555     0.00002719     0.00008364     0.00006580     0.00005121    -0.00001207     0.00010376
   MO  71     0.00029067    -0.00003639     0.00005744    -0.00040196    -0.00010226    -0.00008966    -0.00009456     0.00004490
   MO  72    -0.00002668    -0.00018293     0.00001839     0.00006584    -0.00000377    -0.00002094     0.00001319    -0.00000500
   MO  73     0.00047826    -0.00001353     0.00001697     0.00025730     0.00008350    -0.00009679    -0.00011278    -0.00003554
   MO  74     0.00002754     0.00020772    -0.00006052     0.00007512     0.00000157    -0.00004285    -0.00005505     0.00002362
   MO  75     0.00024672    -0.00003891    -0.00005110    -0.00002345    -0.00014888    -0.00014281     0.00051229    -0.00019763
   MO  76     0.00007472    -0.00001352     0.00000421    -0.00002362     0.00003988     0.00002778    -0.00001828    -0.00000374
   MO  77     0.00005372    -0.00004100     0.00001083    -0.00005416     0.00006495    -0.00000186    -0.00012013     0.00003616
   MO  78    -0.00016201    -0.00016479     0.00010758     0.00002049     0.00011673     0.00021188    -0.00036326     0.00018326
   MO  79    -0.00007212    -0.00004028     0.00002727     0.00001024    -0.00001056    -0.00047304    -0.00016740    -0.00002957
   MO  80    -0.00020265    -0.00006836     0.00005503    -0.00001768     0.00001670     0.00010421    -0.00014660    -0.00023834
   MO  81    -0.00001844    -0.00000534    -0.00000724     0.00026643     0.00027641    -0.00001258     0.00010033    -0.00009315
   MO  82     0.00013048     0.00013091    -0.00003249     0.00000254    -0.00006113     0.00002904    -0.00006946    -0.00030032

                MO  65         MO  66         MO  67         MO  68         MO  69         MO  70         MO  71         MO  72
   MO   3    -0.00066438     0.00015286    -0.00036911     0.00108475     0.00026332    -0.00017248     0.00213811    -0.00021647
   MO   4    -0.00014923     0.00028668     0.00036012    -0.00286239     0.00300120     0.00138240    -0.00038915    -0.00070179
   MO   5     0.00028990     0.00001282     0.00007842    -0.00113199     0.00001223    -0.00010005    -0.00031367    -0.00022199
   MO   6     0.00019179    -0.00029086     0.00002141     0.00066852    -0.00571194    -0.00348182     0.00035169     0.00094624
   MO   7     0.00031275     0.00003567     0.00018016    -0.00199743     0.00601470     0.00024813    -0.00023894     0.00013937
   MO   8    -0.00101424     0.00092337     0.00039753     0.00176012    -0.00005475     0.00065019    -0.00477995     0.00060444
   MO   9     0.00044097    -0.00094462     0.00055940     0.00000727    -0.00000439     0.00003298    -0.00002664     0.00005524
   MO  10    -0.00037934     0.00091554     0.00150114    -0.00061204    -0.00009566    -0.00011442     0.00003875     0.00004680
   MO  11    -0.00009584     0.00008513     0.00002943    -0.00010615     0.00071112     0.00032577    -0.00028119    -0.00002466
   MO  12    -0.00001473    -0.00004128    -0.00008955     0.00027051    -0.00058257    -0.00002173    -0.00008829     0.00006517
   MO  13     0.00003013    -0.00008088    -0.00013098     0.00005941    -0.00035812     0.00006126     0.00008528     0.00000293
   MO  14    -0.00001350     0.00006060     0.00008859    -0.00012591     0.00054657     0.00008769     0.00002932    -0.00000574
   MO  15    -0.00010347    -0.00000765    -0.00013005     0.00071517     0.00009713     0.00024218    -0.00038023     0.00005257
   MO  16    -0.00021564     0.00035864     0.00046182     0.00035492    -0.00057326     0.00018461    -0.00021656     0.00004189
   MO  17    -0.00005642     0.00007509    -0.00005105    -0.00036089     0.00099990     0.00026590    -0.00061705     0.00004622
   MO  18     0.00007979    -0.00007411    -0.00020148    -0.00003174    -0.00056512    -0.00012872     0.00000681     0.00002514
   MO  19     0.00011482    -0.00005924     0.00012889    -0.00023660     0.00021139    -0.00006552     0.00029075    -0.00000946
   MO  20     0.00010192    -0.00013012    -0.00000555     0.00018205    -0.00072385     0.00008761    -0.00025071     0.00011399
   MO  21     0.00008542     0.00000757     0.00011699    -0.00049625     0.00091614    -0.00029308    -0.00003513    -0.00000094
   MO  22    -0.00014613     0.00011964     0.00003739     0.00013447    -0.00016124    -0.00047470     0.00005922     0.00010360
   MO  23    -0.00002054     0.00006221     0.00013443     0.00023681     0.00067111     0.00040973     0.00045254    -0.00012360
   MO  24    -0.00000944     0.00007182     0.00010001    -0.00031876     0.00247853    -0.00009692    -0.00008259     0.00001836
   MO  25     0.00012646    -0.00006964     0.00009357    -0.00000604     0.00022408    -0.00010719     0.00005968    -0.00014055
   MO  26    -0.00001346    -0.00002215    -0.00011009    -0.00000841     0.00016564     0.00063861     0.00006678     0.00000360
   MO  27    -0.00003326     0.00003069     0.00007597     0.00007087    -0.00014036     0.00028446    -0.00002784    -0.00000312
   MO  28    -0.00056052     0.00026225    -0.00029080    -0.00011322     0.00009463    -0.00008776     0.00001678    -0.00005128
   MO  29     0.00010861     0.00003119    -0.00016131     0.00022187     0.00007814    -0.00002695     0.00015790     0.00002908
   MO  30     0.00037321     0.00050437    -0.00016270    -0.00016679    -0.00003093    -0.00005473     0.00012254     0.00000572
   MO  31    -0.00007098     0.00041998    -0.00041641     0.00019192    -0.00002270     0.00011352    -0.00019560    -0.00000898
   MO  32     0.00032146    -0.00040968    -0.00011513    -0.00012957     0.00007642    -0.00002764     0.00008037    -0.00002749
   MO  33    -0.00003563     0.00016681     0.00006011    -0.00017055     0.00025425    -0.00013502     0.00004803    -0.00002757
   MO  34    -0.00000805    -0.00007399    -0.00021404    -0.00058135    -0.00002694    -0.00014160     0.00024116    -0.00003716
   MO  35     0.00007338    -0.00008469    -0.00010195     0.00015665    -0.00019095     0.00009230    -0.00016373     0.00002653
   MO  36     0.00008405    -0.00003772    -0.00010654    -0.00067250     0.00028886    -0.00005483    -0.00021496     0.00000393
   MO  37     0.00013093    -0.00011685    -0.00010249    -0.00065900     0.00025019     0.00006014     0.00006825    -0.00000989
   MO  38     0.00005501    -0.00004269     0.00001430     0.00046172    -0.00128296    -0.00034165     0.00053988     0.00000822
   MO  39     0.00007798    -0.00002058     0.00003364    -0.00031170     0.00132081     0.00009313     0.00028747    -0.00006625
   MO  40     0.00003439    -0.00007147    -0.00012351    -0.00020245    -0.00047276    -0.00011092    -0.00029918     0.00019518
   MO  41     0.00008575    -0.00005496    -0.00003613    -0.00023832     0.00021583     0.00007371    -0.00003167    -0.00002059
   MO  42    -0.00000347    -0.00003361    -0.00007634    -0.00012964    -0.00028269    -0.00007160    -0.00010257    -0.00021235
   MO  43    -0.00007497    -0.00001716    -0.00003035    -0.00005906    -0.00012497    -0.00010227    -0.00026565    -0.00000083
   MO  44     0.00001789    -0.00000806    -0.00003219    -0.00026744     0.00211768    -0.00021912    -0.00008209     0.00002222
   MO  45     0.00020614    -0.00017051    -0.00003131    -0.00044954    -0.00018102    -0.00068501     0.00010161     0.00000347
   MO  46    -0.00055489     0.00024454    -0.00001380    -0.00016108    -0.00014823    -0.00001467     0.00004515    -0.00000146
   MO  47     0.00047490    -0.00050744     0.00044259    -0.00033834    -0.00010498     0.00015707     0.00024932    -0.00005048
   MO  48    -0.00046589    -0.00044685     0.00009314    -0.00024461     0.00007744    -0.00029148     0.00012445    -0.00001024
   MO  49    -0.00023565    -0.00067174     0.00016090     0.00024852    -0.00018334     0.00040145     0.00000257    -0.00004386
   MO  50     0.00009261     0.00014839    -0.00053481    -0.00039379    -0.00040681     0.00033833     0.00036542    -0.00007315
   MO  51    -0.00002061    -0.00006128    -0.00003070    -0.00026483     0.00037373    -0.00033136    -0.00010879     0.00000320
   MO  52     0.00002912     0.00006381     0.00014514     0.00034501     0.00051026     0.00036234     0.00041649    -0.00007107
   MO  53     0.00001264     0.00005313     0.00005106    -0.00018181     0.00035157    -0.00069370     0.00016432    -0.00004343
   MO  54    -0.00001088     0.00000401    -0.00007446    -0.00037335     0.00099650     0.00034058    -0.00014203     0.00001604
   MO  55     0.00040333    -0.00038385     0.00029508    -0.00000802     0.00004177     0.00000470    -0.00000796     0.00001143
   MO  56     0.00002059    -0.00000918     0.00000125     0.00003388    -0.00010002    -0.00063812     0.00008709     0.00016818
   MO  57     0.00001459    -0.00000634    -0.00002371     0.00008950     0.00009509     0.00032020     0.00029067    -0.00002668
   MO  58    -0.00001806     0.00002118     0.00002478     0.00003376     0.00004931    -0.00015555    -0.00003639    -0.00018293
   MO  59    -0.00005182     0.00011234     0.00017750    -0.00010774    -0.00006963     0.00002719     0.00005744     0.00001839
   MO  60    -0.00004372     0.00002996     0.00002711     0.00005275    -0.00001923     0.00008364    -0.00040196     0.00006584
   MO  61     0.00003099    -0.00000459    -0.00001651    -0.00004186    -0.00005237     0.00006580    -0.00010226    -0.00000377
   MO  62    -0.00001205    -0.00000572    -0.00000773    -0.00006410    -0.00006703     0.00005121    -0.00008966    -0.00002094
   MO  63    -0.00004183     0.00004323     0.00006396    -0.00003642     0.00028247    -0.00001207    -0.00009456     0.00001319
   MO  64    -0.00001431    -0.00000395    -0.00001405     0.00011443    -0.00031393     0.00010376     0.00004490    -0.00000500
   MO  65     0.00079016    -0.00001032    -0.00003342    -0.00006006    -0.00001590    -0.00000447     0.00003663     0.00000886
   MO  66    -0.00001032     0.00107708    -0.00021208     0.00003629     0.00003293     0.00000435     0.00000236     0.00000757
   MO  67    -0.00003342    -0.00021208     0.00080254     0.00005394     0.00003560    -0.00000024     0.00003005    -0.00000392
   MO  68    -0.00006006     0.00003629     0.00005394     0.00093524    -0.00019859    -0.00004072     0.00009158     0.00000163
   MO  69    -0.00001590     0.00003293     0.00003560    -0.00019859     0.00199390     0.00011688    -0.00001933    -0.00002390
   MO  70    -0.00000447     0.00000435    -0.00000024    -0.00004072     0.00011688     0.00114401     0.00004140    -0.00006834
   MO  71     0.00003663     0.00000236     0.00003005     0.00009158    -0.00001933     0.00004140     0.00063513    -0.00003944
   MO  72     0.00000886     0.00000757    -0.00000392     0.00000163    -0.00002390    -0.00006834    -0.00003944     0.00026231
   MO  73     0.00003049    -0.00000843    -0.00001116    -0.00002876     0.00003026    -0.00006580     0.00001687     0.00002477
   MO  74     0.00000841    -0.00001126    -0.00001303     0.00005616    -0.00010185    -0.00005666     0.00000438     0.00005489
   MO  75     0.00001120     0.00000865     0.00000651    -0.00008170    -0.00047676    -0.00014775    -0.00003104    -0.00002140
   MO  76    -0.00075485    -0.00028138     0.00002200    -0.00001635    -0.00002505     0.00000420    -0.00003425    -0.00001760
   MO  77    -0.00022013     0.00080109    -0.00052354    -0.00001460    -0.00007460     0.00003457     0.00000400     0.00000053
   MO  78     0.00007216    -0.00015691     0.00008777     0.00003625    -0.00032064     0.00018781     0.00002285    -0.00003535
   MO  79    -0.00000357    -0.00002593     0.00001106     0.00002168    -0.00003231     0.00002099     0.00000706     0.00000280
   MO  80     0.00002550    -0.00001346     0.00000242    -0.00008858    -0.00031349     0.00016820     0.00000559    -0.00002971
   MO  81    -0.00002348     0.00003476     0.00004023     0.00010578    -0.00004656     0.00000545     0.00004191    -0.00000911
   MO  82     0.00001684    -0.00001964    -0.00001966     0.00006279    -0.00035391    -0.00006223    -0.00000557     0.00001300

                MO  73         MO  74         MO  75         MO  76         MO  77         MO  78         MO  79         MO  80
   MO   3    -0.00158573    -0.00059117     0.00003261     0.00005302     0.00003069    -0.00001090     0.00021153     0.00014471
   MO   4     0.00204165    -0.00613993     0.00682067     0.00028705     0.00013886    -0.00023442     0.00057011     0.00327767
   MO   5     0.00033429     0.00022670    -0.00006552     0.00000474     0.00004042    -0.00011150     0.00179102    -0.00020062
   MO   6    -0.00058308     0.00693192    -0.00309175    -0.00005970    -0.00068024    -0.00350517     0.00098956     0.00722462
   MO   7     0.00020438     0.00087138     0.00062259    -0.00002714    -0.00054775    -0.00313048     0.00059341     0.00470132
   MO   8     0.00304395     0.00048569     0.00006515    -0.00007980     0.00012581     0.00048915    -0.00040402    -0.00049368
   MO   9     0.00002377    -0.00001316     0.00013768    -0.00257049    -0.00396100     0.00078421     0.00018984     0.00004885
   MO  10     0.00002229     0.00005522    -0.00003153    -0.00003012    -0.00003594    -0.00001616     0.00002194     0.00001421
   MO  11     0.00022949    -0.00046386     0.00045851     0.00000096    -0.00001301    -0.00012431     0.00001356     0.00016205
   MO  12    -0.00011468     0.00037808    -0.00038976    -0.00001569     0.00002122     0.00016043    -0.00006949    -0.00027513
   MO  13    -0.00012793     0.00008717    -0.00007446     0.00001662     0.00005016     0.00009922     0.00004018    -0.00022897
   MO  14    -0.00006361     0.00000001    -0.00028660     0.00000661     0.00001414    -0.00011390     0.00003979     0.00006696
   MO  15     0.00014665     0.00002796    -0.00010456    -0.00009543    -0.00011627     0.00006134    -0.00002627    -0.00010267
   MO  16     0.00002760     0.00008727     0.00006159     0.00000762     0.00007214     0.00014585    -0.00005338    -0.00022412
   MO  17     0.00039519    -0.00026621     0.00013857     0.00005814     0.00009294    -0.00022787     0.00004587     0.00007284
   MO  18    -0.00017660     0.00039306    -0.00048445     0.00015945     0.00030512     0.00013510    -0.00006626    -0.00012161
   MO  19    -0.00010063    -0.00001325    -0.00013240    -0.00016089    -0.00025650     0.00001084     0.00008538     0.00002502
   MO  20     0.00000276     0.00046579    -0.00047598    -0.00030851    -0.00036410     0.00032377    -0.00005245    -0.00034241
   MO  21     0.00010958    -0.00027985     0.00017032    -0.00017798    -0.00025101    -0.00020435     0.00001108     0.00048349
   MO  22    -0.00013256     0.00025725    -0.00030395     0.00015191     0.00018481    -0.00008303     0.00009565    -0.00002867
   MO  23    -0.00010617    -0.00013741    -0.00027178     0.00004919     0.00006197     0.00002220    -0.00018626    -0.00023665
   MO  24     0.00010162    -0.00003653    -0.00097882    -0.00002275    -0.00008788    -0.00053114     0.00007194     0.00041237
   MO  25    -0.00024203    -0.00013174     0.00008137    -0.00031562    -0.00035253    -0.00003239     0.00023519    -0.00002871
   MO  26    -0.00049688     0.00006639    -0.00080353    -0.00000319     0.00008173     0.00029259     0.00011134    -0.00010973
   MO  27    -0.00011964     0.00061436    -0.00075905     0.00004072     0.00001727     0.00000751    -0.00019479    -0.00025686
   MO  28    -0.00001672    -0.00011211     0.00006959     0.00080531     0.00081223    -0.00018909     0.00002731     0.00007471
   MO  29    -0.00016038     0.00029710    -0.00000693    -0.00016119     0.00003370    -0.00038669     0.00008697    -0.00030690
   MO  30     0.00002474    -0.00007977     0.00003605    -0.00060710     0.00047024    -0.00000855     0.00002879     0.00010670
   MO  31     0.00007473     0.00009422    -0.00020135    -0.00001433     0.00068754    -0.00002317    -0.00011411    -0.00005017
   MO  32    -0.00004810     0.00006572    -0.00030928    -0.00018814    -0.00032976     0.00022227    -0.00014010     0.00015111
   MO  33    -0.00014852    -0.00015448     0.00034195    -0.00003002     0.00004211    -0.00018805    -0.00028335     0.00006372
   MO  34    -0.00009497    -0.00005665     0.00015575     0.00012337     0.00013585    -0.00009918     0.00008994     0.00008916
   MO  35     0.00007662     0.00014342    -0.00020496    -0.00004436    -0.00000012     0.00011913    -0.00006030    -0.00014048
   MO  36    -0.00007566    -0.00012220     0.00009605    -0.00001591     0.00002313    -0.00008068     0.00010825     0.00013641
   MO  37     0.00011737    -0.00008325     0.00004555     0.00001550    -0.00002874    -0.00005079     0.00001809     0.00007956
   MO  38    -0.00032597     0.00053753    -0.00047449    -0.00001082     0.00006958     0.00028347    -0.00004115    -0.00024888
   MO  39    -0.00000009    -0.00027216     0.00002679    -0.00003126    -0.00006781    -0.00029256    -0.00005530     0.00016342
   MO  40     0.00015449     0.00036337    -0.00012193     0.00003858     0.00001108    -0.00000077    -0.00002973     0.00000777
   MO  41     0.00031892    -0.00021970     0.00026883     0.00000287    -0.00001151    -0.00005172    -0.00005351     0.00042081
   MO  42     0.00009246     0.00002304     0.00027339     0.00005583     0.00001091    -0.00001370    -0.00005131     0.00003349
   MO  43    -0.00045055     0.00000592    -0.00009845     0.00004076    -0.00005672     0.00000123     0.00020186    -0.00000545
   MO  44     0.00004079     0.00032864    -0.00140336    -0.00001351    -0.00005561    -0.00033979    -0.00004155    -0.00010637
   MO  45     0.00014372    -0.00041562     0.00085781    -0.00010590    -0.00016044    -0.00003707    -0.00015439    -0.00013567
   MO  46    -0.00003230    -0.00016951    -0.00012614     0.00076937     0.00047222     0.00021500     0.00000937     0.00002241
   MO  47    -0.00010142     0.00015869    -0.00018812    -0.00035661    -0.00082145     0.00022676    -0.00006231    -0.00005748
   MO  48    -0.00012471     0.00021074     0.00051152     0.00105373    -0.00045255    -0.00062407    -0.00007876    -0.00017778
   MO  49    -0.00000203    -0.00015076    -0.00060184     0.00078676    -0.00053193     0.00079115     0.00003916     0.00010660
   MO  50    -0.00034223     0.00013924    -0.00033722    -0.00007586     0.00061471     0.00019374    -0.00009995    -0.00019023
   MO  51     0.00004503    -0.00032074     0.00023436    -0.00000618    -0.00009261    -0.00022272     0.00100175     0.00035268
   MO  52    -0.00037261    -0.00000275    -0.00066838    -0.00014108    -0.00006310     0.00005724     0.00009515     0.00004813
   MO  53     0.00044142    -0.00056081     0.00061930     0.00004807     0.00007900     0.00000577    -0.00014512     0.00095699
   MO  54     0.00003067    -0.00003902    -0.00055618     0.00010096     0.00007743    -0.00018911    -0.00015158     0.00046309
   MO  55     0.00000423    -0.00001210     0.00002238    -0.00061366    -0.00098248     0.00018651     0.00004247     0.00003433
   MO  56     0.00008727     0.00066777    -0.00052782     0.00001689    -0.00001351    -0.00024166    -0.00006249    -0.00002896
   MO  57     0.00047826     0.00002754     0.00024672     0.00007472     0.00005372    -0.00016201    -0.00007212    -0.00020265
   MO  58    -0.00001353     0.00020772    -0.00003891    -0.00001352    -0.00004100    -0.00016479    -0.00004028    -0.00006836
   MO  59     0.00001697    -0.00006052    -0.00005110     0.00000421     0.00001083     0.00010758     0.00002727     0.00005503
   MO  60     0.00025730     0.00007512    -0.00002345    -0.00002362    -0.00005416     0.00002049     0.00001024    -0.00001768
   MO  61     0.00008350     0.00000157    -0.00014888     0.00003988     0.00006495     0.00011673    -0.00001056     0.00001670
   MO  62    -0.00009679    -0.00004285    -0.00014281     0.00002778    -0.00000186     0.00021188    -0.00047304     0.00010421
   MO  63    -0.00011278    -0.00005505     0.00051229    -0.00001828    -0.00012013    -0.00036326    -0.00016740    -0.00014660
   MO  64    -0.00003554     0.00002362    -0.00019763    -0.00000374     0.00003616     0.00018326    -0.00002957    -0.00023834
   MO  65     0.00003049     0.00000841     0.00001120    -0.00075485    -0.00022013     0.00007216    -0.00000357     0.00002550
   MO  66    -0.00000843    -0.00001126     0.00000865    -0.00028138     0.00080109    -0.00015691    -0.00002593    -0.00001346
   MO  67    -0.00001116    -0.00001303     0.00000651     0.00002200    -0.00052354     0.00008777     0.00001106     0.00000242
   MO  68    -0.00002876     0.00005616    -0.00008170    -0.00001635    -0.00001460     0.00003625     0.00002168    -0.00008858
   MO  69     0.00003026    -0.00010185    -0.00047676    -0.00002505    -0.00007460    -0.00032064    -0.00003231    -0.00031349
   MO  70    -0.00006580    -0.00005666    -0.00014775     0.00000420     0.00003457     0.00018781     0.00002099     0.00016820
   MO  71     0.00001687     0.00000438    -0.00003104    -0.00003425     0.00000400     0.00002285     0.00000706     0.00000559
   MO  72     0.00002477     0.00005489    -0.00002140    -0.00001760     0.00000053    -0.00003535     0.00000280    -0.00002971
   MO  73     0.00067193    -0.00004399     0.00010750    -0.00000379     0.00001558     0.00002430     0.00004233     0.00004239
   MO  74    -0.00004399     0.00071487    -0.00019697    -0.00001484    -0.00006788    -0.00035516    -0.00004967    -0.00022246
   MO  75     0.00010750    -0.00019697     0.00129293    -0.00002231    -0.00008797    -0.00043795    -0.00001003     0.00009903
   MO  76    -0.00000379    -0.00001484    -0.00002231     0.00170517    -0.00000554    -0.00001860    -0.00002101     0.00000723
   MO  77     0.00001558    -0.00006788    -0.00008797    -0.00000554     0.00156262    -0.00008926    -0.00003928     0.00004710
   MO  78     0.00002430    -0.00035516    -0.00043795    -0.00001860    -0.00008926     0.00109238     0.00000478     0.00025384
   MO  79     0.00004233    -0.00004967    -0.00001003    -0.00002101    -0.00003928     0.00000478     0.00062645     0.00003580
   MO  80     0.00004239    -0.00022246     0.00009903     0.00000723     0.00004710     0.00025384     0.00003580     0.00081264
   MO  81    -0.00001156    -0.00000786    -0.00000432    -0.00000245     0.00000381     0.00001106     0.00002573     0.00002199
   MO  82    -0.00002774     0.00009946     0.00003052    -0.00000016    -0.00001387    -0.00009065    -0.00000266     0.00012480

                MO  81         MO  82
   MO   3    -0.00203597     0.00015437
   MO   4    -0.00018475    -0.00464156
   MO   5    -0.00013124     0.00008660
   MO   6    -0.00011553    -0.00512600
   MO   7     0.00018338     0.00006750
   MO   8     0.00271851     0.00023947
   MO   9     0.00006522    -0.00003445
   MO  10    -0.00000263    -0.00002315
   MO  11    -0.00025186    -0.00012824
   MO  12    -0.00026530     0.00024244
   MO  13    -0.00006034    -0.00000355
   MO  14     0.00002512     0.00032619
   MO  15     0.00008495     0.00014628
   MO  16     0.00006614     0.00009770
   MO  17    -0.00062456    -0.00002773
   MO  18    -0.00004926     0.00013385
   MO  19     0.00031450    -0.00004855
   MO  20    -0.00048242     0.00026192
   MO  21    -0.00003066    -0.00031913
   MO  22     0.00007806    -0.00022115
   MO  23     0.00047954     0.00011623
   MO  24     0.00008353     0.00061661
   MO  25     0.00016061    -0.00008729
   MO  26     0.00004430     0.00043818
   MO  27    -0.00015509     0.00048884
   MO  28    -0.00016212    -0.00007268
   MO  29     0.00011865     0.00025034
   MO  30    -0.00016523    -0.00005728
   MO  31     0.00012655     0.00003443
   MO  32    -0.00027106    -0.00003472
   MO  33     0.00012555    -0.00007452
   MO  34    -0.00007967    -0.00008292
   MO  35    -0.00011370     0.00008676
   MO  36    -0.00006590    -0.00012274
   MO  37     0.00000783    -0.00010348
   MO  38     0.00048847     0.00034125
   MO  39     0.00027573    -0.00006182
   MO  40    -0.00033414     0.00020965
   MO  41    -0.00004932     0.00025194
   MO  42    -0.00009078     0.00002588
   MO  43     0.00001818     0.00023532
   MO  44    -0.00008362     0.00035987
   MO  45     0.00000687    -0.00014709
   MO  46     0.00003343    -0.00009248
   MO  47    -0.00000387     0.00003988
   MO  48     0.00012503     0.00003518
   MO  49    -0.00004443    -0.00001680
   MO  50     0.00006459     0.00006445
   MO  51     0.00002169    -0.00017923
   MO  52     0.00071492     0.00062908
   MO  53     0.00013853    -0.00082119
   MO  54    -0.00023862     0.00094270
   MO  55    -0.00003796     0.00002147
   MO  56     0.00000190     0.00018600
   MO  57    -0.00001844     0.00013048
   MO  58    -0.00000534     0.00013091
   MO  59    -0.00000724    -0.00003249
   MO  60     0.00026643     0.00000254
   MO  61     0.00027641    -0.00006113
   MO  62    -0.00001258     0.00002904
   MO  63     0.00010033    -0.00006946
   MO  64    -0.00009315    -0.00030032
   MO  65    -0.00002348     0.00001684
   MO  66     0.00003476    -0.00001964
   MO  67     0.00004023    -0.00001966
   MO  68     0.00010578     0.00006279
   MO  69    -0.00004656    -0.00035391
   MO  70     0.00000545    -0.00006223
   MO  71     0.00004191    -0.00000557
   MO  72    -0.00000911     0.00001300
   MO  73    -0.00001156    -0.00002774
   MO  74    -0.00000786     0.00009946
   MO  75    -0.00000432     0.00003052
   MO  76    -0.00000245    -0.00000016
   MO  77     0.00000381    -0.00001387
   MO  78     0.00001106    -0.00009065
   MO  79     0.00002573    -0.00000266
   MO  80     0.00002199     0.00012480
   MO  81     0.00049514     0.00003786
   MO  82     0.00003786     0.00092574

Natural orbital populations,block 1
              MO     1       MO     2       MO     3       MO     4       MO     5       MO     6       MO     7       MO     8
  occ(*)=     2.00000000     2.00000000     1.97562335     1.97488189     1.97090425     1.96126663     1.95180684     1.93525571
              MO     9       MO    10       MO    11       MO    12       MO    13       MO    14       MO    15       MO    16
  occ(*)=     0.06041302     0.02418005     0.02391549     0.01694112     0.01623076     0.01370331     0.01244800     0.00666060
              MO    17       MO    18       MO    19       MO    20       MO    21       MO    22       MO    23       MO    24
  occ(*)=     0.00614544     0.00566696     0.00476421     0.00462206     0.00454276     0.00450041     0.00432507     0.00413365
              MO    25       MO    26       MO    27       MO    28       MO    29       MO    30       MO    31       MO    32
  occ(*)=     0.00243112     0.00192741     0.00164898     0.00119195     0.00102734     0.00077743     0.00073273     0.00070874
              MO    33       MO    34       MO    35       MO    36       MO    37       MO    38       MO    39       MO    40
  occ(*)=     0.00068246     0.00057481     0.00055423     0.00050790     0.00047608     0.00041950     0.00035825     0.00033099
              MO    41       MO    42       MO    43       MO    44       MO    45       MO    46       MO    47       MO    48
  occ(*)=     0.00027675     0.00024949     0.00023895     0.00022590     0.00022456     0.00021470     0.00017216     0.00014498
              MO    49       MO    50       MO    51       MO    52       MO    53       MO    54       MO    55       MO    56
  occ(*)=     0.00010558     0.00009871     0.00009228     0.00007839     0.00006849     0.00006802     0.00006410     0.00005120
              MO    57       MO    58       MO    59       MO    60       MO    61       MO    62       MO    63       MO    64
  occ(*)=     0.00004798     0.00004655     0.00003792     0.00003419     0.00002931     0.00002449     0.00002055     0.00001667
              MO    65       MO    66       MO    67       MO    68       MO    69       MO    70       MO    71       MO    72
  occ(*)=     0.00001331     0.00001131     0.00001062     0.00000974     0.00000924     0.00000592     0.00000530     0.00000526
              MO    73       MO    74       MO    75       MO    76       MO    77       MO    78       MO    79       MO    80
  occ(*)=     0.00000433     0.00000347     0.00000298     0.00000200     0.00000128     0.00000121     0.00000030     0.00000018
              MO    81       MO    82       MO
  occ(*)=     0.00000009     0.00000005


 total number of electrons =   16.0000000000

 test slabel:                     1
  a  


          Mulliken population analysis


  NOTE: For HERMIT use spherical harmonics basis sets !!!
 

                        a   partial gross atomic populations
   ao class       1a         2a         3a         4a         5a         6a  
     2_ s       1.993633   0.000003   0.000205   1.018164   0.031372   0.055369
     2_ p      -0.000000   0.000004   0.000431   0.135388   0.659191   0.245931
     2_ d       0.000001   0.000000  -0.000008  -0.007880   0.000124   0.000278
     4_ s      -0.000001   1.998050   0.601824   0.000168   0.001139   0.087026
     4_ p      -0.000006   0.000060   0.514694   0.001615   0.001780   1.123056
     4_ d      -0.000000   0.000001   0.000954   0.000032  -0.000098   0.017278
     5_ s       0.002536   0.000001   0.000124   0.610189   0.476791  -0.002728
     5_ p       0.000781  -0.000000   0.000028   0.028157  -0.039936   0.005815
     6_ s      -0.000002   0.000417   0.173462   0.000235   0.000539   0.293027
     6_ p      -0.000000   0.000272   0.021964   0.000029   0.000023   0.017004
     7_ s       0.000000   0.001029   0.623025   0.000221  -0.000785   0.094860
     7_ p       0.000000   0.000175   0.038383   0.000034  -0.000119   0.016519
     8_ s       0.002268  -0.000012   0.000637   0.202227   0.847969  -0.000900
     8_ p       0.000789   0.000000  -0.000100  -0.013697  -0.007086   0.008729
 
   ao class       7a         8a         9a        10a        11a        12a  
     2_ s       0.203144   0.008465  -0.000000   0.000005   0.000000   0.003265
     2_ p       1.308982   0.064539   0.000006   0.022951   0.000009   0.004023
     2_ d       0.000784  -0.000122  -0.000000   0.001113   0.000001   0.001270
     4_ s       0.012216   1.115423   0.000001   0.000010   0.000207   0.000002
     4_ p       0.275012   0.726803   0.059211   0.000022   0.011357   0.000029
     4_ d       0.003952   0.008365   0.000557   0.000001   0.001301   0.000008
     5_ s      -0.005253   0.000373   0.000000  -0.000001  -0.000000   0.002950
     5_ p       0.030164   0.000767   0.000000   0.000175  -0.000000  -0.000281
     6_ s       0.076192   0.023069   0.000001   0.000003   0.009803   0.000003
     6_ p       0.002690   0.000295   0.000348  -0.000001  -0.000351   0.000002
     7_ s       0.015362  -0.012793   0.000000   0.000001   0.001266   0.000037
     7_ p       0.004996  -0.003761   0.000288   0.000000   0.000321  -0.000002
     8_ s      -0.005077   0.001243   0.000000   0.000001   0.000000   0.005815
     8_ p       0.028644   0.002588  -0.000000  -0.000100  -0.000000  -0.000179
 
   ao class      13a        14a        15a        16a        17a        18a  
     2_ s       0.000008   0.001703   0.000145   0.000086   0.003358   0.000187
     2_ p       0.000038   0.009659   0.006419   0.003082   0.000324   0.000087
     2_ d       0.000010   0.000706   0.001746   0.003362   0.002097   0.000258
     4_ s       0.000326   0.000131   0.000002   0.000395   0.000117   0.004009
     4_ p       0.008434   0.000103   0.000005   0.000028   0.000022   0.000519
     4_ d       0.000682   0.000024   0.000001  -0.000000   0.000003   0.000493
     5_ s       0.000010   0.000377   0.002818  -0.000360  -0.000704   0.000012
     5_ p       0.000001   0.000214   0.000218   0.000689   0.000154   0.000041
     6_ s       0.000701  -0.000016   0.000001  -0.000015  -0.000004  -0.000208
     6_ p       0.000306  -0.000004   0.000000   0.000001  -0.000001   0.000037
     7_ s       0.006192   0.000032  -0.000001  -0.000001   0.000006   0.000113
     7_ p      -0.000502  -0.000004  -0.000000   0.000001  -0.000000   0.000112
     8_ s       0.000021   0.000466   0.000680  -0.001152   0.000363  -0.000006
     8_ p       0.000004   0.000310   0.000414   0.000543   0.000410   0.000015
 
   ao class      19a        20a        21a        22a        23a        24a  
     2_ s       0.000000   0.000001  -0.000000   0.000009   0.000002   0.000000
     2_ p       0.000046  -0.000001   0.000000   0.000013  -0.000004   0.000001
     2_ d       0.003964   0.003187   0.000000   0.000133   0.000007   0.000000
     4_ s       0.000001   0.000001   0.000000   0.000215   0.000162   0.000000
     4_ p       0.000005   0.000004   0.000004   0.000446   0.003019   0.000015
     4_ d       0.000048   0.000072   0.004076   0.002929   0.001820   0.003141
     5_ s      -0.000000   0.000000  -0.000000   0.000001  -0.000001   0.000000
     5_ p       0.000624   0.000323   0.000000   0.000033   0.000001  -0.000000
     6_ s       0.000002   0.000001   0.000000   0.000054  -0.000341  -0.000002
     6_ p       0.000004   0.000006   0.000138   0.000197   0.000338   0.000504
     7_ s      -0.000002   0.000002  -0.000001   0.000006  -0.000813  -0.000002
     7_ p       0.000008   0.000011   0.000325   0.000438   0.000135   0.000476
     8_ s       0.000000   0.000001   0.000000   0.000009   0.000001   0.000000
     8_ p       0.000063   0.001014  -0.000000   0.000018  -0.000002  -0.000000
 
   ao class      25a        26a        27a        28a        29a        30a  
     2_ s       0.000518   0.000013   0.000001   0.000000   0.000190   0.000001
     2_ p       0.000146   0.002488   0.000023  -0.000000   0.000059   0.000011
     2_ d       0.001347   0.000241   0.000005   0.000001   0.000455   0.000021
     4_ s       0.000001   0.000008   0.000085  -0.000000  -0.000001   0.000048
     4_ p      -0.000002  -0.000002   0.000279   0.000613   0.000002  -0.000012
     4_ d       0.000002   0.000018   0.000942   0.000467   0.000000   0.000299
     5_ s       0.000191  -0.000989  -0.000008  -0.000001  -0.000111   0.000000
     5_ p       0.000059   0.000453   0.000005   0.000000   0.000229   0.000000
     6_ s      -0.000001   0.000001   0.000113   0.000081  -0.000001   0.000006
     6_ p       0.000000   0.000000   0.000066  -0.000056   0.000000   0.000195
     7_ s       0.000001  -0.000003   0.000076   0.000065   0.000001  -0.000017
     7_ p      -0.000000   0.000001   0.000069   0.000024   0.000000   0.000224
     8_ s       0.000103  -0.000767  -0.000009  -0.000000  -0.000062  -0.000003
     8_ p       0.000068   0.000466   0.000003  -0.000001   0.000265   0.000006
 
   ao class      31a        32a        33a        34a        35a        36a  
     2_ s       0.000077   0.000020   0.000005   0.000000   0.000002   0.000001
     2_ p       0.000059   0.000060   0.000033   0.000053   0.000103   0.000003
     2_ d       0.000387   0.000315   0.000121   0.000057   0.000196   0.000002
     4_ s      -0.000000   0.000005   0.000001   0.000000  -0.000000  -0.000034
     4_ p       0.000093   0.000089   0.000323   0.000000   0.000005   0.000045
     4_ d       0.000007   0.000014   0.000025   0.000000   0.000002   0.000072
     5_ s       0.000038   0.000001  -0.000006   0.000000   0.000093  -0.000000
     5_ p      -0.000016   0.000051   0.000061   0.000220   0.000028   0.000000
     6_ s      -0.000000  -0.000000   0.000000  -0.000000   0.000003   0.000351
     6_ p       0.000019   0.000025   0.000061   0.000000   0.000001   0.000146
     7_ s       0.000000  -0.000000   0.000001   0.000000  -0.000000  -0.000046
     7_ p       0.000010   0.000014   0.000023   0.000000   0.000000  -0.000031
     8_ s       0.000003  -0.000003  -0.000001  -0.000000   0.000149  -0.000001
     8_ p       0.000057   0.000120   0.000035   0.000244  -0.000029  -0.000000
 
   ao class      37a        38a        39a        40a        41a        42a  
     2_ s       0.000001   0.000000   0.000019   0.000014  -0.000000   0.000012
     2_ p       0.000000   0.000001   0.000082   0.000093   0.000000   0.000023
     2_ d       0.000000   0.000000   0.000008   0.000005   0.000069   0.000199
     4_ s       0.000055   0.000000   0.000070   0.000049   0.000001   0.000007
     4_ p       0.000038   0.000008   0.000082   0.000058   0.000001   0.000005
     4_ d       0.000149   0.000037   0.000047   0.000021   0.000000   0.000003
     5_ s      -0.000001   0.000000  -0.000000   0.000012   0.000000  -0.000091
     5_ p       0.000000  -0.000000   0.000018   0.000028   0.000106   0.000049
     6_ s      -0.000105   0.000000   0.000037   0.000024   0.000000   0.000004
     6_ p      -0.000008   0.000141  -0.000008  -0.000005   0.000000   0.000000
     7_ s       0.000132  -0.000000  -0.000024  -0.000021  -0.000000  -0.000002
     7_ p       0.000214   0.000232   0.000006   0.000008   0.000001   0.000003
     8_ s       0.000001   0.000000  -0.000002   0.000015   0.000000   0.000020
     8_ p      -0.000001   0.000000   0.000024   0.000030   0.000099   0.000018
 
   ao class      43a        44a        45a        46a        47a        48a  
     2_ s       0.000005   0.000001   0.000004   0.000006   0.000027  -0.000001
     2_ p       0.000064   0.000010   0.000061   0.000038   0.000024   0.000001
     2_ d       0.000045   0.000008   0.000036   0.000103   0.000062   0.000001
     4_ s       0.000001   0.000000   0.000002   0.000002  -0.000000   0.000002
     4_ p      -0.000011   0.000015  -0.000001   0.000004   0.000002  -0.000006
     4_ d       0.000033   0.000038   0.000023   0.000004   0.000001   0.000124
     5_ s      -0.000003  -0.000001   0.000003   0.000006   0.000026   0.000000
     5_ p       0.000005   0.000002   0.000009  -0.000001   0.000045   0.000000
     6_ s       0.000014   0.000003  -0.000001  -0.000001  -0.000002   0.000001
     6_ p       0.000030   0.000077   0.000036   0.000010   0.000001   0.000019
     7_ s       0.000012   0.000003   0.000017   0.000001   0.000002   0.000006
     7_ p       0.000041   0.000064   0.000036   0.000006   0.000001  -0.000002
     8_ s      -0.000010   0.000001   0.000001   0.000020  -0.000018  -0.000002
     8_ p       0.000010   0.000004   0.000001   0.000016   0.000003   0.000002
 
   ao class      49a        50a        51a        52a        53a        54a  
     2_ s      -0.000000  -0.000010  -0.000001  -0.000000  -0.000008  -0.000000
     2_ p       0.000006  -0.000044  -0.000001  -0.000001  -0.000025  -0.000000
     2_ d       0.000089   0.000011   0.000003   0.000000   0.000074   0.000001
     4_ s       0.000001   0.000019  -0.000001   0.000000   0.000000   0.000000
     4_ p      -0.000002   0.000016  -0.000041   0.000019  -0.000001   0.000001
     4_ d       0.000003   0.000003   0.000045   0.000046   0.000002   0.000072
     5_ s      -0.000000   0.000000   0.000001   0.000000   0.000070   0.000001
     5_ p      -0.000002  -0.000003   0.000000  -0.000000  -0.000035  -0.000000
     6_ s       0.000003   0.000002   0.000037   0.000001   0.000001   0.000000
     6_ p       0.000001   0.000003  -0.000000   0.000006  -0.000001  -0.000003
     7_ s       0.000000  -0.000001   0.000050   0.000000   0.000001   0.000000
     7_ p      -0.000000   0.000002  -0.000008   0.000005   0.000001  -0.000003
     8_ s       0.000002   0.000122   0.000009   0.000002   0.000004   0.000000
     8_ p       0.000006  -0.000022  -0.000001  -0.000001  -0.000015  -0.000000
 
   ao class      55a        56a        57a        58a        59a        60a  
     2_ s      -0.000049  -0.000009   0.000000   0.000001   0.000004   0.000000
     2_ p      -0.000014  -0.000001  -0.000000   0.000000   0.000002  -0.000001
     2_ d       0.000043   0.000020   0.000000   0.000002   0.000004  -0.000002
     4_ s       0.000001   0.000004  -0.000005  -0.000001   0.000001   0.000001
     4_ p       0.000000   0.000003  -0.000002  -0.000000   0.000001   0.000007
     4_ d      -0.000000   0.000001   0.000040   0.000012   0.000004   0.000009
     5_ s       0.000119  -0.000002   0.000001   0.000002   0.000012   0.000001
     5_ p      -0.000031  -0.000002  -0.000001  -0.000001  -0.000003   0.000003
     6_ s      -0.000000   0.000001   0.000000   0.000019   0.000002   0.000000
     6_ p       0.000000   0.000006   0.000003   0.000001  -0.000000   0.000001
     7_ s       0.000000  -0.000000   0.000010   0.000010   0.000000   0.000000
     7_ p       0.000001   0.000007   0.000001   0.000000   0.000001   0.000001
     8_ s      -0.000002   0.000026   0.000000  -0.000000   0.000014   0.000000
     8_ p      -0.000005  -0.000001  -0.000000   0.000001  -0.000006   0.000014
 
   ao class      61a        62a        63a        64a        65a        66a  
     2_ s       0.000000  -0.000000  -0.000000  -0.000002   0.000000  -0.000000
     2_ p      -0.000000  -0.000001  -0.000000  -0.000000  -0.000003   0.000000
     2_ d      -0.000004  -0.000000  -0.000012   0.000002  -0.000000   0.000000
     4_ s       0.000000  -0.000026  -0.000001  -0.000000  -0.000000  -0.000001
     4_ p       0.000004  -0.000006   0.000000  -0.000000   0.000000  -0.000000
     4_ d       0.000006   0.000000   0.000000   0.000000  -0.000000  -0.000000
     5_ s       0.000000  -0.000000   0.000000   0.000003  -0.000002  -0.000000
     5_ p       0.000011   0.000000   0.000013  -0.000001   0.000010   0.000001
     6_ s       0.000000   0.000017   0.000001  -0.000000   0.000001  -0.000001
     6_ p       0.000001   0.000001  -0.000000  -0.000000   0.000000   0.000005
     7_ s      -0.000000   0.000042  -0.000000   0.000000   0.000000   0.000001
     7_ p       0.000001  -0.000004  -0.000000   0.000000   0.000001   0.000007
     8_ s      -0.000000   0.000000   0.000000  -0.000007  -0.000001   0.000000
     8_ p       0.000010   0.000000   0.000019   0.000022   0.000007   0.000000
 
   ao class      67a        68a        69a        70a        71a        72a  
     2_ s       0.000000   0.000000  -0.000000   0.000000   0.000000  -0.000000
     2_ p       0.000001  -0.000001  -0.000001   0.000003   0.000000   0.000002
     2_ d       0.000000  -0.000000  -0.000002  -0.000000  -0.000000  -0.000000
     4_ s      -0.000000   0.000000  -0.000000   0.000000  -0.000000   0.000000
     4_ p       0.000000   0.000000   0.000000   0.000002   0.000000   0.000001
     4_ d      -0.000000  -0.000000  -0.000000  -0.000000  -0.000000  -0.000000
     5_ s       0.000000  -0.000001  -0.000003   0.000000   0.000000  -0.000000
     5_ p       0.000005   0.000004   0.000013   0.000001   0.000000   0.000000
     6_ s      -0.000001  -0.000001  -0.000000  -0.000002   0.000000  -0.000003
     6_ p       0.000003   0.000004   0.000000   0.000002   0.000004   0.000005
     7_ s      -0.000001  -0.000001  -0.000000   0.000000  -0.000000  -0.000001
     7_ p       0.000001   0.000005   0.000001  -0.000000   0.000001   0.000002
     8_ s      -0.000000  -0.000000   0.000001  -0.000000  -0.000000  -0.000000
     8_ p       0.000002   0.000000   0.000002   0.000000   0.000000   0.000000
 
   ao class      73a        74a        75a        76a        77a        78a  
     2_ s      -0.000000   0.000000   0.000000   0.000001   0.000002   0.000000
     2_ p       0.000000   0.000000   0.000000   0.000003   0.000002   0.000000
     2_ d      -0.000000  -0.000000  -0.000000  -0.000000  -0.000000  -0.000000
     4_ s       0.000000   0.000000   0.000000  -0.000000  -0.000000   0.000002
     4_ p       0.000002   0.000001   0.000003  -0.000000   0.000000   0.000001
     4_ d      -0.000000  -0.000000  -0.000000  -0.000000  -0.000000  -0.000000
     5_ s      -0.000000   0.000000   0.000000   0.000001  -0.000004  -0.000000
     5_ p       0.000000  -0.000000   0.000000   0.000000   0.000001   0.000000
     6_ s       0.000000   0.000000  -0.000000  -0.000000   0.000000  -0.000000
     6_ p       0.000000   0.000000   0.000000   0.000000  -0.000000   0.000000
     7_ s      -0.000000  -0.000001  -0.000000   0.000000  -0.000000  -0.000002
     7_ p       0.000002   0.000003   0.000000   0.000000   0.000000   0.000000
     8_ s      -0.000000  -0.000001  -0.000000  -0.000003   0.000001  -0.000000
     8_ p       0.000000   0.000000  -0.000000   0.000000   0.000000  -0.000000
 
   ao class      79a        80a        81a        82a  
     2_ s      -0.000000   0.000000  -0.000000  -0.000000
     2_ p      -0.000000  -0.000000  -0.000000  -0.000000
     2_ d      -0.000000  -0.000000  -0.000000  -0.000000
     4_ s       0.000000   0.000000  -0.000000  -0.000000
     4_ p      -0.000000  -0.000000  -0.000000  -0.000000
     4_ d      -0.000000  -0.000000  -0.000000  -0.000000
     5_ s      -0.000000   0.000000   0.000000   0.000000
     5_ p       0.000000   0.000000  -0.000000  -0.000000
     6_ s       0.000000   0.000000   0.000000   0.000000
     6_ p      -0.000000  -0.000000  -0.000000  -0.000000
     7_ s       0.000000   0.000000   0.000000   0.000000
     7_ p      -0.000000  -0.000000  -0.000000  -0.000000
     8_ s       0.000001   0.000000  -0.000000  -0.000000
     8_ p      -0.000000   0.000000  -0.000000  -0.000000


                        gross atomic populations
     ao            2_         4_         5_         6_         7_         8_
      s         3.319972   3.821723   1.086491   0.577524   0.728066   1.054157
      p         2.464480   2.727874   0.029302   0.044599   0.058808   0.023849
      d         0.014946   0.048208   0.000000   0.000000   0.000000   0.000000
    total       5.799398   6.597805   1.115793   0.622123   0.786874   1.078006
 

 Total number of electrons:   16.00000000

 accstate=                     4
 accpdens=                     4
logrecs(*)=   1   2logvrecs(*)=   1   2
 item #                     4 suffix=:.trd1to2:
 computing final density
 =========== Executing IN-CORE method ==========
1e-transition density (   1,   2)
--------------------------------------------------------------------------------
  1e-transition density (root #    1 -> root #   2)
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
  1e-transition density (root #    1 -> root #   2)
--------------------------------------------------------------------------------
================================================================================
   DYZ=     242  DYX=     503  DYW=     590
   D0Z=     147  D0Y=     928  D0X=     875  D0W=    1184
  DDZI=     189 DDYI=     692 DDXI=     675 DDWI=     866
  DDZE=       0 DDYE=     104 DDXE=     105 DDWE=     140
================================================================================
2e-transition density (   1,   2)
--------------------------------------------------------------------------------
  2e-transition density (root #    1 -> root #   2)
--------------------------------------------------------------------------------
 call to prpd23 iwx=1
 call to prpd23 iwx=2
 starting prpd24 .... 
d2temp:fourex offset=     0 N= 59550 init,final=   T   Fmode=   1
d2temp:fourex offset= 59550 N= 59292 init,final=   F   Fmode=   1
d2temp:fourex offset=118842 N= 59913 init,final=   F   Fmode=   1
d2temp:fourex offset=178755 N= 59058 init,final=   F   Fmode=   1
d2temp:fourex offset=237813 N= 59667 init,final=   F   Fmode=   1
d2temp:fourex offset=297480 N= 59730 init,final=   F   Fmode=   1
d2temp:fourex offset=357210 N= 59871 init,final=   F   Fmode=   1
d2temp:fourex offset=417081 N= 59544 init,final=   F   Fmode=   1
d2temp:fourex offset=476625 N= 59835 init,final=   F   Fmode=   1
d2temp:fourex offset=536460 N= 59460 init,final=   F   Fmode=   1
d2temp:fourex offset=595920 N= 59343 init,final=   F   Fmode=   1
d2temp:fourex offset=655263 N= 58947 init,final=   F   Fmode=   1
d2temp:fourex offset=714210 N= 58800 init,final=   F   Fmode=   1
d2temp:fourex offset=773010 N= 58800 init,final=   F   Fmode=   1
d2temp:fourex offset=831810 N= 58845 init,final=   F   Fmode=   1
d2temp:fourex offset=890655 N= 59778 init,final=   F   Fmode=   1
d2temp:fourex offset=950433 N= 57120 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58182 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59445 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59955 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58695 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58620 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58140 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56676 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59124 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56793 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59799 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 55743 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59250 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58128 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57915 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 55545 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58545 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 55800 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58632 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57285 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59235 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59400 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58692 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57828 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57120 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59007 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56460 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56628 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59040 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57795 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59907 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58548 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 54660 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59535 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57525 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58692 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59958 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56793 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 54714 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 60000 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56403 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56658 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59787 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59556 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56082 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56250 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59892 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58140 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 55935 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58860 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N=     0 init,final=   F   Tmode=   1
 starting prpd24 .... 
d2temp:fourex offset=     0 N= 59550 init,final=   T   Fmode=   1
d2temp:fourex offset= 59550 N= 59292 init,final=   F   Fmode=   1
d2temp:fourex offset=118842 N= 59913 init,final=   F   Fmode=   1
d2temp:fourex offset=178755 N= 59058 init,final=   F   Fmode=   1
d2temp:fourex offset=237813 N= 59667 init,final=   F   Fmode=   1
d2temp:fourex offset=297480 N= 59730 init,final=   F   Fmode=   1
d2temp:fourex offset=357210 N= 59871 init,final=   F   Fmode=   1
d2temp:fourex offset=417081 N= 59544 init,final=   F   Fmode=   1
d2temp:fourex offset=476625 N= 59835 init,final=   F   Fmode=   1
d2temp:fourex offset=536460 N= 59460 init,final=   F   Fmode=   1
d2temp:fourex offset=595920 N= 59343 init,final=   F   Fmode=   1
d2temp:fourex offset=655263 N= 58947 init,final=   F   Fmode=   1
d2temp:fourex offset=714210 N= 58800 init,final=   F   Fmode=   1
d2temp:fourex offset=773010 N= 58800 init,final=   F   Fmode=   1
d2temp:fourex offset=831810 N= 58845 init,final=   F   Fmode=   1
d2temp:fourex offset=890655 N= 59778 init,final=   F   Fmode=   1
d2temp:fourex offset=950433 N= 57120 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58182 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59445 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59955 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58695 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58620 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58140 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56676 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59124 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56793 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59799 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 55743 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59250 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58128 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57915 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 55545 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58545 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 55800 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58632 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57285 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59235 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59400 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58692 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57828 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57120 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59007 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56460 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56628 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59040 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57795 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59907 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58548 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 54660 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59535 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57525 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58692 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59958 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56793 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 54714 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 60000 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56403 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56658 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59787 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59556 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56082 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56250 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59892 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58140 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 55935 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58860 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N=     0 init,final=   F   Tmode=   1
================================================================================
   DYZ=     242  DYX=    7351  DYW=    8696
   D0Z=     147  D0Y=    8560  D0X=    7763  D0W=   10896
  DDZI=     735 DDYI=    2660 DDXI=    2520 DDWI=    3140
  DDZE=       0 DDYE=       0 DDXE=       0 DDWE=       0
================================================================================
  entering dump_fourext .... 
d2temp:fourex offset=     1 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset= 60001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=120001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=180001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=240001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=300001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=360001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=420001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=480001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=540001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=600001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=660001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=720001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=780001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=840001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=900001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=960001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
 maximum diagonal element=  2.061512650103098E-003
 accstate=                     4
 accpdens=                     4
logrecs(*)=   1   3logvrecs(*)=   1   3
 item #                     5 suffix=:.trd1to3:
 computing final density
 =========== Executing IN-CORE method ==========
1e-transition density (   1,   3)
--------------------------------------------------------------------------------
  1e-transition density (root #    1 -> root #   3)
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
  1e-transition density (root #    1 -> root #   3)
--------------------------------------------------------------------------------
================================================================================
   DYZ=     242  DYX=     503  DYW=     590
   D0Z=     147  D0Y=     928  D0X=     875  D0W=    1184
  DDZI=     189 DDYI=     692 DDXI=     675 DDWI=     866
  DDZE=       0 DDYE=     104 DDXE=     105 DDWE=     140
================================================================================
2e-transition density (   1,   3)
--------------------------------------------------------------------------------
  2e-transition density (root #    1 -> root #   3)
--------------------------------------------------------------------------------
 call to prpd23 iwx=1
 call to prpd23 iwx=2
 starting prpd24 .... 
d2temp:fourex offset=     0 N= 59550 init,final=   T   Fmode=   1
d2temp:fourex offset= 59550 N= 59292 init,final=   F   Fmode=   1
d2temp:fourex offset=118842 N= 59913 init,final=   F   Fmode=   1
d2temp:fourex offset=178755 N= 59058 init,final=   F   Fmode=   1
d2temp:fourex offset=237813 N= 59667 init,final=   F   Fmode=   1
d2temp:fourex offset=297480 N= 59730 init,final=   F   Fmode=   1
d2temp:fourex offset=357210 N= 59871 init,final=   F   Fmode=   1
d2temp:fourex offset=417081 N= 59544 init,final=   F   Fmode=   1
d2temp:fourex offset=476625 N= 59835 init,final=   F   Fmode=   1
d2temp:fourex offset=536460 N= 59460 init,final=   F   Fmode=   1
d2temp:fourex offset=595920 N= 59343 init,final=   F   Fmode=   1
d2temp:fourex offset=655263 N= 58947 init,final=   F   Fmode=   1
d2temp:fourex offset=714210 N= 58800 init,final=   F   Fmode=   1
d2temp:fourex offset=773010 N= 58800 init,final=   F   Fmode=   1
d2temp:fourex offset=831810 N= 58845 init,final=   F   Fmode=   1
d2temp:fourex offset=890655 N= 59778 init,final=   F   Fmode=   1
d2temp:fourex offset=950433 N= 57120 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58182 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59445 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59955 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58695 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58620 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58140 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56676 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59124 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56793 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59799 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 55743 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59250 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58128 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57915 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 55545 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58545 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 55800 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58632 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57285 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59235 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59400 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58692 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57828 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57120 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59007 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56460 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56628 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59040 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57795 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59907 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58548 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 54660 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59535 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57525 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58692 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59958 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56793 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 54714 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 60000 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56403 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56658 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59787 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59556 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56082 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56250 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59892 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58140 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 55935 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58860 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N=     0 init,final=   F   Tmode=   1
 starting prpd24 .... 
d2temp:fourex offset=     0 N= 59550 init,final=   T   Fmode=   1
d2temp:fourex offset= 59550 N= 59292 init,final=   F   Fmode=   1
d2temp:fourex offset=118842 N= 59913 init,final=   F   Fmode=   1
d2temp:fourex offset=178755 N= 59058 init,final=   F   Fmode=   1
d2temp:fourex offset=237813 N= 59667 init,final=   F   Fmode=   1
d2temp:fourex offset=297480 N= 59730 init,final=   F   Fmode=   1
d2temp:fourex offset=357210 N= 59871 init,final=   F   Fmode=   1
d2temp:fourex offset=417081 N= 59544 init,final=   F   Fmode=   1
d2temp:fourex offset=476625 N= 59835 init,final=   F   Fmode=   1
d2temp:fourex offset=536460 N= 59460 init,final=   F   Fmode=   1
d2temp:fourex offset=595920 N= 59343 init,final=   F   Fmode=   1
d2temp:fourex offset=655263 N= 58947 init,final=   F   Fmode=   1
d2temp:fourex offset=714210 N= 58800 init,final=   F   Fmode=   1
d2temp:fourex offset=773010 N= 58800 init,final=   F   Fmode=   1
d2temp:fourex offset=831810 N= 58845 init,final=   F   Fmode=   1
d2temp:fourex offset=890655 N= 59778 init,final=   F   Fmode=   1
d2temp:fourex offset=950433 N= 57120 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58182 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59445 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59955 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58695 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58620 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58140 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56676 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59124 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56793 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59799 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 55743 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59250 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58128 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57915 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 55545 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58545 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 55800 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58632 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57285 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59235 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59400 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58692 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57828 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57120 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59007 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56460 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56628 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59040 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57795 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59907 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58548 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 54660 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59535 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57525 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58692 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59958 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56793 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 54714 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 60000 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56403 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56658 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59787 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59556 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56082 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56250 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59892 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58140 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 55935 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58860 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N=     0 init,final=   F   Tmode=   1
================================================================================
   DYZ=     242  DYX=    7351  DYW=    8696
   D0Z=     147  D0Y=    8560  D0X=    7763  D0W=   10896
  DDZI=     735 DDYI=    2660 DDXI=    2520 DDWI=    3140
  DDZE=       0 DDYE=       0 DDXE=       0 DDWE=       0
================================================================================
  entering dump_fourext .... 
d2temp:fourex offset=     1 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset= 60001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=120001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=180001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=240001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=300001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=360001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=420001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=480001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=540001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=600001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=660001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=720001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=780001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=840001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=900001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=960001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
 maximum diagonal element=  3.068245202642905E-002
 accstate=                     4
 accpdens=                     4
logrecs(*)=   2   3logvrecs(*)=   2   3
 item #                     6 suffix=:.trd2to3:
 computing final density
 =========== Executing IN-CORE method ==========
1e-transition density (   2,   3)
--------------------------------------------------------------------------------
  1e-transition density (root #    2 -> root #   3)
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
  1e-transition density (root #    2 -> root #   3)
--------------------------------------------------------------------------------
================================================================================
   DYZ=     242  DYX=     503  DYW=     590
   D0Z=     147  D0Y=     928  D0X=     875  D0W=    1184
  DDZI=     189 DDYI=     692 DDXI=     675 DDWI=     866
  DDZE=       0 DDYE=     104 DDXE=     105 DDWE=     140
================================================================================
2e-transition density (   2,   3)
--------------------------------------------------------------------------------
  2e-transition density (root #    2 -> root #   3)
--------------------------------------------------------------------------------
 call to prpd23 iwx=1
 call to prpd23 iwx=2
 starting prpd24 .... 
d2temp:fourex offset=     0 N= 59550 init,final=   T   Fmode=   1
d2temp:fourex offset= 59550 N= 59292 init,final=   F   Fmode=   1
d2temp:fourex offset=118842 N= 59913 init,final=   F   Fmode=   1
d2temp:fourex offset=178755 N= 59058 init,final=   F   Fmode=   1
d2temp:fourex offset=237813 N= 59667 init,final=   F   Fmode=   1
d2temp:fourex offset=297480 N= 59730 init,final=   F   Fmode=   1
d2temp:fourex offset=357210 N= 59871 init,final=   F   Fmode=   1
d2temp:fourex offset=417081 N= 59544 init,final=   F   Fmode=   1
d2temp:fourex offset=476625 N= 59835 init,final=   F   Fmode=   1
d2temp:fourex offset=536460 N= 59460 init,final=   F   Fmode=   1
d2temp:fourex offset=595920 N= 59343 init,final=   F   Fmode=   1
d2temp:fourex offset=655263 N= 58947 init,final=   F   Fmode=   1
d2temp:fourex offset=714210 N= 58800 init,final=   F   Fmode=   1
d2temp:fourex offset=773010 N= 58800 init,final=   F   Fmode=   1
d2temp:fourex offset=831810 N= 58845 init,final=   F   Fmode=   1
d2temp:fourex offset=890655 N= 59778 init,final=   F   Fmode=   1
d2temp:fourex offset=950433 N= 57120 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58182 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59445 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59955 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58695 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58620 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58140 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56676 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59124 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56793 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59799 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 55743 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59250 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58128 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57915 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 55545 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58545 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 55800 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58632 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57285 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59235 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59400 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58692 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57828 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57120 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59007 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56460 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56628 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59040 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57795 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59907 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58548 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 54660 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59535 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57525 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58692 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59958 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56793 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 54714 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 60000 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56403 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56658 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59787 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59556 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56082 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56250 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59892 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58140 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 55935 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58860 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N=     0 init,final=   F   Tmode=   1
 starting prpd24 .... 
d2temp:fourex offset=     0 N= 59550 init,final=   T   Fmode=   1
d2temp:fourex offset= 59550 N= 59292 init,final=   F   Fmode=   1
d2temp:fourex offset=118842 N= 59913 init,final=   F   Fmode=   1
d2temp:fourex offset=178755 N= 59058 init,final=   F   Fmode=   1
d2temp:fourex offset=237813 N= 59667 init,final=   F   Fmode=   1
d2temp:fourex offset=297480 N= 59730 init,final=   F   Fmode=   1
d2temp:fourex offset=357210 N= 59871 init,final=   F   Fmode=   1
d2temp:fourex offset=417081 N= 59544 init,final=   F   Fmode=   1
d2temp:fourex offset=476625 N= 59835 init,final=   F   Fmode=   1
d2temp:fourex offset=536460 N= 59460 init,final=   F   Fmode=   1
d2temp:fourex offset=595920 N= 59343 init,final=   F   Fmode=   1
d2temp:fourex offset=655263 N= 58947 init,final=   F   Fmode=   1
d2temp:fourex offset=714210 N= 58800 init,final=   F   Fmode=   1
d2temp:fourex offset=773010 N= 58800 init,final=   F   Fmode=   1
d2temp:fourex offset=831810 N= 58845 init,final=   F   Fmode=   1
d2temp:fourex offset=890655 N= 59778 init,final=   F   Fmode=   1
d2temp:fourex offset=950433 N= 57120 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58182 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59445 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59955 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58695 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58620 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58140 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56676 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59124 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56793 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59799 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 55743 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59250 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58128 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57915 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 55545 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58545 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 55800 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58632 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57285 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59235 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59400 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58692 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57828 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57120 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59007 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56460 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56628 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59040 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57795 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59907 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58548 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 54660 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59535 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 57525 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58692 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59958 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56793 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 54714 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 60000 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56403 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56658 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59787 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59556 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56082 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 56250 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 59892 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58140 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 55935 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N= 58860 init,final=   F   Fmode=   1
d2temp:fourex offset=****** N=     0 init,final=   F   Tmode=   1
================================================================================
   DYZ=     242  DYX=    7351  DYW=    8696
   D0Z=     147  D0Y=    8560  D0X=    7763  D0W=   10896
  DDZI=     735 DDYI=    2660 DDXI=    2520 DDWI=    3140
  DDZE=       0 DDYE=       0 DDXE=       0 DDWE=       0
================================================================================
  entering dump_fourext .... 
d2temp:fourex offset=     1 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset= 60001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=120001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=180001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=240001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=300001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=360001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=420001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=480001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=540001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=600001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=660001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=720001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=780001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=840001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=900001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=960001 N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
d2temp:fourex offset=****** N= 60000 init,final=   T   Fmode=   2
 maximum diagonal element=  6.334476548383650E-002
 DA ...
