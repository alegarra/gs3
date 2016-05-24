C----- SUBROUTINE GENQMD
C****************************************************************          1.
C****************************************************************          2.
C**********    GENQMD ..... QUOT MIN DEGREE ORDERING    *********          3.
C****************************************************************          4.
C****************************************************************          5.
C                                                                          6.
C     PURPOSE - THIS ROUTINE IMPLEMENTS THE MINIMUM DEGREE                 7.
C        ALGORITHM.  IT MAKES USE OF THE IMPLICIT REPRESENT-               8.
C        ATION OF THE ELIMINATION GRAPHS BY QUOTIENT GRAPHS,               9.
C        AND THE NOTION OF INDISTINGUISHABLE NODES.                       10.
C        CAUTION - THE ADJACENCY VECTOR ADJNCY WILL BE                    11.
C        DESTROYED.                                                       12.
C                                                                         13.
C     INPUT PARAMETERS -                                                  14.
C        NEQNS - NUMBER OF EQUATIONS.                                     15.
C        (XADJ, ADJNCY) - THE ADJACENCY STRUCTURE.                        16.
C                                                                         17.
C     OUTPUT PARAMETERS -                                                 18.
C        PERM - THE MINIMUM DEGREE ORDERING.                              19.
C        INVP - THE INVERSE OF PERM.                                      20.
C                                                                         21.
C     WORKING PARAMETERS -                                                22.
C        DEG - THE DEGREE VECTOR. DEG(I) IS NEGATIVE MEANS                23.
C               NODE I HAS BEEN NUMBERED.                                 24.
C        MARKER - A MARKER VECTOR, WHERE MARKER(I) IS                     25.
C               NEGATIVE MEANS NODE I HAS BEEN MERGED WITH                26.
C               ANOTHER NODE AND THUS CAN BE IGNORED.                     27.
C        RCHSET - VECTOR USED FOR THE REACHABLE SET.                      28.
C        NBRHD - VECTOR USED FOR THE NEIGHBORHOOD SET.                    29.
C        QSIZE - VECTOR USED TO STORE THE SIZE OF                         30.
C               INDISTINGUISHABLE SUPERNODES.                             31.
C        QLINK - VECTOR TO STORE INDISTINGUISHABLE NODES,                 32.
C               I, QLINK(I), QLINK(QLINK(I)) ... ARE THE                  33.
C               MEMBERS OF THE SUPERNODE REPRESENTED BY I.                34.
C                                                                         35.
C     PROGRAM SUBROUTINES -                                               36.
C        QMDRCH, QMDQT, QMDUPD.                                           37.
C                                                                         38.
C****************************************************************         39.
C                                                                         40.
C                                                                         41.
      SUBROUTINE  GENQMD ( NEQNS, XADJ, ADJNCY, PERM, INVP, DEG,          42.
     1                     MARKER, RCHSET, NBRHD, QSIZE, QLINK,           43.
     1                     NOFSUB )                                       44.
C                                                                         45.
C****************************************************************         46.
C                                                                         47.
         INTEGER ADJNCY(1), PERM(1), INVP(1), DEG(1), MARKER(1),          48.
     1           RCHSET(1), NBRHD(1), QSIZE(1), QLINK(1)                  49.
         INTEGER XADJ(1), INODE, IP, IRCH, J, MINDEG, NDEG,               50.
     1           NEQNS, NHDSZE, NODE, NOFSUB, NP, NUM, NUMP1,             51.
     1           NXNODE, RCHSZE, SEARCH, THRESH                           52.
C                                                                         53.
C****************************************************************         54.
C                                                                         55.
C        -----------------------------------------------------            56.
C        INITIALIZE DEGREE VECTOR AND OTHER WORKING VARIABLES.            57.
C        -----------------------------------------------------            58.
         MINDEG = NEQNS                                                   59.
         NOFSUB = 0                                                       60.
         DO 100 NODE = 1, NEQNS                                           61.
            PERM(NODE) = NODE                                             62.
            INVP(NODE) = NODE                                             63.
            MARKER(NODE) = 0                                              64.
            QSIZE(NODE)  = 1                                              65.
            QLINK(NODE)  = 0                                              66.
            NDEG = XADJ(NODE+1) - XADJ(NODE)                              67.
            DEG(NODE) = NDEG                                              68.
            IF ( NDEG .LT. MINDEG )  MINDEG = NDEG                        69.
  100    CONTINUE                                                         70.
         NUM = 0                                                          71.
C        -----------------------------------------------------            72.
C        PERFORM THRESHOLD SEARCH TO GET A NODE OF MIN DEGREE.            73.
C        VARIABLE SEARCH POINTS TO WHERE SEARCH SHOULD START.             74.
C        -----------------------------------------------------            75.
  200    SEARCH = 1                                                       76.
            THRESH = MINDEG                                               77.
            MINDEG = NEQNS                                                78.
  300       NUMP1 = NUM + 1                                               79.
               IF ( NUMP1 .GT. SEARCH )  SEARCH = NUMP1                   80.
               DO 400 J = SEARCH, NEQNS                                   81.
                  NODE = PERM(J)                                          82.
                  IF ( MARKER(NODE) .LT. 0 )  GOTO 400                    83.
                     NDEG = DEG(NODE)                                     84.
                     IF ( NDEG .LE. THRESH )  GO TO 500                   85.
                     IF ( NDEG .LT. MINDEG )  MINDEG =  NDEG              86.
  400          CONTINUE                                                   87.
            GO TO 200                                                     88.
C           ---------------------------------------------------           89.
C           NODE HAS MINIMUM DEGREE. FIND ITS REACHABLE SETS BY           90.
C           CALLING QMDRCH.                                               91.
C           ---------------------------------------------------           92.
  500       SEARCH = J                                                    93.
            NOFSUB = NOFSUB + DEG(NODE)                                   94.
            MARKER(NODE) = 1                                              95.
            CALL QMDRCH (NODE, XADJ, ADJNCY, DEG, MARKER,                 96.
     1                   RCHSZE, RCHSET, NHDSZE, NBRHD )                  97.
C           ------------------------------------------------              98.
C           ELIMINATE ALL NODES INDISTINGUISHABLE FROM NODE.              99.
C           THEY ARE GIVEN BY NODE, QLINK(NODE), ....                    100.
C           ------------------------------------------------             101.
            NXNODE = NODE                                                102.
  600       NUM = NUM + 1                                                103.
               NP  = INVP(NXNODE)                                        104.
               IP  = PERM(NUM)                                           105.
               PERM(NP) = IP                                             106.
               INVP(IP) = NP                                             107.
               PERM(NUM) = NXNODE                                        108.
               INVP(NXNODE) = NUM                                        109.
               DEG(NXNODE) = - 1                                         110.
               NXNODE = QLINK(NXNODE)                                    111.
            IF (NXNODE .GT. 0) GOTO 600                                  112.
C                                                                        113.
            IF ( RCHSZE .LE. 0 )  GO TO 800                              114.
C              ------------------------------------------------          115.
C              UPDATE THE DEGREES OF THE NODES IN THE REACHABLE          116.
C              SET AND IDENTIFY INDISTINGUISHABLE NODES.                 117.
C              ------------------------------------------------          118.
               CALL  QMDUPD ( XADJ, ADJNCY, RCHSZE, RCHSET, DEG,         119.
     1                        QSIZE, QLINK, MARKER, RCHSET(RCHSZE+1),    120.
     1                        NBRHD(NHDSZE+1) )                          121.
C              -------------------------------------------               122.
C              RESET MARKER VALUE OF NODES IN REACH SET.                 123.
C              UPDATE THRESHOLD VALUE FOR CYCLIC SEARCH.                 124.
C              ALSO CALL QMDQT TO FORM NEW QUOTIENT GRAPH.               125.
C              -------------------------------------------               126.
               MARKER(NODE) = 0                                          127.
               DO 700 IRCH = 1, RCHSZE                                   128.
                  INODE = RCHSET(IRCH)                                   129.
                  IF ( MARKER(INODE) .LT. 0 )  GOTO 700                  130.
                     MARKER(INODE) = 0                                   131.
                     NDEG = DEG(INODE)                                   132.
                     IF ( NDEG .LT. MINDEG )  MINDEG = NDEG              133.
                     IF ( NDEG .GT. THRESH )  GOTO 700                   134.
                        MINDEG = THRESH                                  135.
                        THRESH = NDEG                                    136.
                        SEARCH = INVP(INODE)                             137.
  700          CONTINUE                                                  138.
               IF ( NHDSZE .GT. 0 )  CALL  QMDQT ( NODE, XADJ,           139.
     1            ADJNCY, MARKER, RCHSZE, RCHSET, NBRHD )                140.
  800    IF ( NUM .LT. NEQNS )  GO TO 300                                141.
         RETURN                                                          142.
      END                                                                143.
C----- SUBROUTINE QMDUPD
C****************************************************************          1.
C****************************************************************          2.
C**********     QMDUPD ..... QUOT MIN DEG UPDATE      ***********          3.
C****************************************************************          4.
C****************************************************************          5.
C                                                                          6.
C     PURPOSE - THIS ROUTINE PERFORMS DEGREE UPDATE FOR A SET              7.
C        OF NODES IN THE MINIMUM DEGREE ALGORITHM.                         8.
C                                                                          9.
C     INPUT PARAMETERS -                                                  10.
C        (XADJ, ADJNCY) - THE ADJACENCY STRUCTURE.                        11.
C        (NLIST, LIST) - THE LIST OF NODES WHOSE DEGREE HAS TO            12.
C               BE UPDATED.                                               13.
C                                                                         14.
C     UPDATED PARAMETERS -                                                15.
C        DEG - THE DEGREE VECTOR.                                         16.
C        QSIZE - SIZE OF INDISTINGUISHABLE SUPERNODES.                    17.
C        QLINK - LINKED LIST FOR INDISTINGUISHABLE NODES.                 18.
C        MARKER - USED TO MARK THOSE NODES IN REACH/NBRHD SETS.           19.
C                                                                         20.
C     WORKING PARAMETERS -                                                21.
C        RCHSET - THE REACHABLE SET.                                      22.
C        NBRHD -  THE NEIGHBORHOOD SET.                                   23.
C                                                                         24.
C     PROGRAM SUBROUTINES -                                               25.
C        QMDMRG.                                                          26.
C                                                                         27.
C****************************************************************         28.
C                                                                         29.
      SUBROUTINE  QMDUPD ( XADJ, ADJNCY, NLIST, LIST, DEG,                30.
     1                     QSIZE, QLINK, MARKER, RCHSET, NBRHD )          31.
C                                                                         32.
C****************************************************************         33.
C                                                                         34.
         INTEGER  ADJNCY(1), LIST(1), DEG(1), MARKER(1),                  35.
     1            RCHSET(1), NBRHD(1), QSIZE(1), QLINK(1)                 36.
         INTEGER  XADJ(1), DEG0, DEG1, IL, INHD, INODE, IRCH,             37.
     1            J, JSTRT, JSTOP, MARK, NABOR, NHDSZE, NLIST,            38.
     1            NODE, RCHSZE, ROOT                                      39.
C                                                                         40.
C****************************************************************         41.
C                                                                         42.
C        ------------------------------------------------                 43.
C        FIND ALL ELIMINATED SUPERNODES THAT ARE ADJACENT                 44.
C        TO SOME NODES IN THE GIVEN LIST. PUT THEM INTO                   45.
C        (NHDSZE, NBRHD). DEG0 CONTAINS THE NUMBER OF                     46.
C        NODES IN THE LIST.                                               47.
C        ------------------------------------------------                 48.
         IF ( NLIST .LE. 0 )  RETURN                                      49.
         DEG0 = 0                                                         50.
         NHDSZE = 0                                                       51.
         DO 200 IL = 1, NLIST                                             52.
            NODE = LIST(IL)                                               53.
            DEG0 = DEG0 + QSIZE(NODE)                                     54.
            JSTRT = XADJ(NODE)                                            55.
            JSTOP = XADJ(NODE+1) - 1                                      56.
            DO 100 J = JSTRT, JSTOP                                       57.
               NABOR = ADJNCY(J)                                          58.
               IF ( MARKER(NABOR) .NE. 0  .OR.                            59.
     1              DEG(NABOR) .GE. 0 )  GO TO 100                        60.
                  MARKER(NABOR) = - 1                                     61.
                  NHDSZE = NHDSZE + 1                                     62.
                  NBRHD(NHDSZE) = NABOR                                   63.
  100       CONTINUE                                                      64.
  200    CONTINUE                                                         65.
C        --------------------------------------------                     66.
C        MERGE INDISTINGUISHABLE NODES IN THE LIST BY                     67.
C        CALLING THE SUBROUTINE QMDMRG.                                   68.
C        --------------------------------------------                     69.
         IF ( NHDSZE .GT. 0 )                                             70.
     1      CALL  QMDMRG ( XADJ, ADJNCY, DEG, QSIZE, QLINK,               71.
     1                     MARKER, DEG0, NHDSZE, NBRHD, RCHSET,           72.
     1                     NBRHD(NHDSZE+1) )                              73.
C        ----------------------------------------------------             74.
C        FIND THE NEW DEGREES OF THE NODES THAT HAVE NOT BEEN             75.
C        MERGED.                                                          76.
C        ----------------------------------------------------             77.
         DO 600 IL = 1, NLIST                                             78.
            NODE = LIST(IL)                                               79.
            MARK = MARKER(NODE)                                           80.
            IF ( MARK .GT. 1  .OR.  MARK .LT. 0 )  GO TO 600              81.
               MARKER(NODE) = 2                                           82.
               CALL  QMDRCH ( NODE, XADJ, ADJNCY, DEG, MARKER,            83.
     1                        RCHSZE, RCHSET, NHDSZE, NBRHD )             84.
               DEG1 = DEG0                                                85.
               IF ( RCHSZE .LE. 0 )  GO TO 400                            86.
                  DO 300 IRCH = 1, RCHSZE                                 87.
                     INODE = RCHSET(IRCH)                                 88.
                     DEG1 = DEG1 + QSIZE(INODE)                           89.
                     MARKER(INODE) = 0                                    90.
  300             CONTINUE                                                91.
  400          DEG(NODE) = DEG1 - 1                                       92.
               IF ( NHDSZE .LE. 0 )  GO TO 600                            93.
                  DO 500 INHD = 1, NHDSZE                                 94.
                     INODE = NBRHD(INHD)                                  95.
                     MARKER(INODE) = 0                                    96.
  500             CONTINUE                                                97.
  600    CONTINUE                                                         98.
         RETURN                                                           99.
      END                                                                100.
C----- SUBROUTINE QMDQT
C*************************************************************             1.
C*************************************************************             2.
C*******     QMDQT  ..... QUOT MIN DEG QUOT TRANSFORM  *******             3.
C*************************************************************             4.
C*************************************************************             5.
C                                                                          6.
C     PURPOSE - THIS SUBROUTINE PERFORMS THE QUOTIENT GRAPH                7.
C        TRANSFORMATION AFTER A NODE HAS BEEN ELIMINATED.                  8.
C                                                                          9.
C     INPUT PARAMETERS -                                                  10.
C        ROOT - THE NODE JUST ELIMINATED. IT BECOMES THE                  11.
C               REPRESENTATIVE OF THE NEW SUPERNODE.                      12.
C        (XADJ, ADJNCY) - THE ADJACENCY STRUCTURE.                        13.
C        (RCHSZE, RCHSET) - THE REACHABLE SET OF ROOT IN THE              14.
C               OLD QUOTIENT GRAPH.                                       15.
C        NBRHD - THE NEIGHBORHOOD SET WHICH WILL BE MERGED                16.
C               WITH ROOT TO FORM THE NEW SUPERNODE.                      17.
C        MARKER - THE MARKER VECTOR.                                      18.
C                                                                         19.
C     UPDATED PARAMETER -                                                 20.
C        ADJNCY - BECOMES THE ADJNCY OF THE QUOTIENT GRAPH.               21.
C                                                                         22.
C*************************************************************            23.
C                                                                         24.
      SUBROUTINE  QMDQT ( ROOT, XADJ, ADJNCY, MARKER,                     25.
     1                    RCHSZE, RCHSET, NBRHD )                         26.
C                                                                         27.
C*************************************************************            28.
C                                                                         29.
         INTEGER ADJNCY(1), MARKER(1), RCHSET(1), NBRHD(1)                30.
         INTEGER XADJ(1), INHD, IRCH, J, JSTRT, JSTOP, LINK,              31.
     1           NABOR, NODE, RCHSZE, ROOT                                32.
C                                                                         33.
C*************************************************************            34.
C                                                                         35.
         IRCH = 0                                                         36.
         INHD = 0                                                         37.
         NODE = ROOT                                                      38.
  100    JSTRT = XADJ(NODE)                                               39.
         JSTOP = XADJ(NODE+1) - 2                                         40.
         IF ( JSTOP .LT. JSTRT )  GO TO 300                               41.
C           ------------------------------------------------              42.
C           PLACE REACH NODES INTO THE ADJACENT LIST OF NODE              43.
C           ------------------------------------------------              44.
            DO 200 J = JSTRT, JSTOP                                       45.
               IRCH = IRCH + 1                                            46.
               ADJNCY(J) = RCHSET(IRCH)                                   47.
               IF ( IRCH .GE. RCHSZE )  GOTO 400                          48.
  200       CONTINUE                                                      49.
C        ----------------------------------------------                   50.
C        LINK TO OTHER SPACE PROVIDED BY THE NBRHD SET.                   51.
C        ----------------------------------------------                   52.
  300    LINK = ADJNCY(JSTOP+1)                                           53.
         NODE = - LINK                                                    54.
         IF ( LINK .LT. 0 )  GOTO 100                                     55.
            INHD = INHD + 1                                               56.
            NODE = NBRHD(INHD)                                            57.
            ADJNCY(JSTOP+1) = - NODE                                      58.
            GO TO 100                                                     59.
C        -------------------------------------------------------          60.
C        ALL REACHABLE NODES HAVE BEEN SAVED.  END THE ADJ LIST.          61.
C        ADD ROOT TO THE NBR LIST OF EACH NODE IN THE REACH SET.          62.
C        -------------------------------------------------------          63.
  400    ADJNCY(J+1) = 0                                                  64.
         DO 600 IRCH = 1, RCHSZE                                          65.
            NODE = RCHSET(IRCH)                                           66.
            IF ( MARKER(NODE) .LT. 0 )  GOTO 600                          67.
               JSTRT = XADJ(NODE)                                         68.
               JSTOP = XADJ(NODE+1) - 1                                   69.
               DO 500 J = JSTRT, JSTOP                                    70.
                  NABOR = ADJNCY(J)                                       71.
                  IF ( MARKER(NABOR) .GE. 0 ) GO TO 500                   72.
                     ADJNCY(J) = ROOT                                     73.
                     GOTO 600                                             74.
  500          CONTINUE                                                   75.
  600    CONTINUE                                                         76.
         RETURN                                                           77.
      END                                                                 78.
C----- SUBROUTINE QMDMRG
C****************************************************************          1.
C****************************************************************          2.
C**********     QMDMRG ..... QUOT MIN DEG MERGE       ***********          3.
C****************************************************************          4.
C****************************************************************          5.
C                                                                          6.
C     PURPOSE - THIS ROUTINE MERGES INDISTINGUISHABLE NODES IN             7.
C               THE MINIMUM DEGREE ORDERING ALGORITHM.                     8.
C               IT ALSO COMPUTES THE NEW DEGREES OF THESE                  9.
C               NEW SUPERNODES.                                           10.
C                                                                         11.
C     INPUT PARAMETERS -                                                  12.
C        (XADJ, ADJNCY) - THE ADJACENCY STRUCTURE.                        13.
C        DEG0 - THE NUMBER OF NODES IN THE GIVEN SET.                     14.
C        (NHDSZE, NBRHD) - THE SET OF ELIMINATED SUPERNODES               15.
C               ADJACENT TO SOME NODES IN THE SET.                        16.
C                                                                         17.
C     UPDATED PARAMETERS -                                                18.
C        DEG - THE DEGREE VECTOR.                                         19.
C        QSIZE - SIZE OF INDISTINGUISHABLE NODES.                         20.
C        QLINK - LINKED LIST FOR INDISTINGUISHABLE NODES.                 21.
C        MARKER - THE GIVEN SET IS GIVEN BY THOSE NODES WITH              22.
C               MARKER VALUE SET TO 1.  THOSE NODES WITH DEGREE           23.
C               UPDATED WILL HAVE MARKER VALUE SET TO 2.                  24.
C                                                                         25.
C     WORKING PARAMETERS -                                                26.
C        RCHSET - THE REACHABLE SET.                                      27.
C        OVRLP -  TEMP VECTOR TO STORE THE INTERSECTION OF TWO            28.
C               REACHABLE SETS.                                           29.
C                                                                         30.
C****************************************************************         31.
C                                                                         32.
      SUBROUTINE  QMDMRG ( XADJ, ADJNCY, DEG, QSIZE, QLINK,               33.
     1                     MARKER, DEG0, NHDSZE, NBRHD, RCHSET,           34.
     1                     OVRLP )                                        35.
C                                                                         36.
C****************************************************************         37.
C                                                                         38.
         INTEGER  ADJNCY(1), DEG(1), QSIZE(1), QLINK(1),                  39.
     1            MARKER(1), RCHSET(1), NBRHD(1), OVRLP(1)                40.
         INTEGER  XADJ(1), DEG0, DEG1, HEAD, INHD, IOV, IRCH,             41.
     1            J, JSTRT, JSTOP, LINK, LNODE, MARK, MRGSZE,             42.
     1            NABOR, NHDSZE, NODE, NOVRLP, RCHSZE, ROOT               43.
C                                                                         44.
C****************************************************************         45.
C                                                                         46.
C        ------------------                                               47.
C        INITIALIZATION ...                                               48.
C        ------------------                                               49.
         IF ( NHDSZE .LE. 0 )  RETURN                                     50.
         DO 100 INHD = 1, NHDSZE                                          51.
            ROOT = NBRHD(INHD)                                            52.
            MARKER(ROOT) = 0                                              53.
  100    CONTINUE                                                         54.
C        -------------------------------------------------                55.
C        LOOP THROUGH EACH ELIMINATED SUPERNODE IN THE SET                56.
C        (NHDSZE, NBRHD).                                                 57.
C        -------------------------------------------------                58.
         DO 1400 INHD = 1, NHDSZE                                         59.
            ROOT = NBRHD(INHD)                                            60.
            MARKER(ROOT) = - 1                                            61.
            RCHSZE = 0                                                    62.
            NOVRLP = 0                                                    63.
            DEG1   = 0                                                    64.
  200       JSTRT  = XADJ(ROOT)                                           65.
            JSTOP  = XADJ(ROOT+1) - 1                                     66.
C           ----------------------------------------------                67.
C           DETERMINE THE REACHABLE SET AND ITS INTERSECT-                68.
C           ION WITH THE INPUT REACHABLE SET.                             69.
C           ----------------------------------------------                70.
            DO 600 J = JSTRT, JSTOP                                       71.
               NABOR = ADJNCY(J)                                          72.
               ROOT  = - NABOR                                            73.
               IF (NABOR)  200, 700, 300                                  74.
C                                                                         75.
  300          MARK = MARKER(NABOR)                                       76.
               IF ( MARK ) 600, 400, 500                                  77.
  400             RCHSZE = RCHSZE + 1                                     78.
                  RCHSET(RCHSZE) = NABOR                                  79.
                  DEG1 = DEG1 + QSIZE(NABOR)                              80.
                  MARKER(NABOR) = 1                                       81.
                  GOTO 600                                                82.
  500          IF ( MARK .GT. 1 )  GOTO 600                               83.
                  NOVRLP = NOVRLP + 1                                     84.
                  OVRLP(NOVRLP) = NABOR                                   85.
                  MARKER(NABOR) = 2                                       86.
  600       CONTINUE                                                      87.
C           --------------------------------------------                  88.
C           FROM THE OVERLAPPED SET, DETERMINE THE NODES                  89.
C           THAT CAN BE MERGED TOGETHER.                                  90.
C           --------------------------------------------                  91.
  700       HEAD = 0                                                      92.
            MRGSZE = 0                                                    93.
            DO 1100 IOV = 1, NOVRLP                                       94.
               NODE = OVRLP(IOV)                                          95.
               JSTRT = XADJ(NODE)                                         96.
               JSTOP = XADJ(NODE+1) - 1                                   97.
               DO 800 J = JSTRT, JSTOP                                    98.
                  NABOR = ADJNCY(J)                                       99.
                  IF ( MARKER(NABOR) .NE. 0 )  GOTO 800                  100.
                     MARKER(NODE) = 1                                    101.
                     GOTO 1100                                           102.
  800          CONTINUE                                                  103.
C              -----------------------------------------                 104.
C              NODE BELONGS TO THE NEW MERGED SUPERNODE.                 105.
C              UPDATE THE VECTORS QLINK AND QSIZE.                       106.
C              -----------------------------------------                 107.
               MRGSZE = MRGSZE + QSIZE(NODE)                             108.
               MARKER(NODE) = - 1                                        109.
               LNODE = NODE                                              110.
  900          LINK  = QLINK(LNODE)                                      111.
               IF ( LINK .LE. 0 )  GOTO 1000                             112.
                  LNODE = LINK                                           113.
                  GOTO 900                                               114.
 1000          QLINK(LNODE) = HEAD                                       115.
               HEAD = NODE                                               116.
 1100       CONTINUE                                                     117.
            IF ( HEAD .LE. 0 )  GOTO 1200                                118.
               QSIZE(HEAD) = MRGSZE                                      119.
               DEG(HEAD) = DEG0 + DEG1 - 1                               120.
               MARKER(HEAD) = 2                                          121.
C           --------------------                                         122.
C           RESET MARKER VALUES.                                         123.
C           --------------------                                         124.
 1200       ROOT = NBRHD(INHD)                                           125.
            MARKER(ROOT) = 0                                             126.
            IF ( RCHSZE .LE. 0 )  GOTO 1400                              127.
               DO 1300 IRCH = 1, RCHSZE                                  128.
                  NODE = RCHSET(IRCH)                                    129.
                  MARKER(NODE) = 0                                       130.
 1300          CONTINUE                                                  131.
 1400    CONTINUE                                                        132.
         RETURN                                                          133.
      END                                                                134.
C----- SUBROUTINE QMDRCH
C***************************************************************           1.
C***************************************************************           2.
C*********     QMDRCH ..... QUOT MIN DEG REACH SET    **********           3.
C***************************************************************           4.
C***************************************************************           5.
C                                                                          6.
C     PURPOSE - THIS SUBROUTINE DETERMINES THE REACHABLE SET OF            7.
C        A NODE THROUGH A GIVEN SUBSET.  THE ADJACENCY STRUCTURE           8.
C        IS ASSUMED TO BE STORED IN A QUOTIENT GRAPH FORMAT.               9.
C                                                                         10.
C     INPUT PARAMETERS -                                                  11.
C        ROOT - THE GIVEN NODE NOT IN THE SUBSET.                         12.
C        (XADJ, ADJNCY) - THE ADJACENCY STRUCTURE PAIR.                   13.
C        DEG - THE DEGREE VECTOR.  DEG(I) LT 0 MEANS THE NODE             14.
C               BELONGS TO THE GIVEN SUBSET.                              15.
C                                                                         16.
C     OUTPUT PARAMETERS -                                                 17.
C        (RCHSZE, RCHSET) - THE REACHABLE SET.                            18.
C        (NHDSZE, NBRHD) - THE NEIGHBORHOOD SET.                          19.
C                                                                         20.
C     UPDATED PARAMETERS -                                                21.
C        MARKER - THE MARKER VECTOR FOR REACH AND NBRHD SETS.             22.
C               GT 0 MEANS THE NODE IS IN REACH SET.                      23.
C               LT 0 MEANS THE NODE HAS BEEN MERGED WITH                  24.
C               OTHERS IN THE QUOTIENT OR IT IS IN NBRHD SET.             25.
C                                                                         26.
C***************************************************************          27.
C                                                                         28.
      SUBROUTINE  QMDRCH ( ROOT, XADJ, ADJNCY, DEG, MARKER,               29.
     1                     RCHSZE, RCHSET, NHDSZE, NBRHD )                30.
C                                                                         31.
C***************************************************************          32.
C                                                                         33.
         INTEGER ADJNCY(1), DEG(1), MARKER(1),                            34.
     1           RCHSET(1), NBRHD(1)                                      35.
         INTEGER XADJ(1), I, ISTRT, ISTOP, J, JSTRT, JSTOP,               36.
     1           NABOR, NHDSZE, NODE, RCHSZE, ROOT                        37.
C                                                                         38.
C***************************************************************          39.
C                                                                         40.
C        -----------------------------------------                        41.
C        LOOP THROUGH THE NEIGHBORS OF ROOT IN THE                        42.
C        QUOTIENT GRAPH.                                                  43.
C        -----------------------------------------                        44.
         NHDSZE = 0                                                       45.
         RCHSZE = 0                                                       46.
         ISTRT = XADJ(ROOT)                                               47.
         ISTOP = XADJ(ROOT+1) - 1                                         48.
         IF ( ISTOP .LT. ISTRT )  RETURN                                  49.
            DO 600 I = ISTRT, ISTOP                                       50.
               NABOR =  ADJNCY(I)                                         51.
               IF ( NABOR .EQ. 0 ) RETURN                                 52.
               IF ( MARKER(NABOR) .NE. 0 )  GO TO 600                     53.
                  IF ( DEG(NABOR) .LT. 0 )     GO TO 200                  54.
C                    -------------------------------------                55.
C                    INCLUDE NABOR INTO THE REACHABLE SET.                56.
C                    -------------------------------------                57.
                     RCHSZE = RCHSZE + 1                                  58.
                     RCHSET(RCHSZE) = NABOR                               59.
                     MARKER(NABOR) = 1                                    60.
                     GO TO 600                                            61.
C                 -------------------------------------                   62.
C                 NABOR HAS BEEN ELIMINATED. FIND NODES                   63.
C                 REACHABLE FROM IT.                                      64.
C                 -------------------------------------                   65.
  200             MARKER(NABOR) = -1                                      66.
                  NHDSZE = NHDSZE +  1                                    67.
                  NBRHD(NHDSZE) = NABOR                                   68.
  300             JSTRT = XADJ(NABOR)                                     69.
                  JSTOP = XADJ(NABOR+1) - 1                               70.
                  DO 500 J = JSTRT, JSTOP                                 71.
                     NODE = ADJNCY(J)                                     72.
                     NABOR = - NODE                                       73.
                     IF (NODE) 300, 600, 400                              74.
  400                IF ( MARKER(NODE) .NE. 0 )  GO TO 500                75.
                        RCHSZE = RCHSZE + 1                               76.
                        RCHSET(RCHSZE) = NODE                             77.
                        MARKER(NODE) = 1                                  78.
  500             CONTINUE                                                79.
  600       CONTINUE                                                      80.
            RETURN                                                        81.
      END                                                                 82.
C----- SUBROUTINE SMBFCT
C****************************************************************          1.
C****************************************************************          2.
C*********     SMBFCT ..... SYMBOLIC FACTORIZATION       ********          3.
C****************************************************************          4.
C****************************************************************          5.
C                                                                          6.
C     PURPOSE - THIS ROUTINE PERFORMS SYMBOLIC FACTORIZATION               7.
C        ON A PERMUTED LINEAR SYSTEM AND IT ALSO SETS UP THE               8.
C        COMPRESSED DATA STRUCTURE FOR THE SYSTEM.                         9.
C                                                                         10.
C     INPUT PARAMETERS -                                                  11.
C        NEQNS - NUMBER OF EQUATIONS.                                     12.
C        (XADJ, ADJNCY) - THE ADJACENCY STRUCTURE.                        13.
C        (PERM, INVP) - THE PERMUTATION VECTOR AND ITS INVERSE.           14.
C                                                                         15.
C     UPDATED PARAMETERS -                                                16.
C        MAXSUB - SIZE OF THE SUBSCRIPT ARRAY NZSUB.  ON RETURN,          17.
C               IT CONTAINS THE NUMBER OF SUBSCRIPTS USED                 18.
C                                                                         19.
C     OUTPUT PARAMETERS -                                                 20.
C        XLNZ - INDEX INTO THE NONZERO STORAGE VECTOR LNZ.                21.
C        (XNZSUB, NZSUB) - THE COMPRESSED SUBSCRIPT VECTORS.              22.
C        MAXLNZ - THE NUMBER OF NONZEROS FOUND.                           23.
C        FLAG - ERROR FLAG.  POSITIVE VALUE INDICATES THAT.               24.
C               NZSUB ARRAY IS TOO SMALL.                                 25.
C                                                                         26.
C     WORKING PARAMETERS -                                                27.
C        MRGLNK - A VECTOR OF SIZE NEQNS.  AT THE KTH STEP,               28.
C               MRGLNK(K), MRGLNK(MRGLNK(K)) , .........                  29.
C               IS A LIST CONTAINING ALL THOSE COLUMNS L(*,J)             30.
C               WITH J LESS THAN K, SUCH THAT ITS FIRST OFF-              31.
C               DIAGONAL NONZERO IS L(K,J).  THUS, THE                    32.
C               NONZERO STRUCTURE OF COLUMN L(*,K) CAN BE FOUND           33.
C               BY MERGING THAT OF SUCH COLUMNS L(*,J) WITH               34.
C               THE STRUCTURE OF A(*,K).                                  35.
C        RCHLNK - A VECTOR OF SIZE NEQNS.  IT IS USED TO ACCUMULATE       36.
C               THE STRUCTURE OF EACH COLUMN L(*,K).  AT THE              37.
C               END OF THE KTH STEP,                                      38.
C                   RCHLNK(K), RCHLNK(RCHLNK(K)), ........                39.
C               IS THE LIST OF POSITIONS OF NONZEROS IN COLUMN K          40.
C               OF THE FACTOR L.                                          41.
C        MARKER  - AN INTEGER VECTOR OF LENGTH NEQNS. IT IS USED          42.
C               TO TEST IF MASS SYMBOLIC ELIMINATION CAN BE               43.
C               PERFORMED.  THAT IS, IT IS USED TO CHECK WHETHER          44.
C               THE STRUCTURE OF THE CURRENT COLUMN K BEING               45.
C               PROCESSED IS COMPLETELY DETERMINED BY THE SINGLE          46.
C               COLUMN MRGLNK(K).                                         47.
C                                                                         48.
C****************************************************************         49.
C                                                                         50.
      SUBROUTINE  SMBFCT ( NEQNS, XADJ, ADJNCY, PERM, INVP,               51.
     1                     XLNZ, MAXLNZ, XNZSUB, NZSUB, MAXSUB,           52.
     1                     RCHLNK, MRGLNK, MARKER, FLAG )                 53.
C                                                                         54.
C****************************************************************         55.
C                                                                         56.
         INTEGER ADJNCY(1), INVP(1), MRGLNK(1), NZSUB(1),                 57.
     1           PERM(1), RCHLNK(1), MARKER(1)                            58.
         INTEGER XADJ(1), XLNZ(1), XNZSUB(1),                             59.
     1           FLAG, I, INZ, J, JSTOP, JSTRT, K, KNZ,                   60.
     1           KXSUB, MRGK, LMAX, M, MAXLNZ, MAXSUB,                    61.
     1           NABOR, NEQNS, NODE, NP1, NZBEG, NZEND,                   62.
     1           RCHM, MRKFLG                                             63.
C                                                                         64.
C****************************************************************         65.
C                                                                         66.
C       ------------------                                                67.
C       INITIALIZATION ...                                                68.
C       ------------------                                                69.
        NZBEG = 1                                                         70.
        NZEND = 0                                                         71.
        XLNZ(1) = 1                                                       72.
        DO 100 K = 1, NEQNS                                               73.
           MRGLNK(K) = 0                                                  74.
           MARKER(K) = 0                                                  75.
  100   CONTINUE                                                          76.
C       --------------------------------------------------                77.
C       FOR EACH COLUMN ......... .  KNZ COUNTS THE NUMBER                78.
C       OF NONZEROS IN COLUMN K ACCUMULATED IN RCHLNK.                    79.
C       --------------------------------------------------                80.
        NP1 = NEQNS + 1                                                   81.
        DO 1500  K = 1, NEQNS                                             82.
           KNZ = 0                                                        83.
           MRGK = MRGLNK(K)                                               84.
           MRKFLG = 0                                                     85.
           MARKER(K) = K                                                  86.
           IF (MRGK .NE. 0 ) MARKER(K) = MARKER(MRGK)                     87.
           XNZSUB(K) = NZEND                                              88.
           NODE = PERM(K)                                                 89.
           JSTRT = XADJ(NODE)                                             90.
           JSTOP = XADJ(NODE+1) - 1                                       91.
           IF (JSTRT.GT.JSTOP)  GO TO 1500                                92.
C          -------------------------------------------                    93.
C          USE RCHLNK TO LINK THROUGH THE STRUCTURE OF                    94.
C          A(*,K) BELOW DIAGONAL                                          95.
C          -------------------------------------------                    96.
           RCHLNK(K) = NP1                                                97.
           DO 300 J = JSTRT, JSTOP                                        98.
              NABOR = ADJNCY(J)                                           99.
              NABOR = INVP(NABOR)                                        100.
              IF ( NABOR .LE. K )  GO TO 300                             101.
                 RCHM = K                                                102.
  200            M = RCHM                                                103.
                 RCHM = RCHLNK(M)                                        104.
                 IF ( RCHM .LE. NABOR )  GO TO 200                       105.
                    KNZ = KNZ+1                                          106.
                    RCHLNK(M) = NABOR                                    107.
                    RCHLNK(NABOR) = RCHM                                 108.
                    IF ( MARKER(NABOR) .NE. MARKER(K) )  MRKFLG = 1      109.
  300      CONTINUE                                                      110.
C          --------------------------------------                        111.
C          TEST FOR MASS SYMBOLIC ELIMINATION ...                        112.
C          --------------------------------------                        113.
           LMAX = 0                                                      114.
           IF ( MRKFLG .NE. 0 .OR. MRGK .EQ. 0 ) GO TO 350               115.
           IF ( MRGLNK(MRGK) .NE. 0 ) GO TO 350                          116.
           XNZSUB(K) = XNZSUB(MRGK) + 1                                  117.
           KNZ = XLNZ(MRGK+1) - (XLNZ(MRGK) + 1)                         118.
           GO TO 1400                                                    119.
C          -----------------------------------------------               120.
C          LINK THROUGH EACH COLUMN I THAT AFFECTS L(*,K).               121.
C          -----------------------------------------------               122.
  350      I = K                                                         123.
  400      I = MRGLNK(I)                                                 124.
           IF (I.EQ.0)  GO TO 800                                        125.
              INZ = XLNZ(I+1) - (XLNZ(I)+1)                              126.
              JSTRT = XNZSUB(I) +  1                                     127.
              JSTOP = XNZSUB(I) + INZ                                    128.
              IF (INZ.LE.LMAX)  GO TO 500                                129.
                 LMAX = INZ                                              130.
                 XNZSUB(K) = JSTRT                                       131.
C             -----------------------------------------------            132.
C             MERGE STRUCTURE OF L(*,I) IN NZSUB INTO RCHLNK.            133.
C             -----------------------------------------------            134.
  500         RCHM = K                                                   135.
              DO 700 J = JSTRT, JSTOP                                    136.
                 NABOR = NZSUB(J)                                        137.
  600            M = RCHM                                                138.
                 RCHM = RCHLNK(M)                                        139.
                 IF (RCHM.LT.NABOR)  GO TO 600                           140.
                 IF (RCHM.EQ.NABOR)  GO TO 700                           141.
                    KNZ = KNZ+1                                          142.
                    RCHLNK(M) = NABOR                                    143.
                    RCHLNK(NABOR) = RCHM                                 144.
                    RCHM = NABOR                                         145.
  700         CONTINUE                                                   146.
              GO TO 400                                                  147.
C          ------------------------------------------------------        148.
C          CHECK IF SUBSCRIPTS DUPLICATE THOSE OF ANOTHER COLUMN.        149.
C          ------------------------------------------------------        150.
  800      IF (KNZ.EQ.LMAX)  GO TO 1400                                  151.
C             -----------------------------------------------            152.
C             OR IF TAIL OF K-1ST COLUMN MATCHES HEAD OF KTH.            153.
C             -----------------------------------------------            154.
              IF (NZBEG.GT.NZEND)  GO TO 1200                            155.
                 I = RCHLNK(K)                                           156.
                 DO 900 JSTRT=NZBEG,NZEND                                157.
                    IF (NZSUB(JSTRT)-I)  900, 1000, 1200                 158.
  900            CONTINUE                                                159.
                 GO TO 1200                                              160.
 1000            XNZSUB(K) = JSTRT                                       161.
                 DO 1100 J=JSTRT,NZEND                                   162.
                    IF (NZSUB(J).NE.I)  GO TO 1200                       163.
                    I = RCHLNK(I)                                        164.
                    IF (I.GT.NEQNS)  GO TO 1400                          165.
 1100            CONTINUE                                                166.
                 NZEND = JSTRT - 1                                       167.
C             ----------------------------------------                   168.
C             COPY THE STRUCTURE OF L(*,K) FROM RCHLNK                   169.
C             TO THE DATA STRUCTURE (XNZSUB, NZSUB).                     170.
C             ----------------------------------------                   171.
 1200         NZBEG = NZEND +  1                                         172.
              NZEND = NZEND + KNZ                                        173.
              IF (NZEND.GT.MAXSUB)  GO TO 1600                           174.
              I = K                                                      175.
              DO 1300 J=NZBEG,NZEND                                      176.
                 I = RCHLNK(I)                                           177.
                 NZSUB(J) = I                                            178.
                 MARKER(I) = K                                           179.
 1300         CONTINUE                                                   180.
              XNZSUB(K) = NZBEG                                          181.
              MARKER(K) = K                                              182.
C          --------------------------------------------------------      183.
C          UPDATE THE VECTOR MRGLNK.  NOTE COLUMN L(*,K) JUST FOUND      184.
C          IS REQUIRED TO DETERMINE COLUMN L(*,J), WHERE                 185.
C          L(J,K) IS THE FIRST NONZERO IN L(*,K) BELOW DIAGONAL.         186.
C          --------------------------------------------------------      187.
 1400      IF (KNZ.LE.1)  GO TO 1500                                     188.
              KXSUB = XNZSUB(K)                                          189.
              I = NZSUB(KXSUB)                                           190.
              MRGLNK(K) = MRGLNK(I)                                      191.
              MRGLNK(I) = K                                              192.
 1500      XLNZ(K+1) = XLNZ(K) + KNZ                                     193.
        MAXLNZ = XLNZ(NEQNS) - 1                                         194.
        MAXSUB = XNZSUB(NEQNS)                                           195.
        XNZSUB(NEQNS+1) = XNZSUB(NEQNS)                                  196.
        FLAG = 0                                                         197.
        RETURN                                                           198.
C       ----------------------------------------------------             199.
C       ERROR - INSUFFICIENT STORAGE FOR NONZERO SUBSCRIPTS.             200.
C       ----------------------------------------------------             201.
 1600   FLAG = 1                                                         202.
        RETURN                                                           203.
        END                                                              204.
C----- SUBROUTINE GSFCT
C***************************************************************           1.
C***************************************************************           2.
C******     GSFCT ..... GENERAL SPARSE SYMMETRIC FACT     ******           3.
C***************************************************************           4.
C***************************************************************           5.
C                                                                          6.
C     PURPOSE - THIS SUBROUTINE PERFORMS THE SYMMETRIC                     7.
C        FACTORIZATION FOR A GENERAL SPARSE SYSTEM, STORED IN              8.
C        THE COMPRESSED SUBSCRIPT DATA FORMAT.                             9.
C                                                                         10.
C     INPUT PARAMETERS -                                                  11.
C        NEQNS - NUMBER OF EQUATIONS.                                     12.
C        XLNZ - INDEX VECTOR FOR LNZ.  XLNZ(I) POINTS TO THE              13.
C               START OF NONZEROS IN COLUMN I OF FACTOR L.                14.
C        (XNZSUB, NZSUB) - THE COMPRESSED SUBSCRIPT DATA                  15.
C               STRUCTURE FOR FACTOR L.                                   16.
C                                                                         17.
C     UPDATED PARAMETERS -                                                18.
C        LNZ - ON INPUT, CONTAINS NONZEROS OF A, AND ON                   19.
C               RETURN, THE NONZEROS OF L.                                20.
C        DIAG - THE DIAGONAL OF L OVERWRITES THAT OF A.                   21.
C        IFLAG - THE ERROR FLAG.  IT IS SET TO 1 IF A ZERO OR             22.
C               NEGATIVE SQUARE ROOT OCCURS DURING THE                    23.
C               FACTORIZATION.                                            24.
C        OPS   - A DOUBLE PRECISION COMMON PARAMETER THAT IS              25.
C                INCREMENTED BY THE NUMBER OF OPERATIONS                  26.
C                PERFORMED BY THE SUBROUTINE.                             27.
C                                                                         28.
C     WORKING PARAMETERS -                                                29.
C        LINK - AT STEP J, THE LIST IN                                    30.
C                  LINK(J), LINK(LINK(J)), ...........                    31.
C               CONSISTS OF THOSE COLUMNS THAT WILL MODIFY                32.
C               THE COLUMN L(*,J).                                        33.
C        FIRST - TEMPORARY VECTOR TO POINT TO THE FIRST                   34.
C               NONZERO IN EACH COLUMN THAT WILL BE USED                  35.
C               NEXT FOR MODIFICATION.                                    36.
C        TEMP - A TEMPORARY VECTOR TO ACCUMULATE MODIFICATIONS.           37.
C                                                                         38.
C***************************************************************          39.
C                                                                         40.
      SUBROUTINE  GSFCT ( NEQNS, XLNZ, LNZ, XNZSUB, NZSUB, DIAG,          41.
     1                    LINK, FIRST, TEMP, IFLAG )                      42.
C                                                                         43.
C***************************************************************          44.
C                                                                         45.
         DOUBLE PRECISION COUNT, OPS                                      46.
         COMMON  /SPKOPS/ OPS                                             47.
         REAL DIAG(1), LNZ(1), TEMP(1), DIAGJ, LJK                        48.
         INTEGER LINK(1), NZSUB(1)                                        49.
         INTEGER FIRST(1), XLNZ(1), XNZSUB(1),                            50.
     1           I, IFLAG, II, ISTOP, ISTRT, ISUB, J,                     51.
     1           K, KFIRST, NEQNS, NEWK                                   52.
C                                                                         53.
C***************************************************************          54.
C                                                                         55.
C        ------------------------------                                   56.
C        INITIALIZE WORKING VECTORS ...                                   57.
C        ------------------------------                                   58.
         DO 100 I = 1, NEQNS                                              59.
            LINK(I) = 0                                                   60.
            TEMP(I) = 0.0E0                                               61.
  100    CONTINUE                                                         62.
C        --------------------------------------------                     63.
C        COMPUTE COLUMN L(*,J) FOR J = 1,...., NEQNS.                     64.
C        --------------------------------------------                     65.
         DO 600 J = 1, NEQNS                                              66.
C           -------------------------------------------                   67.
C           FOR EACH COLUMN L(*,K) THAT AFFECTS L(*,J).                   68.
C           -------------------------------------------                   69.
            DIAGJ = 0.0E0                                                 70.
            NEWK = LINK(J)                                                71.
  200       K    = NEWK                                                   72.
            IF ( K .EQ. 0 )  GO TO 400                                    73.
               NEWK = LINK(K)                                             74.
C              ---------------------------------------                    75.
C              OUTER PRODUCT MODIFICATION OF L(*,J) BY                    76.
C              L(*,K) STARTING AT FIRST(K) OF L(*,K).                     77.
C              ---------------------------------------                    78.
               KFIRST = FIRST(K)                                          79.
               LJK    = LNZ(KFIRST)                                       80.
               DIAGJ = DIAGJ + LJK*LJK                                    81.
               OPS  = OPS + 1.0D0                                         82.
               ISTRT = KFIRST + 1                                         83.
               ISTOP = XLNZ(K+1) - 1                                      84.
               IF ( ISTOP .LT. ISTRT )  GO TO 200                         85.
C                 ------------------------------------------              86.
C                 BEFORE MODIFICATION, UPDATE VECTORS FIRST,              87.
C                 AND LINK FOR FUTURE MODIFICATION STEPS.                 88.
C                 ------------------------------------------              89.
                  FIRST(K) = ISTRT                                        90.
                  I = XNZSUB(K) + (KFIRST-XLNZ(K)) + 1                    91.
                  ISUB = NZSUB(I)                                         92.
                  LINK(K) = LINK(ISUB)                                    93.
                  LINK(ISUB) = K                                          94.
C                 ---------------------------------------                 95.
C                 THE ACTUAL MOD IS SAVED IN VECTOR TEMP.                 96.
C                 ---------------------------------------                 97.
                  DO 300 II = ISTRT, ISTOP                                98.
                     ISUB = NZSUB(I)                                      99.
                     TEMP(ISUB) = TEMP(ISUB) + LNZ(II)*LJK               100.
                     I = I + 1                                           101.
  300             CONTINUE                                               102.
                  COUNT = ISTOP - ISTRT + 1                              103.
                  OPS  = OPS + COUNT                                     104.
               GO TO 200                                                 105.
C           ----------------------------------------------               106.
C           APPLY THE MODIFICATIONS ACCUMULATED IN TEMP TO               107.
C           COLUMN L(*,J).                                               108.
C           ----------------------------------------------               109.
  400       DIAGJ = DIAG(J) - DIAGJ                                      110.
            IF ( DIAGJ .LE. 0.0E0 )  GO TO 700                           111.
            DIAGJ = SQRT(DIAGJ)                                          112.
            DIAG(J) = DIAGJ                                              113.
            ISTRT = XLNZ(J)                                              114.
            ISTOP = XLNZ(J+1) - 1                                        115.
            IF ( ISTOP .LT. ISTRT )  GO TO 600                           116.
               FIRST(J) = ISTRT                                          117.
               I = XNZSUB(J)                                             118.
               ISUB = NZSUB(I)                                           119.
               LINK(J) = LINK(ISUB)                                      120.
               LINK(ISUB) = J                                            121.
               DO 500 II = ISTRT, ISTOP                                  122.
                  ISUB = NZSUB(I)                                        123.
                  LNZ(II) = ( LNZ(II)-TEMP(ISUB) ) / DIAGJ               124.
                  TEMP(ISUB) = 0.0E0                                     125.
                  I = I + 1                                              126.
  500          CONTINUE                                                  127.
               COUNT = ISTOP - ISTRT + 1                                 128.
               OPS  = OPS + COUNT                                        129.
  600    CONTINUE                                                        130.
         RETURN                                                          131.
C        ------------------------------------------------------          132.
C        ERROR - ZERO OR NEGATIVE SQUARE ROOT IN FACTORIZATION.          133.
C        ------------------------------------------------------          134.
  700    IFLAG = 1                                                       135.
         RETURN                                                          136.
      END                                                                137.

C***************************************************************
C***************************************************************
C****     GENMMD ..... MULTIPLE MINIMUM EXTERNAL DEGREE     ****
C***************************************************************
C***************************************************************
C
C     AUTHOR - JOSEPH W.H. LIU
C              DEPT OF COMPUTER SCIENCE, YORK UNIVERSITY.
C
C     PURPOSE - THIS ROUTINE IMPLEMENTS THE MINIMUM DEGREE
C        ALGORITHM.  IT MAKES USE OF THE IMPLICIT REPRESENTATION
C        OF ELIMINATION GRAPHS BY QUOTIENT GRAPHS, AND THE
C        NOTION OF INDISTINGUISHABLE NODES.  IT ALSO IMPLEMENTS
C        THE MODIFICATIONS BY MULTIPLE ELIMINATION AND MINIMUM
C        EXTERNAL DEGREE.
C        ---------------------------------------------
C        CAUTION - THE ADJACENCY VECTOR ADJNCY WILL BE
C        DESTROYED.
C        ---------------------------------------------
C
C     INPUT PARAMETERS -
C        NEQNS  - NUMBER OF EQUATIONS.
C        (XADJ,ADJNCY) - THE ADJACENCY STRUCTURE.
C        DELTA  - TOLERANCE VALUE FOR MULTIPLE ELIMINATION.
C        MAXINT - MAXIMUM MACHINE REPRESENTABLE (SHORT) INTEGER
C                 (ANY SMALLER ESTIMATE WILL DO) FOR MARKING
C                 NODES.
C
C     OUTPUT PARAMETERS -
C        PERM   - THE MINIMUM DEGREE ORDERING.
C        INVP   - THE INVERSE OF PERM.
C        NOFSUB - AN UPPER BOUND ON THE NUMBER OF NONZERO
C                 SUBSCRIPTS FOR THE COMPRESSED STORAGE SCHEME.
C
C     WORKING PARAMETERS -
C        DHEAD  - VECTOR FOR HEAD OF DEGREE LISTS.
C        INVP   - USED TEMPORARILY FOR DEGREE FORWARD LINK.
C        PERM   - USED TEMPORARILY FOR DEGREE BACKWARD LINK.
C        QSIZE  - VECTOR FOR SIZE OF SUPERNODES.
C        LLIST  - VECTOR FOR TEMPORARY LINKED LISTS.
C        MARKER - A TEMPORARY MARKER VECTOR.
C
C     PROGRAM SUBROUTINES -
C        MMDELM, MMDINT, MMDNUM, MMDUPD.
C
C***************************************************************
C
      SUBROUTINE  GENMMD ( NEQNS, XADJ, ADJNCY, INVP, PERM,
     1                     DELTA, DHEAD, QSIZE, LLIST, MARKER,
     1                     MAXINT, NOFSUB )
C
C***************************************************************
C
C         INTEGER*2  ADJNCY(1), DHEAD(1) , INVP(1)  , LLIST(1) ,
         INTEGER*4  ADJNCY(1), DHEAD(1) , INVP(1)  , LLIST(1) ,
     1              MARKER(1), PERM(1)  , QSIZE(1)
         INTEGER*4  XADJ(1)
         INTEGER*4  DELTA , EHEAD , I     , MAXINT, MDEG  ,
     1              MDLMT , MDNODE, NEQNS , NEXTMD, NOFSUB,
     1              NUM, TAG
C
C***************************************************************
C
         IF  ( NEQNS .LE. 0 )  RETURN
C
C        ------------------------------------------------
C        INITIALIZATION FOR THE MINIMUM DEGREE ALGORITHM.
C        ------------------------------------------------
         NOFSUB = 0
         CALL  MMDINT ( NEQNS, XADJ, ADJNCY, DHEAD, INVP, PERM,
     1                  QSIZE, LLIST, MARKER )
C
C        ----------------------------------------------
C        NUM COUNTS THE NUMBER OF ORDERED NODES PLUS 1.
C        ----------------------------------------------
         NUM = 1
C
C        -----------------------------
C        ELIMINATE ALL ISOLATED NODES.
C        -----------------------------
         NEXTMD = DHEAD(1)
  100    CONTINUE
             IF  ( NEXTMD .LE. 0 )  GO TO 200
                 MDNODE = NEXTMD
                 NEXTMD = INVP(MDNODE)
                 MARKER(MDNODE) = MAXINT
                 INVP(MDNODE) = - NUM
                 NUM = NUM + 1
                 GO TO 100
C
  200    CONTINUE
C        ----------------------------------------
C        SEARCH FOR NODE OF THE MINIMUM DEGREE.
C        MDEG IS THE CURRENT MINIMUM DEGREE;
C        TAG IS USED TO FACILITATE MARKING NODES.
C        ----------------------------------------
         IF  ( NUM .GT. NEQNS )  GO TO 1000
         TAG = 1
         DHEAD(1) = 0
         MDEG = 2
  300    CONTINUE
             IF  ( DHEAD(MDEG) .GT. 0 )  GO TO 400
                 MDEG = MDEG + 1
                 GO TO 300
  400        CONTINUE
C            -------------------------------------------------
C            USE VALUE OF DELTA TO SET UP MDLMT, WHICH GOVERNS
C            WHEN A DEGREE UPDATE IS TO BE PERFORMED.
C            -------------------------------------------------
             MDLMT = MDEG + DELTA
             EHEAD = 0
C
  500        CONTINUE
                 MDNODE = DHEAD(MDEG)
                 IF  ( MDNODE .GT. 0 )  GO TO 600
                     MDEG = MDEG + 1
                     IF  ( MDEG .GT. MDLMT )  GO TO 900
                         GO TO 500
  600            CONTINUE
C                ----------------------------------------
C                REMOVE MDNODE FROM THE DEGREE STRUCTURE.
C                ----------------------------------------
                 NEXTMD = INVP(MDNODE)
                 DHEAD(MDEG) = NEXTMD
                 IF  ( NEXTMD .GT. 0 )  PERM(NEXTMD) = - MDEG
                 INVP(MDNODE) = - NUM
                 NOFSUB = NOFSUB + MDEG + QSIZE(MDNODE) - 2
                 IF  ( NUM+QSIZE(MDNODE) .GT. NEQNS )  GO TO 1000
C                ----------------------------------------------
C                ELIMINATE MDNODE AND PERFORM QUOTIENT GRAPH
C                TRANSFORMATION.  RESET TAG VALUE IF NECESSARY.
C                ----------------------------------------------
                 TAG = TAG + 1
                 IF  ( TAG .LT. MAXINT )  GO TO 800
                     TAG = 1
                     DO  700  I = 1, NEQNS
                         IF  ( MARKER(I) .LT. MAXINT )  MARKER(I) = 0
  700                CONTINUE
  800            CONTINUE
                 CALL  MMDELM ( MDNODE, XADJ, ADJNCY, DHEAD, INVP,
     1                          PERM, QSIZE, LLIST, MARKER, MAXINT,
     1                          TAG )
                 NUM = NUM + QSIZE(MDNODE)
                 LLIST(MDNODE) = EHEAD
                 EHEAD = MDNODE
                 IF  ( DELTA .GE. 0 )  GO TO 500
  900        CONTINUE
C            -------------------------------------------
C            UPDATE DEGREES OF THE NODES INVOLVED IN THE
C            MINIMUM DEGREE NODES ELIMINATION.
C            -------------------------------------------
             IF  ( NUM .GT. NEQNS )  GO TO 1000
             CALL  MMDUPD ( EHEAD, NEQNS, XADJ, ADJNCY, DELTA, MDEG,
     1                      DHEAD, INVP, PERM, QSIZE, LLIST, MARKER,
     1                      MAXINT, TAG )
             GO TO 300
C
 1000    CONTINUE
         CALL  MMDNUM ( NEQNS, PERM, INVP, QSIZE )
         RETURN
C
      END
C***************************************************************
C***************************************************************
C***     MMDINT ..... MULT MINIMUM DEGREE INITIALIZATION     ***
C***************************************************************
C***************************************************************
C
C     AUTHOR - JOSEPH W.H. LIU
C              DEPT OF COMPUTER SCIENCE, YORK UNIVERSITY.
C
C     PURPOSE - THIS ROUTINE PERFORMS INITIALIZATION FOR THE
C        MULTIPLE ELIMINATION VERSION OF THE MINIMUM DEGREE
C        ALGORITHM.
C
C     INPUT PARAMETERS -
C        NEQNS  - NUMBER OF EQUATIONS.
C        (XADJ,ADJNCY) - ADJACENCY STRUCTURE.
C
C     OUTPUT PARAMETERS -
C        (DHEAD,DFORW,DBAKW) - DEGREE DOUBLY LINKED STRUCTURE.
C        QSIZE  - SIZE OF SUPERNODE (INITIALIZED TO ONE).
C        LLIST  - LINKED LIST.
C        MARKER - MARKER VECTOR.
C
C***************************************************************
C
      SUBROUTINE  MMDINT ( NEQNS, XADJ, ADJNCY, DHEAD, DFORW,
     1                     DBAKW, QSIZE, LLIST, MARKER )
C
C***************************************************************
C
C         INTEGER*2  ADJNCY(1), DBAKW(1) , DFORW(1) , DHEAD(1) ,
         INTEGER*4  ADJNCY(1), DBAKW(1) , DFORW(1) , DHEAD(1) ,
     1              LLIST(1) , MARKER(1), QSIZE(1)
         INTEGER*4  XADJ(1)
         INTEGER*4  FNODE , NDEG  , NEQNS , NODE
C
C***************************************************************
C
         DO  100  NODE = 1, NEQNS
             DHEAD(NODE) = 0
             QSIZE(NODE) = 1
             MARKER(NODE) = 0
             LLIST(NODE) = 0
  100    CONTINUE
C        ------------------------------------------
C        INITIALIZE THE DEGREE DOUBLY LINKED LISTS.
C        ------------------------------------------
         DO  200  NODE = 1, NEQNS
             NDEG = XADJ(NODE+1) - XADJ(NODE) + 1
             FNODE = DHEAD(NDEG)
             DFORW(NODE) = FNODE
             DHEAD(NDEG) = NODE
             IF  ( FNODE .GT. 0 )  DBAKW(FNODE) = NODE
             DBAKW(NODE) = - NDEG
  200    CONTINUE
         RETURN
C
      END
C***************************************************************
C***************************************************************
C**     MMDELM ..... MULTIPLE MINIMUM DEGREE ELIMINATION     ***
C***************************************************************
C***************************************************************
C
C     AUTHOR - JOSEPH W.H. LIU
C              DEPT OF COMPUTER SCIENCE, YORK UNIVERSITY.
C
C     PURPOSE - THIS ROUTINE ELIMINATES THE NODE MDNODE OF
C        MINIMUM DEGREE FROM THE ADJACENCY STRUCTURE, WHICH
C        IS STORED IN THE QUOTIENT GRAPH FORMAT.  IT ALSO
C        TRANSFORMS THE QUOTIENT GRAPH REPRESENTATION OF THE
C        ELIMINATION GRAPH.
C
C     INPUT PARAMETERS -
C        MDNODE - NODE OF MINIMUM DEGREE.
C        MAXINT - ESTIMATE OF MAXIMUM REPRESENTABLE (SHORT)
C                 INTEGER.
C        TAG    - TAG VALUE.
C
C     UPDATED PARAMETERS -
C        (XADJ,ADJNCY) - UPDATED ADJACENCY STRUCTURE.
C        (DHEAD,DFORW,DBAKW) - DEGREE DOUBLY LINKED STRUCTURE.
C        QSIZE  - SIZE OF SUPERNODE.
C        MARKER - MARKER VECTOR.
C        LLIST  - TEMPORARY LINKED LIST OF ELIMINATED NABORS.
C
C***************************************************************
C
      SUBROUTINE  MMDELM ( MDNODE, XADJ, ADJNCY, DHEAD, DFORW,
     1                     DBAKW, QSIZE, LLIST, MARKER, MAXINT,
     1                     TAG )
C
C***************************************************************
C
C         INTEGER*2  ADJNCY(1), DBAKW(1) , DFORW(1) , DHEAD(1) ,
         INTEGER*4  ADJNCY(1), DBAKW(1) , DFORW(1) , DHEAD(1) ,
     1              LLIST(1) , MARKER(1), QSIZE(1)
         INTEGER*4  XADJ(1)
         INTEGER*4  ELMNT , I     , ISTOP , ISTRT , J     ,
     1              JSTOP , JSTRT , LINK  , MAXINT, MDNODE,
     1              NABOR , NODE  , NPV   , NQNBRS, NXNODE,
     1              PVNODE, RLMT  , RLOC  , RNODE , TAG   ,
     1              XQNBR
C
C***************************************************************
C
C        -----------------------------------------------
C        FIND REACHABLE SET AND PLACE IN DATA STRUCTURE.
C        -----------------------------------------------
         MARKER(MDNODE) = TAG
         ISTRT = XADJ(MDNODE)
         ISTOP = XADJ(MDNODE+1) - 1
C        -------------------------------------------------------
C        ELMNT POINTS TO THE BEGINNING OF THE LIST OF ELIMINATED
C        NABORS OF MDNODE, AND RLOC GIVES THE STORAGE LOCATION
C        FOR THE NEXT REACHABLE NODE.
C        -------------------------------------------------------
         ELMNT = 0
         RLOC = ISTRT
         RLMT = ISTOP
         DO  200  I = ISTRT, ISTOP
             NABOR = ADJNCY(I)
             IF  ( NABOR .EQ. 0 )  GO TO 300
                 IF  ( MARKER(NABOR) .GE. TAG )  GO TO 200
                     MARKER(NABOR) = TAG
                     IF  ( DFORW(NABOR) .LT. 0 )  GO TO 100
                         ADJNCY(RLOC) = NABOR
                         RLOC = RLOC + 1
                         GO TO 200
  100                CONTINUE
                     LLIST(NABOR) = ELMNT
                     ELMNT = NABOR
  200    CONTINUE
  300    CONTINUE
C            -----------------------------------------------------
C            MERGE WITH REACHABLE NODES FROM GENERALIZED ELEMENTS.
C            -----------------------------------------------------
             IF  ( ELMNT .LE. 0 )  GO TO 1000
                 ADJNCY(RLMT) = - ELMNT
                 LINK = ELMNT
  400            CONTINUE
                     JSTRT = XADJ(LINK)
                     JSTOP = XADJ(LINK+1) - 1
                     DO  800  J = JSTRT, JSTOP
                         NODE = ADJNCY(J)
                         LINK = - NODE
                         IF  ( NODE )  400, 900, 500
  500                    CONTINUE
                         IF  ( MARKER(NODE) .GE. TAG  .OR.
     1                         DFORW(NODE) .LT. 0 )  GO TO 800
                             MARKER(NODE) = TAG
C                            ---------------------------------
C                            USE STORAGE FROM ELIMINATED NODES
C                            IF NECESSARY.
C                            ---------------------------------
  600                        CONTINUE
                                 IF  ( RLOC .LT. RLMT )  GO TO 700
                                     LINK = - ADJNCY(RLMT)
                                     RLOC = XADJ(LINK)
                                     RLMT = XADJ(LINK+1) - 1
                                     GO TO 600
  700                        CONTINUE
                             ADJNCY(RLOC) = NODE
                             RLOC = RLOC + 1
  800                CONTINUE
  900            CONTINUE
                 ELMNT = LLIST(ELMNT)
                 GO TO 300
 1000    CONTINUE
         IF  ( RLOC .LE. RLMT )  ADJNCY(RLOC) = 0
C        --------------------------------------------------------
C        FOR EACH NODE IN THE REACHABLE SET, DO THE FOLLOWING ...
C        --------------------------------------------------------
         LINK = MDNODE
 1100    CONTINUE
             ISTRT = XADJ(LINK)
             ISTOP = XADJ(LINK+1) - 1
             DO  1700  I = ISTRT, ISTOP
                 RNODE = ADJNCY(I)
                 LINK = - RNODE
                 IF  ( RNODE )  1100, 1800, 1200
 1200            CONTINUE
C                --------------------------------------------
C                IF RNODE IS IN THE DEGREE LIST STRUCTURE ...
C                --------------------------------------------
                 PVNODE = DBAKW(RNODE)
                 IF  ( PVNODE .EQ. 0  .OR.
     1                 PVNODE .EQ. (-MAXINT) )  GO TO 1300
C                    -------------------------------------
C                    THEN REMOVE RNODE FROM THE STRUCTURE.
C                    -------------------------------------
                     NXNODE = DFORW(RNODE)
                     IF  ( NXNODE .GT. 0 )  DBAKW(NXNODE) = PVNODE
                     IF  ( PVNODE .GT. 0 )  DFORW(PVNODE) = NXNODE
                     NPV = - PVNODE
                     IF  ( PVNODE .LT. 0 )  DHEAD(NPV) = NXNODE
 1300            CONTINUE
C                ----------------------------------------
C                PURGE INACTIVE QUOTIENT NABORS OF RNODE.
C                ----------------------------------------
                 JSTRT = XADJ(RNODE)
                 JSTOP = XADJ(RNODE+1) - 1
                 XQNBR = JSTRT
                 DO  1400  J = JSTRT, JSTOP
                     NABOR = ADJNCY(J)
                     IF  ( NABOR .EQ. 0 )  GO TO 1500
                         IF  ( MARKER(NABOR) .GE. TAG )  GO TO 1400
                             ADJNCY(XQNBR) = NABOR
                             XQNBR = XQNBR + 1
 1400            CONTINUE
 1500            CONTINUE
C                ----------------------------------------
C                IF NO ACTIVE NABOR AFTER THE PURGING ...
C                ----------------------------------------
                 NQNBRS = XQNBR - JSTRT
                 IF  ( NQNBRS .GT. 0 )  GO TO 1600
C                    -----------------------------
C                    THEN MERGE RNODE WITH MDNODE.
C                    -----------------------------
                     QSIZE(MDNODE) = QSIZE(MDNODE) + QSIZE(RNODE)
                     QSIZE(RNODE) = 0
                     MARKER(RNODE) = MAXINT
                     DFORW(RNODE) = - MDNODE
                     DBAKW(RNODE) = - MAXINT
                     GO TO 1700
 1600            CONTINUE
C                --------------------------------------
C                ELSE FLAG RNODE FOR DEGREE UPDATE, AND
C                ADD MDNODE AS A NABOR OF RNODE.
C                --------------------------------------
                 DFORW(RNODE) = NQNBRS + 1
                 DBAKW(RNODE) = 0
                 ADJNCY(XQNBR) = MDNODE
                 XQNBR = XQNBR + 1
                 IF  ( XQNBR .LE. JSTOP )  ADJNCY(XQNBR) = 0
C
 1700        CONTINUE
 1800    CONTINUE
         RETURN
C
      END
C***************************************************************
C***************************************************************
C*****     MMDUPD ..... MULTIPLE MINIMUM DEGREE UPDATE     *****
C***************************************************************
C***************************************************************
C
C     AUTHOR - JOSEPH W.H. LIU
C              DEPT OF COMPUTER SCIENCE, YORK UNIVERSITY.
C
C     PURPOSE - THIS ROUTINE UPDATES THE DEGREES OF NODES
C        AFTER A MULTIPLE ELIMINATION STEP.
C
C     INPUT PARAMETERS -
C        EHEAD  - THE BEGINNING OF THE LIST OF ELIMINATED
C                 NODES (I.E., NEWLY FORMED ELEMENTS).
C        NEQNS  - NUMBER OF EQUATIONS.
C        (XADJ,ADJNCY) - ADJACENCY STRUCTURE.
C        DELTA  - TOLERANCE VALUE FOR MULTIPLE ELIMINATION.
C        MAXINT - MAXIMUM MACHINE REPRESENTABLE (SHORT)
C                 INTEGER.
C
C     UPDATED PARAMETERS -
C        MDEG   - NEW MINIMUM DEGREE AFTER DEGREE UPDATE.
C        (DHEAD,DFORW,DBAKW) - DEGREE DOUBLY LINKED STRUCTURE.
C        QSIZE  - SIZE OF SUPERNODE.
C        LLIST  - WORKING LINKED LIST.
C        MARKER - MARKER VECTOR FOR DEGREE UPDATE.
C        TAG    - TAG VALUE.
C
C***************************************************************
C
      SUBROUTINE  MMDUPD ( EHEAD, NEQNS, XADJ, ADJNCY, DELTA,
     1                     MDEG, DHEAD, DFORW, DBAKW, QSIZE,
     1                     LLIST, MARKER, MAXINT, TAG )
C
C***************************************************************
C
C         INTEGER*2  ADJNCY(1), DBAKW(1) , DFORW(1) , DHEAD(1) ,
         INTEGER*4  ADJNCY(1), DBAKW(1) , DFORW(1) , DHEAD(1) ,
     1              LLIST(1) , MARKER(1), QSIZE(1)
         INTEGER*4  XADJ(1)
         INTEGER*4  DEG   , DEG0  , DELTA , EHEAD , ELMNT ,
     1              ENODE , FNODE , I     , IQ2   , ISTOP ,
     1              ISTRT , J     , JSTOP , JSTRT , LINK  ,
     1              MAXINT, MDEG  , MDEG0 , MTAG  , NABOR ,
     1              NEQNS , NODE  , Q2HEAD, QXHEAD, TAG
C
C***************************************************************
C
         MDEG0 = MDEG + DELTA
         ELMNT = EHEAD
  100    CONTINUE
C            -------------------------------------------------------
C            FOR EACH OF THE NEWLY FORMED ELEMENT, DO THE FOLLOWING.
C            (RESET TAG VALUE IF NECESSARY.)
C            -------------------------------------------------------
             IF  ( ELMNT .LE. 0 )  RETURN
             MTAG = TAG + MDEG0
             IF  ( MTAG .LT. MAXINT )  GO TO 300
                 TAG = 1
                 DO  200  I = 1, NEQNS
                     IF  ( MARKER(I) .LT. MAXINT )  MARKER(I) = 0
  200            CONTINUE
                 MTAG = TAG + MDEG0
  300        CONTINUE
C            ---------------------------------------------
C            CREATE TWO LINKED LISTS FROM NODES ASSOCIATED
C            WITH ELMNT: ONE WITH TWO NABORS (Q2HEAD) IN
C            ADJACENCY STRUCTURE, AND THE OTHER WITH MORE
C            THAN TWO NABORS (QXHEAD).  ALSO COMPUTE DEG0,
C            NUMBER OF NODES IN THIS ELEMENT.
C            ---------------------------------------------
             Q2HEAD = 0
             QXHEAD = 0
             DEG0 = 0
             LINK = ELMNT
  400        CONTINUE
                 ISTRT = XADJ(LINK)
                 ISTOP = XADJ(LINK+1) - 1
                 DO  700  I = ISTRT, ISTOP
                     ENODE = ADJNCY(I)
                     LINK = - ENODE
                     IF  ( ENODE )  400, 800, 500
C
  500                CONTINUE
                     IF  ( QSIZE(ENODE) .EQ. 0 )  GO TO 700
                         DEG0 = DEG0 + QSIZE(ENODE)
                         MARKER(ENODE) = MTAG
C                        ----------------------------------
C                        IF ENODE REQUIRES A DEGREE UPDATE,
C                        THEN DO THE FOLLOWING.
C                        ----------------------------------
                         IF  ( DBAKW(ENODE) .NE. 0 )  GO TO 700
C                            ---------------------------------------
C                            PLACE EITHER IN QXHEAD OR Q2HEAD LISTS.
C                            ---------------------------------------
                             IF  ( DFORW(ENODE) .EQ. 2 )  GO TO 600
                                 LLIST(ENODE) = QXHEAD
                                 QXHEAD = ENODE
                                 GO TO 700
  600                        CONTINUE
                             LLIST(ENODE) = Q2HEAD
                             Q2HEAD = ENODE
  700            CONTINUE
  800        CONTINUE
C            --------------------------------------------
C            FOR EACH ENODE IN Q2 LIST, DO THE FOLLOWING.
C            --------------------------------------------
             ENODE = Q2HEAD
             IQ2 = 1
  900        CONTINUE
                 IF  ( ENODE .LE. 0 )  GO TO 1500
                 IF  ( DBAKW(ENODE) .NE. 0 )  GO TO 2200
                     TAG = TAG + 1
                     DEG = DEG0
C                    ------------------------------------------
C                    IDENTIFY THE OTHER ADJACENT ELEMENT NABOR.
C                    ------------------------------------------
                     ISTRT = XADJ(ENODE)
                     NABOR = ADJNCY(ISTRT)
                     IF  ( NABOR .EQ. ELMNT )  NABOR = ADJNCY(ISTRT+1)
C                    ------------------------------------------------
C                    IF NABOR IS UNELIMINATED, INCREASE DEGREE COUNT.
C                    ------------------------------------------------
                     LINK = NABOR
                     IF  ( DFORW(NABOR) .LT. 0 )  GO TO 1000
                         DEG = DEG + QSIZE(NABOR)
                         GO TO 2100
 1000                CONTINUE
C                        --------------------------------------------
C                        OTHERWISE, FOR EACH NODE IN THE 2ND ELEMENT,
C                        DO THE FOLLOWING.
C                        --------------------------------------------
                         ISTRT = XADJ(LINK)
                         ISTOP = XADJ(LINK+1) - 1
                         DO  1400  I = ISTRT, ISTOP
                             NODE = ADJNCY(I)
                             LINK = - NODE
                             IF  ( NODE .EQ. ENODE )  GO TO 1400
                             IF  ( NODE )  1000, 2100, 1100
C
 1100                        CONTINUE
                             IF  ( QSIZE(NODE) .EQ. 0 )  GO TO 1400
                             IF  ( MARKER(NODE) .GE. TAG )  GO TO 1200
C                                -------------------------------------
C                                CASE WHEN NODE IS NOT YET CONSIDERED.
C                                -------------------------------------
                                 MARKER(NODE) = TAG
                                 DEG = DEG + QSIZE(NODE)
                                 GO TO 1400
 1200                        CONTINUE
C                            ----------------------------------------
C                            CASE WHEN NODE IS INDISTINGUISHABLE FROM
C                            ENODE.  MERGE THEM INTO A NEW SUPERNODE.
C                            ----------------------------------------
                             IF  ( DBAKW(NODE) .NE. 0 )  GO TO 1400
                             IF  ( DFORW(NODE) .NE. 2 )  GO TO 1300
                                 QSIZE(ENODE) = QSIZE(ENODE) +
     1                                          QSIZE(NODE)
                                 QSIZE(NODE) = 0
                                 MARKER(NODE) = MAXINT
                                 DFORW(NODE) = - ENODE
                                 DBAKW(NODE) = - MAXINT
                                 GO TO 1400
 1300                        CONTINUE
C                            --------------------------------------
C                            CASE WHEN NODE IS OUTMATCHED BY ENODE.
C                            --------------------------------------
                             IF  ( DBAKW(NODE) .EQ.0 )
     1                             DBAKW(NODE) = - MAXINT
 1400                    CONTINUE
                         GO TO 2100
 1500            CONTINUE
C                ------------------------------------------------
C                FOR EACH ENODE IN THE QX LIST, DO THE FOLLOWING.
C                ------------------------------------------------
                 ENODE = QXHEAD
                 IQ2 = 0
 1600            CONTINUE
                     IF  ( ENODE .LE. 0 )  GO TO 2300
                     IF  ( DBAKW(ENODE) .NE. 0 )  GO TO 2200
                         TAG = TAG + 1
                         DEG = DEG0
C                        ---------------------------------
C                        FOR EACH UNMARKED NABOR OF ENODE,
C                        DO THE FOLLOWING.
C                        ---------------------------------
                         ISTRT = XADJ(ENODE)
                         ISTOP = XADJ(ENODE+1) - 1
                         DO  2000  I = ISTRT, ISTOP
                             NABOR = ADJNCY(I)
                             IF  ( NABOR .EQ. 0 )  GO TO 2100
                             IF  ( MARKER(NABOR) .GE. TAG )  GO TO 2000
                                 MARKER(NABOR) = TAG
                                 LINK = NABOR
C                                ------------------------------
C                                IF UNELIMINATED, INCLUDE IT IN
C                                DEG COUNT.
C                                ------------------------------
                                 IF  ( DFORW(NABOR) .LT. 0 )  GO TO 1700
                                     DEG = DEG + QSIZE(NABOR)
                                     GO TO 2000
 1700                            CONTINUE
C                                    -------------------------------
C                                    IF ELIMINATED, INCLUDE UNMARKED
C                                    NODES IN THIS ELEMENT INTO THE
C                                    DEGREE COUNT.
C                                    -------------------------------
                                     JSTRT = XADJ(LINK)
                                     JSTOP = XADJ(LINK+1) - 1
                                     DO  1900  J = JSTRT, JSTOP
                                         NODE = ADJNCY(J)
                                         LINK = - NODE
                                         IF  ( NODE )  1700, 2000, 1800
C
 1800                                    CONTINUE
                                         IF  ( MARKER(NODE) .GE. TAG )
     1                                         GO TO 1900
                                             MARKER(NODE) = TAG
                                             DEG = DEG + QSIZE(NODE)
 1900                                CONTINUE
 2000                    CONTINUE
 2100                CONTINUE
C                    -------------------------------------------
C                    UPDATE EXTERNAL DEGREE OF ENODE IN DEGREE
C                    STRUCTURE, AND MDEG (MIN DEG) IF NECESSARY.
C                    -------------------------------------------
                     DEG = DEG - QSIZE(ENODE) + 1
                     FNODE = DHEAD(DEG)
                     DFORW(ENODE) = FNODE
                     DBAKW(ENODE) = - DEG
                     IF  ( FNODE .GT. 0 )  DBAKW(FNODE) = ENODE
                     DHEAD(DEG) = ENODE
                     IF  ( DEG .LT. MDEG )  MDEG = DEG
 2200                CONTINUE
C                    ----------------------------------
C                    GET NEXT ENODE IN CURRENT ELEMENT.
C                    ----------------------------------
                     ENODE = LLIST(ENODE)
                     IF  ( IQ2 .EQ. 1 )  GO TO 900
                         GO TO 1600
 2300        CONTINUE
C            -----------------------------
C            GET NEXT ELEMENT IN THE LIST.
C            -----------------------------
             TAG = MTAG
             ELMNT = LLIST(ELMNT)
             GO TO 100
C
      END
C***************************************************************
C***************************************************************
C*****     MMDNUM ..... MULTI MINIMUM DEGREE NUMBERING     *****
C***************************************************************
C***************************************************************
C
C     AUTHOR - JOSEPH W.H. LIU
C              DEPT OF COMPUTER SCIENCE, YORK UNIVERSITY.
C
C     PURPOSE - THIS ROUTINE PERFORMS THE FINAL STEP IN
C        PRODUCING THE PERMUTATION AND INVERSE PERMUTATION
C        VECTORS IN THE MULTIPLE ELIMINATION VERSION OF THE
C        MINIMUM DEGREE ORDERING ALGORITHM.
C
C     INPUT PARAMETERS -
C        NEQNS  - NUMBER OF EQUATIONS.
C        QSIZE  - SIZE OF SUPERNODES AT ELIMINATION.
C
C     UPDATED PARAMETERS -
C        INVP   - INVERSE PERMUTATION VECTOR.  ON INPUT,
C                 IF QSIZE(NODE)=0, THEN NODE HAS BEEN MERGED
C                 INTO THE NODE -INVP(NODE); OTHERWISE,
C                 -INVP(NODE) IS ITS INVERSE LABELLING.
C
C     OUTPUT PARAMETERS -
C        PERM   - THE PERMUTATION VECTOR.
C
C***************************************************************
C
      SUBROUTINE  MMDNUM ( NEQNS, PERM, INVP, QSIZE )
C
C***************************************************************
C
C         INTEGER*2  INVP(1)  , PERM(1)  , QSIZE(1)
         INTEGER*4  INVP(1)  , PERM(1)  , QSIZE(1)
         INTEGER*4  FATHER, NEQNS , NEXTF , NODE  , NQSIZE,
     1              NUM   , ROOT
C
C***************************************************************
C
         DO  100  NODE = 1, NEQNS
             NQSIZE = QSIZE(NODE)
             IF  ( NQSIZE .LE. 0 )  PERM(NODE) = INVP(NODE)
             IF  ( NQSIZE .GT. 0 )  PERM(NODE) = - INVP(NODE)
  100    CONTINUE
C        ------------------------------------------------------
C        FOR EACH NODE WHICH HAS BEEN MERGED, DO THE FOLLOWING.
C        ------------------------------------------------------
         DO  500  NODE = 1, NEQNS
             IF  ( PERM(NODE) .GT. 0 )  GO TO 500
C                -----------------------------------------
C                TRACE THE MERGED TREE UNTIL ONE WHICH HAS
C                NOT BEEN MERGED, CALL IT ROOT.
C                -----------------------------------------
                 FATHER = NODE
  200            CONTINUE
                     IF  ( PERM(FATHER) .GT. 0 )  GO TO 300
                         FATHER = - PERM(FATHER)
                         GO TO 200
  300            CONTINUE
C                -----------------------
C                NUMBER NODE AFTER ROOT.
C                -----------------------
                 ROOT = FATHER
                 NUM = PERM(ROOT) + 1
                 INVP(NODE) = - NUM
                 PERM(ROOT) = NUM
C                ------------------------
C                SHORTEN THE MERGED TREE.
C                ------------------------
                 FATHER = NODE
  400            CONTINUE
                     NEXTF = - PERM(FATHER)
                     IF  ( NEXTF .LE. 0 )  GO TO 500
                         PERM(FATHER) = - ROOT
                         FATHER = NEXTF
                         GO TO 400
  500    CONTINUE
C        ----------------------
C        READY TO COMPUTE PERM.
C        ----------------------
         DO  600  NODE = 1, NEQNS
             NUM = - INVP(NODE)
             INVP(NODE) = NUM
             PERM(NUM) = NODE
  600    CONTINUE
         RETURN
C
      END

