C***********************************************************************
C_TITL: PARTFCOM
C
C_DESC:	Partition fucntion file.
C
C_ARGS: See definitions below.
C
C_CALL: No calls.
C
C_HIST: 
C***********************************************************************

      INTEGER MPARTEXTRA,MTABEXTRA,IREADEXTRA,NPARTEXTRA
      PARAMETER (MPARTEXTRA=20,MTABEXTRA=100)
C MPARTEXTRA: maximum number of gases that have tabulated partition functions
C defined.
C MTABEXTRA: maximum number temperatures for the partition functions defined
C for each gas.
      INTEGER IDEXTRA(MPARTEXTRA),ISOEXTRA(MPARTEXTRA)
      INTEGER NTABEXTRA(MPARTEXTRA)
      REAL TEMPEXTRA(MPARTEXTRA,MTABEXTRA)
      REAL PARTFEXTRA(MPARTEXTRA,MTABEXTRA)

      COMMON /PARTFEXTRA/IREADEXTRA,NPARTEXTRA,IDEXTRA,ISOEXTRA,
     &  NTABEXTRA,TEMPEXTRA,PARTFEXTRA

