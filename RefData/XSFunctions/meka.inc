
c Program parameters and arrays...

      INTEGER NUMION
      PARAMETER(NUMION=207)

      INTEGER NREC, NAMAX
      PARAMETER(NREC=23, NAMAX=NREC*256)

      INTEGER NSTO
      PARAMETER (NSTO=3000)

      INTEGER L_CAFE
      PARAMETER (L_CAFE=300)

      INTEGER ne, ncafe

c       ne        from Fillin
c       ncafe     from Filcaf

      INTEGER NU, NG
      PARAMETER (NU=25, NG=25)

c       from Filrec...

      REAL arec(13, NUMION)

c       from Filion...

      REAL aion(NAMAX) 

c       from Fillin

      INTEGER   nrl(NSTO)
      INTEGER   lp(9,NSTO)
      character(8) trans(NSTO)
      REAL        rp(14,NSTO)

c       from Filcaf

      INTEGER   idnr(L_CAFE)
      REAL        cdrcafe(L_CAFE)

c       from Filspe

      REAL p(6,NUMION)

c       from Filkar

      REAL ga(NG,NU)

c  Common blocks for all this data

      COMMON /rmekdat/ arec, aion, rp, cdrcafe, p, ga
      COMMON /imekdat/ ne, ncafe
      COMMON /i2mekdat/ nrl, lp, idnr
      COMMON /cmekdat/ trans
