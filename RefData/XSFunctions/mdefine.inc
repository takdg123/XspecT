C mdefine.inc 

c MAXMDEF              i :  maximum number of models
c MAXPAR               i :  maximum number of parameters for each model
c MAXLXP               i :  maximum length of mexpr
c MAXLMN               i :  maximum length of mname  
c mnam(MAXMDEF)        c :  model name array 
c mexpr (MAXMDEF)      c :  model expression (preprocessed)
c mtype(MAXMDEF)       i :  model type (defined in MODEL/mtype.inc) 
c mnpar(MAXMDEF)       i :  number of parameters 
c memin(MAXMDEF)       r :  min energy
c memax(MAXMDEF)       r :  max energy
c NMDEF                i :  number of models already defined
c MID(MAXMDEF)         i :  index to model 
c mqlc(MAXMDEF)        l :  for convolution model only, true, if position
c                           dependent, (ie .e special parameters are used) 
c
c mpval(MAXPAR,MAXMDEF)    r :  value of parameter/model  
c mpsoft1(MAXPAR,MAXMDEF)  r :  soft low limit of a parameter/model 
c mpsoft2(MAXPAR,MAXMDEF)  r :  soft high  limit of a parameter/model 
c mphard1(MAXPAR,MAXMDEF)  r :  hard low limit of a parameter/model 
c mphard2(MAXPAR,MAXMDEF)  r :  hard high limit of a parameter/model 
c mpdelt(MAXPAR,MAXMDEF)   r :  step delta of parameter/model
c mpsigm(MAXPAR,MAXMDEF)   r :  sigma of parameter/model
c mpnam(MAXPAR,MAXMDEF)    C :  name of the parameter/model
c mpunt(MAXPAR,MAXMDEF)    C :  unit of the parameter/model

c munit(MAXPAR,MAXMDEF) c :  parameter unit 
c

       INTEGER MAXMDEF, MAXPAR, MAXLXP, MAXLMN

       PARAMETER (MAXMDEF=50, MAXPAR=10, MAXLXP=512, MAXLMN=8) 
       CHARACTER  MNAME(MAXMDEF)*(MAXLMN), MEXPR(MAXMDEF)*(MAXLXP),   
     &          mpnam(MAXPAR,MAXMDEF)*(MAXLMN), 
     &          mpunt(MAXPAR,MAXMDEF)*(MAXLMN)

       LOGICAL mqlc(MAXMDEF)

       REAL    memin(MAXMDEF),memax(MAXMDEF), mpval(MAXPAR,MAXMDEF),
     &    mpsoft1(MAXPAR,MAXMDEF),mpsoft2(MAXPAR,MAXMDEF),
     &    mphard1(MAXPAR,MAXMDEF),mphard2(MAXPAR,MAXMDEF),
     &    mpdelt(MAXPAR,MAXMDEF), mpsigm(MAXPAR,MAXMDEF) 
       
       INTEGER mnpar(MAXMDEF), MTYPE(MAXMDEF)   
       INTEGER NMDEF, MID(MAXMDEF)
       COMMON  /mmdef1/mname, mexpr
       COMMON  /mmdef2/memax,memin, nmdef,mid, mtype, mqlc
       COMMON  /mpdef1/mpnam,mpunt
       COMMON  /mpdef2/mnpar, mpval,mpsoft1,mpsoft2,mphard1,mphard2,
     &                mpdelt,mpsigm  
c
c---------------------------------------------------------------------
c Model evaluation related data:
c 
c  STAP connection:  W(I) ---> MEMR(PTR(I))
c
c  MAXNWK  = max number of working pointers 
c  NWK    = number of working pointers actually allocated
c  N(I)   = dimension of the I-th vector 
c            (could be p1, p2, e, and working vectors)
c  PTR(I) = pointer of I-th vector  (I=0, MAXPAR+NWK)
c
c XSPEC: defined names:   E and P1, P2 ....
c ID table:  0 = E, e
C            1 = P1, p1     (first model parameter)
C            2 = P2, p2     (second model parameter)
C            ..........
C           10 = P10, p10
C
C   J, .. MAXCOL =  working vectors
C            ...
c
C   0
c
      INTEGER MAXNWK, MAXCOL

      PARAMETER (MAXNWK = 20, MAXCOL=MAXPAR+MAXNWK)        
      INTEGER N(0:MAXCOL), PTR(0:MAXCOL), NWK 
      COMMON /XSTAP/N, PTR, NWK


