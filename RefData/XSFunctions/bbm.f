      REAL FUNCTION BBM(E, T)
C
C   MODIFIED BLACK BODY
C
      REAL x, e, t
      X = E/T
      BBM = X**1.5*EXP(-X)
      BBM = BBM/(1-EXP(-X))**0.5
      BBM = BBM/E
      RETURN
      END
