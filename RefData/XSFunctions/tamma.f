      SUBROUTINE TAMMA(X, GAM)

      DOUBLE PRECISION SUM(6), A(5)
      DOUBLE PRECISION Z, GAM, X
      INTEGER J, N

      DATA A /-0.5748646D0,  0.9512363D0, -0.6998588D0,
     &         0.4245549D0, -0.1010678D0 /

      SUM(1) = 1.0D0
      Z = DABS(X)
      IF (DABS(X).GT.1.D0) Z = DABS(X-DINT(X))
      IF (X.LT.0.D0) Z = 1.D0 - Z
      DO N = 1, 5
         SUM(N+1) = SUM(N) + A(N)*Z**N
      ENDDO
      GAM = SUM(6)/Z
      IF (X.GT.1.D0) THEN
         N = INT(X)
         DO J = 1, N
            GAM = GAM*(X-J)
         ENDDO
      ELSEIF (X.LT.0.D0) THEN
         N = INT(DABS(X)) + 1
         DO J = 1, N
            GAM = GAM/(X+(J-1.D0))
         ENDDO
      ENDIF

      RETURN
      END
