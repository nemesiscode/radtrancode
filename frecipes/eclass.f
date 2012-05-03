      SUBROUTINE ECLASS(NF,N,LISTA,LISTB,M)
      DIMENSION NF(N),LISTA(M),LISTB(M)
      DO 11 K=1,N
        NF(K)=K
11    CONTINUE
      DO 12 L=1,M
        J=LISTA(L)
1       IF(NF(J).NE.J)THEN
          J=NF(J)
        GOTO 1
        ENDIF
        K=LISTB(L)
2       IF(NF(K).NE.K)THEN
          K=NF(K)
        GOTO 2
        ENDIF
        IF(J.NE.K)NF(J)=K
12    CONTINUE
      DO 13 J=1,N
3       IF(NF(J).NE.NF(NF(J)))THEN
          NF(J)=NF(NF(J))
        GOTO 3
        ENDIF
13    CONTINUE
      RETURN
      END
