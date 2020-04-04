## cca eqs
  ##Pg = G(G'G)^-1G'

  ## X = H'(G(G'G)^-1G')H
  ## W = H'H

library(GSVD)

data(wine)
X <- scale(wine$objective)
Y <- scale(wine$subjective)
R <- t(X) %*% Y

cancor_res <- cancor(X, Y)

 cca.res <- gsvd(
     DAT=(crossprod(X) %^% -1) %*% R %*% (crossprod(Y) %^% -1),
     LW=crossprod(X),
     RW=crossprod(Y)
 )
 cca.res$lx <- (X %*% cca.res$p)
 cca.res$ly <- (Y %*% cca.res$q)


 # takane_cca.res <- gsvd(
 #   ((t(X) %*% X) %^% (-1)) %*% (t(X) %*% Y) %*% ((t(Y) %*% Y) %^% (-1)),
 #   t(X) %*% X,
 #   t(Y) %*% Y
 # )
 # mapply("/", cca.res, takane_cca.res)

 cca.res_gpls <- gplssvd(
        X = X %^% (-1),
        Y = Y %^% (-1),
        XRW=crossprod(X),
        YRW=crossprod(Y)
    )


 ### this is actually the solution. not what Takane lists in his paper.
 this_mat <-
   ((t(X) %*% X) %^% (-1)) %*% (t(X) %*% Y) %*% ((t(Y) %*% Y) %^% (-1))
 takane_test <- geigen(
   this_mat %*% ((t(Y) %*% Y) %^% (1)) %*% t(this_mat),
   (t(X) %*% X)
)

 this_other_mat <-
   ((t(Y) %*% Y) %^% (-1)) %*% (t(Y) %*% X) %*% ((t(X) %*% X) %^% (-1))
 takane_other_test <- geigen(
   this_other_mat %*% ((t(X) %*% X) %^% (1)) %*% t(this_other_mat),
   (t(Y) %*% Y)
 )



 ### alright so *now* I can actually perform tests between the new and old functions...
 cca.res_2 <- gsvd_2(
   (crossprod(X) %^% -1) %*% R %*% (crossprod(Y) %^% -1),
   LW=crossprod(X),
   RW=crossprod(Y)
 )
 cca.res_2$lx <- (X %*% cca.res$p)
 cca.res_2$ly <- (Y %*% cca.res$q)

 mapply("/", cca.res, cca.res_2)

 cca.res_gpls_2 <- gplssvd_2(
   X = X %^% (-1),
   Y = Y %^% (-1),
   XRW=crossprod(X),
   YRW=crossprod(Y)
 )

 mapply("/", cca.res_gpls, cca.res_gpls_2)

 takane_test_2 <- geigen_2(
   this_mat %*% ((t(Y) %*% Y) %^% (1)) %*% t(this_mat),
   (t(X) %*% X)
 )

 mapply("/", takane_test, takane_test_2)

 takane_other_test_2 <- geigen_2(
   this_other_mat %*% ((t(X) %*% X) %^% (1)) %*% t(this_other_mat),
   (t(Y) %*% Y)
 )

 mapply("/", takane_other_test, takane_other_test_2)




#  svd(
#    ((Y %*% ((t(Y) %*% Y)%^%(-1)) %*% t(Y) )) %*% ((X %*% ((t(X) %*% X)%^%(-1)) %*% t(X) ))
#    )
#
#  Ryy <- crossprod(Y)
#  Ryx <- crossprod(Y,X)
#  Rxx <- crossprod(X)
#  Rxy <- crossprod(X,Y)
#
# eigen( (Ryy %^% (-1)) %*% Ryx %*% (Rxx %^% (-1)) %*% Rxy )




 #
 # takane_cca.evd_Y <- geigen(
 #   ((t(Y) %*% (X %*% ((t(X) %*% X)%^%(-1)) %*% t(X) ) %*% Y)),
 #   ((t(Y) %*% Y) %^% (-1))
 # )
 #
 # takane_cca.evd_X <- geigen(
 #   t(X) %*% (Y %*% ((t(Y) %*% Y)%^%(-1)) %*% t(Y) ) %*% X,
 #   ((t(X) %*% X) %^% (-1))
 # )
 #
