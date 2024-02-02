# we need to solve this equation
#  N*(a/b-log(1+a/b))=x for a
# where
# a = effective error of the minimal model
# b = true variance of the latent variable of interest (which is subject to a test)
# x = non-centrality estimate based on von Oertzen & Brandmaier (2013) appendix
#
#
# result is obtained using Lambert's W function (also known as productlog)
#
# a = b*-W(-exp(-1 - x/N))




ncp_minmodel <- function(var_true, var_eff, N) {
  N*(var_true/var_eff - log(1+var_true/var_eff))
}

# according to B.6 in von Oertzen & Brandmaier (2013)
compute_ncp_mix <- function(var_true, eff1, eff2, N1, N2) {
  e1<-ncp_minmodel(var_true, eff1, N1)
  e2 <- ncp_minmodel(var_true, eff2, N2)
  e1+e2
}

cmp_efferr_mix <- function (
    var_true = 1,
    var_eff1 = 28,
    var_eff2 = 18,
    N1 = 500,
    N2 = 500
){
  N <- N1+N2
  ncp_mix <- compute_ncp_mix(var_true, var_eff1, var_eff2, N1, N2)
  
  ncp_min <- ncp_minmodel(var_true, 21.4, N)


  var_true/(-lamW::lambertWm1(-exp(-ncp_mix/N-1))-1)

}

 