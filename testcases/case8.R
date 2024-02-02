# compute mixture of effective errors in a multi-group SEM
# compute it either based on multiple indicator theorem
# or exact NCP mixture
errs<-sapply(seq(50,950,10), function(x, e1=300, e2=10) {
  ce <- cmp_efferr_mix(var_true = 0.08,e1,e2,x,1000-x)
  N <- 1000
  N1 <- x
  N2 <- 1000-x
  ce2 <- N/(N1/e1+N2/e2)
  return(c(ce,ce2))
})

# the MI theoreom approximation gets worse the smaller the true 
# variance is in relation to the effective errors (and with smaller N)

# plot differences
library(tidyverse)
t(errs) %>% as_tibble() %>% rename(MI=V1, NCP=V2) %>% add_column(x=1:ncol(errs)) %>% 
  pivot_longer(1:2) %>% ggplot(aes(y=value,x=x,group=name,color=name ))+geom_line()
