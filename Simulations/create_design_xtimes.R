require(tidyverse)

# input
scene = "heritability_50_sparsity_01_lr_01"
folder = paste0("Simulations_",scene)
p = 100000

# set Design festlegen

SNP_batch = unique(as.integer(c(0.01/100*p,0.1/100*p,1/100*p,1/20*p)))

#prozentual von SNP_batch
m_batch=c(1)

learning_rate=c(0.1)

b_max=c(2,10)

design=data.frame(p_batch=rep(SNP_batch,each=length(m_batch)*length(learning_rate)*length(b_max)),
                  m_batch=rep(rep(m_batch,each=length(learning_rate)*length(b_max)),times=length(SNP_batch)),
                  learning_rate=rep(learning_rate,each=length(b_max),times=length(SNP_batch)*length(m_batch)),
                  b_max=rep(b_max,times=length(learning_rate)*length(SNP_batch)*length(m_batch)))

#design <- design %>% mutate(b_max=unlist(lapply(100/SNP_batch,function(x){max(10,x)})))

design$m_batch=as.integer(design$p_batch*design$m_batch)

design <- design %>% filter(m_batch > 0)

save(design,file=paste0(folder,"/designmatrix.RData"))
