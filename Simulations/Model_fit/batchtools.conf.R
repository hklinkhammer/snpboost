# cluster.functions <- makeClusterFunctionsInteractive()

# cluster.functions <- makeClusterFunctionsSSH(
#   
#   makeSSHWorker(nodename="xcat",
#                 ncpus=16, 
#                 max.jobs=10)
# 
# )
# 
# debug = FALSE

# cluster.functions <- makeClusterFunctionsTORQUE(template = findTemplateFile("torque"),
#                                                 scheduler.latency = 10, fs.latency = 65)
# 
# debug = TRUE 

cluster.functions <- makeClusterFunctionsSlurm(template = "slurm", array.jobs = TRUE,
                                               scheduler.latency = 1, fs.latency = 65, nodename = "localhost")