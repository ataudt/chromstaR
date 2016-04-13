# # ====================================================
# # Estimate the RAM consumption of the multivariate HMM
# # ====================================================
# estimate.RAM.consumption <- function(num.states, num.bins, num.mods) {
# # scalefactoralpha T
# # scalealpha N,T
# # scalebeta N,T
# # densities N,T
# # gamma N,T
# # A N,N
# # sumxi N,N
# # gammaold N,T
# # multiO Nmod,T int
# # binary_states N,Nmod bool
# 
# # double 8bytes
# # int 4bytes
# # bool 4bytes
#     ram <- 8 * (num.bins + num.bins*num.states *5 + num.states^2 *2) + 4 * (num.mods*num.bins) + 4 * (num.states*num.mods)
#     ram.Mb <- as.integer(ram / 2^20)
#     return(ram.Mb)
# }
