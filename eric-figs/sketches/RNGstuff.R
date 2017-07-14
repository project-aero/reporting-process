## Random Number Generator stuff for reproducibility

## To Do: make this into a function and put in helpers.R

# Save current random number generator state to file
if (!exists("./data/.Random.seed")){set.seed(1)}
RNGversion <- as.character(getRversion())
RNGkind <- RNGkind()
RNGseed <- .Random.seed
currentseedfile <- format(Sys.time(), "RNGstate%Y%m%d_%H%M%S.Rdata")
save("RNGversion", "RNGkind","RNGseed", file=currentseedfile)
print(currentseedfile)

# Restore random number generator state from file
load("./data/RNGstate20170602_151509.Rdata")
RNGversion(RNGversion) # enforces RNG version
do.call("RNGkind",as.list(RNGkind)) # enforces RNG kind
.Random.seed <- RNGseed # restore seed

# enforce RNG kind and RNG version
RNGkind("Mersenne-Twister", "Inversion")
RNGversion("3.3.2")


# possible structre of functions (note need to assign to global environment):

# save_rng <- function(savefile=tempfile()) {
#     if (exists(".Random.seed"))  {
#         oldseed <- get(".Random.seed", .GlobalEnv)
#     } else stop("don't know how to save before set.seed() or r*** call")
#     oldRNGkind <- RNGkind()
#     save("oldseed","oldRNGkind",file=savefile)
#     invisible(savefile)
# }
# 
# restore_rng <- function(savefile) {
#     load(savefile)
#     do.call("RNGkind",as.list(oldRNGkind))  ## must be first!
#     assign(".Random.seed", oldseed, .GlobalEnv)
# }
