
# constants tRNA ---------------------------------------------------------------

TRNA_FEATURES <- c(
  "tRNA_length",
  "tRNA_type",
  "tRNA_anticodon",
  "tRNA_seq",
  "tRNA_str",
  "tRNA_CCA.end"
)

TRNA_STRUCTURES <- c(
  "acceptorStem",
  "Dprime5",
  "DStem",
  "Dloop",
  "Dprime3",
  "anticodonStem",
  "anticodonLoop",
  "variableLoop",
  "TStem",
  "Tloop",
  "discriminator"
)

tRNAStructureFunctionList <- list(
  acceptorStem = ".getAcceptorStem",
  Dprime5 = ".getDprime5",
  DStem = ".getDstem",
  Dloop = ".getDloop",
  Dprime3 = ".getDprime3",
  anticodonStem = ".getAnticodonStem",
  anticodonLoop = ".getAnticodonLoop",
  variableLoop = ".getVariableLoop",
  TStem = ".getTstem",
  Tloop = ".getTloop",
  discriminator = ".getDiscriminator")


TRNA_STRUCTURE_ORDER <- c("acceptorStem.prime5",
                          "Dprime5",
                          "DStem.prime5",
                          "Dloop",
                          "DStem.prime3",
                          "Dprime3",
                          "anticodonStem.prime5",
                          "anticodonLoop",
                          "anticodonStem.prime3",
                          "variableLoop",
                          "TStem.prime5",
                          "Tloop",
                          "TStem.prime3",
                          "acceptorStem.prime3",
                          "discriminator")