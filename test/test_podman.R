library(Zorn)
bascetInstance.default <- getBascetPodmanImage(storeAt="~/Downloads/", verbose = TRUE) #, forceInstall = TRUE)
TestBascetInstance(bascetInstance.default)

BascetCellNames(
  bascetRoot,
  "abricate",
  verbose=FALSE
)
