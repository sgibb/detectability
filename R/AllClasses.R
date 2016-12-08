.AAIndex <- setClass("AAIndex",
                     contains="eSet",
                     prototype=prototype(new("VersionedBiobase",
                                             versions=c(Biobase::classVersion("eSet"),
                                             AAIndex="0.1"))))
