#' @importFrom Biobase assayDataElement
#' @noRd
setMethod("exprs", signature(object="AAIndex"),
          function(object) assayDataElement(object, "exprs"))
