setMethod("exprs", signature(object="AAIndex"),
          function(object) assayDataElement(object, "exprs"))
