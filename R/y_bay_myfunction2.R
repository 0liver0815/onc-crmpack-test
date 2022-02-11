##' A Reference Class to represent sequentially updated reporting objects.
##' @name Report
##' @field object The object from which to report
##' @field df the data frame to which columns are sequentially added
##' @field dfNames the names to which strings are sequentially added
ReportBay <-
  setRefClass("ReportBay",
              fields =
                list(object = "ANY",
                     df = "data.frame",
                     dfNames = "character"),
              methods = list(
                dfSave =
                  function(res, name) {
                    df <<- cbind(df, res)
                    dfNames <<- c(dfNames, name)
                    return(res)
                  },
                report =
                  function(slotName,
                           description,
                           percent=TRUE,
                           digits=0,
                           quantiles=c(0, 0.1, 0.9, 1),
                           subset=NULL,
                           doSum=FALSE) {
                    vals <- slot(object, name=slotName)
                    if(! is.null(subset))
                      vals <- vals[subset,]  
                    if(doSum)
                      vals <- apply(vals, 2, sum)  
                    if(percent)
                    {
                      unit <- " %"
                      vals <- vals * 100
                    } else {
                      unit <- ""
                    }
                    
                    res <- paste(round(mean(vals), digits),
                                 unit,
                                 " (",
                                 paste(round(quantile(vals,
                                                      quantiles,
                                                      na.rm=TRUE),
                                             digits),
                                       unit,
                                       collapse=", ",
                                       sep=""),
                                 ")",
                                 sep="")
                    
                    ## print result to the buffer
                    cat(description, ":",
                        "mean",
                        dfSave(res, slotName),
                        "\n")
                  }))