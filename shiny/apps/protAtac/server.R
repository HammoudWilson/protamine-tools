#----------------------------------------------------------------------
# appServer() is called in session context and thus has access to:
#   input, output, session objects
#   values returned from app step modules
#----------------------------------------------------------------------

# objects instantiated here are available to all appStep modules in a session

# session cache objects
binsCache          <- new_dataCache('binsCache')
insertSizesCache   <- new_dataCache('insertSizesCache')
scoresCache        <- new_dataCache('scoresCache')
gcrzCache          <- new_dataCache('gcrzCache')
segmentationCache  <- new_dataCache('segmentationCache')

# track and other reactives
binsWorkingSourceId   <- reactiveVal(NULL)
scoresWorkingSourceId <- reactiveVal(NULL)
gcrzWorkingSourceId   <- reactiveVal(NULL)

# appServer() is called after all modules are instantiated
appServer <- function(){
    # objects instantiated here are available to this app step only
}
