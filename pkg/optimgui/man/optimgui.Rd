\name{optimgui}
\alias{optimgui}
\docType{package}
\title{
optimgui
}
\description{
A graphical user interface to assist R users to set up problems that require optimization 
tools for their solution. At the time of this documentation (2011-8-6) we restrict tools to
unconstrained or box-constrained function minimization.
}
\details{
The principal function of optimgui is called optimgui(). 
It is loaded with
   library(optimgui) 
or 
   require(optimgui)
and invoked with

   optimgui()

at the R prompt.

Functions from events.R

   Events for "CatalogEntry"   

onAddEntry = function(h, ...)
onDeleteEntry = function(h, ...)
onSelectAll = function(h, ...)
onToggleCell = function(cell, path, param)
onEditCell = function(cell, path, new.text, param)


Events for menu items      
onExit = function(h, ...) dispose(h$action$mainWin);
onOpenCatalog = function(h, ...)
onDefaultEvent = function(h, ...) gmessage("Not implemented yet. :(", "Message");
openRopFile = function(filePath, param)
onOpenBuiltinRopFile = function(h, ...)
onOpenRopFile = function(h, ...)
onSaveRopFile = function(h, ...)
onCloseRopFile = function(h, ...)


Events for notebook widget

onFocusTab = function(h, ...)
onAddTab = function(h, ...)
onDeleteTab = function(h, ...)
onTestFunction = function(h, ...)
onRunCode = function(h, ...)
onSaveOutput = function(h, ...)
onClearOutput = function(h, ...)
toggleShowNote = function(h, ...)
toggleShowCode = function(h, ...)
welcomePageAddEvent = function(welcomePage, param)
catalogPageAddEvent = function(catalogPage, param)
onChooseButton = function(h, ...) # JN

}
\author{
Yixuan Qiu \email{yixuan.qiu@cos.name}}
\keyword{ optimization }
\keyword{ graphical user interface }
\seealso{
 optimx
}
\examples{
cat("The examples for optimgui are built into the GUI itself.\n")
}
