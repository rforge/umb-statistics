# $Id$

#
# RcmdrPlugin.UMB - Extensions of the R Commander for statistical teaching
#                       at the Norwegian University of Life Sciences.
#
# Organisation:
# UMB.Utilities.R  - small functions and look-ups
# UMB.GUI.Data.R   - graphical user interface functions
# UMB.GUI.Statistics.R - graphical user interface functions
# UMB.GUI.Models.R - graphical user interface functions
# UMB.GUI.File.R   - graphical user interface functions
# UMB.GUI.Graphs.R - graphical user interface functions
# UMB.Graphics.R   - plot functions
# UMB.Statistics.R - statistical functions

.onAttach <- function(libname, pkgname){
    if (!interactive()) return()
    Rcmdr <- options()$Rcmdr
    plugins <- Rcmdr$plugins
    if (!pkgname %in% plugins) {
        Rcmdr$plugins <- c(plugins, pkgname)
        options(Rcmdr=Rcmdr)
        if("package:Rcmdr" %in% search()) {
            if(!getRcmdr("autoRestart")) {
                closeCommander(ask=FALSE, ask.save=TRUE)
                Commander()
            }
        }
        else {
            Commander()
        }
    }
	putRcmdr("allFileName",NULL) # Prepare variable on package load
	options(contrasts=c('contr.sum','contr.poly'))
	tkbind(CommanderWindow(), "<Control-e>", onExportE)
	tkbind(CommanderWindow(), "<Control-E>", onExportN)
}
