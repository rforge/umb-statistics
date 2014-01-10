# $Id$

#
# RcmdrPlugin.NMBU - Extensions of the R Commander for statistical teaching
#                       at the Norwegian University of Life Sciences.
#
# Organisation:
# NMBU.Utilities.R  - small functions and look-ups
# NMBU.GUI.Data.R   - graphical user interface functions
# NMBU.GUI.Statistics.R - graphical user interface functions
# NMBU.GUI.Models.R - graphical user interface functions
# NMBU.GUI.File.R   - graphical user interface functions
# NMBU.GUI.Graphs.R - graphical user interface functions
# NMBU.Graphics.R   - plot functions
# NMBU.Statistics.R - statistical functions

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
