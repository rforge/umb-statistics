%% $Id$
\encoding{latin1}
\name{RcmdrPlugin.UMB-package}
\alias{RcmdrPlugin.UMB-package}
\alias{RcmdrPlugin.UMB}
\alias{Dummify}
\alias{auto.import}
\alias{dummify}
\alias{formulaFields2}
\alias{modelFormula2}
\alias{modelFormula3}
\alias{radioButtonsUMB}
\alias{addObservationStatisticsUMB}
\alias{anova_reg_GUI}
\alias{aovP}
\alias{auto.import.GUI}
\alias{backwardDrop}
\alias{backwardForward}
\alias{bestSubsets}
\alias{binomialDistributionPlotUMB}
\alias{boxPlotUMB}
\alias{CIplot}
\alias{clustP}
\alias{confint.mvr}
\alias{Confint.mvr}
\alias{contrastGUI}
\alias{contrastGUI2}
\alias{covarianceMatrix}
\alias{createSequence}
\alias{DA.coef}
\alias{DA.scores}
\alias{daP}
\alias{deleteActiveDataSet}
\alias{deviance_tests}
\alias{deviance_tests_GUI}
\alias{discreteDistributionPlotUMB}
\alias{discriminantAnalysis}
\alias{discriminantPlot}
\alias{dotPlot}
\alias{dotplotGUI}
\alias{dots}
\alias{dummy}
\alias{enterTableUMB}
\alias{export.instructions}
\alias{extend.colnames}
\alias{Fisher}
\alias{Fisher.aov}
\alias{fittedLinePlot}
\alias{forward}
\alias{forwardAdd}
\alias{forwardBackward}
\alias{generalizedLinearModelUMB}
\alias{geomDistributionPlotUMB}
\alias{hclust.list}
\alias{hclust.merge}
\alias{hierarchicalClusterVariable}
\alias{histogram_discrete}
\alias{hyperDistributionPlotUMB}
\alias{if.R}
\alias{independentSamplesTTestUMB}
\alias{independentSamplesZTestUMB}
\alias{invXtXGUI}
\alias{linearModelUMB}
\alias{linePlotUMB}
\alias{listHclustSolutions}
\alias{loadAll}
\alias{logic.to.letters}
\alias{make.colours}
\alias{meanCenter}
\alias{mixed.model}
\alias{mixP}
\alias{mixture.contour}
\alias{mixtureGUI}
\alias{multinomialLogitModelUMB}
\alias{multiWayAnovaUMB}
\alias{negbinomialDistributionPlotUMB}
\alias{newDataSetOSX}
\alias{numericalSummariesUMB}
\alias{numSummaryUMB}
\alias{oneSamplesTTestSum}
\alias{oneSamplesZTestSum}
\alias{oneWayAnovaUMB}
\alias{onExport}
\alias{onExportE}
\alias{onExportN}
\alias{ordinalRegressionModelUMB}
\alias{pairedTTestUMB}
\alias{pairedZTestUMB}
\alias{panel.ci.plot}
\alias{pcaP}
\alias{pieChartUMB}
\alias{plot.Bonferroni}
\alias{plot.Fisher}
\alias{plot.TukeyHSD}
\alias{plotByGroups}
\alias{plotDA}
\alias{plotFitDens}
\alias{plotModelUMB}
\alias{plotprops}
\alias{plsP}
\alias{plsRegressionModel}
\alias{PoissonDistributionPlotUMB}
\alias{postHocGUI}
\alias{predict_CI_PI}
\alias{predict_link_response}
\alias{predictRegressionModel}
\alias{PRESS.GUI}
\alias{principalComponentPlots}
\alias{principalComponentsStat}
\alias{proportionTest}
\alias{RelComp}
\alias{resetGLMUMB}
\alias{resetLinearModelUMB}
\alias{saveAll}
\alias{saveAllAs}
\alias{saveAllWorker}
\alias{sdErr}
\alias{simplex.analysis}
\alias{simplexGUI}
\alias{singleSampleTTestUMB}
\alias{singleSampleZTestUMB}
\alias{sortData}
\alias{spearson}
\alias{stepWise}
\alias{stepWiseBack}
\alias{summaryMultinom}
\alias{summaryOrdinal}
\alias{sumSq}
\alias{tally}
\alias{tally.GUI}
\alias{toClipboard}
\alias{twoProportionTest}
\alias{twoSamplesTTest}
\alias{twoSamplesTTestSum}
\alias{twoSamplesZTest}
\alias{twoSamplesZTestSum}
\alias{twoWayTableUMB}
\alias{updateFactor}
\alias{variablesP}
\docType{package}
\title{
  R Commander statistics plugin (RcmdrPlugin.UMB)
}
\description{
  This package provides an Rcmdr "plug-in" extending some functions of 
  linear models and providing new and extended interfaces to PCA, PLS, LDA,
  QDA, clustering of variables, tests, plots, etc.
}
\details{
\tabular{ll}{
  Package: \tab RcmdrPlugin.UMB\cr
  License: \tab GPL version 2 or newer\cr
  LazyLoad: \tab yes\cr
}
  Run the plugin using the code \code{library(RcmdrPlugin.UMB)} or \code{library(umb)}.\cr
  Statistical and summary functions included:
  \itemize{
  \item \code{dummy(y)}
  \item \code{PRESS(object=NULL) (default is current model)}
  \item \code{PRESS.res(object=NULL) (default is current model)}
  \item \code{R2_pred(object=NULL) (default is current model)}
  \item \code{forward(model, alpha=0.2, full=FALSE)}
  \item \code{backward(model, alpha=0.2, full=FALSE)}
  \item \code{stepWise(model, alpha.enter=0.15, alpha.remove=0.15, full=FALSE)}
  \item \code{stepWiseBack(model, alpha.remove=0.15, alpha.enter=0.15, full=FALSE)}
  \item \code{best.subsets(model, nbest=5, nvmax)}
  \item \code{confint.mvr(object, parm, level=0.95, ...)}
  \item \code{confusion(true, predicted)}
  \item \code{DA.scores(object=NULL) (default is current model)}
  \item \code{plotDA(DAobject=NULL, regions=TRUE, contours=FALSE, resolution=100)}
  \item \code{hclust.merge(object) (default is last clustering)}
  \item \code{mixed.model(formula, random.effects=NULL, data, restrictedModel=FALSE, subset="")}
  \item \code{summary.extra(object)}
  \item \code{anova_reg(lm.object)}
  \item \code{predict_CI_PI(model, data, level)}
  \item \code{prop.test.ordinary(x, n, p = NULL, alternative = c("two.sided", "less", 
    "greater"), conf.level = 0.95, correct = TRUE)}
  \item \code{rmsep(object) (default is current model)}
  }
  Utility functions included:
  \itemize{
  \item \code{aovP()}
  \item \code{clustP()}
  \item \code{daP()}
  \item \code{mixP()}
  \item \code{pcaP()}
  \item \code{plsP()}
  \item \code{variablesP()}
  \item \code{DA.coef()}
  \item \code{hclust.list()}
  \item \code{listHclustSolutions(envir=.GlobalEnv, ...)}
  \item \code{make.colours(object)}
  \item \code{confint.mvr(object, parm, level=0.95, ...)}
  \item \code{dummy(y)}
  \item \code{dummify(y,n,name)}
  \item \code{Dummify(data, main.effects, response)}
  \item \code{fparse(f)}
  \item \code{if.R()}           
  }
  GUI Functions included:
  \itemize{
  \item \code{anova_reg_GUI()}              
  \item \code{backwardDrop()}               
  \item \code{backwardForward()}            
  \item \code{bestSubsets()}                
  \item \code{contrastGUI()}                
  \item \code{contrastGUI2()}               
  \item \code{covarianceMatrix()}           
  \item \code{createSequence()}             
  \item \code{discriminantAnalysis()}       
  \item \code{discriminantPlot()}           
  \item \code{dotplotGUI()}                 
  \item \code{enterTableUMB()}              
  \item \code{fittedLinePlot()}             
  \item \code{forwardAdd()}                 
  \item \code{forwardBackward()}            
  \item \code{hierarchicalClusterVariable()}
  \item \code{linearModelANOVA()}           
  \item \code{meanCenter()}                 
  \item \code{mixtureGUI()}                 
  \item \code{plsRegressionModel()}         
  \item \code{postHocGUI()}                 
  \item \code{predictRegressionModel()}     
  \item \code{PRESS.GUI()}                  
  \item \code{principalComponentPlots()}    
  \item \code{principalComponentsStat()}    
  \item \code{proportionTest()}             
  \item \code{sortData()}                   
  \item \code{twoSamplesTTest()}            
  \item \code{twoWayTableUMB()}
  }
  Graphical Functions included:
  \itemize{
  \item \code{CIplot()}         
  \item \code{dotPlot()}        
  \item \code{dots()}           
  \item \code{mixture.contour()}
  \item \code{panel.ci.plot()}  
  \item \code{plotDA()}        
  }
}
\author{
Kristian Hovde Liland <kristian.liland@umb.no>
}
\keyword{ package }
\examples{
\dontrun{
library(RcmdrPlugin.UMB) # Starts up the R Commander including this plugin.
library(umb) # Simpler startup with automatic update (when available)
}
}