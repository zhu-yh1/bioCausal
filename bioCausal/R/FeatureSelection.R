#' Feature selection using chosen focal variable
#'
#' @param data data matrix (feature * sample)
#' @param focalVariable the interested focal variable used for feature selection
#' @param contr provide contrast for limma if focal variable is binomial
#' @param limmaCutoff number of top selected features when using limma for top differential features
#' @param glmnetCutoff number of top selected features when using multi-elastic net regression for top prediction features
#' @param class class of focal variables, only accepts "binomial" or "gaussian"
#' @param nfolds nfolds for runGnetMulti
#' @param bootstrap bootstrap for runGnetMulti
#' @param alpha alpha for runGnetMulti
#'
#' @return A list for feature selection results.
#'
#' @export
featureSelection = function (data, focalVariable, contr = NULL, limmaCutoff = NULL, glmnetCutoff = NULL, class = NULL, nfolds = 3, bootstrap = T, alpha = 0.9, ...) {
	# limma
	if (class == "binomial") {
		limRes = fitLimma(data, as.factor(focalVariable), contr = contr)
	}
	else if (class == "gaussian")
	{
		limRes = fitLimma(data, focalVariable)
	}

	if(!is.null(limmaCutoff)) {
		limSelected = rownames(topTable(limRes, n=limmaCutoff))
	}
	else {
		limSelected = rownames(data)
	}


	# elastic net
	resGN = runGnetMulti(x = t(data), y = focalVariable, nfolds = nfolds, bootstrap = bootstrap, alpha = alpha, ...)

	if (!is.null(glmnetCutoff)) {
		cutoff = sort(resGN$features, decreasing = T)[glmnetCutoff]
		resGNSelected = rownames(data)[which(resGN$features >= cutoff)]
	}
	else {
		resGNSelected = rownames(data)
	}
	selected = rownames(data)[rownames(data) %in% c(limSelected, resGNSelected)]
	results = list(limRes = limRes, limSelected = limSelected, resGN = resGN, resGNSlected = resGNSelected, selected = selected)
	class(results) = "featureSelectionResults"
	results
}

#' Feature selection using chosen focal variable
#'
#' @param labels labels for the features
#' @param res the feature selection result list returned from featureSelection()
#' @param sign whether to plot signed values
#' @param yvalue pval or adj.pval for y-axis
#'
#' @return A list for feature selection results.
#'
#' @export
plotSelection = function(labels, res, sign = F, yvalue = "pval") {
	col = rep("Unselected", length(labels))
	col[labels %in% res$selected] = c("Selected")

  if (yvalue == "pval") {
    y = -log(max(res$limRes$p.value[res$limSelected,]))
  } else if (yvalue == "adj.pval") {
    y = -log(max(res$limRes$q.value[res$limSelected,]))
  } else {
    stop("Please set yvalue as 'pval' or 'adj.pval'")
  }

	x = sort(res$resGN$features, decreasing = T)[length(res$resGNSlected)]
	if(sign) {
	  if (yvalue == "pval") {
	    pval = -log(res$limRes$p.value) * sign(res$limRes$t)
	  } else {
	    pval = -log(res$limRes$q.value) * sign(res$limRes$t)
	  }

		p=myggscatter(res$resGN$features, pval, line = F,
					  label = labels,
					  col = col,
					  palette = c("#ff0000", "#bebebe"),
					  xlab="Chosen frequency in elastic-net regressions",
					  ylab=paste0("signed -log(", yvalue, ") from Limma"),
					  label.select = res$selected,
					  font.label = c(10, "black"),
					  title = paste0("models > ", x-1, " or abs(-log(", yvalue, ")) > ", round(y, digits=2))) +
			geom_hline(yintercept = c(-round(y, digits=2), round(y, digits = 2)), linetype="dashed", color = "red", size=1) +
			geom_vline(xintercept = x-1, linetype="dashed", color = "red", size=1)
	}
	else {
	  if (yvalue == "pval") {
	    pval = -log(res$limRes$p.value)
	  } else {
	    pval = -log(res$limRes$q.value)
	  }
		p=myggscatter(res$resGN$features, pval, line = F,
					  label = labels,
					  col = col,
					  palette = c("#ff0000", "#bebebe"),
					  xlab="Chosen frequency in elastic-net regressions",
					  ylab=paste0("-log(", yvalue, ") from Limma"),
					  label.select = res$selected,
					  font.label = c(10, "black"),
					  title = paste0("models > ", x-1, " or -log(", yvalue, ") > ", round(y, digits=2))) +
			geom_hline(yintercept = round(y, digits=2), linetype="dashed", color = "red", size=1) +
			geom_vline(xintercept = x-1, linetype="dashed", color = "red", size=1)
	}
	p
}
