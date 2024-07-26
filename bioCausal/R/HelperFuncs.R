#' scale data to mean=0 and sd=1
#'
#' @param x data matrix
#' @param margin margin=1 for rows, margin=2 for columns. default margin=1
#'
#' @return A list for feature selection results.
#'
#' @export
tscale = function(x, margin = 1){
	s = apply(x, margin, sd)
	m = apply(x, margin, mean)
	x = sweep(x, margin, m)
	x = sweep(x, margin, s, "/")
	x
}
#' get common rows for 2 dataframes or matrices
#'
#' @export
commonRows = function(data1, data2){
	intersect(rownames(data1), rownames(data2))
}

#' get common columns for 2 dataframes or matrices
#'
#' @export
commonColumns = function(data1, data2){
	intersect(colnames(data1), colnames(data2))
}

#' Calculate q-value. If reports Q-value error, adjusted p-value with BH adjustment will be returned
#'
#' @export
QV = function(pval){

	x = try(qvalue(pval))

	if(!is.list(x)){
		warning("Q-value error")
		# hist(pval)
		return(p.adjust(pval, method = "BH"))
	}
	else{
		return(x$qvalue)
	}
}

#' Add q-value for limma fit results
#'
#' @export
addQV = function(fit){
	fit$q.value = fit$p.value
	for (i in 1:ncol(fit$p.value)){
		ii = which(!is.na(fit$p.value[,i]))
		fit$q.value[ii,i] = QV(fit$p.value[ii,i])
	}
	return(fit)
}

nonEstimable = function(x){
	x = as.matrix(x)
	p = ncol(x)
	QR = qr(x)
	if (QR$rank < p){
		n = colnames(x)
		if (is.null(n)){
			n = as.character(1:p)
		}
		notest = n[QR$pivot[(QR$rank + 1):p]]
		blank = notest == ""
		if (any(blank)){
			notest[blank] = as.character(((QR$rank + 1):p)[blank])
		}
		return(notest)
	}
	else{
		return(NULL)
	}
}

getNonZeroBetas=function(gres,index=T, i=NULL, se=F,vector=F, intercept=F){
	if(is.null(i)){
		i=getBestIndex(gres, se=se)
		message(paste("index is",i))
	}

	ii=which(gres$glmnet.fit$beta[,i]!=0)
	if(index){
		ii
	}
	else if(vector){
		if(intercept){
			(c(gres$glmnet.fit$a0[i],gres$glmnet.fit$beta[,i, drop=F]))
		}
		else{
			gres$glmnet.fit$beta[,i, drop=F]
		}
	}
	else{
		gres$glmnet.fit$beta[ii,i, drop=F]
	}
}

#' fit limma model
#'
#' @param data input data (feature * samples)
#' @param grp data group information
#' @param contr manually set contrast, default=NULL
#'
#' @export
fitLimma = function(data, grp, contr = NULL, extra = NULL, pathway = NULL, method = "ls", plotHists = F, trend = F, subset = NULL){

	if(is.null(contr)) {
		grp = model.matrix(~1+grp)
		colnames(grp)[2] = "con"
		contr = "con-0"
	}

	if(is.null(dim(grp))){
		mod = model.matrix(~0+grp)
		colnames(mod) = levels(grp)
	}
	else{
		mod = grp
	}

## question: What is the meaning of extra
	if(!is.null(extra)){
		mod = cbind(mod, model.matrix(~0+extra))
	}
	colnames(mod)=make.names(colnames(mod))
	if(!is.null(subset)){
		mod = mod[subset,]
		data = data[,subset]
	}
	ne = nonEstimable(mod)
	if(!is.null(ne)){
		iirm = match(ne, colnames(mod))
		mod = mod[,-iirm]
	}

	fit = lmFit(data, mod, method = method)
	contrastMat = makeContrasts(contrasts = contr, levels = mod)
	fitC = contrasts.fit(fit, contrasts = contrastMat)
	fitC = eBayes(fitC, trend = trend, robust=T)
	fitC = addQV(fitC)
	if(plotHists){
		par(mfrow = c(2, 3))
		for(i in 1:ncol(fitC$p.value)){
			hist(fitC$p.value[,i], main = contr[i])
		}
	}
	if(!is.null(pathway)){
		warning("running pathway analysis")
		pathwayout = list()
		np = names(pathway)
		for(i in 1:length(pathway)){
			warning(paste0("pathway ", np[i]))

			respath = wilcoxGMT(fitC$t[,1], pathway[[i]], simple = T)
			ii = which(respath[,4] < 0.2)
			if(length(ii) > 0){
				oo = order(respath[ii,2])
				pathwayout[[np[i]]] = respath[ii[oo],]
			}
			else{
				warning(paste0("min QV = ", min(respath[,4], na.rm=T), " for pathway ", np[i]))
			}
		}

		fitC$pathway = pathwayout
	}
	fitC$fit = fit
	fitC$mod = mod
	#show(colSums(fitC$q.value<0.2))
	fitC
}


resid=function(dat, lab, useMean=T){
	if (is.null(dim(lab))){
		mod=model.matrix(~1+lab);
	}
	else{
		mod=lab
	}
	ne <- nonEstimable(mod)
	if (!is.null(ne)){
		cat("Coefficients not estimable:", paste(ne, collapse = " "),
			"\n")
		mod=mod[, -match(ne, colnames(mod))]
	}

	n=dim(dat)[2]
	Id=diag(n)
	resid=dat %*% (Id - mod %*% solve(t(mod) %*% mod) %*%
				   	t(mod))
	colnames(resid)=colnames(dat)
	if (useMean){
		resid=sweep(resid,1,apply(dat,1,mean), "+")
	}

	resid
}

#' multi elastic net regression
#'
#' @param x data (feature * samples)
#' @param y label
#' @param family "binomial" or "gaussian" for y
#' @param nrep total number of elastic net regressions, default nrep=100
#' @param alpha elastic net mixing parameter, default alpha=0.9
#' @param bootstrap default bootstrap=T
#' @param seed random seed, default seed=123
#'
#' @export
runGnetMulti = function(x, y, family = NULL, nrep = 100, alpha = 0.9, bootstrap = F, seed = 123, ...) {
	message(paste("running ", nrep, " multi-glmnet"))
	set.seed(123)
	if(is.null(family)) {
		if(length(unique(y)) == 2) {
			family = "binomial"
		}
		else {
			family = "gaussian"
			y = scale(y)
		}
		message(paste("family is set to", family))
	}
	vals = double(nrep)
	features = character()
	if(family == "binomial"){
		tm = "auc"
	}
	else {
		tm = "mse"
		vals2 = double()
	}
	xin = x
	yin = y
	numrm = ceiling(length(y) / nrep)
	for(i in 1:nrep) {
		set.seed(i)
		if(bootstrap) {
			ii = sample(length(yin), length(yin) - numrm)
		}
		else {
			ii = 1:length(y)
		}
		x = xin[ii,]
		y = yin[ii]
		if(family == "binomial") {
			gres = cv.glmnet(y = y, x = x, family = family, type.measure = tm, alpha = alpha, ...)
			if(gres$name == "AUC") {
				vals[i] = max(gres$cvm)
				mi = which.max(gres$cvm)
			}
			else {
				vals[i] = min(gres$cvm)
				mi = which.min(gres$cvm)
			}
		}
		else {
			set.seed(i);
			gres = cv.glmnet(y = y, x = x, family = family, type.measure = tm, alpha = alpha, keep = T, ...)
			tmp1 = t(gres$fit.preval)
			tmp2 = model.matrix(~as.factor(gres$foldid))
			resid(tmp1, tmp2)
			preval=resid(t(gres$fit.preval), model.matrix(~as.factor(gres$foldid)))

			cc = cor(y, t(preval))
			vals[i] = max(cc)
			mi = which.min(gres$cvm)
			vals2[i] = min(gres$cvm)
		}

		features = c(features,names(getNonZeroBetas(gres, i = mi,vector = T)))
	}
	if(family == "gaussian") {
		vals = cbind(vals, vals2)
		colnames(vals) = c("cor", "mse")
	}
	ft = table(features)

	num.sel = double(ncol(xin))
	names(num.sel) = colnames(xin)
	num.sel[names(ft)] = ft

	return(list(vals = vals, features = num.sel))
}

myggscatter = function(x, y, col = NULL, col.name = "color", xlab = "", ylab = "", label = NULL, line = T, num.label = 20, label.select=NULL, ...){
	num.label = min(num.label, length(x))
	df = data.frame(x = x, y = y)
	colnames(df) = c("x","y")
	if(!is.null(label.select) && length(label.select) > 100) {
		label.select = sample(label.select,100)
		warning("More than 100 labels were subsampled to 100")
	}
	if(!is.null(col)) {
		df$color = col
		colnames(df)[3] = col.name
	}
	if(!is.null(label)) {
		df$label = label
		colnames(df)[ncol(df)] = "name"
	}
	if(!is.null(label) && is.null(label.select)) {

		set.seed(0); svdres = rsvd(cbind(scale(jitter(x)), scale(jitter(y))))

		Ur = apply(abs(svdres$u), 2, rank)
		Urm = apply(Ur, 1, max)

		label.select = label[which(Urm >= sort(Urm, T)[num.label])]
		#
	}
	# show(label.select)

	if(is.null(label) && is.null(col)) {
		p = ggscatter(df, x = "x", y = "y", xlab = xlab, ylab = ylab, fill = "grey",...)
	}

	else if(is.null(label)&!is.null(col)) {
		p = ggscatter(df, x = "x", y = "y", color = col.name, xlab = xlab, ylab = ylab, ...)
	}
	else if(!is.null(label)&is.null(col)) {

		p = ggscatter(df, x = "x", y = "y", label = "name", label.select = label.select,
					  xlab = xlab, ylab = ylab, repel = T, color = "grey", font.label = c(12, "plain", "black"), ...)
	}
	else if(!is.null(label)&!is.null(col)) {
		p = ggscatter(df, x = "x", y = "y", label = "name", label.select = label.select,
					  color = col.name, xlab = xlab, ylab = ylab, repel = T,  ...)
	}
	if(line) {
		p = p + geom_abline(slope = 1)
	}
	if (!is.null(col)) {
		if(min(col) < 0){
			p = p + gradient_color(c("blue", "white", "red"))
		}}
	p
}

#' get prediction from multi elastic net regression
#'
#' @param data data (feature * samples)
#' @param y label
#' @param family "binomial" or "gaussian" for y
#' @param nrep total number of elastic net regressions, default nrep=100
#' @param alpha elastic net mixing parameter, default alpha=0.9
#' @param seed random seed, default seed=123
#'
#' @export
getPrediction = function (data, y, nrep = 100, seed = 123, family, alpha=0.9, ...) {
	set.seed(seed)
	multi_y = data.frame("sample" = colnames(data))
	coef_y = data.frame("feature" = rownames(data))
	for (i in 1:100) {
		print(paste0("working on: ", i, "th cv.glmnet..."))
		gres = cv.glmnet(y = y, x = t(data), family = "gaussian", type.measure = "mse", alpha = alpha, keep = T, ...)
		multi_y = cbind(multi_y, gres$fit.preval[, gres$cvm == min(gres$cvm)])
		coef_y = cbind(coef_y, gres$glmnet.fit$beta[, gres$cvm == min(gres$cvm)])
	}
	prediction = apply(multi_y[,-1], 1, mean)
	coef = apply(coef_y[,-1], 1, mean)

	return(list("allPredictions" = multi_y, "allCoef" = coef_y,"predictionRes" = prediction, "averagedCoef" = coef))
}

predictScatter = function(real, predict, ylim = NULL, xlim = NULL, title = NULL, xlab = NULL, ylab=NULL,...) {
	data = data.frame("predict" = predict, "real" = real)
	ylim = max(abs(predict))
	xlim = max(abs(real))
	p = ggscatter(data, "real", "predict",
				  add = "reg.line", conf.int = TRUE, # Add confidence interval
				  cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
				  cor.coeff.args = list(method = "pearson", label.x = -xlim, label.y = 0.6*ylim, label.sep = "\n", size=3),
				  title = title,
				  xlab = xlab,
				  ylab = ylab,...)
	p
}

getTopEachColumn = function(mat, top = 10) {
	topIndex = double()
	for(i in 1:ncol(mat)){
		cutoff = sort(mat[, i], T)[top+1]
		topIndex = c(topIndex, which(mat[, i] > cutoff))
	}
	topIndex = unique(topIndex)
	mat[topIndex,]
}

mypheatmap = function(data, rr = NULL, ...) {
	if (is.null(rr)) rr = max(abs(data), na.rm = T)
	pheatmap(data, breaks = seq(-rr,rr,length.out = 101),color = colorRampPalette(rev(brewer.pal(n = 7, name = "PRGn")))(100), ...)
}
