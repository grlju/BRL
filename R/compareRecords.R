#' Creation of Comparison Data 
#'
#' Create comparison vectors for all pairs of records coming from 
#' two datafiles to be linked.
#'
#' @param df1,df2 two datasets to be linked, of class \code{data.frame}, with rows representing records and columns
#' 				representing fields.  Without loss of generality, 
#' 				\code{df1} is assumed to have no less records than \code{df2}.
#' @param flds	a vector indicating the fields to be used in the linkage.  Either a \code{character} vector, in which case  
#' 				all entries need to be names of columns of \code{df1} and \code{df2}, or a \code{numeric} vector
#' 				indicating the columns in \code{df1} and \code{df2} to be used in the linkage.  If provided as a 
#' 				\code{numeric} vector it is assumed that the columns of \code{df1} and \code{df2} are organized such that 
#' 				it makes sense to compare the columns
#' 				\code{df1[,flds]} and \code{df2[,flds]} in that order.  
#' @param flds1,flds2	vectors indicating the fields of \code{df1} and \code{df2} to be used in the linkage.  
#' 				Either \code{character} vectors, in which case  
#' 				all entries need to be names of columns of \code{df1} and \code{df2}, respectively, or \code{numeric} vectors
#' 				indicating the columns in \code{df1} and \code{df2} to be used in the linkage.  It is assumed that
#'				it makes sense to compare the columns
#' 				\code{df1[,flds1]} and \code{df2[,flds2]} in that order.  These arguments are ignored if \code{flds} is specified.
#'				If none of \code{flds,flds1,flds2} are specified, the columns with the same names in \code{df1} and \code{df2} 
#'				are compared, if any.  
#' @param types	a vector of characters indicating the comparison type per comparison field.  The options
#'				are: \code{"str"} for comparisons based on the string distance normalized to \eqn{[0,1]}, with \eqn{0}  
#'				indicating no disagreement and \eqn{1} indicating maximum disagreement;  
#'				\code{"bi"} for binary comparisons (agreement/disagreement); \code{"nu"} for numeric comparisons computed as 
#'				the absolute difference. 
#' 				The default is \code{"str"}.  Fields compared with the \code{"str"} option are first transformed to \code{character}
#'				class.  Factors with different levels compared using the \code{"bi"} option are transformed to factors with the union 
#' 				of the levels.  Fields compared with the \code{"nu"} option need to be of class \code{numeric}.
#' @param breaks	break points for the comparisons to obtain levels of disagreement.  
#'				It can be a list of length equal to the number of comparison fields, containing one numeric vector with the break 
#'				points for each comparison field, where entries corresponding to comparison type \code{"bi"} are ignored.  
#'				It can also be a named list of length two with elements 'str' and 'nu' 
#'				containing numeric vectors with the break 
#'				points for all string-based and numeric comparisons, respectively.  
#'				Finally, it can be a numeric vector with the break points for all comparison fields of type \code{"str"} and \code{"nu"},
#'				which might be meaningful only if all the non-binary comparisons are of a single type, either \code{"str"} or \code{"nu"}.  
#'				For comparisons based on the normalized string distance, a vector of length \eqn{L} of break 
#'				points for the interval \eqn{[0,1]} leads to \eqn{L+1} levels of disagreement.  Similarly, for comparisons based on the absolute 
#'				difference, the break points are for the interval \eqn{[0,\infty)}.  
#'				The default is \code{breaks=c(0,.25,.5)}, which might be meaningful only for comparisons of type \code{"str"}.
#' @param method method for calculating string distances. See \link[stringdist]{stringdistmatrix} for available methods.
#'        The default is \code{method = "lv"} to maintain consistency with previous versions of BRL. This is different from the 
#'        default method use in stringdist.
#' @param ... use to pass additional arguments to \link[stringdist]{stringdistmatrix} to control the string distance calculation. 
#'  
#' @return a list containing: 
#' \describe{
#'   \item{\code{comparisons}}{
#' 						matrix with \code{n1*n2} rows, where the comparison pattern for record pair \eqn{(i,j)}
#'						appears in row \code{(j-1)*n1+i}, for \eqn{i} in \eqn{{1,\dots,n1}}, and \eqn{j} in \eqn{{1,\dots,n2}}.  
#'						A comparison field with \eqn{L+1} levels of disagreement, 
#'						is represented by \eqn{L+1} columns of TRUE/FALSE indicators.  Missing comparisons are coded as FALSE,
#'						which is justified under an assumption of ignorability of the missing comparisons, see Sadinle (2017).
#'						}
#'   \item{\code{n1,n2}}{the datafile sizes, \code{n1 = nrow(df1)} and \code{n2 = nrow(df2)}.}
#'   \item{\code{nDisagLevs}}{a vector containing the number of levels of
#'					 	disagreement per comparison field.}
#'   \item{\code{compFields}}{a data frame containing the names of the fields in the datafiles used in the comparisons 
#'					 	and the types of comparison.}
#' }
#' 
#' @references Mauricio Sadinle (2017). Bayesian Estimation of Bipartite Matchings for Record Linkage. \emph{Journal of the
#' American Statistical Association} 112(518), 600-612. [\href{https://doi.org/10.1080/01621459.2016.1148612}{Published}] [\href{https://arxiv.org/abs/1601.06630}{arXiv}]
#' 
#' @examples
#' data(twoFiles)
#' 
#' myCompData <- compareRecords(df1, df2, 
#'                              flds=c("gname", "fname", "age", "occup"),
#'                              types=c("str","str","bi","bi"), 
#'                              breaks=c(0,.25,.5))
#' 
#' ## same as 
#' myCompData <- compareRecords(df1, df2, types=c("str","str","bi","bi"))
#' 
#' 
#' ## let's transform 'occup' to numeric to illustrate how to obtain numeric comparisons 
#' df1$occup <- as.numeric(df1$occup)
#' df2$occup <- as.numeric(df2$occup)
#' 
#' ## using different break points for 'str' and 'nu' comparisons 
#' myCompData1 <- compareRecords(df1, df2, 
#'                               flds=c("gname", "fname", "age", "occup"),
#'                               types=c("str","str","bi","nu"), 
#'                               breaks=list(str=c(0,.25,.5), nu=0:3))
#' 
#' ## using different break points for each comparison field
#' myCompData2 <- compareRecords(df1, df2, 
#'                               flds=c("gname", "fname", "age", "occup"),
#'                               types=c("str","str","bi","nu"), 
#'                               breaks=list(c(0,.25,.5), c(0,.2,.4,.6), NULL, 0:3))

compareRecords <- function(df1, df2, flds=NULL, flds1=NULL, flds2=NULL, types=NULL, breaks=c(0,.25,.5), method = "lv", ...){
	
	warn <- FALSE
	
	# control the input
	
	if(!is.data.frame(df1)) stop("'df1' is not a data frame")
	if(!is.data.frame(df2)) stop("'df2' is not a data frame")
	
	if(!is.null(flds)){ # user provides 'flds' fields in 'df1' and 'df2' for linkage 
	  if (!is.character(flds) && !is.numeric(flds))
			stop("'flds' should be of class character or numeric")
		
		if(!is.null(flds1)) warning("Argument 'flds' was provided, 'flds1' is ignored")
		if(!is.null(flds2)) warning("Argument 'flds' was provided, 'flds2' is ignored")
		
		flds1 <- flds2 <- flds
		
		if(is.numeric(flds)){
			fldsCheck1 <- flds %in% seq_len(ncol(df1))
			if(!all(fldsCheck1)) stop("Some numbers in 'flds' are out of the column range for 'df1'")
			fldsCheck2 <- flds %in% seq_len(ncol(df2))
			if(!all(fldsCheck2)) stop("Some numbers in 'flds' are out of the column range for 'df2'")
			warn <- TRUE
		}
		if(is.character(flds)){
			fldsCheck1 <- flds %in% colnames(df1)
			if(!all(fldsCheck1)) stop("Some names in 'flds' are not columns of 'df1'")
			fldsCheck2 <- flds %in% colnames(df2)
			if(!all(fldsCheck2)) stop("Some names in 'flds' are not columns of 'df2'")
			flds1 <- match(flds1, colnames(df1))
			flds2 <- match(flds2, colnames(df2))
		}
	}else{
		if(is.null(flds1) && !is.null(flds2)) stop("Argument 'flds1' also needs to be specified")
		if(is.null(flds2) && !is.null(flds1)) stop("Argument 'flds2' also needs to be specified")
		if(is.null(flds1) && is.null(flds2)){
			flds1 <- flds2 <- intersect(colnames(df1), colnames(df2))
			if(length(flds1)==0)
				stop("'df1' and 'df2' do not have fields with the same names")
			flds1 <- match(flds1, colnames(df1))
			flds2 <- match(flds2, colnames(df2))
		}else{ # user provides 'flds1' and 'flds2' fields in 'df1' and 'df2' for linkage 
		  if (!is.character(flds1) && !is.numeric(flds1))
				stop("'flds1' should be of class character or numeric")
		  if (!is.character(flds2) && !is.numeric(flds2))
				stop("'flds2' should be of class character or numeric")
			if(is.numeric(flds1)){
				flds1Check <- flds1 %in% seq_len(ncol(df1))
				if(!all(flds1Check)) stop("Some numbers in 'flds1' are out of the column range for 'df1'")
			}
			if(is.numeric(flds2)){
				flds2Check <- flds2 %in% seq_len(ncol(df2))
				if(!all(flds2Check)) stop("Some numbers in 'flds2' are out of the column range for 'df2'")
			}
			if(is.character(flds1)){
				flds1Check <- flds1 %in% colnames(df1)
				if(!all(flds1Check)) stop("Some names in 'flds1' are not columns of 'df1'")
				flds1 <- match(flds1, colnames(df1))
			}
			if(is.character(flds2)){
				flds2Check <- flds2 %in% colnames(df2)
				if(!all(flds2Check)) stop("Some names in 'flds2' are not columns of 'df2'")
				flds2 <- match(flds2, colnames(df2))
			}
			if(length(flds1)!=length(flds2)) stop("'flds1' and 'flds2' should have the same length")
		}
	}
	
	n1 <- nrow(df1)
	n2 <- nrow(df2)
	
	if(n2 > n1) stop("'df2' has more records than 'df1'")

	F <- length(flds1)
	
	if(is.null(types)) types <- rep("str",F)
	
	if(F != length(types)) 
		stop("Length of 'types' does not equal number of fields used for comparison")
	
	if((!is.character(types)) || any(!(types %in% c("str","bi","nu")))) 
		stop("'types' should be a character vector of 'str', 'bi', and/or 'nu' indicating string-based, 
		binary, or numeric comparisons")
	
	c1 <- sapply(df1[,flds1], class)
	c2 <- sapply(df2[,flds2], class)
	
	if(any(types=="nu") && any(c1[types=="nu", drop=FALSE] != "numeric"))
		stop(paste("Numeric comparison requested for columns '", paste(colnames(df1)[flds1[types=="nu"]],collapse="' '"), 
			"' in 'df1', but at least one of these is not numeric", sep=""))
	
	if(any(types=="nu") && any(c2[types=="nu", drop=FALSE] != "numeric"))
		stop(paste("Numeric comparison requested for columns '", paste(colnames(df2)[flds2[types=="nu"]],collapse="' '"), 
			"' in 'df2', but at least one of these is not numeric", sep=""))
	
	if(!is.list(breaks) && !is.numeric(breaks))
	  stop("'breaks' should be a list or a numeric vector")
	
	if(is.list(breaks)){ # if user specifies list with breaks 
		if(length(breaks) == F){ # breaks for each comparison field
			for(fld in 1:F){ # check the breaks provided for each field
				if(types[fld]=="str"){
					if((!is.numeric(breaks[[fld]])) | any(breaks[[fld]] < 0 | breaks[[fld]] > 1)) 
						stop(paste0("'types[",fld,"]' specified as 'str', so 'breaks[[",fld,"]]' should be a vector of numbers in [0,1]"))
					breaks[[fld]] <- unique(c(-Inf,breaks[[fld]],Inf))
				}
				if(types[fld]=="nu"){
					if((!is.numeric(breaks[[fld]])) | any(breaks[[fld]] < 0)) 
						stop(paste0("'types[",fld,"]' specified as 'nu', so 'breaks[[",fld,"]]' should be a vector of non-negative numbers"))
					breaks[[fld]] <- unique(c(-Inf,breaks[[fld]],Inf))
				}
			}
		}else{ # if length of list is not F, then list needs to be named with elements 'str' and 'nu'
			if((length(breaks) != 2) || (!all(names(breaks) %in% c("str","nu"))))
				stop("When 'breaks' is specified as a list, it should either have length equal to the number of comparison fields, or be a named list with elements 'str' and 'nu'")
			if((!is.numeric(breaks$str)) | any(breaks$str < 0 | breaks$str > 1)) 
				stop("'breaks$str should be a vector of numbers in [0,1]")
			breaks$str <- unique(c(-Inf,breaks$str,Inf))
			if((!is.numeric(breaks$nu)) | any(breaks$nu < 0)) 
				stop("'breaks$nu should be a vector of non-negative numbers")
			breaks$nu <- unique(c(-Inf,breaks$nu,Inf))
			
			breaks1 <- list()
			breaks1[1:F] <- list(NULL)
			breaks1[types=="nu"] <- list(breaks$nu)
			breaks1[types=="str"] <- list(breaks$str)
			breaks <- breaks1
		}
	}
	
	if(is.numeric(breaks)){ # if user specifies only one vector with breaks for all non-binary comparison fields
		if( any(breaks < 0) ) 
			stop("When 'breaks' is specified as a numeric vector, it should contain non-negative numbers")
		if(all(types %in% c("str","bi")) && any(breaks > 1))
			stop("When 'breaks' is specified as a numeric vector, it should contain numbers in [0,1] if all the non-binary comparisons in 'types' are 'str'")
		
		breaks <- unique(c(-Inf,breaks,Inf))
		breaks1 <- list()
		breaks1[1:F] <- list(breaks)
		breaks <- breaks1
	}
	
	pairInds1 <- rep(1:n1, n2) 
	pairInds2 <- rep(1:n2, each=n1)
	
	comparisons <- list()
	
	for(fld in 1:F){
		if(types[fld]=="bi"){
		# computing agreement levels for binary comparisons
			if(c1[fld] != c2[fld])
				warning(paste("Columns '", colnames(df1)[flds1[fld]], "' in 'df1' and '" , 
					colnames(df2)[flds2[fld]], "' in 'df2' are of different classes", sep=""))
		  if (c1[fld] == "factor" && c2[fld] == "factor"){
				levs1 <- sort(levels(df1[,flds1[fld]]))
				levs2 <- sort(levels(df2[,flds2[fld]]))
				if(!identical(levs1,levs2)){
					warning(paste("Columns '", colnames(df1)[flds1[fld]], "' in 'df1' and '" , 
						colnames(df2)[flds2[fld]], "' in 'df2' are factors with different levels", sep=""))
					unilevs <- union(levs1,levs2)
					df1[,flds1[fld]] <- factor(df1[,flds1[fld]], levels = unilevs)
					df2[,flds2[fld]] <- factor(df2[,flds2[fld]], levels = unilevs)
				}
			}
			same <- df1[pairInds1,flds1[fld]] == df2[pairInds2,flds2[fld]]
			AgrLev <- 1*same
			AgrLev[!same] <- 2
			comparisons[[fld]] <- as.factor(AgrLev)
		}
		if(types[fld]=="nu"){
			absDiff <- abs(df1[pairInds1,flds1[fld]] - df2[pairInds2,flds2[fld]])
			AgrLev <- cut(absDiff, breaks=breaks[[fld]], labels=seq_len(length(breaks[[fld]])-1), include.lowest = TRUE)
			comparisons[[fld]] <- as.factor(AgrLev)			
		}
		if(types[fld]=="str"){
		# computing agreement levels for string based comparisons
			df1[,flds1[fld]] <- as.character(df1[,flds1[fld]])
			df2[,flds2[fld]] <- as.character(df2[,flds2[fld]])
			
			# stringsimmatrix in the stringdist package returns the matrix of string distances normalized between 0 and 1
			# the normalization is 0 for complete agreement and 1 for complete disagreement, so we reverse it. 
			lvd <- 1-as.numeric(stringdist::stringsimmatrix(df1[,flds1[fld]], df2[,flds2[fld]], method = method))
			
			AgrLev <- cut(lvd, breaks=breaks[[fld]], labels=seq_len(length(breaks[[fld]])-1), include.lowest = TRUE)
			comparisons[[fld]] <- as.factor(AgrLev)
		}
	}
	
	nDisagLevs <- sapply(comparisons, FUN=nlevels)
		
	# Compute the binary indicators from the agreement levels
	comparisons <- do.call(cbind, lapply(seq_len(F), function(fld){
	  sapply(seq_len(nDisagLevs[fld]), function(lvl){
	    comparisons[[fld]] == lvl
	  })
	}))
	
	# replacing NAs by FALSE or zeroes is justified by ignorability and CI assumptions (see Sadinle 2017)
	comparisons[is.na(comparisons)] <- FALSE 
	
	df1Fields <- colnames(df1)[flds1]
	df2Fields <- colnames(df2)[flds2]
	compFields <- data.frame(file1=df1Fields, file2=df2Fields, types=types)
	
	if(any(df1Fields != df2Fields) || warn){
		warning(paste("The fields '", paste(df1Fields,collapse="' '"), 
			"' in 'df1' are being compared with the fields '",
			paste(df2Fields,collapse="' '"), "' in 'df2'",  
			sep=""))
	}
	
	out <- list(comparisons=comparisons, n1=n1, n2=n2, nDisagLevs=nDisagLevs, 
				compFields=compFields)
	
	return(out)
}
