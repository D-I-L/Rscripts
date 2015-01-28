
getVar <- function(json_data, name, var, asChar) {
	if (length( json_data[[name]] ) != 0) 
		if(asChar)
			var <- sapply(json_data[[name]], as.character)
		else
            var <- sapply(json_data[[name]], as.numeric)
        else {
		if(asChar)
			var <- as.character(var)
		else
			var <- as.numeric(var)
	}
	var
}

getHeader <- function(dprime, rsq) {
	hdr1 <- list(
    	compare = "&nbsp;", length = 10, order = 2,
    	url = "/page/Overview/display/marker_id/REPLACE/build/NCBI36",
    	name = "marker2", title = "Marker", filter_default = ""
    );
    hdr2 <- list(
    	compare = "&ge;", length = 3, order = 3, url = NULL,
    	name = "D.prime", title = "D&#39", filter_default = dprime
    );
	hdr3 <- list(
        compare = "&ge;", length = 3, order = 4, url = NULL,
        name = "R.squared", title = "r<sup>2</sup>", filter_default = rsq
    );
	hdr <-  list(hdr1, hdr2, hdr3)	
}

run <- function(chromosome, dataset, marker1, marker2=NULL, window_size=1000000, dprime=0, rsq=0.8, display="json") {

	if(length(.GlobalEnv$request.body) == 0){
		# probably GET
		json_data = {}
	} else {
		# POST
		json_data <- fromJSON(rawToChar(.GlobalEnv$request.body))
	}

	window_size = getVar(json_data, 'window_size', window_size, FALSE)
 	chromosome  = getVar(json_data, 'chromosome', chromosome, FALSE)
 	dataset     = getVar(json_data, 'dataset', dataset, TRUE)
 	marker1     = getVar(json_data, 'marker1', marker1, TRUE)
 	marker2     = getVar(json_data, 'marker2', marker2, TRUE)
 	dprime      = getVar(json_data, 'dprime', dprime, FALSE)
 	dprime      = getVar(json_data, 'D.prime', dprime, FALSE)
	rsq         = getVar(json_data, 'rsq', rsq, FALSE)
	rsq         = getVar(json_data, 'R.squared', rsq, FALSE)
	display     = getVar(json_data, 'display', display, TRUE)

	root = root <- Sys.getenv("ROOT")
	data_dir = paste(root,"RDATA",dataset,"",sep='/')

	prefix = "genotypes_chr"
	data_file = paste(data_dir,prefix,chromosome,"_",dataset,".RData",sep="")
	if(!file.exists(data_file)){
		err <- paste("Data not found: ",prefix,chromosome,"_",dataset,".RData",sep="")
	} else {
		load(data_file)
		pos=snps$support[marker1,]$position
		if(is.na(pos)){
			err <- paste("Marker ",marker1," was not found in this dataset",sep="")
		} else {
			if(length(marker2) != 0) {
				#print("Marker vs Marker")
				snp.ids <-c(marker1,marker2)
			} else {
				#print("Region based")
				fSNP.pos <-snps$support[marker1,]$position
	
				if(is.na(fSNP.pos)){
					err <- paste("Marker ",marker1," was not found in this dataset",sep="")
				}
				neighbouring.snp.ind <-which( snps$support$position > (fSNP.pos-window_size) &
						snps$support$position < (fSNP.pos+window_size) &
						!is.na(snps$support$marker_mart) )
   	 			snp.ids <-rownames(snps$support)[ neighbouring.snp.ind ]
			}

			snp_data_subset <-snps$genotypes[,snp.ids]
			#Filter based on some quality scores
			snp_stats       <-col.summary(snp_data_subset)
			snp_data_subset <-snp_data_subset[, snp_stats$MAF >= 0.01 & snp_stats$Call.rate >= 0.90 & snp_stats$z.HWE^2 < 25]

			#Check that snp makes it through quality control
			if (marker1 %in% colnames(snp_data_subset)) {
				#for pairwise between 2 markers we need a different output file recipe.
				if(length(marker2) != 0){
					if(!marker2 %in% colnames(snp_data_subset)){
						msg <- paste("Marker ",marker2," did not pass QC",sep="")
					}
					dataset=paste('pairwise',marker1,marker2,dataset,sep="_")
				}

				ld_data <-ld( snp_data_subset, snps$genotypes[,marker1], stats=c("D.prime","R.squared") )
				ld_data_formatted <- data.frame(  rownames(ld_data$D.prime), ld_data )
				colnames(ld_data_formatted) = c("marker2", names(ld_data))
	
				#remove the instance where marker is compared with self
				ld_data_formatted <- ld_data_formatted[ marker1!=rownames(ld_data_formatted), ]
				ld_data_formatted <- na.omit(ld_data_formatted)
				#format R squared before filtering
				ld_data_formatted$R.squared <- round(ld_data_formatted$R.squared, digits=2)
				ld_data_formatted <- ld_data_formatted[ ld_data_formatted$R.squared>=rsq, ]
				msg <- ld_data_formatted[ ld_data_formatted$D.prime>=dprime, ]
			} else {
				err <- paste("Marker ",marker1," did not pass QC",sep="")
			}
		}
	}
	if(length(display) != 0) {
		if(tolower(display) == 'json') {
			if(exists("err")) {
				err <- toJSON(list(error=err))
			} else {
				# make a list of lists
				lol <- do.call(Map, c(list, msg))
      			hdr <- getHeader(dprime, rsq)	
        		msg <- toJSON(list(ld=lol, header=hdr))
        	}
        } else {
        	msg = otable(msg)
        }
    }
}

