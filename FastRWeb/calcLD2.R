
cmdLineRun <- function() {
	args<-commandArgs(TRUE)
	if(length(args) < 4){
   	 	cat("Error incorrect number of args","\n",sep="")
    		q()
	}else{
   	 	for(i in 1:length(args)){
        		eval(parse(text=args[[i]]))
  		}
	}
	library("snpStats")
	run(data_dir, chromosome, dataset, focalSNP, window_size, args)
}

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

run <- function(data_dir, chromosome, dataset, focalSNP, window_size=1000000, dprime=0, rsq=0.8, display="json", targetSNP=NULL, args=0) {

	if(length(.GlobalEnv$request.body) == 0){
		# probably GET
		json_data = {}
	} else {
		# POST
		json_data <- fromJSON(rawToChar(.GlobalEnv$request.body))
	}

	window_size = getVar(json_data, 'window_size', window_size, FALSE)
	data_dir    = getVar(json_data, 'data_dir', data_dir, TRUE)
 	chromosome  = getVar(json_data, 'chromosome', chromosome, FALSE)
 	dataset     = getVar(json_data, 'dataset', dataset, TRUE)
 	focalSNP    = getVar(json_data, 'focalSNP', focalSNP, TRUE)
 	display     = getVar(json_data, 'display', display, TRUE)
 	targetSNP   = getVar(json_data, 'targetSNP', targetSNP, TRUE)
 	dprime      = getVar(json_data, 'dprime', dprime, FALSE)
	rsq         = getVar(json_data, 'rsq', rsq, FALSE)

	prefix = "genotypes_chr"
	data_file = paste(data_dir,prefix,chromosome,"_",dataset,".RData",sep="")
	if(!file.exists(data_file)){
		msg <- paste("Data not found: ",prefix,chromosome,"_",dataset,".RData",sep="")
	} else {
		load(paste(data_dir,prefix,chromosome,"_",dataset,".RData",sep=""))
		pos=snps$support[focalSNP,]$position
		if(is.na(pos)){
			msg <- paste("Marker ",focalSNP," was not found in this dataset",sep="")
		} else {

			if(length(targetSNP) != 0) {
				#print("Marker vs Marker")
				snp.ids <-c(focalSNP,targetSNP)
			} else {
				#print("Region based")
				fSNP.pos <-snps$support[focalSNP,]$position
	
				if(is.na(fSNP.pos)){
					msg <- paste("Marker ",focalSNP," was not found in this dataset",sep="")
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
			if (focalSNP %in% colnames(snp_data_subset)) {
				#for pairwise between 2 markers we need a different output file recipe.
				if(length(targetSNP) != 0){
					if(!targetSNP %in% colnames(snp_data_subset)){
						msg <- paste("Marker ",targetSNP," did not pass QC",sep="")
					}
					dataset=paste('pairwise',focalSNP,targetSNP,dataset,sep="_")
				}

				ld_data <-ld( snp_data_subset, snps$genotypes[,focalSNP], stats=c("D.prime","R.squared") )
				ld_data_formatted <- data.frame( rep(focalSNP,length(ld_data$D.prime)), rownames(ld_data$D.prime), ld_data )
				colnames(ld_data_formatted) = c("marker1", "marker2", names(ld_data))
	
				#remove the instance where marker is compared with self
				ld_data_formatted <- ld_data_formatted[ focalSNP!=rownames(ld_data_formatted), ]
				ld_data_formatted <- ld_data_formatted[ ld_data_formatted$R.squared>=rsq, ]
				msg <- ld_data_formatted[ ld_data_formatted$D.prime>=dprime, ]
			} else {
				msg <- paste("Marker ",focalSNP," did not pass QC",sep="")
			}
		}
	}
	if(length(display) != 0) {
		if(tolower(display) == 'json') {
        	msg <- toJSON(msg)
        } else {
        	msg = otable(msg)
        }
    }
}

#cmdLineRun()
