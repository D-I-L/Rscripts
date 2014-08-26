
cmdLineRun <- function() {
	args<-commandArgs(TRUE)
	if(length(args) < 6){
   	 	cat("Error incorrect number of args","\n",sep="")
    		q()
	}else{
   	 	for(i in 1:length(args)){
        		eval(parse(text=args[[i]]))
  		}
	}
	library("snpStats")
	run(data_dir, chromosome, dataset, output_dir, error_file, focalSNP, window_size, args)
}

getVar <- function(json_data, name, var, asChar) {
	if (length( json_data[[name]] ) != 0) 
		if(asChar)
			var <- sapply(json_data[[name]], as.character)
		else
                	var <- sapply(json_data[[name]], as.integer)
        else {
		if(asChar)
			var <- as.character(var)
		else
                	var <- as.integer(var)
	}
	var
}

run <- function(data_dir, chromosome, dataset, output_dir, error_file, focalSNP, window_size=1000000, args=0) {

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
 	output_dir  = getVar(json_data, 'output_dir', output_dir, TRUE)
 	error_file  = getVar(json_data, 'error_file', error_file, TRUE)
 	focalSNP    = getVar(json_data, 'focalSNP', focalSNP, TRUE)

	prefix = "genotypes_chr"
	data_file = paste(data_dir,prefix,chromosome,"_",dataset,".RData",sep="")
	if(!file.exists(data_file)){
		file = paste(output_dir,"/",dataset,'.txt',sep="")
		write(paste("Dataset file ",file, "does not exist",sep=""),file)
		q("no",0,FALSE)
		stop(paste("Dataset file ",file, "does not exist",sep=""))
	}
	load(paste(data_dir,prefix,chromosome,"_",dataset,".RData",sep=""))

	pos=snps$support[focalSNP,]$position

	if(is.na(pos)){
		file = paste(output_dir,"/",dataset,'.txt',sep="")
		write(paste("Marker ",focalSNP," was not found in this dataset",sep=""),file)
		q("no", 0,FALSE)
		stop(paste("Marker ",focalSNP," was not found in this dataset",sep=""))
	}

	if(length(grep('targetSNP',args))>0){
		print("Marker vs Marker")
		snp.ids <-c(focalSNP,targetSNP)
	} else {
		print("Region based")
		fSNP.pos <-snps$support[focalSNP,]$position
	
		if(is.na(fSNP.pos)){
			file = paste(output_dir,"/",dataset,'.txt',sep="")
			write(paste("Marker ",focalSNP," was not found in this dataset",sep=""),file);
			q("no", 0,FALSE)
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
		if(length(grep('targetSNP',args))>0){
			if(!targetSNP %in% colnames(snp_data_subset)){
				file = paste(output_dir,"/",dataset,'_',targetSNP,'.txt',sep="")
				write(paste("Marker ",targetSNP," did not pass QC",sep=""),file)
				stop(paste("Marker ",targetSNP," did not pass QC",sep=""))
			}
			dataset=paste('pairwise',focalSNP,targetSNP,dataset,sep="_")
		}
	
		ld_data <-ld( snp_data_subset, snps$genotypes[,focalSNP], stats=c("D.prime","R.squared") )
	
		ld_data_formatted <-data.frame( rep(focalSNP,length(ld_data$D.prime)), rownames(ld_data$D.prime), ld_data )
		colnames(ld_data_formatted) = c("marker1", "marker2", names(ld_data))
	
		#remove the instance where marker is compared with self
		ld_data_formatted <-ld_data_formatted[ focalSNP!=rownames(ld_data_formatted), ]
		ld_data_formatted <-ld_data_formatted[ ld_data_formatted$D.prime>=0, ]

		#colnames(ld_data_formatted)[ colnames(ld_data_formatted)=="D.prime" ]   <-"dprime"
		#colnames(ld_data_formatted)[ colnames(ld_data_formatted)=="R.squared" ] <-"rsq2"
		
		output_file = paste(output_dir,"/",dataset,'.ld',sep="")
		write.table(ld_data_formatted, file=output_file, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
	} else {
		file = paste(output_dir,"/",dataset,'.txt',sep="")
		write(paste("Marker ",focalSNP," did not pass QC",sep=""),file)
	}
        otable(ld_data_formatted, tab = "", tr = "", cs = "</td><td>")
	done()
}

#cmdLineRun()

