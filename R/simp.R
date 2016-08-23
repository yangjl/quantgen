#' \code{Simulate phenotype}
#'
#' simulate QTL effect and phentoype.
#'
#'
#' @param inputdf An input data.frame for fastq files. Must contains fq1, fq2, out (and/or bam).
#' If inputdf contained bam, bwa alignment will be escaped.
#' Additional columns: group (group id), sample (sample id), PL (platform, i.e. illumina),
#' LB (library id), PU (unit, i.e. unit1). These strings (or info) will pass to BWA mem through -R.
#'
#' @param h2 The full path of genome with bwa indexed reference fasta file.
#' @param gatkpwd The absolute path of GenomeAnalysisTK.jar.
#' @param picardpwd The absolute path of picard.jar.
#' @param minscore Minimum score to output, default=5, [bwa 30]. It will pass to bwa mem -T INT.
#'
#' @param markDup Mark Duplicates, default=TRUE.
#' @param addRG Add or replace Read Groups using Picard AddOrReplaceReadGroups, default=FALSE.
#' @param realignInDels Realign Indels, default=FALSE. IF TRUE, a golden indel.vcf file should be provided.
#' @param indels.vcf The full path of indels.vcf.
#' @param recalBases Recalibrate Bases, default=FALSE. IF TRUE, a golden snps.vcf file should be provided.
#' @param dbsnp.vcf The full path of dbsnp.vcf.
#' @param email Your email address that farm will email to once the jobs were done/failed.
#' @param runinfo Parameters specify the array job partition information.
#' A vector of c(FALSE, "bigmemh", "1"): 1) run or not, default=FALSE
#' 2) -p partition name, default=bigmemh and 3) --cpus, default=1.
#' It will pass to \code{set_array_job}.
#'
#' @return return a batch of shell scripts.
#'
#' @examples
#' inputdf <- data.frame(fq1="fq_1.fq", fq2="f1_2.fq", out="mysample",
#'                  group="g1", sample="s1", PL="illumina", LB="lib1", PU="unit1")
#'
#' run_GATK(inputdf,
#'          ref.fa="~/dbcenter/Ecoli/reference/Ecoli_k12_MG1655.fasta",
#'          gatkpwd="$HOME/bin/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar",
#'          picardpwd="$HOME/bin/picard-tools-2.1.1/picard.jar",
#'          markDup=TRUE,
#'          realignInDels=FALSE, indels.vcf="indels.vcf",
#'          recalBases=FALSE, dbsnp.vcf="dbsnp.vcf",
#'          email=NULL, runinfo = c(FALSE, "bigmemh", 1))
#'
#' @export
simp <- function(X, h2, alpha, NQTN, distribution, a2=0){
  n=nrow(X)
  m=ncol(X)
  #Sampling QTN
  QTN.position <- sample(m,NQTN,replace=F)
  SNPQ=as.matrix(X[,QTN.position])
  QTN.position

  #QTN effects
  if(distribution=="normal"){
    addeffect <- rnorm(NQTN,0,1)
  }else{
    addeffect=alpha^(1:NQTN)
  }
  #Simulate phenotype
  effect <- SNPQ%*%addeffect
  effectvar=var(effect)

  #Simulate Interaction
  cp=0*effect
  nint=4
  if(a2 > 0 & NQTN >= nint){
    for(i in nint:nint){
      Int.position=sample(NQTN,i,replace=F)
      cp=apply(SNPQ[,Int.position],1,prod)
    }
    print(dim(cp))
    cpvar=var(cp)
    intvar=(effectvar-a2*effectvar)/a2
    if(is.na(cp[1]))stop("something wrong in simulating interaction")

    if(cpvar!=0){
      print(c(effectvar,intvar,cpvar,var(cp),a2))
      print(dim(cp))
      cp=cp/sqrt(cpvar)
      cp=cp*sqrt(intvar)
      effectvar=effectvar+intvar
    }else{cp=0*effect}
  }

  residualvar <- (effectvar - h2*effectvar)/h2
  residual <- rnorm(n, 0, sqrt(residualvar))
  y <- effect+residual+cp
  return(list(addeffect = addeffect, y=y, add = effect, residual = residual,
              QTN.position=QTN.position, SNPQ=SNPQ,int=cp))
}
