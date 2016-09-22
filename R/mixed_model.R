#' \code{Linear mixed model}
#'
#' Get the fixed effect from a linear mixed model using `nlme` package.
#'
#'
#' @param data data to fit the model. [data.frame]
#' @param model Fixed effects of the model, i.e. KRN ~ Genotype. [model]
#' @param random Random effects of the model, i.e. ~1 | Farm/Rep. [model]
#' @param trait Name of the phenotypic trait. [character]
#'
#' @return return a batch of shell scripts.
#'
#' @examples
#' blue <- mixed_model(data = df, model = KRN ~ Pedigree, random = ~1 | Farm/Rep, trait = "KRN")
#'
#' @export
mixed_model <- function(data = df, model = KRN ~ Pedigree, random = ~1 | Farm/Rep,
                        trait = "KRN") {
  #library("nlme")
  trait <- as.character(model)[2]
  data <- data[!is.na(data[, trait]), ]
  data[, trait] <- as.numeric(as.character(data[, trait]))

  lmeout1 <- lme(model, data = data, random = random)
  ped.hat1 <- lmeout1$coef$fixed
  ped.hat1[-1] <- ped.hat1[-1] + ped.hat1[1]

  fix <- as.character(model)[3]
  names(ped.hat1)[1] <- as.character(data[order(data[, fix]), fix][1])
  names(ped.hat1) <- gsub(fix, "", names(ped.hat1))
  tped <- data.frame(Genotype = names(ped.hat1), trait = ped.hat1)
  names(tped)[2] <- trait
  return(tped)
}

