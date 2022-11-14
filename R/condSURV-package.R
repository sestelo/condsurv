

#' Bladder Cancer Recurrences.
#'
#' bladderCS is a data frame with 8 variables and 85 observations.  Data on
#' recurrences of bladder cancer, used by many people to demonstrate
#' methodology for recurrent event modelling.
#'
#'
#' @name bladderCS
#' @docType data
#' @format  A data frame with 85 observations on the following 8 variables.
#' Below a brief description is given for some of these variables.
#' \describe{
#' \item{t1}{Time to first recurrence/censoring, whichever occurs first.}
#' \item{e1}{Recurrence/censoring indicator (first recurrence=1) for the first
#' time (t1).} \item{t2}{Time to second recurrence/censoring, whichever occurs
#' first.} \item{e2}{Recurrence/censoring indicator (second recurrence=1) for
#' the second time (t2)} \item{t3}{Time to recurrence/censoring, whichever
#' occurs first.} \item{e3}{Recurrence/censoring indicator (third recurrence=1)
#' for the third time (t3)} \item{t4}{Time to fourth recurrence/censoring,
#' whichever occurs first.} \item{e4}{Recurrence/censoring indicator (fourth
#' recurrence=1) for the fourth time (t4)}
#' }
#' @references
#'
#' Byar, D. (1980) Veterans administration study of chemoprophylaxis for
#' recurrent stage i bladder tumors: Comparisons of placebo, pyridoxine and
#' topical thiotepa. Bladder Tumors and Other Topics in Urological Oncology,
#' 18:363--370.
#' @keywords datasets
#' @examples
#'
#'   data(bladderCS)
#'   head(bladderCS)
#'
NULL





#' Chemotherapy for Stage B/C colon cancer.
#'
#' These are data from one of the first successful trials of adjuvant
#' chemotherapy for colon cancer. Levamisole is a low-toxicity compound
#' previously used to treat worm infestations in animals; 5-FU is a moderately
#' toxic (as these things go) chemotherapy agent.
#'
#'
#' @name colonCS
#' @docType data
#' @format  A data frame with 929 observations on the following 15 variables.
#' Below a brief description is given for some of these variables. \describe{
#' \item{time1}{Time to recurrence/censoring/death, whichever occurs first.}
#' \item{event1}{Recurrence/censoring indicator (recurrence=1, alive=0).}
#' \item{Stime}{Time to censoring/death, whichever occurs first.}
#' \item{event}{Death/censoring indicator (death=1, alive=0).}
#' \item{rx}{Treatment - Obs(ervation), Lev(amisole), Lev(amisole)+5-FU.}
#' \item{sex}{Sex indicator (male=1, female=0).} \item{age}{Age in years.}
#' \item{obstruct}{Obstruction of colon by tumour.} \item{perfor}{Perforation
#' of colon.} \item{adhere}{Adherence to nearby organs.} \item{nodes}{Number of
#' lymph nodes with detectable cancer.} \item{differ}{Differentiation of tumour
#' (1=well, 2=moderate, 3=poor).} \item{extent}{Extent of local spread
#' (1=submucosa, 2=muscle, 3=serosa, 4=contiguous structures).}
#' \item{surg}{Time from surgery to registration (0=short, 1=long).}
#' \item{node4}{More than 4 positive lymph nodes.} }
#' @references JA Laurie, CG Moertel, TR Fleming, HS Wieand, JE Leigh, J Rubin,
#' GW McCormack, JB Gerstner, JE Krook and J Malliard. Surgical adjuvant
#' therapy of large-bowel carcinoma: An evaluation of levamisole and the
#' combination of levamisole and fluorouracil: The North Central Cancer
#' Treatment Group and the Mayo Clinic. Journal of Clinical Oncology,
#' 7:1447-1456, 1989.
#'
#' DY Lin. Cox regression analysis of multivariate failure time data: the
#' marginal approach. Statistics in Medicine, 13:2233-2247, 1994.
#'
#' CG Moertel, TR Fleming, JS MacDonald, DG Haller, JA Laurie, PJ Goodman, JS
#' Ungerleider, WA Emerson, DC Tormey, JH Glick, MH Veeder and JA Maillard.
#' Levamisole and fluorouracil for adjuvant therapy of resected colon
#' carcinoma. New England Journal of Medicine, 332:352-358, 1990.
#'
#' CG Moertel, TR Fleming, JS MacDonald, DG Haller, JA Laurie, CM Tangen, JS
#' Ungerleider, WA Emerson, DC Tormey, JH Glick, MH Veeder and JA Maillard.
#' Fluorouracil plus Levamisole as and effective adjuvant therapy after
#' resection of stage II colon carcinoma: a final report. Annals of Internal
#' Medicine, 122:321-326, 1991.
#' @source The study is originally described in Laurie (1989).The main report
#' is found in Moertel (1990). This data set is closest to that of the final
#' report in Moertel (1991). A version of the data with less follow-up time was
#' used in the paper by Lin (1994).
#' @keywords datasets
#' @examples
#'
#' data(colonCS)
#' head(colonCS)
#'
NULL





#' condSURV:
#' c("\\Sexpr[results=rd,stage=build]{tools:::Rd_package_title(\"#1\")}",
#' "condSURV")\Sexpr{tools:::Rd_package_title("condSURV")} %A package for
#' nonparametric estimation of the survival functions for ordered multivariate
#' failure time data.
#'
#' Newly developed methods for the estimation of the conditional survival
#' function are implemented in this package. The condSURV package implements
#' nonparametric and semiparametric estimators for these quantities. The
#' package also implements feasible estimation methods for these quantities
#' conditionally on current or past covariate measures. Other related
#' estimators are also implemented in the package. One of these estimators is
#' the Kaplan-Meier estimator typically assumed to estimate the survival
#' function. A modification of the Kaplan-Meier estimator, based on a
#' preliminary estimation (presmoothing) of the censoring probability for the
#' survival time, given the available information is also implemented.
#'
#'
#' @name condSURV-package
#' @aliases condSURV-package condSURV
#' @docType package
#' @author Luis Meira-Machado and Marta Sestelo.
#'
#' Maintainer: Marta Sestelo, sestelo@@uvigo.es
#' @references L. Meira-Machado, M. Sestelo, and A. Goncalves (2016).
#' Nonparametric estimation of the survival function for ordered multivariate
#' failure time data: a comparative study. Biometrical Journal, 58(3),
#' 623--634.
NULL





#' German Breast Cancer Study Data.
#'
#' gbcsCS is a data frame with 16 variables and 686 observations. Cancer
#' clinical trials are a rich source for examples of applications of methods
#' for the analysis of time to event. Willi Sauerbrei and Patrick Royston have
#' graciously provided us with data obtained from the German Breast Cancer
#' Study Group, which they used to illustrate methods for building prognostic
#' models (Sauerbrei and Royston, 1999). In the main study, a total of 720
#' patients with primary node positive breast cancer were recruited between
#' July 1984, and December 1989, (see Schmoor, Olschweski and Schumacher M.
#' 1996 and Schumacher et al. (1994)).
#'
#'
#' @name gbcsCS
#' @docType data
#' @format  A data frame with 686 observations on the following 16 variables.
#' Below a brief description is given for some of these variables. \describe{
#' \item{rectime}{Time to recurrence/censoring, whichever occurs
#' first.}
#' \item{censrec}{Recurrence/censoring indicator (recurrence=1,
#' alive=0).}
#' \item{survtime}{Time to censoring/death, whichever occurs
#' first.}
#' \item{censdead}{Death/censoring indicator (death=1, alive=0).}
#' \item{age}{Age in years.}
#' \item{size}{Tumour size.}
#' #' \item{deathdate}{}
#' \item{diagdateb}{}
#' \item{estrg_recp}{}
#' \item{grade}{}
#' \item{hormone}{}
#' \item{id}{}
#' \item{nodes}{}
#' \item{prog_recp}{}
#' \item{recdate}{}
#' \item{menopause}{}
#' }
#' @references Schmoor, C., Sauerbrei, W. Bastert, G., Schumacher, M. (2000).
#' Role of Isolated Locoregional Recurrence of Breast Cancer: Results of Four
#' Prospective Studies. Journal of Clinical Oncology, 18(8), 1696-1708.
#'
#' Schumacher, M., Bastert, G., Bojar, H., Hiibner, K., Olschewski, M.,
#' Sauerbrei, W., Schmoor, C., Beyerle, C., Neumann, R.L.A. and Rauschecker,
#' H.F. for the German Breast Cancer Study Group (GBSG) (1994). A randomized 2
#' x 2 trial evaluating hormonal treatment and the duration of chemotherapy in
#' node-positive breast cancer patients. Journal of Clinical Oncology, 12,
#' 2086-2093.
#'
#' Hosmer, D.W. and Lemeshow, S. and May, S. (2008). Applied Survival Analysis:
#' Regression Modeling of Time to Event Data: Second Edition, John Wiley and
#' Sons Inc., New York, NY
#' @keywords datasets
#' @examples
#'
#' data(gbcsCS)
#' head(gbcsCS)
#'
NULL



