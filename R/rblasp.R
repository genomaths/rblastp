## Copyright (C) 2019 Robersy Sanchez <https://genomaths.com/>
## Author: Robersy Sanchez
##
## This program is part 'GCcube' R package (not published yet).
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## 'GCcube' is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

#' @rdname blastp
#' @title BLASTP programs to search protein subjects using a protein query
#' @description This function is a wrapping to use the 
#'     \href{https://www.ncbi.nlm.nih.gov/books/NBK279671/}{NCBI} BLASTP
#'     programs for search protein subjects using a protein query. The function
#'     applies the command line version of the Basic Local Alignment Search Tool
#'     (BLAST), as provided by the
#'     \href{https://www.ncbi.nlm.nih.gov/books/NBK279671/}{NCBI}. The reason
#'     why the function is addressed to only BLASTP is because we are interested
#'     only in the codon sequences and, in general, the sequence alignment based
#'     on DNA sequence breaks the codons. The best accuracy is obtained by
#'     traslating the codon sequence into amino acid sequence.
#'     
#' @details BLAST command line sofware must be previously installed in the local
#'     computer where this function will be called, as explained at  
#'     \href{https://www.ncbi.nlm.nih.gov/books/NBK279671/}{NCBI}. it is worthy
#'     to notice that once created the ncbi-blast database contains six files 
#'     with the following structure:       
#'     \itemize{
#'         \item{dataBaseName.phr}
#'         \item{dataBaseName.pin}
#'         \item{dataBaseName.pog}
#'         \item{ataBaseName.psd}
#'         \item{dataBaseName.psi}
#'         \item{dataBaseName.psq}
#'     }
#' @param query.seq An AAStringSet object 
#'     (\code{link[Biostrings]{XStringSet-class}}) from *Biostrings* package
#'     carrying the query amino acid sequence(s). This obejct can be created 
#'     reading the fasta file of interest with function *readDNAStringSet* from
#'     *Biostrings* package.
#' @param dtb whole path to the database of sequences. For example, it would be
#'     dtb = "/Path_to_database/dataBaseName" (do not add any extension 
#'     after the "dataBaseName", see details). The database must be 
#'     created by function 'blastdbcmd', as specified for 
#'     \href{https://www.ncbi.nlm.nih.gov/books/NBK279671/}{NCBI} blast. Default
#'     is NULL. If not provided, then function *blasp* will try to create it 
#'     using the information provided in parameters *db.fa* and *dir.fa*.
#' @param db.fa,dir.fa Name and directory of the file containing the fasta 
#'     sequence(s) to build a ncbi-blast database if dtb = NULL.
#' @param tmp A directory where to write a temporal file.  
#' @param maxTargetSeqs Integer value to pass to blastp. The maximum number of 
#'     sequences to target in the database. Default is 2, i.e., for each query
#'     sequence, the two top sequences from the database with the minimum 
#'     pairwise alignment *evalue* will be returned. 
#' @param numcode The NCBI genetic code number for translation. By default the
#'     standard genetic code is used.
#' @param num.cores,tasks Parameters for parallel computation using
#'     \code{\link[BiocParallel]{BiocParallel-package}}: the number of cores to
#'     use, i.e. at most how many child processes will be run simultaneously
#'     (see \code{\link[BiocParallel]{bplapply}} and the number of tasks per job
#'     (only for Linux OS).
#' @importFrom seqinr syncodons
#' @importFrom Biostrings writeXStringSet
#' @importFrom BiocParallel MulticoreParam SnowParam bplapply
#' @author Robersy Sanchez <https://genomaths.com>
#' @export
#' @examples
#' # ## ============================= For Linux OS =============================
#' # ## ---------- Download the amino acid sequence from github ----------------
#' # url1 <- paste0("https://github.com/genomaths/seqalignments/raw/master/",
#' #                "Pyrococcus/pep_seq_all/",
#' #                "Pyrococcus_abyssi_ge5.ASM19593v2.pep.all.fa")
#' # 
#' # url2 <- paste0("https://github.com/genomaths/seqalignments/raw/master/",
#' #                "Pyrococcus/pep_seq_all/Pyrococcus_furiosus_dsm_3638.",
#' #                "ASM730v1.pep.all.fa")
#' # 
#' # ## The destination files are in the temporal 'tmp' local folder
#' # dir.create("/tmp/pyroc")
#' # outfile1 <- "/tmp/pyroc/p_abiss.fa"
#' # outfile2 <- "/tmp/pyroc/p_furiosus.fa"
#' # 
#' # download.file(url = url1, destfile = outfile1, method = "wget", quiet=TRUE)
#' # download.file(url = url2, destfile = outfile2, method = "wget", quiet=TRUE)
#' # 
#' ## ---------------- Set the information needed for blastp ------------------
#' # db.fa <- "p_abiss.fa"
#' # dir.fa <- "/tmp/pyroc/"
#' # dir.db <- "/tmp/pyroc/"
#' # tmp <- dir.db
#' # 
#' # file <- "/tmp/pyroc/p_furiosus.fa"
#' # seqs <- readAAStringSet(filepath = file, format = "fasta")
#' # query.seq <- seqs[1:20]
#' # seq.name  <- substr(x = names(query.seq), start = 1, stop = 8)
#' # 
#' # blastp(query.seq=query.seq, dtb = NULL, seq.name = seq.name, db.fa = db.fa,
#' #       dir.fa = dir.fa, tmp = tmp)
#' # 
#' # ## ------------------------- delete 'pyroc' folder  -----------------------
#' # unlink(x = "/tmp/pyroc", recursive = TRUE)


blastp <- function(query.seq, dtb = NULL, seq.name = NULL, db.fa = NULL,
                   dir.fa = NULL, tmp = getwd(), maxTargetSeqs = 2, numcode = 1, 
                   num.cores = 1L, tasks = 0L) {
  if (missing(query.seq)) stop("Provide query sequence(s) in FASTa format") 
  if (is.null(db.fa) && is.null(dtb)) 
    stop("Provide sequence database or a fasta file to create it") 
  
  if (!is.null(dtb)) {
    db <- suppressWarnings(try(system(paste0("blastdbcmd -db ", dtb,
                                             " -info"),
                                      intern = TRUE, ignore.stderr = TRUE), 
                               silent = TRUE))
    if (inherits(db, "try-error") || length(db) == 0) 
      stop("The database provided is not valid")
  } 
  
  if (is.null(dir.fa) && is.null(dtb)) stop("Provide the fasta file directory")
  else {
    if (is.null(dtb) && !is.null(db.fa)) {
      db <- suppressWarnings(try(system(paste0("blastdbcmd -db ", dir.db,
                                               db.fa, ".prot -info"),
                                        intern = TRUE, 
                                        ignore.stderr = TRUE),
                                 silent = TRUE))
      if (!inherits(db, "try-error") && length(db) > 0)
        dtb <-paste0(dir.db, db.fa, ".prot")
    }
  }
  
  if (is.null(dtb) && !is.null(db.fa) && !is.null(dir.fa)) {
    newdb <- paste0("makeblastdb -in ", dir.fa, db.fa, " -dbtype prot ", 
                    "-parse_seqids -out ", dir.db, db.fa, ".prot ",
                    "-title ", db.fa)
    system(newdb)
    dtb <- try(system(paste0("blastdbcmd -db ", dir.db, db.fa, ".prot -info")
                      , intern = TRUE))
    if (inherits(dtb, "try-error") || length(dtb) == 0)
      stop("Database creation failed")
    dtb <-paste0(dir.db, db.fa, ".prot")
  } 
  
  # === Auxiliar function to perform blastp through OS command  ===
  blast <- function(k, query, dtb, tmp, seq.name, maxTargetSeqs) {
    sname <- seq.name[k]
    writeXStringSet(query[k], filepath = paste0(tmp, "tmp", sname, ".fasta"))
    str1 = paste0("blastp -db ", dtb, " -query ", tmp, "tmp", sname, ".fasta")
    str2 = paste0("-out ", tmp, "tmp", sname, 
                  ".txt -outfmt '6 qseqid sseqid pident ",
                  "qcovs bitscore score evalue' -max_target_seqs ", 
                  maxTargetSeqs)
    system(paste( str1, str2, sep = " " ))
    tmp1 <- try( read.delim(paste0(tmp, "tmp", sname,".txt"), 
                            header = FALSE ), silent = TRUE)
    if (!inherits(tmp1, "try-error")) {
      res <- DataFrame(tmp1)
      colnames(res) <- c("qseqid", "sseqid", "pident", "qcovs", 
                         "bitscore", "score", "evalue")
      file.remove(c(paste0(tmp, "tmp", sname,".txt"), 
                    paste0(tmp, "tmp", sname, ".fasta" ) ) )
      return(res)
    } else {
      file.remove(paste0(tmp, "tmp", sname, ".fasta" ))
      res = DataFrame(qseqid = NA, sseqid = NA, pident = NA, 
                      qcovs = NA, bitscore = NA, score = NA, evalue = NA)
    }
    return(res)
  }
  # ------------------------------------------------------------------------- #
  
  if (length(query.seq) > 1) {
    # Set parallel computation
    if (Sys.info()['sysname'] == "Linux") {
      bpparam <- MulticoreParam(workers=num.cores, tasks = tasks)
    } else bpparam <- SnowParam(workers = num.cores, type = "SOCK")
    num.seq <- 1:length(query.seq)
    res <- bplapply(num.seq, blast, query=query.seq, dtb=dtb, tmp=tmp, 
                    seq.name=seq.name, maxTargetSeqs=maxTargetSeqs,
                    BPPARAM = bpparam)
    res <- do.call(rbind, res)
  } else res <- blast(1, query=query.seq, dtb=dtb, tmp=tmp, 
                      seq.name=seq.name, maxTargetSeqs)
  
  return(res)
}