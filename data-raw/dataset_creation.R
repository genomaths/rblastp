## Copyright (C) 2022 Robersy Sanchez <https://genomaths.com/>
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

# ========================================================================== #
#
# ======== Script used to generate the datasets used in the examples ======= #
#
# ========================================================================== #

library(Biostrings)

url1 <- paste0("https://github.com/genomaths/seqalignments/raw/master/",
               "Pyrococcus/pep_seq_all/",
               "Pyrococcus_abyssi_ge5_gca_000195935.ASM19593v2.pep.all.fa.gz")

url2 <- paste0("https://github.com/genomaths/seqalignments/raw/master/",
               "Pyrococcus/pep_seq_all/",
               "Pyrococcus_furiosus_dsm_3638_gca_000007305.ASM730v1.pep.all.fa.gz")

## The destination files are in the temporal 'tmp' local folder
dir.create("/tmp/pyroc")
outfile1 <- "/tmp/pyroc/p_abiss.fa"
outfile2 <- "/tmp/pyroc/p_furiosus.fa"

download.file(url = url1, destfile = outfile1, method = "wget", quiet=TRUE)
download.file(url = url2, destfile = outfile2, method = "wget", quiet=TRUE)

file <- "/tmp/pyroc/p_furiosus.fa"
p_furiosus <- readAAStringSet(filepath = file, format = "fasta")
file <- "/tmp/pyroc/p_abiss.fa"
p_abiss <- readAAStringSet(filepath = file, format = "fasta")

unlink(x = "/tmp/pyroc", recursive = TRUE)

usethis::use_data(p_furiosus, p_abiss, overwrite = TRUE)








