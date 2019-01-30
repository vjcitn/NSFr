#' transform fixed template into program using user parameter settings
#' @param tsvname character(1)
#' @export
doSubs = function(tsvname, clustVblName="Clusters",
    numtrees=100, numthread=1, num_inf_genes = 5,
    num_top_ranked = 3) {
#
# grab template file
    templ = readLines(system.file("python/NSforestTemplate.py", package="NSFr"))
    prog = gsub("%%TSVNAME%%", tsvname, templ)
    prog = gsub("%%CLUSTVBLNAME%%", clustVblName, prog)
    prog = gsub("%%NUMTREES%%", numtrees, prog)
    prog = gsub("%%NUMTHREADS%%", numthread, prog)
    prog = gsub("%%NUM_INF_GENES%%", num_inf_genes, prog)
    prog = gsub("%%NUM_TOP_RANKED%%", num_top_ranked, prog)
    chk = grep("%%", prog)
    stopifnot(length(chk)==0)
    prog
}

# generate the TSV File from an SCE

# invoke program having written output of doSubs to a tempfile
# where you source it  -- run it in a temp folder and get
# all the csv read in and put in a list

# version 2 -- minimize the use of files

