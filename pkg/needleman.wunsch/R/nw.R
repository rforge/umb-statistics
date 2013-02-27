# Recode character to integer vectors
char2int <- function(str){
  lapply(1:length(str), function(s)as.integer(factor(strsplit(str[s],"")[[1]], levels=c("A","C","G","T")))-1L)
}

# Core function calling compiled C++ of Needleman-Wunsch
nw <- function(seq1, seq2, S, g){
  .Call( "nw", seq1, seq2, as.integer(S), as.integer(g), PACKAGE = "needleman.wunsch" )
}

# Core function calling compiled C++ of Needleman-Wunsch (two element vector S)
nws <- function(seq1, seq2, S, g){
  .Call( "nws", seq1, seq2, as.integer(S), as.integer(g), PACKAGE = "needleman.wunsch" )
}

# Core function calling compiled C++ of Needleman-Wunsch (no S)
qnw <- function(seq1, seq2){
  .Call( "qnw", seq1, seq2, PACKAGE = "needleman.wunsch" )
}

# Score satellites and/or sequences
nw.sat.quer <- function(satellites, queries, S, g){
  .Call( "nw_sat_quer", satellites, queries, as.integer(S), as.integer(g), as.integer(max(unlist(lapply(satellites,length)))), as.integer(max(unlist(lapply(queries,length)))), PACKAGE = "needleman.wunsch" )
}

# Score satellites and/or sequences (two element vector S)
nws.sat.quer <- function(satellites, queries, S, g){
  .Call( "nws_sat_quer", satellites, queries, as.integer(S), as.integer(g), as.integer(max(unlist(lapply(satellites,length)))), as.integer(max(unlist(lapply(queries,length)))), PACKAGE = "needleman.wunsch" )
}

# Score satellites and/or sequences (no S)
qnw.sat.quer <- function(satellites, queries){
  .Call( "qnw_sat_quer", satellites, queries, as.integer(max(unlist(lapply(queries,length)))), PACKAGE = "needleman.wunsch" )
}

# Semi-global score calling compiled C++ of alternative Needleman-Wunsch
nw.sg <- function(short, long, S, g){
  .Call( "nw_sg", short, long, as.integer(S), as.integer(g), PACKAGE = "needleman.wunsch" )
}

# Semi-global score calling compiled C++ of alternative Needleman-Wunsch (two element vector S)
nws.sg <- function(short, long, S, g){
  .Call( "nws_sg", short, long, as.integer(S), as.integer(g), PACKAGE = "needleman.wunsch" )
}

# Semi-global score calling compiled C++ of alternative Needleman-Wunsch (no S)
qnw.sg <- function(short, long){
  .Call( "qnw_sg", short, long, PACKAGE = "needleman.wunsch" )
}

# Semi-global score (multiple comparisons)
nw.sg.mult <- function(short, long, S, g){
  .Call( "nw_sg_mult", short, long, as.integer(S), as.integer(g), as.integer(max(unlist(lapply(short,length)))), as.integer(max(unlist(lapply(long,length)))), PACKAGE = "needleman.wunsch" )
}

# Semi-global score (multiple comparisons) (two element vector S)
nws.sg.mult <- function(short, long, S, g){
  .Call( "nws_sg_mult", short, long, as.integer(S), as.integer(g), as.integer(max(unlist(lapply(short,length)))), as.integer(max(unlist(lapply(long,length)))), PACKAGE = "needleman.wunsch" )
}

# Semi-global score (multiple comparisons) (no S)
qnw.sg.mult <- function(short, long){
  .Call( "qnw_sg_mult", short, long, as.integer(max(unlist(lapply(long,length)))), PACKAGE = "needleman.wunsch" )
}
