
labeler <- function(bin_num, bin_size) {
  i = 1
  l <- c(paste("0-", bin_size, sep = ""))

  while (i < bin_num + 1 ) {
    i = i+1
    l <- append(l, paste((i-1)*bin_size, "-", i*bin_size, sep = "" ))
  }
  return(l)
}

get_bin_label <- function(x, bin_num, bin_size) {
  b <- labeler(bin_num, bin_size)
  return(b[x])
}

get_bin_label_V <- Vectorize(get_bin_label, vectorize.args = c('x','bin_num', 'bin_size'))
