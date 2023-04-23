#' @importFrom XICOR "calculateXI"
#' @importFrom stats "var"

choose_p <- function(test_data,cov_e){
  data = test_data
  xi = c()
  res = c()
  d = 0
  C = 0
  n = nrow(data)
  p = ncol(data)-2

  for( i in 1:p){
    w = data[,i]
    var_e = cov_e[i,i]
    var_w = var(w)
    x = mean(w) + as.numeric( (var_w-var_e)*solve(var_w) )*(w-mean(w))
    C = calculateXI(x,  data$Y)
    xi = c(xi , C)
  }

  var.list = colnames(data)[1:p]
  names(xi) <- c( var.list )

  d = round(n/log(n),0)
  selected_Var = sort(xi,decreasing = TRUE)[1:d]
  return(selected_Var)
}
