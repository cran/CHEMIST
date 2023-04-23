#' @importFrom XICOR "calculateXI"
#' @importFrom stats "var"

choose_p_D <- function(test_data,cov_e){
  data = test_data
  xi = c()
  res = c()
  d = 0
  C = 0
  n = nrow(data)
  p = ncol(data)-2

  index0 <- which(data$D == 0)
  index1 <- which(data$D == 1)
  alpha = length(index0)/n
  for( i in 1:p){
    w_0 = data[index0,i]
    var_e = cov_e[i,i]
    var_w_0 = var(w_0)
    #x_0 = mean(w_0) + as.numeric( (var_w_0-var_e)*solve(var_w_0) ) *(w_0-mean(w_0))

    w_1 = data[index1,i]
    var_e = cov_e[i,i]
    var_w_1 = var(w_1)
    #x_1 = mean(w_1) + as.numeric( (var_w_1-var_e)*solve(var_w_1) ) *(w_1-mean(w_1))

    x_0 = data[index0,i]
    x_1 = data[index1,i]
    C = alpha*calculateXI(x_0 , data$Y[index0] ) +
      (1-alpha)*calculateXI(x_1 , data$Y[index1])
    xi = c(xi , C )
  }

  var.list = colnames(data)[1:p]
  names(xi) <- c( var.list )

  d = round(n/log(n),0)
  selected_Var = sort(xi,decreasing = TRUE)[1:d]
  return(selected_Var)
}

