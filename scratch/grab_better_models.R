pp <- structure(list(
  foo = c(1L, 2L, 3L, 1L, 3L),
  group = c(1, 1, 1, 2, 2)),
  .Names = c("foo", "group"),
  row.names = c(2L, 3L, 1L, 4L, 5L),
  class = "data.frame"
)

pp$x2 <- rnorm(5)
pp$x3 <- as.factor(c(1,0,1,0,1))

pp_form <- as.formula(~as.factor(foo) + x2 + x3)

## good
full_model_frame <- model.frame(pp_form, data=pp)
full_model_matrix <- model.matrix(pp_form, data=pp)

factor_vars <- rep(NA, NCOL(full_model_frame))
for (ii in 1:NCOL(full_model_frame)){
  factor_vars[ii] <- is.factor(full_model_frame[,ii])
}

##good
model.matrix(pp_form, data=pp[1:3,])

##fails
model.matrix(pp_form, data=pp[4:5,])

##good?
x_levels_vec <- levels(as.factor(pp$foo))
x_levels_list1 <- list(
  foo = x_levels_vec
)
x_levels_list2 <- list(
  `as.factor(foo)` = x_levels_vec
)

model.matrix(
  pp_form,
  data=pp[4:5,],
  xlev = x_levels_list1
)

model.matrix(
  pp_form,
  data=pp[4:5,],
  xlev = x_levels_list2
)

pmf <- model.frame(
  pp_form,
  data=pp[4:5,],
  xlev = x_levels_list2
)





###

dc_mf <- attr(attr(full_model_frame, "terms"), "dataClasses")
dcf <- dc_mf=="factor"

x_levels_list <- list()
ii = 1
for (var in 1:length(dc_mf)){
  if (dc_mf[var]=="factor"){
    these_levels <- levels(
      full_model_frame[,names(dc_mf)[var]])
    x_levels_list[[ii]] <- these_levels
    names(x_levels_list)[[ii]] <- names(dc_mf)[var]
    ii <- ii+1
  }


}

model.matrix(
  pp_form,
  data = pp[4, ],
  xlev = x_levels_list
)



myfun <- function(dfm, x_levels, form){

  model.matrix(form, data = dfm, xlev = x_levels)


}


myfun(dfm=pp[4,], x_levels = x_levels_list, form = pp_form)
