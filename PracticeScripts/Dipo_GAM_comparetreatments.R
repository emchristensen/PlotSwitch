


# stuff I was playing around with -- at end of Diop_GAM_analysis.Rmd


# extract and plot smooth for treatment effect for each treatment separately

pd = plot(m1)
df1 = data.frame(n = pd[[1]]$fit, treatment = rep('CC'), numericdate = pd[[1]]$x, se = pd[[1]]$se)
df2 = data.frame(n = pd[[2]]$fit, treatment = rep('EC'), numericdate = pd[[2]]$x, se = pd[[2]]$se)
df3 = data.frame(n = pd[[3]]$fit, treatment = rep('XC'), numericdate = pd[[3]]$x, se = pd[[3]]$se)

df = rbind(df1,df2,df3)
df$lower = df$n - 2 * df$se
df$upper = df$n + 2 * df$se


# smooth of EC treatment
ggplot(df[df$treatment=='EC',],aes(x=numericdate,y=n,colour=treatment)) +
  geom_line() +
  geom_vline(xintercept=16.53) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) 

# smooth of XC treatment
ggplot(df[df$treatment=='XC',],aes(x=numericdate,y=n,colour=treatment)) +
  geom_line() +
  geom_vline(xintercept=16.53) +
  geom_ribbon(aes(ymin=lower,ymax=upper), alpha = .2)
