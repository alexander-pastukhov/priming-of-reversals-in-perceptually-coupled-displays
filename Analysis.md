Priming of reversals and stability in perceptually coupled multistable
displays
================
Alexander (Sasha) Pastukhov
20 March 2018

This document contains the entire analysis code to generate descriptive
statistics and figures for the manuscript.

# Preprocessing

First, we merge files of individual observers into a single
table

``` r
reports <- tibble(filename= list.files(path='Data', pattern= 'csv$')) %>%
  rowwise() %>%
  do(read_csv2(file.path('Data', .$filename))) %>%
  
  # converting Python True/False to R boolean
  mutate(Yellow= Yellow == 'True') %>%
  
  # making observer a factor
  ungroup() %>%
  mutate(Observer = as.factor(Observer))
```

Next, we create additional variables that merge `left`/`right` triggers
and reports into a `single` trigger and report.

``` r
reports <- reports %>%
  # exogenous trigger
  mutate(Trigger = fct_recode(Switch, unperturbed= 'neither', 
                                      'single trigger'= 'left', 
                                      'single trigger'= 'right', 
                                      'double trigger'= 'both'), 
         Trigger = fct_relevel(Trigger, 'unperturbed', 'single trigger', 'double trigger')) %>%
  # responses
  mutate(ResponseMerged = fct_recode(Response, 'both stable'= 'neither', 
                              'single switch'= 'left', 
                              'single switch'= 'right', 
                              'both switched'= 'both'), 
         ResponseMerged = fct_relevel(ResponseMerged, 'both stable', 'single switch', 'both switched'))
```

We also recode `Gap`from numeric values into a factor into a new
`Layout` variable for plotting.

``` r
reports <- reports %>%
  mutate(Layout = fct_recode(as.character(Gap), overlap = '-0.15', touching= '0', gap= '0.15'), 
         Layout = fct_relevel(Layout, 'overlap', 'touching', 'gap'))
```

# Probability of the report given the exogenous trigger type, the spatial layout, and timing of the switch

This generates *Figure X* and the associated *Table X*.

First, we plot group averages grouped by the *response type*. However,
the proportion is computed relative to the *trigger type* and,
therefore, values with each plot do not add up to 1. Rather, values for
the same *trigger type* across all plots add up to 1.

``` r
observer_average_report <- reports %>%
    # counting each type of report for all available combinations of exogenous trigger and layouts  
    dplyr::group_by(Observer, Layout, Trigger, ResponseMerged) %>%
    dplyr::summarize(response_count= n()) %>%
    dplyr::ungroup() %>%
  
    # filling out missing combinations with zero report counts
    tidyr::complete(Observer, Layout, Trigger, ResponseMerged, fill= list(response_count = 0)) %>%

    # computing proportion of given report types for all combinations of exogenous trigger and layouts
    # first, at the level of individual observers
    dplyr::group_by(Observer, Layout, Trigger) %>%
    dplyr::mutate(report_proportion= response_count/sum(response_count)) %>%
    dplyr::ungroup()

group_average_report <- observer_average_report %>%
    # then, at the group level
    dplyr::group_by(Layout, Trigger, ResponseMerged) %>%
    dplyr::summarize(Preport= mean(report_proportion),
                     SERRreport= sd(report_proportion)/sqrt(n()-1))
    
ggplot(data= group_average_report,
       aes(x= Layout, y= Preport, ymin= Preport-SERRreport, ymax= Preport+SERRreport,
           color= Trigger, fill= Trigger, linetype= Trigger,
           shape= Trigger, group= Trigger))+
  geom_line() +
  geom_point(size=3) +
  geom_errorbar(width=0.1, linetype= 'solid') +
  xlab('Layout') +
  ylab('P(report | layout, exogenous trigger)')+
  # theme(legend.position='none') +
  facet_grid(. ~ ResponseMerged)+
  scale_shape_manual(values=c(21, 22, 23))
```

![](Analysis_files/figure-gfm/Report%20probability%20plot-1.png)<!-- -->

``` r
# ggsave('Trigger type and layout.pdf', width= 20, height= 12, units = 'cm', useDingbats= FALSE)
```

Next we analyze how the spatial layout and the type of the exogenous
trigger (single or double) alters the probability of **each** report.
Please note that we are **excluding** unperturbed trials from the
analysis as our focus is on single/double trigger
difference.

``` r
effect_on_reports <- function(observer_reports, response_type, matching_trigger){
  observer_averages <- observer_reports %>%
      # counting each type of report for all available combinations of exogenous trigger and layouts  
      dplyr::group_by(Observer, Layout, Gap, Trigger, ResponseMerged, SwitchTime) %>%
      dplyr::summarize(response_count= n()) %>%
      dplyr::ungroup() %>%
    
      # filling out missing combinations with zero report counts
      tidyr::complete(Observer, Layout, Gap, Trigger, ResponseMerged, SwitchTime, fill= list(response_count = 0)) %>%
  
      # computing proportion of given report types for all combinations of exogenous trigger and layouts
      # first, at the level of individual observers
      dplyr::group_by(Observer, Layout, Gap, Trigger, SwitchTime) %>%
      dplyr::mutate(report_proportion= response_count/sum(response_count), 
                    report_proportion= ifelse(is.na(report_proportion), 0, report_proportion), 
                    report_logit= car::logit(report_proportion, adjust= 0.025)) %>%
      dplyr::ungroup()  
  
  selected_reports <- observer_averages %>%
      dplyr::filter((ResponseMerged==response_type) & (Trigger!='unperturbed')) %>%
      mutate(Trigger= fct_drop(Trigger))

  if (!is.null(matching_trigger)){
    selected_reports <- selected_reports %>%
      mutate(Trigger= fct_relevel(Trigger, matching_trigger))
  }
  
  # performing multi-level linear mixed modelling
  `null model` <- lme4::lmer(report_logit ~ 1 + (1|Observer), data= selected_reports, REML= FALSE)
  `+ trigger time` <- update(`null model`, .~. + SwitchTime)
  `+ layout` <- update(`+ trigger time`, .~. + Gap)
  `+ trigger type` <- update(`+ layout`, .~. + Trigger)
  `+ type x layout` <- update(`+ trigger type`, .~. + Trigger*Gap)
  lmANOVA <- anova(`null model`, `+ trigger time`, `+ layout`, `+ trigger type`, `+ type x layout`)
  
  # computing effect size for each model
  effect_size= t(sapply(list(`null model`, `+ trigger time`, `+ layout`, `+ trigger type`, `+ type x layout`), r.squaredGLMM))
  colnames(effect_size) <- c('Rm2', 'Rc2')
  lmANOVA <- cbind(data.frame(model= rownames(lmANOVA)), lmANOVA, effect_size)

  # performing Bayesian ANOVA 
   bfANOVA <- extractBF(selected_reports %>%
    ungroup() %>%
    mutate(SwitchTime= as.factor(SwitchTime), 
           iObserver= as.factor(as.numeric(Observer))) %>%
    anovaBF(report_logit ~ SwitchTime + Layout + Trigger +  iObserver, data= ., whichRandom = 'iObserver'))
   
  tested_LMM_Models <- c('SwitchTime + iObserver', 'Layout + SwitchTime + iObserver', 'Layout + Trigger + SwitchTime + iObserver', 'Layout + Trigger + Layout:Trigger + SwitchTime + iObserver')   
  
  bfANOVA$model <- rownames(bfANOVA)
  bfTable<- rbind(data.frame(bf= c(1), error= 1), bfANOVA %>% filter(model %in% tested_LMM_Models) %>% select(bf, error)) %>%
    # turning comparison into "with the previous model", as for LMM
    mutate(bf_prev= bf/lag(bf))

  list('ANOVA'= cbind(data.frame(Response= rep(response_type, nrow(lmANOVA))), lmANOVA, bfTable), 
       'Layout model'= `+ layout`)
}
```

First, we used linear mixed models for each type of the response and
apply multiple comparisons correction, before reporting
    them.

``` r
both_stable_results <- effect_on_reports(reports, 'both stable', NULL)
```

    ## Warning: 'r.squaredGLMM' now calculates a revised statistic. See the help
    ## page.

    ## Warning: data coerced from tibble to data frame

``` r
single_switch_results <- effect_on_reports(reports, 'single switch', 'single trigger')
```

    ## Warning: data coerced from tibble to data frame

``` r
both_switched_results <- effect_on_reports(reports, 'both switched', 'single trigger')
```

    ## Warning: data coerced from tibble to data frame

``` r
all_response_results <- rbind(both_stable_results[['ANOVA']], single_switch_results[['ANOVA']], both_switched_results[['ANOVA']])
all_response_results$p_adj <- p.adjust(all_response_results$`Pr(>Chisq)`)
```

## Both spheres remained stable

``` r
all_response_results %>%
  filter(Response == 'both stable') %>%
  select(-Response) %>%
  kable(., col.names = c('model', 'Df', 'AIC', 'BIC', 'logLik', 'deviance', '$\\chi^2$', '$\\chi^2$ df', '$p$', '$R_m^2$', '$R_c^2$', 'BF', 'BF error', 'BF(previous)', '$p_{adj}$'), 
        digits = c(0, 0, 1, 1, 1, 1, 1, 0, 4, 2, 2, 1, 2, 1, 4))
```

| model            | Df |    AIC |    BIC |   logLik | deviance | \(\chi^2\) | \(\chi^2\) df |  \(p\) | \(R_m^2\) | \(R_c^2\) | BF | BF error | BF(previous) | \(p_{adj}\) |
| :--------------- | -: | -----: | -----: | -------: | -------: | ---------: | ------------: | -----: | --------: | --------: | -: | -------: | -----------: | ----------: |
| null model       |  3 | 4202.2 | 4217.2 | \-2098.1 |   4196.2 |         NA |            NA |     NA |      0.00 |      0.05 |  1 |     1.00 |           NA |          NA |
| \+ trigger time  |  4 | 4203.4 | 4223.4 | \-2097.7 |   4195.4 |        0.8 |             1 | 0.3682 |      0.00 |      0.05 |  0 |     0.01 |          0.0 |      1.0000 |
| \+ layout        |  5 | 4204.8 | 4229.7 | \-2097.4 |   4194.8 |        0.7 |             1 | 0.4155 |      0.00 |      0.05 |  0 |     0.01 |          0.0 |      1.0000 |
| \+ trigger type  |  6 | 4194.7 | 4224.6 | \-2091.4 |   4182.7 |       12.0 |             1 | 0.0005 |      0.01 |      0.07 |  0 |     0.02 |         26.8 |      0.0052 |
| \+ type x layout |  7 | 4196.7 | 4231.6 | \-2091.4 |   4182.7 |        0.0 |             1 | 0.9966 |      0.01 |      0.07 |  0 |     0.02 |          0.0 |      1.0000 |

## Only one sphere reversed

``` r
all_response_results %>%
  filter(Response == 'single switch') %>%
  select(-Response) %>%
  kable(., col.names = c('model', 'Df', 'AIC', 'BIC', 'logLik', 'deviance', '$\\chi^2$', '$\\chi^2$ df', '$p$', '$R_m^2$', '$R_c^2$', 'BF', 'BF error', 'BF(previous)', '$p_{adj}$'), 
        digits = c(0, 0, 1, 1, 1, 1, 1, 0, 4, 2, 2, 1, 2, 1, 4))
```

| model            | Df |    AIC |    BIC |   logLik | deviance | \(\chi^2\) | \(\chi^2\) df |  \(p\) | \(R_m^2\) | \(R_c^2\) |        BF | BF error | BF(previous) | \(p_{adj}\) |
| :--------------- | -: | -----: | -----: | -------: | -------: | ---------: | ------------: | -----: | --------: | --------: | --------: | -------: | -----------: | ----------: |
| null model       |  3 | 3651.5 | 3666.5 | \-1822.8 |   3645.5 |         NA |            NA |     NA |      0.00 |      0.01 |         1 |     1.00 |           NA |          NA |
| \+ trigger time  |  4 | 3653.5 | 3673.4 | \-1822.8 |   3645.5 |        0.0 |             1 | 0.9044 |      0.00 |      0.01 |         0 |     0.01 | 0.000000e+00 |      1.0000 |
| \+ layout        |  5 | 3650.6 | 3675.5 | \-1820.3 |   3640.6 |        4.9 |             1 | 0.0267 |      0.00 |      0.01 |         0 |     0.01 | 1.000000e-01 |      0.2405 |
| \+ trigger type  |  6 | 3589.5 | 3619.4 | \-1788.7 |   3577.5 |       63.1 |             1 | 0.0000 |      0.06 |      0.07 | 377319516 |     0.02 | 3.040537e+12 |      0.0000 |
| \+ type x layout |  7 | 3591.5 | 3626.3 | \-1788.7 |   3577.5 |        0.0 |             1 | 0.8545 |      0.06 |      0.07 |  10427261 |     0.05 | 0.000000e+00 |      1.0000 |

``` r
single_switch_Gap <- summary(single_switch_results[['Layout model']])$coefficients['Gap', ]
```

For the `single switch` response the dependence between the response
probability and the distance between objects was positive (0.7168045±
0.3231647, *t*= 2.2180781). In other words, large distance increased the
likelihood of a single object switch.

## Both spheres reversed the rotation

``` r
all_response_results %>%
  filter(Response == 'both switched') %>%
  select(-Response) %>%
  kable(., col.names = c('model', 'Df', 'AIC', 'BIC', 'logLik', 'deviance', '$\\chi^2$', '$\\chi^2$ df', '$p$', '$R_m^2$', '$R_c^2$', 'BF', 'BF error', 'BF(previous)', '$p_{adj}$'), 
        digits = c(0, 0, 1, 1, 1, 1, 1, 0, 4, 2, 2, 1, 2, 1, 4))
```

| model            | Df |    AIC |    BIC |   logLik | deviance | \(\chi^2\) | \(\chi^2\) df |  \(p\) | \(R_m^2\) | \(R_c^2\) |   BF | BF error | BF(previous) | \(p_{adj}\) |
| :--------------- | -: | -----: | -----: | -------: | -------: | ---------: | ------------: | -----: | --------: | --------: | ---: | -------: | -----------: | ----------: |
| null model       |  3 | 4311.7 | 4326.6 | \-2152.8 |   4305.7 |         NA |            NA |     NA |      0.00 |      0.03 |  1.0 |     1.00 |           NA |          NA |
| \+ trigger time  |  4 | 4313.1 | 4333.0 | \-2152.5 |   4305.1 |        0.6 |             1 | 0.4265 |      0.00 |      0.03 |  0.0 |     0.01 |          0.0 |       1.000 |
| \+ layout        |  5 | 4311.6 | 4336.5 | \-2150.8 |   4301.6 |        3.5 |             1 | 0.0617 |      0.00 |      0.03 |  0.0 |     0.01 |          0.1 |       0.494 |
| \+ trigger type  |  6 | 4280.7 | 4310.6 | \-2134.4 |   4268.7 |       32.9 |             1 | 0.0000 |      0.03 |      0.06 | 89.7 |     0.02 |     836729.8 |       0.000 |
| \+ type x layout |  7 | 4282.6 | 4317.5 | \-2134.3 |   4268.6 |        0.1 |             1 | 0.8110 |      0.03 |      0.06 |  2.0 |     0.02 |          0.0 |       1.000 |

``` r
both_switched_Gap <- summary(both_switched_results[['Layout model']])$coefficients['Gap', ]
```

When both objects switched, the dependence between the response and the
between-objects distance was negative (-0.8173491± 0.4371662, *t*=
-1.8696529). In other words, large distance increased the likelihood of
a single object switch.

# Effect of the previous response, a.k.a. Priming

``` r
# adding information about the prior Response
reports_with_history <- reports %>%
  mutate(PreviousResponse_merged = lag(ResponseMerged), 
         PreviousResponse = lag(Response), 
         PreviousLayout = lag(Layout)) %>%
  
  # we will ignore first trial in each block, as they lack history
  filter(!is.na(PreviousResponse))
  
observer_report_given_history <- reports_with_history %>%
  dplyr::group_by(Observer, Layout,  PreviousResponse_merged, ResponseMerged) %>%
  dplyr::summarize(response_count= n()) %>%
  dplyr::ungroup() %>%  

  # filling out missing combinations with zero report counts
  tidyr::complete(Observer, Layout, PreviousResponse_merged, ResponseMerged, fill= list(response_count = 0)) %>% 
  
  # computing proportion of given report types for all combinations of prior Response and layouts
  dplyr::group_by(Observer, Layout, PreviousResponse_merged) %>%
  dplyr::mutate(report_proportion= response_count/sum(response_count)) %>%
  dplyr::ungroup()

group_report_given_history <- observer_report_given_history %>%
    dplyr::group_by(Observer, Layout, PreviousResponse_merged) %>%
    dplyr::mutate(report_proportion= response_count/sum(response_count)) %>%
    # then, at the group level
    dplyr::group_by(Layout, PreviousResponse_merged, ResponseMerged) %>%
    dplyr::summarize(Preport= mean(report_proportion),
                     SERRreport= sd(report_proportion)/sqrt(n()-1)) %>%
    dplyr::rename('Previous Response' = PreviousResponse_merged)

ggplot(data= group_report_given_history,
       aes(x= Layout, y= Preport, ymin= Preport-SERRreport, ymax= Preport+SERRreport,
           color= `Previous Response`, fill= `Previous Response`, linetype= `Previous Response`,
           shape= `Previous Response`, group= `Previous Response`))+
  geom_line() +
  geom_point(size=3) +
  geom_errorbar(width=0.1, linetype= 'solid') +
  xlab('Layout') +
  ylab('P(report | layout, prior Response)')+
  # theme(legend.position='none') +
  facet_grid(. ~ ResponseMerged)+
  scale_shape_manual(values=c(21, 22, 23))
```

![](Analysis_files/figure-gfm/Number%20of%20prior%20switches-1.png)<!-- -->

``` r
# ggsave('Priming.pdf', width= 20, height= 12, units = 'cm', useDingbats= FALSE)
```

``` r
priming_on_reports <- function(observer_reports, response_type){
  observer_averages <- observer_reports %>%
      # counting each type of report for all available combinations of exogenous trigger and layouts  
      dplyr::group_by(Observer, Gap, PreviousResponse_merged, ResponseMerged) %>%
      dplyr::summarize(response_count= n()) %>%
      dplyr::ungroup() %>%
    
      # filling out missing combinations with zero report counts
      tidyr::complete(Observer, Gap, PreviousResponse_merged, ResponseMerged, fill= list(response_count = 0)) %>%
  
      # computing proportion of given report types for all combinations of exogenous trigger and layouts
      # first, at the level of individual observers
      dplyr::group_by(Observer, Gap, PreviousResponse_merged) %>%
      dplyr::mutate(report_proportion= response_count/sum(response_count), 
                    report_proportion= ifelse(is.na(report_proportion), 0, report_proportion), 
                    report_logit= car::logit(report_proportion, adjust= 0.025)) %>%
      dplyr::ungroup()  
  
  selected_reports <- observer_averages %>%
      dplyr::filter((ResponseMerged==response_type)) %>%
      mutate(PreviousResponse_merged= fct_relevel(PreviousResponse_merged, response_type))
    

  # using lmerTest to compare other responses to the target one
  sorted_reports <- selected_reports %>%
    mutate(PreviousResponse_merged= fct_relevel(PreviousResponse_merged, response_type))
  lmTestModel <- lmerTest::lmer(report_logit ~ Gap + PreviousResponse_merged + (1|Observer), data= sorted_reports)
  
  lmTest <- data.frame(summary(lmTestModel)$coefficients)
  model <- rownames(lmTest)
  model[3:4] <- levels(sorted_reports$PreviousResponse_merged)[-1]
  lmTest <- cbind(model, data.frame(Response= rep(response_type, nrow(lmTest))), lmTest)
  rownames(lmTest) <- NULL


  list('lmTest'= lmTest,
       'lmTestModel'= lmTestModel)
}

# modelling each response type
both_stable_results <- priming_on_reports(reports_with_history, 'both stable')
single_switch_results <- priming_on_reports(reports_with_history, 'single switch')
both_switched_results <- priming_on_reports(reports_with_history, 'both switched')

# adjusting p-values for multiple comparisons 
all_priming_results <- rbind(both_stable_results[['lmTest']], single_switch_results[['lmTest']], both_switched_results[['lmTest']])
all_priming_results$p_adj <- p.adjust(all_priming_results$`Pr...t..`, method = 'holm')
```

## Priming for `both stable` response

``` r
all_priming_results %>%
  filter(Response == 'both stable') %>%
  select(-Response) %>%
  kable(., col.names = c('Model', 'Estimate', 'Std. Error', 'df', '$t$', '$p$', '$p_{adj}$'), 
        digits = c(0, 4, 2, 1, 1, 4, 4))
```

| Model         | Estimate | Std. Error |   df | \(t\) |  \(p\) | \(p_{adj}\) |
| :------------ | -------: | ---------: | ---: | ----: | -----: | ----------: |
| (Intercept)   |   0.3763 |       0.30 | 13.6 |   1.3 | 0.2292 |      0.4188 |
| Gap           | \-0.6940 |       0.55 | 93.0 | \-1.3 | 0.2094 |      0.4188 |
| single switch | \-0.8014 |       0.16 | 93.0 | \-4.9 | 0.0000 |      0.0001 |
| both switched | \-0.5029 |       0.16 | 93.0 | \-3.1 | 0.0030 |      0.0148 |

## Priming for `single switch` response

``` r
all_priming_results %>%
  filter(Response == 'single switch') %>%
  select(-Response) %>%
  kable(., col.names = c('Model', 'Estimate', 'Std. Error', 'df', '$t$', '$p$', '$p_{adj}$'), 
        digits = c(0, 4, 2, 1, 1, 4, 4))
```

| Model         | Estimate | Std. Error |   df | \(t\) |  \(p\) | \(p_{adj}\) |
| :------------ | -------: | ---------: | ---: | ----: | -----: | ----------: |
| (Intercept)   | \-0.5776 |       0.20 | 16.3 | \-2.9 | 0.0114 |      0.0456 |
| Gap           |   1.8932 |       0.50 | 93.0 |   3.8 | 0.0003 |      0.0023 |
| both stable   | \-0.6601 |       0.15 | 93.0 | \-4.4 | 0.0000 |      0.0003 |
| both switched | \-0.5680 |       0.15 | 93.0 | \-3.8 | 0.0003 |      0.0023 |

## Priming for `both switched` response

``` r
all_priming_results %>%
  filter(Response == 'both switched') %>%
  select(-Response) %>%
  kable(. , col.names = c('Model', 'Estimate', 'Std. Error', 'df', '$t$', '$p$', '$p_{adj}$'), 
        digits = c(0, 4, 2, 1, 1, 4, 4))
```

| Model         | Estimate | Std. Error |   df | \(t\) |  \(p\) | \(p_{adj}\) |
| :------------ | -------: | ---------: | ---: | ----: | -----: | ----------: |
| (Intercept)   | \-1.1193 |       0.23 | 13.9 | \-4.8 | 0.0003 |      0.0023 |
| Gap           | \-1.7708 |       0.45 | 93.0 | \-3.9 | 0.0002 |      0.0017 |
| both stable   | \-0.4289 |       0.14 | 93.0 | \-3.2 | 0.0021 |      0.0127 |
| single switch | \-0.2932 |       0.14 | 93.0 | \-2.2 | 0.0332 |      0.0995 |

# Location-specificity of priming for single switch only

Here, we are looking only at the `single switch` responses that were
preceded also by `single switch` responses. Specifically, we are
interested the layout would influence whether the second reponse is more
likely to appear on the same side as the previous one.

``` r
single_reports <- reports_with_history %>%
  filter(ResponseMerged == 'single switch' & PreviousResponse %in% c('left', 'right')) %>%
  mutate(ResponseSide= ifelse(Response== PreviousResponse, 'same', 'opposite')) %>%
  group_by(Observer, Layout) %>%
  summarize(Psame_side= mean(ResponseSide == 'same')) %>%
  ungroup() %>%
  tidyr::complete(Observer, Layout, fill= list(Psame_side = 0))

# plotting 
single_reports %>%
  mutate(Layout= fct_relevel(Layout, 'overlap', 'touching', 'gap')) %>%
ggplot(aes(x= Layout, y= Psame_side, group= Layout, color= Observer))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.2))+
  theme(legend.position = 'none') + 
  ylab('P(same side)')
```

![](Analysis_files/figure-gfm/Side-specificity%20for%20single%20switch%20priming-1.png)<!-- -->

``` r
# ggsave('Priming-Side.pdf', width= 16, height= 16, units = 'cm', useDingbats= FALSE)

# statistical comparison of gap condition versus two other layouts
single_reports <- single_reports %>%
  ungroup() %>%
  mutate(Layout= fct_relevel(Layout, 'gap'), 
         Plogit= car::logit(Psame_side, adjust = 0.025))

summary(lmerTest::lmer(Plogit ~ Layout + (1|Observer), data= single_reports))$coefficients %>%
  kable(., col.names = c('Estimate', 'Std. Error', 'df', '$t$', '$p$'), 
        digits = c(4, 2, 1, 1, 4))
```

|                | Estimate | Std. Error |   df | \(t\) |  \(p\) |
| -------------- | -------: | ---------: | ---: | ----: | -----: |
| (Intercept)    |   0.9314 |       0.43 | 29.4 |   2.2 | 0.0377 |
| Layoutoverlap  | \-1.0995 |       0.53 | 22.0 | \-2.1 | 0.0481 |
| Layouttouching | \-1.2748 |       0.53 | 22.0 | \-2.4 | 0.0239 |
