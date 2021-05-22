# descr -------------------------------------------------------------------
project. <- 'LoA'
out.     <- 'result.i.csv, fits.i'
desc.    <- { 
  'result.i / a table with all possible time-to-loss results; 
   fits / a list of survival-fits 
   result
   (including time to diagnosis)'
}
# last change
# 2019-10-01
# 2020-10-18 T25FW
# require -----------------------------------------------------------------

require(ReIns)
source("../.DATA other/scores.lists.R")
# scores. <- read_rds('DATA derived/scores.rds')
demo.   <- .dd.FA('demo') %>% filter(study == 'FACOMS')

# select parameters and severity groups -----------------------------------

# .dd.FA('scafi_cat') %>% 
#   filter(paramcd == 'w25.i', study == 'FACOMS') %>% 
#   ungroup %>% select(reason) %>% table(exclude = F)

w25 <- .dd.FA('scafi') %>% 
  filter( paramcd == 'w25.i') %>% filter(study == 'FACOMS') %>% 
  mutate(aval = as.numeric(unable)) %>% 
  select(-unable)

fds <- .dd.FA('steps') %>% 
  filter( study == 'FACOMS') %>% 
  select( study, sjid, avisitn, fds.act) %>% 
  rename(aval = fds.act) %>% 
  mutate(paramcd = 'fds') %>%
  left_join ( .ds.FACOMS('fsatax') %>% select(study, sjid, avisitn, adt)) %>% 
  filter(!is.na(aval))

dt. <- bind_rows(
  .dd.FA('fars'),
  fds,
  w25
  ) %>% 
  .add.time(tm = 'agedur', keepadt = T) %>% 
  .add.demo() %>% 
  left_join(demo. %>% select(sjid, diag)) %>% 
  filter( paramcd %in% c('fds','w25.i', .l.mFARS)) %>%
  # filter( paramcd %in% c('fds')) %>%
  # filter(!(paramcd %in% c('fane1','fane6','fane5','fane4','fane3b'))) %>%
  # mutate( time = dur) %>% 
  mutate(time = age) %>% 
  mutate( aval = floor(aval)) %>% 
  droplevels() %>%
  group_by( sjid, paramcd)

table(dt.$paramcd)

# scores for diagnosis parameter ------------------------------------------
# 2020 October
# symp to diag only makes sense with duration
# otherwise you get age of diagnosis (!) -> tried now

dt.sympdiag <- dt. %>%
  filter( paramcd == 'fane7') %>% # could be any parameter (I used fane7 previously)
  mutate( paramcd  = 'diagnosis' ) %>%
  group_by(sjid) %>%
  filter( dur     == min(dur)) %>%
  ungroup %>%
  # mutate( time     = diag-symp ) %>% # otherwise you get ages in the result!
  mutate( time     = diag ) %>% # yes!
  mutate( age      = diag ) %>%
  mutate( time     = ifelse(is.na(time), 0, time)) %>% # unknown diag needs a zero
  mutate( dur      =  0 ) %>%
  mutate( avisitn  = -1 ) %>%
  mutate( aval     =  1 ) # all are diagnosed by now!

# dublicate those with diag ( without will be left censored / truncated )

dt.sympdiag.0 <- dt.sympdiag %>%
  filter(!is.na(diag)) %>%
  mutate(adt     = adt - 1) %>% # adt is used to calcluate visit count later, removes one day
  mutate(aval    = 0)

dt.sympdiag %<>% 
  bind_rows(dt.sympdiag.0)

dt.sympdiag %<>% # group_by(sjid) %>% filter(symp<15) %>% 
  arrange(sjid) %>% 
  filter(!(time < 0 & aval == 0)) # removes the above added double lines for diag before symp

rm(dt.sympdiag.0)

dt.sympdiag %<>% 
  mutate(time = ifelse(time <  0, 0    , time)) %>% # set negative times to 0 (or they get filtertered from tte below)
  mutate(time = ifelse(time == 0, 0.001, time))     # make diag == symp visible in plots

# add to scores. ----------------------------------------------------------

dt. %<>%
  bind_rows(dt.sympdiag)

.ct(dt.sympdiag)
rm(dt.sympdiag)

# bind_rows(scores.sympdiag, scores.sympdiag.0) %>% filter(sjid == 1) %>% arrange(adt)

# analysis ----------------------------------------------------------------

dt. %<>% # two visits with the same date
  filter(!(sjid == 54   & avisitn == 4 & paramcd == 'w25.i')) %>% 
  filter(!(sjid == 4362 & avisitn == 6 & paramcd == 'fds')) %>% 
  filter(!is.na(aval))

# dt. %>% filter(sjid == 4900, paramcd == 'w25.i')

dt. %<>%
  group_by( sjid, paramcd) %>%
  arrange ( sjid, paramcd, time) %>% 
  # filter(sjid == 118, paramcd == 'fds')
  # filter(paramcd == 'diagnosis') %>% 
  # filter(symp<15) %>% 
  mutate  ( viscount = rank(unique(adt)))

# dt. %<>% 
#   filter(age>8)

# dt. %>% filter(paramcd == 'diagnosis') %>% 
#   filter(sjid > 200)

params   <- levels(as.factor(dt.$paramcd)) 

tte.i    <- as_tibble(data.frame())
result.i <- as_tibble(data.frame())

# multiple fits (this takes time) -----------------------------------------
# loop through parameters in scores ----

fits.i  <- list()

# for debugg
# par <- 'fane7'; thr <- 1 ; bx <-  "<15y"
# par <- 'diagnosis'; thr <- 1 ; bx <-  "<15y"

dt. %<>% 
  mutate( sever = factor('all'))
  # mutate( sever = sev.o)
  # mutate( sever = cut( symp, c( 0,    14, 24, 75), include.lowest = T, labels = c(         '<15y', '15-24y', '>24y')))

dt. %<>% 
  mutate( time = age)



for(par in levels(factor(dt.$paramcd))) { 
  
  # print(paste("---------- parameter:", par))
  # print(paste("---", par, "---"))
  
  dt <- dt. %>% filter(paramcd == par)
  
  # loop through parameters levels ----  
  for(thr in as.numeric(levels(factor(dt$aval)))) {  
    
    # print(paste("---", thr, "----------"))

    # build tte for interval censoring
    tte <- dt %>% #filter(!(sjid %in% c(1,10,103,104,105))) %>% 
      arrange(sjid, viscount) %>%
      group_by(sjid) %>% 
      mutate(left_censor = ifelse(viscount == min(viscount) & aval>=thr, 1,0)) %>% 
      mutate(left_censor = max(left_censor)) %>% 
      # filter(left_censor<1) %>%
      mutate(event  = ifelse(aval >= thr, 1, 0)) %>%
      mutate(pevent = max(event)) %>%
      group_by(sjid) %>%
      mutate(ev.vis  =  ifelse(event == 1, viscount
                               , 100)) %>%
      mutate(ev.vis  =  min(ev.vis, na.rm = T)) %>%
      filter(viscount <= ev.vis) %>% 
      group_by(sjid, event) %>%
      mutate(
        slct = ifelse(event == 0, max(time), min(time)),
        tpnt = ifelse(event == 0, "t", "t2")
      ) %>% 
      group_by(sjid) %>%
      select(study, sjid, symp, sever, left_censor, pevent, slct, tpnt) %>%
      # select(sjid, Q7_AGE, enrolled_stage, pevent, slct, tpnt) %>%
      unique() %>%
      spread(tpnt, slct)
    
    # add potentially missing columns t or t2
    {if(!("t"  %in% names( tte ))) tte$t  <- -Inf}
    {if(!("t2" %in% names( tte ))) tte$t2 <- NA}
    
    # censoring / currently only Turnbull is used. Only Diagnosis has event = 1 ( event at time )
    tte %<>% #filter(sjid == 280) %>% 
      mutate(
        event = ifelse( left_censor == 1, 2   , 3     ), # when enroll.non-amb. -> event == 2 (left censored) otherwise 3 (interval)
        event = ifelse( pevent      == 0, 0   , event ), # if no event          -> event == 0 (right censored)
        event = ifelse( t == t2         , 1   , event ), # put this here for diagnosis - but I don't think I ever use it below 2019-10-01
        event = ifelse( is.na(t2)       , 0   , event ), # this is to catch right cesored becoming Inf with the line before 2019-10-01
        t2    = ifelse( is.na(t2)       , Inf , t2    ), # necessary for some calls
        t     = ifelse( is.na(t)        , -Inf, t )      # t needs to be -Inf for interval to work (not NA). 
      ) %>%
      mutate( censor = ifelse(event == 3        , 0     , 1      )) %>% 
      mutate( censor = ifelse(event == 1        , 0     , censor )) %>% 
      mutate( censor = ifelse(is.infinite(t)    , 1     , censor )) # censor column only for Turnbull call
    
    tte %<>% mutate(paramcd = par, thr = thr) %>% 
      mutate( grpx = sever) %>%
      droplevels %>%
      filter()
    
    labels = levels(tte$grpx)
    
    tte %<>%
      # filter( is.finite(t)) %>%                           # left censored lines will be removed before the model 
      mutate( t = ifelse(is.infinite(t), 0, t)) %>%         # ReIns Turnbull censroring
      mutate( L = t, R = t2) %>% ungroup                    # %>% select(symp, t, t2, censor)
    
    # fit
    
    fits   <- list() ; i = 1
    stable <- tbl_df(data.frame())

    # Run Models ---- 

    for (bx in labels) { 
      
      f_n <- paste(par,thr,bx,sep=".")
      
      # print(paste("--- group:", bx))

      # ONE dataset ----

      tte.bx <- tte %>% 
        filter(grpx == bx) %>% 
        filter(t    <= t2)
      
      # 2020-02-11
      # this next line ensures that it's only tried to run tte if events occur (3) or (1)? 
      # when there is a 1, there also has to be a zero, so I added this condition
      if( (3 %in% tte.bx$event) | ((0 %in% tte.bx$event) & (1 %in% tte.bx$event)) ){ 
      # if( (3 %in% tte.bx$event) | (1 %in% tte.bx$event) ){ 
        # not all left or right censored or no lines: length(unique(tte.x$event))!=1 & nrow(tte.x)>0
        # replaced with "at least one event"
        
        # left censored model goes first
        print(paste("--- ", par, " --- ", thr ," --- ", bx, " --- running model, left censored"))
        
        # m.fit      <- Turnbull( L = tte.bx$t, R = tte.bx$t2, censored = tte.bx$censor, x = seq(0,1), trunclower = -Inf, truncupper = Inf)  
        m.fit      <- Turnbull( L = tte.bx$L, R = tte.bx$R, censored = tte.bx$censor, x = seq(0,1), trunclower = -Inf, truncupper = Inf)  
        # fits[i]    <- m.fit["fit"]
        f_nx        <- paste(f_n,"lcn",sep=".")
        fits.i[f_nx] <- m.fit["fit"]
        i = i + 1
        result  <-  t(summary(m.fit$fit)$table) %>% tbl_df() %>% mutate(strata = bx) 
        result$q25 <- unname(quantile(m.fit$fit, conf.int = F)[1])
        result$q75 <- unname(quantile(m.fit$fit, conf.int = F)[3])
        result %<>% mutate(censoring = "lcn")
        stable  <-  bind_rows(stable, result)
        rm(result)
        
        # remove left censored observations / but keep them for diagnosis # not so sure anymore
        # if( par != 'diagnosis' ){
        tte.bx %<>% 
          filter(L != 0)
        # }
        
        # right censored model second
        print(paste("--- ", par, " --- ", thr ," --- ", bx, " --- running model, only observed events"))
        
        m.fit  <- Turnbull( L = tte.bx$t, R = tte.bx$t2, censored = tte.bx$censor, x = seq(0,1), trunclower = -Inf, truncupper = Inf)  
        # fits[i]  <- m.fit["fit"]
        f_nx        <- paste(f_n,"ooe",sep=".")
        fits.i[f_nx] <- m.fit["fit"]
        i = i + 1
        result   <-  t(summary(m.fit$fit)$table) %>% tbl_df() %>% mutate(strata = bx) 
        result$q25 <- unname(quantile(m.fit$fit, conf.int = F)[1])
        result$q75 <- unname(quantile(m.fit$fit, conf.int = F)[3])
        result  %<>% mutate(censoring = "ooe")
        stable   <-  bind_rows(stable, result)

        rm(result, m.fit)
      }

      else {
        
        print(paste("--- ", par, " --- ", thr ," --- ", bx, " --- empty frame"))
        
        result = tbl_df( data.frame( a = c(NA), b = c(NA), c = c(NA), d = c(NA), e = c(NA), f = c(NA),
                                     g = c(NA), h = c(NA) ,i = c(NA))) %>% mutate(strata = bx)
        names(result) <-  c("records", "n.max", "n.start", "events", "*rmean", "*se(rmean)","median", "0.95LCL", "0.95UCL","strata")
        result %<>% mutate(censoring = "lcn")
        stable %<>% bind_rows(result)
        rm(result)
      }
    }
    
    rm(tte.bx, labels)
    
    result.i  %<>% bind_rows(stable %>% mutate(paramcd = par, thr = thr))
    
    rm(stable)
  }
} 

rm ( par, thr, bx, i, dt, dt., tte.i, tte, fits)
saveRDS ( result.i, "DATA derived/result.i.age.all.rds")
.wt     ( result.i, "DATA derived/result.i.age.all.csv", sep=',' )
saveRDS ( fits.i  , "DATA derived/fits.i.age.all.rds"  )

# -------------------------------------------------------------------------


# dm ----------------------------------------------------------------------
# output options
# .wt(table1, paste('code/', project., '.' , out., '.', gsub(":| ",".", Sys.time()),".csv", sep=""), sep = ',')

# read_pptx("templates/template.pptx") %>%
#   add_slide(layout = '21', master = 'CR') %>%
#   ph_with_text(type = "title", str = 'decline rates by gaa1') %>%
#   ph_with_vg(print(last_plot(),newpage = F), type = "body", index = 1) %>%
#   print(target = paste(project.,'/', out., '.', gsub(":| ",".", Sys.time()),".pptx", sep=""))


rm(desc., f_n, f_nx, params)
