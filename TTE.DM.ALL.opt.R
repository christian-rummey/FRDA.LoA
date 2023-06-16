# setup, data -------------------------------------------------------------

require(ReIns)

# parameter ---------------------------------------------------------------

# pars  <- c('fane2a','fds')

title  <-  c('Steps to LoA')
pars   <-  c(.l.FARS.E[c(2,3,4,5,6,7,9)])

# data --------------------------------------------------------------------

demo.   <- .dd('demo.l') %>% filter(study == 'FACOMS')

fds <- .dd('steps') %>% 
  .add.time(tm = 'age', keepadt = T) %>% 
  filter( study == 'FACOMS') %>% 
  ungroup %>% 
  select( study, sjid, adt, avisitn, fds.act) %>% 
  rename(aval = fds.act) %>% 
  mutate(paramcd = 'fds') %>%
  filter(!is.na(aval))

dt. <- bind_rows(
  # w25,
  .dd('adl'), 
  .dd('fars'),
  fds
  ) %>% 
  filter(paramcd %in% pars)%>%
  select(-fpf, -hpf) 

dt. %<>%
  # filter( aval == 6) %>% 
  left_join(.rt('../DATA other/scales.txt') %>% select(paramcd, maxscore)) %>% 
  mutate( event = 0 ) %>% 
  mutate( event = ifelse( aval >= maxscore            , 1, event )) %>% 
  mutate( event = ifelse( paramcd == 'fds' & aval == 5, 1, event )) %>% 
  .add.time() %>% 
  left_join( demo. %>% select(sjid, sev.o, symp, pm)) %>% filter(!is.na(symp)) %>% 
  mutate(dur = age-symp) %>% 
  select( study, sjid, symp, sev.o, age, dur, paramcd, aval, event )

dt. %<>% 
  droplevels() %>%
  group_by( sjid, paramcd)

dt. %<>% filter(!is.na(aval))

# correct fds/fsatax duplicate date entries -------------------------------
  # filter(!(sjid ==  231 & avisitn == 13 & adt == as.Date('2021-07-09'))) %>% 
  # filter(!(sjid ==  160 & avisitn == 10 & adt == as.Date('2020-05-11'))) %>% 
  # filter(!(sjid == 4620 & avisitn ==  4 & adt == as.Date('2021-06-14'))) %>% 
  # filter(!(sjid == 4235 & avisitn == 11 )) %>% 
  # filter(!(sjid == 4362 & avisitn ==  5 )) %>% 

dt. %<>% unique()

with(dt., table(event, paramcd, exclude = F))

rm(demo., fds)

# add visit counts  -------------------------------------------------------

dt. %<>%
  arrange ( sjid, age ) %>% 
  mutate  ( viscount = rank(unique(age)))

# groups, no time --------------------------------------------------------

dt. %<>% 
  mutate( aval  = event        ) %>%
  # mutate( grp   = factor('all')) %>% 
  mutate( grp   = sev.o        )


# multiple fits (this takes time) -----------------------------------------

tte.i    <- as_tibble(data.frame())
result.i <- as_tibble(data.frame())
fits.i   <- list()

par = 'fds'
thr = '5'

for(par in levels(factor(dt.$paramcd))) { 
  
  for(tm in c('age', 'dur') ) {
    
    if(tm == 'age') {
      dt <- dt. %>% filter(paramcd == par) %>% mutate(time = age)
    } else {
      dt <- dt. %>% filter(paramcd == par) %>% mutate(time = dur)
    }
    
    # loop through parameters levels ----
    for(thr in as.numeric(levels(factor(dt$aval)))) {  
      
      # print(paste("---", thr, "----------"))
  
      # build tte for interval censoring
      tte <- dt %>% 
        # filter((sjid %in% c(10,103,104,105))) %>% 
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
        select(study, sjid, symp, grp, left_censor, pevent, slct, tpnt) %>%
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
        mutate( grpx = grp) %>%
        droplevels %>%
        filter()
      
      labels = levels(tte$grpx)
      
      tte %<>%
        # filter( is.finite(t)) %>%                           # left censored lines will be removed before the model 
        mutate( t = ifelse(is.infinite(t), 0, t)) %>%         # ReIns Turnbull censroring
        mutate( L = t, R = t2) %>% ungroup                    # %>% select(symp, t, t2, censor)
      
      # fit
      
      fits   <- list() ; i = 1
      stable <- tibble::as_tibble(data.frame())
  
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
          result  <-  t(summary(m.fit$fit)$table) %>% tibble::as_tibble() %>% mutate(strata = bx) 
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
          result   <-  t(summary(m.fit$fit)$table) %>% tibble::as_tibble() %>% mutate(strata = bx) 
          result$q25 <- unname(quantile(m.fit$fit, conf.int = F)[1])
          result$q75 <- unname(quantile(m.fit$fit, conf.int = F)[3])
          result  %<>% mutate(time = tm, censoring = "ooe")
          stable   <-  bind_rows(stable, result)
  
          rm(result, m.fit)
        }
  
        else {
          
          print(paste("--- ", par, " --- ", thr ," --- ", bx, " --- empty frame"))
          
          result = tibble::as_tibble( data.frame( a = c(NA), b = c(NA), c = c(NA), d = c(NA), e = c(NA), f = c(NA),
                                       g = c(NA), h = c(NA) ,i = c(NA))) %>% mutate(strata = bx)
          names(result) <-  c("records", "n.max", "n.start", "events", "*rmean", "*se(rmean)","median", "0.95LCL", "0.95UCL","strata")
          result %<>% mutate(time = tm, censoring = "lcn")
          stable %<>% bind_rows(result)
          rm(result)
        }
      }
      
      rm(tte.bx, labels)
      
      result.i  %<>% bind_rows(stable %>% mutate(time = tm, paramcd = par, thr = thr))
      
      rm(stable)
      
    }
  }
} 

.wds   ( result.i, paste0('Results', '/', format(Sys.time(), "%Y-%m-%d %H%M%S"),' ', title, 'result'))
saveRDS( fits.i  , paste0('Results', '/', format(Sys.time(), "%Y-%m-%d %H%M%S"),' ', title, 'fits'))
rm ( par, thr, bx, i, dt, dt., tte.i, tte, fits)

# -------------------------------------------------------------------------