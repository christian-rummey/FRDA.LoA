# descr -------------------------------------------------------------------
# require -----------------------------------------------------------------
require(survminer)

# result.i <- read_rds("DATA derived/TTE.ADL.EFACTS.age.result.txt.rds")
# fits.i   <- read_rds("DATA derived/TTE.ADL.EFACTS.age.fits.rds")

# p0 LoA ------------------------------------------------------------------

with(result.i, table(paramcd, thr, censoring))

# grid. <- c(
#   'a7.walk' , c(5),
#   'w25.iu'  , c(1),
#   'fds', c(5)    
# )

res.tmp <- result.i %>%
  filter(thr != 0) %>%
  filter(paramcd == 'fds') %>% 
  # filter(
  #     # ( paramcd == 'fane7'  & thr %in% c( 5 ) ) |
  #     # ( paramcd == 'w25.i'  & thr %in% c( 1 ) ) |
  #   ( paramcd == 'fds'    & thr %in% c( 1,2,3, 4 ) ) |
  #   ( paramcd == 'fds'    & thr %in% c( 5 ) )
  # ) %>%
  filter(censoring == 'lcn') %>%
  # filter(strata    == '<15y') %>%
  mutate(fn        =  paste(paramcd, thr,strata,censoring, sep='.')) %>%
  arrange(median) %>% 
  droplevels()

p0 <- ggsurvplot_combine(fits.i[res.tmp$fn],
                         xlim = c(0,41),
                         break.time.by = 5,
                         # legend = 'none',
                         # legend.labs = res.tmp$strata,
                         censor = F,
                         # palette = c('#2171B5', '#6BAED6','black', '#BDD7E7', 'grey' ),
                         # xlab = 'Disease Duration [years]',
                         # ylab = 'proportion retaining function',
                         # legend.title = 'Strata',
                         risk.table       = T, risk.table.fontsize = .1, tables.height = 0.3,
                         fun = 'pct'
)



vs <- c(res.tmp %>% filter(paramcd %in% pars))$median

p0$plot <- p0$plot +
  geom_hline(yintercept = 50, linetype = 'dashed')+
  geom_vline(xintercept = vs, linetype = 'dashed')

# p1 <15y -----------------------------------------------------------------

res.tmp <- result.i %>%
  filter(thr != 0) %>%
  # filter(paramcd %in% c("dis_stage", .l.FARS.E)) %>%
  # filter(paramcd %in% c("dis_stage")) %>%
  filter(
    # (paramcd == 'dis_stage' ) |
    (paramcd == 'fane2a'  & thr == 4) |
      # (paramcd == 'fane1'   & thr == 4) |
      (paramcd == 'fane2b'  & thr == 4) |
      (paramcd == 'fane3a'  & thr == 4) |
      (paramcd == 'diagnosis'  & thr == 1) |
      (paramcd == 'fane3b'  & thr == 4) |
      (paramcd == 'fane4'  & thr == 4) |
      (paramcd == 'fane5'  & thr == 4) |
      # (paramcd == 'fane6'  & thr == 3) |
      (paramcd == 'fane7'  & thr == 5) 
      # paramcd == 'fanc1' & thr %in% c(5,6,7,8)
    ) %>%
  filter(censoring == 'lcn') %>%
  filter(strata    == '<15y') %>%
  mutate(fn        =  paste(paramcd, thr,strata,censoring, sep='.')) %>%
  arrange(median) %>% 
  droplevels()

p1 <- ggsurvplot_combine(fits.i[res.tmp$fn],
                        xlim = c(0,20),
                        break.time.by = 2,
                        legend = 'none',
                        legend.labs = c('E5','E4','E3B','Diagnosis','E2B','E3A','E2A','LoA/E7'),
                        censor = F,
                        palette = c(RColorBrewer::brewer.pal(9, "Greens")[3:5],'blue',RColorBrewer::brewer.pal(9, "Greens")[6:8],"red"),
                        # palette = c(RColorBrewer::brewer.pal(9, "Greens")[3:8],"red"),
                        xlab = 'Disease Duration [years]',
                        ylab = 'Proportion of subjects with\nabilities retained',
                        legend.title = 'LoA / Stance Items',
                        risk.table       = T, risk.table.fontsize = .1, tables.height = 0.4,
                        fun = 'pct'
)



vs <- c(res.tmp %>% filter(paramcd %in% c('fane3a','fane2b','fane2a','fane7')))$median

p1$plot <- p1$plot +
  geom_hline(yintercept = 50, linetype = 'dashed')+
  geom_vline(xintercept = vs, linetype = 'dashed')

# p2 15-24y -----------------------------------------------------------------

res.tmp <- result.i %>%
  filter(thr != 0) %>%
  # filter(paramcd %in% c("dis_stage", .l.FARS.E)) %>%
  # filter(paramcd %in% c("dis_stage")) %>%
  filter(
    # (paramcd == 'dis_stage' ) |
    (paramcd == 'fane2a'  & thr == 4) |
      # (paramcd == 'fane1'   & thr == 4) |
      (paramcd == 'fane2b'  & thr == 4) |
      (paramcd == 'fane3a'  & thr == 4) |
      (paramcd == 'diagnosis'  & thr == 1) |
      (paramcd == 'fane3b'  & thr == 4) |
      (paramcd == 'fane4'  & thr == 4) |
      (paramcd == 'fane5'  & thr == 4) |
      # (paramcd == 'fane6'  & thr == 3) |
      (paramcd == 'fane7'  & thr == 5) 
    # paramcd == 'fanc1' & thr %in% c(5,6,7,8)
  ) %>%
  filter(censoring == 'lcn') %>%
  filter(strata    == '15-24y') %>%
  mutate(fn        =  paste(paramcd, thr,strata,censoring, sep='.')) %>%
  arrange(median) %>% 
  droplevels()

p2 <- ggsurvplot_combine(fits.i[res.tmp$fn],
                         xlim = c(0,20),
                         break.time.by = 2,
                         legend = 'none',
                         legend.labs = c('E5','E4','E3B','Diagnosis','E2B','E3A','E2A','LoA/E7'),
                         censor = F,
                         palette = c(RColorBrewer::brewer.pal(9, "Greens")[3:5],'blue',RColorBrewer::brewer.pal(9, "Greens")[6:8],"red"),
                         # palette = c(RColorBrewer::brewer.pal(9, "Greens")[3:8],"red"),
                         xlab = 'Disease Duration [years]',
                         ylab = 'Proportion of subjects with\nabilities retained',
                         legend.title = 'LoA / Stance Items',
                         risk.table       = T, risk.table.fontsize = .1, tables.height = 0.4,
                         fun = 'pct'
)

vs <- c(res.tmp %>% filter(paramcd %in% c('fane3a','fane2b','fane2a','fane7')))$median

p2$plot <- p2$plot +
  geom_hline(yintercept = 50, linetype = 'dashed')+
  geom_vline(xintercept = vs, linetype = 'dashed')

# p2 >24y -----------------------------------------------------------------

res.tmp <- result.i %>%
  filter(thr != 0) %>%
  # filter(paramcd %in% c("dis_stage", .l.FARS.E)) %>%
  # filter(paramcd %in% c("dis_stage")) %>%
  filter(
    # (paramcd == 'dis_stage' ) |
    (paramcd == 'fane2a'  & thr == 4) |
      # (paramcd == 'fane1'   & thr == 4) |
      (paramcd == 'fane2b'  & thr == 4) |
      (paramcd == 'fane3a'  & thr == 4) |
      (paramcd == 'diagnosis'  & thr == 1) |
      (paramcd == 'fane3b'  & thr == 4) |
      (paramcd == 'fane4'  & thr == 4) |
      (paramcd == 'fane5'  & thr == 4) |
      # (paramcd == 'fane6'  & thr == 3) |
      (paramcd == 'fane7'  & thr == 5) 
    # paramcd == 'fanc1' & thr %in% c(5,6,7,8)
  ) %>%
  filter(censoring == 'lcn') %>%
  filter(strata    == '>24y') %>%
  mutate(fn        =  paste(paramcd, thr,strata,censoring, sep='.')) %>%
  arrange(median) %>% 
  droplevels()

p3 <- ggsurvplot_combine(fits.i[res.tmp$fn],
                         xlim = c(0,20),
                         break.time.by = 2,
                         legend = 'none',
                         legend.labs = c('E5','E4','E3B','Diagnosis','E2B','E3A','E2A','LoA/E7'),
                         censor = F,
                         palette = c(RColorBrewer::brewer.pal(9, "Greens")[3:5],'blue',RColorBrewer::brewer.pal(9, "Greens")[6:8],"red"),
                         # palette = c(RColorBrewer::brewer.pal(9, "Greens")[3:8],"red"),
                         xlab = 'Disease Duration [years]',
                         ylab = 'Proportion of subjects with\nabilities retained',
                         legend.title = 'LoA / Stance Items',
                         risk.table       = T, risk.table.fontsize = .1, tables.height = 0.4,
                         fun = 'pct'
)

vs <- c(res.tmp %>% filter(paramcd %in% c('fane3a','fane2b','fane2a','fane7')))$median

p3$plot <- p3$plot +
  geom_hline(yintercept = 50, linetype = 'dashed')+
  geom_vline(xintercept = vs, linetype = 'dashed')

# output ------------------------------------------------------------------

title <- NA

read_pptx("../Templates/CR.template.pptx") %>%
  add_slide   ( layout = "F", master = "CR") %>%
  ph_with(value = format(Sys.time()), location = ph_location_type(type = "dt")) %>% 
  ph_with_vg  ( print(p0, newpage = F) , location = ph_location(left = 0.97, top = 1.314, width = 16/2.54, height = 14.4/2.54) ) %>%
  add_slide   ( layout = "F", master = "CR") %>%
  ph_with(value = format(Sys.time()), location = ph_location_type(type = "dt")) %>% 
  ph_with_vg  ( print(p1, newpage = F) , location = ph_location(left = 0.97, top = 1.314, width = 16/2.54, height = 17.0/2.54) ) %>%
  add_slide   ( layout = "F", master = "CR") %>%
  ph_with(value = format(Sys.time()), location = ph_location_type(type = "dt")) %>% 
  ph_with_vg  ( print(p2, newpage = F) , location = ph_location(left = 0.97, top = 1.314, width = 16/2.54, height = 17.0/2.54) ) %>%
  add_slide   ( layout = "F", master = "CR") %>%
  ph_with(value = format(Sys.time()), location = ph_location_type(type = "dt")) %>% 
  ph_with_vg  ( print(p3, newpage = F) , location = ph_location(left = 0.97, top = 1.314, width = 16/2.54, height = 17.0/2.54) ) %>%
  print       ( target = paste('code/', project., '.', out., '.', gsub(":","-", Sys.time()), " - ", title,".pptx", sep="") )

.sp
