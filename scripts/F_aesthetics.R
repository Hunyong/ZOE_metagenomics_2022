### color palette
# data.frame(genus = coef.plot[[1]]$data$genus, color = ggplot_build(coef.plot[[1]])$data[[1]]$colour) %>% unique %>% arrange(genus)
col.genera = 
  c(Actinomyces = "#F8766D", Arcanobacterium = "#ED8141", Bradyrhizobium = "#E08B00", 
    Campylobacter = "#CF9400", Candida = "#BB9D00", Cardiobacterium = "#A3A500", 
    Centipeda = "#85AD00", Lachnoanaerobaculum = "#5BB300", Lachnospiraceae = "#00B81F", 
    Lancefieldella = "#00BC59", Leptotrichia = "#00BF7D", Megasphaera = "#00C19C", 
    Mitsuokella = "#00C0B8", Mogibacterium = "#00BDD0", Neisseria = "#00B8E5", 
    Olsenella = "#00B0F6", Oribacterium = "#00A5FF", Prevotella = "#7997FF", 
    Scardovia = "#AC88FF", Selenomonas = "#CF78FF", Slackia = "#E76BF3", 
    Solobacterium = "#F763E0", Stomatobaculum = "#FF61C9", 
    Streptococcus = "#FF65AE", Veillonella = "#FF6C91")


col.species = 
  c(`Leptotrichia wadei` = "#00BF7D", 
    `Prevotella salivae` = "#7997FF", 
    `Selenomonas sputigena` = "#CF78FF", 
    `Streptococcus mutans` = "#FF65AE")


# data.frame(size = ggplot_build(coef.plot[[1]])$data[[1]]$size, logq = coef.plot[[1]]$data$zoe1.q %>% {-log(.)}) %>% arrange(size)
size.log.q = c(0, 4)
size.log.p = c(0, 6)
  