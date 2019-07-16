## Code to generate the plots used in my thesis.
## Requires that the code in thesis_analysis.R has been run and that the resulting datasets and functions
## are in the global environment.


## Plot constants and save function ----

theme_spacing <- theme(axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)),
                       axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)))
lib_theme <- theme_light() + theme(text=element_text(family="Linux Libertine G", size=12))

save_plot <- function(plotname, plotwidth, plotheight){
  path <- paste0("r-analysis/plots/", deparse(substitute(plotname)), ".pdf")
  ggsave(filename=path, plot=plotname,
         device=cairo_pdf,
         width=plotwidth, height=plotheight, units = "in")
}

## Fig. 3: Scatterplot of types for each suffix by year ----

er_doc_info_norm_typ <- data.frame(types = er_doc_info_norm$types, 
                                   types.norm = er_doc_info_norm$types.norm, 
                                   date = er_doc_info_norm$date)
er_doc_info_norm_typ_tall <- melt(er_doc_info_norm_typ, id="date")

hkeit_doc_info_norm_typ <- data.frame(types = hkeit_doc_info_norm$types, 
                                      types.norm = hkeit_doc_info_norm$types.norm, 
                                      date = hkeit_doc_info_norm$date)
hkeit_doc_info_norm_typ_tall <- melt(hkeit_doc_info_norm_typ, id="date")

ung_doc_info_norm_typ <- data.frame(types = ung_doc_info_norm$types, 
                                    types.norm = ung_doc_info_norm$types.norm, 
                                    date = ung_doc_info_norm$date)
ung_doc_info_norm_typ_tall <- melt(ung_doc_info_norm_typ, id="date")

sfxs_doc_info_norm_typ <- dplyr::bind_rows(list(er=er_doc_info_norm_typ_tall, 
                                                hkeit=hkeit_doc_info_norm_typ_tall, 
                                                ung=ung_doc_info_norm_typ_tall),
                                           .id = 'sfx')

plot_raw_types_year <- ggplot(subset(sfxs_doc_info_norm_typ, variable="types"),
       aes(date, value, colour=sfx)) +
  geom_point() +
  geom_smooth(se=F) +
  labs(x="Year", y="Types (V)") +
  lib_theme + 
  # scale_colour_brewer(name="Suffix", labels=c("-er", "-heit/-keit", "-ung"),
                      # palette = "Set1") +
  scale_colour_grey(name="Suffix", labels=c("-er", "-heit/-keit", "-ung"))

save_plot(plot_raw_types_year, 4, 3)


## Fig. 4: Scatterplot of text length for each suffix by year ----

plot_tokens_texts_year <- ggplot(ung_doc_info_norm, aes(date, corpus.actual)) +
  geom_point() +
  geom_smooth(se=F, colour='darkgrey') +
  scale_y_continuous(trans="log10") +
  labs(x="Year", y="Text length (tokens on dipl)") +
  # scale_colour_brewer(palette = "Set1") +
  scale_colour_grey() +
  theme_spacing + 
  lib_theme

save_plot(plot_tokens_texts_year, 3.5, 2.5)


## Fig. 5: Bar plot of new types per 100-year subperiod ----

ung_sp_roa_info <- data.frame(tokens = ung_sp_info$tokens,
                              new.types = ung_sp_roa$new.types,
                              subset = ung_sp_info$subset)
ung_sp_roa_info_tall <- melt(ung_sp_roa_info, id="subset")

er_sp_roa_info <- data.frame(tokens = er_sp_info$tokens,
                             new.types = er_sp_roa$new.types,
                             subset = er_sp_info$subset)
er_sp_roa_info_tall <- melt(er_sp_roa_info, id="subset")

hkeit_sp_roa_info <- data.frame(tokens = hkeit_sp_info$tokens,
                                new.types = hkeit_sp_roa$new.types,
                                subset = hkeit_sp_info$subset)
hkeit_sp_roa_info_tall <- melt(hkeit_sp_roa_info, id="subset")

sfxs_sp_roa <- dplyr::bind_rows(list(er=er_sp_roa_info_tall, 
                                     hkeit=hkeit_sp_roa_info_tall, 
                                     ung=ung_sp_roa_info_tall),
                                .id = 'sfx')


plot_newtypes_subperiod <- ggplot(subset(sfxs_sp_roa, sfxs_sp_roa$variable=="new.types"),
       aes(subset, value, fill=sfx)) +
  geom_bar(stat='identity', position='dodge') +
  labs(x="Subperiod", y="New types") +
  scale_fill_grey(name="Suffix",
                  labels=c("-er", "-heit/-keit", "-ung")) +
  theme_spacing +
  lib_theme

save_plot(plot_newtypes_subperiod, 5, 3)


## Fig. 6: Faceted scatterplots of actual and normalised type counts ----

sfxs_doc_info_norm <- dplyr::bind_rows(list(er=er_doc_info_norm,
                                            hkeit=hkeit_doc_info_norm,
                                            ung=ung_doc_info_norm),
                                       .id= 'sfx')

sfxs_doc_info_norm_tog <- data.frame(types = sfxs_doc_info_norm$types,
                                     types.norm = sfxs_doc_info_norm$types.norm,
                                     corpus.actual = sfxs_doc_info_norm$corpus.actual,
                                     sfx = sfxs_doc_info_norm$sfx)
sfxs_doc_info_norm_tog_tall <- melt(sfxs_doc_info_norm_tog, id.vars = c("sfx", "corpus.actual"))

collabels <- c(er = "-er", hkeit = "-heit/-keit", ung="-ung")
rowlabels <- c(types = "Actual types", types.norm = "Normalised types")


plot_faceted_types_norm_textlen <- ggplot(sfxs_doc_info_norm_tog_tall, aes(x=corpus.actual, y=value, colour=sfx)) +
  geom_smooth(se=F) +
  geom_point() +
  facet_grid(variable ~ sfx, labeller = labeller(.rows = rowlabels, .cols = collabels)) +
  lib_theme +
  # scale_colour_brewer(palette = "Set1") +
  scale_colour_grey() +
  labs(x="Text length", y = "Types") +
  theme(legend.position = 'none',
        strip.text.x = element_text(colour="black"),
        strip.text.y = element_text(colour="black"))

save_plot(plot_faceted_types_norm_textlen, 6, 4)


## Fig. 7: Dependency of normalised type counts on N ----

max_corpus_size <- 10000
norm_type_df <- data.frame(matrix(nrow=max_corpus_size, ncol=2))
names(norm_type_df) <- c("corpus.size", "norm.types")
norm_type_df$corpus.size <- 1:max_corpus_size
norm_type_df$norm.types <- normalise(100, norm_type_df$corpus.size, 500)

plot_normed_type_dependency <- ggplot(norm_type_df, aes(corpus.size, norm.types)) +
  geom_line() +
  theme_classic() +
  theme(axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        text=element_text(family="Linux Libertine G", size=12),
        axis.line = element_line(arrow = arrow(length = unit(0.1, "inches")))) +
  coord_cartesian(ylim = c(0, 150), xlim = c(0, 9500)) +
  labs(x="Original subcorpus size", y="Normalised type count") +
  theme_spacing

save_plot(plot_normed_type_dependency, 3.5, 3)


## Fig. 8 and Fig. 18: Monte Carlo plots for potential productivity ----

## -er
plot_mc_pp_er_nosps <- ggplot(er_ppmc_stats, aes(n.seen, mean)) +
  geom_errorbar(data=er_ppmc_stats, aes(ymin=lower, ymax=upper), colour="darkgrey") +
  geom_line() + geom_point() +
  # scale_x_continuous(trans="log10") +
  labs(title="-er", x="N", y=paste0("Potential productivity (2000 iterations)")) +
  lib_theme +
  # scale_colour_brewer(palette="Set1") +
  scale_colour_grey() +
  theme(legend.position = "none") +
  theme_spacing

plot_mc_pp_er_wsps <- plot_mc_pp_er_nosps +
  geom_point(data=er_sp_info, aes(x=tokens, y=pp, colour="sfx")) +
  geom_text_repel(data=er_sp_info,
                  family="Linux Libertine G",
                  aes(x=tokens, y=pp, label=er_sp_info$subset, colour="sfx"),
                  box.padding = 0.5,
                  force=2)

save_plot(plot_mc_pp_er_nosps, 2.5, 4)
save_plot(plot_mc_pp_er_wsps, 2.5, 4)


## -hkeit

plot_mc_pp_hkeit_nosps <- ggplot(hkeit_ppmc_stats, aes(n.seen, mean)) +
  geom_errorbar(data=hkeit_ppmc_stats, aes(ymin=lower, ymax=upper), colour="darkgrey") +
  geom_line() + geom_point() +
  # scale_x_continuous(trans="log10") +
  labs(title="-heit/-keit", x="N", y=paste0("Potential productivity (2000 iterations)")) +
  lib_theme +
  # scale_colour_brewer(palette="Set1") +
  scale_colour_grey() +
  theme(legend.position = "none") +
  theme_spacing

plot_mc_pp_hkeit_wsps <- plot_mc_pp_hkeit_nosps +
  geom_point(data=hkeit_sp_info, aes(x=tokens, y=pp, colour="sfx")) +
  geom_text_repel(data=hkeit_sp_info,
                  family="Linux Libertine G",
                  aes(x=tokens, y=pp, label=hkeit_sp_info$subset, colour="sfx"),
                  box.padding = 0.5,
                  force=2)

save_plot(plot_mc_pp_hkeit_nosps, 2.5, 4)
save_plot(plot_mc_pp_hkeit_wsps, 2.5, 4)


## -ung

plot_mc_pp_ung_nosps <- ggplot(ung_ppmc_stats, aes(n.seen, mean)) +
  geom_errorbar(data=ung_ppmc_stats, aes(ymin=lower, ymax=upper), colour="darkgrey") +
  geom_line() + geom_point() +
  # scale_x_continuous(trans="log10") +
  labs(title="-ung", x="N", y=paste0("Potential productivity (2000 iterations)")) +
  lib_theme +
  # scale_colour_brewer(palette="Set1") +
  scale_colour_grey() +
  theme(legend.position = "none") +
  theme_spacing

plot_mc_pp_ung_wsps <- plot_mc_pp_ung_nosps +
  geom_point(data=ung_sp_info, aes(x=tokens, y=pp, colour="sfx")) +
  geom_text_repel(data=ung_sp_info,
                  family="Linux Libertine G",
                  aes(x=tokens, y=pp, label=ung_sp_info$subset, colour="sfx"),
                  box.padding = 0.5,
                  force=2)

save_plot(plot_mc_pp_ung_nosps, 2.5, 4)
save_plot(plot_mc_pp_ung_wsps, 2.5, 4)


## Fig. 9: VGC for first twenty -er derivations in RIDGES ----

er_20_vgc <- data.frame(matrix(nrow=20, ncol=2))
names(er_20_vgc) <- c("N", "V")
er_20_vgc$N <- 1:20
er_20_vgc$V <- c(1, 2, 2, 3, 4, 4, 4, 5, 5, 6, 7, 8, 9, 9, 10, 11, 11, 12, 13, 14)

plot_er_twenty_vgc <- ggplot(er_20_vgc, aes(N, V)) +
  geom_line() +
  lib_theme +
  scale_y_continuous(breaks=seq(2, 14, 2)) +
  scale_x_continuous(breaks=seq(2, 20, 2)) +
  theme_spacing

save_plot(plot_er_twenty_vgc, 3, 4)


## Fig. 10: Idealised productive and unproductive VGCs ----

make_ideal_vgc <- function(wdlist){
  freq_list <- as.data.frame(table(wdlist))
  freq_spc <- as.data.frame(table(freq_list$Freq))
  
  # ## Rename the headers the way zipfR wants them.
  names(freq_spc) <- c("m", "Vm")
  # 
  # ## Write dataframe to external .spc file.
  write.table(freq_spc, "r-analysis/data.spc", sep="\t", row.names=FALSE)
  
  # ## Read data.spc back in with read.spc; now the spectrum is ready for use with zipfR.
  ideal_prod_spc <- read.spc("r-analysis/data.spc")
  
  ideal_prod_spc_vgc <- vgc.interp(ideal_prod_spc, 1:N(ideal_prod_spc))
  
  return(ideal_prod_spc_vgc)
}

f <- 6

prod_list <- c(rep("a", round2(f/1, 0)), 
               rep("b", round2(f/1.1, 0)), 
               rep("c", round2(f/1.2, 0)), 
               rep("d", round2(f/1.5, 0)), 
               rep("e", round2(f/2, 0)), 
               rep("f", round2(f/3, 0)), 
               rep("g", round2(f/5, 0)), 
               rep("h", round2(f/8, 0)), 
               rep("i", round2(f/12, 0)), 
               rep("j", round2(f/17, 0)))

unprod_list <- c(rep("a", round2(f/1, 0)), 
                 rep("b", round2(f/.9, 0)), 
                 rep("c", round2(f/.8, 0)), 
                 rep("d", round2(f/.7, 0)), 
                 rep("e", round2(f/.6, 0)), 
                 rep("f", round2(f/.5, 0)), 
                 rep("g", round2(f/.4, 0)), 
                 rep("h", round2(f/.3, 0)), 
                 rep("i", round2(f/.2, 0)), 
                 rep("j", round2(f/.1, 0)))


ideal_prod_spc_vgc <- make_ideal_vgc(prod_list)
ideal_unprod_spc_vgc <- make_ideal_vgc(unprod_list)

plot_ideal_prod_vgc <- ggplot(ideal_prod_spc_vgc, aes(N, V)) + 
  geom_line() +
  theme_classic() +
  theme(axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        text=element_text(family="Linux Libertine G", size=12),
        axis.line = element_line(arrow = arrow(length = unit(0.05, "inches")))) +
  ggtitle("productive")


plot_ideal_unprod_vgc <- ggplot(ideal_unprod_spc_vgc, aes(N, V)) + 
  geom_line() +
  theme_classic() +
  coord_cartesian(ylim=c(0, 18)) +
  theme(axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        text=element_text(family="Linux Libertine G", size=12),
        axis.line = element_line(arrow = arrow(length = unit(0.05, "inches")))) +
  ggtitle("unproductive")

save_plot(plot_ideal_prod_vgc, 2, 2)
save_plot(plot_ideal_unprod_vgc, 2, 2)


## Fig. 11: Empirical VGCs for each suffix per subperiod ----

plot_er_emp_vgc_sps <- ggplot(er_sp_emp_vgc_tall, aes(N, V, colour=period)) +
  geom_line() +
  ggtitle("-er") +
  lib_theme +
  # scale_colour_brewer(palette = "Set1",
                      # name="Subperiod") +
  scale_colour_grey(name="Subperiod") +
  theme_spacing

plot_hkeit_emp_vgc_sps <- ggplot(hkeit_sp_emp_vgc_tall, aes(N, V, colour=period)) +
  geom_line() +
  ggtitle("-heit/-keit") +
  lib_theme +
  # scale_colour_brewer(palette = "Set1",
  # name="Subperiod") +
  scale_colour_grey(name="Subperiod") +
  theme_spacing

plot_ung_emp_vgc_sps <- ggplot(ung_sp_emp_vgc_tall, aes(N, V, colour=period)) +
  geom_line() +
  ggtitle("-ung") +
  lib_theme +
  # scale_colour_brewer(palette = "Set1",
  # name="Subperiod") +
  scale_colour_grey(name="Subperiod") +
  theme_spacing

save_plot(plot_er_emp_vgc_sps, 3, 3)
save_plot(plot_hkeit_emp_vgc_sps, 3, 3)
save_plot(plot_ung_emp_vgc_sps, 3, 3)


## Fig. 12: Empirical and interpolated VGCs for full samples ----

## Empirical VGCs:
sfxs_conc_emp_vgc <- dplyr::bind_rows(list(er=er_conc_emp_vgc, 
                                           hkeit=hkeit_conc_emp_vgc, 
                                           ung=ung_conc_emp_vgc),
                                      .id = 'sfx')

plot_sfxs_conc_emp_vgc <- ggplot(sfxs_conc_emp_vgc, aes(N, V, colour=sfx)) +
  geom_line() +
  ggtitle("Empirical VGCs") +
  lib_theme +
  scale_color_grey(name="Suffix", labels=c("-er", "-heit/-keit", "-ung")) +
  theme_spacing

save_plot(plot_sfxs_conc_emp_vgc, 4, 3)


## Interpolated VGCs:
sfxs_conc_bin_vgc <- dplyr::bind_rows(list(er=er_conc_bin_vgc, 
                                           hkeit=hkeit_conc_bin_vgc, 
                                           ung=ung_conc_bin_vgc),
                                      .id = 'sfx')

plot_sfxs_conc_bin_vgc <- ggplot(sfxs_conc_bin_vgc, aes(N, V, colour=sfx)) +
  geom_line() +
  ggtitle("Interpolated VGCs") +
  lib_theme +
  scale_color_grey(name="Suffix", labels=c("-er", "-heit/-keit", "-ung")) +
  theme_spacing +
  labs(y="E[V]")

save_plot(plot_sfxs_conc_bin_vgc, 4, 3)


## Fig. 13: Two curves, same slope ----

## This one is a Big Mess because a lot of this code came from the internet 
## and was cobbled together in gross ways. Not going to invest the time to 
## clean it up.

x <- seq(0,40)
y <- dnorm(seq(0,40), mean=25, sd=5)

## Getting deriv of prod curve ready.

spl <- smooth.spline(ideal_prod_spc_vgc$N, ideal_prod_spc_vgc$V, spar=0.3)
newx <- seq(min(ideal_prod_spc_vgc$N), max(ideal_prod_spc_vgc$N), 0.1)
pred <- predict(spl, x=newx, deriv=0)

newx <- 10
pred0 <- predict(spl, x=newx, deriv=0)
pred1 <- predict(spl, x=newx, deriv=1)
yint <- pred0$y - (pred1$y*newx)
xint <- -yint/pred1$y

mod <- lm(yint + pred1$y*x ~ x)
abline(mod)
# coef(mod)[2] ## The slope of the tangent line at this point.


## Getting deriv of unprod curve ready.

test1 <- ideal_unprod_spc_vgc
test1$V <- ideal_unprod_spc_vgc$V * 1.0625

spl2 <- smooth.spline(test1$N, test1$V, spar=0.3)
newx2 <- seq(min(test1$N), max(test1$N), 0.1)
pred2 <- predict(spl2, x=newx2, deriv=0)

newx2 <- 10
pred02 <- predict(spl2, x=newx2, deriv=0)
pred12 <- predict(spl2, x=newx2, deriv=1)
yint2 <- pred02$y - (pred1$y*newx2)
xint2 <- -yint2/pred1$y

mod2 <- lm(yint2 + pred12$y*x ~ x)
abline(mod2, col=3)
# print(coef(mod2)[2]) ## The slope of the tangent line at this point.


# Prep for plotting

curve2 <- test1[1:28, ]
curve1 <- ideal_prod_spc_vgc

curves_df <- data.frame(N = curve1$N,
                        V1 = curve1$V,
                        V2 = curve2$V)

plot_same_pp <- ggplot(curves_df, aes(x=N)) +
  geom_segment(aes(y=0, yend=5.9, x=10, xend=10), linetype='dotted') +
  
  geom_line(aes(y=V1), colour="#cccccc") +
  geom_abline(slope=coef(mod)[2], intercept=coef(mod)[1], colour="#cccccc", linetype="dashed") +
  geom_point(aes(x=10, y=curves_df$V1[10]), colour="#cccccc") +
  
  geom_line(aes(y=V2), colour="#333333") +
  geom_abline(slope=coef(mod2)[2], intercept=coef(mod2)[1], colour="#333333", linetype="dashed") +
  geom_point(aes(x=10, y=curves_df$V2[10]), colour="#333333") +
  
  theme_classic() +
  theme_spacing +
  theme(axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        text=element_text(family="Linux Libertine G", size=12),
        axis.line = element_line(arrow = arrow(length = unit(0.1, "inches")))) +
  labs(y="V") +
  coord_cartesian(xlim=c(1, 15), ylim=c(1, 7))

save_plot(plot_same_pp, 4, 3)


## Fig. 14: S for two VGCS ----

s_comp_vgc_df <- ideal_unprod_spc_vgc
s_comp_vgc_df$V2 <- 2*s_comp_vgc_df$V
s_comp_vgc_df$V2[1] <- s_comp_vgc_df$V[1]

s_comp_vgc_df_tall <- melt(s_comp_vgc_df, id.vars=c("N"))

plot_comparing_s <- ggplot(subset(s_comp_vgc_df_tall, s_comp_vgc_df_tall$N <= 100), aes(N, value, colour=variable)) + 
  geom_line() +
  theme_classic() +
  theme(axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        text=element_text(family="Linux Libertine G", size=12),
        axis.line = element_line(arrow = arrow(length = unit(0.1, "inches")))) +
  scale_colour_grey() +
  labs(y="E[V]") + 
  theme(legend.position = "none",
        axis.text = element_blank()) +
  geom_hline(yintercept=10.001, linetype='dashed') +
  geom_hline(yintercept=20.001, linetype='dashed') +
  annotate("text", x = 2, y = 21, label = expression(S[1]), family="Linux Libertine G") +
  annotate("text", x = 2, y = 11, label = expression(S[2]), family="Linux Libertine G")

save_plot(plot_comparing_s, 4, 3)


## Fig. 15: Type frequency plots for each suffix ----

er_top12_freq <- head(er_conc_freqlist[order(er_conc_freqlist$frequency, decreasing=TRUE),], 12)
hkeit_top12_freq <- head(hkeit_conc_freqlist[order(hkeit_conc_freqlist$frequency, decreasing=TRUE),], 12)
ung_top12_freq <- head(ung_conc_freqlist[order(ung_conc_freqlist$frequency, decreasing=TRUE),], 12)

plot_er_freqlist <- ggplot(data=er_top12_freq, aes(x=lemma,y=frequency)) +
  geom_bar(position="dodge",stat="identity") +
  coord_flip() +
  ggtitle("-er") +
  labs(x=element_blank(), y=element_blank()) +
  lib_theme + theme_spacing +
  scale_x_discrete(limits=rev(er_top12_freq$lemma),
                   labels=c("Kräutler", "Gräber", "Forscher", "Fehler", "Römer", "Italiener",
                            "Blümler", "Botaniker", "Gärtner", "Träger", "Apotheker", "Leser"))

plot_hkeit_freqlist <- ggplot(data=hkeit_top12_freq, aes(x=lemma,y=frequency)) +
  geom_bar(position="dodge",stat="identity") +
  coord_flip() +
  ggtitle("-heit/-keit") +
  labs(x=element_blank(), y=element_blank()) +
  lib_theme + theme_spacing +
  scale_x_discrete(limits=rev(hkeit_top12_freq$lemma),
                   labels=c("Reinigkeit", "Sonderheit", "Keuschheit", "Wahrheit", "Bitterkeit", "Gesundheit", 
                            "Engbrüstigkeit", "Gelegenheit", "Beschaffenheit", "Ähnlichkeit", "Feuchtigkeit", "Krankheit"))

plot_ung_freqlist <- ggplot(data=ung_top12_freq, aes(x=lemma,y=frequency)) +
  geom_bar(position="dodge",stat="identity") +
  coord_flip() +
  ggtitle("-ung") +
  labs(x=element_blank(), y=element_blank()) +
  lib_theme + theme_spacing +
  scale_x_discrete(limits=rev(ung_top12_freq$lemma),
                   labels=c("Beobachtung", "Entwicklung", "Entzündung", "Befruchtung", "Abteilung", "Meinung",
                            "Zusammensetzung", "Verstopfung", "Bildung", "Benennung", "Wirkung", "Gattung"))

save_plot(plot_er_freqlist, 2, 4)
save_plot(plot_hkeit_freqlist, 2.25, 4)
save_plot(plot_ung_freqlist, 2.25, 4)


## Fig. 16 and Fig. 18: Monte Carlo plots for S ----

## -er:
plot_mc_s_er_nosps <- ggplot(er_mc_approx_stats_no_outliers, aes(n.seen, mean)) +
  geom_errorbar(data=er_mc_approx_stats_no_outliers, aes(ymin=lower, ymax=upper), colour="darkgrey") +
  geom_line() + geom_point() +
  coord_cartesian(ylim = c(0, 900)) +
  # scale_x_continuous(trans="log10") +
  labs(title="-er", x="N", y=paste0("S (", min(er_mc_approx_overview$data.pts), " iterations)")) +
  lib_theme + theme_spacing +
  # scale_colour_brewer(palette="Set1") +
  scale_colour_grey() +
  theme(legend.position = "none")

plot_mc_s_er_wsps <- plot_mc_s_er_nosps +
  geom_point(data=er_sp_info, aes(x=tokens, y=S.approx, colour="sfx")) +
  geom_text_repel(data=er_sp_info,
                  aes(x=tokens, y=S.approx, label=er_sp_info$subset, colour="sfx"),
                  box.padding = 0.5,
                  force=2,
                  family="Linux Libertine G")

save_plot(plot_mc_s_er_nosps, 2.5, 4)
save_plot(plot_mc_s_er_wsps, 2.5, 4)


## -hkeit:
plot_mc_s_hkeit_nosps <- ggplot(hkeit_mc_approx_stats_no_outliers, aes(n.seen, mean)) +
  geom_errorbar(data=hkeit_mc_approx_stats_no_outliers, aes(ymin=lower, ymax=upper), colour="darkgrey") +
  geom_line() + geom_point() +
  coord_cartesian(ylim = c(0, 400)) +
  # scale_x_continuous(trans="log10") +
  labs(title="-heit/-keit", x="N", y=paste0("S (", min(hkeit_mc_approx_overview$data.pts), " iterations)")) +
  lib_theme + theme_spacing +
  # scale_colour_brewer(palette="Set1") +
  scale_colour_grey() +
  theme(legend.position = "none")

plot_mc_s_hkeit_wsps <- plot_mc_s_hkeit_nosps +
  geom_point(data=hkeit_sp_info, aes(x=tokens, y=S.approx, colour="ung")) +
  geom_text_repel(data=hkeit_sp_info,
                  aes(x=tokens, y=S.approx, label=hkeit_sp_info$subset, colour="ung"),
                  box.padding = 0.5,
                  force=2,
                  family="Linux Libertine G")

save_plot(plot_mc_s_hkeit_nosps, 2.5, 4)
save_plot(plot_mc_s_hkeit_wsps, 2.5, 4)


## -ung:

plot_mc_s_ung_nosps <- ggplot(ung_mc_approx_stats_no_outliers, aes(n.seen, mean)) +
  geom_errorbar(data=ung_mc_approx_stats_no_outliers, aes(ymin=lower, ymax=upper), colour="darkgrey") +
  geom_line() + geom_point() +
  coord_cartesian(ylim = c(0, 900)) +
  # scale_x_continuous(trans="log10") +
  labs(title="-ung", x="N", y=paste0("S (", min(ung_mc_approx_overview$data.pts), " iterations)")) +
  lib_theme + theme_spacing +
  # scale_colour_brewer(palette="Set1") +
  scale_colour_grey() +
  theme(legend.position = "none")

plot_mc_s_ung_wsps <- plot_mc_s_ung_nosps +
  geom_point(data=ung_sp_info, aes(x=tokens, y=S.approx, colour="ung")) +
  geom_text_repel(data=subset(ung_sp_info, S.approx < 450), 
                  aes(x=tokens, y=S.approx, label=subset(ung_sp_info, S.approx < 450)$subset, colour="ung"),
                  box.padding = 0.5,
                  force=5,
                  nudge_y = 10,
                  family = "Linux Libertine G"
                  # ylim = c(NA, 400)
  ) +
  geom_text_repel(data=subset(ung_sp_info, S.approx > 450), 
                  aes(x=tokens, y=S.approx, label=subset(ung_sp_info, S.approx > 450)$subset, colour="ung"),
                  # point.padding = 0.75,
                  box.padding = 0.5,
                  force=2,
                  nudge_x = 0,
                  ylim = c(800, NA),
                  family = "Linux Libertine G"
  )

save_plot(plot_mc_s_ung_nosps, 2.5, 4)
save_plot(plot_mc_s_ung_wsps, 2.5, 4)


## Fig. 17: Monte Carlo for hkeit from DECOW ----

plot_mc_s_hkeit_decow <- ggplot(hkeitDECOW_mc_approx_stats_no_outliers, aes(n.seen, mean)) +
  geom_errorbar(data=hkeitDECOW_mc_approx_stats_no_outliers, aes(ymin=lower, ymax=upper), colour="darkgrey") +
  geom_line() + geom_point() +
  # coord_cartesian(ylim = c(0, 900)) +
  # scale_x_continuous(trans="log10") +
  labs(x="N", y=paste0("S (", min(hkeitDECOW_mc_approx_overview$data.pts), " iterations)")) +
  lib_theme + theme_spacing +
  theme(legend.position = "none")

save_plot(plot_mc_s_hkeit_decow, 3.5, 3.5)
