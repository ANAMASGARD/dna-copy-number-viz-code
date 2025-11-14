## Modified from https://github.com/tdhock/animint/blob/master/inst/examples/intreg.R
## For issue #40: DNA copy number large margin visualization

library(animint2)
data(intreg)

signal.colors <- c(estimate="#0adb0a", latent="#0098ef")
breakpoint.colors <- c("1breakpoint"="#ff7d7d", "0breakpoints"='#f6f4bf')
model.linetypes <- c(margin="dotted",limit="dashed",regression="solid")
intreg$annotations$logratio <- max(intreg$sig$log)

## Regions with linetype indicating errors.
breaks.by.signal <- split(intreg$breaks, intreg$breaks$signal)
anns.by.signal <- split(intreg$ann, intreg$ann$signal)
error.regions.list <- list()
for(signal in names(breaks.by.signal)){
  signal.breaks <- breaks.by.signal[[signal]]
  signal.anns <- anns.by.signal[[signal]]
  signal.anns$target.breaks <-
    ifelse(signal.anns$annotation=="1breakpoint", 1, 0)
  for(model.i in 1:20){
    model.breaks <- subset(signal.breaks, segments==model.i)
    signal.anns$breaks <- NA
    for(region.i in 1:nrow(signal.anns)){
      region <- signal.anns[region.i, ]
      after.start <- region$first.base < model.breaks$base
      before.end <- model.breaks$base < region$last.base
      signal.anns$breaks[region.i] <- sum(after.start & before.end)
    }
    signal.anns$error.type <- with(signal.anns, {
      ifelse(breaks < target.breaks, "false negative",
             ifelse(target.breaks < breaks, "false positive", "correct"))
    })
    error.regions.list[[paste(model.i, signal)]] <-
      data.frame(segments=model.i, signal.anns)
  }
}
error.regions <- do.call(rbind, error.regions.list)

reg <- subset(intreg$model, line=="regression")
slope <- with(reg, (min.L-max.L)/(min.feature-max.feature))
intreg$intervals$pred.L <-
  slope * (intreg$intervals$feature - reg$min.feature) + reg$min.L

## Main visualization with error annotations
intreg.errors <- 
  animint(
       title="Max-margin interval regression for supervised segmentation model selection",
       source="https://github.com/ANAMASGARD/dna-copy-number-viz-code/blob/main/intreg.R",
       
       signal=ggplot()+
       theme_bw()+
       theme_animint(height=300, width=800, last_in_row=TRUE)+       
       scale_x_continuous("position on chromosome (mega base pairs)",
                          breaks=c(100,200))+
       ylab("noisy copy number logratio signal")+
       geom_tallrect(aes(xmin=first.base/1e6, xmax=last.base/1e6,
                         fill=annotation,
                         linetype=error.type),
                     showSelected=c("signal", "segments"),
                     color="black",
                     alpha=0.5,
                     data=error.regions)+
       scale_linetype_manual("error type", values=c(correct=0,
                                   "false negative"=3,
                                   "false positive"=1))+
       guides(linetype=guide_legend(override.aes=list(fill="white")))+       
       scale_fill_manual(values=breakpoint.colors)+
       geom_blank(aes(first.base/1e6, logratio+2/8), data=intreg$ann)+
       geom_point(aes(base/1e6, logratio),
                  showSelected="signal",
                  data=intreg$sig)+
       geom_segment(aes(first.base/1e6, mean, xend=last.base/1e6, yend=mean),
                    showSelected=c("signal", "segments"),
                    data=intreg$seg, colour=signal.colors[["estimate"]])+
       geom_vline(aes(xintercept=base/1e6),
                  showSelected=c("signal", "segments"),
                  colour=signal.colors[["estimate"]],
                  linetype="dashed",
                  data=intreg$breaks),
       
       penalty=ggplot()+
         theme_bw()+
         theme_animint(height=500, width=800)+
       geom_tallrect(aes(xmin=min.L, xmax=max.L),
                     showSelected="signal",
                     clickSelects="segments",
                     data=intreg$selection,
                     alpha=1/2)+
       ylab("")+
       geom_vline(aes(xintercept=pred.L),
                  showSelected="signal",
                  color="violet",
                  data=intreg$intervals)+
       geom_segment(aes(min.L, feature, xend=max.L, yend=feature),
                    clickSelects="signal",
                    size=6,
                    data=data.frame(intreg$intervals, what="regression"))+
       geom_segment(aes(min.L, min.feature, xend=max.L, yend=max.feature,
                        linetype=line),
                    colour="violet",
                    size=3,
                    data=data.frame(intreg$model, what="regression"))+
       scale_linetype_manual(values=model.linetypes)+
       geom_segment(aes(min.L, cost, xend=max.L, yend=cost),
                    showSelected="signal",
                    data=data.frame(intreg$selection, what="incorrect labels"))+
       geom_segment(aes(min.L, segments, xend=max.L, yend=segments),
                    showSelected="signal",
                    data=data.frame(intreg$selection, what="segments"))+
       xlab("penalty value log(lambda)")+
       facet_grid(what~., scales="free")
  )

print("Visualization created successfully!")
print(paste("Title:", intreg.errors$title))
print(paste("Source:", intreg.errors$source))

