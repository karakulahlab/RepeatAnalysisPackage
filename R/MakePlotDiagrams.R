
drawPlotESx_UpGene<-function(esx){ #to plot diagram to Enrichment Score table for uploaded genes with using plotly library
  p <- plot_ly(
    type = 'scatter',
    x = esx$r,
    y = esx$Repeats,
    text = paste("Repeats: ", esx$Repeats,
                 "<br>r: ",esx$r,
                 "<br>g: ", esx$g,
                 "<br>R: ",esx$R,
                 "<br>G: ", esx$G,
                 "<br>ESx: ", esx$ESx,
                 "<br>pValue: ",esx$pValue,
                 "<br>pAdjValue: ",esx$pAdjust),
    hoverinfo = 'text',
    mode = 'markers',
    transforms = list(
      list(
        type = 'groupby',
        groups = esx$r,
        styles = list(
          list(target = 1, value = list(marker =list(color = 'blue'))),
          list(target = 3, value = list(marker =list(color = 'red'))),
          list(target = 5, value = list(marker =list(color = 'yellow'))),
          list(target = 7, value = list(marker =list(color = 'ping'))),
          list(target = 13, value = list(marker =list(color = 'green'))),
          list(target = 27, value = list(marker =list(color = 'purple')))
        ))))
  return(p)
}

drawPlotESx_UpCluster<-function(esx){ #to plot diagram to Enrichment Score table for uploaded genes and its modules(clustered genes) with using plotly library
  p <- plot_ly(
    type = 'scatter',
    x = esx$r,
    y = esx$Repeats,
    text = paste("Repeats: ", esx$Repeats,
                 "Modules: ", esx$Modules,
                 "<br>r: ",esx$r,
                 "<br>g: ", esx$g,
                 "<br>R: ",esx$R,
                 "<br>G: ", esx$G,
                 "<br>ESx: ", esx$ESx,
                 "<br>pValue: ",esx$pValue,
                 "<br>pAdjValue: ",esx$pAdjust),
    hoverinfo = 'text',
    mode = 'markers',
    transforms = list(
      list(
        type = 'groupby',
        groups = esx$Modules,
        styles = list(
          list(target = "black", value = list(marker =list(color = 'black'))),
          list(target = "blue", value = list(marker =list(color = 'blue'))),
          list(target = "brown", value = list(marker =list(color = 'brown'))),
          list(target = "cyan", value = list(marker =list(color = 'cyan'))),
          list(target = "green", value = list(marker =list(color = 'green'))),
          list(target = "greenyellow", value = list(marker =list(color = 'greenyellow'))),
          list(target = "grey", value = list(marker =list(color = 'grey'))),
          list(target = "grey60", value = list(marker =list(color = 'grey60'))),
          list(target = "lightcyan", value = list(marker =list(color = 'lightcyan'))),
          list(target = "lightgreen", value = list(marker =list(color = 'lightgreen'))),
          list(target = "lightyellow", value = list(marker =list(color = 'lightgyellow'))),
          list(target = "magenta", value = list(marker =list(color = 'magenta'))),
          list(target = "midnightblue", value = list(marker =list(color = 'midnightblue'))),
          list(target = "pink", value = list(marker =list(color = 'pink'))),
          list(target = "purple", value = list(marker =list(color = 'purple'))),
          list(target = "royalblue", value = list(marker =list(color = 'royalblue'))),
          list(target = "tan", value = list(marker =list(color = 'tan'))),
          list(target = "yellow", value = list(marker =list(color = 'yellow')))
        ))))
  return(p)
}
