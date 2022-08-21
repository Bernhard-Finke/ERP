FROM dukegcb/openshift-shiny-verse:3.6.2
RUN R -e "install.packages(c('shiny','dplyr','tidyr','ggplot2','RSQLite','UpSetR',c('igraph', version='2.1.6'),'ggrepel','corrplot','DT','randomcoloR', 'plotly'), dependencies = T)"
ADD ./BioErp /srv/code
