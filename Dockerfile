FROM dukegcb/openshift-shiny-verse:4.1.2
RUN R -e "install.packages(c('shiny','dplyr','tidyr','ggplot2','RSQLite','UpSetR','igraph','ggrepel','corrplot','DT','randomcoloR', 'plotly'), dependencies = T)"
ADD ./BioErp /srv/code
