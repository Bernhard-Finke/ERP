FROM dukegcb/openshift-shiny-verse:4.1.2
RUN sudo apt-get update -y && sudo apt-get install -y libglpk-dev
RUN R -e "install.packages(c('shiny','dplyr','tidyr','ggplot2','RSQLite','UpSetR',c('igraph', version='2.1.6'),'ggrepel','corrplot','DT','randomcoloR', 'plotly'), dependencies = T)"
ADD ./BioErp /srv/code
