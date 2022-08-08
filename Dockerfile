FROM dukegcb/openshift-shiny-verse:3.6.3
RUN install2.r -e shiny dplyr tidyr ggplot2 RSQLite ggdag UpSetR igraph ggrepel corrplot DT randomcoloR
ADD ./BioErp /srv/code
