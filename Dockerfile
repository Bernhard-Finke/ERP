FROM dukegcb/openshift-shiny-verse:4.1.2
RUN install2.r -e shiny dplyr tidyr ggplot2 RSQLite ggdag UpSetR igraph ggrepel corrplot DT randomcoloR
ADD ./BioErp /srv/code
