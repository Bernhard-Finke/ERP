FROM dukegcb/openshift-shiny-verse:4.0.2
RUN apt-get update && apt-get install -y \ libglpk40
RUN install2.r -e shiny dplyr tidyr ggplot2 RSQLite ggdag UpSetR igraph ggrepel corrplot DT randomcoloR
ADD ./BioErp /srv/code
