FROM dukegcb/openshift-shiny-verse:4.0.2
RUN install2.r -e shiny \
 dplyr \
 tidyr \
 ggplot2 \
 RSQLite \
 UpSetR \
 igraph \
 ggrepel \
 corrplot \
 DT \
 randomcoloR\
ADD ./BioErp /srv/code
