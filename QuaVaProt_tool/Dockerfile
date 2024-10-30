# Use the rocker/shiny image as a base
FROM rocker/shiny:latest

# Install required R packages
RUN R -e "install.packages(c('shiny', 'shinythemes', 'shinyWidgets', 'dplyr', 'DT', 'bslib', 'httr', 'jsonlite', 'BiocManager', 'bioseq', 'rentrez'))"
RUN R -e "BiocManager::install('UniProt.ws')"

# Copy the Shiny app to the image
COPY shiny-tool /srv/shiny-server/

# Expose port 3838
EXPOSE 3838

# Start the Shiny server
CMD ["/usr/bin/shiny-server"]