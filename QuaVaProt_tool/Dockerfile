# Use the rocker/shiny image as a base
FROM rocker/shiny:latest

# Install required R packages
RUN R -e "install.packages(c('shiny', 'shinythemes', 'shinyWidgets', 'dplyr', 'DT', 'httr', 'jsonlite', 'BiocManager', 'bioseq', 'rentrez', 'devtools', 'callr', 'shinycssloaders'))"
RUN R -e "BiocManager::install('UniProt.ws')"
RUN R -e "require(devtools); install_version('bslib', version = '0.5.1', repos = 'http://cran.us.r-project.org')"

# Copy the Shiny app to the image
COPY shiny-tool /srv/shiny-server/

# Create a data folder inside the container with proper permissions
RUN mkdir -p /srv/shiny-server/data \
    && chown -R shiny:shiny /srv/shiny-server/data \
    && chmod -R 755 /srv/shiny-server/data

#add conf to disable timeouts
COPY shiny-server.conf /etc/shiny-server/shiny-server.conf

# Expose port 3838
EXPOSE 3838

# Start the Shiny server

CMD ["/usr/bin/shiny-server"]