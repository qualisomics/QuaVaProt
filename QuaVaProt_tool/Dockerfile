# Use the rocker/shiny image as a base
FROM rocker/shiny:latest

# Install required R packages
RUN R -e "install.packages(c('shiny', 'shinythemes', 'shinyWidgets', 'dplyr', 'DT', 'httr', 'jsonlite', 'BiocManager', 'bioseq', 'rentrez', 'devtools', 'callr', 'shinycssloaders'))"
RUN R -e "BiocManager::install('UniProt.ws')"
RUN R -e "require(devtools); install_version('bslib', version = '0.5.1', repos = 'http://cran.us.r-project.org')"

# Set working directory
WORKDIR /srv/shiny-server/

# Copy the Shiny app to the image
COPY shiny-tool /srv/shiny-server/

# Create a data folder inside the container with proper permissions
RUN mkdir -p /srv/shiny-server/data && \
    chown -R shiny:shiny /srv/shiny-server/ && \
    chmod -R 775 /srv/shiny-server/

# Ensure the Shiny user owns the app directory
RUN chown -R shiny:shiny /srv/shiny-server/data

#add conf to disable timeouts
COPY shiny-server.conf /etc/shiny-server/shiny-server.conf

# Switch to non-root user
USER shiny

# Expose port 3838
EXPOSE 3838

# Start the Shiny server

CMD ["/usr/bin/shiny-server"]