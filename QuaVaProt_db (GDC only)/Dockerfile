# Use the rocker/shiny image as a base
FROM rocker/shiny:latest

# Install required R packages
RUN R -e "install.packages(c('shiny', 'shinythemes', 'shinyWidgets', 'dplyr', 'DT', 'plotly', 'seqinr', 'httr', 'jsonlite', 'DBI', 'RSQLite', 'shinyjs', 'devtools'))"
RUN R -e "require(devtools); install_version('bslib', version = '0.5.1', repos = 'http://cran.us.r-project.org')"

# Copy the Shiny app to the image
COPY shiny-app /srv/shiny-server/

# Expose port 3838
EXPOSE 3838

# Start the Shiny server
CMD ["/usr/bin/shiny-server"]
