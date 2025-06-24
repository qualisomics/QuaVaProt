# Use the rocker/shiny image as a base
FROM rocker/shiny:latest

# Install required R packages
RUN R -e "install.packages(c('devtools'))"
RUN R -e "devtools::install_version('shiny', version = '1.7.5', repos = 'http://cran.us.r-project.org')"
RUN R -e "devtools::install_version('shinythemes', version = '1.2.0', repos = 'http://cran.us.r-project.org')"
RUN R -e "devtools::install_version('shinyWidgets', version = '0.8.0', repos = 'http://cran.us.r-project.org')"
RUN R -e "devtools::install_version('dplyr', version = '1.1.2', repos = 'http://cran.us.r-project.org')"
RUN R -e "devtools::install_version('DT', version = '0.29', repos = 'http://cran.us.r-project.org')"
RUN R -e "devtools::install_version('plotly', version = '4.10.2', repos = 'http://cran.us.r-project.org')"
RUN R -e "devtools::install_version('seqinr', version = '4.2-30', repos = 'http://cran.us.r-project.org')"
RUN R -e "devtools::install_version('httr', version = '1.4.7', repos = 'http://cran.us.r-project.org')"
RUN R -e "devtools::install_version('jsonlite', version = '1.8.7', repos = 'http://cran.us.r-project.org')"
RUN R -e "devtools::install_version('shinyjs', version = '2.1.0', repos = 'http://cran.us.r-project.org')"
RUN R -e "devtools::install_version('bslib', version = '0.5.1', repos = 'http://cran.us.r-project.org')"

# Copy the Shiny app to the image
COPY shiny-app /srv/shiny-server/

# Expose port 3838
EXPOSE 3838

# Start the Shiny server
CMD ["/usr/bin/shiny-server"]
