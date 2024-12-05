# Install necessary system dependencies for R packages
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev

# Install required R packages
RUN R -e "install.packages(c('shiny', 'shinythemes', 'shinyWidgets', 'dplyr', 'DT', 'plotly', 'bslib', 'seqinr', 'httr', 'jsonlite', 'shinyjs'))"

# Copy the Shiny app to the image
COPY shiny-app /srv/shiny-server/

# Set proper permissions for the shiny app
RUN chmod -R 755 /srv/shiny-server

# Expose port 3838
EXPOSE 3838

# Start the Shiny server
CMD ["/usr/bin/shiny-server"]
