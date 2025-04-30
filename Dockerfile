# Base image
FROM rocker/shiny:4.3.1

# General updates
# System dependencies for CRAN packages
RUN apt-get update && apt-get install -y \
    libudunits2-dev \
    libgdal-dev \
    libgeos-dev \
    libproj-dev \
    libxml2-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    libfontconfig1-dev \
    libmagick++-dev \
    git \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Create working directory
WORKDIR /srv/shiny-server/

RUN rm -rf /srv/shiny-server/*

# Install the required packages
# Recreate the R environment using renv package
RUN Rscript -e 'install.packages(c("sf", "units", "xml2", "tidyverse", "viridis", "smoothr", "RImageJROI", "utils" ,"codetools", "purrr"), dependencies = TRUE)'

# Copy app code
COPY ./app/app.R /srv/shiny-server/app.R
COPY ./app/www/* /srv/shiny-server/www/
COPY ./app/R/* /srv/shiny-server/R/


# Ensure that the expected user is present in the container
RUN if id shiny &>/dev/null && [ "$(id -u shiny)" -ne 999 ]; then \
        userdel -r shiny; \
        id -u 999 &>/dev/null && userdel -r "$(id -un 999)"; \
    fi; \
    useradd -u 999 -m -s /bin/bash shiny; \
    chown -R shiny:shiny /srv/shiny-server/ /var/lib/shiny-server/ /var/log/shiny-server/

# Other settings
USER shiny
# Open ports (Shiny default)
EXPOSE 3838

# Startup command
CMD ["shiny-server"]
