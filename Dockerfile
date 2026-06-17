# Base image from Rocker with tidyverse (which includes many system dependencies)
FROM rocker/tidyverse:4.3.1

# Install additional system dependencies commonly needed for bioinformatics
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    git \
    make \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    && rm -rf /var/lib/apt/lists/*

# Install renv specifically
RUN R -e "install.packages('renv', repos = c(CRAN = 'https://cloud.r-project.org'))"

# Set the working directory inside the container
WORKDIR /project

# Copy renv configuration files to leverage Docker layer caching
COPY renv.lock ./
COPY renv/activate.R renv/
COPY .Rprofile ./

# Restore the precise R environment used in the study
RUN R -e "renv::restore()"

# Copy the rest of the project files
COPY . .

# Set default command to an interactive bash shell for researchers to run the pipeline
CMD ["bash"]
