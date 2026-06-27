# Dockerfile to build a Docker container for NEMESIS
# 
# ~Shubham K - 08/04/2026 - initial release
# ______________________________________________________________________________________________


# Using an older version of Ubuntu that supports older gfortran
FROM ubuntu:16.04

ENV DEBIAN_FRONTEND=noninteractive

# -o pipefail helps catching errors during build
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

# Install development tools and compilers from Ubuntu 16.04
RUN apt-get update && apt-get install -y --no-install-recommends \
    bash \
    build-essential \
    findutils \
    gcc-5 \
    gfortran-5 \
    sed \
    tcsh \
    && rm -rf /var/lib/apt/lists/*

# Env variables required to compile NEMESIS
ENV RADREPO=/app
ENV RADSRC=/app/radtran
ENV BIN=/app/bin
ENV OBJ=/app/obj
ENV LIB=/app/lib
ENV PATH=/app/bin:${PATH}

# Flags
ENV CC=gcc-5
ENV FCOMP=gfortran-5
ENV FC=gfortran-5
ENV STATIC_FLAG="-mcmodel=large -O3 -fno-pic -no-pie"
ENV ARFLAGS=rvU
ENV FCFLAGS1=
ENV FCFLAGS1_FOVGREG=-ffree-form

WORKDIR /app

RUN mkdir -p "${BIN}" "${OBJ}" "${LIB}"

# Copy all files from root dir to /app
COPY . .

# Check the ISYS flag - only as a precaution
RUN sed -i 's/ISYS=1/ISYS=4/g' "${RADSRC}/rtm_util/isys.f" || true

# Older instruction - this line is now redundant, but kept here for backward compatibility
RUN find "${RADSRC}" -type f -name "*.f" -exec sed -i "s/READONLY/ACTION='READ'/g" {} + || true

# In all the makefiles, the line below adds the U flag to avoid timestamp issues
RUN find . -type f \( -name "makefile" -o -name "Makefile" \) -exec sed -i 's/ar rv /ar rvU /g' {} + || true

# The block below is where the actual compilation of NEMESIS happens. Here, makeradtranlib and makeradtranbin 
# are replaced by their individual commands. This helps identify any compilation issues. 'touch *.f' is used
# to update timestamps and force recompilation to avoid stale dependencies. 
RUN set -eux && \
    ulimit -s unlimited && \
    cd /app/FOVgreg && make clean && make lib && \
    cd /app/frecipes && make clean && make lib && \
    cd /app/radtran/rtm_util && make clean && touch *.f && make lib && \
    cd /app/radtran/spec_data && make clean && touch *.f && make lib && \
    cd /app/radtran/path && make clean && touch *.f && make lib && \
    cd /app/radtran/scatter && make clean && touch *.f && make lib && \
    cd /app/radtran/radtran && make clean && touch *.f && make lib && \
    cd /app/radtran/ciatable && make clean && touch *.f && make lib && \
    cd /app/radtran/cirsrad && make clean && touch *.f && make lib && \
    cd /app/radtran/cirsradg && make clean && touch *.f && make lib && \
    cd /app/radtran/matrices && make clean && touch *.f && make lib && \
    cd /app/nemesis && make clean && make lib && \
    cd /app/radtran/rtm_util && make clean && make bin && \
    cd /app/radtran/spec_data && make clean && make bin && \
    cd /app/radtran/path && make clean && make bin && \
    cd /app/radtran/scatter && make clean && make bin && \
    cd /app/radtran/radtran && make clean && make bin && \
    cd /app/radtran/ciatable && make clean && make bin && \
    cd /app/radtran/cirsrad && make clean && make bin && \
    cd /app/radtran/cirsradg && make clean && make bin && \
    cd /app/nemesis && make clean && make bin

COPY entrypoint.sh /usr/local/bin/entrypoint.sh
RUN chmod +x /usr/local/bin/entrypoint.sh

# Call entrypoint.sh when 'docker run' is used. 
ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
CMD ["/bin/bash"]
