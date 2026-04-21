FROM ubuntu:16.04

ENV DEBIAN_FRONTEND=noninteractive

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

RUN sed -i \
    -e 's|http://archive.ubuntu.com/ubuntu|http://old-releases.ubuntu.com/ubuntu|g' \
    -e 's|http://security.ubuntu.com/ubuntu|http://old-releases.ubuntu.com/ubuntu|g' \
    /etc/apt/sources.list && \
    printf 'Acquire::Check-Valid-Until "false";\n' > /etc/apt/apt.conf.d/99no-check-valid && \
    apt-get update && apt-get install -y --no-install-recommends \
    bash \
    build-essential \
    findutils \
    gcc-5 \
    gfortran-5 \
    sed \
    tcsh \
    && rm -rf /var/lib/apt/lists/*

ENV RADREPO=/app
ENV RADSRC=/app/radtran
ENV BIN=/app/bin
ENV OBJ=/app/obj
ENV LIB=/app/lib
ENV PATH=/app/bin:${PATH}

ENV CC=gcc-5
ENV FCOMP=gfortran-5
ENV FC=gfortran-5
ENV STATIC_FLAG="-mcmodel=large -O3 -fno-pic -no-pie"
ENV ARFLAGS=rvU
ENV FCFLAGS1=
ENV FCFLAGS1_FOVGREG=-ffree-form

WORKDIR /app

RUN mkdir -p "${BIN}" "${OBJ}" "${LIB}"

COPY . .

RUN sed -i 's/ISYS=1/ISYS=4/g' "${RADSRC}/rtm_util/isys.f" || true
RUN find "${RADSRC}" -type f -name "*.f" -exec sed -i "s/READONLY/ACTION='READ'/g" {} + || true
RUN find . -type f \( -name "makefile" -o -name "Makefile" \) -exec sed -i 's/ar rv /ar rvU /g' {} + || true

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

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
CMD ["/bin/bash"]
