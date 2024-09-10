FROM gitlab-registry.internal.sanger.ac.uk/dermatlas/analysis-methods/maf:latest as base
RUN mv /opt/repo /opt/maf

USER root

# Set the top level environment variables
ENV \
    DATA_DIRECTORY="/data" \
    OPT_DIRECTORY="/opt" \
    USER_NAME="admin" \
    USER_DIRECTORY="/home/admin" \
    LC_ALL="en_US.UTF-8" \
    LANG="en_US.UTF-8" 
    
# Set next environment variables that interpolate the top level environment
# variables
ENV \
    USER_BASHRC="${USER_DIRECTORY:?}/.bashrc" \
    USER_BIN_DIRECTORY="${USER_DIRECTORY:?}/.local/bin" \
    SSH_DIR="${USER_DIRECTORY:?}/.ssh" \
    PROJECT_DIRECTORY="${OPT_DIRECTORY:?}/repo" \
    RENV_DIRECTORY="${OPT_DIRECTORY:?}/renv" \
    RENV_PATHS_CACHE="${OPT_DIRECTORY:?}/renv-cache" \
    LOGGING_DIRECTORY="${DATA_DIRECTORY:?}/logs" 
# Set the environment variables for the versions of the software 
ENV \
    RENV_PATHS_LIBRARY="${RENV_DIRECTORY:?}/library"


# Copy the rest of the files into the image (except the .Dockerignore files)
COPY --chown="${USER_NAME}:${USER_NAME}" . .

RUN chown -R "${USER_NAME:?}:${USER_NAME:?}" "${PROJECT_DIRECTORY:?}" "${DATA_DIRECTORY:?}" "${USER_DIRECTORY:?}" "${OPT_DIRECTORY:?}" "/opt/maf"
# # Switch to the non-root user
USER "${USER_NAME:?}"


WORKDIR $PROJECT_DIRECTORY
CMD ["/bin/bash"]
