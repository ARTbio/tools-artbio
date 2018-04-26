#!/usr/bin/env bash

set -e

if [ -s changed_tools_chunk.list ]; then
    planemo test --conda_dependency_resolution --conda_auto_install --conda_channels iuc,bioconda,conda-forge,defaults --galaxy_branch "$GALAXY_RELEASE" --galaxy_source "$GALAXY_REPO" $(cat changed_tools_chunk.list)
elif [ -s changed_repositories_chunk.list ]; then
    while read -r DIR; do
        planemo test --conda_dependency_resolution --conda_auto_install --conda_channels iuc,bioconda,conda-forge,defaults --galaxy_branch "$GALAXY_RELEASE" --galaxy_source "$GALAXY_REPO" "$DIR"
done < changed_repositories_chunk.list
fi

