set -e
if [ -s changed_tools_chunk.list ]; then
    planemo test --conda_dependency_resolution --conda_auto_install --conda_channels iuc,conda-forge,bioconda,defaul$
elif [ -s changed_repositories_chunk.list ]; then
    while read -r DIR; do
        planemo test --conda_dependency_resolution --conda_auto_install --conda_channels iuc,conda-forge,bioconda,de$
    done < changed_repositories_chunk.list
fi
