#!/usr/bin/env bash
set -e
cd "$TRAVIS_BUILD_DIR" && flake8 --exclude=.git,deprecated,helper_scripts,python_modules,unstable .
while read -r DIR;
    do planemo shed_lint --tools --ensure_metadata --urls --report_level warn --fail_level error --recursive "$DIR";
done < changed_repositories.list
cat changed_repositories.list

