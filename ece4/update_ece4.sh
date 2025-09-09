#!/bin/bash

# commands to update the local main branch to the upstream
# use with caution

git checkout main
git submodule foreach --recursive git fetch --all
git fetch --all
git merge upstream main
git submodule update --recursive
