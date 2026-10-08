#!/bin/bash

# commands to update the local main branch to the upstream
# use with caution

cd /hpcperm/${user}/ec-earth-4-fork
git checkout main
git submodule foreach --recursive git fetch --all
git fetch --all
git merge upstream main
git submodule update --recursive
