#!/bin/bash

# commands to update the local main branch to the upstream
# use with caution

cd /lus/h2resw01/hpcperm/ecme3497/ec-earth-4-fork
git checkout main
git submodule foreach --recursive git fetch --all
git fetch --all
git merge upstream main
git submodule update --recursive
