#!/bin/bash

# commands to update the local main branch to the upstream
# use with caution

cd /lus/h2resw01/hpcperm/ccpd/ecearth4/revisions/main/
git checkout main
git submodule foreach --recursive git fetch --all
git fetch --all
git merge upstream main
git submodule update --recursive

# show changes to the runtime templates introduced by the merge (generate-job.py may need updates)
git --no-pager log --oneline --stat ORIG_HEAD..HEAD -- scripts/runtime/*-example.yml
