#!/bin/bash
# Create tgz for each year for multiple subdirectories (restart, log, ...)
#
# USAGE:
#   ./pack_year.sh BASEDIR START END [--delete]
#
# EXAMPLE:
#   ./pack_year.sh /data 100 100 --delete

BASEDIR=$1
START=$2
END=$3
CLEANUP=$4    # --delete

PAD=3

# list of folders
DATASETS=("restart" "log")

if [ -z "$BASEDIR" ] || [ -z "$START" ] || [ -z "$END" ]; then
    echo "Usage: $0 BASEDIR START END [--delete]"
    exit 1
fi

if [ ! -d "$BASEDIR" ]; then
    echo "Error: BASEDIR does not exist: $BASEDIR"
    exit 1
fi

for (( i=START; i<=END; i++ )); do
    YEAR=$(printf "%0${PAD}d" $i)
    STATUS=0
    FOUND=0

    for DS in "${DATASETS[@]}"; do
        SRC="${BASEDIR}/${DS}/${YEAR}"
        TAR="${BASEDIR}/${DS}/${YEAR}.tgz"

        if [ ! -d "$SRC" ]; then
            continue
        fi

        FOUND=1

        if [ -f "$TAR" ]; then
            echo "Skip ${DS}/${YEAR} (tgz exists)"
            continue
        fi

        echo "Compress ${DS}/${YEAR} → ${TAR}"
        tar -czf "$TAR" -C "${BASEDIR}/${DS}" "$YEAR" || STATUS=1
    done

    if [ $FOUND -eq 0 ]; then
        echo "Skip $YEAR (not found)"
        continue
    fi

    if [ $STATUS -eq 0 ] && [ "$CLEANUP" = "--delete" ]; then
        for DS in "${DATASETS[@]}"; do
            rm -rf "${BASEDIR}/${DS}/${YEAR}"
        done
        echo "$YEAR deleted"
    fi
done

echo "Compression done"
