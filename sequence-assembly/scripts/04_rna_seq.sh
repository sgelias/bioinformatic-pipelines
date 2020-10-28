#!/bin/bash
#
# NAME: RNAseq analysis pipeline
# PLAT: Linux
# DESCRIPTION: ...

# Parametrization example from:
# https://stackoverflow.com/questions/7069682/how-to-get-arguments-with-flags-in-bash

while test $# -gt 0; do
    case "$1" in
    -h | --help)
        echo "$package - attempt to capture frames"
        echo " "
        echo "$package [options] application [arguments]"
        echo " "
        echo "options:"
        echo "-h, --help                show brief help"
        echo "-a, --action=ACTION       specify an action to use"
        echo "-o, --output-dir=DIR      specify a directory to store output in"
        exit 0
        ;;
    -a)
        shift
        if test $# -gt 0; then
            export PROCESS=$1
        else
            echo "no process specified"
            exit 1
        fi
        shift
        ;;
    --action*)
        export PROCESS=$(echo $1 | sed -e 's/^[^=]*=//g')
        shift
        ;;
    -o)
        shift
        if test $# -gt 0; then
            export OUTPUT=$1
        else
            echo "no output dir specified"
            exit 1
        fi
        shift
        ;;
    --output-dir*)
        export OUTPUT=$(echo $1 | sed -e 's/^[^=]*=//g')
        shift
        ;;
    *)
        break
        ;;
    esac
done
