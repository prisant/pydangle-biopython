#!/bin/bash
# jsonl2txt.sh — Convert pydangle JSONL output to colon-separated text
#
# Produces output matching pydangle's native CSV format:
#   file:model:chain:resnum:ins:resname:field1:field2:...
#
# Handles any set of measurement fields dynamically. Null values
# render as __?__ (matching pydangle's CSV convention). Resnum is
# right-justified to 4 characters.
#
# Skips metadata records (lines with "_meta" field) if present.
#
# Usage:
#   pydangle-biopython -o jsonl -c "phi; psi; omega" file.pdb | ./jsonl2txt.sh
#   ./jsonl2txt.sh input.jsonl
#   ./jsonl2txt.sh input.jsonl output.txt
#
# Prerequisites:
#   jq (https://jqlang.github.io/jq/)

set -euo pipefail

JQ_FILTER='
  select(has("_meta") | not) |
  [
    to_entries[] |
    if .key == "resnum" then
      (("    " + (.value|tostring))[-4:])
    elif .value == null then
      "__?__"
    else
      (.value | tostring)
    end
  ] | join(":")
'

# Check if data is coming from a pipe (stdin)
if [ -p /dev/stdin ]; then
    jq -r "$JQ_FILTER"
    exit 0
fi

if [ "$#" -lt 1 ] || [ "$#" -gt 2 ]; then
    echo "Usage: $0 <input.jsonl> [output.txt]" >&2
    echo "       or pipe data into this script." >&2
    exit 1
fi

INPUT_FILE=$1

if [ "$#" -eq 1 ]; then
    jq -r "$JQ_FILTER" "$INPUT_FILE"
else
    OUTPUT_FILE=$2
    jq -r "$JQ_FILTER" "$INPUT_FILE" > "$OUTPUT_FILE"
    echo "Converted: $INPUT_FILE -> $OUTPUT_FILE" >&2
fi
