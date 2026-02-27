#!/usr/bin/env bash
#
# generate_changelog.sh
#
# Generates CHANGELOG.md from git tags and commit history.
# Called automatically by the post-commit hook, or run manually:
#   ./scripts/generate_changelog.sh
#
# Commit messages starting with "chore:" or "wip:" are excluded.
# Tag your releases with: git tag -a v0.1.0 -m "First release"

set -euo pipefail

REPO_ROOT="$(git rev-parse --show-toplevel)"
CHANGELOG="$REPO_ROOT/CHANGELOG.md"

# Header
cat > "$CHANGELOG" << 'HEADER'
# Changelog

All notable changes to this project will be documented in this file.
This file is auto-generated from git history after each commit.

HEADER

# Get tags sorted by version (newest first)
mapfile -t TAG_ARRAY < <(git tag --sort=-version:refname 2>/dev/null || true)

# Unreleased commits (since last tag, or all if no tags)
if [ ${#TAG_ARRAY[@]} -gt 0 ]; then
    UNRELEASED=$(git log "${TAG_ARRAY[0]}"..HEAD --pretty=format:"- %s (%h)" \
        --no-merges | grep -vE "^- (chore|wip):" || true)
else
    UNRELEASED=$(git log --pretty=format:"- %s (%h)" \
        --no-merges | grep -vE "^- (chore|wip):" || true)
fi

if [ -n "$UNRELEASED" ]; then
    echo "## Unreleased" >> "$CHANGELOG"
    echo "" >> "$CHANGELOG"
    echo "$UNRELEASED" >> "$CHANGELOG"
    echo "" >> "$CHANGELOG"
fi

# Tagged releases — iterate newest to oldest
for i in "${!TAG_ARRAY[@]}"; do
    TAG="${TAG_ARRAY[$i]}"
    TAG_DATE=$(git log -1 --format="%as" "$TAG")
    echo "## $TAG — $TAG_DATE" >> "$CHANGELOG"
    echo "" >> "$CHANGELOG"

    # Next index is the older tag
    NEXT_INDEX=$((i + 1))

    if [ $NEXT_INDEX -lt ${#TAG_ARRAY[@]} ]; then
        # Commits between the older tag and this tag
        OLDER_TAG="${TAG_ARRAY[$NEXT_INDEX]}"
        COMMITS=$(git log "$OLDER_TAG".."$TAG" --pretty=format:"- %s (%h)" \
            --no-merges | grep -vE "^- (chore|wip):" || true)
    else
        # Oldest tag — show all commits up to and including this tag
        COMMITS=$(git log "$TAG" --pretty=format:"- %s (%h)" \
            --no-merges | grep -vE "^- (chore|wip):" || true)
    fi

    if [ -n "$COMMITS" ]; then
        echo "$COMMITS" >> "$CHANGELOG"
    else
        echo "- No notable changes" >> "$CHANGELOG"
    fi

    echo "" >> "$CHANGELOG"
done
