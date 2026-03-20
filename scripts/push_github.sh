#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

REMOTE="${1:-origin}"
BRANCH="${2:-$(git rev-parse --abbrev-ref HEAD)}"
PUSH_TAGS="${PUSH_TAGS:-0}"

REMOTE_URL="$(git remote get-url "$REMOTE")"

echo "Repository: $ROOT"
echo "Remote:     $REMOTE"
echo "URL:        $REMOTE_URL"
echo "Branch:     $BRANCH"
echo
git status --short --branch
echo

if [[ "$REMOTE_URL" == https://github.com/* ]]; then
  cat <<'EOF'
Note:
  This remote uses HTTPS. GitHub no longer accepts account passwords for push.
  Use either:
    1. a GitHub PAT, or
    2. an SSH remote such as git@github.com:<owner>/<repo>.git
EOF
  echo
fi

git push -u "$REMOTE" "$BRANCH"

if [[ "$PUSH_TAGS" == "1" ]]; then
  git push "$REMOTE" --tags
fi
