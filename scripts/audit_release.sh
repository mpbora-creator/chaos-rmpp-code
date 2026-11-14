#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

required=(README.md LICENSE CITATION.cff .zenodo.json Project.toml Manifest.toml)
for path in "${required[@]}"; do
  test -s "$path" || { printf 'Missing required file: %s\n' "$path" >&2; exit 1; }
done

if find . -type f -size +95M -not -path './.git/*' | grep -q .; then
  printf 'Files larger than 95 MB must not be committed to GitHub:\n' >&2
  find . -type f -size +95M -not -path './.git/*' >&2
  exit 1
fi

if rg --hidden -n '/Users/|/home/' \
  --glob '!build/**' --glob '!diagnostics/**' --glob '!fig/**' \
  --glob '!.git/**' --glob '!scripts/audit_release.sh' .; then
  printf 'Remove the machine-specific paths reported above before release.\n' >&2
  exit 1
fi

if rg -n 'OWNER/REPOSITORY|zenodo\.XXXXXXX' CITATION.cff README.md; then
  if test "${ALLOW_RELEASE_PLACEHOLDERS:-0}" != 1; then
    printf 'Replace the repository/DOI placeholders reported above before release.\n' >&2
    exit 1
  fi
  printf 'Warning: release metadata still contains approved preparation placeholders.\n'
fi

printf 'Release audit passed.\n'
