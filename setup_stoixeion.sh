#!/usr/bin/env bash
# Prepare the exact upstream Stoixeion revision used by run_stoixeion_miniscope.m.
set -euo pipefail

repo_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
target=${1:-"${repo_dir}/external/Stoixeion"}
upstream=https://github.com/luiscareid/Stoixeion.git
revision=be6f547cfb6d2b49c1c0ec3b6331d9ea98f08491

if [[ ! -d "${target}/.git" ]]; then
  mkdir -p "$(dirname "${target}")"
  git clone "${upstream}" "${target}"
fi

if [[ $(git -C "${target}" rev-parse HEAD) != "${revision}" ]]; then
  if [[ -n $(git -C "${target}" status --porcelain) ]]; then
    echo "Stoixeion has local changes; refusing to switch revisions: ${target}" >&2
    exit 1
  fi
  git -C "${target}" fetch origin "${revision}"
  git -C "${target}" checkout --detach "${revision}"
fi

patch_file="${repo_dir}/stoixeion_exports.patch"
if git -C "${target}" apply --reverse --check "${patch_file}"; then
  echo "Stoixeion patch is already applied: ${target}"
elif git -C "${target}" apply --check "${patch_file}"; then
  git -C "${target}" apply "${patch_file}"
  echo "Stoixeion is ready: ${target}"
else
  echo "Stoixeion does not match the pinned revision and patch: ${target}" >&2
  exit 1
fi
