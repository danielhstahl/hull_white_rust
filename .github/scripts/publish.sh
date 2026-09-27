#!/usr/bin/env bash
#
# Publish this crate to crates.io, with the tolerated failure case narrowed.
#
# The workflow step this replaced was
#
#     cargo publish --token ${{ secrets.CARGO_TOKEN }} --allow-dirty || true
#
# and `|| true` made *every* publish failure invisible: an empty or expired token, a
# rejected manifest, a version that was already uploaded, a compile error -- the
# deploy job reported green regardless, so "we never shipped 0.9.1" was
# indistinguishable from "we shipped it".  Here everything is fatal except one
# case, stated below, which is surfaced as a warning rather than passed silently.
#
# Usage:
#   bash .github/scripts/publish.sh              # publish (needs CARGO_TOKEN)
#   bash .github/scripts/publish.sh --dry-run    # every check + `cargo publish --dry-run`, no upload
#
# Env:
#   CARGO_TOKEN  crates.io API token.  Required unless --dry-run.

# `-e`: an unexpected failure ends the run red.  `-u`: an unset variable is a bug,
# not an empty string.  `-o pipefail`: a failure inside `foo | bar` is not lost.
set -euo pipefail

dry_run=0
if [ "${1:-}" = "--dry-run" ]; then
  dry_run=1
fi

# --- credential -------------------------------------------------------------
# An empty secret used to be swallowed along with everything else.  Say so.
if [ "$dry_run" -eq 0 ] && [ -z "${CARGO_TOKEN:-}" ]; then
  echo "::error::secrets.CARGO_TOKEN is not set; cannot publish." >&2
  exit 1
fi

# --- what is cargo about to package? ----------------------------------------
# Unified Package Identifier: `...#<name>@<version>`.  Splitting that avoids a `jq`
# dependency in the deploy path -- cargo already produced the structured value.
pkgid=$(cargo pkgid)
namever=${pkgid##*#}
pkg=${namever%@*}
version=${namever##*@}
if [ -z "$pkg" ] || [ -z "$version" ] || [ "$pkg" = "$namever" ]; then
  echo "::error::could not read the package name/version from 'cargo pkgid': ${pkgid}" >&2
  exit 1
fi
echo "Package: ${pkg} ${version}"

# --- the tree we are packaging ----------------------------------------------
# `--allow-dirty` is retained below: generated artefacts (`target/`, and the
# gitignored `Cargo.lock`) sit next to the manifest, and cargo should not refuse to
# ship over that.  What it is *not* allowed to hide is a modified tracked file, so
# assert the tracked tree is clean before relying on the flag.  If this step fires,
# some earlier step in the workflow wrote to source control instead of to target/.
tracked_dirty=$(git status --porcelain --untracked-files=no)
if [ -n "$tracked_dirty" ]; then
  echo "::error::tracked files were modified during the build; refusing to publish:" >&2
  printf '%s\n' "$tracked_dirty" >&2
  exit 1
fi

# --- the one tolerated case -------------------------------------------------
# Already on crates.io at this exact version => this push did not bump
# `package.version`.  crates.io rejects the re-upload and retrying will never make
# it pass, so skip -- but loudly, as a warning annotation on the run.  "Skipped,
# not shipped" has to be a state someone can see; it must not be the default.
#
# crates.io returns 404 for an unknown (name, version) pair and requires a
# User-Agent on API calls.
status=$(curl -sS -o /dev/null -w '%{http_code}' \
  -A "hull-white-ci (https://github.com/danielhstahl/hull_white_rust)" \
  "https://crates.io/api/v1/crates/${pkg}/${version}")

if [ "$status" = "200" ]; then
  echo "::warning::${pkg} ${version} is already on crates.io -- skipping publish. Bump package.version to ship this commit."
  exit 0
fi

echo "${pkg} ${version} is not on crates.io (HTTP ${status}); proceeding."

# --- go ---------------------------------------------------------------------
if [ "$dry_run" -eq 1 ]; then
  # Verifies that the package builds from the packaged tarball and that the manifest
  # is publishable, without uploading and without a token.
  echo "Dry run: cargo publish --dry-run --allow-dirty (nothing is uploaded)"
  cargo publish --dry-run --allow-dirty
  echo "::notice::Dry run complete for ${pkg} ${version}; nothing uploaded."
  exit 0
fi

cargo publish --token "${CARGO_TOKEN}" --allow-dirty
echo "::notice::Published ${pkg} ${version} to crates.io."
