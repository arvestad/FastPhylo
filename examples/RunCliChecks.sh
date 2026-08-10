#!/bin/bash
# CLI-level regression checks for apps migrated off gengetopt onto
# CLI11 (see gengetopt_migration_plan.md). RunExamples.sh only
# exercises real flag combinations that produce output to diff -
# --help/--version content and error-path exit codes aren't covered
# by anything else, so each newly migrated app gets a check_app call
# here. fnj is the first (Phase B).

set -u

if [ -z "${FASTPHYLOPATH:-}" ]; then
    echo "Set FASTPHYLOPATH to the directory with the project binaries!"
    exit 1
fi

export PATH="${FASTPHYLOPATH}:$PATH"

failures=0

check() {
    local description="$1"
    shift
    if "$@" >/dev/null 2>&1; then
        echo "  OK: $description"
    else
        echo "  FAIL: $description"
        failures=$((failures + 1))
    fi
}

# $1 = app name, remaining args = long option names ("outfile", not
# "--outfile") expected to appear in --help output. --version is
# checked against "$app <anything>" rather than an exact string -
# PACKAGE_VERSION is no longer a fixed number to hardcode here (it is
# git-describe-derived for an ordinary build, or tag-derived for a
# release build - see github_actions_release_builds_plan.md's Phase E),
# just a format/plumbing check that the app name made it through.
check_app() {
    local app="$1"
    shift

    echo "== $app =="

    check "$app --version exits 0" "$app" --version
    local actual_version
    actual_version="$("$app" --version 2>&1)"
    case "$actual_version" in
        "$app "*)
            echo "  OK: $app --version starts with '$app '"
            ;;
        *)
            echo "  FAIL: $app --version printed '$actual_version', expected it to start with '$app '"
            failures=$((failures + 1))
            ;;
    esac

    check "$app --help exits 0" "$app" --help
    for opt in "$@"; do
        if "$app" --help 2>&1 | grep -q -- "--$opt"; then
            echo "  OK: $app --help mentions --$opt"
        else
            echo "  FAIL: $app --help does not mention --$opt"
            failures=$((failures + 1))
        fi
    done

    if ! "$app" -I bogus </dev/null >/dev/null 2>&1; then
        echo "  OK: $app -I bogus is rejected (nonzero exit)"
    else
        echo "  FAIL: $app -I bogus was accepted"
        failures=$((failures + 1))
    fi

    if ! "$app" --this-flag-does-not-exist </dev/null >/dev/null 2>&1; then
        echo "  OK: $app --this-flag-does-not-exist is rejected (nonzero exit)"
    else
        echo "  FAIL: $app --this-flag-does-not-exist was accepted"
        failures=$((failures + 1))
    fi
}

check_app fnj outfile input-format output-format print-counts \
    analyze-run-number method dm-per-run number-of-runs bootstraps \
    print-relaxng-input print-relaxng-output

check_app fastdist outfile input-format memory-efficient \
    output-format distance-function bootstraps no-incl-orig seed \
    no-ambiguities no-ambig-resolve no-transprob ambiguity-frequency-model \
    tstvratio pyrtvratio no-tstvratio fixfactor number-of-runs \
    print-relaxng-input print-relaxng-output

check_app fastprot outfile input-format memory-efficient \
    output-format bootstraps no-incl-orig seed distance-function model-file \
    remove-indels maximum-likelihood sd pfam speed print-relaxng-input \
    print-relaxng-output

# fastprot_mpi is out of scope for this script - only built with
# -DBUILD_WITH_MPI=ON, which needs an MPI installation this repo's CI
# doesn't provide (see gengetopt_migration_plan.md's Phase D note).

echo
if [ "$failures" -eq 0 ]; then
    echo "All CLI checks passed."
    exit 0
else
    echo "$failures CLI check(s) failed."
    exit 1
fi
