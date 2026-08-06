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

# $1 = app name, $2 = expected exact "--version" output, remaining
# args = long option names ("outfile", not "--outfile") expected to
# appear in --help output.
check_app() {
    local app="$1"
    local expected_version="$2"
    shift 2

    echo "== $app =="

    check "$app --version exits 0" "$app" --version
    if [ "$("$app" --version 2>&1)" = "$expected_version" ]; then
        echo "  OK: $app --version prints '$expected_version'"
    else
        echo "  FAIL: $app --version printed '$("$app" --version 2>&1)', expected '$expected_version'"
        failures=$((failures + 1))
    fi

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

check_app fnj "fnj 1.0.10" outfile input-format output-format print-counts \
    analyze-run-number method dm-per-run number-of-runs bootstraps \
    print-relaxng-input print-relaxng-output

check_app fastdist "fastdist 1.0.10" outfile input-format memory-efficient \
    output-format distance-function bootstraps no-incl-orig seed \
    no-ambiguities no-ambig-resolve no-transprob ambiguity-frequency-model \
    tstvratio pyrtvratio no-tstvratio fixfactor number-of-runs \
    print-relaxng-input print-relaxng-output

echo
if [ "$failures" -eq 0 ]; then
    echo "All CLI checks passed."
    exit 0
else
    echo "$failures CLI check(s) failed."
    exit 1
fi
