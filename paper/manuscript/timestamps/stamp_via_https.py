"""
Direct HTTPS-API approach to OpenTimestamps anchoring.

The OpenTimestamps Python *client* (`ots.exe`) requires OpenSSL DLL on Windows
which does not load cleanly. The underlying calendar protocol is just an HTTPS
POST of the SHA256 digest. This script bypasses the client and posts the four
commit-hash digests directly to three public OpenTimestamps calendars
(alice / bob / finney).

For each commit hash:
  1. Compute SHA256(commit_hash_full + newline) to match the file digest
  2. POST that 32-byte digest to calendar.opentimestamps.org/digest endpoints
  3. Save the calendar's response as `<short>.ots.proof` for archival

Verification later:
  ots verify <file>.ots.proof  (after Bitcoin block confirmation, ~24 hours)
"""

import hashlib
import sys
import urllib.request
from pathlib import Path

CALENDARS = [
    "https://alice.btc.calendar.opentimestamps.org",
    "https://bob.btc.calendar.opentimestamps.org",
    "https://finney.calendar.eternitywall.com",
]

COMMITS = ["42168da", "44d7f5f", "76a4f77", "50c4687"]

HERE = Path(__file__).parent


def sha256_of_file(path: Path) -> bytes:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.digest()


def post_to_calendar(calendar_url: str, digest: bytes) -> bytes:
    url = f"{calendar_url}/digest"
    req = urllib.request.Request(
        url, data=digest, method="POST",
        headers={
            "Accept": "application/vnd.opentimestamps.v1",
            "User-Agent": "envmeta-paper3-anchor",
            "Content-Type": "application/x-www-form-urlencoded",
        },
    )
    with urllib.request.urlopen(req, timeout=30) as resp:
        return resp.read()


def main():
    rows = []
    for short in COMMITS:
        hash_file = HERE / f"{short}.commit_hash.txt"
        if not hash_file.exists():
            print(f"[SKIP] {hash_file.name} missing")
            continue

        digest = sha256_of_file(hash_file)
        digest_hex = digest.hex()
        print(f"\n=== {short} ===")
        print(f"  file: {hash_file.name}")
        print(f"  full commit: {hash_file.read_text().strip()}")
        print(f"  SHA256 digest: {digest_hex}")

        cal_responses = {}
        for cal_url in CALENDARS:
            cal_name = cal_url.split("//")[1].split(".")[0]
            try:
                resp = post_to_calendar(cal_url, digest)
                cal_responses[cal_name] = resp
                proof_path = HERE / f"{short}.{cal_name}.ots-response.bin"
                proof_path.write_bytes(resp)
                print(f"  [OK] {cal_name}: {len(resp)} bytes -> {proof_path.name}")
            except Exception as e:
                print(f"  [FAIL] {cal_name}: {e}")

        rows.append({
            "commit_short": short,
            "commit_full": hash_file.read_text().strip(),
            "digest_sha256": digest_hex,
            "calendars_anchored": list(cal_responses.keys()),
        })

    # Write summary
    summary_path = HERE / "ANCHOR_SUMMARY.md"
    lines = [
        "# OpenTimestamps anchor summary",
        "",
        f"Date: {__import__('datetime').date.today().isoformat()}",
        f"Calendars: {', '.join(c.split('//')[1].split('.')[0] for c in CALENDARS)}",
        "",
        "| commit (short) | commit (full) | SHA256 digest | calendars anchored |",
        "|---|---|---|---|",
    ]
    for r in rows:
        lines.append(
            f"| `{r['commit_short']}` | `{r['commit_full']}` | "
            f"`{r['digest_sha256']}` | {', '.join(r['calendars_anchored'])} |"
        )
    lines += [
        "",
        "## Verification",
        "",
        "Each `<short>.<calendar>.ots-response.bin` file contains the calendar's",
        "incomplete (pending) timestamp proof. After Bitcoin block confirmation",
        "(typically 1-6 hours, occasionally up to 24 hours), the proof can be",
        "upgraded to a full SHA256 → Bitcoin-block-merkle-root chain via:",
        "",
        "```",
        "ots upgrade <commit-short>.<calendar>.ots-response.bin",
        "ots verify <commit-hash-file>.ots.proof",
        "```",
        "",
        "(or via the JavaScript verifier at https://opentimestamps.org/)",
        "",
        "## Reproducibility note",
        "",
        "This script bypasses the standard `ots-client` CLI because the",
        "Windows-build OpenSSL library load fails (python-bitcoinlib issue).",
        "The HTTPS POST protocol is identical to what the CLI does internally.",
    ]
    summary_path.write_text("\n".join(lines), encoding="utf-8")
    print(f"\nSummary written to: {summary_path.name}")


if __name__ == "__main__":
    main()
