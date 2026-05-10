# OpenTimestamps anchor summary

Date: 2026-05-10
Calendars: alice, bob, finney

| commit (short) | commit (full) | SHA256 digest | calendars anchored |
|---|---|---|---|
| `42168da` | `42168da153f9c64a51966d1c93a9f596a80c834c` | `fadde64521e5224c7bc943f4caf5aee59aa33afc8a7ec0a3bf2cd8609db569f2` | alice, bob, finney |
| `44d7f5f` | `44d7f5f14769b7693650111305fe7d30ccfe2c27` | `4667511a1f9c909806faa72cc603df6c98ea7d5addbf087f2aa81cc8d5fb8a8d` | alice, bob, finney |
| `76a4f77` | `76a4f775b74bd5134e628c0e3ec6077c04cc40ae` | `c4117b4250d987043ea555e9fa80050b99859cfc79b8064e117b24054efe8a83` | alice, bob, finney |
| `50c4687` | `50c4687031277cff9fa8957b0429e8d85d1a8149` | `bc2a82a09edf7dbe18048f900b734ce68e0683a8baf52c5bafe22106bdf0b2bb` | alice, bob, finney |

## Verification

Each `<short>.<calendar>.ots-response.bin` file contains the calendar's
incomplete (pending) timestamp proof. After Bitcoin block confirmation
(typically 1-6 hours, occasionally up to 24 hours), the proof can be
upgraded to a full SHA256 → Bitcoin-block-merkle-root chain via:

```
ots upgrade <commit-short>.<calendar>.ots-response.bin
ots verify <commit-hash-file>.ots.proof
```

(or via the JavaScript verifier at https://opentimestamps.org/)

## Reproducibility note

This script bypasses the standard `ots-client` CLI because the
Windows-build OpenSSL library load fails (python-bitcoinlib issue).
The HTTPS POST protocol is identical to what the CLI does internally.