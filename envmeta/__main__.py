"""EnvMeta CLI dispatcher.

用法：
    python -m envmeta kb-build --elements arsenic,nitrogen ...

子命令：
    kb-build      Build elements.json from KEGG snapshot + seed KO list
"""
from __future__ import annotations

import sys


def _kb_build(argv: list[str]) -> int:
    from envmeta.tools.kb_builder import main as kb_main
    return kb_main(argv)


_COMMANDS = {
    "kb-build": _kb_build,
}


def main() -> int:
    argv = sys.argv[1:]
    if not argv or argv[0] in ("-h", "--help"):
        print("envmeta CLI - subcommands:")
        for cmd in _COMMANDS:
            print(f"  {cmd}")
        print("\nUse: python -m envmeta <subcommand> --help")
        return 0
    cmd = argv[0]
    if cmd not in _COMMANDS:
        print(f"Unknown subcommand: {cmd}", file=sys.stderr)
        print(f"Available: {', '.join(_COMMANDS)}", file=sys.stderr)
        return 2
    return _COMMANDS[cmd](argv[1:])


if __name__ == "__main__":
    raise SystemExit(main())
