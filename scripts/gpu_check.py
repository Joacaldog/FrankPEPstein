#!/usr/bin/env python3
"""Report what the GPU backends would do on this machine, without running the pipeline.

Run it on the machine you intend to use:

    python3 scripts/gpu_check.py
    FRANKPEP_RESCORE_BACKEND=consensus python3 scripts/gpu_check.py

The second form shows whether an opt-in would actually take effect or fall back to CPU, which
is the question worth answering before launching a long run.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from gpu_backends import describe_environment  # noqa: E402


def main():
    for line in describe_environment():
        print(line)


if __name__ == "__main__":
    main()
