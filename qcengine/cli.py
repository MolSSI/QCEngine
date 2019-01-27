"""
Provides a CLI for QCEngine
"""

import argparse
import json

from . import compute

__all__ = ["main"]


def parse_args():
    parser = argparse.ArgumentParser(description='A CLI for the QCEngine.')

    parser.add_argument("program", type=str, default=None, help="The program to call")
    parser.add_argument("data", type=str, default=None, help="A JSON blob of the task description.")

    return vars(parser.parse_args())


def main(args=None):

    # Grab CLI args if not present
    if args is None:
        args = parse_args()

    ret = compute(json.loads(args["data"]), args["program"])
    print(json.dumps(ret))


if __name__ == '__main__':
    main()
