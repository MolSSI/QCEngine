"""
Provides a CLI for QCEngine
"""

import argparse
import json
import sys

from . import compute, compute_procedure  # run and run-procedure
from . import __version__, list_all_procedures, list_all_programs, list_available_procedures, \
    list_available_programs, get_procedure, get_program  # info
from .config import global_repr  # info

__all__ = ["main"]


def parse_args():
    parser = argparse.ArgumentParser(description='A CLI for the QCEngine.')
    subparsers = parser.add_subparsers(dest="command")

    info = subparsers.add_parser('info', help="Print information about QCEngine setup, version, and environment.")
    info.add_argument('--version', action='store_true', help="Print version of qcengine and qcelemental.")
    info.add_argument('--programs', action='store_true', help="Print detected and supported programs.")
    info.add_argument('--procedures', action='store_true', help="Print detected and supported procedures.")
    info.add_argument('--config', action='store_true', help="Print host, compute, and job configuration.")
    info.add_argument('--all', action='store_true', help="Print all available information.")

    run = subparsers.add_parser('run', help="Run a program on a given task. Output is printed as a JSON blob.")
    run.add_argument('program', type=str, help="The program to run.")
    run.add_argument('data', type=str, help="A JSON blob of the task description. "
                                            "If '-', data will be read from STDIN.")

    run_procedure = subparsers.add_parser('run-procedure', help="Run a procedure on a given task. "
                                                                "Output is printed as a JSON blob.")
    run_procedure.add_argument('procedure', type=str, help="The procedure to run.")
    run_procedure.add_argument('data', type=str, help="A JSON blob of the task description. "
                                                      "If '-', data will be read from STDIN.")

    args = vars(parser.parse_args())
    if args["command"] is None:
        parser.print_help(sys.stderr)
        exit(1)

    # Print usage and exit if no info options are provided
    if args["command"] == "info" and not any([v for v in args.values() if v != "info"]):
        info.print_help()
        exit(1)

    return args

def info_cli(args):

    def info_version():
        print(">>> Version information:")
        print(f"QCEngine version:    {__version__}")
        try:
            import qcelemental
            print(f"QCElemental version: {qcelemental.__version__}")
        except ImportError as e:
            print(e)
        print()

    def info_programs():
        print(">>> Program information")
        all_progs = list_all_programs()
        avail_progs = list_available_programs()
        print("Available programs:")
        for prog_name in sorted(avail_progs):
            version = get_program(prog_name).get_version()
            if version is None:
                version = "???"
            print(f"{prog_name} v{version}")

        print()
        print("Other supported programs:")
        print(" ".join(sorted(all_progs - avail_progs)))
        print()

    def info_procedures():
        print(">>> Procedure information")
        all_procs = list_all_procedures()
        avail_procs = list_available_procedures()
        print("Available procedures:")
        for proc_name in sorted(avail_procs):
            version = get_procedure(proc_name).get_version()
            if version is None:
                version = "???"
            print(f"{proc_name} v{version}")

        print()
        print("Other supported procedures:")
        print(" ".join(sorted(all_procs - avail_procs)))
        print()

    if args["version"] or args["all"]:
        info_version()
    if args["programs"] or args["all"]:
        info_programs()
    if args["procedures"] or args["all"]:
        info_procedures()
    if args["config"] or args["all"]:
        print(">>> Configuration information")
        print()
        print(global_repr())


def main(args=None):
    # Grab CLI args if not present
    if args is None:
        args = parse_args()
    command = args.pop("command")
    if command == "info":
        info_cli(args)
    elif command == "run":
        ret = compute(json.loads(args["data"] if args["data"] != "-" else sys.stdin.read()), args["program"])
        print(ret.json())
    elif command == "run-procedure":
        ret = compute_procedure(json.loads(args["data"] if args["data"] != "-" else sys.stdin.read()), args["procedure"])
        print(ret.json())


if __name__ == '__main__':
    main()
