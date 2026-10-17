import os
import sys
import time

import pytest

from qcengine import util
from qcengine.exceptions import InputError
from qcengine.testing import schema_versions2


def test_model_wrapper(schema_versions2):
    models, _, _ = schema_versions2

    with pytest.raises(InputError):
        util.model_wrapper({"bad": "yup"}, models.AtomicInput)


def test_compute_wrapper_capture():

    oldout = id(sys.stdout)
    with util.compute_wrapper(capture_output=True):

        assert id(sys.stdout) != id(oldout)

    assert id(sys.stdout) == oldout


def test_compute_wrapper_capture_exception():

    oldout = id(sys.stdout)
    with util.compute_wrapper(capture_output=True) as metadata:

        assert id(sys.stdout) != id(oldout)

        raise KeyError("Hello there!")

    assert id(sys.stdout) == oldout

    assert metadata["success"] is False
    assert metadata["error_type"] == "unknown_error"


def test_terminate():

    t = time.time()
    with util.popen(["sleep", "30"]) as proc:

        util.terminate_process(proc["proc"])

    assert (time.time() - t) < 1


def test_tmpdir():

    with util.temporary_directory(child="this") as tmpdir:
        assert str(tmpdir).split(os.path.sep)[-1] == "this"

    with util.temporary_directory() as parentdir:
        with util.temporary_directory(parent=parentdir, child="this") as tmpdir:
            assert str(tmpdir).split(os.path.sep)[-1] == "this"
            assert str(tmpdir).rsplit(os.path.sep, 1)[0] == str(parentdir)

    with util.temporary_directory(suffix="this") as tmpdir:
        assert str(tmpdir).split(os.path.sep)[-1].endswith("this")


def test_disk_files():

    infiles = {"thing1": "hello", "thing2": "world", "other": "everyone"}
    outfiles = {"thing*": None, "other": None}
    with util.temporary_directory(suffix="this") as tmpdir:
        with util.disk_files(infiles=infiles, outfiles=outfiles, cwd=tmpdir):
            pass

    assert outfiles.keys() == {"thing*", "other"}
    assert outfiles["thing*"]["thing1"] == "hello"
    assert outfiles["other"] == "everyone"


def test_popen_tee_output(capsys):
    # Test without passing
    with util.popen(["echo", "hello"]) as proc:
        proc["proc"].wait()
    assert proc["stdout"].strip() == "hello"

    # Test with passing
    with util.popen(["echo", "hello"], pass_output_forward=True) as proc:
        proc["proc"].wait()
    assert proc["stdout"] == "hello\n"
    captured = capsys.readouterr()
    assert captured.out == "hello\n"


def test_popen_non_utf8_output():
    """Test that popen gracefully handles non-UTF-8 bytes in subprocess output.

    This is a regression test for Windows encoding issues where subprocess output
    may contain non-UTF-8 bytes (e.g., raw 0xa7 byte for section symbol instead of
    UTF-8 encoded 0xc2 0xa7). The popen function should use errors='replace' to
    handle these gracefully without raising UnicodeDecodeError.
    """
    # Create a Python script that writes non-UTF-8 bytes to stdout and stderr
    script = r"""
import sys
# Write valid UTF-8
sys.stdout.write("Valid UTF-8: ")
# Write invalid byte (raw 0xa7, which is not valid UTF-8 on its own)
sys.stdout.buffer.write(b"\xa7")
sys.stdout.write(" end\n")

# Also test stderr
sys.stderr.write("Stderr with invalid: ")
sys.stderr.buffer.write(b"\xa7")
sys.stderr.write("\n")
"""

    # Execute the script through popen - it should NOT raise UnicodeDecodeError
    with util.popen([sys.executable, "-c", script]) as proc:
        proc["proc"].wait()

    # Verify that the function completed successfully
    assert proc["proc"].returncode == 0

    # Verify that stdout and stderr are valid strings (not crashed with UnicodeDecodeError)
    assert isinstance(proc["stdout"], str)
    assert isinstance(proc["stderr"], str)

    # Verify that valid UTF-8 portions are preserved
    assert "Valid UTF-8:" in proc["stdout"]
    assert "end" in proc["stdout"]
    assert "Stderr with invalid:" in proc["stderr"]

    # Verify that invalid bytes are replaced with the replacement character
    # (U+FFFD appears as '?' when displayed, but in the string it's the actual replacement char)
    assert "\ufffd" in proc["stdout"], "Non-UTF-8 byte should be replaced with U+FFFD"
    assert "\ufffd" in proc["stderr"], "Non-UTF-8 byte in stderr should be replaced with U+FFFD"


def test_popen_non_utf8_output_section_symbol():
    """Test handling of section symbol (§, U+00A7) as raw byte 0xa7.

    This specific test addresses the Windows issue where Psi4's n-body output
    contains section symbols that may be written as raw bytes instead of UTF-8.
    """
    # Create a script that outputs the section symbol as raw byte (as Windows might do)
    script = r"""
import sys
sys.stdout.write("N-Body: ")
# Write section symbol as raw byte instead of UTF-8 (0xc2 0xa7)
sys.stdout.buffer.write(b"\xa7")
sys.stdout.write("A_(2)@(2) total energy\n")
"""

    # Execute - should not crash
    with util.popen([sys.executable, "-c", script]) as proc:
        proc["proc"].wait()

    assert proc["proc"].returncode == 0
    assert isinstance(proc["stdout"], str)

    # The output should be recoverable and contain the replacement character
    assert "N-Body:" in proc["stdout"]
    assert "total energy" in proc["stdout"]
    assert "\ufffd" in proc["stdout"], "Section symbol as raw byte should be replaced"
