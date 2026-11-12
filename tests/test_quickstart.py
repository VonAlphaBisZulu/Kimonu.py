"""The README example, examples/quickstart.py: it runs, prints what the README shows, and is the README's code."""
import os
import re
import subprocess
import sys

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir)
SCRIPT = os.path.join(ROOT, "examples", "quickstart.py")


def _readme_blocks():
    """The first python block of README.md and the plain block that follows it, its output."""
    text = open(os.path.join(ROOT, "README.md"), encoding="utf-8").read()
    code, output = re.search(r"```python\n(.*?)```\s*```\n(.*?)```", text, re.S).groups()
    return code, output


def test_quickstart_is_the_readme_example():
    code, _ = _readme_blocks()
    script = open(SCRIPT, encoding="utf-8").read()
    assert script.split('"""\n', 1)[1] == code


def test_quickstart_prints_the_readme_output():
    _, output = _readme_blocks()
    env = dict(os.environ, PYTHONPATH=os.pathsep.join([os.path.join(ROOT, "src"), os.environ.get("PYTHONPATH", "")]))
    run = subprocess.run([sys.executable, SCRIPT], capture_output=True, text=True, env=env, timeout=600)
    assert run.returncode == 0, run.stderr
    assert output in run.stdout                         # a commercial solver may add lines about its license
