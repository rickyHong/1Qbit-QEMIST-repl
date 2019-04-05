# QEMIST documentation

The documentation is built with `sphinx` and `sphinx-apidoc`. To build the
documentation run
```bash
sphinx-apidoc  -o refman/ .. ../setup.py ../qemist/tests/
```
in this directory to generate the reference manual, the run
```bash
make html
```
to build the HTML documentation. The output files will be found in `_build/html`.
