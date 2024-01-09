This documentation was generated using Sphinx. The following steps can be taken to do so:

```
sudo apt-get install python3-sphinx
pip install sphinx_rtd_theme

mkdir docs
cd docs
sphinx-quickstart --no-batchfile --sep --ext-autodoc --ext-viewcode --ext-todo
```

Add the following in `source/conf.py`:
```
import os
import sys
sys.path.insert(0, os.path.abspath('../../'))
```
Also set `html_theme = 'sphinx_rtd_theme'`.

```
cd ..
sphinx-apidoc -o docs/source . setup.py # excludes setup.py
```

Include `modules` in `docs/source/index.rst`

```
cd docs
make html
```
